/*
   Copyright (c) 2009-2016, Jack Poulson.
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

#include "./Process.hpp"
#include "./LowerSolve/Forward.hpp"
#include "./LowerSolve/Backward.hpp"
#include "./LowerMultiply/Forward.hpp"
#include "./LowerMultiply/Backward.hpp"

namespace El {

// Prototypes
namespace ldl {

template<typename Field>
void ChangeFrontType
( Front<Field>& front, LDLFrontType type, bool recurse=true );
template<typename Field>
void DiagonalSolve
( const NodeInfo& info,
  const Front<Field>& front,
        MatrixNode<Field>& B );
template<typename Field>
void DiagonalScale
( const NodeInfo& info,
  const Front<Field>& front,
        MatrixNode<Field>& B );

} // namespace ldl

template<typename Field>
SparseLDLFactorization<Field>::SparseLDLFactorization()
{ }

template<typename Field>
void SparseLDLFactorization<Field>::Initialize
( const SparseMatrix<Field>& A,
        bool hermitian,
  const BisectCtrl& bisectCtrl )
{
    EL_DEBUG_CSE
    info_.reset( new ldl::NodeInfo );
    separator_.reset( new ldl::Separator );
    ldl::NestedDissection
    ( A.LockedGraph(), map_, *separator_, *info_, bisectCtrl );
    InvertMap( map_, inverseMap_ );
    front_.reset( new ldl::Front<Field>(A,map_,*info_,hermitian) );

    initialized_ = true;
    factored_ = false;
}

template<typename Field>
void SparseLDLFactorization<Field>::Initialize2DGridGraph
( Int gridDim0,
  Int gridDim1,
  const SparseMatrix<Field>& A,
        bool hermitian,
  const BisectCtrl& bisectCtrl )
{
    EL_DEBUG_CSE
    info_.reset( new ldl::NodeInfo );
    separator_.reset( new ldl::Separator );
    ldl::NaturalNestedDissection
    ( gridDim0, gridDim1, 1, A.LockedGraph(),
      map_, *separator_, *info_, bisectCtrl.cutoff );
    InvertMap( map_, inverseMap_ );
    front_.reset( new ldl::Front<Field>(A,map_,*info_,hermitian) );

    initialized_ = true;
    factored_ = false;
}

template<typename Field>
void SparseLDLFactorization<Field>::Initialize3DGridGraph
( Int gridDim0,
  Int gridDim1,
  Int gridDim2,
  const SparseMatrix<Field>& A,
        bool hermitian,
  const BisectCtrl& bisectCtrl )
{
    EL_DEBUG_CSE
    info_.reset( new ldl::NodeInfo );
    separator_.reset( new ldl::Separator );
    ldl::NaturalNestedDissection
    ( gridDim0, gridDim1, gridDim2, A.LockedGraph(),
      map_, *separator_, *info_, bisectCtrl.cutoff );
    InvertMap( map_, inverseMap_ );
    front_.reset( new ldl::Front<Field>(A,map_,*info_,hermitian) );

    initialized_ = true;
    factored_ = false;
}

template<typename Field>
void SparseLDLFactorization<Field>::Factor( LDLFrontType frontType )
{
    EL_DEBUG_CSE
    if( !initialized_ )
        LogicError("Must initialize before calling 'Factor()'");
    // We make use of the following rather than checking 'factored_' since it
    // is sometimes useful to directly manipulate the fronts.
    if( !Unfactored(front_->type) )
        LogicError("Fronts are already marked as factored");
    
    // Convert from 1D to 2D if necessary
    ChangeFrontType( SYMM_2D );
    
    // Perform the initial factorization
    ldl::Process( *info_, *front_, InitialFactorType(frontType) );
    factored_ = true;
    
    // Convert the fronts from the initial factorization to the requested form
    ChangeFrontType( frontType );
}

template<typename Field>
void SparseLDLFactorization<Field>::ChangeFrontType( LDLFrontType frontType )
{
    EL_DEBUG_CSE
    if( !initialized_ )
        LogicError("Must initialize before calling 'ChangeFrontType()'");
    ldl::ChangeFrontType( *front_, frontType );
}

template<typename Field>
void SparseLDLFactorization<Field>::ChangeNonzeroValues
( const SparseMatrix<Field>& ANew )
{
    EL_DEBUG_CSE
    if( !initialized_ )
        LogicError("Must initialize before calling 'ChangeNonzeroValues()'");
    front_->Pull( ANew, map_, *info_ );
    factored_ = false;
}

template<typename Field>
void SparseLDLFactorization<Field>::Solve( Matrix<Field>& B ) const
{
    EL_DEBUG_CSE
    if( !factored_ )
        LogicError("Must call Factor() before Solve()");
    // TODO(poulson): Reuse a single MatrixNode data structure?
    ldl::MatrixNode<Field> BNodal( inverseMap_, *info_, B );
    Solve( BNodal );
    BNodal.Push( inverseMap_, *info_, B );
}

template<typename Field>
void SparseLDLFactorization<Field>::Solve( ldl::MatrixNode<Field>& B ) const
{
    EL_DEBUG_CSE
    if( !factored_ )
        LogicError("Must call Factor() before Solve()");
    const Orientation orientation = front_->isHermitian ? ADJOINT : TRANSPOSE;
    if( BlockFactorization(front_->type) )
    {   
        // Solve against block diagonal factor, L D
        SolveAgainstL( NORMAL, B );
        // Solve against the (conjugate-)transpose of the block unit diagonal L
        SolveAgainstL( orientation, B );
    }
    else
    {   
        // Solve against unit diagonal L 
        SolveAgainstL( NORMAL, B );
        // Solve against diagonal
        SolveAgainstD( NORMAL, B );
        // Solve against the (conjugate-)transpose of the unit diagonal L
        SolveAgainstL( orientation, B );
    }
}

template<typename Field>
void SparseLDLFactorization<Field>::SolveWithIterativeRefinement
( const SparseMatrix<Field>& A,
        Matrix<Field>& B,
  const Base<Field>& minReductionFactor,
        Int maxRefineIts ) const
{
    EL_DEBUG_CSE
    if( !factored_ )
        LogicError("Must call Factor() before SolveWithIterativeRefinement()");
    // TODO(poulson): Generalize this implementation
    if( B.Width() > 1 )
        LogicError("Iterative Refinement currently only supported for one RHS");
    auto BOrig = B;

    // Compute the initial guess
    // =========================
    Matrix<Field> X( B );
    Solve( X );

    Int refineIt = 0;
    if( maxRefineIts > 0 )
    {
        Matrix<Field> dX, XCand;
        Multiply( NORMAL, Field(-1), A, X, Field(1), B );
        Base<Field> errorNorm = FrobeniusNorm( B );
        for( ; refineIt<maxRefineIts; ++refineIt )
        {
            // Compute the proposed update to the solution
            // -------------------------------------------
            dX = B;
            Solve( dX );
            XCand = X;
            XCand += dX;

            // If the proposed update lowers the residual, accept it
            // -----------------------------------------------------
            B = BOrig;
            Multiply( NORMAL, Field(-1), A, XCand, Field(1), B );
            Base<Field> newErrorNorm = FrobeniusNorm( B );
            if( minReductionFactor*newErrorNorm < errorNorm )
            {
                X = XCand;
                errorNorm = newErrorNorm;
            }
            else if( newErrorNorm < errorNorm )
            {
                X = XCand;
                errorNorm = newErrorNorm;
                break;
            }
            else
                break;
        }
    }
    // Store the final result
    // ======================
    B = X;
}

template<typename Field>
void SparseLDLFactorization<Field>::SolveAgainstL
( Orientation orientation, Matrix<Field>& B ) const
{
    EL_DEBUG_CSE
    if( !factored_ )
        LogicError("Must call Factor() before SolveAgainstL()");
    ldl::MatrixNode<Field> BNodal( inverseMap_, *info_, B );
    SolveAgainstL( orientation, BNodal );
    BNodal.Push( inverseMap_, *info_, B );
}

template<typename Field>
void SparseLDLFactorization<Field>::SolveAgainstL
( Orientation orientation, ldl::MatrixNode<Field>& B ) const
{
    EL_DEBUG_CSE
    if( !factored_ )
        LogicError("Must call Factor() before SolveAgainstL()");
    if( orientation == NORMAL )
        ldl::LowerForwardSolve( *info_, *front_, B );
    else
        ldl::LowerBackwardSolve( *info_, *front_, B, orientation==ADJOINT );
}

template<typename Field>
void SparseLDLFactorization<Field>::MultiplyWithL
( Orientation orientation, Matrix<Field>& B ) const
{
    EL_DEBUG_CSE
    if( !factored_ )
        LogicError("Must call Factor() before MultiplyWithL()");
    ldl::MatrixNode<Field> BNodal( inverseMap_, *info_, B );
    MultiplyWithL( orientation, BNodal );
    BNodal.Push( inverseMap_, *info_, B );
}

template<typename Field>
void SparseLDLFactorization<Field>::MultiplyWithL
( Orientation orientation, ldl::MatrixNode<Field>& B ) const
{
    EL_DEBUG_CSE
    if( !factored_ )
        LogicError("Must call Factor() before MultiplyWithL()");
    if( orientation == NORMAL )
        ldl::LowerForwardMultiply( *info_, *front_, B );
    else
        ldl::LowerBackwardMultiply( *info_, *front_, B, orientation==ADJOINT );
}

template<typename Field>
void SparseLDLFactorization<Field>::SolveAgainstD
( Orientation orientation, Matrix<Field>& B ) const
{
    EL_DEBUG_CSE
    if( !factored_ )
        LogicError("Must call Factor() before SolveAgainstD()");
    ldl::MatrixNode<Field> BNodal( inverseMap_, *info_, B );
    SolveAgainstD( orientation, BNodal );
    BNodal.Push( inverseMap_, *info_, B );
}

template<typename Field>
void SparseLDLFactorization<Field>::SolveAgainstD
( Orientation orientation, ldl::MatrixNode<Field>& B ) const
{
    EL_DEBUG_CSE
    if( !factored_ )
        LogicError("Must call Factor() before SolveAgainstD()");
    if( orientation == NORMAL )
    {
        ldl::DiagonalSolve( *info_, *front_, B );
    }
    else
    {
        LogicError("Conjugated solve with D not yet supported");
    }
}

template<typename Field>
void SparseLDLFactorization<Field>::MultiplyWithD
( Orientation orientation, Matrix<Field>& B ) const
{
    EL_DEBUG_CSE
    if( !factored_ )
        LogicError("Must call Factor() before MultiplyWithD()");
    ldl::MatrixNode<Field> BNodal( inverseMap_, *info_, B );
    MultiplyWithD( orientation, BNodal );
    BNodal.Push( inverseMap_, *info_, B );
}

template<typename Field>
void SparseLDLFactorization<Field>::MultiplyWithD
( Orientation orientation, ldl::MatrixNode<Field>& B ) const
{
    EL_DEBUG_CSE
    if( !factored_ )
        LogicError("Must call Factor() before MultiplyWithD()");
    if( orientation == NORMAL )
    {
        ldl::DiagonalScale( *info_, *front_, B );
    }
    else
    {
        LogicError("Conjugated multiplication with D not yet supported");
    }
}

template<typename Field>
bool SparseLDLFactorization<Field>::Factored() const
{ return factored_; }

template<typename Field>
Int SparseLDLFactorization<Field>::NumEntries() const
{
    EL_DEBUG_CSE
    if( !initialized_ )
        LogicError("Must initialize before calling 'NumEntries()'");
    return front_->NumEntries();
}

template<typename Field>
Int SparseLDLFactorization<Field>::NumTopLeftEntries() const
{
    EL_DEBUG_CSE
    if( !initialized_ )
        LogicError("Must initialize before calling 'NumTopLeftEntries()'");
    return front_->NumTopLeftEntries();
}

template<typename Field>
Int SparseLDLFactorization<Field>::NumBottomLeftEntries() const
{
    EL_DEBUG_CSE
    if( !initialized_ )
        LogicError("Must initialize before calling 'NumBottomLeftEntries()'");
    return front_->NumBottomLeftEntries();
}

template<typename Field>
double SparseLDLFactorization<Field>::FactorGFlops() const
{
    EL_DEBUG_CSE
    if( !initialized_ )
        LogicError("Must initialize before calling 'FactorGFlops()'");
    return front_->FactorGFlops();
}

template<typename Field>
double SparseLDLFactorization<Field>::SolveGFlops( Int numRHS ) const
{
    EL_DEBUG_CSE
    if( !initialized_ )
        LogicError("Must initialize before calling 'SolveGFlops()'");
    return front_->SolveGFlops( numRHS );
}

template<typename Field>
ldl::Front<Field>& SparseLDLFactorization<Field>::Front()
{
    EL_DEBUG_CSE
    return *front_;
}

template<typename Field>
const ldl::Front<Field>& SparseLDLFactorization<Field>::Front() const
{
    EL_DEBUG_CSE
    return *front_;
}

template<typename Field>
ldl::NodeInfo& SparseLDLFactorization<Field>::NodeInfo()
{
    EL_DEBUG_CSE
    return *info_;
}

template<typename Field>
const ldl::NodeInfo& SparseLDLFactorization<Field>::NodeInfo() const
{
    EL_DEBUG_CSE
    return *info_;
}

template<typename Field>
ldl::Separator& SparseLDLFactorization<Field>::Separator()
{
    EL_DEBUG_CSE
    return *separator_;
}

template<typename Field>
const ldl::Separator& SparseLDLFactorization<Field>::Separator() const
{
    EL_DEBUG_CSE
    return *separator_;
}

template<typename Field>
vector<Int>& SparseLDLFactorization<Field>::Map()
{
    EL_DEBUG_CSE
    return map_;
}

template<typename Field>
const vector<Int>& SparseLDLFactorization<Field>::Map() const
{
    EL_DEBUG_CSE
    return map_;
}

template<typename Field>
vector<Int>& SparseLDLFactorization<Field>::InverseMap()
{
    EL_DEBUG_CSE
    return inverseMap_;
}

template<typename Field>
const vector<Int>& SparseLDLFactorization<Field>::InverseMap() const
{
    EL_DEBUG_CSE
    return inverseMap_;
}

#define PROTO(Field) template class SparseLDLFactorization<Field>;

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
