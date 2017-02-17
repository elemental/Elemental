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
( DistFront<Field>& front, LDLFrontType type, bool recurse=true );
template<typename Field>
void DiagonalSolve
( const DistNodeInfo& info,
  const DistFront<Field>& front,
        DistMultiVecNode<Field>& B );
template<typename Field>
void DiagonalSolve
( const DistNodeInfo& info,
  const DistFront<Field>& front,
        DistMatrixNode<Field>& B );
template<typename Field>
void DiagonalScale
( const DistNodeInfo& info,
  const DistFront<Field>& front,
        DistMultiVecNode<Field>& B );
template<typename Field>
void DiagonalScale
( const DistNodeInfo& info,
  const DistFront<Field>& front,
        DistMatrixNode<Field>& B );

} // namespace ldl

template<typename Field>
DistSparseLDLFactorization<Field>::DistSparseLDLFactorization()
{ }

template<typename Field>
void DistSparseLDLFactorization<Field>::Initialize
( const DistSparseMatrix<Field>& A,
        bool hermitian,
  const BisectCtrl& bisectCtrl )
{
    EL_DEBUG_CSE
    info_.reset( new ldl::DistNodeInfo(A.Grid()) );
    separator_.reset( new ldl::DistSeparator );
    ldl::NestedDissection
    ( A.LockedDistGraph(), map_, *separator_, *info_, bisectCtrl );
    InvertMap( map_, inverseMap_ );
    front_.reset
    ( new ldl::DistFront<Field>(A,map_,*separator_,*info_,hermitian) );

    initialized_ = true;
    factored_ = false;
}

template<typename Field>
void DistSparseLDLFactorization<Field>::Initialize2DGridGraph
( Int gridDim0,
  Int gridDim1,
  const DistSparseMatrix<Field>& A,
        bool hermitian,
  const BisectCtrl& bisectCtrl )
{
    EL_DEBUG_CSE
    info_.reset( new ldl::DistNodeInfo(A.Grid()) );
    separator_.reset( new ldl::DistSeparator );
    ldl::NaturalNestedDissection
    ( gridDim0, gridDim1, 1, A.LockedDistGraph(),
      map_, *separator_, *info_, bisectCtrl.cutoff );
    InvertMap( map_, inverseMap_ );
    front_.reset
    ( new ldl::DistFront<Field>(A,map_,*separator_,*info_,hermitian) );

    initialized_ = true;
    factored_ = false;
}

template<typename Field>
void DistSparseLDLFactorization<Field>::Initialize3DGridGraph
( Int gridDim0,
  Int gridDim1,
  Int gridDim2,
  const DistSparseMatrix<Field>& A,
        bool hermitian,
  const BisectCtrl& bisectCtrl )
{
    EL_DEBUG_CSE
    info_.reset( new ldl::DistNodeInfo(A.Grid()) );
    separator_.reset( new ldl::DistSeparator );
    ldl::NaturalNestedDissection
    ( gridDim0, gridDim1, gridDim2, A.LockedDistGraph(),
      map_, *separator_, *info_, bisectCtrl.cutoff );
    InvertMap( map_, inverseMap_ );
    front_.reset
    ( new ldl::DistFront<Field>(A,map_,*separator_,*info_,hermitian) );

    initialized_ = true;
    factored_ = false;
}

template<typename Field>
void DistSparseLDLFactorization<Field>::Factor( LDLFrontType frontType )
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
void DistSparseLDLFactorization<Field>::ChangeFrontType
( LDLFrontType frontType )
{
    EL_DEBUG_CSE
    if( !initialized_ )
        LogicError("Must initialize before calling 'ChangeFrontType()'");
    ldl::ChangeFrontType( *front_, frontType );
}

template<typename Field>
void DistSparseLDLFactorization<Field>::ChangeNonzeroValues
( const DistSparseMatrix<Field>& ANew )
{
    EL_DEBUG_CSE
    if( !initialized_ )
        LogicError("Must initialize before calling 'ChangeNonzeroValues()'");
    if( !formedPullMetadata_ )
    {
        ANew.MappedSources( map_, mappedSources_ );
        ANew.MappedTargets( map_, mappedTargets_, columnOffsets_ );
        formedPullMetadata_ = true;
    }
    front_->Pull
    ( ANew, map_, *separator_, *info_,
      mappedSources_, mappedTargets_, columnOffsets_ );
    factored_ = false;
}

template<typename Field>
void DistSparseLDLFactorization<Field>::Solve( DistMultiVec<Field>& B ) const
{
    EL_DEBUG_CSE
    if( !factored_ )
        LogicError("Must call Factor() before Solve()");
    if( FrontIs1D(front_->type) )
    {
        ldl::DistMultiVecNode<Field> BNodal( inverseMap_, *info_, B );
        Solve( BNodal );
        BNodal.Push( inverseMap_, *info_, B );
    }
    else
    {
        ldl::DistMatrixNode<Field> BNodal( inverseMap_, *info_, B );
        Solve( BNodal );
        BNodal.Push( inverseMap_, *info_, B );
    }
}

template<typename Field>
void DistSparseLDLFactorization<Field>::Solve
( ldl::DistMultiVecNode<Field>& B ) const
{
    EL_DEBUG_CSE
    if( !factored_ )
        LogicError("Must call Factor() before Solve()");

    // TODO(poulson): Only perform the switch if there are a sufficient
    // number of right-hand sides?
    /*
    if( !FrontIs1D(front_->type) )
    {
        // TODO(poulson): Add warning?
        ldl::DistMatrixNode<Field> BMat( B );
        Solve( BMat );
        B = BMat;
        return;
    }
    */

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
void DistSparseLDLFactorization<Field>::Solve
( ldl::DistMatrixNode<Field>& B ) const
{
    EL_DEBUG_CSE
    if( !factored_ )
        LogicError("Must call Factor() before Solve()");

    if( FrontIs1D(front_->type) )
    {
        // TODO: Add warning?
        ldl::DistMultiVecNode<Field> BMV( B );
        Solve( BMV );
        B = BMV;
        return;
    }

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
void DistSparseLDLFactorization<Field>::SolveWithIterativeRefinement
( const DistSparseMatrix<Field>& A,
        DistMultiVec<Field>& B,
  const Base<Field>& minReductionFactor,
        Int maxRefineIts ) const
{
    EL_DEBUG_CSE
    if( !factored_ )
        LogicError("Must call Factor() before SolveWithIterativeRefinement()");
    // TODO(poulson): Generalize this implementation
    if( B.Width() > 1 )
        LogicError("Iterative Refinement currently only supported for one RHS");
    const Grid& grid = B.Grid();

    DistMultiVec<Field> BOrig(grid);
    BOrig = B;

    // Compute the initial guess
    // =========================
    DistMultiVec<Field> X(B);
    Solve( X );

    Int refineIt = 0;
    if( maxRefineIts > 0 )
    {
        DistMultiVec<Field> dX(grid), XCand(grid);
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
void DistSparseLDLFactorization<Field>::SolveAgainstL
( Orientation orientation, DistMultiVec<Field>& B ) const
{
    EL_DEBUG_CSE
    if( !factored_ )
        LogicError("Must call Factor() before SolveAgainstL()");
    if( FrontIs1D(front_->type) )
    {
        ldl::DistMultiVecNode<Field> BNodal( inverseMap_, *info_, B );
        SolveAgainstL( orientation, BNodal );
        BNodal.Push( inverseMap_, *info_, B );
    }
    else
    {
        ldl::DistMatrixNode<Field> BNodal( inverseMap_, *info_, B );
        SolveAgainstL( orientation, BNodal );
        BNodal.Push( inverseMap_, *info_, B );
    }
}

template<typename Field>
void DistSparseLDLFactorization<Field>::SolveAgainstL
( Orientation orientation, ldl::DistMultiVecNode<Field>& B ) const
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
void DistSparseLDLFactorization<Field>::SolveAgainstL
( Orientation orientation, ldl::DistMatrixNode<Field>& B ) const
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
void DistSparseLDLFactorization<Field>::MultiplyWithL
( Orientation orientation, DistMultiVec<Field>& B ) const
{
    EL_DEBUG_CSE
    if( !factored_ )
        LogicError("Must call Factor() before MultiplyWithL()");
    if( FrontIs1D(front_->type) )
    {
        ldl::DistMultiVecNode<Field> BNodal( inverseMap_, *info_, B );
        MultiplyWithL( orientation, BNodal );
        BNodal.Push( inverseMap_, *info_, B );
    }
    else
    {
        ldl::DistMatrixNode<Field> BNodal( inverseMap_, *info_, B );
        MultiplyWithL( orientation, BNodal );
        BNodal.Push( inverseMap_, *info_, B );
    }
}

template<typename Field>
void DistSparseLDLFactorization<Field>::MultiplyWithL
( Orientation orientation, ldl::DistMultiVecNode<Field>& B ) const
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
void DistSparseLDLFactorization<Field>::MultiplyWithL
( Orientation orientation, ldl::DistMatrixNode<Field>& B ) const
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
void DistSparseLDLFactorization<Field>::SolveAgainstD
( Orientation orientation, DistMultiVec<Field>& B ) const
{
    EL_DEBUG_CSE
    if( !factored_ )
        LogicError("Must call Factor() before SolveAgainstD()");
    if( FrontIs1D(front_->type) )
    {
        ldl::DistMultiVecNode<Field> BNodal( inverseMap_, *info_, B );
        SolveAgainstD( orientation, BNodal );
        BNodal.Push( inverseMap_, *info_, B );
    }
    else
    {
        ldl::DistMatrixNode<Field> BNodal( inverseMap_, *info_, B );
        SolveAgainstD( orientation, BNodal );
        BNodal.Push( inverseMap_, *info_, B );
    }
}

template<typename Field>
void DistSparseLDLFactorization<Field>::SolveAgainstD
( Orientation orientation, ldl::DistMultiVecNode<Field>& B ) const
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
        LogicError("Conjugated SolveAgainstD() not yet supported");
    }
}

template<typename Field>
void DistSparseLDLFactorization<Field>::SolveAgainstD
( Orientation orientation, ldl::DistMatrixNode<Field>& B ) const
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
        LogicError("Conjugated SolveAgainstD() not yet supported");
    }
}

template<typename Field>
void DistSparseLDLFactorization<Field>::MultiplyWithD
( Orientation orientation, DistMultiVec<Field>& B ) const
{
    EL_DEBUG_CSE
    if( !factored_ )
        LogicError("Must call Factor() before MultiplyWithD()");
    if( FrontIs1D(front_->type) )
    {
        ldl::DistMultiVecNode<Field> BNodal( inverseMap_, *info_, B );
        MultiplyWithD( orientation, BNodal );
        BNodal.Push( inverseMap_, *info_, B );
    }
    else
    {
        ldl::DistMatrixNode<Field> BNodal( inverseMap_, *info_, B );
        MultiplyWithD( orientation, BNodal );
        BNodal.Push( inverseMap_, *info_, B );
    }
}

template<typename Field>
void DistSparseLDLFactorization<Field>::MultiplyWithD
( Orientation orientation, ldl::DistMultiVecNode<Field>& B ) const
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
        LogicError("Conjugated MultiplyWithD() not yet supported");
    }
}

template<typename Field>
void DistSparseLDLFactorization<Field>::MultiplyWithD
( Orientation orientation, ldl::DistMatrixNode<Field>& B ) const
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
        LogicError("Conjugated MultiplyWithD() not yet supported");
    }
}

template<typename Field>
bool DistSparseLDLFactorization<Field>::Factored() const
{ return factored_; }

template<typename Field>
Int DistSparseLDLFactorization<Field>::NumLocalEntries() const
{
    EL_DEBUG_CSE
    if( !initialized_ )
        LogicError("Must initialize before calling 'NumLocalEntries()'");
    return front_->NumLocalEntries();
}

template<typename Field>
Int DistSparseLDLFactorization<Field>::NumTopLeftLocalEntries() const
{
    EL_DEBUG_CSE
    if( !initialized_ )
        LogicError("Must initialize before calling 'NumTopLeftLocalEntries()'");
    return front_->NumTopLeftLocalEntries();
}

template<typename Field>
Int DistSparseLDLFactorization<Field>::NumBottomLeftLocalEntries() const
{
    EL_DEBUG_CSE
    if( !initialized_ )
        LogicError
        ("Must initialize before calling 'NumBottomLeftLocalEntries()'");
    return front_->NumBottomLeftLocalEntries();
}

template<typename Field>
double DistSparseLDLFactorization<Field>::LocalFactorGFlops
( bool selectiveInversion ) const
{
    EL_DEBUG_CSE
    if( !initialized_ )
        LogicError("Must initialize before calling 'LocalFactorGFlops()'");
    return front_->LocalFactorGFlops( selectiveInversion );
}

template<typename Field>
double DistSparseLDLFactorization<Field>::LocalSolveGFlops( Int numRHS ) const
{
    EL_DEBUG_CSE
    if( !initialized_ )
        LogicError("Must initialize before calling 'LocalSolveGFlops()'");
    return front_->LocalSolveGFlops( numRHS );
}

template<typename Field>
ldl::DistFront<Field>& DistSparseLDLFactorization<Field>::Front()
{
    EL_DEBUG_CSE
    return *front_;
}

template<typename Field>
const ldl::DistFront<Field>& DistSparseLDLFactorization<Field>::Front() const
{
    EL_DEBUG_CSE
    return *front_;
}

template<typename Field>
ldl::DistNodeInfo& DistSparseLDLFactorization<Field>::NodeInfo()
{
    EL_DEBUG_CSE
    return *info_;
}

template<typename Field>
const ldl::DistNodeInfo& DistSparseLDLFactorization<Field>::NodeInfo() const
{
    EL_DEBUG_CSE
    return *info_;
}

template<typename Field>
ldl::DistSeparator& DistSparseLDLFactorization<Field>::Separator()
{
    EL_DEBUG_CSE
    return *separator_;
}

template<typename Field>
const ldl::DistSeparator& DistSparseLDLFactorization<Field>::Separator() const
{
    EL_DEBUG_CSE
    return *separator_;
}

template<typename Field>
DistMap& DistSparseLDLFactorization<Field>::Map()
{
    EL_DEBUG_CSE
    return map_;
}

template<typename Field>
const DistMap& DistSparseLDLFactorization<Field>::Map() const
{
    EL_DEBUG_CSE
    return map_;
}

template<typename Field>
DistMap& DistSparseLDLFactorization<Field>::InverseMap()
{
    EL_DEBUG_CSE
    return inverseMap_;
}

template<typename Field>
const DistMap& DistSparseLDLFactorization<Field>::InverseMap() const
{
    EL_DEBUG_CSE
    return inverseMap_;
}

template<typename Field>
ldl::DistMultiVecNodeMeta&
DistSparseLDLFactorization<Field>::DistMultiVecNodeMeta() const
{
    EL_DEBUG_CSE
    return dmvMeta_;
}

#define PROTO(Field) template class DistSparseLDLFactorization<Field>;

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
