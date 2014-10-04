/*
   Copyright (c) 2009-2014, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, Stanford University, and the
   Georgia Insitute of Technology.
   All rights reserved.
 
   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_SPARSEDIRECT_NUMERIC_HPP
#define EL_SPARSEDIRECT_NUMERIC_HPP

namespace El {

// Forward declaration
template<typename T> class DistNodalMultiVec;

// For handling a matrix distributed in a [MC,MR] manner over each node
// of the elimination tree
template<typename T>
class DistNodalMatrix
{
public:
    std::vector<Matrix<T>> localNodes;
    std::vector<DistMatrix<T>> distNodes;

    DistNodalMatrix();
    DistNodalMatrix
    ( const DistMap& inverseMap, const DistSymmInfo& info,
      const DistMultiVec<T>& X );
    DistNodalMatrix( const DistNodalMultiVec<T>& X );

    const DistNodalMatrix<T>& operator=( const DistNodalMultiVec<T>& X );

    void Pull
    ( const DistMap& inverseMap, const DistSymmInfo& info,
      const DistMultiVec<T>& X );
    void Push
    ( const DistMap& inverseMap, const DistSymmInfo& info,
            DistMultiVec<T>& X ) const;

    int Height() const;
    int Width() const;

    mutable std::vector<MatrixCommMeta> commMetas;
    void ComputeCommMetas( const DistSymmInfo& info ) const;
private:
    int height_, width_;
};

// For handling a set of vectors distributed in a [VC,* ] manner over each node
// of the elimination tree
template<typename T>
class DistNodalMultiVec
{
public:
    std::vector<Matrix<T>> localNodes;
    std::vector<DistMatrix<T,VC,STAR>> distNodes;

    DistNodalMultiVec();
    DistNodalMultiVec
    ( const DistMap& inverseMap, const DistSymmInfo& info,
      const DistMultiVec<T>& X );
    DistNodalMultiVec( const DistNodalMatrix<T>& X );

    const DistNodalMultiVec<T>& operator=( const DistNodalMatrix<T>& X );

    void Pull
    ( const DistMap& inverseMap, const DistSymmInfo& info,
      const DistMultiVec<T>& X );
    void Push
    ( const DistMap& inverseMap, const DistSymmInfo& info,
            DistMultiVec<T>& X ) const;

    int Height() const;
    int Width() const;

    int LocalHeight() const;

    void UpdateHeight();
    void UpdateWidth();

private:
    int height_, width_;
};

enum SymmFrontType
{
  SYMM_1D,                SYMM_2D,
  LDL_1D,                 LDL_2D,
  LDL_SELINV_1D,          LDL_SELINV_2D,
  LDL_INTRAPIV_1D,        LDL_INTRAPIV_2D,
  LDL_INTRAPIV_SELINV_1D, LDL_INTRAPIV_SELINV_2D,
  BLOCK_LDL_1D,           BLOCK_LDL_2D,
  BLOCK_LDL_INTRAPIV_1D,  BLOCK_LDL_INTRAPIV_2D
};

inline bool Unfactored( SymmFrontType type )
{ return type == SYMM_1D || type == SYMM_2D; }

inline bool FrontsAre1d( SymmFrontType type )
{
    return type == SYMM_1D                ||
           type == LDL_1D                 ||
           type == LDL_SELINV_1D          ||
           type == LDL_INTRAPIV_1D        ||
           type == LDL_INTRAPIV_SELINV_1D ||
           type == BLOCK_LDL_1D           ||
           type == BLOCK_LDL_INTRAPIV_1D;
}

inline bool BlockFactorization( SymmFrontType type )
{
    return type == BLOCK_LDL_1D ||
           type == BLOCK_LDL_2D ||
           type == BLOCK_LDL_INTRAPIV_1D ||
           type == BLOCK_LDL_INTRAPIV_2D;
}

inline bool SelInvFactorization( SymmFrontType type )
{
    return type == LDL_SELINV_1D ||
           type == LDL_SELINV_2D ||
           type == LDL_INTRAPIV_SELINV_1D ||
           type == LDL_INTRAPIV_SELINV_2D;
}

inline bool PivotedFactorization( SymmFrontType type )
{
    return type == LDL_INTRAPIV_1D ||
           type == LDL_INTRAPIV_2D ||
           type == LDL_INTRAPIV_SELINV_1D ||
           type == LDL_INTRAPIV_SELINV_2D ||
           type == BLOCK_LDL_INTRAPIV_1D  ||
           type == BLOCK_LDL_INTRAPIV_2D;
}

inline SymmFrontType ConvertTo2d( SymmFrontType type )
{
    DEBUG_ONLY(CallStackEntry cse("ConvertTo2d"))
    SymmFrontType newType;
    switch( type )
    {
    case SYMM_1D:
    case SYMM_2D:                newType = SYMM_2D;                break;
    case LDL_1D:
    case LDL_2D:                 newType = LDL_2D;                 break;
    case LDL_SELINV_1D:
    case LDL_SELINV_2D:          newType = LDL_SELINV_2D;          break;
    case LDL_INTRAPIV_1D:
    case LDL_INTRAPIV_2D:        newType = LDL_INTRAPIV_2D;        break;
    case LDL_INTRAPIV_SELINV_1D:
    case LDL_INTRAPIV_SELINV_2D: newType = LDL_INTRAPIV_SELINV_2D; break;
    case BLOCK_LDL_1D:
    case BLOCK_LDL_2D:           newType = BLOCK_LDL_2D;           break;
    case BLOCK_LDL_INTRAPIV_1D:
    case BLOCK_LDL_INTRAPIV_2D:  newType = BLOCK_LDL_INTRAPIV_2D;  break;
    default: LogicError("Invalid front type");
    }
    return newType;
}

inline SymmFrontType ConvertTo1d( SymmFrontType type )
{
    DEBUG_ONLY(CallStackEntry cse("ConvertTo1d"))
    SymmFrontType newType;
    switch( type )
    {
    case SYMM_1D:
    case SYMM_2D:                newType = SYMM_1D;                break;
    case LDL_1D:
    case LDL_2D:                 newType = LDL_1D;                 break;
    case LDL_SELINV_1D:
    case LDL_SELINV_2D:          newType = LDL_SELINV_1D;          break;
    case LDL_INTRAPIV_1D:
    case LDL_INTRAPIV_2D:        newType = LDL_INTRAPIV_1D;        break;
    case LDL_INTRAPIV_SELINV_1D:
    case LDL_INTRAPIV_SELINV_2D: newType = LDL_INTRAPIV_SELINV_1D; break;
    case BLOCK_LDL_1D:
    case BLOCK_LDL_2D:           newType = BLOCK_LDL_1D;           break;
    case BLOCK_LDL_INTRAPIV_1D:
    case BLOCK_LDL_INTRAPIV_2D:  newType = BLOCK_LDL_INTRAPIV_1D;  break;
    default: LogicError("Invalid front type");
    }
    return newType;
}

inline SymmFrontType AppendSelInv( SymmFrontType type )
{
    DEBUG_ONLY(CallStackEntry cse("AppendSelInv"))
    SymmFrontType newType;
    switch( type )
    {
    case LDL_1D:          newType = LDL_SELINV_1D; break;
    case LDL_2D:          newType = LDL_SELINV_2D; break;
    case LDL_INTRAPIV_1D: newType = LDL_INTRAPIV_SELINV_1D; break;
    case LDL_INTRAPIV_2D: newType = LDL_INTRAPIV_SELINV_2D; break;
    default: LogicError("Sel-inv does not make sense for this type");
    }
    return newType;
}

inline SymmFrontType RemoveSelInv( SymmFrontType type )
{
    DEBUG_ONLY(CallStackEntry cse("RemoveSelInv"))
    SymmFrontType newType;
    switch( type )
    {
    case LDL_SELINV_1D: newType = LDL_1D; break;
    case LDL_SELINV_2D: newType = LDL_2D; break;
    case LDL_INTRAPIV_SELINV_1D: newType = LDL_INTRAPIV_1D; break;
    case LDL_INTRAPIV_SELINV_2D: newType = LDL_INTRAPIV_2D; break;
    default: LogicError("This type did not involve selective inversion");
    }
    return newType;
}

inline SymmFrontType InitialFactorType( SymmFrontType type )
{
    if( Unfactored(type) )
        LogicError("Front type does not require factorization");
    if( BlockFactorization(type) )
        return ConvertTo2d(type);
    else if( PivotedFactorization(type) )
        return LDL_INTRAPIV_2D;
    else
        return LDL_2D;
}

// Only keep track of the left and bottom-right piece of the fronts
// (with the bottom-right piece stored in workspace) since only the left side
// needs to be kept after the factorization is complete.

template<typename T>
struct SymmFront
{
    Matrix<T> frontL;

    Matrix<T> diag;
    Matrix<T> subdiag;
    Matrix<Int> piv;

    mutable Matrix<T> work;
};

template<typename T>
struct DistSymmFront
{
    // The 'frontType' member variable of the parent 'DistSymmFrontTree' 
    // determines which of the following fronts is active.
    //
    // Split each front into a left and right piece such that the right piece
    // is not needed after the factorization (and can be freed).

    DistMatrix<T,VC,STAR> front1dL;
    DistMatrix<T> front2dL;

    DistMatrix<T,VC,STAR> diag1d;
    DistMatrix<T,VC,STAR> subdiag1d;
    DistMatrix<Int,VC,STAR> piv;

    mutable DistMatrix<T,VC,STAR> work1d;
    mutable DistMatrix<T> work2d;
};

template<typename T>
struct DistSymmFrontTree
{
    bool isHermitian;
    SymmFrontType frontType;
    std::vector<SymmFront<T>> localFronts;
    std::vector<DistSymmFront<T>> distFronts;

    DistSymmFrontTree();

    DistSymmFrontTree
    ( const DistSparseMatrix<T>& A,
      const DistMap& reordering,
      const DistSeparatorTree& sepTree,
      const DistSymmInfo& info,
      bool conjugate=false );

    void Initialize
    ( const DistSparseMatrix<T>& A,
      const DistMap& reordering,
      const DistSeparatorTree& sepTree,
      const DistSymmInfo& info,
      bool conjugate=false );

    void TopLeftMemoryInfo
    ( double& numLocalEntries, double& minLocalEntries, double& maxLocalEntries,
      double& numGlobalEntries ) const;

    void BottomLeftMemoryInfo
    ( double& numLocalEntries, double& minLocalEntries, double& maxLocalEntries,
      double& numGlobalEntries ) const;

    void MemoryInfo
    ( double& numLocalEntries, double& minLocalEntries, double& maxLocalEntries,
      double& numGlobalEntries ) const;

    void FactorizationWork
    ( double& numLocalFlops, double& minLocalFlops, double& maxLocalFlops,
      double& numGlobalFlops, bool selInv=false ) const;

    void SolveWork
    ( double& numLocalFlops, double& minLocalFlops, double& maxLocalFlops,
      double& numGlobalFlops, int numRhs=1 ) const;
};

template<typename F>
void ChangeFrontType( DistSymmFrontTree<F>& L, SymmFrontType frontType );

template<typename F>
void DiagonalSolve
( const DistSymmInfo& info, const DistSymmFrontTree<F>& L,
  DistNodalMultiVec<F>& X );
template<typename F>
void DiagonalSolve
( const DistSymmInfo& info, const DistSymmFrontTree<F>& L,
  DistNodalMatrix<F>& X );

// All fronts of L are required to be initialized to the expansions of the 
// original sparse matrix before calling LDL.
template<typename F>
void LDL
( DistSymmInfo& info, DistSymmFrontTree<F>& L,
  SymmFrontType newFrontType=LDL_2D );

// TODO: Decide if this interface needs to be exposed
template<typename T>
void InitializeDistLeaf( const DistSymmInfo& info, DistSymmFrontTree<T>& L );

template<typename T>
void LowerMultiply
( Orientation orientation, int diagOffset,
  const DistSymmInfo& info, const DistSymmFrontTree<T>& L,
  DistNodalMultiVec<T>& X );

template<typename F>
void LowerSolve
( Orientation orientation, const DistSymmInfo& info,
  const DistSymmFrontTree<F>& L, DistNodalMultiVec<F>& X );
template<typename F>
void LowerSolve
( Orientation orientation, const DistSymmInfo& info,
  const DistSymmFrontTree<F>& L, DistNodalMatrix<F>& X );

template<typename F>
void Solve
( const DistSymmInfo& info,
  const DistSymmFrontTree<F>& L, DistNodalMultiVec<F>& X );
template<typename F>
void Solve
( const DistSymmInfo& info,
  const DistSymmFrontTree<F>& L, DistNodalMatrix<F>& X );

template<typename F>
void SymmetricSolve
( const DistSparseMatrix<F>& A, DistMultiVec<F>& X,
  bool conjugate=false,
  bool sequential=true, int numDistSeps=1, int numSeqSeps=1, int cutoff=128 );
template<typename F>
void HermitianSolve
( const DistSparseMatrix<F>& A, DistMultiVec<F>& X,
  bool sequential=true, int numDistSeps=1, int numSeqSeps=1, int cutoff=128 );

} // namespace El

#endif // ifndef EL_SPARSEDIRECT_NUMERIC_HPP
