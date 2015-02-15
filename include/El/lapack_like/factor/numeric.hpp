/*
   Copyright 2009-2011, Jack Poulson.
   All rights reserved.

   Copyright 2011-2012, Jack Poulson, Lexing Ying, and 
   The University of Texas at Austin.
   All rights reserved.

   Copyright 2013, Jack Poulson, Lexing Ying, and Stanford University.
   All rights reserved.

   Copyright 2013-2014, Jack Poulson and The Georgia Institute of Technology.
   All rights reserved.

   Copyright 2014-2015, Jack Poulson and Stanford University.
   All rights reserved.
   
   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_FACTOR_NUMERIC_HPP
#define EL_FACTOR_NUMERIC_HPP

namespace El {

template<typename T>
struct DistMatrixNode;
template<typename T>
struct DistMultiVecNode;

template<typename T>
struct MatrixNode
{
    Matrix<T> matrix;
    Matrix<T> work;

    MatrixNode<T>* parent;
    vector<MatrixNode<T>*> children;
    DistMatrixNode<T>* duplicateMat;
    DistMultiVecNode<T>* duplicateMV;

    MatrixNode( MatrixNode<T>* parentNode=nullptr );
    MatrixNode( DistMatrixNode<T>* dupNode );
    MatrixNode( DistMultiVecNode<T>* dupNode );

    MatrixNode
    ( const vector<Int>& invMap, const SymmNodeInfo& info, const Matrix<T>& X );

    void Pull
    ( const vector<Int>& invMap, const SymmNodeInfo& info, const Matrix<T>& X );
    void Push
    ( const vector<Int>& invMap, const SymmNodeInfo& info, Matrix<T>& X ) const;

    Int Height() const;
};

struct MultiVecCommMeta
{
    Int localOff, localSize;
    vector<int> numChildSendInds;
    vector<vector<Int>> childRecvInds;

    void Empty()
    {
        SwapClear( numChildSendInds );
        SwapClear( childRecvInds );
    }
};

// For handling a set of vectors distributed in a [VC,* ] manner over each node
// of the elimination tree
template<typename T>
struct DistMultiVecNode
{   
    DistMatrix<T,VC,STAR> matrix;
    DistMatrix<T,VC,STAR> work;

    DistMultiVecNode<T>* parent;
    DistMultiVecNode<T>* child;
    MatrixNode<T>* duplicate;

    mutable MultiVecCommMeta commMeta; 

    DistMultiVecNode( DistMultiVecNode<T>* parentNode=nullptr );

    DistMultiVecNode
    ( const DistMap& invMap, const DistSymmNodeInfo& info,
      const DistMultiVec<T>& X );

    DistMultiVecNode( const DistMatrixNode<T>& X );

    const DistMultiVecNode<T>& operator=( const DistMatrixNode<T>& X );

    Int LocalHeight() const;

    void Pull
    ( const DistMap& invMap, const DistSymmNodeInfo& info,
      const DistMultiVec<T>& X );
    void Push
    ( const DistMap& invMap, const DistSymmNodeInfo& info,
            DistMultiVec<T>& X ) const;

    void ComputeCommMeta( const DistSymmNodeInfo& info ) const;
};

struct MatrixCommMeta
{
    vector<int> numChildSendInds;
    vector<vector<Int>> childRecvInds;

    void Empty()
    {
        SwapClear( numChildSendInds );
        SwapClear( childRecvInds );
    }
};

// For handling a matrix distributed in a [MC,MR] manner over each node
// of the elimination tree
template<typename T>
struct DistMatrixNode
{
    DistMatrix<T> matrix;
    DistMatrix<T> work;

    DistMatrixNode<T>* parent;
    DistMatrixNode<T>* child;
    MatrixNode<T>* duplicate;

    mutable MatrixCommMeta commMeta; 

    DistMatrixNode( DistMatrixNode<T>* parentNode=nullptr );

    DistMatrixNode
    ( const DistMap& invMap, const DistSymmNodeInfo& info,
      const DistMultiVec<T>& X );

    DistMatrixNode( const DistMultiVecNode<T>& X );

    const DistMatrixNode<T>& operator=( const DistMultiVecNode<T>& X );

    void Pull
    ( const DistMap& invMap, const DistSymmNodeInfo& info,
      const DistMultiVec<T>& X );
    void Push
    ( const DistMap& invMap, const DistSymmNodeInfo& info,
            DistMultiVec<T>& X ) const;

    void ComputeCommMeta( const DistSymmNodeInfo& info ) const;
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

inline bool FrontIs1D( SymmFrontType type )
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

inline SymmFrontType ConvertTo2D( SymmFrontType type )
{
    DEBUG_ONLY(CallStackEntry cse("ConvertTo2D"))
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

inline SymmFrontType ConvertTo1D( SymmFrontType type )
{
    DEBUG_ONLY(CallStackEntry cse("ConvertTo1D"))
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
        return ConvertTo2D(type);
    else if( PivotedFactorization(type) )
        return LDL_INTRAPIV_2D;
    else
        return LDL_2D;
}

// Only keep track of the left and bottom-right piece of the fronts
// (with the bottom-right piece stored in workspace) since only the left side
// needs to be kept after the factorization is complete.

template<typename F>
struct DistSymmFront;

template<typename F>
struct SymmFront
{
    bool isHermitian;
    SymmFrontType type;

    Matrix<F> L;

    Matrix<F> diag;
    Matrix<F> subdiag;
    Matrix<Int> piv;

    Matrix<F> work;

    SymmFront<F>* parent;
    vector<SymmFront<F>*> children;
    DistSymmFront<F>* duplicate;

    SymmFront( SymmFront<F>* parentNode=nullptr );
    SymmFront( DistSymmFront<F>* dupNode );
    SymmFront
    ( const SparseMatrix<F>& A, const vector<Int>& reordering,
      const SymmNodeInfo& rootInfo, bool conjugate=true );

    void Pull
    ( const SparseMatrix<F>& A, const vector<Int>& reordering, 
      const SymmNodeInfo& rootInfo, bool conjugate=true );
    void Push
    ( SparseMatrix<F>& A, const vector<Int>& reordering, 
      const SymmNodeInfo& rootInfo ) const;

    void Unpack( SparseMatrix<F>& A, const SymmNodeInfo& rootInfo ) const;

    // NOTE: These are temporarily not functioning
    double NumEntries() const;
    double NumTopLeftEntries() const;
    double NumBottomLeftEntries() const;
    double FactorGFlops( bool selInv=false ) const;
    double SolveGFlops( Int numRHS=1 ) const;
};

struct FactorCommMeta
{
    vector<int> numChildSendInds;
    // This information does not necessarily have to be kept and can be
    // computed from the above information (albeit somewhat expensively).
    mutable vector<vector<Int>> childRecvInds;

    void EmptyChildRecvIndices() const
    { SwapClear(childRecvInds); }

    void Empty()
    {
        SwapClear( numChildSendInds );
        EmptyChildRecvIndices();
    }
};
void ComputeFactRecvInds( const DistSymmNodeInfo& info );

template<typename F>
struct DistSymmFront
{
    bool isHermitian;
    SymmFrontType type;

    // Split each front into a left and right piece such that the right piece
    // is not needed after the factorization (and can be freed).

    // When this node is a duplicate of a sequential node, L1D or L2D will be 
    // attached to the sequential L matrix of the duplicate

    DistMatrix<F,VC,STAR> L1D;
    DistMatrix<F> L2D;

    DistMatrix<F,VC,STAR> diag;
    DistMatrix<F,VC,STAR> subdiag;
    DistMatrix<Int,VC,STAR> piv;

    DistMatrix<F> work;
    mutable FactorCommMeta commMeta;

    DistSymmFront<F>* parent;
    DistSymmFront<F>* child;
    SymmFront<F>* duplicate;

    DistSymmFront( DistSymmFront<F>* parentNode=nullptr );

    DistSymmFront
    ( const DistSparseMatrix<F>& A,
      const DistMap& reordering,
      const DistSeparator& rootSep,
      const DistSymmNodeInfo& info,
      bool conjugate=false );

    void Pull
    ( const DistSparseMatrix<F>& A,
      const DistMap& reordering,
      const DistSeparator& rootSep,
      const DistSymmNodeInfo& info,
      bool conjugate=false );
    void Push
    ( DistSparseMatrix<F>& A, const DistMap& reordering, 
      const DistSeparator& rootSep, const DistSymmNodeInfo& rootInfo ) const;

    void Unpack
    ( DistSparseMatrix<F>& A, 
      const DistSeparator& rootSep, const DistSymmNodeInfo& rootInfo ) const;

    // NOTE: These are temporarily not functioning
    double NumLocalEntries() const;
    double NumTopLeftLocalEntries() const;
    double NumBottomLeftLocalEntries() const;
    double LocalFactorGFlops( bool selInv=false ) const;
    double LocalSolveGFlops( Int numRHS=1 ) const;

    void ComputeRecvInds( const DistSymmNodeInfo& info ) const;
    void ComputeCommMeta
    ( const DistSymmNodeInfo& info, bool computeRecvInds ) const;
};

template<typename F>
void ChangeFrontType
( DistSymmFront<F>& front, SymmFrontType type, bool recurse=true );

// TODO: Move into BLAS1?
template<typename F>
void DiagonalScale
( const SymmNodeInfo& info, const SymmFront<F>& front,
  MatrixNode<F>& X );
template<typename F>
void DiagonalScale
( const DistSymmNodeInfo& info, const DistSymmFront<F>& front,
  DistMultiVecNode<F>& X );
template<typename F>
void DiagonalScale
( const DistSymmNodeInfo& info, const DistSymmFront<F>& L,
  DistMatrixNode<F>& X );

// TODO: Move into BLAS1?
template<typename F>
void DiagonalSolve
( const SymmNodeInfo& info, const SymmFront<F>& front,
  MatrixNode<F>& X );
template<typename F>
void DiagonalSolve
( const DistSymmNodeInfo& info, const DistSymmFront<F>& front,
  DistMultiVecNode<F>& X );
template<typename F>
void DiagonalSolve
( const DistSymmNodeInfo& info, const DistSymmFront<F>& L,
  DistMatrixNode<F>& X );

// TODO: Rename to Trsm and move to BLAS3?
template<typename F>
void LowerSolve
( Orientation orientation, const SymmNodeInfo& info,
  const SymmFront<F>& L, MatrixNode<F>& X );
template<typename F>
void LowerSolve
( Orientation orientation, const DistSymmNodeInfo& info,
  const DistSymmFront<F>& L, DistMultiVecNode<F>& X );
template<typename F>
void LowerSolve
( Orientation orientation, const DistSymmNodeInfo& info,
  const DistSymmFront<F>& L, DistMatrixNode<F>& X );

// TODO: Rename to Trmm and move to BLAS3?
template<typename F>
void LowerMultiply
( Orientation orientation, const SymmNodeInfo& info,
  const SymmFront<F>& L, MatrixNode<F>& X );
template<typename F>
void LowerMultiply
( Orientation orientation, const DistSymmNodeInfo& info,
  const DistSymmFront<F>& L, DistMultiVecNode<F>& X );
template<typename F>
void LowerMultiply
( Orientation orientation, const DistSymmNodeInfo& info,
  const DistSymmFront<F>& L, DistMatrixNode<F>& X );

} // namespace El

#endif // ifndef EL_FACTOR_NUMERIC_HPP
