/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   Copyright (c) 2012 Jack Poulson, Lexing Ying, and
   The University of Texas at Austin.
   All rights reserved.

   Copyright (c) 2013 Jack Poulson, Lexing Ying, and Stanford University.
   All rights reserved.

   Copyright (c) 2014 Jack Poulson and The Georgia Institute of Technology.
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_FACTOR_LDL_SPARSE_NUMERIC_HPP
#define EL_FACTOR_LDL_SPARSE_NUMERIC_HPP

#define EL_SUITESPARSE_NO_SCALAR_FUNCS
#include <ElSuiteSparse/ldl.hpp>

#include <El/lapack_like/factor/ldl/sparse/symbolic.hpp>

namespace El {

enum LDLFrontType
{
  SYMM_1D,                SYMM_2D,
  LDL_1D,                 LDL_2D,
  LDL_SELINV_1D,          LDL_SELINV_2D,
  LDL_INTRAPIV_1D,        LDL_INTRAPIV_2D,
  LDL_INTRAPIV_SELINV_1D, LDL_INTRAPIV_SELINV_2D,
  BLOCK_LDL_1D,           BLOCK_LDL_2D,
  BLOCK_LDL_INTRAPIV_1D,  BLOCK_LDL_INTRAPIV_2D
};

bool Unfactored( LDLFrontType type );
bool FrontIs1D( LDLFrontType type );
bool BlockFactorization( LDLFrontType type );
bool SelInvFactorization( LDLFrontType type );
bool PivotedFactorization( LDLFrontType type );
LDLFrontType ConvertTo2D( LDLFrontType type );
LDLFrontType ConvertTo1D( LDLFrontType type );
LDLFrontType AppendSelInv( LDLFrontType type );
LDLFrontType RemoveSelInv( LDLFrontType type );
LDLFrontType InitialFactorType( LDLFrontType type );

namespace ldl {

template<typename T>
struct DistMatrixNode;
template<typename T>
struct DistMultiVecNode;

template<typename T>
struct MatrixNode
{
    Matrix<T> matrix;
    Matrix<T> work;

    // An observing pointer to the parent (should it exist).
    MatrixNode<T>* parent=nullptr;

    // An observing pointer to the equivalent (trivially) 2D distributed node
    // (should it exist).
    DistMatrixNode<T>* duplicateMat=nullptr;

    // An observing pointer to the equivalent (trivially) 1D distributed node
    // (should it exist).
    DistMultiVecNode<T>* duplicateMV=nullptr;

    // Unique pointers to the children.
    vector<unique_ptr<MatrixNode<T>>> children;

    MatrixNode( MatrixNode<T>* parentNode=nullptr );

    MatrixNode( DistMatrixNode<T>* dupNode );

    MatrixNode( DistMultiVecNode<T>* dupNode );

    MatrixNode
    ( const vector<Int>& invMap,
      const NodeInfo& info,
      const Matrix<T>& X );

    ~MatrixNode();

    const MatrixNode<T>& operator=( const MatrixNode<T>& X );

    void Pull
    ( const vector<Int>& invMap,
      const NodeInfo& info,
      const Matrix<T>& X );
    void Push
    ( const vector<Int>& invMap,
      const NodeInfo& info,
            Matrix<T>& X ) const;

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

struct DistMultiVecNodeMeta
{
    vector<Int> sendInds;
    vector<Int> recvInds;
    vector<int> mappedOwners;
    vector<int> sendSizes;
    vector<int> sendOffs;
    vector<int> recvSizes;
    vector<int> recvOffs;

    template<typename T>
    void Initialize
    ( const DistMultiVecNode<T>& XNode,
      const DistNodeInfo& info,
      const DistMap& invMap,
      const DistMultiVec<T>& X );
};

// For handling a set of vectors distributed in a [VC,* ] manner over each node
// of the elimination tree
template<typename T>
struct DistMultiVecNode
{
    DistMatrix<T,VC,STAR> matrix;
    DistMatrix<T,VC,STAR> work;

    // An observing pointer to the parent (should it exist).
    DistMultiVecNode<T>* parent=nullptr;

    // A unique pointer to the equivalent sequential node (should it exist).
    unique_ptr<MatrixNode<T>> duplicate;

    // A unique pointer to the child node shared by this process
    // (should it exist).
    unique_ptr<DistMultiVecNode<T>> child;

    mutable MultiVecCommMeta commMeta;

    DistMultiVecNode( DistMultiVecNode<T>* parentNode=nullptr );

    DistMultiVecNode
    ( const DistMap& invMap,
      const DistNodeInfo& info,
      const DistMultiVec<T>& X );

    DistMultiVecNode( const DistMatrixNode<T>& X );

    ~DistMultiVecNode();

    const DistMultiVecNode<T>& operator=( const DistMatrixNode<T>& X );

    void Pull
    ( const DistMap& invMap,
      const DistNodeInfo& info,
      const DistMultiVec<T>& X );
    void Pull
    ( const DistMap& invMap,
      const DistNodeInfo& info,
      const DistMultiVec<T>& X,
            DistMultiVecNodeMeta& meta );
    void Push
    ( const DistMap& invMap,
      const DistNodeInfo& info,
            DistMultiVec<T>& X ) const;
    void Push
    ( const DistMap& invMap,
      const DistNodeInfo& info,
            DistMultiVec<T>& X,
            DistMultiVecNodeMeta& meta ) const;

    Int LocalHeight() const;
    void ComputeCommMeta( const DistNodeInfo& info ) const;
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

    // An observing pointer to the parent node (should it exist).
    DistMatrixNode<T>* parent=nullptr;

    // A unique pointer to the equivalent sequential node (should it exist).
    unique_ptr<MatrixNode<T>> duplicate;

    // A unique pointer to the child node shared by this process
    // (should it exist).
    unique_ptr<DistMatrixNode<T>> child;

    mutable MatrixCommMeta commMeta;

    DistMatrixNode( DistMatrixNode<T>* parentNode=nullptr );

    DistMatrixNode
    ( const DistMap& invMap, const DistNodeInfo& info,
      const DistMultiVec<T>& X );

    DistMatrixNode( const DistMultiVecNode<T>& X );

    ~DistMatrixNode();

    const DistMatrixNode<T>& operator=( const DistMultiVecNode<T>& X );

    void Pull
    ( const DistMap& invMap, const DistNodeInfo& info,
      const DistMultiVec<T>& X );
    void Push
    ( const DistMap& invMap, const DistNodeInfo& info,
            DistMultiVec<T>& X ) const;

    void ComputeCommMeta( const DistNodeInfo& info ) const;
};

// Only keep track of the left and bottom-right piece of the fronts
// (with the bottom-right piece stored in workspace) since only the left side
// needs to be kept after the factorization is complete.

template<typename Field>
struct DistFront;

template<typename Field>
struct Front
{
    bool isHermitian;
    bool sparseLeaf=false;
    LDLFrontType type;

    Matrix<Field> LDense;
    SparseMatrix<Field> LSparse;

    Matrix<Field> diag;
    Matrix<Field> subdiag;
    Permutation p;

    Matrix<Field> workDense;
    SparseMatrix<Field> workSparse;

    // An observing pointer for the parent front (should it exist).
    Front<Field>* parent=nullptr;

    // An observing pointer for the duplicate distributed front
    // (should it exist).
    DistFront<Field>* duplicate=nullptr;

    // Unique pointers to the child fronts (should they exist).
    vector<unique_ptr<Front<Field>>> children;

    Front( Front<Field>* parentNode=nullptr );

    Front( DistFront<Field>* dupNode );

    Front
    ( const SparseMatrix<Field>& A,
      const vector<Int>& reordering,
      const NodeInfo& rootInfo,
      bool hermitian=true );

    ~Front();

    void Pull
    ( const SparseMatrix<Field>& A,
      const vector<Int>& reordering,
      const NodeInfo& rootInfo,
      bool hermitian=true );
    void PullUpdate
    ( const SparseMatrix<Field>& A,
      const vector<Int>& reordering,
      const NodeInfo& rootInfo );

    void Push
    (       SparseMatrix<Field>& A,
      const vector<Int>& reordering,
      const NodeInfo& rootInfo ) const;

    void Unpack( SparseMatrix<Field>& A, const NodeInfo& rootInfo ) const;

    const Front<Field>& operator=( const Front<Field>& front );

    Int Height() const;
    Int NumEntries() const;
    Int NumTopLeftEntries() const;
    Int NumBottomLeftEntries() const;
    double FactorGFlops() const;
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
void ComputeFactRecvInds( const DistNodeInfo& info );

template<typename Field>
struct DistFront
{
    bool isHermitian;
    LDLFrontType type;

    // Split each front into a left and right piece such that the right piece
    // is not needed after the factorization (and can be freed).

    // When this node is a duplicate of a sequential node, L1D or L2D will be
    // attached to the sequential L matrix of the duplicate

    DistMatrix<Field,VC,STAR> L1D;
    DistMatrix<Field> L2D;

    DistMatrix<Field,VC,STAR> diag;
    DistMatrix<Field,VC,STAR> subdiag;
    DistPermutation p;

    DistMatrix<Field> work;
    mutable FactorCommMeta commMeta;

    // An observing pointer to the parent node (should it exist).
    DistFront<Field>* parent=nullptr;

    // A unique pointer to the sequential duplicate front (should it exist).
    unique_ptr<Front<Field>> duplicate;

    // A unique pointer to the distributed child front shared by this process
    // (should it exist).
    unique_ptr<DistFront<Field>> child;

    DistFront( DistFront<Field>* parentNode=nullptr );

    DistFront
    ( const DistSparseMatrix<Field>& A,
      const DistMap& reordering,
      const DistSeparator& rootSep,
      const DistNodeInfo& info,
      bool hermitian=false );

    ~DistFront();

    void Pull
    ( const DistSparseMatrix<Field>& A,
      const DistMap& reordering,
      const DistSeparator& rootSep,
      const DistNodeInfo& info,
      bool hermitian=false );
    // Allow the reuse of mapped{Sources,Targets}, which are expensive to form
    void Pull
    ( const DistSparseMatrix<Field>& A,
      const DistMap& reordering,
      const DistSeparator& rootSep,
      const DistNodeInfo& info,
            vector<Int>& mappedSources,
            vector<Int>& mappedTargets,
            vector<Int>& colOffs,
      bool hermitian=false );

    void PullUpdate
    ( const DistSparseMatrix<Field>& A,
      const DistMap& reordering,
      const DistSeparator& rootSep,
      const DistNodeInfo& info );
    // Allow the reuse of mapped{Sources,Targets}, which are expensive to form
    void PullUpdate
    ( const DistSparseMatrix<Field>& A,
      const DistMap& reordering,
      const DistSeparator& rootSep,
      const DistNodeInfo& info,
            vector<Int>& mappedSources,
            vector<Int>& mappedTargets,
            vector<Int>& colOffs );

    // NOTE: This routine is not yet functioning
    void Push
    ( DistSparseMatrix<Field>& A, const DistMap& reordering,
      const DistSeparator& rootSep, const DistNodeInfo& rootInfo ) const;

    void Unpack
    ( DistSparseMatrix<Field>& A,
      const DistSeparator& rootSep, const DistNodeInfo& rootInfo ) const;

    const DistFront<Field>& operator=( const DistFront<Field>& front );

    Int NumLocalEntries() const;
    Int NumTopLeftLocalEntries() const;
    Int NumBottomLeftLocalEntries() const;
    double LocalFactorGFlops( bool selInv=false ) const;
    double LocalSolveGFlops( Int numRHS=1 ) const;

    void ComputeRecvInds( const DistNodeInfo& info ) const;
    void ComputeCommMeta
    ( const DistNodeInfo& info, bool computeRecvInds ) const;
};

template<typename Field>
void ChangeFrontType
( Front<Field>& front, LDLFrontType type, bool recurse=true );
template<typename Field>
void ChangeFrontType
( DistFront<Field>& front, LDLFrontType type, bool recurse=true );

template<typename Field>
void DiagonalScale
( const NodeInfo& info, const Front<Field>& front, MatrixNode<Field>& X );
template<typename Field>
void DiagonalScale
( const DistNodeInfo& info,
  const DistFront<Field>& L,
        DistMultiVecNode<Field>& X );
template<typename Field>
void DiagonalScale
( const DistNodeInfo& info,
  const DistFront<Field>& L,
        DistMatrixNode<Field>& X );

template<typename Field>
void DiagonalSolve
( const NodeInfo& info,
  const Front<Field>& front,
        MatrixNode<Field>& X );
template<typename Field>
void DiagonalSolve
( const DistNodeInfo& info,
  const DistFront<Field>& front,
        DistMultiVecNode<Field>& X );
template<typename Field>
void DiagonalSolve
( const DistNodeInfo& info,
  const DistFront<Field>& L,
        DistMatrixNode<Field>& X );

template<typename Field>
void LowerSolve
( Orientation orientation, const NodeInfo& info,
  const Front<Field>& L, MatrixNode<Field>& X );
template<typename Field>
void LowerSolve
( Orientation orientation, const DistNodeInfo& info,
  const DistFront<Field>& L, DistMultiVecNode<Field>& X );
template<typename Field>
void LowerSolve
( Orientation orientation, const DistNodeInfo& info,
  const DistFront<Field>& L, DistMatrixNode<Field>& X );

template<typename Field>
void LowerMultiply
( Orientation orientation, const NodeInfo& info,
  const Front<Field>& L, MatrixNode<Field>& X );
template<typename Field>
void LowerMultiply
( Orientation orientation, const DistNodeInfo& info,
  const DistFront<Field>& L, DistMultiVecNode<Field>& X );
template<typename Field>
void LowerMultiply
( Orientation orientation, const DistNodeInfo& info,
  const DistFront<Field>& L, DistMatrixNode<Field>& X );

} // namespace ldl

template<typename Field>
class SparseLDLFactorization
{
public:
    SparseLDLFactorization();

    void Initialize
    ( const SparseMatrix<Field>& A,
            bool hermitian=true,
      const BisectCtrl& bisectCtrl=BisectCtrl() );

    // TODO(poulson): Combine these two interfaces
    void Initialize2DGridGraph
    ( Int gridDim0,
      Int gridDim1,
      const SparseMatrix<Field>& A,
            bool hermitian=true,
      const BisectCtrl& bisectCtrl=BisectCtrl() );
    void Initialize3DGridGraph
    ( Int gridDim0,
      Int gridDim1,
      Int gridDim2,
      const SparseMatrix<Field>& A,
            bool hermitian=true,
      const BisectCtrl& bisectCtrl=BisectCtrl() );

    void Factor( LDLFrontType frontType=LDL_2D );

    void ChangeNonzeroValues( const SparseMatrix<Field>& ANew );

    // Overwrite 'B' with the solution to 'A X = B'.
    void Solve( Matrix<Field>& B ) const;

    // Overwrite 'B' with the solution to 'A X = B' using Iterative Refinement.
    void SolveWithIterativeRefinement
    ( const SparseMatrix<Field>& A,
            Matrix<Field>& B,
      const Base<Field>& relTolRefine,
      Int maxRefineIts ) const;

    ldl::Front<Field>& Front();
    const ldl::Front<Field>& Front() const;

    ldl::NodeInfo& NodeInfo();
    const ldl::NodeInfo& NodeInfo() const;

    ldl::Separator& Separator();
    const ldl::Separator& Separator() const;

    vector<Int>& Map();
    const vector<Int>& Map() const;

    vector<Int>& InverseMap();
    const vector<Int>& InverseMap() const;
    
private:
    unique_ptr<ldl::Front<Field>> front_;
    unique_ptr<ldl::NodeInfo> info_;
    unique_ptr<ldl::Separator> separator_;

    vector<Int> map_, inverseMap_;
};

template<typename Field>
class DistSparseLDLFactorization
{
public:
    DistSparseLDLFactorization();

    void Initialize
    ( const DistSparseMatrix<Field>& A,
            bool hermitian=true,
      const BisectCtrl& bisectCtrl=BisectCtrl() );

    // TODO(poulson): Combine these two interfaces
    void Initialize2DGridGraph
    ( Int gridDim0,
      Int gridDim1,
      const DistSparseMatrix<Field>& A,
            bool hermitian=true,
      const BisectCtrl& bisectCtrl=BisectCtrl() );
    void Initialize3DGridGraph
    ( Int gridDim0,
      Int gridDim1,
      Int gridDim2,
      const DistSparseMatrix<Field>& A,
            bool hermitian=true,
      const BisectCtrl& bisectCtrl=BisectCtrl() );

    void Factor( LDLFrontType frontType=LDL_2D );

    void ChangeNonzeroValues( const DistSparseMatrix<Field>& ANew );

    // Overwrite 'B' with the solution to 'A X = B'.
    void Solve( DistMultiVec<Field>& B ) const;

    // Overwrite 'B' with the solution to 'A X = B' using Iterative Refinement.
    void SolveWithIterativeRefinement
    ( const DistSparseMatrix<Field>& A,
            DistMultiVec<Field>& B,
      const Base<Field>& relTolRefine,
      Int maxRefineIts ) const;

    ldl::DistFront<Field>& Front();
    const ldl::DistFront<Field>& Front() const;

    ldl::DistNodeInfo& NodeInfo();
    const ldl::DistNodeInfo& NodeInfo() const;

    ldl::DistSeparator& Separator();
    const ldl::DistSeparator& Separator() const;

    DistMap& Map();
    const DistMap& Map() const;

    DistMap& InverseMap();
    const DistMap& InverseMap() const;

    ldl::DistMultiVecNodeMeta& DistMultiVecNodeMeta() const;

private:
    unique_ptr<ldl::DistFront<Field>> front_;
    unique_ptr<ldl::DistNodeInfo> info_;
    unique_ptr<ldl::DistSeparator> separator_;

    DistMap map_, inverseMap_;

    // Metadata for repeated calls to DistFront<Field>::Pull
    mutable bool formedPullMetadata_=false;
    mutable vector<Int> mappedSources_, mappedTargets_, columnOffsets_;

    // Metadata for future use.
    mutable ldl::DistMultiVecNodeMeta dmvMeta_;
};

} // namespace El

#endif // ifndef EL_FACTOR_LDL_SPARSE_NUMERIC_HPP
