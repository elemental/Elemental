/*
   Copyright (c) 2009-2012, Jack Poulson, Lexing Ying, and
   The University of Texas at Austin.
   All rights reserved.

   Copyright (c) 2013, Jack Poulson, Lexing Ying, and Stanford University.
   All rights reserved.

   Copyright (c) 2013-2014, Jack Poulson and
   The Georgia Institute of Technology.
   All rights reserved.

   Copyright (c) 2014-2015, Jack Poulson and Stanford University.
   All rights reserved.

   Copyright (c) 2016, Jack Poulson.
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {
namespace ldl {

template<typename Field>
Front<Field>::Front( Front<Field>* parentNode )
: parent(parentNode)
{
    if( parentNode != nullptr )
    {
        isHermitian = parentNode->isHermitian;
        type = parentNode->type;
    }
}

template<typename Field>
Front<Field>::Front( DistFront<Field>* dupNode )
: duplicate(dupNode)
{
    isHermitian = dupNode->isHermitian;
    type = dupNode->type;
}

template<typename Field>
Front<Field>::Front
( const SparseMatrix<Field>& A,
  const vector<Int>& reordering,
  const NodeInfo& info,
  bool conjugate )
{
    EL_DEBUG_CSE
    Pull( A, reordering, info, conjugate );
}

template<typename Field>
Front<Field>::~Front() { }

template<typename Field>
void Front<Field>::Pull
( const SparseMatrix<Field>& A,
  const vector<Int>& reordering,
  const NodeInfo& rootInfo,
  bool conjugate )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( A.Height() != (Int)reordering.size() )
          LogicError("Mapping was not the right size");
    )
    type = SYMM_2D;
    isHermitian = conjugate;

    // Invert the reordering
    const Int n = reordering.size();
    EL_DEBUG_ONLY(
      if( A.Height() != n || A.Width() != n )
          LogicError("Expected A to be square and consistent with reordering");
    )
    vector<Int> invReorder(n);
    for( Int j=0; j<n; ++j )
        invReorder[reordering[j]] = j;

    function<void(const NodeInfo&,Front<Field>&)> pull =
      [&]( const NodeInfo& node, Front<Field>& front )
      {
        // Delete any existing children
        SwapClear( front.children );

        const Int numChildren = node.children.size();
        front.children.resize( numChildren );
        for( Int c=0; c<numChildren; ++c )
        {
            front.children[c].reset( new Front<Field>(&front) );
            pull( *node.children[c], *front.children[c] );
        }
        // Mark this node as a sparse leaf if it does not have any children
        if( numChildren == 0 )
            front.sparseLeaf = true;

        const Int lowerSize = node.lowerStruct.size();
        const Field* AValBuf = A.LockedValueBuffer();
        const Int* AColBuf = A.LockedTargetBuffer();
        const Int* AOffsetBuf = A.LockedOffsetBuffer();
        if( front.sparseLeaf )
        {
            front.workSparse.Empty();
            Zeros( front.workSparse, node.size, node.size );
            Zeros( front.LDense, lowerSize, node.size );

            // Count the number of sparse entries to queue into the top-left
            Int numEntriesTopLeft = 0;
            for( Int t=0; t<node.size; ++t )
            {
                const Int j = invReorder[node.off+t];
                const Int rowOff = AOffsetBuf[j];
                const Int numConn = AOffsetBuf[j+1] - rowOff;
                for( Int k=0; k<numConn; ++k )
                {
                    const Int iOrig = AColBuf[rowOff+k];
                    const Int i = reordering[iOrig];

                    if( i < node.off+t )
                        continue;
                    else if( i < node.off+node.size )
                        ++numEntriesTopLeft;
                }
            }
            front.workSparse.Reserve( numEntriesTopLeft );

            for( Int t=0; t<node.size; ++t )
            {
                const Int j = invReorder[node.off+t];
                const Int rowOff = AOffsetBuf[j];
                const Int numConn = AOffsetBuf[j+1] - rowOff;
                for( Int k=0; k<numConn; ++k )
                {
                    const Int iOrig = AColBuf[rowOff+k];
                    const Int i = reordering[iOrig];

                    const Field transVal = AValBuf[rowOff+k];
                    const Field value = conjugate ? Conj(transVal) : transVal;

                    if( i < node.off+t )
                        continue;
                    else if( i < node.off+node.size )
                    {
                        // Since SuiteSparse makes use of column-major ordering,
                        // and Elemental uses row-major ordering of its sparse
                        // matrices, we are implicitly storing the transpose.
                        front.workSparse.QueueUpdate( i-node.off, t, transVal );
                    }
                    else
                    {
                        const Int origOff = Find( node.origLowerStruct, i );
                        const Int row = node.origLowerRelInds[origOff];
                        EL_DEBUG_ONLY(
                          if( row < t )
                              LogicError("Tried to touch upper triangle");
                        )
                        front.LDense(row-node.size,t) = value;
                    }
                }
            }
            front.workSparse.ProcessQueues();
            MakeSymmetric( LOWER, front.workSparse, front.isHermitian );
        }
        else
        {
            Zeros( front.LDense, node.size+lowerSize, node.size );
            for( Int t=0; t<node.size; ++t )
            {
                const Int j = invReorder[node.off+t];
                const Int rowOff = AOffsetBuf[j];
                const Int numConn = AOffsetBuf[j+1] - rowOff;
                for( Int k=0; k<numConn; ++k )
                {
                    const Int iOrig = AColBuf[rowOff+k];
                    const Int i = reordering[iOrig];

                    const Field transVal = AValBuf[rowOff+k];
                    const Field value = conjugate ? Conj(transVal) : transVal;

                    if( i < node.off+t )
                        continue;
                    else if( i < node.off+node.size )
                    {
                        front.LDense(i-node.off,t) = value;
                    }
                    else
                    {
                        const Int origOff = Find( node.origLowerStruct, i );
                        const Int row = node.origLowerRelInds[origOff];
                        EL_DEBUG_ONLY(
                          if( row < t )
                              LogicError("Tried to touch upper triangle");
                        )
                        front.LDense(row,t) = value;
                    }
                }
            }
        }
      };
    pull( rootInfo, *this );
}

template<typename Field>
void Front<Field>::PullUpdate
( const SparseMatrix<Field>& A,
  const vector<Int>& reordering,
  const NodeInfo& rootInfo )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( A.Height() != (Int)reordering.size() )
          LogicError("Mapping was not the right size");
    )

    // Invert the reordering
    const Int n = reordering.size();
    EL_DEBUG_ONLY(
      if( A.Height() != n || A.Width() != n )
          LogicError("Expected A to be square and consistent with reordering");
    )
    vector<Int> invReorder(n);
    for( Int j=0; j<n; ++j )
        invReorder[reordering[j]] = j;

    function<void(const NodeInfo&,Front<Field>&)> pull =
      [&]( const NodeInfo& node, Front<Field>& front )
      {
        const Int numChildren = node.children.size();
        for( Int c=0; c<numChildren; ++c )
            pull( *node.children[c], *front.children[c] );

        const Field* AValBuf = A.LockedValueBuffer();
        const Int* AColBuf = A.LockedTargetBuffer();
        const Int* AOffsetBuf = A.LockedOffsetBuffer();

        if( front.sparseLeaf )
        {
            LogicError("Sparse leaves not yet handled in Front::PullUpdate");
        }
        else
        {
            for( Int t=0; t<node.size; ++t )
            {
                const Int j = invReorder[node.off+t];
                const Int rowOff = AOffsetBuf[j];
                const Int numConn = AOffsetBuf[j+1] - rowOff;
                for( Int k=0; k<numConn; ++k )
                {
                    const Int iOrig = AColBuf[rowOff+k];
                    const Int i = reordering[iOrig];

                    const Field transVal = AValBuf[rowOff+k];
                    const Field value = isHermitian ? Conj(transVal) : transVal;

                    if( i < node.off+t )
                        continue;
                    else if( i < node.off+node.size )
                    {
                        front.LDense(i-node.off,t) += value;
                    }
                    else
                    {
                        const Int origOff = Find( node.origLowerStruct, i );
                        const Int row = node.origLowerRelInds[origOff];
                        EL_DEBUG_ONLY(
                          if( row < t )
                              LogicError("Tried to touch upper triangle");
                        )
                        front.LDense(row,t) += value;
                    }
                }
            }
        }
      };
    pull( rootInfo, *this );
}

// TODO(poulson): Use lower-level access
template<typename Field>
void Front<Field>::Push
( SparseMatrix<Field>& A,
  const vector<Int>& reordering,
  const NodeInfo& rootInfo ) const
{
    EL_DEBUG_CSE

    // Invert the reordering
    const Int n = reordering.size();
    vector<Int> invReorder( n );
    for( Int j=0; j<n; ++j )
        invReorder[reordering[j]] = j;

    Zeros( A, n, n );

    // Reserve space for the lower triangle
    // TODO(poulson): Compensate for sparse leaves
    Int numLower = 0;
    function<void(const Front<Field>&)> countLower =
      [&]( const Front<Field>& front )
      {
          for( const auto& child : front.children )
              countLower( *child );
          const Int nodeSize = front.LDense.Width();
          const Int structSize = front.Height() - nodeSize;
          numLower += (nodeSize*(nodeSize+1))/2 + nodeSize*structSize;
      };
    countLower( *this );
    A.Reserve( numLower );

    function<void(const NodeInfo&,const Front<Field>&)>
      push =
      [&]( const NodeInfo& node, const Front<Field>& front )
      {
        const Int numChildren = node.children.size();
        for( Int c=0; c<numChildren; ++c )
            push( *node.children[c], *front.children[c] );

        const Int lowerSize = node.lowerStruct.size();
        if( front.sparseLeaf )
        {
            // Push in the diagonal block
            const Int numEntries = front.LSparse.NumEntries();
            if( numEntries == 0 )
            {
                // We have not yet factored, so use the original sparse matrix
                // stored in front.workSparse
                const Int numWorkEntries = front.workSparse.NumEntries();
                for( Int e=0; e<numWorkEntries; ++e )
                {
                    const Int iSparse = front.workSparse.Row(e);
                    const Int jSparse = front.workSparse.Col(e);
                    const Field value = front.workSparse.Value(e);
                    if( iSparse < jSparse || value == Field(0) )
                        continue;

                    const Int j = invReorder[jSparse + node.off];
                    A.QueueUpdate( invReorder[iSparse+node.off], j, value );
                }
            }
            else
            {
                // We have already factored, so use front.LSparse
                for( Int e=0; e<numEntries; ++e )
                {
                    const Int iSparse = front.LSparse.Row(e);
                    const Int jSparse = front.LSparse.Col(e);
                    const Field value = front.LSparse.Value(e);
                    if( iSparse < jSparse || value == Field(0) )
                        continue;

                    const Int j = invReorder[jSparse + node.off];
                    A.QueueUpdate( invReorder[iSparse+node.off], j, value );
                }
            }

            // Push in the lower connectivity
            for( Int t=0; t<node.size; ++t )
            {
                const Int j = invReorder[node.off+t];
                for( Int s=0; s<lowerSize; ++s )
                {
                    const Int i = invReorder[node.lowerStruct[s]];
                    const Field value = front.LDense(s,t);
                    if( value != Field(0) )
                        A.QueueUpdate( i, j, value );
                }
            }
        }
        else
        {
            for( Int t=0; t<node.size; ++t )
            {
                const Int j = invReorder[node.off+t];

                // Push in the lower triangle of the diagonal block
                for( Int s=t; s<node.size; ++s )
                {
                    const Int i = invReorder[node.off+s];
                    const Field value = front.LDense(s,t);
                    if( value != Field(0) )
                        A.QueueUpdate( i, j, value );
                }

                // Push in the connectivity
                for( Int s=0; s<lowerSize; ++s )
                {
                    const Int i = invReorder[node.lowerStruct[s]];
                    const Field value = front.LDense(s+node.size,t);
                    if( value != Field(0) )
                        A.QueueUpdate( i, j, value );
                }
            }
        }
      };
    push( rootInfo, *this );
    A.ProcessQueues();
    MakeSymmetric( LOWER, A, isHermitian );
}

// TODO(poulson): Use lower-level access
template<typename Field>
void Front<Field>::Unpack
( SparseMatrix<Field>& A, const NodeInfo& rootInfo ) const
{
    EL_DEBUG_CSE
    const Int n = rootInfo.off + rootInfo.size;
    Zeros( A, n, n );

    // Reserve space for the lower triangle
    // TODO(poulson): Compensate for sparse leaves
    Int numLower = 0;
    function<void(const Front<Field>&)> countLower =
      [&]( const Front<Field>& front )
      {
          for( const auto& child : front.children )
              countLower( *child );
          const Int nodeSize = front.LDense.Width();
          const Int structSize = front.Height() - nodeSize;
          numLower += (nodeSize*(nodeSize+1))/2 + nodeSize*structSize;
      };
    countLower( *this );
    A.Reserve( numLower );

    function<void(const NodeInfo&,const Front<Field>&)> push =
      [&]( const NodeInfo& node, const Front<Field>& front )
      {
        const Int numChildren = node.children.size();
        for( Int c=0; c<numChildren; ++c )
        {
            PushIndent();
            push( *node.children[c], *front.children[c] );
            PopIndent();
        }

        const Int lowerSize = node.lowerStruct.size();
        if( front.sparseLeaf )
        {
            // Push in the diagonal block
            const Int numEntries = front.LSparse.NumEntries();
            if( numEntries == 0 )
            {
                // We have not yet factored, so use the original sparse matrix
                // stored in front.workSparse
                const Int numWorkEntries = front.workSparse.NumEntries();
                for( Int e=0; e<numWorkEntries; ++e )
                {
                    const Int iSparse = front.workSparse.Row(e);
                    const Int jSparse = front.workSparse.Col(e);
                    const Field value = front.workSparse.Value(e);
                    if( iSparse < jSparse || value == Field(0) )
                        continue;

                    const Int j = jSparse + node.off;
                    A.QueueUpdate( iSparse+node.off, j, value );
                }
            }
            else
            {
                // We have already factored, so use front.LSparse
                for( Int e=0; e<numEntries; ++e )
                {
                    const Int iSparse = front.LSparse.Row(e);
                    const Int jSparse = front.LSparse.Col(e);
                    const Field value = front.LSparse.Value(e);
                    if( iSparse < jSparse || value == Field(0) )
                        continue;

                    const Int j = jSparse + node.off;
                    A.QueueUpdate( iSparse+node.off, j, value );
                }
            }

            // Push in the lower connectivity
            for( Int t=0; t<node.size; ++t )
            {
                const Int j = node.off+t;
                for( Int s=0; s<lowerSize; ++s )
                {
                    const Int i = node.lowerStruct[s];
                    const Field value = front.LDense(s,t);
                    if( value != Field(0) )
                        A.QueueUpdate( i, j, value );
                }
            }
        }
        else
        {
            for( Int t=0; t<node.size; ++t )
            {
                const Int j = node.off+t;

                // Push in the lower triangle of the diagonal block
                for( Int s=t; s<node.size; ++s )
                {
                    const Int i = node.off+s;
                    const Field value = front.LDense(s,t);
                    if( value != Field(0) )
                        A.QueueUpdate( i, j, value );
                }

                // Push in the connectivity
                for( Int s=0; s<lowerSize; ++s )
                {
                    const Int i = node.lowerStruct[s];
                    const Field value = front.LDense(s+node.size,t);
                    if( value != Field(0) )
                        A.QueueUpdate( i, j, value );
                }
            }
        }
      };
    push( rootInfo, *this );
    A.ProcessQueues();
}

template<typename Field>
const Front<Field>& Front<Field>::operator=( const Front<Field>& front )
{
    EL_DEBUG_CSE
    isHermitian = front.isHermitian;
    sparseLeaf = front.sparseLeaf;
    type = front.type;
    LDense = front.LDense;
    LSparse = front.LSparse;
    diag = front.diag;
    subdiag = front.subdiag;
    p = front.p;
    workDense = front.workDense;
    workSparse = front.workSparse;
    // Do not copy parent...

    // Delete any existing children
    SwapClear( children );

    const int numChildren = front.children.size();
    children.resize( numChildren );
    for( int c=0; c<numChildren; ++c )
    {
        children[c].reset( new Front<Field>(this) );
        *children[c] = *front.children[c];
    }
    return *this;
}

template<typename Field>
Int Front<Field>::Height() const
{ return sparseLeaf ? LDense.Height()+LDense.Width() : LDense.Height(); }

template<typename Field>
Int Front<Field>::NumEntries() const
{
    EL_DEBUG_CSE
    Int numEntries = 0;
    function<void(const Front<Field>&)> count =
      [&]( const Front<Field>& front )
      {
        for( const auto& child : front.children )
            count( *child );

        if( front.sparseLeaf )
        {
            // Count the diagonal block
            const Int numSparseEntries = front.LSparse.NumEntries();
            if( numSparseEntries == 0 )
            {
                // The matrix has not yet been factored
                numEntries += front.workSparse.NumEntries();
            }
            else
            {
                // The matrix has already been factored
                numEntries += numSparseEntries;
            }

            // Count the connectivity
            numEntries += front.LDense.Height() * front.LDense.Width();
        }
        else
        {
            // Add in L
            numEntries += front.LDense.Height() * front.LDense.Width();
        }
        // Add in the workspace for the Schur complement
        numEntries += front.workDense.Height()*front.workDense.Width();
      };
    count( *this );
    return numEntries;
}

template<typename Field>
Int Front<Field>::NumTopLeftEntries() const
{
    EL_DEBUG_CSE
    Int numEntries = 0;
    function<void(const Front<Field>&)> count =
      [&]( const Front<Field>& front )
      {
        for( const auto& child : front.children )
            count( *child );
        if( front.sparseLeaf )
        {
            // Count the diagonal block
            const Int numSparseEntries = front.LSparse.NumEntries();
            if( numSparseEntries == 0 )
            {
                // The matrix has not yet been factored
                numEntries += front.workSparse.NumEntries();
            }
            else
            {
                // The matrix has already been factored
                numEntries += numSparseEntries;
            }
        }
        else
        {
            const Int n = front.LDense.Width();
            numEntries += n*n;
        }
      };
    count( *this );
    return numEntries;
}

template<typename Field>
Int Front<Field>::NumBottomLeftEntries() const
{
    EL_DEBUG_CSE
    Int numEntries = 0;
    function<void(const Front<Field>&)> count =
      [&]( const Front<Field>& front )
      {
        for( const auto& child : front.children )
            count( *child );
        const Int m = front.LDense.Height();
        const Int n = front.LDense.Width();
        if( front.sparseLeaf )
        {
            numEntries += m*n;
        }
        else
        {
            numEntries += (m-n)*n;
        }
      };
    count( *this );
    return numEntries;
}

template<typename Field>
double Front<Field>::FactorGFlops() const
{
    EL_DEBUG_CSE
    double gflops = 0.;
    function<void(const Front<Field>&)> count =
      [&]( const Front<Field>& front )
      {
        for( const auto& child : front.children )
            count( *child );
        const double m = front.LDense.Height();
        const double n = front.LDense.Width();
        double realFrontFlops=0;
        if( front.sparseLeaf )
        {
            if( Unfactored(front.type) )
                LogicError("Matrix has not yet been factored");

            // Count the flops from the sparse factorization
            const Int* offsetBuf = front.LSparse.LockedOffsetBuffer();
            for( Int j=0; j<n; ++j )
            {
                const Int nnz = offsetBuf[j+1]-offsetBuf[j];
                realFrontFlops += nnz*(nnz+2.);
            }
            // Count the flops from the dense trsm
            realFrontFlops += m*n;
            // Count the flops from the Schur-complement update
            realFrontFlops += m*m*n;
        }
        else
        {
            realFrontFlops = (n*n*n/3) + (m-n)*n + (m-n)*(m-n)*n;
        }
        gflops += (IsComplex<Field>::value ? 4*realFrontFlops
                                       : realFrontFlops)/1.e9;
      };
    count( *this );
    return gflops;
}

template<typename Field>
double Front<Field>::SolveGFlops( Int numRHS ) const
{
    EL_DEBUG_CSE
    double gflops = 0.;
    function<void(const Front<Field>&)> count =
      [&]( const Front<Field>& front )
      {
        for( const auto& child : front.children )
            count( *child );
        const double m = front.LDense.Height();
        const double n = front.LDense.Width();
        double realFrontFlops = 0;
        if( front.sparseLeaf )
        {
            if( Unfactored(front.type) )
                LogicError("Matrix has not yet been factored");

            const double numEntries = front.LSparse.NumEntries();
            realFrontFlops = (numEntries+m*n)*numRHS;
        }
        else
        {
            realFrontFlops = m*n*numRHS;
        }
        gflops += (IsComplex<Field>::value ? 4*realFrontFlops
                                       : realFrontFlops)/1.e9;
      };
    count( *this );
    return gflops;
}

#define PROTO(Field) template struct Front<Field>;
#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace ldl
} // namespace El
