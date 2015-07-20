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
   
   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {
namespace ldl {

template<typename F>
Front<F>::Front( Front<F>* parentNode )
: parent(parentNode), duplicate(nullptr)
{ 
    if( parentNode != nullptr )
    {
        isHermitian = parentNode->isHermitian;
        type = parentNode->type;
    }
}

template<typename F>
Front<F>::Front( DistFront<F>* dupNode )
: parent(nullptr), duplicate(dupNode)
{
    isHermitian = dupNode->isHermitian;
    type = dupNode->type;
}

template<typename F>
Front<F>::Front
( const SparseMatrix<F>& A, 
  const vector<Int>& reordering,
  const NodeInfo& info,
  bool conjugate )
: parent(nullptr), duplicate(nullptr)
{
    DEBUG_ONLY(CSE cse("Front::Front"))
    Pull( A, reordering, info, conjugate );
}

template<typename F>
void Front<F>::Pull
( const SparseMatrix<F>& A, 
  const vector<Int>& reordering,
  const NodeInfo& rootInfo,
  bool conjugate )
{
    DEBUG_ONLY(
      CSE cse("Front::Pull");
      if( A.Height() != (Int)reordering.size() )
          LogicError("Mapping was not the right size");
    )
    type = SYMM_2D;
    isHermitian = conjugate;

    // Invert the reordering
    const Int n = reordering.size();
    DEBUG_ONLY(
      if( A.Height() != n || A.Width() != n )
          LogicError("Expected A to be square and consistent with reordering");
    )
    vector<Int> invReorder(n);
    for( Int j=0; j<n; ++j )
        invReorder[reordering[j]] = j;
    
    function<void(const NodeInfo&,Front<F>&)> pull = 
      [&]( const NodeInfo& node, Front<F>& front )
      {
        // Delete any existing children
        for( Front<F>* child : front.children )
            delete child;

        const Int numChildren = node.children.size();
        front.children.resize( numChildren );
        for( Int c=0; c<numChildren; ++c )
        {
            front.children[c] = new Front<F>(&front);
            pull( *node.children[c], *front.children[c] );
        }

        const Int lowerSize = node.lowerStruct.size();
        Zeros( front.L, node.size+lowerSize, node.size );

        for( Int t=0; t<node.size; ++t )
        {
            const Int j = invReorder[node.off+t];
            const Int numConn = A.NumConnections( j );
            const Int entryOff = A.RowOffset( j );
            for( Int k=0; k<numConn; ++k )
            {
                const Int iOrig = A.Col( entryOff+k );
                const Int i = reordering[iOrig];
                const F transVal = A.Value( entryOff+k );
                const F value = ( conjugate ? Conj(transVal) : transVal );

                if( i < node.off+t )
                    continue;
                else if( i < node.off+node.size )
                {
                    front.L.Set( i-node.off, t, value );
                }
                else
                {
                    const Int origOff = Find( node.origLowerStruct, i );
                    const Int row = node.origLowerRelInds[origOff];
                    DEBUG_ONLY(
                        if( row < t )
                            LogicError("Tried to touch upper triangle");
                    )
                    front.L.Set( row, t, value );
                }
            }
        }
      };
    pull( rootInfo, *this );
}

template<typename F>
void Front<F>::PullUpdate
( const SparseMatrix<F>& A, 
  const vector<Int>& reordering,
  const NodeInfo& rootInfo )
{
    DEBUG_ONLY(
      CSE cse("Front::PullUpdate");
      if( A.Height() != (Int)reordering.size() )
          LogicError("Mapping was not the right size");
    )

    // Invert the reordering
    const Int n = reordering.size();
    DEBUG_ONLY(
      if( A.Height() != n || A.Width() != n )
          LogicError("Expected A to be square and consistent with reordering");
    )
    vector<Int> invReorder(n);
    for( Int j=0; j<n; ++j )
        invReorder[reordering[j]] = j;
    
    function<void(const NodeInfo&,Front<F>&)> pull = 
      [&]( const NodeInfo& node, Front<F>& front )
      {
        const Int numChildren = node.children.size();
        for( Int c=0; c<numChildren; ++c )
            pull( *node.children[c], *front.children[c] );

        for( Int t=0; t<node.size; ++t )
        {
            const Int j = invReorder[node.off+t];
            const Int numConn = A.NumConnections( j );
            const Int entryOff = A.RowOffset( j );
            for( Int k=0; k<numConn; ++k )
            {
                const Int iOrig = A.Col( entryOff+k );
                const Int i = reordering[iOrig];
                const F transVal = A.Value( entryOff+k );
                const F value = ( isHermitian ? Conj(transVal) : transVal );

                if( i < node.off+t )
                    continue;
                else if( i < node.off+node.size )
                {
                    front.L.Update( i-node.off, t, value );
                }
                else
                {
                    const Int origOff = Find( node.origLowerStruct, i );
                    const Int row = node.origLowerRelInds[origOff];
                    DEBUG_ONLY(
                        if( row < t )
                            LogicError("Tried to touch upper triangle");
                    )
                    front.L.Update( row, t, value );
                }
            }
        }
      };
    pull( rootInfo, *this );
}


template<typename F>
void Front<F>::Push
( SparseMatrix<F>& A, 
  const vector<Int>& reordering,
  const NodeInfo& rootInfo ) const
{
    DEBUG_ONLY(CSE cse("Front::Push"))

    // Invert the reordering
    const Int n = reordering.size();
    vector<Int> invReorder( n ); 
    for( Int j=0; j<n; ++j )
        invReorder[reordering[j]] = j;

    Zeros( A, n, n );

    // Reserve space for the lower triangle
    Int numLower = 0;
    function<void(const Front<F>&)> countLower = 
      [&]( const Front<F>& front )
      {
          for( const Front<F>* child : front.children )
              countLower( *child );
          const Int nodeSize = front.L.Width();
          const Int structSize = front.L.Height() - nodeSize;
          numLower += (nodeSize*(nodeSize+1))/2 + nodeSize*structSize;
      };
    countLower( *this );
    A.Reserve( numLower );

    function<void(const NodeInfo&,const Front<F>&)> 
      push = 
      [&]( const NodeInfo& node, const Front<F>& front )
      {
        const Int numChildren = node.children.size();
        for( Int c=0; c<numChildren; ++c )
            push( *node.children[c], *front.children[c] );

        const Int lowerSize = node.lowerStruct.size();
        for( Int t=0; t<node.size; ++t )
        {
            const Int j = invReorder[node.off+t];

            // Push in the lower triangle of the diagonal block
            for( Int s=t; s<node.size; ++s )
            {
                const Int i = invReorder[node.off+s];
                A.QueueUpdate( i, j, front.L.Get(s,t) );
            }

            // Push in the connectivity 
            for( Int s=0; s<lowerSize; ++s )
            {
                const Int i = invReorder[node.lowerStruct[s]];
                A.QueueUpdate( i, j, front.L.Get(s+node.size,t) );
            }
        }
      };
    push( rootInfo, *this );
    A.ProcessQueues();
    MakeSymmetric( LOWER, A, isHermitian );
}

template<typename F>
void Front<F>::Unpack( SparseMatrix<F>& A, const NodeInfo& rootInfo ) const
{
    DEBUG_ONLY(CSE cse("Front::Push"))
    const Int n = rootInfo.off + rootInfo.size;
    Zeros( A, n, n );

    // Reserve space for the lower triangle
    Int numLower = 0;
    function<void(const Front<F>&)> countLower = 
      [&]( const Front<F>& front )
      {
          for( const Front<F>* child : front.children )
              countLower( *child );
          const Int nodeSize = front.L.Width();
          const Int structSize = front.L.Height() - nodeSize;
          numLower += (nodeSize*(nodeSize+1))/2 + nodeSize*structSize;
      };
    countLower( *this );
    A.Reserve( numLower );

    function<void(const NodeInfo&,const Front<F>&)> push = 
      [&]( const NodeInfo& node, const Front<F>& front )
      {
        const Int numChildren = node.children.size();
        for( Int c=0; c<numChildren; ++c )
            push( *node.children[c], *front.children[c] );

        const Int lowerSize = node.lowerStruct.size();
        for( Int t=0; t<node.size; ++t )
        {
            const Int j = node.off+t;

            // Push in the lower triangle of the diagonal block
            for( Int s=t; s<node.size; ++s )
            {
                const Int i = node.off+s;
                A.QueueUpdate( i, j, front.L.Get(s,t) );
            }

            // Push in the connectivity 
            for( Int s=0; s<lowerSize; ++s )
            {
                const Int i = node.lowerStruct[s];
                A.QueueUpdate( i, j, front.L.Get(s+node.size,t) );
            }
        }
      };
    push( rootInfo, *this );
    A.ProcessQueues();
}

template<typename F>
const Front<F>& Front<F>::operator=( const Front<F>& front )
{
    DEBUG_ONLY(CSE cse("Front::operator="))
    isHermitian = front.isHermitian;
    type = front.type;
    L = front.L;
    diag = front.diag;
    subdiag = front.subdiag;
    piv = front.piv;
    work = front.work;
    // Do not copy parent...
    const int numChildren = front.children.size();
    children.resize( numChildren );
    for( int c=0; c<numChildren; ++c )
    {
        children[c] = new Front<F>(this);
        *children[c] = *front.children[c];
    }
    return *this;
}

template<typename F>
Int Front<F>::NumEntries() const
{
    DEBUG_ONLY(CSE cse("Front::NumEntries"))
    Int numEntries = 0;
    function<void(const Front<F>&)> count =
      [&]( const Front<F>& front )
      {
        for( auto* child : front.children )
            count( *child );

        // Add in L
        numEntries += front.L.Height() * front.L.Width();
 
        // Add in the workspace
        numEntries += front.work.Height()*front.work.Width(); 
      };
    count( *this );
    return numEntries;
}

template<typename F>
Int Front<F>::NumTopLeftEntries() const
{
    DEBUG_ONLY(CSE cse("Front::NumTopLeftEntries"))
    Int numEntries = 0;
    function<void(const Front<F>&)> count =
      [&]( const Front<F>& front )
      {
        for( auto* child : front.children )
            count( *child );
        const Int n = front.L.Width();
        numEntries += n*n;
      };
    count( *this );
    return numEntries;
}

template<typename F>
Int Front<F>::NumBottomLeftEntries() const
{
    DEBUG_ONLY(CSE cse("Front::NumBottomLeftEntries"))
    Int numEntries = 0;
    function<void(const Front<F>&)> count =
      [&]( const Front<F>& front )
      {
        for( auto* child : front.children )
            count( *child );
        const Int m = front.L.Height();
        const Int n = front.L.Width();
        numEntries += (m-n)*n;
      };
    count( *this );
    return numEntries;
}

template<typename F>
double Front<F>::FactorGFlops() const
{
    DEBUG_ONLY(CSE cse("DistFront::FactorGFlops"))
    double gflops = 0.;
    function<void(const Front<F>&)> count =
      [&]( const Front<F>& front )
      {
        for( auto* child : front.children )
            count( *child );
        const double m = front.L.Height();
        const double n = front.L.Width();
        const double realFrontFlops = (n*n*n/3) + (m-n)*n + (m-n)*(m-n)*n;
        gflops += (IsComplex<F>::val ? 4*realFrontFlops : realFrontFlops)/1.e9;
      };
    count( *this );
    return gflops;
}

template<typename F>
double Front<F>::SolveGFlops( Int numRHS ) const
{
    DEBUG_ONLY(CSE cse("Front::SolveGFlops"))
    double gflops = 0.;
    function<void(const Front<F>&)> count =
      [&]( const Front<F>& front )
      {
        for( auto* child : front.children )
            count( *child );
        const double m = front.L.Height();
        const double n = front.L.Width();
        const double realFrontFlops = m*n*numRHS;
        gflops += (IsComplex<F>::val ? 4*realFrontFlops : realFrontFlops)/1.e9;
      };
    count( *this );
    return gflops;
}

#define PROTO(F) template struct Front<F>;
#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace ldl
} // namespace El
