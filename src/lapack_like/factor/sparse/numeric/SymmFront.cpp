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

template<typename F>
SymmFront<F>::SymmFront( SymmFront<F>* parentNode )
: parent(parentNode), duplicate(nullptr)
{ }

template<typename F>
SymmFront<F>::SymmFront( DistSymmFront<F>* dupNode )
: parent(nullptr), duplicate(dupNode)
{
    isHermitian = dupNode->isHermitian;
    type = dupNode->type;
}

template<typename F>
SymmFront<F>::SymmFront
( const SparseMatrix<F>& A, 
  const vector<Int>& reordering,
  const SymmNodeInfo& info,
  bool conjugate )
: parent(nullptr), duplicate(nullptr)
{
    DEBUG_ONLY(CallStackEntry cse("SymmFront::SymmFront"))
    Pull( A, reordering, info, conjugate );
}

template<typename F>
void SymmFront<F>::Pull
( const SparseMatrix<F>& A, 
  const vector<Int>& reordering,
  const SymmNodeInfo& rootInfo,
  bool conjugate )
{
    DEBUG_ONLY(
        CallStackEntry cse("SymmFront::Pull");
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
    
    function<void(const SymmNodeInfo&,SymmFront<F>&)> pull = 
      [&]( const SymmNodeInfo& node, SymmFront<F>& front )
      {
        // Delete any existing children
        for( SymmFront<F>* child : front.children )
            delete child;

        const Int numChildren = node.children.size();
        front.children.resize( numChildren );
        for( Int c=0; c<numChildren; ++c )
        {
            front.children[c] = new SymmFront<F>(&front);
            pull( *node.children[c], *front.children[c] );
        }

        const Int lowerSize = node.lowerStruct.size();
        Zeros( front.L, node.size+lowerSize, node.size );

        for( Int t=0; t<node.size; ++t )
        {
            const Int j = invReorder[node.off+t];
            const Int numConn = A.NumConnections( j );
            const Int entryOff = A.EntryOffset( j );
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
void SymmFront<F>::Push
( SparseMatrix<F>& A, 
  const vector<Int>& reordering,
  const SymmNodeInfo& rootInfo ) const
{
    DEBUG_ONLY(CallStackEntry cse("SymmFront::Push"))

    // Invert the reordering
    const Int n = reordering.size();
    vector<Int> invReorder( n ); 
    for( Int j=0; j<n; ++j )
        invReorder[reordering[j]] = j;

    Zeros( A, n, n );

    // Reserve space for the lower triangle
    Int numLower = 0;
    function<void(const SymmFront<F>&)> countLower = 
      [&]( const SymmFront<F>& front )
      {
          for( const SymmFront<F>* child : front.children )
              countLower( *child );
          const Int nodeSize = front.L.Width();
          const Int structSize = front.L.Height() - nodeSize;
          numLower += (nodeSize*(nodeSize+1))/2 + nodeSize*structSize;
      };
    countLower( *this );
    A.Reserve( numLower );

    function<void(const SymmNodeInfo&,const SymmFront<F>&)> 
      push = 
      [&]( const SymmNodeInfo& node, const SymmFront<F>& front )
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
    A.MakeConsistent();
    MakeSymmetric( LOWER, A, isHermitian );
}

template<typename F>
void SymmFront<F>::Unpack
( SparseMatrix<F>& A, const SymmNodeInfo& rootInfo ) const
{
    DEBUG_ONLY(CallStackEntry cse("SymmFront::Push"))
    const Int n = rootInfo.off + rootInfo.size;
    Zeros( A, n, n );

    // Reserve space for the lower triangle
    Int numLower = 0;
    function<void(const SymmFront<F>&)> countLower = 
      [&]( const SymmFront<F>& front )
      {
          for( const SymmFront<F>* child : front.children )
              countLower( *child );
          const Int nodeSize = front.L.Width();
          const Int structSize = front.L.Height() - nodeSize;
          numLower += (nodeSize*(nodeSize+1))/2 + nodeSize*structSize;
      };
    countLower( *this );
    A.Reserve( numLower );

    function<void(const SymmNodeInfo&,const SymmFront<F>&)> 
      push = 
      [&]( const SymmNodeInfo& node, const SymmFront<F>& front )
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
    A.MakeConsistent();
}

template<typename F>
Int SymmFront<F>::NumEntries() const
{
    DEBUG_ONLY(CallStackEntry cse("SymmFront::NumEntries"))
    Int numEntries = 0;
    function<void(const SymmFront<F>&)> count =
      [&]( const SymmFront<F>& front )
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
Int SymmFront<F>::NumTopLeftEntries() const
{
    DEBUG_ONLY(CallStackEntry cse("SymmFront::NumTopLeftEntries"))
    Int numEntries = 0;
    function<void(const SymmFront<F>&)> count =
      [&]( const SymmFront<F>& front )
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
Int SymmFront<F>::NumBottomLeftEntries() const
{
    DEBUG_ONLY(CallStackEntry cse("SymmFront::NumBottomLeftEntries"))
    Int numEntries = 0;
    function<void(const SymmFront<F>&)> count =
      [&]( const SymmFront<F>& front )
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
double SymmFront<F>::FactorGFlops() const
{
    DEBUG_ONLY(CallStackEntry cse("DistSymmFront::FactorGFlops"))
    double gflops = 0.;
    function<void(const SymmFront<F>&)> count =
      [&]( const SymmFront<F>& front )
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
double SymmFront<F>::SolveGFlops( Int numRHS ) const
{
    DEBUG_ONLY(CallStackEntry cse("SymmFront::SolveGFlops"))
    double gflops = 0.;
    function<void(const SymmFront<F>&)> count =
      [&]( const SymmFront<F>& front )
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

#define PROTO(F) template class SymmFront<F>;
#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
