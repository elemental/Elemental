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
#include "El.hpp"

namespace El {

template<typename T>
DistMultiVecNode<T>::DistMultiVecNode( DistMultiVecNode<T>* parentNode )
: parent(parentNode), child(nullptr), duplicate(nullptr)
{ }

template<typename T>
DistMultiVecNode<T>::DistMultiVecNode
( const DistMap& invMap, const DistSymmNodeInfo& info,
  const DistMultiVec<T>& X )
: parent(nullptr), child(nullptr), duplicate(nullptr)
{
    DEBUG_ONLY(CallStackEntry cse("DistMultiVecNode::DistMultiVecNode"))
    Pull( invMap, info, X );
}

template<typename T>
DistMultiVecNode<T>::DistMultiVecNode( const DistMatrixNode<T>& X )
: parent(nullptr), child(nullptr), duplicate(nullptr)
{
    DEBUG_ONLY(CallStackEntry cse("DistMultiVecNode::DistMultiVecNode"))
    *this = X;
}

template<typename T>
const DistMultiVecNode<T>&
DistMultiVecNode<T>::operator=( const DistMatrixNode<T>& X )
{
    DEBUG_ONLY(CallStackEntry cse("DistMultiVecNode::operator="))

    if( X.child == nullptr )
    {
        delete duplicate; 
        duplicate = new MatrixNode<T>(this);
        *duplicate = *X.duplicate;

        matrix.Attach( X.matrix.Grid(), duplicate->matrix );

        return *this;
    }

    matrix.SetGrid( X.matrix.Grid() );
    matrix = X.matrix;

    delete child;
    child = new DistMultiVecNode<T>(this);
    *child = *X.child;

    return *this;
}

template<typename T>
void DistMultiVecNode<T>::Pull
( const DistMap& invMap, const DistSymmNodeInfo& info,
  const DistMultiVec<T>& X )
{
    DEBUG_ONLY(CallStackEntry cse("DistMultiVecNode::Pull"))
    const Int width = X.Width();

    // Count the number of indices assigned to this process
    // (and build the subtree)
    int numRecvInds = 0;
    function<void(const SymmNodeInfo&,MatrixNode<T>&)> localCount =
      [&]( const SymmNodeInfo& node, MatrixNode<T>& XNode )
      {
        const Int numChildren = node.children.size();
        XNode.children.resize(numChildren);
        for( Int c=0; c<numChildren; ++c )
        {
            XNode.children[c] = new MatrixNode<T>(&XNode);
            localCount( *node.children[c], *XNode.children[c] );
        }

        XNode.matrix.Resize( node.size, width );
        numRecvInds += node.size;
      };
    function<void(const DistSymmNodeInfo&,DistMultiVecNode<T>&)> count =
      [&]( const DistSymmNodeInfo& node, DistMultiVecNode<T>& XNode )
      {
        if( node.child == nullptr )
        {
            delete XNode.duplicate;
            XNode.duplicate = new MatrixNode<T>(&XNode);
            localCount( *node.duplicate, *XNode.duplicate );

            XNode.matrix.Attach( *node.grid, XNode.duplicate->matrix );

            return;
        }
        XNode.child = new DistMultiVecNode<T>(&XNode);
        count( *node.child, *XNode.child );

        XNode.commMeta.Empty();
        XNode.matrix.SetGrid( *node.grid );
        XNode.matrix.Resize( node.size, width );
        numRecvInds += XNode.matrix.LocalHeight();
      };
    count( info, *this );
   
    // Pack the indices for mapping to the original ordering
    Int off = 0;
    vector<Int> mappedInds( numRecvInds );
    function<void(const SymmNodeInfo&,MatrixNode<T>&)> localPack =
      [&]( const SymmNodeInfo& node, MatrixNode<T>& XNode )
      {
        const Int numChildren = node.children.size();
        for( Int c=0; c<numChildren; ++c )
            localPack( *node.children[c], *XNode.children[c] );

        for( Int t=0; t<node.size; ++t )
            mappedInds[off++] = node.off+t;
      };
    function<void(const DistSymmNodeInfo&,DistMultiVecNode<T>&)> pack =
      [&]( const DistSymmNodeInfo& node, DistMultiVecNode<T>& XNode )
      {
        if( node.child == nullptr )
        {
            localPack( *node.duplicate, *XNode.duplicate );
            return;
        }
        pack( *node.child, *XNode.child );

        for( Int tLoc=0; tLoc<XNode.matrix.LocalHeight(); ++tLoc )
        {
            const Int t = XNode.matrix.GlobalRow(tLoc);
            mappedInds[off++] = node.off+t;
        }
      };
    pack( info, *this );

    // Convert the indices to the original ordering
    invMap.Translate( mappedInds );

    // Figure out how many entries each process owns that we need
    mpi::Comm comm = X.Comm();
    const int commSize = mpi::Size( comm );
    vector<int> recvSizes( commSize, 0 );
    for( int s=0; s<numRecvInds; ++s )
        ++recvSizes[ X.RowOwner( mappedInds[s] ) ];
    vector<int> recvOffs;
    Scan( recvSizes, recvOffs );
    vector<Int> recvInds( numRecvInds );
    auto offs = recvOffs;
    for( int s=0; s<numRecvInds; ++s )
    {
        const Int i = mappedInds[s];
        recvInds[ offs[X.RowOwner(i)]++ ] = i;
    }

    // Coordinate for the coming AllToAll to exchange the indices of X
    vector<int> sendSizes(commSize);
    mpi::AllToAll( recvSizes.data(), 1, sendSizes.data(), 1, comm );
    vector<int> sendOffs;
    const int numSendInds = Scan( sendSizes, sendOffs );

    // Request the indices
    vector<Int> sendInds( numSendInds );
    mpi::AllToAll
    ( recvInds.data(), recvSizes.data(), recvOffs.data(),
      sendInds.data(), sendSizes.data(), sendOffs.data(), comm );

    // Fulfill the requests
    vector<T> sendVals( numSendInds*width );
    const Int firstLocalRow = X.FirstLocalRow();
    for( Int s=0; s<numSendInds; ++s )
        for( Int j=0; j<width; ++j )
            sendVals[s*width+j] = X.GetLocal( sendInds[s]-firstLocalRow, j );

    // Reply with the values
    vector<T> recvVals( numRecvInds*width );
    for( int q=0; q<commSize; ++q )
    {
        sendSizes[q] *= width;
        sendOffs[q] *= width;
        recvSizes[q] *= width;
        recvOffs[q] *= width;
    }
    mpi::AllToAll
    ( sendVals.data(), sendSizes.data(), sendOffs.data(),
      recvVals.data(), recvSizes.data(), recvOffs.data(), comm );
    SwapClear( sendVals );
    SwapClear( sendSizes );
    SwapClear( sendOffs );

    // Unpack the values
    off = 0;
    offs = recvOffs;
    function<void(const SymmNodeInfo&,MatrixNode<T>&)> localUnpack =
      [&]( const SymmNodeInfo& node, MatrixNode<T>& XNode )
      {
        const Int numChildren = node.children.size();
        for( Int c=0; c<numChildren; ++c )
            localUnpack( *node.children[c], *XNode.children[c] );

        for( Int t=0; t<node.size; ++t )
        {
            const Int i = mappedInds[off++];
            const int q = X.RowOwner(i);
            for( Int j=0; j<width; ++j )
                XNode.matrix.Set( t, j, recvVals[offs[q]++] );
        }
      };
    function<void(const DistSymmNodeInfo&,DistMultiVecNode<T>&)> unpack =
      [&]( const DistSymmNodeInfo& node, DistMultiVecNode<T>& XNode )
      {
        if( node.child == nullptr )
        {
            localUnpack( *node.duplicate, *XNode.duplicate );
            return;
        }
        unpack( *node.child, *XNode.child );

        const Int localHeight = XNode.matrix.LocalHeight();
        for( Int tLoc=0; tLoc<localHeight; ++tLoc )
        {
            const Int i = mappedInds[off++];
            const int q = X.RowOwner(i);
            for( Int j=0; j<width; ++j )
                XNode.matrix.SetLocal( tLoc, j, recvVals[offs[q]++] );
        }
      };
    unpack( info, *this );
    DEBUG_ONLY(
      if( off != numRecvInds )
          LogicError("Unpacked wrong number of indices");
    )
}

template<typename T>
Int DistMultiVecNode<T>::LocalHeight() const
{
    DEBUG_ONLY(CallStackEntry cse("DistMultiVecNode::LocalHeight"))
    Int localHeight = 0;
    function<void(const DistMultiVecNode<T>&)> count = 
      [&]( const DistMultiVecNode<T>& node )
      {
        if( node.child == nullptr )
        {
            localHeight += node.duplicate->Height();
            return;
        }
        count( *node.child );
        localHeight += node.matrix.LocalHeight();
      };
    count( *this );
    return localHeight;
}

template<typename T>
void DistMultiVecNode<T>::Push
( const DistMap& invMap, const DistSymmNodeInfo& info,
        DistMultiVec<T>& X ) const
{
    DEBUG_ONLY(CallStackEntry cse("DistMultiVecNode::Push"))

    mpi::Comm comm = info.comm;
    const Int height = info.size + info.off;
    const Int width = matrix.Width();
    X.SetComm( comm );
    X.Resize( height, width );

    // Fill the set of indices that we need to map to the original ordering
    const int numSendInds = LocalHeight();
    vector<Int> mappedInds( numSendInds );
    {
        Int off=0;

        function<void(const SymmNodeInfo&,MatrixNode<T>&)> localPack = 
          [&]( const SymmNodeInfo& node, MatrixNode<T>& XNode )
          {
            const Int numChildren = node.children.size();
            for( Int c=0; c<numChildren; ++c )
                localPack( *node.children[c], *XNode.children[c] );

            for( Int t=0; t<node.size; ++t )
                mappedInds[off++] = node.off+t;
          };
        function<void(const DistSymmNodeInfo&,const DistMultiVecNode<T>&)> 
          pack = 
          [&]( const DistSymmNodeInfo& node, const DistMultiVecNode<T>& XNode )
          {
            if( node.child == nullptr )
            {
                localPack( *node.duplicate, *XNode.duplicate );
                return;
            }
            pack( *node.child, *XNode.child );

            const Int nodeHeight = XNode.matrix.Height();
            const Int colShift = XNode.matrix.ColShift();
            const Int colStride = XNode.matrix.ColStride();
            for( Int t=colShift; t<nodeHeight; t+=colStride )
                mappedInds[off++] = node.off+t;
          };

        pack( info, *this );
    }

    // Convert the indices to the original ordering
    invMap.Translate( mappedInds );

    // Figure out how many indices each process owns that we need to send
    const int commSize = mpi::Size( comm );
    vector<int> sendSizes(commSize,0);
    for( Int s=0; s<numSendInds; ++s )
        ++sendSizes[ X.RowOwner(mappedInds[s]) ];
    vector<int> sendOffs;
    Scan( sendSizes, sendOffs );

    // Pack the send indices and values
    vector<T> sendVals( numSendInds*width );
    vector<Int> sendInds( numSendInds );
    {
        Int off=0;
        auto offs = sendOffs;

        function<void(const SymmNodeInfo&,const MatrixNode<T>&)> localPack = 
          [&]( const SymmNodeInfo& node, const MatrixNode<T>& XNode )
          {
            const Int numChildren = node.children.size();
            for( Int c=0; c<numChildren; ++c )
                localPack( *node.children[c], *XNode.children[c] );

            for( Int t=0; t<node.size; ++t )
            {
                const Int i = mappedInds[off++];
                const int q = X.RowOwner(i);
                for( Int j=0; j<width; ++j )
                    sendVals[offs[q]*width+j] = XNode.matrix.Get(t,j);    
                sendInds[offs[q]++] = i;
            }
          };
        function<void(const DistSymmNodeInfo&,const DistMultiVecNode<T>&)> 
          pack =
          [&]( const DistSymmNodeInfo& node, const DistMultiVecNode<T>& XNode )
          {
            if( node.child == nullptr )
            {
                localPack( *node.duplicate, *XNode.duplicate );
                return;
            }
            pack( *node.child, *XNode.child );
            
            const Int localHeight = XNode.matrix.LocalHeight();
            for( Int tLoc=0; tLoc<localHeight; ++tLoc )
            {
                const Int i = mappedInds[off++];
                const int q = X.RowOwner(i);
                for( Int j=0; j<width; ++j )
                    sendVals[offs[q]*width+j] = XNode.matrix.GetLocal(tLoc,j);
                sendInds[offs[q]++] = i;
            }
          };
        pack( info, *this );
    }

    // Coordinate for the coming AllToAll to exchange the indices of x
    vector<int> recvSizes(commSize);
    mpi::AllToAll( sendSizes.data(), 1, recvSizes.data(), 1, comm );
    vector<int> recvOffs;
    const int numRecvInds = Scan( recvSizes, recvOffs );
    DEBUG_ONLY(
      if( numRecvInds != X.LocalHeight() )
          LogicError("numRecvInds was not equal to local height");
    )

    // Send the indices
    vector<Int> recvInds( numRecvInds );
    mpi::AllToAll
    ( sendInds.data(), sendSizes.data(), sendOffs.data(),
      recvInds.data(), recvSizes.data(), recvOffs.data(), comm );

    // Send the values
    vector<T> recvVals( numRecvInds*width );
    for( int q=0; q<commSize; ++q )
    {
        sendSizes[q] *= width;
        sendOffs[q] *= width;
        recvSizes[q] *= width;
        recvOffs[q] *= width;
    }
    mpi::AllToAll
    ( sendVals.data(), sendSizes.data(), sendOffs.data(),
      recvVals.data(), recvSizes.data(), recvOffs.data(), comm );
    SwapClear( sendVals );
    SwapClear( sendSizes );
    SwapClear( sendOffs );

    // Unpack the values
    const Int firstLocalRow = X.FirstLocalRow();
    for( Int s=0; s<numRecvInds; ++s )
    {
        const Int i = recvInds[s];
        const Int iLoc = i - firstLocalRow;
        for( Int j=0; j<width; ++j )
            X.SetLocal( iLoc, j, recvVals[s*width+j] );
    }
}

template<typename T>
void DistMultiVecNode<T>::ComputeCommMeta( const DistSymmNodeInfo& info ) const
{
    DEBUG_ONLY(CallStackEntry cse("DistMultiVecNode::ComputeCommMeta"))
    if( commMeta.numChildSendInds.size() != 0 )
        return;

    commMeta.Empty();
    if( child == nullptr )
    {
        commMeta.localOff = info.duplicate->myOff;
        commMeta.localSize = info.duplicate->size;
        return;
    }

    // This is currently assumed (and will eventually be lifted)
    const Int numChildren = 2;
    vector<int> teamSizes(numChildren), teamOffs(numChildren);

    const int teamSize = mpi::Size( info.comm );
    const int teamRank = mpi::Rank( info.comm );
    const auto& childNode = *info.child;
    const int childTeamSize = mpi::Size( childNode.comm );
    const int childTeamRank = mpi::Rank( childNode.comm );
    const bool inFirstTeam = ( childTeamRank == teamRank );
    const bool leftIsFirst = ( childNode.onLeft==inFirstTeam );
    teamSizes[0] =
        ( childNode.onLeft ? childTeamSize : teamSize-childTeamSize );
    teamSizes[1] = teamSize - teamSizes[0];
    teamOffs[0] = ( leftIsFirst ? 0 : teamSizes[1] );
    teamOffs[1] = ( leftIsFirst ? teamSizes[0] : 0 );

    const auto& myRelInds =
        ( childNode.onLeft ? info.childRelInds[0] : info.childRelInds[1] );

    // Fill numChildSendInds
    commMeta.numChildSendInds.resize( teamSize );
    El::MemZero( commMeta.numChildSendInds.data(), teamSize );
    const Int updateSize = childNode.lowerStruct.size();
    // TODO: Use 'matrix' or 'work' instead?
    {
        const Int align = childNode.size % childTeamSize;
        const Int shift = Shift( childTeamRank, align, childTeamSize );
        const Int localHeight = Length( updateSize, shift, childTeamSize );
        for( Int iChildLoc=0; iChildLoc<localHeight; ++iChildLoc )
        {
            const Int iChild = shift + iChildLoc*childTeamSize;
            const int destRank = myRelInds[iChild] % teamSize;
            ++commMeta.numChildSendInds[destRank];
        }
    }

    // Compute the solve recv indices
    commMeta.childRecvInds.resize( teamSize );
    for( Int c=0; c<numChildren; ++c )
    {
        const Int numInds = info.childRelInds[c].size();
        vector<Int> inds;
        for( Int i=0; i<numInds; ++i )
            if( work.IsLocalRow( info.childRelInds[c][i] ) )
                inds.push_back( i );

        const Int numSolveInds = inds.size();
        for( Int iPre=0; iPre<numSolveInds; ++iPre )
        {
            const Int iChild = inds[iPre];
            const Int i = info.childRelInds[c][iChild];
            const Int iLoc = (i-teamRank) / teamSize;
            const int childRank = (info.childSizes[c]+iChild) % teamSizes[c];
            const int frontRank = teamOffs[c] + childRank;
            commMeta.childRecvInds[frontRank].push_back(iLoc);
        }
    }
    commMeta.localOff = child->commMeta.localOff + child->commMeta.localSize;
    commMeta.localSize = matrix.LocalHeight();
}

#define PROTO(T) template class DistMultiVecNode<T>;
#include "El/macros/Instantiate.h"

} // namespace El
