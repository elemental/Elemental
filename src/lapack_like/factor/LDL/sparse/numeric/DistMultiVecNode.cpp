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
namespace ldl {

template<typename T>
DistMultiVecNode<T>::DistMultiVecNode( DistMultiVecNode<T>* parentNode )
: parent(parentNode), child(nullptr), duplicate(nullptr)
{ }

template<typename T>
DistMultiVecNode<T>::DistMultiVecNode
( const DistMap& invMap,
  const DistNodeInfo& info,
  const DistMultiVec<T>& X )
: parent(nullptr), child(nullptr), duplicate(nullptr)
{
    DEBUG_ONLY(CSE cse("DistMultiVecNode::DistMultiVecNode"))
    Pull( invMap, info, X );
}

template<typename T>
DistMultiVecNode<T>::DistMultiVecNode( const DistMatrixNode<T>& X )
: parent(nullptr), child(nullptr), duplicate(nullptr)
{
    DEBUG_ONLY(CSE cse("DistMultiVecNode::DistMultiVecNode"))
    *this = X;
}

template<typename T>
DistMultiVecNode<T>::~DistMultiVecNode()
{
    delete child;
    delete duplicate;
}

template<typename T>
const DistMultiVecNode<T>&
DistMultiVecNode<T>::operator=( const DistMatrixNode<T>& X )
{
    DEBUG_ONLY(CSE cse("DistMultiVecNode::operator="))

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

// TODO: Consider passing in a pre-built origOwners
template<typename T>
void DistMultiVecNode<T>::Pull
( const DistMap& invMap,
  const DistNodeInfo& info,
  const DistMultiVec<T>& X )
{
    DEBUG_ONLY(CSE cse("DistMultiVecNode::Pull"))
    DistMultiVecNodeMeta meta;
    Pull( invMap, info, X, meta );
}

template<typename T>
void DistMultiVecNode<T>::Pull
( const DistMap& invMap,
  const DistNodeInfo& info,
  const DistMultiVec<T>& X,
        DistMultiVecNodeMeta& meta )
{
    DEBUG_ONLY(CSE cse("DistMultiVecNode::Pull"))
    const Int width = X.Width();
    const T* XBuf = X.LockedMatrix().LockedBuffer();
    const Int XLDim = X.LockedMatrix().LDim();
    const Int firstLocalRow = X.FirstLocalRow();
    mpi::Comm comm = X.Comm();
    const int commSize = mpi::Size( comm );

    function<void(const NodeInfo&,MatrixNode<T>&)> localInit =
      [&]( const NodeInfo& node, MatrixNode<T>& XNode )
      {
        const Int numChildren = node.children.size();
        if( XNode.children.size() != node.children.size() )
        {
            for( auto* childNode : XNode.children )
                delete childNode;
            XNode.children.resize(numChildren);
            for( Int c=0; c<numChildren; ++c )
                XNode.children[c] = new MatrixNode<T>(&XNode);
        }

        for( Int c=0; c<numChildren; ++c )
            localInit( *node.children[c], *XNode.children[c] );

        XNode.matrix.Resize( node.size, width );
      };
    function<void(const DistNodeInfo&,DistMultiVecNode<T>&)> init =
      [&]( const DistNodeInfo& node, DistMultiVecNode<T>& XNode )
      {
        if( node.child == nullptr )
        {
            DEBUG_ONLY(
              if( XNode.child != nullptr )
                  LogicError("Child should have been a nullptr");
            )
            if( XNode.duplicate == nullptr )
                XNode.duplicate = new MatrixNode<T>(&XNode);
            localInit( *node.duplicate, *XNode.duplicate );

            XNode.matrix.Attach( *node.grid, XNode.duplicate->matrix );
            return;
        }

        DEBUG_ONLY(
          if( XNode.duplicate != nullptr )
              LogicError("Duplicate should have been a nullptr");
        )
        if( XNode.child == nullptr )
            XNode.child = new DistMultiVecNode<T>(&XNode);
        init( *node.child, *XNode.child );

        if( XNode.matrix.Grid() != *node.grid ||
            XNode.matrix.Height() != node.size ||
            XNode.matrix.Width() != width )
        {
            XNode.commMeta.Empty();
            XNode.matrix.SetGrid( *node.grid );
            XNode.matrix.Resize( node.size, width );
        }
      };
    init( info, *this );

    const Int numRecvInds = LocalHeight();
    if( meta.sendInds.size() == 0 || meta.recvInds.size() == 0 ||
        meta.mappedOwners.size() == 0 ||
        meta.sendSizes.size() == 0 || meta.sendOffs.size() == 0 ||
        meta.recvSizes.size() == 0 || meta.recvOffs.size() == 0 )
    {
        // Pack the indices for mapping to the original ordering
        Int off = 0;
        vector<Int> mappedInds( numRecvInds );
        function<void(const NodeInfo&,MatrixNode<T>&)> localPack =
          [&]( const NodeInfo& node, MatrixNode<T>& XNode )
          {
            const Int numChildren = node.children.size();
            for( Int c=0; c<numChildren; ++c )
                localPack( *node.children[c], *XNode.children[c] );

            for( Int t=0; t<node.size; ++t )
                mappedInds[off++] = node.off+t;
          };
        function<void(const DistNodeInfo&,DistMultiVecNode<T>&)> pack =
          [&]( const DistNodeInfo& node, DistMultiVecNode<T>& XNode )
          {
            if( node.child == nullptr )
            {
                localPack( *node.duplicate, *XNode.duplicate );
                return;
            }
            pack( *node.child, *XNode.child );

            const Int XNodeLocalHeight = XNode.matrix.LocalHeight();
            for( Int tLoc=0; tLoc<XNodeLocalHeight; ++tLoc )
            {
                const Int t = XNode.matrix.GlobalRow(tLoc);
                mappedInds[off++] = node.off+t;
            }
          };
        pack( info, *this );

        // Convert the indices to the original ordering
        // TODO: Consider passing in a pre-built origOwners
        invMap.Translate( mappedInds );

        meta.mappedOwners.resize( numRecvInds );
        for( int s=0; s<numRecvInds; ++s )
            meta.mappedOwners[s] = X.RowOwner(mappedInds[s]);

        // Figure out how many entries each process owns that we need
        meta.recvSizes.resize( commSize ); 
        for( int q=0; q<commSize; ++q )
            meta.recvSizes[q] = 0;
        for( int s=0; s<numRecvInds; ++s )
            ++meta.recvSizes[meta.mappedOwners[s]];
        Scan( meta.recvSizes, meta.recvOffs );
        meta.recvInds.resize( numRecvInds );
        auto offs = meta.recvOffs;
        for( int s=0; s<numRecvInds; ++s )
            meta.recvInds[offs[meta.mappedOwners[s]]++] = mappedInds[s];

        // Coordinate for the coming AllToAll to exchange the indices of X
        meta.sendSizes.resize( commSize );
        mpi::AllToAll
        ( meta.recvSizes.data(), 1, meta.sendSizes.data(), 1, comm );
        const int numSendInds = Scan( meta.sendSizes, meta.sendOffs );

        // Request the indices
        meta.sendInds.resize( numSendInds );
        mpi::AllToAll
        ( meta.recvInds.data(), meta.recvSizes.data(), meta.recvOffs.data(),
          meta.sendInds.data(), meta.sendSizes.data(), meta.sendOffs.data(),
          comm );
    }

    // Fulfill the requests
    const Int numSendInds = meta.sendInds.size();
    vector<T> sendVals( numSendInds*width );
    for( Int s=0; s<numSendInds; ++s )
        for( Int j=0; j<width; ++j )
            sendVals[s*width+j] = XBuf[meta.sendInds[s]-firstLocalRow+j*XLDim];

    // Reply with the values
    vector<T> recvVals( numRecvInds*width );
    for( int q=0; q<commSize; ++q )
    {
        meta.sendSizes[q] *= width;
        meta.sendOffs[q] *= width;
        meta.recvSizes[q] *= width;
        meta.recvOffs[q] *= width;
    }
    mpi::AllToAll
    ( sendVals.data(), meta.sendSizes.data(), meta.sendOffs.data(),
      recvVals.data(), meta.recvSizes.data(), meta.recvOffs.data(), comm );
    SwapClear( sendVals );

    // Unpack the values
    Int off = 0;
    auto offs = meta.recvOffs;
    function<void(const NodeInfo&,MatrixNode<T>&)> localUnpack =
      [&]( const NodeInfo& node, MatrixNode<T>& XNode )
      {
        const Int numChildren = node.children.size();
        for( Int c=0; c<numChildren; ++c )
            localUnpack( *node.children[c], *XNode.children[c] );

        T* XNodeBuf = XNode.matrix.Buffer();
        const Int XNodeLDim = XNode.matrix.LDim();
        for( Int t=0; t<node.size; ++t )
        {
            const int q = meta.mappedOwners[off++];
            for( Int j=0; j<width; ++j )
                XNodeBuf[t+j*XNodeLDim] = recvVals[offs[q]++];
        }
      };
    function<void(const DistNodeInfo&,DistMultiVecNode<T>&)> unpack =
      [&]( const DistNodeInfo& node, DistMultiVecNode<T>& XNode )
      {
        if( node.child == nullptr )
        {
            localUnpack( *node.duplicate, *XNode.duplicate );
            return;
        }
        unpack( *node.child, *XNode.child );

        const Int localHeight = XNode.matrix.LocalHeight();
        T* XNodeBuf = XNode.matrix.Buffer();
        const Int XNodeLDim = XNode.matrix.LDim();
        for( Int tLoc=0; tLoc<localHeight; ++tLoc )
        {
            const int q = meta.mappedOwners[off++];
            for( Int j=0; j<width; ++j )
                XNodeBuf[tLoc+j*XNodeLDim] = recvVals[offs[q]++];
        }
      };
    unpack( info, *this );
    DEBUG_ONLY(
      if( off != numRecvInds )
          LogicError("Unpacked wrong number of indices");
    )

    // Exit with the sizes/offsets calibrated to a single right-hand side
    for( int q=0; q<commSize; ++q )
    {
        meta.sendSizes[q] /= width;
        meta.sendOffs[q] /= width;
        meta.recvSizes[q] /= width;
        meta.recvOffs[q] /= width;
    }
}

template<typename T>
void DistMultiVecNode<T>::Push
( const DistMap& invMap,
  const DistNodeInfo& info,
        DistMultiVec<T>& X ) const
{
    DEBUG_ONLY(CSE cse("DistMultiVecNode::Push"))
    DistMultiVecNodeMeta meta;
    Push( invMap, info, X, meta );
}

template<typename T>
void DistMultiVecNode<T>::Push
( const DistMap& invMap,
  const DistNodeInfo& info,
        DistMultiVec<T>& X,
        DistMultiVecNodeMeta& meta ) const
{
    DEBUG_ONLY(CSE cse("DistMultiVecNode::Push"))
    const Int height = info.size + info.off;
    const Int width = matrix.Width();

    mpi::Comm comm = info.comm;
    const int commSize = mpi::Size( comm );
    X.SetComm( comm );
    X.Resize( height, width );
    const Int firstLocalRow = X.FirstLocalRow();
    T* XBuf = X.Matrix().Buffer();
    const Int XLDim = X.Matrix().LDim();

    // Fill the set of indices that we need to map to the original ordering
    const int numSendInds = LocalHeight();

    if( meta.sendInds.size() == 0 || meta.recvInds.size() == 0 ||
        meta.mappedOwners.size() == 0 ||
        meta.sendSizes.size() == 0 || meta.sendOffs.size() == 0 ||
        meta.recvSizes.size() == 0 || meta.recvOffs.size() == 0 )
    {
        std::vector<Int> mappedInds( numSendInds );
        {
            Int off=0;

            function<void(const NodeInfo&,MatrixNode<T>&)> localPack = 
              [&]( const NodeInfo& node, MatrixNode<T>& XNode )
              {
                const Int numChildren = node.children.size();
                for( Int c=0; c<numChildren; ++c )
                    localPack( *node.children[c], *XNode.children[c] );

                for( Int t=0; t<node.size; ++t )
                    mappedInds[off++] = node.off+t;
              };
            function<void(const DistNodeInfo&,const DistMultiVecNode<T>&)> 
              pack = 
              [&]( const DistNodeInfo& node, const DistMultiVecNode<T>& XNode )
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
        // TODO: Consider passing in a pre-built origOwners
        invMap.Translate( mappedInds );

        meta.mappedOwners.resize( numSendInds );
        for( Int s=0; s<numSendInds; ++s )
            meta.mappedOwners[s] = X.RowOwner(mappedInds[s]);

        // Figure out how many indices each process owns that we need to send
        meta.sendSizes.resize( commSize );
        for( int q=0; q<commSize; ++q )
            meta.sendSizes[q] = 0;
        for( Int s=0; s<numSendInds; ++s )
            ++meta.sendSizes[meta.mappedOwners[s]];
        Scan( meta.sendSizes, meta.sendOffs );

        // Pack the send indices
        meta.sendInds.resize( numSendInds );
        {
            Int off=0;
            auto offs = meta.sendOffs;

            function<void(const NodeInfo&,const MatrixNode<T>&)> localPackInds =
              [&]( const NodeInfo& node, const MatrixNode<T>& XNode )
              {
                const Int numChildren = node.children.size();
                for( Int c=0; c<numChildren; ++c )
                    localPackInds( *node.children[c], *XNode.children[c] );

                const T* XNodeBuf = XNode.matrix.LockedBuffer();
                const Int XNodeLDim = XNode.matrix.LDim();
                for( Int t=0; t<node.size; ++t )
                {
                    const Int i = mappedInds[off];
                    const int q = meta.mappedOwners[off++];
                    meta.sendInds[offs[q]++] = i;
                }
              };
            function<void(const DistNodeInfo&,const DistMultiVecNode<T>&)> 
              packInds =
              [&]( const DistNodeInfo& node, const DistMultiVecNode<T>& XNode )
              {
                if( node.child == nullptr )
                {
                    localPackInds( *node.duplicate, *XNode.duplicate );
                    return;
                }
                packInds( *node.child, *XNode.child );
            
                const Int localHeight = XNode.matrix.LocalHeight();
                const T* XNodeBuf = XNode.matrix.LockedBuffer();
                const Int XNodeLDim = XNode.matrix.LDim();
                for( Int tLoc=0; tLoc<localHeight; ++tLoc )
                {
                    const Int i = mappedInds[off];
                    const int q = meta.mappedOwners[off++];
                    meta.sendInds[offs[q]++] = i;
                }
              };
            packInds( info, *this );
        }

        // Coordinate for the coming AllToAll to exchange the indices of x
        meta.recvSizes.resize( commSize );
        mpi::AllToAll
        ( meta.sendSizes.data(), 1, meta.recvSizes.data(), 1, comm );
        const int numRecvInds = Scan( meta.recvSizes, meta.recvOffs );
        DEBUG_ONLY(
          if( numRecvInds != X.LocalHeight() )
              LogicError("numRecvInds was not equal to local height");
        )

        // Send the indices
        meta.recvInds.resize( numRecvInds );
        mpi::AllToAll
        ( meta.sendInds.data(), meta.sendSizes.data(), meta.sendOffs.data(),
          meta.recvInds.data(), meta.recvSizes.data(), meta.recvOffs.data(),
          comm );
    }

    // Pack the send values
    const Int numRecvInds = meta.recvInds.size();
    vector<T> sendVals( numSendInds*width );
    {
        Int off=0;
        auto offs = meta.sendOffs;
        function<void(const NodeInfo&,const MatrixNode<T>&)> localPack = 
          [&]( const NodeInfo& node, const MatrixNode<T>& XNode )
          {
            const Int numChildren = node.children.size();
            for( Int c=0; c<numChildren; ++c )
                localPack( *node.children[c], *XNode.children[c] );

            const T* XNodeBuf = XNode.matrix.LockedBuffer();
            const Int XNodeLDim = XNode.matrix.LDim();
            for( Int t=0; t<node.size; ++t )
            {
                const int q = meta.mappedOwners[off++];
                for( Int j=0; j<width; ++j )
                    sendVals[offs[q]*width+j] = XNodeBuf[t+j*XNodeLDim];
                offs[q]++;
            }
          };
        function<void(const DistNodeInfo&,const DistMultiVecNode<T>&)> 
          pack =
          [&]( const DistNodeInfo& node, const DistMultiVecNode<T>& XNode )
          {
            if( node.child == nullptr )
            {
                localPack( *node.duplicate, *XNode.duplicate );
                return;
            }
            pack( *node.child, *XNode.child );
            
            const Int localHeight = XNode.matrix.LocalHeight();
            const T* XNodeBuf = XNode.matrix.LockedBuffer();
            const Int XNodeLDim = XNode.matrix.LDim();
            for( Int tLoc=0; tLoc<localHeight; ++tLoc )
            {
                const int q = meta.mappedOwners[off++];
                for( Int j=0; j<width; ++j )
                    sendVals[offs[q]*width+j] = XNodeBuf[tLoc+j*XNodeLDim];
                offs[q]++;
            }
          };
        pack( info, *this );
    }

    // Send the values
    vector<T> recvVals( numRecvInds*width );
    for( int q=0; q<commSize; ++q )
    {
        meta.sendSizes[q] *= width;
        meta.sendOffs[q] *= width;
        meta.recvSizes[q] *= width;
        meta.recvOffs[q] *= width;
    }
    mpi::AllToAll
    ( sendVals.data(), meta.sendSizes.data(), meta.sendOffs.data(),
      recvVals.data(), meta.recvSizes.data(), meta.recvOffs.data(), comm );
    SwapClear( sendVals );

    // Unpack the values
    for( Int s=0; s<numRecvInds; ++s )
    {
        const Int i = meta.recvInds[s];
        const Int iLoc = i - firstLocalRow;
        for( Int j=0; j<width; ++j )
            XBuf[iLoc+j*XLDim] = recvVals[s*width+j];
    }

    for( int q=0; q<commSize; ++q )
    {
        meta.sendSizes[q] /= width;
        meta.sendOffs[q] /= width;
        meta.recvSizes[q] /= width;
        meta.recvOffs[q] /= width;
    }
}

template<typename T>
Int DistMultiVecNode<T>::LocalHeight() const
{
    DEBUG_ONLY(CSE cse("DistMultiVecNode::LocalHeight"))
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
void DistMultiVecNode<T>::ComputeCommMeta( const DistNodeInfo& info ) const
{
    DEBUG_ONLY(CSE cse("DistMultiVecNode::ComputeCommMeta"))
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

#define PROTO(T) template struct DistMultiVecNode<T>;
#include "El/macros/Instantiate.h"

} // namespace ldl
} // namespace El
