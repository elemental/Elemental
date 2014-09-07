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
#ifndef EL_SPARSEDIRECT_NUMERIC_LOWERMULTIPLY_DIST_HPP
#define EL_SPARSEDIRECT_NUMERIC_LOWERMULTIPLY_DIST_HPP

namespace El {

template<typename T> 
inline void DistLowerMultiplyNormal
( int diagOff, const DistSymmInfo& info, 
  const DistSymmFrontTree<T>& L, DistNodalMultiVec<T>& X )
{
    DEBUG_ONLY(CallStackEntry cse("DistLowerMultiplyNormal"))
    const int numDistNodes = info.distNodes.size();
    const int width = X.Width();
    if( L.frontType != SYMM_1D && L.frontType != LDL_1D )
        LogicError("This multiply mode is not yet implemented");

    // Copy the information from the local portion into the distributed leaf
    const SymmFront<T>& localRootFront = L.localFronts.back();
    const DistSymmFront<T>& distLeafFront = L.distFronts[0];
    distLeafFront.work1d.LockedAttach
    ( localRootFront.work.Height(), localRootFront.work.Width(), 
      distLeafFront.front1dL.Grid(), 0, 0,
      localRootFront.work.LockedBuffer(), localRootFront.work.LDim() );
    
    // Perform the distributed portion of the forward multiply
    for( int s=1; s<numDistNodes; ++s )
    {
        const DistSymmNodeInfo& childNode = info.distNodes[s-1];
        const DistSymmNodeInfo& node = info.distNodes[s];
        const DistSymmFront<T>& childFront = L.distFronts[s-1];
        const DistSymmFront<T>& front = L.distFronts[s];
        const Grid& childGrid = childFront.front1dL.Grid();
        const Grid& grid = front.front1dL.Grid();
        mpi::Comm comm = grid.VCComm();
        mpi::Comm childComm = childGrid.VCComm();
        const int commSize = mpi::Size( comm );
        const int childCommSize = mpi::Size( childComm );

        // Set up a workspace
        DistMatrix<T,VC,STAR>& W = front.work1d;
        W.SetGrid( grid );
        W.Resize( front.front1dL.Height(), width );
        DistMatrix<T,VC,STAR> WT(grid), WB(grid);
        PartitionDown( W, WT, WB, node.size );
        WT = X.distNodes[s-1];
        Zero( WB );

        // Now that the right-hand side is set up, perform the multiply
        FrontLowerMultiply( NORMAL, diagOff, front.front1dL, W );

        // Pack our child's update
        const MultiVecCommMeta& commMeta = node.multiVecMeta;
        DistMatrix<T,VC,STAR>& childW = childFront.work1d;
        const int updateSize = childW.Height()-childNode.size;
        auto childUpdate = 
            LockedView( childW, childNode.size, 0, updateSize, width );
        int sendBufferSize = 0;
        std::vector<int> sendCounts(commSize), sendDispls(commSize);
        for( int proc=0; proc<commSize; ++proc )
        {
            const int sendSize = commMeta.numChildSendInds[proc]*width;
            sendCounts[proc] = sendSize;
            sendDispls[proc] = sendBufferSize;
            sendBufferSize += sendSize;
        }
        std::vector<T> sendBuffer( sendBufferSize );

        const std::vector<int>& myChildRelInds = 
            ( childNode.onLeft ? node.leftRelInds : node.rightRelInds );
        const int colShift = childUpdate.ColShift();
        const int localHeight = childUpdate.LocalHeight();
        std::vector<int> packOffs = sendDispls;
        for( int iChildLoc=0; iChildLoc<localHeight; ++iChildLoc )
        {
            const int iChild = colShift + iChildLoc*childCommSize;
            const int destRank = myChildRelInds[iChild] % commSize;
            for( int jChild=0; jChild<width; ++jChild )
                sendBuffer[packOffs[destRank]++] = 
                    childUpdate.GetLocal(iChildLoc,jChild);
        }
        SwapClear( packOffs );
        childW.Empty();
        if( s == 1 )
            L.localFronts.back().work.Empty();

        // Set up the receive buffer
        int recvBufferSize = 0;
        std::vector<int> recvCounts(commSize), recvDispls(commSize);
        for( int proc=0; proc<commSize; ++proc )
        {
            const int recvSize = commMeta.childRecvInds[proc].size()*width;
            recvCounts[proc] = recvSize;
            recvDispls[proc] = recvBufferSize;
            recvBufferSize += recvSize;
        }
        std::vector<T> recvBuffer( recvBufferSize );
        DEBUG_ONLY(VerifySendsAndRecvs( sendCounts, recvCounts, comm ))

        // AllToAll to send and receive the child updates
        SparseAllToAll
        ( sendBuffer, sendCounts, sendDispls,
          recvBuffer, recvCounts, recvDispls, comm );
        SwapClear( sendBuffer );
        SwapClear( sendCounts );
        SwapClear( sendDispls );

        // Unpack the child updates (with an Axpy)
        for( int proc=0; proc<commSize; ++proc )
        {
            const T* recvVals = &recvBuffer[recvDispls[proc]];
            const std::vector<int>& recvInds = commMeta.childRecvInds[proc];
            const int numRecvInds = recvInds.size();
            for( int k=0; k<numRecvInds; ++k )
            {
                const int iFrontLoc = recvInds[k];
                const T* recvRow = &recvVals[k*width];
                T* WRow = W.Buffer( iFrontLoc, 0 );
                const int WLDim = W.LDim();
                for( int jFront=0; jFront<width; ++jFront )
                    WRow[jFront*WLDim] += recvRow[jFront];
            }
        }
        SwapClear( recvBuffer );
        SwapClear( recvCounts );
        SwapClear( recvDispls );

        // Store this node's portion of the result
        X.distNodes[s-1] = WT;
    }
    L.localFronts.back().work.Empty();
    L.distFronts.back().work1d.Empty();
}

template<typename T> 
inline void DistLowerMultiplyTranspose
( int diagOff, const DistSymmInfo& info, 
  const DistSymmFrontTree<T>& L, DistNodalMultiVec<T>& X, 
  bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("DistLowerMultiplyTranspose"))
    const int numDistNodes = info.distNodes.size();
    const int width = X.Width();
    const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );
    if( L.frontType != SYMM_1D && L.frontType != LDL_1D )
        LogicError("This multiply mode is not yet implemented");

    // Directly operate on the root separator's portion of the right-hand sides
    const DistSymmNodeInfo& rootNode = info.distNodes.back();
    const SymmFront<T>& localRootFront = L.localFronts.back();
    if( numDistNodes == 1 )
    {
        Matrix<T>& XRoot = X.localNodes.back();
        localRootFront.work = XRoot;
        FrontLowerMultiply
        ( orientation, diagOff, localRootFront.frontL, XRoot );
    }
    else
    {
        const DistSymmFront<T>& rootFront = L.distFronts.back();
        DistMatrix<T,VC,STAR>& XRoot = X.distNodes.back();
        rootFront.work1d = XRoot;
        FrontLowerMultiply( orientation, diagOff, rootFront.front1dL, XRoot );
    }

    for( int s=numDistNodes-2; s>=0; --s )
    {
        const DistSymmNodeInfo& parentNode = info.distNodes[s+1];
        const DistSymmNodeInfo& node = info.distNodes[s];
        const DistSymmFront<T>& parentFront = L.distFronts[s+1];
        const DistSymmFront<T>& front = L.distFronts[s];
        const Grid& grid = front.front1dL.Grid();
        const Grid& parentGrid = parentFront.front1dL.Grid();
        mpi::Comm comm = grid.VCComm(); 
        mpi::Comm parentComm = parentGrid.VCComm();
        const int commSize = mpi::Size( comm );
        const int parentCommSize = mpi::Size( parentComm );

        // Set up a copy of the RHS in our workspace.
        DistMatrix<T,VC,STAR>& W = front.work1d;
        W.SetGrid( grid );
        W.Resize( front.front1dL.Height(), width );
        DistMatrix<T,VC,STAR> WT(grid), WB(grid);
        PartitionDown( W, WT, WB, node.size );
        Matrix<T>& localXT = 
          ( s>0 ? X.distNodes[s-1].Matrix() : X.localNodes.back() );
        WT.Matrix() = localXT;

        //
        // Set the bottom from the parent's workspace
        //

        // Pack the relevant portions of the parent's RHS's
        // (which are stored in 'work1d')
        const MultiVecCommMeta& commMeta = parentNode.multiVecMeta;
        int sendBufferSize = 0;
        std::vector<int> sendCounts(parentCommSize), sendDispls(parentCommSize);
        for( int proc=0; proc<parentCommSize; ++proc )
        {
            const int sendSize = commMeta.childRecvInds[proc].size()*width;
            sendCounts[proc] = sendSize;
            sendDispls[proc] = sendBufferSize;
            sendBufferSize += sendSize;
        }
        std::vector<T> sendBuffer( sendBufferSize );

        DistMatrix<T,VC,STAR>& parentWork = parentFront.work1d;
        for( int proc=0; proc<parentCommSize; ++proc )
        {
            T* sendVals = &sendBuffer[sendDispls[proc]];
            const std::vector<int>& recvInds = commMeta.childRecvInds[proc];
            const int numRecvInds = recvInds.size();
            for( int k=0; k<numRecvInds; ++k )
            {
                const int iFrontLoc = recvInds[k];
                T* packedRow = &sendVals[k*width];
                const T* workRow = parentWork.LockedBuffer( iFrontLoc, 0 );
                const int workLDim = parentWork.LDim();
                for( int jFront=0; jFront<width; ++jFront )
                    packedRow[jFront] = workRow[jFront*workLDim];
            }
        }
        parentWork.Empty();

        // Set up the receive buffer
        int recvBufferSize = 0;
        std::vector<int> recvCounts(parentCommSize), recvDispls(parentCommSize);
        for( int proc=0; proc<parentCommSize; ++proc )
        {
            const int recvSize = commMeta.numChildSendInds[proc]*width;
            recvCounts[proc] = recvSize;
            recvDispls[proc] = recvBufferSize;
            recvBufferSize += recvSize;
        }
        std::vector<T> recvBuffer( recvBufferSize );
        DEBUG_ONLY(VerifySendsAndRecvs( sendCounts, recvCounts, parentComm ))

        // AllToAll to send and recv parent updates
        SparseAllToAll
        ( sendBuffer, sendCounts, sendDispls,
          recvBuffer, recvCounts, recvDispls, parentComm );
        SwapClear( sendBuffer );
        SwapClear( sendCounts );
        SwapClear( sendDispls );

        // Unpack the updates using the send approach from the forward solve
        const std::vector<int>& myRelInds = 
            ( node.onLeft ? parentNode.leftRelInds : parentNode.rightRelInds );
        const int colShift = WB.ColShift();
        const int localHeight = WB.LocalHeight();
        for( int iUpdateLoc=0; iUpdateLoc<localHeight; ++iUpdateLoc )
        {
            const int iUpdate = colShift + iUpdateLoc*commSize;
            const int startRank = myRelInds[iUpdate] % parentCommSize;
            const T* recvBuf = &recvBuffer[recvDispls[startRank]];
            for( int jUpdate=0; jUpdate<width; ++jUpdate )
                WB.SetLocal(iUpdateLoc,jUpdate,recvBuf[jUpdate]);
            recvDispls[startRank] += width;
        }
        SwapClear( recvBuffer );
        SwapClear( recvCounts );
        SwapClear( recvDispls );

        // Make a copy of the unmodified RHS
        DistMatrix<T,VC,STAR> XNode( W );

        // Perform the multiply for this front
        if( s > 0 )
            FrontLowerMultiply( orientation, diagOff, front.front1dL, XNode );
        else
        {
            localRootFront.work = W.Matrix();
            FrontLowerMultiply
            ( orientation, diagOff, localRootFront.frontL, XNode.Matrix() );
        }

        // Store the supernode portion of the result
        DistMatrix<T,VC,STAR> XNodeT(grid), XNodeB(grid);
        PartitionDown( XNode, XNodeT, XNodeB, node.size );
        localXT = XNodeT.Matrix();
        XNode.Empty();
    }
}

} // namespace El

#endif // ifndef EL_SPARSEDIRECT_NUMERIC_LOWERMULTIPLY_DIST_HPP
