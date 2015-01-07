/*
   Copyright (c) 2009-2015, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, Stanford University, and the
   Georgia Insitute of Technology.
   All rights reserved.
 
   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_SPARSEDIRECT_NUMERIC_LOWERSOLVE_DIST_HPP
#define EL_SPARSEDIRECT_NUMERIC_LOWERSOLVE_DIST_HPP

namespace El {

template<typename F> 
inline void DistLowerForwardSolve
( const DistSymmInfo& info, 
  const DistSymmFrontTree<F>& L, DistNodalMultiVec<F>& X )
{
    DEBUG_ONLY(CallStackEntry cse("DistLowerForwardSolve"))
    const Int numDistNodes = info.distNodes.size();
    const Int width = X.Width();
    const SymmFrontType frontType = L.frontType;
    const bool frontsAre1d = FrontsAre1d( frontType );
    if( Unfactored(frontType) )
        LogicError("Nonsensical front type for solve");
    if( frontType == LDL_2D || frontType == LDL_INTRAPIV_2D )
        LogicError("2D non-inverted solves supported by this routine");
    if( frontType == BLOCK_LDL_1D || frontType == BLOCK_LDL_INTRAPIV_1D )
        LogicError("1D block solves not yet supported by this routine");
    const bool blocked = BlockFactorization( frontType );

    // Copy the information from the local portion into the distributed leaf
    const SymmFront<F>& localRootFront = L.localFronts.back();
    const DistSymmFront<F>& distLeafFront = L.distFronts[0];
    const Grid& leafGrid = ( frontsAre1d ? distLeafFront.front1dL.Grid() 
                                         : distLeafFront.front2dL.Grid() );
    distLeafFront.work1d.LockedAttach
    ( localRootFront.work.Height(), localRootFront.work.Width(), leafGrid, 0, 0,
      localRootFront.work );
    
    // Perform the distributed portion of the forward solve
    for( Int s=1; s<numDistNodes; ++s )
    {
        const DistSymmNodeInfo& node = info.distNodes[s];
        const DistSymmFront<F>& front = L.distFronts[s];
        const Grid& grid = ( frontsAre1d ? front.front1dL.Grid()
                                         : front.front2dL.Grid() );
        mpi::Comm comm = grid.VCComm();
        const int commSize = mpi::Size( comm );

        const DistSymmNodeInfo& childNode = info.distNodes[s-1];
        const DistSymmFront<F>& childFront = L.distFronts[s-1];
        const Grid& childGrid = ( frontsAre1d ? childFront.front1dL.Grid()
                                              : childFront.front2dL.Grid() );
        mpi::Comm childComm = childGrid.VCComm();
        const int childCommSize = mpi::Size( childComm );

        // Set up a workspace
        const Int frontHeight = ( frontsAre1d ? front.front1dL.Height()
                                              : front.front2dL.Height() );
        DistMatrix<F,VC,STAR>& W = front.work1d;
        W.SetGrid( grid );
        W.Resize( frontHeight, width );
        DistMatrix<F,VC,STAR> WT(grid), WB(grid);
        PartitionDown( W, WT, WB, node.size );
        WT = X.distNodes[s-1];
        Zero( WB );

        // Pack our child's update
        const MultiVecCommMeta& commMeta = node.multiVecMeta;
        DistMatrix<F,VC,STAR>& childW = childFront.work1d;
        const Int updateSize = childW.Height()-childNode.size;
        DistMatrix<F,VC,STAR> childUpdate( childW.Grid() );
        LockedView( childUpdate, childW, childNode.size, 0, updateSize, width );
        int sendBufferSize = 0;
        std::vector<int> sendCounts(commSize), sendDispls(commSize);
        for( int proc=0; proc<commSize; ++proc )
        {
            const int sendSize = commMeta.numChildSendInds[proc]*width;
            sendCounts[proc] = sendSize;
            sendDispls[proc] = sendBufferSize;
            sendBufferSize += sendSize;
        }
        std::vector<F> sendBuffer( sendBufferSize );

        const bool onLeft = childNode.onLeft;
        const std::vector<Int>& myChildRelInds = 
            ( onLeft ? node.leftRelInds : node.rightRelInds );
        const int colShift = childUpdate.ColShift();
        const Int localHeight = childUpdate.LocalHeight();
        auto packOffs = sendDispls;
        for( Int iChildLoc=0; iChildLoc<localHeight; ++iChildLoc )
        {
            const Int iChild = colShift + iChildLoc*childCommSize;
            const int destRank = myChildRelInds[iChild] % commSize;
            for( Int jChild=0; jChild<width; ++jChild )
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
        std::vector<F> recvBuffer( recvBufferSize );
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
            const F* recvVals = &recvBuffer[recvDispls[proc]];
            const std::vector<Int>& recvInds = commMeta.childRecvInds[proc];
            for( unsigned k=0; k<recvInds.size(); ++k )
                blas::Axpy
                ( width, F(1), 
                  &recvVals[k*width], 1, W.Buffer(recvInds[k],0), W.LDim() );
        }
        SwapClear( recvBuffer );
        SwapClear( recvCounts );
        SwapClear( recvDispls );

        // Now that the RHS is set up, perform this node's solve
        if( frontType == LDL_1D )
            FrontLowerForwardSolve( front.front1dL, W );
        else if( frontType == LDL_SELINV_1D )
            FrontFastLowerForwardSolve( front.front1dL, W );
        else if( frontType == LDL_SELINV_2D )
            FrontFastLowerForwardSolve( front.front2dL, W );
        else if( frontType == LDL_INTRAPIV_1D )
            FrontIntraPivLowerForwardSolve( front.front1dL, front.piv, W );
        else if( frontType == LDL_INTRAPIV_SELINV_1D )
            FrontFastIntraPivLowerForwardSolve( front.front1dL, front.piv, W );
        else if( frontType == LDL_INTRAPIV_SELINV_2D )
            FrontFastIntraPivLowerForwardSolve( front.front2dL, front.piv, W );
        else if( blocked )
            FrontBlockLowerForwardSolve( front.front2dL, W );
        else
            LogicError("Unsupported front type");

        // Store this node's portion of the result
        X.distNodes[s-1] = WT;
    }
    L.localFronts.back().work.Empty();
    L.distFronts.back().work1d.Empty();
}

template<typename F> 
inline void DistLowerForwardSolve
( const DistSymmInfo& info, 
  const DistSymmFrontTree<F>& L, DistNodalMatrix<F>& X )
{
    DEBUG_ONLY(CallStackEntry cse("DistLowerForwardSolve"))
    const Int numDistNodes = info.distNodes.size();
    const Int width = X.Width();
    const SymmFrontType frontType = L.frontType;
    if( Unfactored(frontType) )
        LogicError("Nonsensical front type for solve");
    if( FrontsAre1d(frontType) )
        LogicError("1d solves not yet implemented");
    const bool blocked = BlockFactorization( frontType );

    const bool computeCommMetas = ( X.commMetas.size() == 0 );
    if( computeCommMetas )
        X.ComputeCommMetas( info );

    // Copy the information from the local portion into the distributed leaf
    const SymmFront<F>& localRootFront = L.localFronts.back();
    const DistSymmFront<F>& distLeafFront = L.distFronts[0];
    const Grid& leafGrid = distLeafFront.front2dL.Grid();
    distLeafFront.work2d.LockedAttach
    ( localRootFront.work.Height(), localRootFront.work.Width(), leafGrid, 0, 0,
      localRootFront.work );
    
    // Perform the distributed portion of the forward solve
    for( Int s=1; s<numDistNodes; ++s )
    {
        const DistSymmNodeInfo& node = info.distNodes[s];
        const DistSymmFront<F>& front = L.distFronts[s];
        const Grid& grid = front.front2dL.Grid();
        const int gridHeight = grid.Height();
        const int gridWidth = grid.Width();
        mpi::Comm comm = grid.VCComm();
        const int commSize = mpi::Size( comm );

        const DistSymmNodeInfo& childNode = info.distNodes[s-1];
        const DistSymmFront<F>& childFront = L.distFronts[s-1];
        const Grid& childGrid = childFront.front2dL.Grid();
        const int childGridHeight = childGrid.Height();
        const int childGridWidth = childGrid.Width();

        // Set up a workspace
        const Int frontHeight = front.front2dL.Height();
        DistMatrix<F>& W = front.work2d;
        W.SetGrid( grid );
        W.Resize( frontHeight, width );
        DistMatrix<F> WT(grid), WB(grid);
        PartitionDown( W, WT, WB, node.size );
        WT = X.distNodes[s-1];
        Zero( WB );

        // Pack our child's update
        const MatrixCommMeta& commMeta = X.commMetas[s-1];
        DistMatrix<F>& childW = childFront.work2d;
        const Int updateSize = childW.Height()-childNode.size;
        DistMatrix<F> childUpdate( childW.Grid() );
        LockedView( childUpdate, childW, childNode.size, 0, updateSize, width );
        int sendBufferSize = 0;
        std::vector<int> sendCounts(commSize), sendDispls(commSize);
        for( int proc=0; proc<commSize; ++proc )
        {
            const int sendSize = commMeta.numChildSendInds[proc];
            sendCounts[proc] = sendSize;
            sendDispls[proc] = sendBufferSize;
            sendBufferSize += sendSize;
        }
        std::vector<F> sendBuffer( sendBufferSize );

        // Pack send data
        const bool onLeft = childNode.onLeft;
        const std::vector<Int>& myChildRelInds = 
            ( onLeft ? node.leftRelInds : node.rightRelInds );
        const int colShift = childUpdate.ColShift();
        const int rowShift = childUpdate.RowShift();
        const Int localWidth = childUpdate.LocalWidth();
        const Int localHeight = childUpdate.LocalHeight();
        auto packOffs = sendDispls;
        for( Int iChildLoc=0; iChildLoc<localHeight; ++iChildLoc )
        {
            const Int iChild = colShift + iChildLoc*childGridHeight;
            const int destRow = myChildRelInds[iChild] % gridHeight;
            for( Int jChildLoc=0; jChildLoc<localWidth; ++jChildLoc )
            {
                const Int jChild = rowShift + jChildLoc*childGridWidth;
                const int destCol = jChild % gridWidth;
                const int destRank = destRow + destCol*gridHeight;
                sendBuffer[packOffs[destRank]++] = 
                    childUpdate.GetLocal(iChildLoc,jChildLoc);
            }
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
            const int recvSize = commMeta.childRecvInds[proc].size()/2;
            recvCounts[proc] = recvSize;
            recvDispls[proc] = recvBufferSize;
            recvBufferSize += recvSize;
        }
        std::vector<F> recvBuffer( recvBufferSize );
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
            const F* recvVals = &recvBuffer[recvDispls[proc]];
            const std::vector<Int>& recvInds = commMeta.childRecvInds[proc];
            for( unsigned k=0; k<recvInds.size()/2; ++k )
            {
                const Int iFrontLoc = recvInds[2*k+0];
                const Int jLoc = recvInds[2*k+1];
                W.UpdateLocal( iFrontLoc, jLoc, recvVals[k] );
            }
        }
        SwapClear( recvBuffer );
        SwapClear( recvCounts );
        SwapClear( recvDispls );

        // Now that the RHS is set up, perform this node's solve
        if( frontType == LDL_2D )
            FrontLowerForwardSolve( front.front2dL, W );
        else if( frontType == LDL_SELINV_2D )
            FrontFastLowerForwardSolve( front.front2dL, W );
        else if( frontType == LDL_INTRAPIV_2D )
            FrontIntraPivLowerForwardSolve( front.front2dL, front.piv, W );
        else if( frontType == LDL_INTRAPIV_SELINV_2D )
            FrontFastIntraPivLowerForwardSolve( front.front2dL, front.piv, W );
        else if( blocked )
            FrontBlockLowerForwardSolve( front.front2dL, W );
        else
            LogicError("Unsupported front type");

        // Store this node's portion of the result
        X.distNodes[s-1] = WT;
    }
    L.localFronts.back().work.Empty();
    L.distFronts.back().work2d.Empty();
}

template<typename F>
inline void DistLowerBackwardSolve
( const DistSymmInfo& info, 
  const DistSymmFrontTree<F>& L, DistNodalMultiVec<F>& X, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("DistLowerBackwardSolve"))
    const int numDistNodes = info.distNodes.size();
    const Int width = X.Width();
    const SymmFrontType frontType = L.frontType;
    if( Unfactored(frontType) )
        LogicError("Nonsensical front type for solve");
    if( frontType == LDL_2D || frontType == LDL_INTRAPIV_2D )
        LogicError("2D non-inverted solves supported by this routine");
    if( frontType == BLOCK_LDL_1D || frontType == BLOCK_LDL_INTRAPIV_1D )
        LogicError("1D block solves not yet supported by this routine");
    const bool frontsAre1d = FrontsAre1d( frontType );
    const bool blocked = BlockFactorization( frontType );
    const bool pivoted = PivotedFactorization( frontType );

    // Directly operate on the root separator's portion of the right-hand sides
    const SymmFront<F>& localRootFront = L.localFronts.back();
    if( numDistNodes == 1 )
    {
        View( localRootFront.work, X.localNodes.back() );
        if( blocked )
            FrontBlockLowerBackwardSolve
            ( localRootFront.frontL, localRootFront.work, conjugate );
        else if( pivoted )
            FrontIntraPivLowerBackwardSolve
            ( localRootFront.frontL, localRootFront.piv, localRootFront.work, 
              conjugate );
        else
            FrontLowerBackwardSolve
            ( localRootFront.frontL, localRootFront.work, conjugate );
    }
    else
    {
        const DistSymmFront<F>& rootFront = L.distFronts.back();
        View( rootFront.work1d, X.distNodes.back() );
        auto& W = rootFront.work1d;
        if( frontType == LDL_1D )
            FrontLowerBackwardSolve( rootFront.front1dL, W, conjugate );
        else if( frontType == LDL_SELINV_1D )
            FrontFastLowerBackwardSolve( rootFront.front1dL, W, conjugate );
        else if( frontType == LDL_SELINV_2D )
            FrontFastLowerBackwardSolve( rootFront.front2dL, W, conjugate );
        else if( frontType == LDL_INTRAPIV_1D )
            FrontIntraPivLowerBackwardSolve
            ( rootFront.front1dL, rootFront.piv, W, conjugate );
        else if( frontType == LDL_INTRAPIV_SELINV_1D )
            FrontFastIntraPivLowerBackwardSolve
            ( rootFront.front1dL, rootFront.piv, W, conjugate );
        else if( frontType == LDL_INTRAPIV_SELINV_2D )
            FrontFastIntraPivLowerBackwardSolve
            ( rootFront.front2dL, rootFront.piv, W, conjugate );
        else if( blocked )
            FrontBlockLowerBackwardSolve( rootFront.front2dL, W, conjugate );
        else
            LogicError("Unsupported front type");
    }

    for( Int s=numDistNodes-2; s>=0; --s )
    {
        const DistSymmNodeInfo& parentNode = info.distNodes[s+1];
        const DistSymmNodeInfo& node = info.distNodes[s];
        const DistSymmFront<F>& parentFront = L.distFronts[s+1];
        const DistSymmFront<F>& front = L.distFronts[s];
        const Grid& grid = ( frontsAre1d ? front.front1dL.Grid() 
                                         : front.front2dL.Grid() );
        const Grid& parentGrid = ( frontsAre1d ? parentFront.front1dL.Grid()
                                               : parentFront.front2dL.Grid() );
        mpi::Comm comm = grid.VCComm(); 
        mpi::Comm parentComm = parentGrid.VCComm();
        const int commSize = mpi::Size( comm );
        const int parentCommSize = mpi::Size( parentComm );
        const Int frontHeight = ( frontsAre1d ? front.front1dL.Height()
                                              : front.front2dL.Height() );

        // Set up a workspace
        DistMatrix<F,VC,STAR>& W = front.work1d;
        W.SetGrid( grid );
        W.Resize( frontHeight, width );
        DistMatrix<F,VC,STAR> WT(grid), WB(grid);
        PartitionDown( W, WT, WB, node.size );
        Matrix<F>& XT = 
          ( s>0 ? X.distNodes[s-1].Matrix() : X.localNodes.back() );
        WT.Matrix() = XT;

        //
        // Set the bottom from the parent
        //

        // Pack the updates using the recv approach from the forward solve
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
        std::vector<F> sendBuffer( sendBufferSize );

        DistMatrix<F,VC,STAR>& parentWork = parentFront.work1d;
        for( int proc=0; proc<parentCommSize; ++proc )
        {
            F* sendVals = &sendBuffer[sendDispls[proc]];
            const std::vector<Int>& recvInds = commMeta.childRecvInds[proc];
            for( unsigned k=0; k<recvInds.size(); ++k )
                StridedMemCopy
                ( &sendVals[k*width],                     1,
                  parentWork.LockedBuffer(recvInds[k],0), parentWork.LDim(),
                  width );
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
        std::vector<F> recvBuffer( recvBufferSize );
        DEBUG_ONLY(VerifySendsAndRecvs( sendCounts, recvCounts, parentComm ))

        // AllToAll to send and recv parent updates
        SparseAllToAll
        ( sendBuffer, sendCounts, sendDispls,
          recvBuffer, recvCounts, recvDispls, parentComm );
        SwapClear( sendBuffer );
        SwapClear( sendCounts );
        SwapClear( sendDispls );

        // Unpack the updates using the send approach from the forward solve
        const bool onLeft = node.onLeft;
        const std::vector<Int>& myRelInds = 
            ( onLeft ? parentNode.leftRelInds : parentNode.rightRelInds );
        const int colShift = WB.ColShift();
        const Int localHeight = WB.LocalHeight();
        for( Int iUpdateLoc=0; iUpdateLoc<localHeight; ++iUpdateLoc )
        {
            const Int iUpdate = colShift + iUpdateLoc*commSize;
            const int startRank = myRelInds[iUpdate] % parentCommSize;
            const F* recvBuf = &recvBuffer[recvDispls[startRank]];
            for( Int j=0; j<width; ++j )
                WB.SetLocal(iUpdateLoc,j,recvBuf[j]);
            recvDispls[startRank] += width;
        }
        SwapClear( recvBuffer );
        SwapClear( recvCounts );
        SwapClear( recvDispls );

        // Call the custom node backward solve
        if( s > 0 )
        {
            if( frontType == LDL_1D )
                FrontLowerBackwardSolve( front.front1dL, W, conjugate );
            else if( frontType == LDL_SELINV_1D )
                FrontFastLowerBackwardSolve( front.front1dL, W, conjugate );
            else if( frontType == LDL_SELINV_2D )
                FrontFastLowerBackwardSolve( front.front2dL, W, conjugate );
            else if( frontType == LDL_INTRAPIV_1D )
                FrontIntraPivLowerBackwardSolve
                ( front.front1dL, front.piv, W, conjugate );
            else if( frontType == LDL_INTRAPIV_SELINV_1D )
                FrontFastIntraPivLowerBackwardSolve
                ( front.front1dL, front.piv, W, conjugate );
            else if( frontType == LDL_INTRAPIV_SELINV_2D )
                FrontFastIntraPivLowerBackwardSolve
                ( front.front2dL, front.piv, W, conjugate );
            else if( blocked )
                FrontBlockLowerBackwardSolve( front.front2dL, W, conjugate );
            else
                LogicError("Unsupported front type");
        }
        else
        {
            View( localRootFront.work, W.Matrix() );
            if( blocked )
                FrontBlockLowerBackwardSolve
                ( localRootFront.frontL, localRootFront.work, conjugate );
            else if( pivoted )
                FrontIntraPivLowerBackwardSolve
                ( localRootFront.frontL, localRootFront.piv, 
                  localRootFront.work, conjugate );
            else
                FrontLowerBackwardSolve
                ( localRootFront.frontL, localRootFront.work, conjugate );
        }

        // Store this node's portion of the result
        XT = WT.Matrix();
    }
}

template<typename F>
inline void DistLowerBackwardSolve
( const DistSymmInfo& info, 
  const DistSymmFrontTree<F>& L, DistNodalMatrix<F>& X, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("DistLowerBackwardSolve"))
    const int numDistNodes = info.distNodes.size();
    const Int width = X.Width();
    const SymmFrontType frontType = L.frontType;
    if( Unfactored(frontType) )
        LogicError("Nonsensical front type for solve");
    if( FrontsAre1d(frontType) )
        LogicError("1d solves not yet implemented");
    const bool blocked = BlockFactorization( frontType );
    const bool pivoted = PivotedFactorization( frontType );

    // Directly operate on the root separator's portion of the right-hand sides
    const SymmFront<F>& localRootFront = L.localFronts.back();
    if( numDistNodes == 1 )
    {
        View( localRootFront.work, X.localNodes.back() );
        if( blocked )
            FrontBlockLowerBackwardSolve
            ( localRootFront.frontL, localRootFront.work, conjugate );
        else if( pivoted )
            FrontIntraPivLowerBackwardSolve
            ( localRootFront.frontL, localRootFront.piv, localRootFront.work,
              conjugate );
        else
            FrontLowerBackwardSolve
            ( localRootFront.frontL, localRootFront.work, conjugate );
    }
    else
    {
        const DistSymmFront<F>& rootFront = L.distFronts.back();
        View( rootFront.work2d, X.distNodes.back() );

        auto& frontL = rootFront.front2dL;
        auto& piv = rootFront.piv;
        auto& W = rootFront.work2d;
        if( frontType == LDL_2D )
            FrontLowerBackwardSolve( frontL, W, conjugate );
        else if( frontType == LDL_SELINV_2D )
            FrontFastLowerBackwardSolve( frontL, W, conjugate );
        else if( frontType == LDL_INTRAPIV_2D )
            FrontIntraPivLowerBackwardSolve( frontL, piv, W, conjugate );
        else if( frontType == LDL_INTRAPIV_SELINV_2D )
            FrontFastIntraPivLowerBackwardSolve( frontL, piv, W, conjugate );
        else if( blocked )
            FrontBlockLowerBackwardSolve( frontL, W, conjugate );
        else
            LogicError("Unsupported front type");
    }

    for( Int s=numDistNodes-2; s>=0; --s )
    {
        const DistSymmNodeInfo& parentNode = info.distNodes[s+1];
        const DistSymmNodeInfo& node = info.distNodes[s];
        const DistSymmFront<F>& parentFront = L.distFronts[s+1];
        const DistSymmFront<F>& front = L.distFronts[s];
        const Grid& grid = front.front2dL.Grid();
        const int gridHeight = grid.Height(); 
        const int gridWidth = grid.Width();
        const Grid& parentGrid = parentFront.front2dL.Grid();
        const int parentGridHeight = parentGrid.Height();
        const int parentGridWidth = parentGrid.Width();
        mpi::Comm parentComm = parentGrid.VCComm();
        const int parentCommSize = mpi::Size( parentComm );
        const Int frontHeight = front.front2dL.Height();

        // Set up a workspace
        DistMatrix<F>& W = front.work2d;
        W.SetGrid( grid );
        W.Resize( frontHeight, width );
        DistMatrix<F> WT(grid), WB(grid);
        PartitionDown( W, WT, WB, node.size );
        Matrix<F>& XT = 
          ( s>0 ? X.distNodes[s-1].Matrix() : X.localNodes.back() );
        WT.Matrix() = XT;

        //
        // Set the bottom from the parent
        //

        // Pack the updates using the recv approach from the forward solve
        const MatrixCommMeta& commMeta = X.commMetas[s];
        int sendBufferSize = 0;
        std::vector<int> sendCounts(parentCommSize), sendDispls(parentCommSize);
        for( int proc=0; proc<parentCommSize; ++proc )
        {
            // childRecvInds contains pairs of indices, but we will send one
            // floating-point value per pair
            const int sendSize = commMeta.childRecvInds[proc].size()/2;
            sendCounts[proc] = sendSize;
            sendDispls[proc] = sendBufferSize;
            sendBufferSize += sendSize;
        }
        std::vector<F> sendBuffer( sendBufferSize );

        DistMatrix<F>& parentWork = parentFront.work2d;
        for( int proc=0; proc<parentCommSize; ++proc )
        {
            F* sendVals = &sendBuffer[sendDispls[proc]];
            const std::vector<Int>& recvInds = commMeta.childRecvInds[proc];
            for( unsigned k=0; k<recvInds.size()/2; ++k )
            {
                const Int iFrontLoc = recvInds[2*k+0];
                const Int jLoc = recvInds[2*k+1];
                sendVals[k] = parentWork.GetLocal( iFrontLoc, jLoc );
            }
        }
        parentWork.Empty();

        // Set up the receive buffer
        int recvBufferSize = 0;
        std::vector<int> recvCounts(parentCommSize), recvDispls(parentCommSize);
        for( int proc=0; proc<parentCommSize; ++proc )
        {
            const int recvSize = commMeta.numChildSendInds[proc];
            recvCounts[proc] = recvSize;
            recvDispls[proc] = recvBufferSize;
            recvBufferSize += recvSize;
        }
        std::vector<F> recvBuffer( recvBufferSize );
        DEBUG_ONLY(VerifySendsAndRecvs( sendCounts, recvCounts, parentComm ))

        // AllToAll to send and recv parent updates
        SparseAllToAll
        ( sendBuffer, sendCounts, sendDispls,
          recvBuffer, recvCounts, recvDispls, parentComm );
        SwapClear( sendBuffer );
        SwapClear( sendCounts );
        SwapClear( sendDispls );

        // Unpack the updates using the send approach from the forward solve
        const bool onLeft = node.onLeft;
        const std::vector<Int>& myRelInds = 
            ( onLeft ? parentNode.leftRelInds : parentNode.rightRelInds );
        const int colShift = WB.ColShift();
        const int rowShift = WB.RowShift();
        const Int localHeight = WB.LocalHeight();
        const Int localWidth = WB.LocalWidth();
        for( Int iUpdateLoc=0; iUpdateLoc<localHeight; ++iUpdateLoc )
        {
            const Int iUpdate = colShift + iUpdateLoc*gridHeight;
            const int startRow = myRelInds[iUpdate] % parentGridHeight;
            for( int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const Int jUpdate = rowShift + jLoc*gridWidth;
                const int startCol = jUpdate % parentGridWidth;
                const int startRank = startRow + startCol*parentGridHeight;
                WB.SetLocal
                ( iUpdateLoc, jLoc, recvBuffer[recvDispls[startRank]++] );
            }
        }
        SwapClear( recvBuffer );
        SwapClear( recvCounts );
        SwapClear( recvDispls );

        // Call the custom node backward solve
        if( s > 0 )
        {
            auto& frontL = front.front2dL;
            auto& piv = front.piv;
            if( frontType == LDL_2D )
                FrontLowerBackwardSolve( frontL, W, conjugate );
            else if( frontType == LDL_SELINV_2D )
                FrontFastLowerBackwardSolve( frontL, W, conjugate );
            else if( frontType == LDL_INTRAPIV_2D )
                FrontIntraPivLowerBackwardSolve( frontL, piv, W, conjugate );
            else if( frontType == LDL_INTRAPIV_SELINV_2D )
                FrontFastIntraPivLowerBackwardSolve
                ( frontL, piv, W, conjugate );
            else if( blocked )
                FrontBlockLowerBackwardSolve( frontL, W, conjugate );
            else
                LogicError("Unsupported front type");
        }
        else
        {
            View( localRootFront.work, W.Matrix() );
            if( blocked )
                FrontBlockLowerBackwardSolve
                ( localRootFront.frontL, localRootFront.work, conjugate );
            else if( pivoted )
                FrontIntraPivLowerBackwardSolve
                ( localRootFront.frontL, localRootFront.piv, 
                  localRootFront.work, conjugate );
            else
                FrontLowerBackwardSolve
                ( localRootFront.frontL, localRootFront.work, conjugate );
        }

        // Store this node's portion of the result
        XT = WT.Matrix();
    }
}

} // namespace El

#endif // ifndef EL_SPARSEDIRECT_NUMERIC_LOWERSOLVE_DIST_HPP
