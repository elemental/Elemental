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
#ifndef EL_FACTOR_LDL_NUMERIC_LOWERSOLVE_FORWARD_HPP
#define EL_FACTOR_LDL_NUMERIC_LOWERSOLVE_FORWARD_HPP

#include "./FrontForward.hpp"

namespace El {
namespace ldl {

template<typename F> 
void LowerForwardSolve
( const NodeInfo& info, 
  const Front<F>& front,
        MatrixNode<F>& X )
{
    EL_DEBUG_CSE

    const Int numChildren = info.children.size();
    for( Int c=0; c<numChildren; ++c )
        LowerForwardSolve
        ( *info.children[c], *front.children[c], *X.children[c] );

    // Set up a workspace
    // TODO: Only set up a workspace if there is not a parent 
    //       (or a duplicate's parent)
    auto& W = X.work;
    const Int numRHS = X.matrix.Width();
    W.Resize( front.Height(), numRHS );
    auto WT = W( IR(0,info.size), ALL );
    auto WB = W( IR(info.size,END), ALL );
    WT = X.matrix;
    Zero( WB );

    // Update using the children (if they exist)
    for( Int c=0; c<numChildren; ++c )
    {
        auto& childW = X.children[c]->work;
        const Int childSize = info.children[c]->size;
        const Int childHeight = childW.Height();
        const Int childUSize = childHeight-childSize;

        auto childU = childW( IR(childSize,childHeight), IR(0,numRHS) );
        for( Int iChild=0; iChild<childUSize; ++iChild )
        {
            const Int iFront = info.childRelInds[c][iChild]; 
            for( Int j=0; j<numRHS; ++j )
                W(iFront,j) += childU(iChild,j);
        }
        childW.Empty();
    }

    // Solve against this front
    FrontLowerForwardSolve( front, W );

    // Store this node's portion of the result
    X.matrix = WT;
}

template<typename F>
void LowerForwardSolve
( const DistNodeInfo& info,
  const DistFront<F>& front,
        DistMultiVecNode<F>& X )
{
    EL_DEBUG_CSE

    const bool frontIs1D = FrontIs1D( front.type );
    const Grid& grid = ( frontIs1D ? front.L1D.Grid() : front.L2D.Grid() );
    if( front.duplicate != nullptr )
    {
        LowerForwardSolve( *info.duplicate, *front.duplicate, *X.duplicate );
        X.work.LockedAttach( grid, X.duplicate->work );
        return;
    }

    const auto& childInfo = *info.child;
    const auto& childFront = *front.child;
    EL_DEBUG_ONLY(
      if( FrontIs1D(front.type) != FrontIs1D(childFront.type) )
          LogicError("Incompatible front type mixture");
    )

    LowerForwardSolve( childInfo, childFront, *X.child );

    // Set up a workspace
    // TODO: Only set up a workspace if there is a parent
    const Int numRHS = X.matrix.Width();
    const Int frontHeight =
      ( frontIs1D ? front.L1D.Height() : front.L2D.Height() );
    auto& W = X.work;
    W.SetGrid( grid );
    W.Resize( frontHeight, numRHS );
    auto WT = W( IR(0,info.size), ALL );
    auto WB = W( IR(info.size,END), ALL );
    WT = X.matrix;
    Zero( WB );
    mpi::Comm comm = W.DistComm();
    const int commSize = mpi::Size( comm );

    // Compute the metadata for transmitting child updates
    X.ComputeCommMeta( info );
    auto& childW = X.child->work;
    auto childU = childW( IR(childInfo.size,childW.Height()), IR(0,numRHS) );
    vector<int> sendSizes(commSize), recvSizes(commSize);
    for( int q=0; q<commSize; ++q )
    {
        sendSizes[q] = X.commMeta.numChildSendInds[q]*numRHS;
        recvSizes[q] = X.commMeta.childRecvInds[q].size()*numRHS;
    }
    EL_DEBUG_ONLY(VerifySendsAndRecvs( sendSizes, recvSizes, comm ))
    vector<int> sendOffs, recvOffs;
    const int sendBufSize = Scan( sendSizes, sendOffs );
    const int recvBufSize = Scan( recvSizes, recvOffs );

    // Pack our child's update
    vector<F> sendBuf( sendBufSize );
    const Int myChild = ( childInfo.onLeft ? 0 : 1 );
    auto packOffs = sendOffs;
    const Int localHeight = childU.LocalHeight();
    const Matrix<F>& childULoc = childU.LockedMatrix();
    for( Int iChildLoc=0; iChildLoc<localHeight; ++iChildLoc )
    {
        const Int iChild = childU.GlobalRow(iChildLoc);
        const Int q = W.RowOwner( info.childRelInds[myChild][iChild] );
        for( Int j=0; j<numRHS; ++j )
            sendBuf[packOffs[q]++] = childULoc(iChildLoc,j);
    }
    SwapClear( packOffs );
    childW.Empty();
    if( X.child->duplicate != nullptr )
        X.child->duplicate->work.Empty();

    // AllToAll to send and receive the child updates
    vector<F> recvBuf( recvBufSize );
    SparseAllToAll
    ( sendBuf, sendSizes, sendOffs,
      recvBuf, recvSizes, recvOffs, comm );
    SwapClear( sendBuf );
    SwapClear( sendSizes );
    SwapClear( sendOffs );

    // Unpack the child updates
    Matrix<F>& WLoc = W.Matrix();
    for( int q=0; q<commSize; ++q )
    {
        const F* recvVals = &recvBuf[recvOffs[q]];
        const auto& recvInds = X.commMeta.childRecvInds[q];
        const Int numRecvInds = recvInds.size();
        for( Int k=0; k<numRecvInds; ++k )
            for( Int j=0; j<numRHS; ++j )
                WLoc(recvInds[k],j) += recvVals[k*numRHS+j]; 
    }
    SwapClear( recvBuf );
    SwapClear( recvSizes );
    SwapClear( recvOffs );

    // Now that the RHS is set up, perform this node's solve
    FrontLowerForwardSolve( front, W );

    // Unpack the workspace
    X.matrix = WT;
}

template<typename F>
void LowerForwardSolve
( const DistNodeInfo& info,
  const DistFront<F>& front,
        DistMatrixNode<F>& X )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( FrontIs1D(front.type) )
          LogicError("Front was not 2D");
    )

    const Grid& grid = front.L2D.Grid();
    if( front.duplicate != nullptr )
    {
        LowerForwardSolve( *info.duplicate, *front.duplicate, *X.duplicate );
        X.work.LockedAttach( grid, X.duplicate->work );
        return;
    }

    const auto& childInfo = *info.child;
    const auto& childFront = *front.child;
    EL_DEBUG_ONLY(
      if( FrontIs1D(front.type) != FrontIs1D(childFront.type) )
          LogicError("Incompatible front type mixture");
    )

    LowerForwardSolve( childInfo, childFront, *X.child );

    // Set up a workspace
    // TODO: Only set up a workspace if there is a parent
    const Int numRHS = X.matrix.Width();
    const Int frontHeight = front.L2D.Height();
    auto& W = X.work;
    W.SetGrid( grid );
    W.Align( 0, 0 );
    W.Resize( frontHeight, numRHS );
    auto WT = W( IR(0,info.size), ALL );
    auto WB = W( IR(info.size,END), ALL );
    WT = X.matrix;
    Zero( WB );
    mpi::Comm comm = W.DistComm();
    const int commSize = mpi::Size( comm );

    // Compute the metadata
    X.ComputeCommMeta( info );
    auto& childW = X.child->work;
    auto childU = childW( IR(childInfo.size,childW.Height()), IR(0,numRHS) );
    vector<int> sendSizes(commSize), recvSizes(commSize);
    for( int q=0; q<commSize; ++q )
    {
        sendSizes[q] = X.commMeta.numChildSendInds[q];
        recvSizes[q] = X.commMeta.childRecvInds[q].size()/2;
    }
    EL_DEBUG_ONLY(VerifySendsAndRecvs( sendSizes, recvSizes, comm ))
    vector<int> sendOffs, recvOffs;
    const int sendBufSize = Scan( sendSizes, sendOffs );
    const int recvBufSize = Scan( recvSizes, recvOffs );

    // Pack send data
    vector<F> sendBuf( sendBufSize );
    const Int myChild = ( childInfo.onLeft ? 0 : 1 );
    auto packOffs = sendOffs;
    const Int localHeight = childU.LocalHeight();
    const Int localWidth = childU.LocalWidth();
    const Matrix<F>& childULoc = childU.LockedMatrix();
    for( Int iChildLoc=0; iChildLoc<localHeight; ++iChildLoc )
    {
        const Int iChild = childU.GlobalRow(iChildLoc);
        const Int iParent = info.childRelInds[myChild][iChild];
        for( Int jChildLoc=0; jChildLoc<localWidth; ++jChildLoc )
        {
            const Int j = childU.GlobalCol(jChildLoc);
            const int q = W.Owner( iParent, j );
            EL_DEBUG_ONLY(
              if( packOffs[q] >= sendBufSize )
                  LogicError("packOffs[",q,"]=",packOffs[q]," >= ",sendBufSize);
            )
            sendBuf[packOffs[q]++] = childULoc(iChildLoc,jChildLoc);
        }
    }
    SwapClear( packOffs );
    childW.Empty();
    if( X.child->duplicate != nullptr )
        X.child->duplicate->work.Empty();

    // AllToAll to send and receive the child updates
    vector<F> recvBuf( recvBufSize );
    SparseAllToAll
    ( sendBuf, sendSizes, sendOffs,
      recvBuf, recvSizes, recvOffs, comm );
    SwapClear( sendBuf );
    SwapClear( sendSizes );
    SwapClear( sendOffs );

    // Unpack the child updates (with an Axpy)
    Matrix<F>& WLoc = W.Matrix();
    for( int q=0; q<commSize; ++q )
    {
        const auto& recvInds = X.commMeta.childRecvInds[q];
        const Int numRecvInds = recvInds.size();
        for( Int k=0; k<numRecvInds/2; ++k )
        {
            const Int iLoc = recvInds[2*k+0];
            const Int jLoc = recvInds[2*k+1];
            WLoc(iLoc,jLoc) += recvBuf[recvOffs[q]+k];
        }
    }
    SwapClear( recvBuf );
    SwapClear( recvSizes );
    SwapClear( recvOffs );

    // Now that the RHS is set up, perform this node's solve
    FrontLowerForwardSolve( front, W );

    // Store this node's portion of the result
    X.matrix = WT;
}

} // namespace ldl
} // namespace El

#endif // EL_FACTOR_LDL_NUMERIC_LOWERSOLVE_FORWARD_HPP
