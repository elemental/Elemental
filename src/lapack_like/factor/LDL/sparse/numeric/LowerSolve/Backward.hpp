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
#ifndef EL_FACTOR_LDL_NUMERIC_LOWERSOLVE_BACKWARD_HPP
#define EL_FACTOR_LDL_NUMERIC_LOWERSOLVE_BACKWARD_HPP

#include "./FrontBackward.hpp"

namespace El {
namespace ldl {

template<typename F> 
inline void LowerBackwardSolve
( const NodeInfo& info, 
  const Front<F>& front,
        MatrixNode<F>& X, bool conjugate )
{
    EL_DEBUG_CSE

    auto* dupMV = X.duplicateMV;
    auto* dupMat = X.duplicateMat;
    const bool haveParent = X.parent != nullptr;
    bool haveDupMVParent = dupMV != nullptr && dupMV->parent != nullptr;
    bool haveDupMatParent = dupMat != nullptr && dupMat->parent != nullptr;
    // Ugh. Surely this can be simplified.
    auto& W = 
      (haveParent ? X.work 
                  : (haveDupMVParent ? dupMV->work.Matrix() 
                                     : (haveDupMatParent ? dupMat->work.Matrix()
                                                         : X.matrix)));

    FrontLowerBackwardSolve( front, W, conjugate );

    const Int numRHS = X.matrix.Width();
    if( haveParent || haveDupMVParent || haveDupMatParent )
        X.matrix = W( IR(0,info.size), IR(0,numRHS) );

    const Int numChildren = front.children.size();
    for( Int c=0; c<numChildren; ++c )
    {
        // Set up a workspace for the child
        auto& childW = X.children[c]->work;
        childW.Resize( front.children[c]->Height(), numRHS );
        const Int childSize = info.children[c]->size;
        auto childWT = childW( IR(0,childSize), ALL );
        auto childWB = childW( IR(childSize,END), ALL );
        childWT = X.children[c]->matrix;

        // Update the child's workspace
        const Int childUSize = childWB.Height();
        for( Int iChild=0; iChild<childUSize; ++iChild )
        {
            const Int i = info.childRelInds[c][iChild];
            for( Int j=0; j<numRHS; ++j )
                childWB(iChild,j) = W(i,j);
        }
    }
    if( haveParent )
        X.work.Empty();
    else if( haveDupMVParent )
        dupMV->work.Empty();
    else if( haveDupMatParent )
        dupMat->work.Empty();

    for( Int c=0; c<numChildren; ++c )
        LowerBackwardSolve
        ( *info.children[c], *front.children[c], *X.children[c], conjugate );
}

template<typename F>
inline void LowerBackwardSolve
( const DistNodeInfo& info,
  const DistFront<F>& front, DistMultiVecNode<F>& X, bool conjugate )
{
    EL_DEBUG_CSE
    if( front.duplicate != nullptr )
    {
        LowerBackwardSolve
        ( *info.duplicate, *front.duplicate, *X.duplicate, conjugate );
        return;
    }

    const bool haveParent = ( X.parent != nullptr );
    auto& W = ( haveParent ? X.work : X.matrix );
    FrontLowerBackwardSolve( front, W, conjugate );

    const Int numRHS = X.matrix.Width();
    if( haveParent )
        X.matrix = W( IR(0,info.size), IR(0,numRHS) );

    if( front.child == nullptr )
        return;

    // Set up a workspace for our child
    const bool frontIs1D = FrontIs1D( front.type );
    const auto& childFront = *front.child;
    EL_DEBUG_ONLY(
      if( FrontIs1D(front.type) != FrontIs1D(childFront.type) )
          LogicError("Incompatible front type mixture");
    )
    const Grid& childGrid =
      ( frontIs1D ? childFront.L1D.Grid() : childFront.L2D.Grid() );
    const Int childFrontHeight =
      ( frontIs1D ? childFront.L1D.Height() : childFront.L2D.Height() );
    auto& childW = X.child->work;
    childW.SetGrid( childGrid );
    childW.Resize( childFrontHeight, numRHS );
    const Int childSize = info.child->size;
    auto childWT = childW( IR(0,childSize), ALL );
    auto childWB = childW( IR(childSize,END), ALL );
    childWT = X.child->matrix;

    // Compute the metadata for updating the children
    X.ComputeCommMeta( info );
    mpi::Comm comm = W.DistComm();
    const int commSize = mpi::Size( comm );
    vector<int> sendSizes(commSize), recvSizes(commSize);
    for( int q=0; q<commSize; ++q )
    {
        sendSizes[q] = X.commMeta.childRecvInds[q].size()*numRHS;
        recvSizes[q] = X.commMeta.numChildSendInds[q]*numRHS;
    }
    EL_DEBUG_ONLY(VerifySendsAndRecvs( sendSizes, recvSizes, comm ))
    vector<int> sendOffs, recvOffs;
    const int sendBufSize = Scan( sendSizes, sendOffs );
    const int recvBufSize = Scan( recvSizes, recvOffs );

    // Pack the updates for the children
    vector<F> sendBuf( sendBufSize );
    Matrix<F>& WLoc = W.Matrix();
    for( int q=0; q<commSize; ++q )
    {
        F* sendVals = &sendBuf[sendOffs[q]];
        const auto& recvInds = X.commMeta.childRecvInds[q];
        const Int numRecvInds = recvInds.size();
        for( Int k=0; k<numRecvInds; ++k )
        {
            const Int recvInd = recvInds[k];
            for( Int j=0; j<numRHS; ++j )
                sendVals[k*numRHS+j] = WLoc(recvInd,j);
        }
    }
    if( haveParent )
        W.Empty();

    // AllToAll to send and recv parent updates
    vector<F> recvBuf( recvBufSize );
    SparseAllToAll
    ( sendBuf, sendSizes, sendOffs,
      recvBuf, recvSizes, recvOffs, comm );
    SwapClear( sendBuf );
    SwapClear( sendSizes );
    SwapClear( sendOffs );

    // Unpack the updates using the send approach from the forward solve
    const Int myChild = ( info.child->onLeft ? 0 : 1 );
    const Int localHeight = childWB.LocalHeight();
    Matrix<F>& childWBLoc = childWB.Matrix();
    for( Int iUpdateLoc=0; iUpdateLoc<localHeight; ++iUpdateLoc )
    {
        const Int iUpdate = childWB.GlobalRow(iUpdateLoc);
        const int q = W.RowOwner(info.childRelInds[myChild][iUpdate]);
        for( Int j=0; j<numRHS; ++j )
            childWBLoc(iUpdateLoc,j) = recvBuf[recvOffs[q]++];
    }
    SwapClear( recvBuf );
    SwapClear( recvSizes );
    SwapClear( recvOffs );

    LowerBackwardSolve( *info.child, *front.child, *X.child, conjugate );
}

template<typename F>
inline void LowerBackwardSolve
( const DistNodeInfo& info,
  const DistFront<F>& front,
        DistMatrixNode<F>& X, bool conjugate )
{
    EL_DEBUG_CSE
    if( front.duplicate != nullptr )
    {
        LowerBackwardSolve
        ( *info.duplicate, *front.duplicate, *X.duplicate, conjugate );
        return;
    }

    const bool haveParent = ( X.parent != nullptr );
    auto& W = ( haveParent ? X.work : X.matrix );
    FrontLowerBackwardSolve( front, W, conjugate );

    const Int numRHS = X.matrix.Width();
    if( haveParent )
        X.matrix = W( IR(0,info.size), IR(0,numRHS) );

    if( front.child == nullptr )
        return;

    // Set up a workspace for our child
    const auto& childFront = *front.child;
    EL_DEBUG_ONLY(
      if( FrontIs1D(front.type) != FrontIs1D(childFront.type) )
          LogicError("Incompatible front type mixture");
    )
    const Grid& childGrid = childFront.L2D.Grid();
    const Int childFrontHeight = childFront.L2D.Height();
    auto& childW = X.child->work;
    childW.SetGrid( childGrid );
    childW.Align( 0, 0 );
    childW.Resize( childFrontHeight, numRHS );
    const Int childSize = info.child->size;
    auto childWT = childW( IR(0,childSize), ALL );
    auto childWB = childW( IR(childSize,END), ALL );
    childWT = X.child->matrix;

    // Set up the metadata for the child updates
    X.ComputeCommMeta( info );
    mpi::Comm comm = W.DistComm();
    const int commSize = mpi::Size( comm );
    vector<int> sendSizes(commSize), recvSizes(commSize);
    for( int q=0; q<commSize; ++q )
    {
        sendSizes[q] = X.commMeta.childRecvInds[q].size()/2;
        recvSizes[q] = X.commMeta.numChildSendInds[q];
    }
    EL_DEBUG_ONLY(VerifySendsAndRecvs( sendSizes, recvSizes, comm ))
    vector<int> sendOffs, recvOffs;
    const int sendBufSize = Scan( sendSizes, sendOffs );
    const int recvBufSize = Scan( recvSizes, recvOffs );

    // Pack the updates
    vector<F> sendBuf( sendBufSize );
    Matrix<F>& WLoc = W.Matrix();
    for( int q=0; q<commSize; ++q )
    {
        F* sendVals = &sendBuf[sendOffs[q]];
        const auto& recvInds = X.commMeta.childRecvInds[q];
        const Int numRecvInds = recvInds.size();
        for( Int k=0; k<numRecvInds/2; ++k )
        {
            const Int iLoc = recvInds[2*k+0];
            const Int jLoc = recvInds[2*k+1];
            sendVals[k] = WLoc(iLoc,jLoc);
        }
    }
    if( haveParent )
        W.Empty();

    // AllToAll to send and recv parent updates
    vector<F> recvBuf( recvBufSize );
    SparseAllToAll
    ( sendBuf, sendSizes, sendOffs,
      recvBuf, recvSizes, recvOffs, comm );
    SwapClear( sendBuf );
    SwapClear( sendSizes );
    SwapClear( sendOffs );

    // Unpack the updates using the send approach from the forward solve
    const auto& childInfo = *info.child;
    const Int myChild = ( childInfo.onLeft ? 0 : 1 );
    const Int localHeight = childWB.LocalHeight();
    const Int localWidth = childWB.LocalWidth();
    Matrix<F>& childWBLoc = childWB.Matrix();
    for( Int iUpdateLoc=0; iUpdateLoc<localHeight; ++iUpdateLoc )
    {
        const Int iUpdate = childWB.GlobalRow(iUpdateLoc);
        const Int i = info.childRelInds[myChild][iUpdate];
        for( int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const Int j = childWB.GlobalCol(jLoc);
            const int q = W.Owner( i, j );
            childWBLoc(iUpdateLoc,jLoc) = recvBuf[recvOffs[q]++];
        }
    }
    SwapClear( recvBuf );
    SwapClear( recvSizes );
    SwapClear( recvOffs );

    LowerBackwardSolve( *info.child, *front.child, *X.child, conjugate );
}

} // namespace ldl
} // namespace El

#endif // EL_FACTOR_LDL_NUMERIC_LOWERSOLVE_BACKWARD_HPP
