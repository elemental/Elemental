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
#pragma once
#ifndef EL_SPARSEDIRECT_NUMERIC_LOWERSOLVE_BACKWARD_HPP
#define EL_SPARSEDIRECT_NUMERIC_LOWERSOLVE_BACKWARD_HPP

#include "./FrontBackward.hpp"

namespace El {

template<typename F> 
inline void LowerBackwardSolve
( const SymmNodeInfo& info, 
  const SymmFront<F>& front, MatrixNode<F>& X, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("LowerBackwardSolve"))

    auto* dup = front.duplicate;
    const bool haveParent = front.parent != nullptr;
    const bool have1DParent = 
      dup != nullptr && FrontIs1D(dup->type) && dup->parent != nullptr;
    const bool have2DParent =
      dup != nullptr && !FrontIs1D(dup->type) && dup->parent != nullptr;
    auto& W = 
      ( haveParent ? front.work 
                   : (have1DParent ? dup->work1D.Matrix() 
                                   : (have2DParent ? dup->work2D.Matrix() 
                                                   : X.matrix ) ) );

    FrontLowerBackwardSolve( front, W, conjugate );

    const Int numRHS = X.matrix.Width();
    if( haveParent || have1DParent || have2DParent )
        X.matrix = W( IR(0,info.size), IR(0,numRHS) );

    const Int numChildren = front.children.size();
    for( Int c=0; c<numChildren; ++c )
    {
        const auto& childFront = *front.children[c];
        auto& childW = childFront.work;

        // Set up a workspace for the child
        childW.Resize( front.children[c]->L.Height(), numRHS );
        Matrix<F> childWT, childWB; 
        PartitionDown( childW, childWT, childWB, info.children[c]->size );
        childWT = X.children[c]->matrix;

        // Update the child's workspace
        const Int childUSize = childWB.Height();
        for( Int iChild=0; iChild<childUSize; ++iChild )
        {
            const Int i = info.childRelInds[c][iChild];
            for( Int j=0; j<numRHS; ++j )
                childWB.Set( iChild, j, W.Get(i,j) );
        }
    }
    if( haveParent )
        front.work.Empty();
    else if( have1DParent )
        dup->work1D.Empty();    
    else if( have2DParent )
        dup->work2D.Empty();

    for( Int c=0; c<numChildren; ++c )
        LowerBackwardSolve
        ( *info.children[c], *front.children[c], *X.children[c], conjugate );
}

template<typename F>
inline void LowerBackwardSolve
( const DistSymmNodeInfo& info,
  const DistSymmFront<F>& front, DistMultiVecNode<F>& X, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("LowerBackwardSolve"))
    if( front.duplicate != nullptr )
    {
        LowerBackwardSolve
        ( *info.duplicate, *front.duplicate, *X.duplicate, conjugate );
        return;
    }

    const bool haveParent = ( front.parent != nullptr );
    auto& W = ( haveParent ? front.work1D : X.matrix );
    FrontLowerBackwardSolve( front, W, conjugate );

    const Int numRHS = X.matrix.Width();
    if( haveParent )
        X.matrix = W( IR(0,info.size), IR(0,numRHS) );

    if( front.child == nullptr )
        return;
    if( front.type != front.child->type )
        LogicError("Expected front types to match");

    // Set up a workspace for our child
    const bool frontIs1D = FrontIs1D( front.type );

    const auto& childFront = *front.child;
    const Grid& childGrid =
        ( frontIs1D ? childFront.L1D.Grid() : childFront.L2D.Grid() );
    const Int childFrontHeight =
        ( frontIs1D ? childFront.L1D.Height() : childFront.L2D.Height() );
    auto& childW = childFront.work1D;
    childW.SetGrid( childGrid );
    childW.Resize( childFrontHeight, numRHS );
    DistMatrix<F,VC,STAR> childWT(childGrid), childWB(childGrid);
    PartitionDown( childW, childWT, childWB, info.child->size );
    childWT = X.child->matrix;

    // Compute the metadata for updating the children
    mpi::Comm comm = W.DistComm();
    const int commSize = mpi::Size( comm );
    const auto& commMeta = info.multiVecMeta;
    vector<int> sendSizes(commSize), recvSizes(commSize);
    for( int q=0; q<commSize; ++q )
    {
        sendSizes[q] = commMeta.childRecvInds[q].size()*numRHS;
        recvSizes[q] = commMeta.numChildSendInds[q]*numRHS;
    }
    vector<int> sendOffs, recvOffs;
    const int sendBufSize = Scan( sendSizes, sendOffs );
    const int recvBufSize = Scan( recvSizes, recvOffs );

    // Pack the updates for the children
    vector<F> sendBuf( sendBufSize );
    for( int q=0; q<commSize; ++q )
    {
        F* sendVals = &sendBuf[sendOffs[q]];
        const auto& recvInds = commMeta.childRecvInds[q];
        for( unsigned k=0; k<recvInds.size(); ++k )
            StridedMemCopy
            ( &sendVals[k*numRHS],           1,
              W.LockedBuffer(recvInds[k],0), W.LDim(),
              numRHS );
    }
    if( haveParent )
        W.Empty();

    // AllToAll to send and recv parent updates
    vector<F> recvBuf( recvBufSize );
    DEBUG_ONLY(VerifySendsAndRecvs( sendSizes, recvSizes, comm ))
    SparseAllToAll
    ( sendBuf, sendSizes, sendOffs,
      recvBuf, recvSizes, recvOffs, comm );
    SwapClear( sendBuf );
    SwapClear( sendSizes );
    SwapClear( sendOffs );

    // Unpack the updates using the send approach from the forward solve
    const auto& childRelInds =
        ( info.child->onLeft ? info.childRelInds[0] : info.childRelInds[1] );
    const Int localHeight = childWB.LocalHeight();
    for( Int iUpdateLoc=0; iUpdateLoc<localHeight; ++iUpdateLoc )
    {
        const Int iUpdate = childWB.GlobalRow(iUpdateLoc);
        const int q = W.RowOwner(childRelInds[iUpdate]);
        for( Int j=0; j<numRHS; ++j )
            childWB.SetLocal( iUpdateLoc, j, recvBuf[recvOffs[q]++] );
    }
    SwapClear( recvBuf );
    SwapClear( recvSizes );
    SwapClear( recvOffs );

    LowerBackwardSolve( *info.child, *front.child, *X.child, conjugate );
}

template<typename F>
inline void LowerBackwardSolve
( const DistSymmNodeInfo& info,
  const DistSymmFront<F>& front, DistMatrixNode<F>& X, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("LowerBackwardSolve"))
    if( front.duplicate != nullptr )
    {
        LowerBackwardSolve
        ( *info.duplicate, *front.duplicate, *X.duplicate, conjugate );
        return;
    }

    const bool haveParent = ( front.parent != nullptr );
    auto& W = ( haveParent ? front.work2D : X.matrix );
    FrontLowerBackwardSolve( front, W, conjugate );

    const Int numRHS = X.matrix.Width();
    if( haveParent )
        X.matrix = W( IR(0,info.size), IR(0,numRHS) );

    if( front.child == nullptr )
        return;
    if( front.type != front.child->type )
        LogicError("Expected front types to match");

    // Set up a workspace for our child
    const auto& childFront = *front.child;
    const Grid& childGrid = childFront.L2D.Grid();
    const Int childFrontHeight = childFront.L2D.Height();
    auto& childW = childFront.work2D;
    childW.SetGrid( childGrid );
    childW.Resize( childFrontHeight, numRHS );
    DistMatrix<F> childWT(childGrid), childWB(childGrid);
    PartitionDown( childW, childWT, childWB, info.child->size );
    childWT = X.child->matrix;

    // Set up the metadata for the child updates
    mpi::Comm comm = W.DistComm();
    const int commSize = mpi::Size( comm );
    const auto& commMeta = X.commMeta;
    vector<int> sendSizes(commSize), recvSizes(commSize);
    for( int q=0; q<commSize; ++q )
    {
        sendSizes[q] = commMeta.childRecvInds[q].size()/2;
        recvSizes[q] = commMeta.numChildSendInds[q];
    }
    vector<int> sendOffs, recvOffs;
    const int sendBufSize = Scan( sendSizes, sendOffs );
    const int recvBufSize = Scan( recvSizes, recvOffs );

    // Pack the updates
    vector<F> sendBuf( sendBufSize );
    for( int q=0; q<commSize; ++q )
    {
        F* sendVals = &sendBuf[sendOffs[q]];
        const auto& recvInds = commMeta.childRecvInds[q];
        for( unsigned k=0; k<recvInds.size()/2; ++k )
        {
            const Int iLoc = recvInds[2*k+0];
            const Int jLoc = recvInds[2*k+1];
            sendVals[k] = W.GetLocal( iLoc, jLoc );
        }
    }
    if( haveParent )
        W.Empty();

    // AllToAll to send and recv parent updates
    vector<F> recvBuf( recvBufSize );
    DEBUG_ONLY(VerifySendsAndRecvs( sendSizes, recvSizes, comm ))
    SparseAllToAll
    ( sendBuf, sendSizes, sendOffs,
      recvBuf, recvSizes, recvOffs, comm );
    SwapClear( sendBuf );
    SwapClear( sendSizes );
    SwapClear( sendOffs );

    // Unpack the updates using the send approach from the forward solve
    const auto& childInfo = *info.child;
    const auto& childRelInds =
        ( childInfo.onLeft ? info.childRelInds[0] : info.childRelInds[1] );
    const Int localHeight = childWB.LocalHeight();
    const Int localWidth = childWB.LocalWidth();
    for( Int iUpdateLoc=0; iUpdateLoc<localHeight; ++iUpdateLoc )
    {
        const Int iUpdate = childWB.GlobalRow(iUpdateLoc);
        for( int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const Int j = childWB.GlobalCol(jLoc);
            const int q = W.Owner( childRelInds[iUpdate], j );
            childWB.SetLocal( iUpdateLoc, jLoc, recvBuf[recvOffs[q]++] );
        }
    }
    SwapClear( recvBuf );
    SwapClear( recvSizes );
    SwapClear( recvOffs );

    LowerBackwardSolve( *info.child, *front.child, *X.child, conjugate );
}

} // namespace El

#endif // ifndef EL_SPARSEDIRECT_NUMERIC_LOWERSOLVE_BACKWARD_HPP
