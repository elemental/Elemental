/*
   Copyright (c) 2009-2015, Jack Poulson.
   All rights reserved.
 
   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_FACTOR_REGQSDLDL_PROCESSDISTTREE_HPP
#define EL_FACTOR_REGQSDLDL_PROCESSDISTTREE_HPP

#include "./ProcessFront.hpp"

namespace El {
namespace reg_qsd_ldl {

template<typename F>
inline void
Process
( const SymmNodeInfo& info, SymmFront<F>& front,
  Base<F> pivTol, 
  const MatrixNode<Base<F>>& regCand, 
        MatrixNode<Base<F>>& reg,
  bool aPriori,
  SymmFrontType factorType )
{
    DEBUG_ONLY(CallStackEntry cse("reg_qsd_ldl::Process"))

    const int updateSize = info.lowerStruct.size();
    auto& FL = front.L;
    auto& FBR = front.work;
    FBR.Empty();
    DEBUG_ONLY(
      if( FL.Height() != info.size+updateSize || FL.Width() != info.size )
          LogicError("Front was not the proper size");
    )

    // Process children and add in their updates
    Zeros( FBR, updateSize, updateSize );
    const int numChildren = info.children.size();
    for( Int c=0; c<numChildren; ++c )
    {
        Process
        ( *info.children[c], *front.children[c], 
          pivTol, *regCand.children[c], *reg.children[c], aPriori,
          factorType );

        auto& childU = front.children[c]->work;
        const int childUSize = childU.Height();
        for( int jChild=0; jChild<childUSize; ++jChild )
        {
            const int j = info.childRelInds[c][jChild];
            for( int iChild=jChild; iChild<childUSize; ++iChild )
            {
                const int i = info.childRelInds[c][iChild];
                const F value = childU.Get(iChild,jChild);
                if( j < info.size )
                    FL.Update( i, j, value );
                else
                    FBR.Update( i-info.size, j-info.size, value );
            }
        }
        childU.Empty();
    }

    ProcessFront
    ( front, pivTol, regCand.matrix, reg.matrix, aPriori, factorType );
}

template<typename F>
inline void
Process
( const DistSymmNodeInfo& info, DistSymmFront<F>& front,
  Base<F> pivTol, 
  const DistMultiVecNode<Base<F>>& regCand, 
        DistMultiVecNode<Base<F>>& reg,
  bool aPriori,
  SymmFrontType factorType )
{
    DEBUG_ONLY(CallStackEntry cse("reg_qsd_ldl::Process"))

    // Switch to a sequential algorithm if possible
    if( front.duplicate != nullptr )
    {
        const Grid& grid = *info.grid;
        auto& frontDup = *front.duplicate;

        Process
        ( *info.duplicate, frontDup, 
          pivTol, *regCand.duplicate, *reg.duplicate, aPriori,
          factorType );

        // Pull the relevant information from the duplicate
        front.type = frontDup.type;
        front.work.LockedAttach( grid, frontDup.work );
        if( !BlockFactorization(front.type) )
        {
            front.diag.LockedAttach( grid, frontDup.diag );
            if( PivotedFactorization(front.type) )
            {
                front.piv.LockedAttach( grid, frontDup.piv );
                front.subdiag.LockedAttach( grid, frontDup.subdiag );
            }
        }
        return;
    }

    const auto& childInfo = *info.child;
    auto& childFront = *front.child;
    Process
    ( childInfo, childFront, pivTol, *regCand.child, *reg.child, aPriori,
      factorType );

    const Int updateSize = info.lowerStruct.size();
    front.work.Empty();
    DEBUG_ONLY(
      if( front.L2D.Height() != info.size+updateSize ||
          front.L2D.Width() != info.size )
          LogicError("Front was not the proper size");
    )

    // Compute the metadata for sharing child updates
    front.ComputeCommMeta( info, true );
    mpi::Comm comm = front.L2D.DistComm();
    const unsigned commSize = mpi::Size( comm );
    const auto& childU = childFront.work;
    vector<int> sendSizes(commSize), recvSizes(commSize);
    for( int q=0; q<commSize; ++q )
    {
        sendSizes[q] = front.commMeta.numChildSendInds[q];
        recvSizes[q] = front.commMeta.childRecvInds[q].size()/2;
    }
    DEBUG_ONLY(VerifySendsAndRecvs( sendSizes, recvSizes, comm ))
    vector<int> sendOffs, recvOffs;
    const int sendBufSize = Scan( sendSizes, sendOffs );
    const int recvBufSize = Scan( recvSizes, recvOffs );

    // Pack the updates
    vector<F> sendBuf( sendBufSize );
    const Int myChild = ( childInfo.onLeft ? 0 : 1 );
    auto offs = sendOffs;
    const Int updateLocHeight = childU.LocalHeight();
    const Int updateLocWidth = childU.LocalWidth();
    for( Int jChildLoc=0; jChildLoc<updateLocWidth; ++jChildLoc )
    {
        const Int jChild = childU.GlobalCol(jChildLoc);
        const Int j = info.childRelInds[myChild][jChild];
        const Int iChildOff = childU.LocalRowOffset( jChild );
        for( Int iChildLoc=iChildOff; iChildLoc<updateLocHeight; ++iChildLoc )
        {
            const Int iChild = childU.GlobalRow(iChildLoc);
            const Int i = info.childRelInds[myChild][iChild];
            const int q = front.L2D.Owner( i, j );
            sendBuf[offs[q]++] = childU.GetLocal(iChildLoc,jChildLoc);
        }
    }
    DEBUG_ONLY(
      for( int q=0; q<commSize; ++q )
      {
          if( offs[q]-sendOffs[q] != front.commMeta.numChildSendInds[q] )
              LogicError("Error in packing stage");
      }
    )
    SwapClear( offs );
    childFront.work.Empty();
    if( childFront.duplicate != nullptr )
        childFront.duplicate->work.Empty();

    // AllToAll to send and receive the child updates
    vector<F> recvBuf( recvBufSize );
    SparseAllToAll
    ( sendBuf, sendSizes, sendOffs,
      recvBuf, recvSizes, recvOffs, comm );
    SwapClear( sendBuf );
    SwapClear( sendSizes );
    SwapClear( sendOffs );

    // Unpack the child udpates (with an Axpy)
    auto& FL = front.L2D;
    auto FTL = FL( IR(0,info.size), IR(0,info.size) );
    auto& FBR = front.work;
    const Int topLocHeight = FTL.LocalHeight();
    const Int leftLocWidth = FTL.LocalWidth();
    FBR.SetGrid( FTL.Grid() );
    FBR.Align( FTL.RowOwner(info.size), FTL.ColOwner(info.size) );
    Zeros( FBR, updateSize, updateSize );
    for( int q=0; q<commSize; ++q )
    {
        const Int numRecvIndPairs = front.commMeta.childRecvInds[q].size()/2;
        for( Int k=0; k<numRecvIndPairs; ++k )
        {
            const Int iLoc = front.commMeta.childRecvInds[q][2*k+0];
            const Int jLoc = front.commMeta.childRecvInds[q][2*k+1];
            const F value = recvBuf[recvOffs[q]+k];
            if( jLoc < leftLocWidth )
                FL.UpdateLocal( iLoc, jLoc, value );
            else
                FBR.UpdateLocal( iLoc-topLocHeight, jLoc-leftLocWidth, value );
        }
    }
    SwapClear( recvBuf );
    SwapClear( recvSizes );
    SwapClear( recvOffs );

    ProcessFront
    ( front, pivTol, regCand.matrix, reg.matrix, aPriori, factorType );
}

} // namespace reg_qsd_ldl
} // namespace El

#endif // ifndef EL_FACTOR_REGQSDLDL_PROCESSDISTTREE_HPP
