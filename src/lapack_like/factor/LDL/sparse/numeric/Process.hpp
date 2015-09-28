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
#ifndef EL_LDL_PROCESS_HPP
#define EL_LDL_PROCESS_HPP

#include "ElSuiteSparse/ldl.hpp"
#include "./ProcessFront.hpp"

namespace El {
namespace ldl {

template<typename F> 
inline void 
Process( const NodeInfo& info, Front<F>& front, LDLFrontType factorType )
{
    DEBUG_ONLY(CSE cse("ldl::Process"))
    const int updateSize = info.lowerStruct.size();
    auto& FBR = front.workDense;
    FBR.Empty();
    Zeros( FBR, updateSize, updateSize );

    if( front.sparseLeaf )
    {
        front.type = factorType;
        const Int m = front.LDense.Height();
        const Int n = front.LDense.Width();
        const Int numEntries = info.LOffsets.back();
        const Int numSources = info.LOffsets.size()-1;

        // TODO: Add support for pivoting here
        if( PivotedFactorization(factorType) )
            Zeros( front.subdiag, n-1, 1 );

        Zeros( front.LSparse, numSources, numSources );
        front.LSparse.ForceNumEntries( numEntries );
        F* LValBuf = front.LSparse.ValueBuffer();
        Int* LRowBuf = front.LSparse.SourceBuffer();
        Int* LColBuf = front.LSparse.TargetBuffer();
        Int* LOffsetBuf = front.LSparse.OffsetBuffer();

        for( Int i=0; i<numSources; ++i )
        {
            const Int iStart = info.LOffsets[i];
            const Int iEnd = info.LOffsets[i+1];
            LOffsetBuf[i] = iStart;
            for( Int e=iStart; e<iEnd; ++e )
                LRowBuf[e] = i;
        }
        LOffsetBuf[numSources] = info.LOffsets[numSources];
        front.diag.Resize( numSources, 1 );

        // Factor the transpose of L
        // TODO: Reuse these workspaces
        vector<Int> LNnz(numSources), pattern(numSources), flag(numSources);
        vector<F> y(numSources);
        suite_sparse::ldl::Numeric
        ( numSources,
          front.workSparse.LockedOffsetBuffer(),
          front.workSparse.LockedTargetBuffer(),
          front.workSparse.LockedValueBuffer(),
          LOffsetBuf,
          info.LParents.data(),
          LNnz.data(),
          LColBuf,
          LValBuf,
          front.diag.Buffer(),
          y.data(),
          pattern.data(),
          flag.data(),
          (const Int*)nullptr,
          (const Int*)nullptr,
          front.isHermitian );
        front.LSparse.ForceConsistency();

        // Solve against L_{TL}^T from the right
        bool onLeft = false;
        suite_sparse::ldl::LTSolveMulti
        ( onLeft, m, n, front.LDense.Buffer(), front.LDense.LDim(),
          LOffsetBuf, LColBuf, LValBuf, front.isHermitian );

        // Save a copy of ABL
        auto ABLCopy = front.LDense;

        // Solve against the diagonal
        suite_sparse::ldl::DSolveMulti
        ( onLeft, m, n, front.LDense.Buffer(), front.LDense.LDim(),
          front.diag.Buffer() );

        // Form the Schur complement
        Orientation orientation = ( front.isHermitian ? ADJOINT : TRANSPOSE );
        Trrk
        ( LOWER, NORMAL, orientation,
          F(-1), front.LDense, ABLCopy, F(0), front.workDense );
    }
    else
    {
        auto& FL = front.LDense;
        DEBUG_ONLY(
          if( FL.Height() != info.size+updateSize || FL.Width() != info.size )
              LogicError("Front was not the proper size");
        )

        // Process children and add in their updates
        const int numChildren = info.children.size();
        for( Int c=0; c<numChildren; ++c )
        {
            Process( *info.children[c], *front.children[c], factorType );

            auto& childU = front.children[c]->workDense;
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
        ProcessFront( front, factorType );
    }
}

template<typename F>
inline void
Process
( const DistNodeInfo& info, DistFront<F>& front, LDLFrontType factorType )
{
    DEBUG_ONLY(CSE cse("ldl::Process"))

    // Switch to a sequential algorithm if possible
    if( front.duplicate != nullptr )
    {
        const Grid& grid = *info.grid;
        auto& frontDup = *front.duplicate;

        Process( *info.duplicate, frontDup, factorType );

        // Pull the relevant information up from the duplicate
        front.type = frontDup.type;
        front.work.LockedAttach( grid, frontDup.workDense );
        if( !BlockFactorization(factorType) )
        {
            front.diag.LockedAttach( grid, frontDup.diag );
            if( PivotedFactorization(factorType) )
            {
                front.piv.LockedAttach( grid, frontDup.piv );
                front.subdiag.LockedAttach( grid, frontDup.subdiag );
            }
        }

        return;
    }

    const auto& childInfo = *info.child;
    auto& childFront = *front.child;
    Process( childInfo, childFront, factorType );

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
    const int commSize = mpi::Size( comm );
    const auto& childU = childFront.work;
    vector<int> sendSizes(commSize), recvSizes(commSize);
    for( int q=0; q<commSize; ++q )
    {
        sendSizes[q] = front.commMeta.numChildSendInds[q];
        recvSizes[q] = front.commMeta.childRecvInds[q].size()/2;
    }
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
        childFront.duplicate->workDense.Empty();

    // AllToAll to send and receive the child updates
    vector<F> recvBuf( recvBufSize );
    DEBUG_ONLY(VerifySendsAndRecvs( sendSizes, recvSizes, comm ))
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

    ProcessFront( front, factorType );
}

} // namespace ldl
} // namespace El

#endif // ifndef EL_LDL_PROCESS_HPP
