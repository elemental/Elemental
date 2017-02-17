/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_BLAS_COPY_GENERALPURPOSE_HPP
#define EL_BLAS_COPY_GENERALPURPOSE_HPP

namespace El {
namespace copy {

template<typename S,typename T,typename=EnableIf<CanCast<S,T>>>
void Helper
( const AbstractDistMatrix<S>& A,
        AbstractDistMatrix<T>& B )
{
    EL_DEBUG_CSE

    // TODO: Decide whether S or T should be used as the transmission type
    //       based upon which is smaller. Transmit S by default.
    const Int height = A.Height();
    const Int width = A.Width();
    const Grid& g = B.Grid();
    B.Resize( height, width );
    Zero( B );
    const bool BPartic = B.Participating();
    const int BRoot = B.Root();

    const bool includeViewers = (A.Grid() != B.Grid());

    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();
    auto& ALoc = A.LockedMatrix();
    auto& BLoc = B.Matrix();

    // TODO: Break into smaller pieces to avoid excessive memory usage?
    vector<Entry<S>> remoteEntries;
    vector<int> distOwners;
    if( A.RedundantRank() == 0 )
    {
        const bool noRedundant = B.RedundantSize() == 1;
        const int colStride = B.ColStride();
        const int rowRank = B.RowRank();
        const int colRank = B.ColRank();

        vector<Int> globalRows(localHeight), localRows(localHeight);
        vector<int> ownerRows(localHeight);
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = A.GlobalRow(iLoc);
            const int ownerRow = B.RowOwner(i);
            globalRows[iLoc] = i;
            ownerRows[iLoc] = ownerRow;
            localRows[iLoc] = B.LocalRow(i,ownerRow);
        }

        remoteEntries.reserve( localHeight*localWidth );
        distOwners.reserve( localHeight*localWidth );
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const Int j = A.GlobalCol(jLoc);
            const int ownerCol = B.ColOwner(j);
            const Int localCol = B.LocalCol(j,ownerCol);
            const bool isLocalCol = ( BPartic && ownerCol == rowRank );
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
            {
                const int ownerRow = ownerRows[iLoc];
                const Int localRow = localRows[iLoc];
                const bool isLocalRow = ( BPartic && ownerRow == colRank );
                const S& alpha = ALoc(iLoc,jLoc);
                if( noRedundant && isLocalRow && isLocalCol )
                {
                    BLoc(localRow,localCol) = Caster<S,T>::Cast(alpha);
                }
                else
                {
                    remoteEntries.push_back
                    ( Entry<S>{localRow,localCol,alpha} );
                    distOwners.push_back( ownerRow + colStride*ownerCol );
                }
            }
        }
    }

    // We will first push to redundant rank 0 of B
    const int redundantRootB = 0;

    // Compute the metadata
    // ====================
    const Int totalSend = remoteEntries.size();
    mpi::Comm comm;
    vector<int> sendCounts, owners(totalSend);
    if( includeViewers )
    {
        comm = g.ViewingComm();
        const int viewingSize = mpi::Size( g.ViewingComm() );
        const int distBSize = mpi::Size( B.DistComm() );

        vector<int> distBToViewing(distBSize);
        for( int distBRank=0; distBRank<distBSize; ++distBRank )
        {
            const int vcOwner =
              g.CoordsToVC
              (B.ColDist(),B.RowDist(),distBRank,BRoot,redundantRootB);
            distBToViewing[distBRank] = g.VCToViewing(vcOwner);
        }

        sendCounts.resize(viewingSize,0);
        for( Int k=0; k<totalSend; ++k )
        {
            owners[k] = distBToViewing[distOwners[k]];
            ++sendCounts[owners[k]];
        }
    }
    else
    {
        if( !g.InGrid() )
            return;
        comm = g.VCComm();

        const int distBSize = mpi::Size( B.DistComm() );
        vector<int> distBToVC(distBSize);
        for( int distBRank=0; distBRank<distBSize; ++distBRank )
        {
            distBToVC[distBRank] =
              g.CoordsToVC
              (B.ColDist(),B.RowDist(),distBRank,BRoot,redundantRootB);
        }

        const int vcSize = mpi::Size( g.VCComm() );
        sendCounts.resize(vcSize,0);
        for( Int k=0; k<totalSend; ++k )
        {
            owners[k] = distBToVC[distOwners[k]];
            ++sendCounts[owners[k]];
        }
    }
    SwapClear( distOwners );

    // Pack the data
    // =============
    vector<int> sendOffs;
    Scan( sendCounts, sendOffs );
    vector<Entry<S>> sendBuf;
    FastResize( sendBuf, totalSend );
    auto offs = sendOffs;
    for( Int k=0; k<totalSend; ++k )
        sendBuf[offs[owners[k]]++] = remoteEntries[k];
    SwapClear( remoteEntries );
    SwapClear( owners );

    // Exchange and unpack the data
    // ============================
    auto recvBuf = mpi::AllToAll( sendBuf, sendCounts, sendOffs, comm );
    if( BPartic )
    {
        if( B.RedundantRank() == redundantRootB )
        {
            Int recvBufSize = recvBuf.size();
            for( Int k=0; k<recvBufSize; ++k )
            {
                const auto& entry = recvBuf[k];
                BLoc(entry.i,entry.j) = Caster<S,T>::Cast(entry.value);
            }
        }
        El::Broadcast( B, B.RedundantComm(), redundantRootB );
    }
}

template<typename S,typename T,typename>
void GeneralPurpose
( const AbstractDistMatrix<S>& A,
        AbstractDistMatrix<T>& B )
{
    EL_DEBUG_CSE

    if( A.Grid().Size() == 1 && B.Grid().Size() == 1 )
    {
        B.Resize( A.Height(), A.Width() );
        Copy( A.LockedMatrix(), B.Matrix() );
        return;
    }

    Helper( A, B );
}

template<typename T,typename>
void GeneralPurpose
( const AbstractDistMatrix<T>& A,
        AbstractDistMatrix<T>& B )
{
    EL_DEBUG_CSE

    const Int height = A.Height();
    const Int width = A.Width();

    if( A.Grid().Size() == 1 && B.Grid().Size() == 1 )
    {
        B.Resize( height, width );
        Copy( A.LockedMatrix(), B.Matrix() );
        return;
    }

#ifdef EL_HAVE_SCALAPACK
    const bool useBLACSRedist = true;
    if( useBLACSRedist &&
        A.ColDist() == MC && A.RowDist() == MR &&
        B.ColDist() == MC && B.RowDist() == MR &&
        A.ColCut() == 0 && A.RowCut() == 0 &&
        B.ColCut() == 0 && B.RowCut() == 0 )
    {
        B.Resize( height, width );
        const int contextA = blacs::Context( A );
        auto descA = FillDesc( A );
        auto descB = FillDesc( B );

        // This appears to be noticeably faster than the current
        // Elemental-native scheme which also transmits metadata
        //
        // Hmmm...should there be some type of check to ensure that the
        // entire set of processes in contextA encompasses the entire set of
        // processes used for A and B?
        blacs::Redistribute
        ( A.Height(), A.Width(),
          A.LockedBuffer(), descA.data(),
          B.Buffer(),       descB.data(), contextA );

        return;
    }
#endif

    Helper( A, B );
}

} // namespace copy
} // namespace El

#endif // ifndef EL_BLAS_COPY_GENERALPURPOSE_HPP
