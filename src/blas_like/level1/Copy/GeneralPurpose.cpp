/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

#include "El/blas_like/level1/copy_internal.hpp"

namespace El {
namespace copy {

template<typename S,typename T,typename>
void GeneralPurpose
( const AbstractDistMatrix<S>& A,
        AbstractDistMatrix<T>& B ) 
{
    DEBUG_ONLY(CSE cse("copy::GeneralPurpose"))
    // TODO: Decide whether S or T should be used as the transmission type
    //       based upon which is smaller. Transmit S by default.
    const Int height = A.Height();
    const Int width = A.Width();
    const Grid& g = B.Grid();
    const Dist colDist=B.ColDist(), rowDist=B.RowDist();
    const int root = B.Root();
    B.Resize( height, width );

    if( A.Grid().Size() == 1 && B.Grid().Size() == 1 )
    {
        Copy( A.LockedMatrix(), B.Matrix() );
        return;
    }
    const bool includeViewers = (A.Grid() != B.Grid());

    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();

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
            ownerRows[iLoc] = ownerRow;
            const bool isLocalRow = ( ownerRow == colRank );
            globalRows[iLoc] = i;
            localRows[iLoc] =
              ( isLocalRow ? B.LocalRow(i)
                           : B.LocalRow(i,ownerRow) );
        }

        remoteEntries.reserve( localHeight*localWidth );
        distOwners.reserve( localHeight*localWidth );
        const S* ABuf = A.LockedBuffer();
        const Int ALDim = A.LDim();
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const Int j = A.GlobalCol(jLoc);
            const int ownerCol = B.ColOwner(j);
            const bool isLocalCol = ( ownerCol == rowRank );
            const Int localCol = 
              ( isLocalCol ? B.LocalCol(j)
                           : B.LocalCol(j,ownerCol) );
            for( Int iLoc=0; iLoc<localHeight; ++iLoc ) 
            {
                const Int i = globalRows[iLoc];
                const int ownerRow = ownerRows[iLoc];
                const bool isLocalRow = ( ownerRow == colRank );
                const Int localRow = localRows[iLoc];
                const S alpha = ABuf[iLoc+jLoc*ALDim];
                if( noRedundant && isLocalRow && isLocalCol )
                {
                    B.SetLocal( localRow, localCol, Caster<S,T>::Cast(alpha) );
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

    // Compute the metadata
    // ====================
    const Int totalSend = remoteEntries.size();
    mpi::Comm comm;
    vector<int> sendCounts, owners(totalSend);
    if( includeViewers )
    {
        comm = g.ViewingComm();
        const int commSize = mpi::Size( comm );

        vector<int> distMap(commSize);
        for( int q=0; q<commSize; ++q )
        {
            const int vcOwner = g.CoordsToVC(colDist,rowDist,q,root);
            distMap[q] = g.VCToViewing(vcOwner);
        }

        sendCounts.resize(commSize,0);
        for( Int k=0; k<totalSend; ++k )
        {
            owners[k] = distMap[distOwners[k]];
            ++sendCounts[owners[k]];
        }
    }
    else
    {
        if( !B.Participating() )
            return;
        comm = g.VCComm();
        const int commSize = mpi::Size( comm );

        vector<int> distMap(commSize);
        for( int q=0; q<commSize; ++q )
            distMap[q] = g.CoordsToVC(colDist,rowDist,q,root);

        sendCounts.resize(commSize,0);
        for( Int k=0; k<totalSend; ++k )
        {
            owners[k] = distMap[distOwners[k]];
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
    Int recvBufSize = recvBuf.size();
    mpi::Broadcast( recvBufSize, 0, B.RedundantComm() );
    FastResize( recvBuf, recvBufSize );
    mpi::Broadcast( recvBuf.data(), recvBufSize, 0, B.RedundantComm() );
    T* BBuf = B.Buffer();
    const Int BLDim = B.LDim();
    for( Int k=0; k<recvBufSize; ++k )
    {
        const auto& entry = recvBuf[k];
        BBuf[entry.i+entry.j*BLDim] = Caster<S,T>::Cast(entry.value);
    }
}

#define CONVERT(S,T) \
  template void GeneralPurpose \
  ( const AbstractDistMatrix<S>& A, \
          AbstractDistMatrix<T>& B );

#define SAME(T) CONVERT(T,T)

#define PROTO_INT(T) SAME(T) 

#define PROTO_REAL(Real) \
  SAME(Real) \
  CONVERT(Int,Real) \
  CONVERT(Real,Complex<Real>)

#define PROTO_COMPLEX(C) \
  SAME(C) \
  CONVERT(Int,C)

#ifdef EL_HAVE_QUAD

#define PROTO_FLOAT \
  PROTO_REAL(float) \
  CONVERT(float,double) \
  CONVERT(float,Quad) \
  CONVERT(float,Complex<double>) \
  CONVERT(float,Complex<Quad>)

#define PROTO_DOUBLE \
  PROTO_REAL(double) \
  CONVERT(double,float) \
  CONVERT(double,Quad) \
  CONVERT(double,Complex<float>) \
  CONVERT(double,Complex<Quad>)

#define PROTO_QUAD \
  PROTO_REAL(Quad) \
  CONVERT(Quad,float) \
  CONVERT(Quad,double) \
  CONVERT(Quad,Complex<float>) \
  CONVERT(Quad,Complex<double>)

#define PROTO_COMPLEX_FLOAT \
  PROTO_COMPLEX(Complex<float>) \
  CONVERT(Complex<float>,Complex<double>) \
  CONVERT(Complex<float>,Complex<Quad>)

#define PROTO_COMPLEX_DOUBLE \
  PROTO_COMPLEX(Complex<double>) \
  CONVERT(Complex<double>,Complex<float>) \
  CONVERT(Complex<double>,Complex<Quad>)

#define PROTO_COMPLEX_QUAD \
  PROTO_COMPLEX(Complex<Quad>) \
  CONVERT(Complex<Quad>,Complex<float>) \
  CONVERT(Complex<Quad>,Complex<double>)

#else

#define PROTO_FLOAT \
  PROTO_REAL(float) \
  CONVERT(float,double) \
  CONVERT(float,Complex<double>)

#define PROTO_DOUBLE \
  PROTO_REAL(double) \
  CONVERT(double,float) \
  CONVERT(double,Complex<float>)

#define PROTO_COMPLEX_FLOAT \
  PROTO_COMPLEX(Complex<float>) \
  CONVERT(Complex<float>,Complex<double>)

#define PROTO_COMPLEX_DOUBLE \
  PROTO_COMPLEX(Complex<double>) \
  CONVERT(Complex<double>,Complex<float>)

#endif

#define EL_ENABLE_QUAD
#include "El/macros/Instantiate.h"

} // namespace copy
} // namespace El
