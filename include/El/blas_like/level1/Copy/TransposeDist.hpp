/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_BLAS_COPY_TRANSPOSEDIST_HPP
#define EL_BLAS_COPY_TRANSPOSEDIST_HPP

namespace El {
namespace copy {

// TODO: Generalize the below implementation
template<typename T,Dist U,Dist V>
void TransposeDist( const DistMatrix<T,U,V>& A, DistMatrix<T,V,U>& B ) 
{
    DEBUG_ONLY(CSE cse("copy::TransposeDist"))
    AssertSameGrids( A, B );

    const Grid& g = B.Grid();
    B.Resize( A.Height(), A.Width() );
    if( !B.Participating() )
        return;

    const Int colStrideA = A.ColStride();
    const Int rowStrideA = A.RowStride();
    const Int distSize = A.DistSize();

    if( A.DistSize() == 1 && B.DistSize() == 1 ) 
    {
        Copy( A.LockedMatrix(), B.Matrix() );
    }
    else if( A.Width() == 1 )
    {
        const Int height = A.Height();
        const Int maxLocalHeight = MaxLength(height,distSize);
        const Int portionSize = mpi::Pad( maxLocalHeight );

        const Int colDiff = Shift(A.DistRank(),A.ColAlign(),distSize) - 
                            Shift(B.DistRank(),B.ColAlign(),distSize);
        const Int sendRankB = Mod( B.DistRank()+colDiff, distSize );
        const Int recvRankA = Mod( A.DistRank()-colDiff, distSize );
        const Int recvRankB = 
            (recvRankA/colStrideA)+rowStrideA*(recvRankA%colStrideA);

        vector<T> buffer;
        FastResize( buffer, (colStrideA+rowStrideA)*portionSize );
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[colStrideA*portionSize];

        if( A.RowRank() == A.RowAlign() )
        {
            // Pack
            // TODO: Use kernel from copy::util
            const Int AColShift = A.ColShift();
            const T* ABuf = A.LockedBuffer();
            EL_PARALLEL_FOR
            for( Int k=0; k<rowStrideA; ++k )
            {
                T* data = &recvBuf[k*portionSize];

                const Int shift = 
                  Shift_(A.ColRank()+colStrideA*k,A.ColAlign(),distSize);
                const Int offset = (shift-AColShift) / colStrideA;
                const Int thisLocalHeight = Length_(height,shift,distSize);

                for( Int iLoc=0; iLoc<thisLocalHeight; ++iLoc )
                    data[iLoc] = ABuf[offset+iLoc*rowStrideA];
            }
        }

        // (e.g., A[VC,STAR] <- A[MC,MR])
        mpi::Scatter
        ( recvBuf, portionSize,
          sendBuf, portionSize, A.RowAlign(), A.RowComm() );

        // (e.g., A[VR,STAR] <- A[VC,STAR])
        mpi::SendRecv
        ( sendBuf, portionSize, sendRankB,
          recvBuf, portionSize, recvRankB, B.DistComm() );

        // (e.g., A[MR,MC] <- A[VR,STAR])
        mpi::Gather
        ( recvBuf, portionSize,
          sendBuf, portionSize, B.RowAlign(), B.RowComm() );

        if( B.RowRank() == B.RowAlign() )
        {
            // Unpack
            // TODO: Use kernel from copy::util
            T* bufB = B.Buffer();
            EL_PARALLEL_FOR
            for( Int k=0; k<colStrideA; ++k )
            {
                const T* data = &sendBuf[k*portionSize];

                const Int shift = 
                  Shift_(B.ColRank()+rowStrideA*k,B.ColAlign(),distSize);
                const Int offset = (shift-B.ColShift()) / rowStrideA;
                const Int thisLocalHeight = Length_(height,shift,distSize);

                for( Int iLoc=0; iLoc<thisLocalHeight; ++iLoc )
                    bufB[offset+iLoc*colStrideA] = data[iLoc];
            }
        }
    }
    else if( A.Height() == 1 )
    {
        const Int width = A.Width();
        const Int maxLocalWidth = MaxLength(width,distSize);
        const Int portionSize = mpi::Pad( maxLocalWidth );

        const Int rowDiff = Shift(B.DistRank(),A.RowAlign(),distSize) -
                            Shift(A.DistRank(),B.RowAlign(),distSize);
        const Int sendRankA = Mod( A.DistRank()+rowDiff, distSize );
        const Int recvRankB = Mod( B.DistRank()-rowDiff, distSize );
        const Int recvRankA = 
            (recvRankB/rowStrideA)+colStrideA*(recvRankB%rowStrideA);

        vector<T> buffer;
        FastResize( buffer, (colStrideA+rowStrideA)*portionSize );
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[rowStrideA*portionSize];

        if( A.ColRank() == A.ColAlign() )
        {
            // Pack
            // TODO: Use kernel from copy::util
            const T* ABuf = A.LockedBuffer();
            EL_PARALLEL_FOR
            for( Int k=0; k<colStrideA; ++k )
            {
                T* data = &recvBuf[k*portionSize];

                const Int shift = 
                  Shift_(A.RowRank()+rowStrideA*k,A.RowAlign(),distSize);
                const Int offset = (shift-A.RowShift()) / rowStrideA;
                const Int thisLocalWidth = Length_(width,shift,distSize);

                for( Int jLoc=0; jLoc<thisLocalWidth; ++jLoc )
                    data[jLoc] = ABuf[(offset+jLoc*colStrideA)*A.LDim()];
            }
        }

        // (e.g., A[STAR,VR] <- A[MC,MR])
        mpi::Scatter
        ( recvBuf, portionSize,
          sendBuf, portionSize, A.ColAlign(), A.ColComm() );

        // A[STAR,VC] <- A[STAR,VR]
        mpi::SendRecv
        ( sendBuf, portionSize, sendRankA,
          recvBuf, portionSize, recvRankA, A.DistComm() );

        // A[MR,MC] <- A[STAR,VC]
        mpi::Gather
        ( recvBuf, portionSize,
          sendBuf, portionSize, B.ColAlign(), B.ColComm() );

        if( B.ColRank() == B.ColAlign() )
        {
            // Unpack
            // TODO: Use kernel from copy::util
            T* bufB = B.Buffer();
            EL_PARALLEL_FOR
            for( Int k=0; k<rowStrideA; ++k )
            {
                const T* data = &sendBuf[k*portionSize];

                const Int shift = 
                  Shift_(B.RowRank()+colStrideA*k,B.RowAlign(),distSize);
                const Int offset = (shift-B.RowShift()) / colStrideA;
                const Int thisLocalWidth = Length_(width,shift,distSize);

                for( Int jLoc=0; jLoc<thisLocalWidth; ++jLoc )
                    bufB[(offset+jLoc*rowStrideA)*B.LDim()] = data[jLoc];
            }
        }
    }
    else
    {
        if( A.Height() >= A.Width() )
        {
            DistMatrix<T,ProductDist<U,V>(),ProductDistPartner<U,V>()>
              A_ProdDistA( A );
            DistMatrix<T,ProductDist<V,U>(),ProductDistPartner<V,U>()>
              A_ProdDistB( g );
            A_ProdDistB.AlignColsWith( B );
            A_ProdDistB = A_ProdDistA;
            A_ProdDistA.Empty();
            B = A_ProdDistB;
        }
        else
        {
            DistMatrix<T,ProductDistPartner<V,U>(),ProductDist<V,U>()>
                A_ProdDistB( A );
            DistMatrix<T,ProductDistPartner<U,V>(),ProductDist<U,V>()>
                A_ProdDistA( g );

            A_ProdDistA.AlignRowsWith( B );
            A_ProdDistA = A_ProdDistB;
            A_ProdDistB.Empty();
            B = A_ProdDistA;
        }
    }
}

} // namespace copy
} // namespace El

#endif // ifndef EL_BLAS_COPY_TRANSPOSEDIST_HPP
