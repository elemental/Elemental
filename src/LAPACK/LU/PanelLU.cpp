/*
   Copyright 2009-2010 Jack Poulson

   This file is part of Elemental.

   Elemental is free software: you can redistribute it and/or modify it under
   the terms of the GNU Lesser General Public License as published by the
   Free Software Foundation; either version 3 of the License, or 
   (at your option) any later version.

   Elemental is distributed in the hope that it will be useful, but 
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with Elemental. If not, see <http://www.gnu.org/licenses/>.
*/
#include "Elemental/LAPACKInternal.hpp"
using namespace std;
using namespace Elemental;
using namespace Elemental::BLAS;
using namespace Elemental::wrappers::MPI;

template<typename T>
void
Elemental::LAPACK::Internal::PanelLU
( DistMatrix<T,Star,Star>& A, 
  DistMatrix<T,VC,  Star>& B, 
  DistMatrix<int,Star,Star>& p, 
  const int pivotOffset        )
{
#ifndef RELEASE
    PushCallStack("LAPACK::Internal::PanelLU");
    if( A.GetGrid() != p.GetGrid() || p.GetGrid() != B.GetGrid() )
        throw "Matrices must be distributed over the same grid.";
    if( A.Width() != B.Width() )
        throw "A and B must be the same width.";
    if( A.Height() != p.Height() || p.Width() != 1 )
        throw "p must be a vector that conforms with A.";
#endif
    const Grid& grid = A.GetGrid();
    const int colShift = B.ColShift();

    // Matrix views
    DistMatrix<T,Star,Star> 
        ATL(grid), ATR(grid),  A00(grid), a01(grid),     A02(grid),  
        ABL(grid), ABR(grid),  a10(grid), alpha11(grid), a12(grid),  
                               A20(grid), a21(grid),     A22(grid);

    DistMatrix<T,VC,Star>
        BL(grid), BR(grid),
        B0(grid), b1(grid), B2(grid);

    DistMatrix<int,Star,Star>
        pT(grid),  p0(grid),
        pB(grid),  psi1(grid),
                   p2(grid);

    Matrix<T> pivotCol, pivotRow, a1;

    const int width = A.Width();
    const int numBytes = (width+1)*sizeof(T)+sizeof(int);
    vector<char> sendData(numBytes);
    vector<char> recvData(numBytes);

    // Extract pointers to send and recv data
    char* sendBuf = &sendData[0];
    char* recvBuf = &recvData[0];

    // Start the algorithm
    PushBlocksizeStack( 1 );
    PartitionDownDiagonal( A, ATL, ATR,
                              ABL, ABR );
    PartitionRight( B, BL, BR );
    PartitionDown( p, pT,
                      pB );
    while( ATL.Height() < A.Height() )
    {
        RepartitionDownDiagonal( ATL, /**/ ATR,  A00, /**/ a01,     A02,
                                /*************/ /**********************/
                                      /**/       a10, /**/ alpha11, a12,
                                 ABL, /**/ ABR,  A20, /**/ a21,     A22 );

        RepartitionRight( BL, /**/ BR,  
                          B0, /**/ b1, B2 );

        RepartitionDown( pT,  p0,
                        /**/ /****/
                              psi1,
                         pB,  p2   );

        //--------------------------------------------------------------------//
        
        // Store the index/value of the pivot candidate in A
        T pivotValue = alpha11.LocalEntry(0,0);
        int pivotIndex = a01.Height();
        for( int i=0; i<a21.Height(); ++i )
        {
            T value = a21.LocalEntry(i,0);
            if( FastAbs(value) > FastAbs(pivotValue) )
            {
                pivotValue = value;
                pivotIndex = a01.Height() + i + 1;
            }
        }

        // Update the pivot candidate to include local data from B
        for( int i=0; i<B.LocalHeight(); ++i )
        {
            T value = b1.LocalEntry(i,0);
            if( FastAbs(value) > FastAbs(pivotValue) )
            {
                pivotValue = value;
                pivotIndex = A.Height() + colShift + i*grid.Size();
            }
        }

        // Fill the send buffer with:
        // [ pivotValue | pivotRow | pivotIndex ]
        if( pivotIndex < A.Height() )
        {
            pivotCol.View( A.LocalMatrix(), 0, A00.Width(), A.Height(), 1 );
            pivotRow.View( A.LocalMatrix(), pivotIndex, 0, 1, A.Width() );
            ((T*)sendBuf)[0] = pivotCol(pivotIndex,0); 
            for( int j=0; j<width; ++j )
                ((T*)sendBuf)[j+1] = pivotRow(0,j);
        }
        else
        {
            const int localIndex = 
                ((pivotIndex-A.Height())-colShift)/grid.Size();

            pivotCol.View( b1.LocalMatrix() );
            pivotRow.View( B.LocalMatrix(), localIndex, 0, 1, width );
            ((T*)sendBuf)[0] = pivotCol(localIndex,0);
            for( int j=0; j<width; ++j )
                ((T*)sendBuf)[j+1] = pivotRow(0,j);
        }
        ((int*)(((T*)sendBuf)+width+1))[0] = pivotIndex;

        // Communicate to establish the pivot information
        AllReduce
        ( sendBuf, recvBuf, numBytes, PivotOp<T>(), grid.VCComm() );

        // Update the pivot vector
        const int maxIndex = ((int*)(((T*)recvBuf)+width+1))[0];
        p.LocalEntry(a01.Height(),0) = maxIndex + pivotOffset;

        // Copy a1 into the pivot row
        a1.View( A.LocalMatrix(), A00.Height(), 0, 1, width );
        if( maxIndex < A.Height() )
        {
            pivotRow.View( A.LocalMatrix(), maxIndex, 0, 1, width );
            for( int j=0; j<width; ++j )
                pivotRow(0,j) = a1(0,j);
        }
        else
        {
            const int ownerRank = 
                (B.ColAlignment()+(maxIndex-A.Height())) % grid.Size();
            if( grid.VCRank() == ownerRank )
            {
                const int localIndex = 
                    ((maxIndex-A.Height())-colShift) / grid.Size();
                pivotRow.View( B.LocalMatrix(), localIndex, 0, 1, width );
                for( int j=0; j<width; ++j )
                    pivotRow(0,j) = a1(0,j);
            }
        }

        // Copy the pivot row into a1
        for( int j=0; j<width; ++j )
            a1(0,j) = ((T*)recvBuf)[j+1];

        // Now we can perform the update
        T alpha11Inv = ((T)1) / alpha11.LocalEntry(0,0);
        Scal( alpha11Inv, a21.LocalMatrix() );
        Scal( alpha11Inv, b1.LocalMatrix()  );
        Geru( (T)-1, a21.LocalMatrix(), a12.LocalMatrix(), A22.LocalMatrix() );
        Geru( (T)-1, b1.LocalMatrix(), a12.LocalMatrix(), B2.LocalMatrix() );

        //--------------------------------------------------------------------//

        SlidePartitionDownDiagonal( ATL, /**/ ATR,  A00, a01,     /**/ A02,
                                         /**/       a10, alpha11, /**/ a12,
                                   /*************/ /**********************/
                                    ABL, /**/ ABR,  A20, a21,     /**/ A22 );

        SlidePartitionRight( BL,     /**/ BR,  
                             B0, b1, /**/ B2 );

        SlidePartitionDown( pT,  p0,
                                 psi1,
                           /**/ /****/
                            pB,  p2   );
    }
    PopBlocksizeStack();

#ifndef RELEASE
    PopCallStack();
#endif
}

template void
Elemental::LAPACK::Internal::PanelLU
( DistMatrix<float,Star,Star>& A, 
  DistMatrix<float,VC,  Star>& B, 
  DistMatrix<int,  Star,Star>& p,
  const int pivotOffset          );

template void
Elemental::LAPACK::Internal::PanelLU
( DistMatrix<double,Star,Star>& A, 
  DistMatrix<double,VC,  Star>& B, 
  DistMatrix<int,   Star,Star>& p,
  const int pivotOffset           );

#ifndef WITHOUT_COMPLEX
template void
Elemental::LAPACK::Internal::PanelLU
( DistMatrix<scomplex,Star,Star>& A,
  DistMatrix<scomplex,VC,  Star>& B,
  DistMatrix<int,     Star,Star>& p,
  const int pivotOffset             );

template void
Elemental::LAPACK::Internal::PanelLU
( DistMatrix<dcomplex,Star,Star>& A,
  DistMatrix<dcomplex,VC,  Star>& B,
  DistMatrix<int,     Star,Star>& p,
  const int pivotOffset             );
#endif

