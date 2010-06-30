/*
   This file is part of Elemental, a library for distributed-memory dense 
   linear algebra.

   Copyright (c) 2009-2010 Jack Poulson <jack.poulson@gmail.com>.
   All rights reserved.

   This file is released under the terms of the license contained in the file
   LICENSE-PURE.
*/
#include "elemental/lapack_internal.hpp"
using namespace std;
using namespace elemental;
using namespace elemental::blas;
using namespace elemental::wrappers::mpi;

template<typename T>
void
elemental::lapack::internal::PanelLU
( DistMatrix<T,Star,Star>& A, 
  DistMatrix<T,VC,  Star>& B, 
  DistMatrix<int,Star,Star>& p, 
  int pivotOffset )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::PanelLU");
    if( A.GetGrid() != p.GetGrid() || p.GetGrid() != B.GetGrid() )
        throw logic_error( "Matrices must be distributed over the same grid." );
    if( A.Width() != B.Width() )
        throw logic_error( "A and B must be the same width." );
    if( A.Height() != p.Height() || p.Width() != 1 )
        throw logic_error( "p must be a vector that conforms with A." );
#endif
    const Grid& g = A.GetGrid();
    const int np = g.Size();
    const int colShift = B.ColShift();
    const int colAlignment = B.ColAlignment();

    // Matrix views
    DistMatrix<T,Star,Star> 
        ATL(g), ATR(g),  A00(g), a01(g),     A02(g),  
        ABL(g), ABR(g),  a10(g), alpha11(g), a12(g),  
                         A20(g), a21(g),     A22(g);

    DistMatrix<T,VC,Star>
        BL(g), BR(g),
        B0(g), b1(g), B2(g);

    DistMatrix<int,Star,Star>
        pT(g),  p0(g),
        pB(g),  psi1(g),
                p2(g);

    const int width = A.Width();
    const int numBytes = (width+1)*sizeof(T)+sizeof(int);
    vector<char> sendData(numBytes);
    vector<char> recvData(numBytes);

    // Extract pointers to send and recv data
    char* sendBuf = &sendData[0];
    char* recvBuf = &recvData[0];

    // Start the algorithm
    PushBlocksizeStack( 1 );
    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    PartitionRight( B, BL, BR, 0 );
    PartitionDown
    ( p, pT,
         pB, 0 );
    while( ATL.Height() < A.Height() )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ a01,     A02,
         /*************/ /**********************/
               /**/       a10, /**/ alpha11, a12,
          ABL, /**/ ABR,  A20, /**/ a21,     A22 );

        RepartitionRight
        ( BL, /**/ BR,  
          B0, /**/ b1, B2 );

        RepartitionDown
        ( pT,  p0,
         /**/ /****/
               psi1,
          pB,  p2 );

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
                pivotIndex = A.Height() + colShift + i*np;
            }
        }

        // Fill the send buffer with:
        // [ pivotValue | pivotRow | pivotIndex ]
        if( pivotIndex < A.Height() )
        {
            ((T*)sendBuf)[0] = A.LocalEntry(pivotIndex,a10.Width());
            for( int j=0; j<width; ++j )
                ((T*)sendBuf)[j+1] = A.LocalEntry(pivotIndex,j);
        }
        else
        {
            const int localIndex = ((pivotIndex-A.Height())-colShift)/np;
            ((T*)sendBuf)[0] = b1.LocalEntry(localIndex,0);
            for( int j=0; j<width; ++j )
                ((T*)sendBuf)[j+1] = B.LocalEntry(localIndex,j);
        }
        ((int*)(((T*)sendBuf)+width+1))[0] = pivotIndex;

        // Communicate to establish the pivot information
        AllReduce
        ( sendBuf, recvBuf, numBytes, PivotOp<T>(), g.VCComm() );

        // Update the pivot vector
        const int maxIndex = ((int*)(((T*)recvBuf)+width+1))[0];
        p.LocalEntry(a01.Height(),0) = maxIndex + pivotOffset;

        // Copy the current row into the pivot row
        if( maxIndex < A.Height() )
        {
            for( int j=0; j<width; ++j )
                A.LocalEntry(maxIndex,j) = A.LocalEntry(A00.Height(),j);
        }
        else
        {
            const int ownerRank = (colAlignment+(maxIndex-A.Height())) % np;
            if( g.VCRank() == ownerRank )
            {
                const int localIndex = ((maxIndex-A.Height())-colShift) / np;
                for( int j=0; j<width; ++j )
                    B.LocalEntry(localIndex,j) = A.LocalEntry(A00.Height(),j);
            }
        }

        // Copy the pivot row into the current row
        for( int j=0; j<width; ++j )
            A.LocalEntry(A00.Height(),j) = ((T*)recvBuf)[j+1];

        // Now we can perform the update
        T alpha11Inv = ((T)1) / alpha11.LocalEntry(0,0);
        Scal( alpha11Inv, a21.LocalMatrix() );
        Scal( alpha11Inv, b1.LocalMatrix()  );
        Geru( (T)-1, a21.LocalMatrix(), a12.LocalMatrix(), A22.LocalMatrix() );
        Geru( (T)-1, b1.LocalMatrix(), a12.LocalMatrix(), B2.LocalMatrix() );

        //--------------------------------------------------------------------//

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, a01,     /**/ A02,
               /**/       a10, alpha11, /**/ a12,
         /*************/ /**********************/
          ABL, /**/ ABR,  A20, a21,     /**/ A22 );

        SlidePartitionRight
        ( BL,     /**/ BR,  
          B0, b1, /**/ B2 );

        SlidePartitionDown
        ( pT,  p0,
               psi1,
         /**/ /****/
          pB,  p2 );
    }
    PopBlocksizeStack();

#ifndef RELEASE
    PopCallStack();
#endif
}

template void
elemental::lapack::internal::PanelLU
( DistMatrix<float,Star,Star>& A, 
  DistMatrix<float,VC,  Star>& B, 
  DistMatrix<int,  Star,Star>& p,
  int pivotOffset );

template void
elemental::lapack::internal::PanelLU
( DistMatrix<double,Star,Star>& A, 
  DistMatrix<double,VC,  Star>& B, 
  DistMatrix<int,   Star,Star>& p,
  int pivotOffset );

#ifndef WITHOUT_COMPLEX
template void
elemental::lapack::internal::PanelLU
( DistMatrix<scomplex,Star,Star>& A,
  DistMatrix<scomplex,VC,  Star>& B,
  DistMatrix<int,     Star,Star>& p,
  int pivotOffset );

template void
elemental::lapack::internal::PanelLU
( DistMatrix<dcomplex,Star,Star>& A,
  DistMatrix<dcomplex,VC,  Star>& B,
  DistMatrix<int,     Star,Star>& p,
  int pivotOffset );
#endif

