/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_GEMM_NN_HPP
#define ELEM_GEMM_NN_HPP

#include ELEM_SCALE_INC

namespace elem {
namespace gemm {

// Cannon's algorithm
template<typename T>
inline void
Cannon_NN
( T alpha, const DistMatrix<T>& A,
           const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C )
{
    DEBUG_ONLY(
        CallStackEntry cse("gemm::Cannon_NN");
        if( A.Grid() != B.Grid() || B.Grid() != C.Grid() )
            LogicError("{A,B,C} must have the same grid");
        if( A.Height() != C.Height() || B.Width() != C.Width() ||
            A.Width() != B.Height() )
            LogicError
            ("Nonconformal matrices:\n",
             DimsString(A,"A"),"\n",DimsString(B,"B"),"\n",DimsString(C,"C"));
    )
    const Grid& g = A.Grid();
    if( g.Height() != g.Width() )
        LogicError("Process grid must be square for Cannon's");
    if( C.ColAlign() != A.ColAlign() || C.RowAlign() != B.RowAlign() )
        LogicError("C is not properly aligned");

    const Int row = g.Row();
    const Int col = g.Col();
    const Int pSqrt = g.Height();
    mpi::Comm rowComm = g.RowComm();
    mpi::Comm colComm = g.ColComm(); 
    if( A.Width() % pSqrt != 0 )
        LogicError("For now, width(A) must be integer multiple of sqrt(p)");

    // Begin by scaling our local portion of C
    Scale( beta, C );

    // Load the initial A and B packages (may want to transpose B...)
    const Int localHeightA = A.LocalHeight();
    const Int localHeightB = B.LocalHeight();
    const Int localWidthA = A.LocalWidth();
    const Int localWidthB = B.LocalWidth();
    Matrix<T> pkgA(localHeightA,localWidthA,localHeightA), 
              pkgB(localHeightB,localWidthB,localHeightB);
    for( Int jLoc=0; jLoc<localWidthA; ++jLoc )
        MemCopy
        ( pkgA.Buffer(0,jLoc), A.LockedBuffer(0,jLoc), localHeightA );
    for( Int jLoc=0; jLoc<localWidthB; ++jLoc )
        MemCopy
        ( pkgB.Buffer(0,jLoc), B.LockedBuffer(0,jLoc), localHeightB );

    // Perform the initial circular shifts so that our A and B packages align
    const Int rowShiftA = A.RowShift();
    const Int colShiftB = B.ColShift();
    const Int leftInitA  = Mod(col-colShiftB,pSqrt);
    const Int rightInitA = Mod(col+colShiftB,pSqrt);
    const Int aboveInitB = Mod(row-rowShiftA,pSqrt);
    const Int belowInitB = Mod(row+rowShiftA,pSqrt);
    const Int pkgSizeA = localHeightA*localWidthA;
    const Int pkgSizeB = localHeightB*localWidthB;
    mpi::SendRecv( pkgA.Buffer(), pkgSizeA, leftInitA, rightInitA, rowComm );
    mpi::SendRecv( pkgB.Buffer(), pkgSizeB, aboveInitB, belowInitB, colComm );

    // Now begin the data flow
    const Int aboveRow = Mod(row-1,pSqrt);
    const Int belowRow = Mod(row+1,pSqrt);
    const Int leftCol  = Mod(col-1,pSqrt);
    const Int rightCol = Mod(col+1,pSqrt);
    for( Int q=0; q<pSqrt; ++q )
    {
        Gemm( NORMAL, NORMAL, alpha, pkgA, pkgB, T(1), C.Matrix() );
        if( q != pSqrt-1 )
        {
            mpi::SendRecv
            ( pkgA.Buffer(), pkgSizeA, leftCol, rightCol, rowComm );
            mpi::SendRecv
            ( pkgB.Buffer(), pkgSizeB, aboveRow, belowRow, colComm );
        }
    }
}

// Normal Normal Gemm that avoids communicating the matrix A
template<typename T>
inline void
SUMMA_NNA
( T alpha, const DistMatrix<T>& A,
           const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C )
{
    DEBUG_ONLY(
        CallStackEntry cse("gemm::SUMMA_NNA");
        if( A.Grid() != B.Grid() || B.Grid() != C.Grid() )
            LogicError("{A,B,C} must have the same grid");
        if( A.Height() != C.Height() || B.Width() != C.Width() ||
            A.Width() != B.Height() )
            LogicError
            ("Nonconformal matrices:\n",
             DimsString(A,"A"),"\n",DimsString(B,"B"),"\n",DimsString(C,"C"));
    )
    const Grid& g = A.Grid();

    // Matrix views
    DistMatrix<T> BL(g), BR(g),
                  B0(g), B1(g), B2(g);
    DistMatrix<T> CL(g), CR(g),
                  C0(g), C1(g), C2(g);

    // Temporary distributions
    DistMatrix<T,VR,STAR> B1_VR_STAR(g);
    DistMatrix<T,STAR,MR> B1Trans_STAR_MR(g);
    DistMatrix<T,MC,STAR> D1_MC_STAR(g);

    B1_VR_STAR.AlignWith( A );
    B1Trans_STAR_MR.AlignWith( A );
    D1_MC_STAR.AlignWith( A );

    // Start the algorithm
    Scale( beta, C );
    LockedPartitionRight( B, BL, BR, 0 );
    PartitionRight( C, CL, CR, 0 );
    while( BR.Width() > 0 )
    {
        LockedRepartitionRight
        ( BL, /**/     BR,
          B0, /**/ B1, B2 );

        RepartitionRight
        ( CL, /**/     CR,
          C0, /**/ C1, C2 );

        //--------------------------------------------------------------------//
        B1_VR_STAR = B1;
        B1_VR_STAR.TransposePartialColAllGather( B1Trans_STAR_MR );

        // D1[MC,*] := alpha A[MC,MR] B1[MR,*]
        LocalGemm( NORMAL, TRANSPOSE, alpha, A, B1Trans_STAR_MR, D1_MC_STAR );

        // C1[MC,MR] += scattered result of D1[MC,*] summed over grid rows
        C1.RowSumScatterUpdate( T(1), D1_MC_STAR );
        //--------------------------------------------------------------------//

        SlideLockedPartitionRight
        ( BL,     /**/ BR,
          B0, B1, /**/ B2 );

        SlidePartitionRight
        ( CL,     /**/ CR,
          C0, C1, /**/ C2 );
    }
}

// Normal Normal Gemm that avoids communicating the matrix B
template<typename T>
inline void 
SUMMA_NNB
( T alpha, const DistMatrix<T>& A,
           const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C )
{
    DEBUG_ONLY(
        CallStackEntry cse("gemm::SUMMA_NNB");
        if( A.Grid() != B.Grid() || B.Grid() != C.Grid() )
            LogicError("{A,B,C} must be distributed over the same grid");
        if( A.Height() != C.Height() || B.Width() != C.Width() ||
            A.Width() != B.Height() )
            LogicError
            ("Nonconformal matrices:\n",
             DimsString(A,"A"),"\n",DimsString(B,"B"),"\n",DimsString(C,"C"));
    )
    const Grid& g = A.Grid();

    // Matrix views
    DistMatrix<T> AT(g),  A0(g),
                  AB(g),  A1(g),
                          A2(g);
    DistMatrix<T> CT(g),  C0(g),
                  CB(g),  C1(g),
                          C2(g);

    // Temporary distributions
    DistMatrix<T,STAR,MC> A1_STAR_MC(g);
    DistMatrix<T,MR,STAR> D1Trans_MR_STAR(g);

    A1_STAR_MC.AlignWith( B );
    D1Trans_MR_STAR.AlignWith( B );

    // Start the algorithm
    Scale( beta, C );
    LockedPartitionDown
    ( A, AT,
         AB, 0 );
    PartitionDown
    ( C, CT,
         CB, 0 );
    while( AB.Height() > 0 )
    {
        LockedRepartitionDown
        ( AT,  A0,
         /**/ /**/
               A1,
          AB,  A2 );

        RepartitionDown
        ( CT,  C0,
         /**/ /**/
               C1,
          CB,  C2 );

        //--------------------------------------------------------------------//
        A1_STAR_MC = A1; // A1[*,MC] <- A1[MC,MR]

        // D1^T[MR,* ] := alpha B^T[MR,MC] A1^T[MC,* ]
        LocalGemm
        ( TRANSPOSE, TRANSPOSE, alpha, B, A1_STAR_MC, D1Trans_MR_STAR );

        C1.TransposeColSumScatterUpdate( T(1), D1Trans_MR_STAR );
        //--------------------------------------------------------------------//

        SlideLockedPartitionDown
        ( AT,  A0,
               A1,
         /**/ /**/
          AB,  A2 );
 
        SlidePartitionDown
        ( CT,  C0,
               C1,
         /**/ /**/
          CB,  C2 );
    }
}                     

// Normal Normal Gemm that avoids communicating the matrix C
template<typename T>
inline void 
SUMMA_NNC
( T alpha, const DistMatrix<T>& A,
           const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C )
{
    DEBUG_ONLY(
        CallStackEntry cse("gemm::SUMMA_NNC");
        if( A.Grid() != B.Grid() || B.Grid() != C.Grid() )
            LogicError("{A,B,C} must be distributed over the same grid");
        if( A.Height() != C.Height() || B.Width() != C.Width() ||
            A.Width() != B.Height() )
            LogicError
            ("Nonconformal matrices:\n",
             DimsString(A,"A"),"\n",DimsString(B,"B"),"\n",DimsString(C,"C"));
    )
    const Grid& g = A.Grid();

    // Matrix views
    DistMatrix<T> AL(g), AR(g),
                  A0(g), A1(g), A2(g);         
    DistMatrix<T> BT(g),  B0(g),
                  BB(g),  B1(g),
                          B2(g);

    // Temporary distributions
    DistMatrix<T,MC,STAR> A1_MC_STAR(g);
    DistMatrix<T,MR,STAR> B1Trans_MR_STAR(g); 

    A1_MC_STAR.AlignWith( C );
    B1Trans_MR_STAR.AlignWith( C );

    // Start the algorithm
    Scale( beta, C );
    LockedPartitionRight( A, AL, AR, 0 ); 
    LockedPartitionDown
    ( B, BT, 
         BB, 0 ); 
    while( AR.Width() > 0 )
    {
        LockedRepartitionRight( AL, /**/ AR,
                                A0, /**/ A1, A2 );

        LockedRepartitionDown( BT,  B0,
                              /**/ /**/
                                    B1, 
                               BB,  B2 );

        //--------------------------------------------------------------------//
        A1_MC_STAR = A1; 
        B1.TransposeColAllGather( B1Trans_MR_STAR );

        // C[MC,MR] += alpha A1[MC,*] (B1^T[MR,*])^T
        //           = alpha A1[MC,*] B1[*,MR]
        LocalGemm
        ( NORMAL, TRANSPOSE, alpha, A1_MC_STAR, B1Trans_MR_STAR, T(1), C );
        //--------------------------------------------------------------------//

        SlideLockedPartitionRight( AL,     /**/ AR,
                                   A0, A1, /**/ A2 );

        SlideLockedPartitionDown( BT,  B0,
                                       B1,
                                 /**/ /**/
                                  BB,  B2 );
    }
}

// Normal Normal Gemm for panel-panel dot products
template<typename T>
inline void 
SUMMA_NNDot
( T alpha, const DistMatrix<T>& A,
           const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C )
{
    DEBUG_ONLY(
        CallStackEntry cse("gemm::SUMMA_NNDot");
        if( A.Grid() != B.Grid() || B.Grid() != C.Grid() )
            LogicError("{A,B,C} must have the same grid");
        if( A.Height() != C.Height() || B.Width() != C.Width() ||
            A.Width() != B.Height() )
            LogicError
            ("Nonconformal matrices:\n",
             DimsString(A,"A"),"\n",DimsString(B,"B"),"\n",DimsString(C,"C"));
    )
    const Grid& g = A.Grid();

    if( A.Height() > B.Width() )
    {
        // Matrix views
        DistMatrix<T> AT(g), AB(g),
                      A0(g), A1(g), A2(g);         
        DistMatrix<T> BL(g),  B0(g),
                      BR(g),  B1(g),
                              B2(g);
        DistMatrix<T> CT(g), C0(g), C1L(g), C1R(g),
                      CB(g), C1(g), C10(g), C11(g), C12(g),
                             C2(g);

        // Temporary distributions
        DistMatrix<T,STAR,VC> A1_STAR_VC(g);
        DistMatrix<T,VC,STAR> B1_VC_STAR(g);
        DistMatrix<T,STAR,STAR> C11_STAR_STAR(g);

        // Star the algorithm
        Scale( beta, C );
        LockedPartitionDown
        ( A, AT,
             AB, 0 );
        PartitionDown
        ( C, CT,
             CB, 0 );
        while( AB.Height() > 0 )
        {
            LockedRepartitionDown
            ( AT,  A0,
             /**/ /**/
                   A1,
              AB,  A2 );

            RepartitionDown
            ( CT,  C0,
             /**/ /**/
                   C1,
              CB,  C2 );

            A1_STAR_VC = A1; 
            B1_VC_STAR.AlignWith( A1_STAR_VC );

            LockedPartitionRight( B, BL, BR, 0 );
            PartitionRight( C1, C1L, C1R, 0 );
            while( BR.Width() > 0 )
            {
                LockedRepartitionRight
                ( BL, /**/ BR,
                  B0, /**/ B1, B2 );

                RepartitionRight
                ( C1L, /**/ C1R,
                  C10, /**/ C11, C12 );

                //------------------------------------------------------------//
                B1_VC_STAR = B1;
                LocalGemm
                ( NORMAL, NORMAL, 
                  alpha, A1_STAR_VC, B1_VC_STAR, C11_STAR_STAR );
                C11.SumScatterUpdate( T(1), C11_STAR_STAR );
                //------------------------------------------------------------//

                SlideLockedPartitionRight
                ( BL,     /**/ BR,
                  B0, B1, /**/ B2 );

                SlidePartitionRight
                ( C1L,      /**/ C1R,
                  C10, C11, /**/ C12 );
            }

            SlideLockedPartitionDown
            ( AT,  A0,
                   A1,
             /**/ /**/
              AB,  A2 );

            SlidePartitionDown
            ( CT,  C0,
                   C1,
             /**/ /**/
              CB,  C2 );
        }
    }
    else
    {
        // Matrix views
        DistMatrix<T> AT(g), AB(g),
                      A0(g), A1(g), A2(g);         
        DistMatrix<T> BL(g),  B0(g),
                      BR(g),  B1(g),
                              B2(g);
        DistMatrix<T> 
            CL(g), CR(g),         C1T(g),  C01(g),
            C0(g), C1(g), C2(g),  C1B(g),  C11(g),
                                           C21(g);

        // Temporary distributions
        DistMatrix<T,STAR,VR> A1_STAR_VR(g);
        DistMatrix<T,VR,STAR> B1_VR_STAR(g);
        DistMatrix<T,STAR,STAR> C11_STAR_STAR(g);

        // Star the algorithm
        Scale( beta, C );
        LockedPartitionRight( B, BL, BR, 0 );
        PartitionRight( C, CL, CR, 0 );
        while( BR.Width() > 0 )
        {
            LockedRepartitionRight
            ( BL, /**/ BR,
              B0, /**/ B1, B2 );

            RepartitionRight
            ( CL, /**/ CR,
              C0, /**/ C1, C2 );

            B1_VR_STAR = B1;
            A1_STAR_VR.AlignWith( B1_VR_STAR );

            LockedPartitionDown
            ( A, AT,
                 AB, 0 );
            PartitionDown
            ( C1, C1T,
                  C1B, 0 );
            while( AB.Height() > 0 )
            {
                LockedRepartitionDown
                ( AT,  A0,
                 /**/ /**/
                       A1,
                  AB,  A2 );

                RepartitionDown
                ( C1T,  C01,
                 /***/ /***/
                        C11,
                  C1B,  C21 );

                //------------------------------------------------------------//
                A1_STAR_VR = A1;
                LocalGemm
                ( NORMAL, NORMAL, 
                  alpha, A1_STAR_VR, B1_VR_STAR, C11_STAR_STAR );
                C11.SumScatterUpdate( T(1), C11_STAR_STAR );
                //------------------------------------------------------------//

                SlideLockedPartitionDown
                ( AT,  A0,
                       A1,
                 /**/ /**/
                  AB,  A2 );

                SlidePartitionDown
                ( C1T,  C01,
                        C11,
                 /***/ /***/
                  C1B,  C21 );
            }

            SlideLockedPartitionRight
            ( BL,     /**/ BR,
              B0, B1, /**/ B2 ); 

            SlidePartitionRight
            ( CL,     /**/ CR,
              C0, C1, /**/ C2 );
        }
    }
}

template<typename T>
inline void
SUMMA_NN
( T alpha, const DistMatrix<T>& A,
           const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C )
{
    DEBUG_ONLY(CallStackEntry cse("gemm::SUMMA_NN"))
    const Int m = C.Height();
    const Int n = C.Width();
    const Int k = A.Width();
    const double weightTowardsC = 2.;
    const double weightAwayFromDot = 10.;

    if( weightAwayFromDot*m <= k && weightAwayFromDot*n <= k )
        SUMMA_NNDot( alpha, A, B, beta, C );
    else if( m <= n && weightTowardsC*m <= k )
        SUMMA_NNB( alpha, A, B, beta, C );    
    else if( n <= m && weightTowardsC*n <= k )
        SUMMA_NNA( alpha, A, B, beta, C );
    else
        SUMMA_NNC( alpha, A, B, beta, C );
}

} // namespace gemm
} // namespace elem

#endif // ifndef ELEM_GEMM_NN_HPP
