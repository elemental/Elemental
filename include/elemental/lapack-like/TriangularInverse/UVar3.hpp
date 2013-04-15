/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_TRIANGULARINVERSE_UVAR3_HPP
#define LAPACK_TRIANGULARINVERSE_UVAR3_HPP

#include "elemental/blas-like/level3/Gemm.hpp"
#include "elemental/blas-like/level3/Trsm.hpp"

namespace elem {
namespace triangular_inverse {

template<typename F>
inline void
UVar3Unb( UnitOrNonUnit diag, Matrix<F>& U )
{
#ifndef RELEASE
    PushCallStack("triangular_inverse::UVar3Unb");
    if( U.Height() != U.Width() )
        throw std::logic_error("Nonsquare matrices cannot be triangular");
#endif
    const int n = U.Height();
    const int ldu = U.LDim();
    F* UBuffer = U.Buffer();
    for( int j=n-1; j>=0; --j )
    {
        const F upsilon = ( diag==NON_UNIT ? UBuffer[j+j*ldu] : F(1) );
        for( int k=0; k<j; ++k )
            UBuffer[k+j*ldu] /= -upsilon;
        blas::Geru
        ( j, n-(j+1), F(1),
          &UBuffer[j*ldu], 1, &UBuffer[j+(j+1)*ldu], ldu, 
          &UBuffer[(j+1)*ldu], ldu );
        if( diag == NON_UNIT )
        {
            for( int k=j+1; k<n; ++k )
                UBuffer[j+k*ldu] /= upsilon;
            UBuffer[j+j*ldu] = F(1) / UBuffer[j+j*ldu];
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
UVar3( UnitOrNonUnit diag, Matrix<F>& U )
{
#ifndef RELEASE
    PushCallStack("triangular_inverse::UVar3");
    if( U.Height() != U.Width() )
        throw std::logic_error("Nonsquare matrices cannot be triangular");
#endif
    // Matrix views
    Matrix<F> 
        UTL, UTR,  U00, U01, U02,
        UBL, UBR,  U10, U11, U12,
                   U20, U21, U22;

    // Start the algorithm
    PartitionUpDiagonal
    ( U, UTL, UTR,
         UBL, UBR, 0 );
    while( UBR.Height() < U.Height() )
    {
        RepartitionUpDiagonal
        ( UTL, /**/ UTR,  U00, U01, /**/ U02,
               /**/       U10, U11, /**/ U12,
         /*************/ /******************/
          UBL, /**/ UBR,  U20, U21, /**/ U22 );

        //--------------------------------------------------------------------//
        Trsm( RIGHT, UPPER, NORMAL, diag, F(-1), U11, U01 );
        Gemm( NORMAL, NORMAL, F(1), U01, U12, F(1), U02 );
        Trsm( LEFT, UPPER, NORMAL, diag, F(1), U11, U12 );
        UVar3Unb( diag, U11 );
        //--------------------------------------------------------------------//

        SlidePartitionUpDiagonal
        ( UTL, /**/ UTR,  U00, /**/ U01, U02,
         /*************/ /******************/
               /**/       U10, /**/ U11, U12,
          UBL, /**/ UBR,  U20, /**/ U21, U22 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
UVar3( UnitOrNonUnit diag, DistMatrix<F>& U )
{
#ifndef RELEASE
    PushCallStack("triangular_inverse::UVar3");
    if( U.Height() != U.Width() )
        throw std::logic_error("Nonsquare matrices cannot be triangular");
#endif
    const Grid& g = U.Grid();

    // Matrix views
    DistMatrix<F> 
        UTL(g), UTR(g),  U00(g), U01(g), U02(g),
        UBL(g), UBR(g),  U10(g), U11(g), U12(g),
                         U20(g), U21(g), U22(g);

    // Temporary distributions
    DistMatrix<F,VC,  STAR> U01_VC_STAR(g);
    DistMatrix<F,STAR,STAR> U11_STAR_STAR(g);
    DistMatrix<F,STAR,VR  > U12_STAR_VR(g);
    DistMatrix<F,STAR,MC  > U01Trans_STAR_MC(g);
    DistMatrix<F,MR,  STAR> U12Trans_MR_STAR(g);

    // Start the algorithm
    PartitionUpDiagonal
    ( U, UTL, UTR,
         UBL, UBR, 0 );
    while( UBR.Height() < U.Height() )
    {
        RepartitionUpDiagonal
        ( UTL, /**/ UTR,  U00, U01, /**/ U02,
               /**/       U10, U11, /**/ U12,
         /*************/ /******************/
          UBL, /**/ UBR,  U20, U21, /**/ U22 );

        U01Trans_STAR_MC.AlignWith( U02 );
        U12Trans_MR_STAR.AlignWith( U02 );
        //--------------------------------------------------------------------//
        U01_VC_STAR = U01;
        U11_STAR_STAR = U11;
        LocalTrsm
        ( RIGHT, UPPER, NORMAL, diag, F(-1), U11_STAR_STAR, U01_VC_STAR );

        // We transpose before the communication to avoid cache-thrashing
        // in the unpacking stage.
        U12Trans_MR_STAR.TransposeFrom( U12 );
        U01Trans_STAR_MC.TransposeFrom( U01_VC_STAR );

        LocalGemm
        ( TRANSPOSE, TRANSPOSE, 
          F(1), U01Trans_STAR_MC, U12Trans_MR_STAR, F(1), U02 );
        U01.TransposeFrom( U01Trans_STAR_MC );

        U12_STAR_VR.TransposeFrom( U12Trans_MR_STAR );
        LocalTrsm
        ( LEFT, UPPER, NORMAL, diag, F(1), U11_STAR_STAR, U12_STAR_VR );
        LocalTriangularInverse( UPPER, diag, U11_STAR_STAR );
        U11 = U11_STAR_STAR;
        U12 = U12_STAR_VR;
        //--------------------------------------------------------------------//
        U01Trans_STAR_MC.FreeAlignments();
        U12Trans_MR_STAR.FreeAlignments();

        SlidePartitionUpDiagonal
        ( UTL, /**/ UTR,  U00, /**/ U01, U02,
         /*************/ /******************/
               /**/       U10, /**/ U11, U12,
          UBL, /**/ UBR,  U20, /**/ U21, U22 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace triangular_inverse
} // namespace elem

#endif // ifndef LAPACK_TRIANGULARINVERSE_UVAR3_HPP
