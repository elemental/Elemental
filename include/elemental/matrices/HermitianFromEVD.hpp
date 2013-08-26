/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_HERMITIANFROMEVD_HPP
#define LAPACK_HERMITIANFROMEVD_HPP

#include "elemental/blas-like/level1/DiagonalScale.hpp"
#include "elemental/blas-like/level1/MakeTrapezoidal.hpp"

namespace elem {

// A :=  Z Omega Z^T, where Omega is diagonal and real-valued

template<typename F>
inline void
HermitianFromEVD
( UpperOrLower uplo,
        Matrix<F>& A,
  const Matrix<BASE(F)>& w,
  const Matrix<F>& Z )
{
#ifndef RELEASE
    CallStackEntry cse("HermitianFromEVD");
#endif
    typedef BASE(F) R;

    Matrix<F> ZL, ZR,
              Z0, Z1, Z2;
    Matrix<R> wT,  w0,
              wB,  w1,
                   w2;

    Matrix<F> Z1Copy, Y1;

    A.ResizeTo( Z.Height(), Z.Height() );
    if( uplo == LOWER )
        MakeTrapezoidal( UPPER, A, 1 );
    else
        MakeTrapezoidal( LOWER, A, -1 );
    LockedPartitionRight( Z, ZL, ZR, 0 );
    LockedPartitionDown
    ( w, wT,
         wB, 0 );
    while( ZL.Width() < Z.Width() )
    {
        LockedRepartitionRight
        ( ZL, /**/ ZR,
          Z0, /**/ Z1, Z2 );
        LockedRepartitionDown
        ( wT,  w0,
         /**/ /**/
               w1,
          wB,  w2 );

        //--------------------------------------------------------------------//
        Y1 = Z1Copy = Z1;
        DiagonalScale( RIGHT, NORMAL, w1, Y1 );
        Trrk( uplo, NORMAL, ADJOINT, F(1), Z1Copy, Y1, F(1), A );
        //--------------------------------------------------------------------//

        SlideLockedPartitionDown
        ( wT,  w0,
               w1,
         /**/ /**/
          wB,  w2 );
        SlideLockedPartitionRight
        ( ZL,     /**/ ZR,
          Z0, Z1, /**/ Z2 );
    }
}

template<typename F>
inline Matrix<F>
HermitianFromEVD
( UpperOrLower uplo,
  const Matrix<BASE(F)>& w,
  const Matrix<F>& Z )
{
    Matrix<F> A;
    HermitianFromEVD( uplo, A, w, Z );
    return A;
}

template<typename F>
inline void
HermitianFromEVD
( UpperOrLower uplo,
        DistMatrix<F>& A,
  const DistMatrix<BASE(F),VR,STAR>& w,
  const DistMatrix<F>& Z )
{
#ifndef RELEASE
    CallStackEntry cse("HermitianFromEVD");
#endif
    const Grid& g = A.Grid();
    typedef BASE(F) R;

    DistMatrix<F> ZL(g), ZR(g),
                  Z0(g), Z1(g), Z2(g);
    DistMatrix<R,VR,STAR> wT(g),  w0(g),
                          wB(g),  w1(g),
                                  w2(g);

    DistMatrix<F,MC,  STAR> Z1_MC_STAR(g);
    DistMatrix<F,VR,  STAR> Z1_VR_STAR(g);
    DistMatrix<F,STAR,MR  > Z1Adj_STAR_MR(g);
    DistMatrix<R,STAR,STAR> w1_STAR_STAR(g);

    A.ResizeTo( Z.Height(), Z.Height() );
    if( uplo == LOWER )
        MakeTrapezoidal( UPPER, A, 1 );
    else
        MakeTrapezoidal( LOWER, A, -1 );
    LockedPartitionRight( Z, ZL, ZR, 0 );
    LockedPartitionDown
    ( w, wT,
         wB, 0 );
    while( ZL.Width() < Z.Width() )
    {
        LockedRepartitionRight
        ( ZL, /**/ ZR,
          Z0, /**/ Z1, Z2 );
        LockedRepartitionDown
        ( wT,  w0,
         /**/ /**/
               w1,
          wB,  w2 );

        Z1_MC_STAR.AlignWith( A );
        Z1_VR_STAR.AlignWith( A );
        Z1Adj_STAR_MR.AlignWith( A );
        //--------------------------------------------------------------------//
        Z1_MC_STAR = Z1;
        Z1_VR_STAR = Z1_MC_STAR;
        w1_STAR_STAR = w1;

        DiagonalScale( RIGHT, NORMAL, w1_STAR_STAR, Z1_VR_STAR );

        Z1Adj_STAR_MR.AdjointFrom( Z1_VR_STAR );
        LocalTrrk( uplo, F(1), Z1_MC_STAR, Z1Adj_STAR_MR, F(1), A );
        //--------------------------------------------------------------------//

        SlideLockedPartitionDown
        ( wT,  w0,
               w1,
         /**/ /**/
          wB,  w2 );
        SlideLockedPartitionRight
        ( ZL,     /**/ ZR,
          Z0, Z1, /**/ Z2 );
    }
}

template<typename F>
inline DistMatrix<F>
HermitianFromEVD
( UpperOrLower uplo,
  const DistMatrix<BASE(F),VR,STAR>& w,
  const DistMatrix<F>& Z )
{
    DistMatrix<F> A(w.Grid());
    HermitianFromEVD( uplo, A, w, Z );
    return A;
}

} // namespace elem

#endif // ifndef LAPACK_HERMITIANFROMEVD_HPP
