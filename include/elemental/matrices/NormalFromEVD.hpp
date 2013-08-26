/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_NORMALFROMEVD_HPP
#define LAPACK_NORMALFROMEVD_HPP

#include "elemental/blas-like/level1/DiagonalScale.hpp"
#include "elemental/blas-like/level3/Gemm.hpp"
#include "elemental/matrices/Zeros.hpp"

namespace elem {

// A :=  Z Omega Z^T, where Omega is complex-valued and diagonal

template<typename R>
inline void
NormalFromEVD
(       Matrix<Complex<R> >& A,
  const Matrix<Complex<R> >& w,
  const Matrix<Complex<R> >& Z )
{
#ifndef RELEASE
    CallStackEntry cse("NormalFromEVD");
#endif
    typedef Complex<R> C;

    Matrix<C> ZL, ZR,
              Z0, Z1, Z2;
    Matrix<C> wT,  w0,
              wB,  w1,
                   w2;

    Matrix<C> Y1, Z1Copy;

    Zeros( A, Z.Height(), Z.Height() );
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
        DiagonalScale( RIGHT, ADJOINT, w1, Y1 );
        Gemm( NORMAL, NORMAL, C(1), Z1Copy, Y1, C(1), A );
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

template<typename R>
inline Matrix<Complex<R> >
NormalFromEVD
( const Matrix<Complex<R> >& w,
  const Matrix<Complex<R> >& Z )
{
    Matrix<Complex<R> > A;
    NormalFromEVD( A, w, Z );
    return A;
}

template<typename R>
inline void
NormalFromEVD
(       DistMatrix<Complex<R> >& A,
  const DistMatrix<Complex<R>,VR,STAR>& w,
  const DistMatrix<Complex<R> >& Z )
{
#ifndef RELEASE
    CallStackEntry cse("NormalFromEVD");
#endif
    const Grid& g = A.Grid();
    typedef Complex<R> C;

    DistMatrix<C> ZL(g), ZR(g),
                  Z0(g), Z1(g), Z2(g);
    DistMatrix<C,VR,STAR> wT(g),  w0(g),
                          wB(g),  w1(g),
                                  w2(g);

    DistMatrix<C,MC,  STAR> Z1_MC_STAR(g);
    DistMatrix<C,VR,  STAR> Z1_VR_STAR(g);
    DistMatrix<C,STAR,MR  > Z1Adj_STAR_MR(g);
    DistMatrix<C,STAR,STAR> w1_STAR_STAR(g);

    Zeros( A, Z.Height(), Z.Height() );
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

        DiagonalScale( RIGHT, ADJOINT, w1_STAR_STAR, Z1_VR_STAR );

        Z1Adj_STAR_MR.AdjointFrom( Z1_VR_STAR );
        LocalGemm( NORMAL, NORMAL, C(1), Z1_MC_STAR, Z1Adj_STAR_MR, C(1), A );
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

template<typename R>
inline DistMatrix<Complex<R> >
NormalFromEVD
( const DistMatrix<Complex<R>,VR,STAR>& w,
  const DistMatrix<Complex<R> >& Z )
{
    DistMatrix<Complex<R> > A( w.Grid() );
    NormalFromEVD( A, w, Z );
    return A;
}

} // namespace elem

#endif // ifndef LAPACK_NORMALFROMEVD_HPP
