/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_HERMITIANFUNCTION_HPP
#define LAPACK_HERMITIANFUNCTION_HPP

#include "elemental/blas-like/level1/DiagonalScale.hpp"
#include "elemental/blas-like/level1/MakeTrapezoidal.hpp"
#include "elemental/blas-like/level3/Gemm.hpp"
#include "elemental/lapack-like/HermitianEig.hpp"
#include "elemental/matrices/Zeros.hpp"

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
    CallStackEntry entry("HermitianFromEVD");
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
inline void
HermitianFromEVD
( UpperOrLower uplo,
        DistMatrix<F>& A,
  const DistMatrix<BASE(F),VR,STAR>& w,
  const DistMatrix<F>& Z )
{
#ifndef RELEASE
    CallStackEntry entry("HermitianFromEVD");
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

// A :=  Z Omega Z^T, where Omega is complex-valued and diagonal

template<typename R>
inline void
NormalFromEVD
(       Matrix<Complex<R> >& A,
  const Matrix<Complex<R> >& w,
  const Matrix<Complex<R> >& Z )
{
#ifndef RELEASE
    CallStackEntry entry("NormalFromEVD");
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
inline void
NormalFromEVD
(       DistMatrix<Complex<R> >& A,
  const DistMatrix<Complex<R>,VR,STAR>& w,
  const DistMatrix<Complex<R> >& Z )
{
#ifndef RELEASE
    CallStackEntry entry("NormalFromEVD");
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

//
// Modify the eigenvalues of A with the real-valued function f, which will 
// therefore result in a Hermitian matrix, which we store in-place.
//

template<typename F,class RealFunctor>
inline void
RealHermitianFunction
( UpperOrLower uplo, Matrix<F>& A, const RealFunctor& f )
{
#ifndef RELEASE
    CallStackEntry entry("RealHermitianFunction");
#endif
    if( A.Height() != A.Width() )
        throw std::logic_error("Hermitian matrices must be square");
    typedef BASE(F) R;

    // Get the EVD of A
    Matrix<R> w;
    Matrix<F> Z;
    HermitianEig( uplo, A, w, Z );

    // Replace w with f(w)
    const int n = w.Height();
    for( int i=0; i<n; ++i )
    {
        const R omega = w.Get(i,0);
        w.Set(i,0,f(omega));
    }

    // A := Z f(Omega) Z^H
    HermitianFromEVD( uplo, A, w, Z );
}

template<typename F,class RealFunctor>
inline void
RealHermitianFunction
( UpperOrLower uplo, DistMatrix<F>& A, const RealFunctor& f )
{
#ifndef RELEASE
    CallStackEntry entry("RealHermitianFunction");
#endif
    EnsurePMRRR();
    if( A.Height() != A.Width() )
        throw std::logic_error("Hermitian matrices must be square");
    typedef BASE(F) R;

    // Get the EVD of A
    const Grid& g = A.Grid();
    DistMatrix<R,VR,STAR> w(g);
    DistMatrix<F> Z(g);
    HermitianEig( uplo, A, w, Z );

    // Replace w with f(w)
    const int numLocalEigs = w.LocalHeight();
    for( int iLoc=0; iLoc<numLocalEigs; ++iLoc )
    {
        const R omega = w.GetLocal(iLoc,0);
        w.SetLocal(iLoc,0,f(omega));
    }

    // A := Z f(Omega) Z^H
    HermitianFromEVD( uplo, A, w, Z ); 
}

//
// Modify the eigenvalues of A with the complex-valued function f, which will
// therefore result in a normal (in general, non-Hermitian) matrix, which we 
// store in-place. At some point a version will be written which takes a real
// symmetric matrix as input and produces a complex normal matrix.
//

template<typename R,class ComplexFunctor>
inline void
ComplexHermitianFunction
( UpperOrLower uplo, Matrix<Complex<R> >& A, const ComplexFunctor& f )
{
#ifndef RELEASE
    CallStackEntry entry("ComplexHermitianFunction");
#endif
    if( A.Height() != A.Width() )
        throw std::logic_error("Hermitian matrices must be square");
    typedef Complex<R> C;

    // Get the EVD of A
    Matrix<R> w;
    Matrix<C> Z;
    HermitianEig( uplo, A, w, Z );

    // Form f(w)
    const int n = w.Height();
    Matrix<C> fw( n, 1 );
    for( int i=0; i<n; ++i )
    {
        const R omega = w.Get(i,0);
        fw.Set(i,0,f(omega));
    }

    // A := Z f(Omega) Z^H
    NormalFromEVD( A, fw, Z );
}

template<typename R,class ComplexFunctor>
inline void
ComplexHermitianFunction
( UpperOrLower uplo, DistMatrix<Complex<R> >& A, const ComplexFunctor& f )
{
#ifndef RELEASE
    CallStackEntry entry("ComplexHermitianFunction");
#endif
    EnsurePMRRR();
    if( A.Height() != A.Width() )
        throw std::logic_error("Hermitian matrices must be square");
    typedef Complex<R> C;

    // Get the EVD of A
    const Grid& g = A.Grid();
    DistMatrix<R,VR,STAR> w(g);
    DistMatrix<C> Z(g);
    HermitianEig( uplo, A, w, Z );

    // Form f(w)
    DistMatrix<C,VR,STAR> fw(g);
    fw.AlignWith( w );
    fw.ResizeTo( w.Height(), 1 );
    const int numLocalEigs = w.LocalHeight();
    for( int iLoc=0; iLoc<numLocalEigs; ++iLoc )
    {
        const R omega = w.GetLocal(iLoc,0);
        fw.SetLocal(iLoc,0,f(omega));
    }

    // A := Z f(Omega) Z^H
    NormalFromEVD( A, fw, Z );
}

} // namespace elem

#endif // ifndef LAPACK_HERMITIANFUNCTION_HPP
