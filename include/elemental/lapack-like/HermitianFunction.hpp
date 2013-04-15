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

namespace elem {

namespace hermitian_function {

// A :=  Z Omega Z^T, where Omega is diagonal and real-valued

template<typename F>
inline void
ReformHermitianMatrix
( UpperOrLower uplo,
        Matrix<F>& A,
  const Matrix<typename Base<F>::type>& w,
  const Matrix<F>& Z )
{
#ifndef RELEASE
    PushCallStack("hermitian_function::ReformHermitianMatrix");
#endif
    typedef typename Base<F>::type R;

    Matrix<F> ZL, ZR,
              Z0, Z1, Z2;
    Matrix<R> wT,  w0,
              wB,  w1,
                   w2;

    Matrix<F> Z1Copy, Y1;

    if( uplo == LOWER )
        MakeTrapezoidal( LEFT, UPPER, 1, A );
    else
        MakeTrapezoidal( LEFT, LOWER, -1, A );
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
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
ReformHermitianMatrix
( UpperOrLower uplo,
        DistMatrix<F>& A,
  const DistMatrix<typename Base<F>::type,VR,STAR>& w,
  const DistMatrix<F>& Z )
{
#ifndef RELEASE
    PushCallStack("hermitian_function::ReformHermitianMatrix");
#endif
    const Grid& g = A.Grid();
    typedef typename Base<F>::type R;

    DistMatrix<F> ZL(g), ZR(g),
                  Z0(g), Z1(g), Z2(g);
    DistMatrix<R,VR,STAR> wT(g),  w0(g),
                          wB(g),  w1(g),
                                  w2(g);

    DistMatrix<F,MC,  STAR> Z1_MC_STAR(g);
    DistMatrix<F,VR,  STAR> Z1_VR_STAR(g);
    DistMatrix<F,STAR,MR  > Z1Adj_STAR_MR(g);
    DistMatrix<R,STAR,STAR> w1_STAR_STAR(g);

    if( uplo == LOWER )
        MakeTrapezoidal( LEFT, UPPER, 1, A );
    else
        MakeTrapezoidal( LEFT, LOWER, -1, A );
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
        Z1Adj_STAR_MR.FreeAlignments();
        Z1_VR_STAR.FreeAlignments();
        Z1_MC_STAR.FreeAlignments();

        SlideLockedPartitionDown
        ( wT,  w0,
               w1,
         /**/ /**/
          wB,  w2 );
        SlideLockedPartitionRight
        ( ZL,     /**/ ZR,
          Z0, Z1, /**/ Z2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

// A :=  Z Omega Z^T, where Omega is complex-valued and diagonal

template<typename R>
inline void
ReformNormalMatrix
(       Matrix<Complex<R> >& A,
  const Matrix<Complex<R> >& w,
  const Matrix<Complex<R> >& Z )
{
#ifndef RELEASE
    PushCallStack("hermitian_function::ReformNormalMatrix");
#endif
    typedef Complex<R> C;

    Matrix<C> ZL, ZR,
              Z0, Z1, Z2;
    Matrix<C> wT,  w0,
              wB,  w1,
                   w2;

    Matrix<C> Y1, Z1Copy;

    Zero( A );
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
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline void
ReformNormalMatrix
(       DistMatrix<Complex<R> >& A,
  const DistMatrix<Complex<R>,VR,STAR>& w,
  const DistMatrix<Complex<R> >& Z )
{
#ifndef RELEASE
    PushCallStack("hermitian_function::ReformNormalMatrix");
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

    Zero( A );
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
        Z1Adj_STAR_MR.FreeAlignments();
        Z1_VR_STAR.FreeAlignments();
        Z1_MC_STAR.FreeAlignments();

        SlideLockedPartitionDown
        ( wT,  w0,
               w1,
         /**/ /**/
          wB,  w2 );
        SlideLockedPartitionRight
        ( ZL,     /**/ ZR,
          Z0, Z1, /**/ Z2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace hermitian_eig

#ifdef HAVE_PMRRR

//
// Modify the eigenvalues of A with the real-valued function f, which will 
// therefore result in a Hermitian matrix, which we store in-place.
//

template<typename F,class RealFunctor>
inline void
RealHermitianFunction
( UpperOrLower uplo, DistMatrix<F>& A, const RealFunctor& f )
{
#ifndef RELEASE
    PushCallStack("RealHermitianFunction");
#endif
    if( A.Height() != A.Width() )
        throw std::logic_error("Hermitian matrices must be square");
    typedef typename Base<F>::type R;

    // Get the EVD of A
    const Grid& g = A.Grid();
    DistMatrix<R,VR,STAR> w(g);
    DistMatrix<F> Z(g);
    HermitianEig( uplo, A, w, Z );

    // Replace w with f(w)
    const int numLocalEigs = w.LocalHeight();
    for( int iLocal=0; iLocal<numLocalEigs; ++iLocal )
    {
        const R omega = w.GetLocal(iLocal,0);
        w.SetLocal(iLocal,0,f(omega));
    }

    // Form the custom outer product, Z Omega Z^T
    hermitian_function::ReformHermitianMatrix( uplo, A, w, Z );
#ifndef RELEASE
    PopCallStack();
#endif
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
( UpperOrLower uplo, DistMatrix<Complex<R> >& A, const ComplexFunctor& f )
{
#ifndef RELEASE
    PushCallStack("ComplexHermitianFunction");
#endif
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
    fw.AlignWith( w.DistData() );
    fw.ResizeTo( w.Height(), 1 );
    const int numLocalEigs = w.LocalHeight();
    for( int iLocal=0; iLocal<numLocalEigs; ++iLocal )
    {
        const R omega = w.GetLocal(iLocal,0);
        fw.SetLocal(iLocal,0,f(omega));
    }

    // Form the custom outer product, Z f(Omega) Z^H
    hermitian_function::ReformNormalMatrix( A, fw, Z );
#ifndef RELEASE
    PopCallStack();
#endif
}

#endif // ifdef HAVE_PMRRR

} // namespace elem

#endif // ifndef LAPACK_HERMITIANFUNCTION_HPP
