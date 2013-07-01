/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_SIGN_HPP
#define LAPACK_SIGN_HPP

#include "elemental/blas-like/level1/Axpy.hpp"
#include "elemental/blas-like/level1/Scale.hpp"
#include "elemental/blas-like/level3/Gemm.hpp"
#include "elemental/lapack-like/HermitianFunction.hpp"
#include "elemental/lapack-like/Inverse.hpp"
#include "elemental/lapack-like/Norm/Frobenius.hpp"
#include "elemental/lapack-like/Norm/One.hpp"
#include "elemental/lapack-like/Determinant.hpp"

// See Chapter 5 of Nicholas J. Higham's "Functions of Matrices: Theory and
// Computation", which is currently available at:
// http://www.siam.org/books/ot104/OT104HighamChapter5.pdf

namespace elem {
namespace sign {

enum Scaling {
    NONE,    
    DETERMINANT,
    FROB_NORM
};

template<typename F>
inline void
NewtonStep( const Matrix<F>& X, Matrix<F>& XNew, Scaling scaling=FROB_NORM )
{
#ifndef RELEASE
    CallStackEntry entry("sign::NewtonStep");
#endif
    typedef BASE(F) R;

    // Calculate mu while forming XNew := inv(X)
    R mu;
    Matrix<int> p;
    XNew = X;
    LU( XNew, p );
    if( scaling == DETERMINANT )
    {
        SafeProduct<F> det = determinant::AfterLUPartialPiv( XNew, p );
        mu = R(1)/Exp(det.kappa);
    }
    inverse::AfterLUPartialPiv( XNew, p );
    if( scaling == FROB_NORM )
        mu = Sqrt( FrobeniusNorm(XNew)/FrobeniusNorm(X) );
    else if( scaling == NONE )
        mu = 1;

    // Overwrite XNew with the new iterate
    const R halfMu = mu/R(2);
    const R halfMuInv = R(1)/(2*mu); 
    Scale( halfMuInv, XNew );
    Axpy( halfMu, X, XNew );
}

template<typename F>
inline void
NewtonStep
( const DistMatrix<F>& X, DistMatrix<F>& XNew, Scaling scaling=FROB_NORM )
{
#ifndef RELEASE
    CallStackEntry entry("sign::NewtonStep");
#endif
    typedef BASE(F) R;

    // Calculate mu while forming B := inv(X)
    R mu;
    DistMatrix<int,VC,STAR> p( X.Grid() );
    XNew = X;
    LU( XNew, p );
    if( scaling == DETERMINANT )
    {
        SafeProduct<F> det = determinant::AfterLUPartialPiv( XNew, p );
        mu = R(1)/Exp(det.kappa);
    }
    inverse::AfterLUPartialPiv( XNew, p );
    if( scaling == FROB_NORM )
        mu = Sqrt( FrobeniusNorm(XNew)/FrobeniusNorm(X) );
    else if( scaling == NONE )
        mu = 1;

    // Overwrite XNew with the new iterate
    const R halfMu = mu/R(2);
    const R halfMuInv = R(1)/(2*mu); 
    Scale( halfMuInv, XNew );
    Axpy( halfMu, X, XNew );
}

template<typename F>
inline void
NewtonSchulzStep( const Matrix<F>& X, Matrix<F>& XTmp, Matrix<F>& XNew )
{
#ifndef RELEASE
    CallStackEntry entry("sign::NewtonSchulzStep");
#endif
    typedef BASE(F) R;
    const int n = X.Height();
 
    // XTmp := 3I - X^2
    Identity( XTmp, n, n );
    Gemm( NORMAL, NORMAL, R(-1), X, X, R(3), XTmp );

    // XNew := 1/2 X XTmp
    Gemm( NORMAL, NORMAL, R(1)/R(2), X, XTmp, XNew );
}

template<typename F>
inline void
NewtonSchulzStep
( const DistMatrix<F>& X, DistMatrix<F>& XTmp, DistMatrix<F>& XNew )
{
#ifndef RELEASE
    CallStackEntry entry("sign::NewtonSchulzStep");
#endif
    typedef BASE(F) R;
    const int n = X.Height();

    // XTmp := 3I - X^2
    Identity( XTmp, n, n );
    Gemm( NORMAL, NORMAL, R(-1), X, X, R(3), XTmp );

    // XNew := 1/2 X XTmp
    Gemm( NORMAL, NORMAL, R(1)/R(2), X, XTmp, XNew );
}

template<typename F>
inline int
Newton
( Matrix<F>& A, Scaling scaling=FROB_NORM, int maxIts=100, BASE(F) tol=0 )
{
#ifndef RELEASE
    CallStackEntry entry("sign::Newton");
#endif
    typedef BASE(F) R;
    Matrix<F> B;
    Matrix<F> *X=&A, *XNew=&B;

    if( tol == R(0) )
        tol = A.Height()*lapack::MachineEpsilon<R>();

    int numIts=0;
    while( numIts < maxIts )
    {
        // Overwrite XNew with the new iterate
        NewtonStep( *X, *XNew, scaling );

        // Use the difference in the iterates to test for convergence
        Axpy( R(-1), *XNew, *X );
        const R oneDiff = OneNorm( *X );
        const R oneNew = OneNorm( *XNew );

        // Ensure that X holds the current iterate and break if possible
        ++numIts;
        std::swap( X, XNew );
        if( oneDiff/oneNew <= tol )
            break;
    }
    if( X != &A )
        A = *XNew;
    return numIts;
}

template<typename F>
inline int
Newton
( DistMatrix<F>& A, Scaling scaling=FROB_NORM, 
  int maxIts=100, BASE(F) tol=0 )
{
#ifndef RELEASE
    CallStackEntry entry("sign::Newton");
#endif
    typedef BASE(F) R;
    DistMatrix<F> B( A.Grid() );
    DistMatrix<F> *X=&A, *XNew=&B;

    if( tol == R(0) )
        tol = A.Height()*lapack::MachineEpsilon<R>();

    int numIts=0;
    while( numIts < maxIts )
    {
        // Overwrite XNew with the new iterate
        NewtonStep( *X, *XNew, scaling );

        // Use the difference in the iterates to test for convergence
        Axpy( R(-1), *XNew, *X );
        const R oneDiff = OneNorm( *X );
        const R oneNew = OneNorm( *XNew );

        // Ensure that X holds the current iterate and break if possible
        ++numIts;
        std::swap( X, XNew );
        if( oneDiff/oneNew < tol )
            break;
    }
    if( X != &A )
        A = *XNew;
    return numIts;
}

// TODO: NewtonSchulzHybrid which estimates when || X^2 - I ||_2 < 1

} // namespace sign

template<typename F>
inline void
Sign( Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("Sign");
#endif
    sign::Newton( A );
}

template<typename F>
inline void
Sign( Matrix<F>& A, Matrix<F>& N )
{
#ifndef RELEASE
    PushCallStack("Sign");
#endif
    Matrix<F> ACopy( A );
    sign::Newton( A );
    Gemm( NORMAL, NORMAL, F(1), A, ACopy, N );
}

template<typename F>
inline void
Sign( DistMatrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("Sign");
#endif
    sign::Newton( A );
}

template<typename F>
inline void
Sign( DistMatrix<F>& A, DistMatrix<F>& N )
{
#ifndef RELEASE
    PushCallStack("Sign");
#endif
    DistMatrix<F> ACopy( A );
    sign::Newton( A );
    Gemm( NORMAL, NORMAL, F(1), A, ACopy, N );
}

// The Hermitian sign decomposition is equivalent to the Hermitian polar
// decomposition... A = (U sgn(Lambda) U') (U sgn(Lambda)Lambda U')
//                    = (U sgn(Lambda) U') (U |Lambda| U')

// Even though sgn(lambda) isn't well-defined when lambda=0, we will extend it
// from the right so that the sign decomposition of a singular Hermitian matrix
// is a polar decomposition (which always exists).

template<typename F>
inline void
HermitianSign( UpperOrLower uplo, Matrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("HermitianSign");
#endif
    typedef BASE(F) R;

    // Get the EVD of A
    Matrix<R> w;
    Matrix<F> Z;
    HermitianEig( uplo, A, w, Z );

    const int n = A.Height();
    for( int i=0; i<n; ++i )
    {
        const R omega = w.Get(i,0);
        if( omega >= 0 )
            w.Set(i,0,R(1));
        else
            w.Set(i,0,R(-1));
    }

    // Reform the Hermitian matrix with the modified eigenvalues
    HermitianFromEVD( uplo, A, w, Z );
}

template<typename F>
inline void
HermitianSign( UpperOrLower uplo, Matrix<F>& A, Matrix<F>& N )
{
#ifndef RELEASE
    CallStackEntry entry("HermitianSign");
#endif
    typedef BASE(F) R;

    // Get the EVD of A
    Matrix<R> w;
    Matrix<F> Z;
    HermitianEig( uplo, A, w, Z );

    const int n = A.Height();
    Matrix<R> wSgn( n, 1 ), wAbs( n, 1 );
    for( int i=0; i<n; ++i )
    {
        const R omega = w.Get(i,0);
        if( omega >= 0 )
        {
            wSgn.Set(i,0,R(1));
            wAbs.Set(i,0,omega);
        }
        else
        {
            wSgn.Set(i,0,R(-1));
            wAbs.Set(i,0,-omega);
        }
    }

    // Form the Hermitian matrices with modified eigenvalues
    HermitianFromEVD( uplo, A, wSgn, Z );
    HermitianFromEVD( uplo, N, wAbs, Z );
}

template<typename F>
inline void
HermitianSign( UpperOrLower uplo, DistMatrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("HermitianSign");
#endif
    EnsurePMRRR();
    typedef BASE(F) R;

    // Get the EVD of A
    const Grid& g = A.Grid();
    DistMatrix<R,VR,STAR> w(g);
    DistMatrix<F> Z(g);
    HermitianEig( uplo, A, w, Z );

    const int numLocalEigs = w.LocalHeight();
    for( int iLoc=0; iLoc<numLocalEigs; ++iLoc )
    {
        const R omega = w.GetLocal(iLoc,0);
        if( omega >= 0 )
            w.SetLocal(iLoc,0,R(1));
        else
            w.SetLocal(iLoc,0,R(-1));
    }

    // Reform the Hermitian matrix with the modified eigenvalues
    HermitianFromEVD( uplo, A, w, Z );
}

template<typename F>
inline void
HermitianSign( UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<F>& N )
{
#ifndef RELEASE
    CallStackEntry entry("HermitianSign");
#endif
    EnsurePMRRR();
    typedef BASE(F) R;

    // Get the EVD of A
    const Grid& g = A.Grid();
    DistMatrix<R,VR,STAR> w(g);
    DistMatrix<F> Z(g);
    HermitianEig( uplo, A, w, Z );

    const int n = A.Height();
    const int numLocalEigs = w.LocalHeight();
    DistMatrix<R,VR,STAR> wSgn(g), wAbs(g);
    wSgn.AlignWith( w );
    wAbs.AlignWith( w );
    wSgn.ResizeTo( n, 1 );
    wAbs.ResizeTo( n, 1 );
    for( int iLoc=0; iLoc<numLocalEigs; ++iLoc )
    {
        const R omega = w.GetLocal(iLoc,0);
        if( omega >= 0 )
        {
            wSgn.SetLocal(iLoc,0,R(1));
            wAbs.SetLocal(iLoc,0,omega);
        }
        else
        {
            wSgn.SetLocal(iLoc,0,R(-1));
            wAbs.SetLocal(iLoc,0,-omega);
        }
    }

    // Form the Hermitian matrix with the modified eigenvalues
    HermitianFromEVD( uplo, A, wSgn, Z );
    HermitianFromEVD( uplo, N, wAbs, Z );
}

} // namespace elem

#endif // ifndef LAPACK_SIGN_HPP
