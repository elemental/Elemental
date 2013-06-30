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

} // namespace elem

#endif // ifndef LAPACK_SIGN_HPP
