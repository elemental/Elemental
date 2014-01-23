/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_SIGN_HPP
#define ELEM_LAPACK_SIGN_HPP

#include "elemental/blas-like/level1/Axpy.hpp"
#include "elemental/blas-like/level1/Scale.hpp"
#include "elemental/blas-like/level3/Gemm.hpp"
#include "elemental/lapack-like/Inverse.hpp"
#include "elemental/lapack-like/Norm/Frobenius.hpp"
#include "elemental/lapack-like/Norm/One.hpp"
#include "elemental/lapack-like/Determinant.hpp"
#include "elemental/matrices/HermitianFromEVD.hpp"

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
    DEBUG_ONLY(CallStackEntry cse("sign::NewtonStep"))
    typedef Base<F> Real;

    // Calculate mu while forming XNew := inv(X)
    Real mu;
    Matrix<Int> p;
    XNew = X;
    LU( XNew, p );
    if( scaling == DETERMINANT )
    {
        SafeProduct<F> det = determinant::AfterLUPartialPiv( XNew, p );
        mu = Real(1)/Exp(det.kappa);
    }
    inverse::AfterLUPartialPiv( XNew, p );
    if( scaling == FROB_NORM )
        mu = Sqrt( FrobeniusNorm(XNew)/FrobeniusNorm(X) );
    else if( scaling == NONE )
        mu = 1;
    else
        LogicError("Scaling case not handled");

    // Overwrite XNew with the new iterate
    const Real halfMu = mu/Real(2);
    const Real halfMuInv = Real(1)/(2*mu); 
    Scale( halfMuInv, XNew );
    Axpy( halfMu, X, XNew );
}

template<typename F>
inline void
NewtonStep
( const DistMatrix<F>& X, DistMatrix<F>& XNew, Scaling scaling=FROB_NORM )
{
    DEBUG_ONLY(CallStackEntry cse("sign::NewtonStep"))
    typedef Base<F> Real;

    // Calculate mu while forming B := inv(X)
    Real mu;
    DistMatrix<Int,VC,STAR> p( X.Grid() );
    XNew = X;
    LU( XNew, p );
    if( scaling == DETERMINANT )
    {
        SafeProduct<F> det = determinant::AfterLUPartialPiv( XNew, p );
        mu = Real(1)/Exp(det.kappa);
    }
    inverse::AfterLUPartialPiv( XNew, p );
    if( scaling == FROB_NORM )
        mu = Sqrt( FrobeniusNorm(XNew)/FrobeniusNorm(X) );
    else if( scaling == NONE )
        mu = 1;
    else
        LogicError("Scaling case not handled");

    // Overwrite XNew with the new iterate
    const Real halfMu = mu/Real(2);
    const Real halfMuInv = Real(1)/(2*mu); 
    Scale( halfMuInv, XNew );
    Axpy( halfMu, X, XNew );
}

template<typename F>
inline void
NewtonSchulzStep( const Matrix<F>& X, Matrix<F>& XTmp, Matrix<F>& XNew )
{
    DEBUG_ONLY(CallStackEntry cse("sign::NewtonSchulzStep"))
    typedef Base<F> Real;
    const Int n = X.Height();
 
    // XTmp := 3I - X^2
    Identity( XTmp, n, n );
    Gemm( NORMAL, NORMAL, Real(-1), X, X, Real(3), XTmp );

    // XNew := 1/2 X XTmp
    Gemm( NORMAL, NORMAL, Real(1)/Real(2), X, XTmp, XNew );
}

template<typename F>
inline void
NewtonSchulzStep
( const DistMatrix<F>& X, DistMatrix<F>& XTmp, DistMatrix<F>& XNew )
{
    DEBUG_ONLY(CallStackEntry cse("sign::NewtonSchulzStep"))
    typedef Base<F> Real;
    const Int n = X.Height();

    // XTmp := 3I - X^2
    Identity( XTmp, n, n );
    Gemm( NORMAL, NORMAL, Real(-1), X, X, Real(3), XTmp );

    // XNew := 1/2 X XTmp
    Gemm( NORMAL, NORMAL, Real(1)/Real(2), X, XTmp, XNew );
}

// Please see Chapter 5 of Higham's 
// "Functions of Matrices: Theory and Computation" for motivation behind
// the different choices of p, which are usually in {0,1,2}
template<typename F>
inline Int
Newton
( Matrix<F>& A, Scaling scaling=FROB_NORM, Int maxIts=100, BASE(F) tol=0, 
  BASE(F) p=1, bool progress=false )
{
    DEBUG_ONLY(CallStackEntry cse("sign::Newton"))
    typedef Base<F> Real;
    if( tol == Real(0) )
        tol = A.Height()*lapack::MachineEpsilon<Real>();

    Int numIts=0;
    Matrix<F> B;
    Matrix<F> *X=&A, *XNew=&B;
    while( numIts < maxIts )
    {
        // Overwrite XNew with the new iterate
        NewtonStep( *X, *XNew, scaling );

        // Use the difference in the iterates to test for convergence
        Axpy( Real(-1), *XNew, *X );
        const Real oneDiff = OneNorm( *X );
        const Real oneNew = OneNorm( *XNew );

        // Ensure that X holds the current iterate and break if possible
        ++numIts;
        std::swap( X, XNew );
        if( progress )
            std::cout << "after " << numIts << " Newton iter's: " 
                      << "oneDiff=" << oneDiff << ", oneNew=" << oneNew 
                      << ", oneDiff/oneNew=" << oneDiff/oneNew << ", tol=" 
                      << tol << std::endl;
        if( oneDiff/oneNew <= Pow(oneNew,p)*tol )
            break;
    }
    if( X != &A )
        A = *X;
    return numIts;
}

template<typename F>
inline Int
Newton
( DistMatrix<F>& A, Scaling scaling=FROB_NORM, 
  Int maxIts=100, BASE(F) tol=0, BASE(F) p=1, bool progress=false )
{
    DEBUG_ONLY(CallStackEntry cse("sign::Newton"))
    typedef Base<F> Real;
    if( tol == Real(0) )
        tol = A.Height()*lapack::MachineEpsilon<Real>();

    Int numIts=0;
    DistMatrix<F> B( A.Grid() );
    DistMatrix<F> *X=&A, *XNew=&B;
    while( numIts < maxIts )
    {
        // Overwrite XNew with the new iterate
        NewtonStep( *X, *XNew, scaling );

        // Use the difference in the iterates to test for convergence
        Axpy( Real(-1), *XNew, *X );
        const Real oneDiff = OneNorm( *X );
        const Real oneNew = OneNorm( *XNew );

        // Ensure that X holds the current iterate and break if possible
        ++numIts;
        std::swap( X, XNew );
        if( progress && A.Grid().Rank() == 0 )
            std::cout << "after " << numIts << " Newton iter's: "
                      << "oneDiff=" << oneDiff << ", oneNew=" << oneNew
                      << ", oneDiff/oneNew=" << oneDiff/oneNew << ", tol=" 
                      << tol << std::endl;
        if( oneDiff/oneNew <= Pow(oneNew,p)*tol )
            break;
    }
    if( X != &A )
        A = *X;
    return numIts;
}

// TODO: NewtonSchulzHybrid which estimates when || X^2 - I ||_2 < 1

} // namespace sign

template<typename F>
inline void
Sign( Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("Sign"))
    sign::Newton( A );
}

template<typename F>
inline void
Sign( Matrix<F>& A, Matrix<F>& N )
{
    DEBUG_ONLY(CallStackEntry cse("Sign"))
    Matrix<F> ACopy( A );
    sign::Newton( A );
    Gemm( NORMAL, NORMAL, F(1), A, ACopy, N );
}

template<typename F>
inline void
Sign( DistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("Sign"))
    sign::Newton( A );
}

template<typename F>
inline void
Sign( DistMatrix<F>& A, DistMatrix<F>& N )
{
    DEBUG_ONLY(CallStackEntry cse("Sign"))
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
    DEBUG_ONLY(CallStackEntry cse("HermitianSign"))
    typedef Base<F> Real;

    // Get the EVD of A
    Matrix<Real> w;
    Matrix<F> Z;
    HermitianEig( uplo, A, w, Z );

    const Int n = A.Height();
    for( Int i=0; i<n; ++i )
    {
        const Real omega = w.Get(i,0);
        if( omega >= 0 )
            w.Set(i,0,Real(1));
        else
            w.Set(i,0,Real(-1));
    }

    // Reform the Hermitian matrix with the modified eigenvalues
    HermitianFromEVD( uplo, A, w, Z );
}

template<typename F>
inline void
HermitianSign( UpperOrLower uplo, Matrix<F>& A, Matrix<F>& N )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianSign"))
    typedef Base<F> Real;

    // Get the EVD of A
    Matrix<Real> w;
    Matrix<F> Z;
    HermitianEig( uplo, A, w, Z );

    const Int n = A.Height();
    Matrix<Real> wSgn( n, 1 ), wAbs( n, 1 );
    for( Int i=0; i<n; ++i )
    {
        const Real omega = w.Get(i,0);
        if( omega >= 0 )
        {
            wSgn.Set(i,0,Real(1));
            wAbs.Set(i,0,omega);
        }
        else
        {
            wSgn.Set(i,0,Real(-1));
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
    DEBUG_ONLY(CallStackEntry cse("HermitianSign"))
    EnsurePMRRR();
    typedef Base<F> Real;

    // Get the EVD of A
    const Grid& g = A.Grid();
    DistMatrix<Real,VR,STAR> w(g);
    DistMatrix<F> Z(g);
    HermitianEig( uplo, A, w, Z );

    const Int numLocalEigs = w.LocalHeight();
    for( Int iLoc=0; iLoc<numLocalEigs; ++iLoc )
    {
        const Real omega = w.GetLocal(iLoc,0);
        if( omega >= 0 )
            w.SetLocal(iLoc,0,Real(1));
        else
            w.SetLocal(iLoc,0,Real(-1));
    }

    // Reform the Hermitian matrix with the modified eigenvalues
    HermitianFromEVD( uplo, A, w, Z );
}

template<typename F>
inline void
HermitianSign( UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<F>& N )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianSign"))
    EnsurePMRRR();
    typedef Base<F> Real;

    // Get the EVD of A
    const Grid& g = A.Grid();
    DistMatrix<Real,VR,STAR> w(g);
    DistMatrix<F> Z(g);
    HermitianEig( uplo, A, w, Z );

    const Int n = A.Height();
    const Int numLocalEigs = w.LocalHeight();
    DistMatrix<Real,VR,STAR> wSgn(g), wAbs(g);
    wSgn.AlignWith( w );
    wAbs.AlignWith( w );
    wSgn.Resize( n, 1 );
    wAbs.Resize( n, 1 );
    for( Int iLoc=0; iLoc<numLocalEigs; ++iLoc )
    {
        const Real omega = w.GetLocal(iLoc,0);
        if( omega >= 0 )
        {
            wSgn.SetLocal(iLoc,0,Real(1));
            wAbs.SetLocal(iLoc,0,omega);
        }
        else
        {
            wSgn.SetLocal(iLoc,0,Real(-1));
            wAbs.SetLocal(iLoc,0,-omega);
        }
    }

    // Form the Hermitian matrix with the modified eigenvalues
    HermitianFromEVD( uplo, A, wSgn, Z );
    HermitianFromEVD( uplo, N, wAbs, Z );
}

} // namespace elem

#endif // ifndef ELEM_LAPACK_SIGN_HPP
