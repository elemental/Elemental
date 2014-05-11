/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_PSEUDOINVERSE_HPP
#define ELEM_PSEUDOINVERSE_HPP

#include ELEM_DIAGONALSCALE_INC
#include ELEM_GEMM_INC
#include ELEM_SVD_INC
#include ELEM_MAXNORM_INC
#include ELEM_HERMITIANFROMEVD_INC

namespace elem {

// Replace A with its pseudoinverse

template<typename F>
inline void
Pseudoinverse( Matrix<F>& A, Base<F> tolerance=0 )
{
    DEBUG_ONLY(CallStackEntry cse("Pseudoinverse"))
    typedef Base<F> R;

    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = Max( m, n );

    // Get the SVD of A
    Matrix<R> s;
    Matrix<F> U, V;
    U = A;
    SVD( U, s, V );

    if( tolerance == R(0) )
    {
        // Set the tolerance equal to k ||A||_2 eps
        const R eps = lapack::MachineEpsilon<R>();
        const R twoNorm = MaxNorm( s );
        tolerance = k*twoNorm*eps;
    }
    // Invert above the tolerance
    const Int numVals = s.Height();
    for( Int i=0; i<numVals; ++i )
    {
        const R sigma = s.Get(i,0);
        if( sigma < tolerance )
            s.Set(i,0,0);
        else
            s.Set(i,0,1/sigma);
    }

    // Scale U with the singular values, U := U Sigma
    DiagonalScale( RIGHT, NORMAL, s, U );

    // Form pinvA = (U Sigma V^H)^H = V (U Sigma)^H
    Gemm( NORMAL, ADJOINT, F(1), V, U, A );
}

template<typename F>
inline void
HermitianPseudoinverse( UpperOrLower uplo, Matrix<F>& A, Base<F> tolerance=0 )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianPseudoinverse"))
    typedef Base<F> R;
    const Int n = A.Height();

    // Get the EVD of A
    Matrix<R> w;
    Matrix<F> Z;
    HermitianEig( uplo, A, w, Z );

    if( tolerance == R(0) )
    {
        // Set the tolerance equal to n ||A||_2 eps
        const R eps = lapack::MachineEpsilon<R>();
        const R twoNorm = MaxNorm( w );
        tolerance = n*twoNorm*eps;
    }
    // Invert above the tolerance
    for( Int i=0; i<n; ++i )
    {
        const R omega = w.Get(i,0);
        if( Abs(omega) < tolerance )
            w.Set(i,0,0);
        else
            w.Set(i,0,1/omega);
    }

    // Form the pseudoinverse
    HermitianFromEVD( uplo, A, w, Z );
}

template<typename F>
inline void
Pseudoinverse( DistMatrix<F>& A, Base<F> tolerance=0 )
{
    DEBUG_ONLY(CallStackEntry cse("Pseudoinverse"))
    typedef Base<F> R;

    const Grid& g = A.Grid();
    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = Max( m, n );

    // Get the SVD of A
    DistMatrix<R,VR,STAR> s(g);
    DistMatrix<F> U(g), V(g);
    U = A;
    SVD( U, s, V );

    if( tolerance == R(0) )
    {
        // Set the tolerance equal to k ||A||_2 eps
        const R eps = lapack::MachineEpsilon<R>();
        const R twoNorm = MaxNorm( s );
        tolerance = k*twoNorm*eps;
    }
    // Invert above the tolerance
    const Int numLocalVals = s.LocalHeight();
    for( Int iLoc=0; iLoc<numLocalVals; ++iLoc )
    {
        const R sigma = s.GetLocal(iLoc,0);
        if( sigma < tolerance )
            s.SetLocal(iLoc,0,0);
        else
            s.SetLocal(iLoc,0,1/sigma);
    }

    // Scale U with the singular values, U := U Sigma
    DiagonalScale( RIGHT, NORMAL, s, U );

    // Form pinvA = (U Sigma V^H)^H = V (U Sigma)^H
    Gemm( NORMAL, ADJOINT, F(1), V, U, A );
}

template<typename F>
inline void
HermitianPseudoinverse
( UpperOrLower uplo, DistMatrix<F>& A, Base<F> tolerance=0 )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianPseudoinverse"))
    typedef Base<F> R;
    const Int n = A.Height();

    // Get the EVD of A
    const Grid& g = A.Grid();
    DistMatrix<R,VR,STAR> w(g);
    DistMatrix<F> Z(g);
    HermitianEig( uplo, A, w, Z );

    if( tolerance == R(0) )
    {
        // Set the tolerance equal to n ||A||_2 eps
        const R eps = lapack::MachineEpsilon<R>();
        const R twoNorm = MaxNorm( w );
        tolerance = n*twoNorm*eps;
    }
    // Invert above the tolerance
    const Int numLocalEigs = w.LocalHeight();
    for( Int iLoc=0; iLoc<numLocalEigs; ++iLoc )
    {
        const R omega = w.GetLocal(iLoc,0);
        if( Abs(omega) < tolerance )
            w.SetLocal(iLoc,0,0);
        else
            w.SetLocal(iLoc,0,1/omega);
    }

    // Form the pseudoinverse
    HermitianFromEVD( uplo, A, w, Z );
}

} // namespace elem

#endif // ifndef ELEM_PSEUDOINVERSE_HPP
