/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_PSEUDOINVERSE_HPP
#define ELEM_LAPACK_PSEUDOINVERSE_HPP

#include "elemental/blas-like/level1/DiagonalScale.hpp"
#include "elemental/blas-like/level3/Gemm.hpp"
#include "elemental/lapack-like/Norm/Max.hpp"
#include "elemental/lapack-like/SVD.hpp"
#include "elemental/matrices/HermitianFromEVD.hpp"

namespace elem {

//
// Replace A with its pseudoinverse
//

template<typename F>
inline void
Pseudoinverse( Matrix<F>& A, BASE(F) tolerance=0 )
{
#ifndef RELEASE
    CallStackEntry entry("Pseudoinverse");
#endif
    typedef BASE(F) R;

    const int m = A.Height();
    const int n = A.Width();
    const int k = std::max(m,n);

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
    const int numVals = s.Height();
    for( int i=0; i<numVals; ++i )
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
HermitianPseudoinverse( UpperOrLower uplo, Matrix<F>& A, BASE(F) tolerance=0 )
{
#ifndef RELEASE
    CallStackEntry entry("HermitianPseudoinverse");
#endif
    typedef BASE(F) R;
    const int n = A.Height();

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
    for( int i=0; i<n; ++i )
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
Pseudoinverse( DistMatrix<F>& A, BASE(F) tolerance=0 )
{
#ifndef RELEASE
    CallStackEntry entry("Pseudoinverse");
#endif
    typedef BASE(F) R;

    const Grid& g = A.Grid();
    const int m = A.Height();
    const int n = A.Width();
    const int k = std::max(m,n);

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
    const int numLocalVals = s.LocalHeight();
    for( int iLoc=0; iLoc<numLocalVals; ++iLoc )
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
( UpperOrLower uplo, DistMatrix<F>& A, BASE(F) tolerance=0 )
{
#ifndef RELEASE
    CallStackEntry entry("HermitianPseudoinverse");
#endif
    EnsurePMRRR();
    typedef BASE(F) R;
    const int n = A.Height();

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
    const int numLocalEigs = w.LocalHeight();
    for( int iLoc=0; iLoc<numLocalEigs; ++iLoc )
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

#endif // ifndef ELEM_LAPACK_PSEUDOINVERSE_HPP
