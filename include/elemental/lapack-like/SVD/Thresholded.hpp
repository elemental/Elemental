/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_SVD_THRESHOLDED_HPP
#define LAPACK_SVD_THRESHOLDED_HPP

#include "elemental/blas-like/level3/Gemm.hpp"
#include "elemental/blas-like/level3/Herk.hpp"
#include "elemental/lapack-like/HermitianEig.hpp"
#include "elemental/lapack-like/Norm/Frobenius.hpp"

namespace elem {
namespace svd {

template<typename F>
inline void
ThresholdedTall
( Matrix<F>& A, Matrix<BASE(F)>& s, Matrix<F>& V, BASE(F) tol=0 )
{
#ifndef RELEASE
    CallStackEntry entry("svd::ThresholdedTall");
    if( A.Height() < A.Width() )
        throw std::logic_error("A must be at least as tall as it is wide");
    if( tol < 0 )
        throw std::logic_error("negative threshold does not make sense");
#endif
    typedef BASE(F) R;
    const int m = A.Height();
    const int n = A.Width();
    const R frobNorm = FrobeniusNorm( A );
    if( tol == R(0) )
    {
        const R eps = lapack::MachineEpsilon<R>();
        tol = m*frobNorm*eps;
    }

    // C := A^H A
    Matrix<F> C;
    Herk( LOWER, ADJOINT, F(1), A, C );

    // [V,Sigma^2] := eig(C), where each sigma > tol
    HermitianEig( LOWER, C, s, V, tol*tol, frobNorm*frobNorm );
    
    // Sigma := sqrt(Sigma^2)
    const int k = s.Height();
    for( int i=0; i<k; ++i )
        s.Set( i, 0, Sqrt(s.Get(i,0)) );

    // Y := A V
    Matrix<F> Y;
    Gemm( NORMAL, NORMAL, F(1), A, V, Y );

    // Set each column of A to be the corresponding normalized column of Y
    // NOTE: A (potentially better) alternative would be to compute the norm of
    //       each column of A and normalize via it, as it might vary slightly
    //       from the corresponding computed singular value.
    A = Y;
    for( int j=0; j<n; ++j )
    {
        const R sigma = s.Get( j, 0 );
        for( int i=0; i<m; ++i )
            A.Set( i, j, A.Get(i,j)/sigma );
    }
}

#ifdef HAVE_PMRRR
template<typename F>
inline void
ThresholdedTall
( DistMatrix<F>& A, DistMatrix<BASE(F),VR,STAR>& s, DistMatrix<F>& V,
  BASE(F) tol=0 )
{
#ifndef RELEASE
    CallStackEntry entry("svd::ThresholdedTall");
    if( A.Height() < A.Width() )
        throw std::logic_error("A must be at least as tall as it is wide");
    if( tol < 0 )
        throw std::logic_error("negative threshold does not make sense");
#endif
    typedef BASE(F) R;
    const Grid& g = A.Grid();
    const int m = A.Height();
    const R frobNorm = FrobeniusNorm( A );
    if( tol == R(0) )
    {
        const R eps = lapack::MachineEpsilon<R>();
        tol = m*frobNorm*eps;
    }

    // C := A^H A
    DistMatrix<F> C( g );
    Herk( LOWER, ADJOINT, F(1), A, C );

    // [V,Sigma^2] := eig(C), where each sigma > tol
    HermitianEig( LOWER, C, s, V, tol*tol, frobNorm*frobNorm );
    
    // Sigma := sqrt(Sigma^2)
    {
        const int localHeight = s.LocalHeight();
        for( int iLoc=0; iLoc<localHeight; ++iLoc )
            s.SetLocal( iLoc, 0, Sqrt(s.GetLocal(iLoc,0)) );
    }

    // Y := A V
    DistMatrix<F> Y( g );
    Gemm( NORMAL, NORMAL, F(1), A, V, Y );

    // Set each column of A to be the corresponding normalized column of Y
    // NOTE: A (potentially better) alternative would be to compute the norm of
    //       each column of A and normalize via it, as it might vary slightly
    //       from the corresponding computed singular value.
    A = Y;
    {
        DistMatrix<R,MR,STAR> s_MR_STAR( g );
        s_MR_STAR.AlignWith( A.DistData() );
        s_MR_STAR = s;
        const int localWidth = A.LocalWidth();
        const int localHeight = A.LocalHeight();
        for( int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const R sigma = s_MR_STAR.GetLocal( jLoc, 0 );
            for( int iLoc=0; iLoc<localHeight; ++iLoc )
                A.SetLocal( iLoc, jLoc, A.GetLocal(iLoc,jLoc)/sigma );
        }
    }
}
#endif // ifdef HAVE_PMRRR

template<typename F>
inline void
ThresholdedWide
( Matrix<F>& A, Matrix<BASE(F)>& s, Matrix<F>& V, BASE(F) tol=0 )
{
#ifndef RELEASE
    CallStackEntry entry("svd::ThresholdedWide");
    if( A.Width() < A.Height() )
        throw std::logic_error("A must be at least as wide as it is tall");
    if( tol < 0 )
        throw std::logic_error("negative threshold does not make sense");
#endif
    typedef BASE(F) R;
    const int n = A.Width();
    const R frobNorm = FrobeniusNorm( A );
    if( tol == R(0) )
    {
        const R eps = lapack::MachineEpsilon<R>();
        tol = n*frobNorm*eps;
    }

    // C := A A^H
    Matrix<F> C;
    Herk( LOWER, NORMAL, F(1), A, C );

    // [U,Sigma^2] := eig(C), where each sigma > tol
    Matrix<F> U;
    HermitianEig( LOWER, C, s, U, tol*tol, frobNorm*frobNorm );
    
    // Sigma := sqrt(Sigma^2)
    const int k = s.Height();
    for( int i=0; i<k; ++i )
        s.Set( i, 0, Sqrt(s.Get(i,0)) );

    // (Sigma V) := A^H U
    Gemm( ADJOINT, NORMAL, F(1), A, U, V );

    // Divide each column of (Sigma V) by sigma
    // NOTE: A (potentially better) alternative would be to compute the norm of
    //       each column of V and normalize via it, as it might vary slightly
    //       from the corresponding computed singular value.
    for( int j=0; j<k; ++j )
    {
        const R sigma = s.Get( j, 0 );
        for( int i=0; i<n; ++i )
            V.Set( i, j, V.Get(i,j)/sigma );
    }
    A = U;
}

#ifdef HAVE_PMRRR
template<typename F>
inline void
ThresholdedWide
( DistMatrix<F>& A, DistMatrix<BASE(F),VR,STAR>& s, DistMatrix<F>& V,
  BASE(F) tol=0 )
{
#ifndef RELEASE
    CallStackEntry entry("svd::ThresholdedWide");
    if( A.Width() < A.Height() )
        throw std::logic_error("A must be at least as wide as it is tall");
    if( tol < 0 )
        throw std::logic_error("negative threshold does not make sense");
#endif
    typedef BASE(F) R;
    const Grid& g = A.Grid();
    const int n = A.Width();
    const R frobNorm = FrobeniusNorm( A );
    if( tol == R(0) )
    {
        const R eps = lapack::MachineEpsilon<R>();
        tol = n*frobNorm*eps;
    }

    // C := A A^H
    DistMatrix<F> C( g );
    Herk( LOWER, NORMAL, F(1), A, C );

    // [U,Sigma^2] := eig(C), where each sigma > tol
    DistMatrix<F> U( g );
    HermitianEig( LOWER, C, s, U, tol*tol, frobNorm*frobNorm );
    
    // Sigma := sqrt(Sigma^2)
    {
        const int localHeight = s.LocalHeight();
        for( int iLoc=0; iLoc<localHeight; ++iLoc )
            s.SetLocal( iLoc, 0, Sqrt(s.GetLocal(iLoc,0)) );
    }

    // (Sigma V) := A^H U
    Gemm( ADJOINT, NORMAL, F(1), A, U, V );

    // Divide each column of (Sigma V) by sigma
    // NOTE: A (potentially better) alternative would be to compute the norm of
    //       each column of V and normalize via it, as it might vary slightly
    //       from the corresponding computed singular value.
    {
        DistMatrix<R,MR,STAR> s_MR_STAR( g );
        s_MR_STAR.AlignWith( V.DistData() );
        s_MR_STAR = s;
        const int localWidth = V.LocalWidth();
        const int localHeight = V.LocalHeight();
        for( int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const R sigma = s_MR_STAR.GetLocal( jLoc, 0 );
            for( int iLoc=0; iLoc<localHeight; ++iLoc )
                V.SetLocal( iLoc, jLoc, V.GetLocal(iLoc,jLoc)/sigma );
        }
    }
    A = U;
}
#endif // ifdef HAVE_PMRRR

template<typename F>
inline void
Thresholded
( Matrix<F>& A, Matrix<BASE(F)>& s, Matrix<F>& V, BASE(F) tol=0 )
{
#ifndef RELEASE
    CallStackEntry entry("svd::Thresholded");
#endif
    if( A.Height() >= A.Width() )
        ThresholdedTall( A, s, V, tol );
    else
        ThresholdedWide( A, s, V, tol );
}

#ifdef HAVE_PMRRR
template<typename F>
inline void
Thresholded
( DistMatrix<F>& A, DistMatrix<BASE(F),VR,STAR>& s, DistMatrix<F>& V,
  BASE(F) tol=0 )
{
#ifndef RELEASE
    CallStackEntry entry("svd::Thresholded");
#endif
    if( A.Height() >= A.Width() )
        ThresholdedTall( A, s, V, tol );
    else
        ThresholdedWide( A, s, V, tol );
}
#endif // ifdef HAVE_PMRRR

} // namespace svd
} // namespace elem

#endif // ifndef LAPACK_SVD_THRESHOLDED_HPP
