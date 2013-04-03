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
#include "elemental/lapack-like/Norm/Frobenius.hpp"
#include "elemental/matrices/Zeros.hpp"

namespace elem {
namespace svd {

#ifdef HAVE_PMRRR

template<typename F>
inline void
ThresholdedTall
( DistMatrix<F>& A,
  DistMatrix<typename Base<F>::type,VR,STAR>& s,
  DistMatrix<F>& V,
  typename Base<F>::type tol=0 )
{
#ifndef RELEASE
    PushCallStack("svd::ThresholdedTall");
    if( A.Height() < A.Width() )
        throw std::logic_error("A must be at least as tall as it is wide");
    if( tol < 0 )
        throw std::logic_error("negative threshold does not make sense");
#endif
    typedef typename Base<F>::type R;
    const Grid& g = A.Grid();
    const int m = A.Height();
    const int n = A.Width();
    const R frobNorm = FrobeniusNorm( A );
    if( tol == R(0) )
    {
        const R eps = lapack::MachineEpsilon<R>();
        tol = m*frobNorm*eps;
    }

    // C := A^H A
    DistMatrix<F> C( g );
    Zeros( n, n, C );
    Herk( LOWER, ADJOINT, F(1), A, F(0), C );

    // [V,Sigma^2] := eig(C), where each sigma > tol
    HermitianEig( LOWER, C, s, V, tol*tol, frobNorm*frobNorm );
    
    // Sigma := sqrt(Sigma^2)
    {
        const int localHeight = s.LocalHeight();
        for( int iLocal=0; iLocal<localHeight; ++iLocal )
            s.SetLocal( iLocal, 0, Sqrt(s.GetLocal(iLocal,0)) );
    }

    // Y := A V
    const int k = V.Width();
    DistMatrix<F> Y( g );
    Zeros( m, k, Y );
    Gemm( NORMAL, NORMAL, F(1), A, V, F(0), Y );

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
        for( int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            const R sigma = s_MR_STAR.GetLocal( jLocal, 0 );
            for( int iLocal=0; iLocal<localHeight; ++iLocal )
                A.SetLocal( iLocal, jLocal, A.GetLocal(iLocal,jLocal)/sigma );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
ThresholdedWide
( DistMatrix<F>& A,
  DistMatrix<typename Base<F>::type,VR,STAR>& s,
  DistMatrix<F>& V,
  typename Base<F>::type tol=0 )
{
#ifndef RELEASE
    PushCallStack("svd::ThresholdedWide");
    if( A.Width() < A.Height() )
        throw std::logic_error("A must be at least as wide as it is tall");
    if( tol < 0 )
        throw std::logic_error("negative threshold does not make sense");
#endif
    typedef typename Base<F>::type R;
    const Grid& g = A.Grid();
    const int m = A.Height();
    const int n = A.Width();
    const R frobNorm = FrobeniusNorm( A );
    if( tol == R(0) )
    {
        const R eps = lapack::MachineEpsilon<R>();
        tol = n*frobNorm*eps;
    }

    // C := A A^H
    DistMatrix<F> C( g );
    Zeros( m, m, C );
    Herk( LOWER, NORMAL, F(1), A, F(0), C );

    // [U,Sigma^2] := eig(C), where each sigma > tol
    DistMatrix<F> U( g );
    HermitianEig( LOWER, C, s, U, tol*tol, frobNorm*frobNorm );
    
    // Sigma := sqrt(Sigma^2)
    {
        const int localHeight = s.LocalHeight();
        for( int iLocal=0; iLocal<localHeight; ++iLocal )
            s.SetLocal( iLocal, 0, Sqrt(s.GetLocal(iLocal,0)) );
    }

    // (Sigma V) := A^H U
    const int k = U.Width();
    Zeros( n, k, V );
    Gemm( ADJOINT, NORMAL, F(1), A, U, F(0), V );

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
        for( int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            const R sigma = s_MR_STAR.GetLocal( jLocal, 0 );
            for( int iLocal=0; iLocal<localHeight; ++iLocal )
                V.SetLocal( iLocal, jLocal, V.GetLocal(iLocal,jLocal)/sigma );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
Thresholded
( DistMatrix<F>& A,
  DistMatrix<typename Base<F>::type,VR,STAR>& s,
  DistMatrix<F>& V,
  typename Base<F>::type tol=0 )
{
#ifndef RELEASE
    PushCallStack("svd::Thresholded");
#endif
    if( A.Height() >= A.Width() )
        ThresholdedTall( A, s, V, tol );
    else
        ThresholdedWide( A, s, V, tol );
#ifndef RELEASE
    PopCallStack();
#endif
}

#endif // ifdef HAVE_PMRRR

} // namespace svd
} // namespace elem

#endif // ifndef LAPACK_SVD_CHAN_HPP
