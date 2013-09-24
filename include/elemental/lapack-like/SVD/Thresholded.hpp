/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_SVD_THRESHOLDED_HPP
#define ELEM_LAPACK_SVD_THRESHOLDED_HPP

#include "elemental/blas-like/level3/Gemm.hpp"
#include "elemental/blas-like/level3/Herk.hpp"
#include "elemental/lapack-like/HermitianEig.hpp"
#include "elemental/lapack-like/Norm/Frobenius.hpp"

namespace elem {
namespace svd {

template<typename F>
inline void
ThresholdedTall( Matrix<F>& A, Matrix<BASE(F)>& s, Matrix<F>& V, BASE(F) tol=0 )
{
#ifndef RELEASE
    CallStackEntry entry("svd::ThresholdedTall");
    if( A.Height() < A.Width() )
        LogicError("A must be at least as tall as it is wide");
    if( tol < 0 )
        LogicError("negative threshold does not make sense");
#endif
    typedef BASE(F) Real;
    const Int m = A.Height();
    const Int n = A.Width();
    const Real frobNorm = FrobeniusNorm( A );
    if( tol == Real(0) )
    {
        const Real eps = lapack::MachineEpsilon<Real>();
        tol = m*frobNorm*eps;
    }

    // C := A^H A
    Matrix<F> C;
    Herk( LOWER, ADJOINT, F(1), A, C );

    // [V,Sigma^2] := eig(C), where each sigma > tol
    HermitianEig( LOWER, C, s, V, tol*tol, frobNorm*frobNorm );
    
    // Sigma := sqrt(Sigma^2)
    const Int k = s.Height();
    for( Int i=0; i<k; ++i )
        s.Set( i, 0, Sqrt(s.Get(i,0)) );

    // Y := A V
    Matrix<F> Y;
    Gemm( NORMAL, NORMAL, F(1), A, V, Y );

    // Set each column of A to be the corresponding normalized column of Y
    // NOTE: A (potentially better) alternative would be to compute the norm of
    //       each column of A and normalize via it, as it might vary slightly
    //       from the corresponding computed singular value.
    A = Y;
    for( Int j=0; j<n; ++j )
    {
        const Real sigma = s.Get( j, 0 );
        for( Int i=0; i<m; ++i )
            A.Set( i, j, A.Get(i,j)/sigma );
    }
}

template<typename F>
inline void
ThresholdedTall
( DistMatrix<F>& A, DistMatrix<BASE(F),VR,STAR>& s, DistMatrix<F>& V,
  BASE(F) tol=0 )
{
#ifndef RELEASE
    CallStackEntry entry("svd::ThresholdedTall");
    if( A.Height() < A.Width() )
        LogicError("A must be at least as tall as it is wide");
    if( tol < 0 )
        LogicError("negative threshold does not make sense");
#endif
    EnsurePMRRR();
    typedef BASE(F) Real;
    const Grid& g = A.Grid();
    const Int m = A.Height();
    const Real frobNorm = FrobeniusNorm( A );
    if( tol == Real(0) )
    {
        const Real eps = lapack::MachineEpsilon<Real>();
        tol = m*frobNorm*eps;
    }

    // C := A^H A
    DistMatrix<F> C(g);
    Herk( LOWER, ADJOINT, F(1), A, C );

    // [V,Sigma^2] := eig(C), where each sigma > tol
    HermitianEig( LOWER, C, s, V, tol*tol, frobNorm*frobNorm );
    
    // Sigma := sqrt(Sigma^2)
    {
        const Int localHeight = s.LocalHeight();
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
            s.SetLocal( iLoc, 0, Sqrt(s.GetLocal(iLoc,0)) );
    }

    // Y := A V
    DistMatrix<F> Y(g);
    Gemm( NORMAL, NORMAL, F(1), A, V, Y );

    // Set each column of A to be the corresponding normalized column of Y
    // NOTE: A (potentially better) alternative would be to compute the norm of
    //       each column of A and normalize via it, as it might vary slightly
    //       from the corresponding computed singular value.
    A = Y;
    {
        DistMatrix<Real,MR,STAR> s_MR_STAR(g);
        s_MR_STAR.AlignWith( A.DistData() );
        s_MR_STAR = s;
        const Int localWidth = A.LocalWidth();
        const Int localHeight = A.LocalHeight();
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const Real sigma = s_MR_STAR.GetLocal( jLoc, 0 );
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                A.SetLocal( iLoc, jLoc, A.GetLocal(iLoc,jLoc)/sigma );
        }
    }
}

template<typename F>
inline void
ThresholdedWide( Matrix<F>& A, Matrix<BASE(F)>& s, Matrix<F>& V, BASE(F) tol=0 )
{
#ifndef RELEASE
    CallStackEntry entry("svd::ThresholdedWide");
    if( A.Width() < A.Height() )
        LogicError("A must be at least as wide as it is tall");
    if( tol < 0 )
        LogicError("negative threshold does not make sense");
#endif
    typedef BASE(F) Real;
    const Int n = A.Width();
    const Real frobNorm = FrobeniusNorm( A );
    if( tol == Real(0) )
    {
        const Real eps = lapack::MachineEpsilon<Real>();
        tol = n*frobNorm*eps;
    }

    // C := A A^H
    Matrix<F> C;
    Herk( LOWER, NORMAL, F(1), A, C );

    // [U,Sigma^2] := eig(C), where each sigma > tol
    Matrix<F> U;
    HermitianEig( LOWER, C, s, U, tol*tol, frobNorm*frobNorm );
    
    // Sigma := sqrt(Sigma^2)
    const Int k = s.Height();
    for( Int i=0; i<k; ++i )
        s.Set( i, 0, Sqrt(s.Get(i,0)) );

    // (Sigma V) := A^H U
    Gemm( ADJOINT, NORMAL, F(1), A, U, V );

    // Divide each column of (Sigma V) by sigma
    // NOTE: A (potentially better) alternative would be to compute the norm of
    //       each column of V and normalize via it, as it might vary slightly
    //       from the corresponding computed singular value.
    for( Int j=0; j<k; ++j )
    {
        const Real sigma = s.Get( j, 0 );
        for( Int i=0; i<n; ++i )
            V.Set( i, j, V.Get(i,j)/sigma );
    }
    A = U;
}

template<typename F>
inline void
ThresholdedWide
( DistMatrix<F>& A, DistMatrix<BASE(F),VR,STAR>& s, DistMatrix<F>& V,
  BASE(F) tol=0 )
{
#ifndef RELEASE
    CallStackEntry entry("svd::ThresholdedWide");
    if( A.Width() < A.Height() )
        LogicError("A must be at least as wide as it is tall");
    if( tol < 0 )
        LogicError("negative threshold does not make sense");
#endif
    EnsurePMRRR();
    typedef BASE(F) Real;
    const Grid& g = A.Grid();
    const Int n = A.Width();
    const Real frobNorm = FrobeniusNorm( A );
    if( tol == Real(0) )
    {
        const Real eps = lapack::MachineEpsilon<Real>();
        tol = n*frobNorm*eps;
    }

    // C := A A^H
    DistMatrix<F> C( g );
    Herk( LOWER, NORMAL, F(1), A, C );

    // [U,Sigma^2] := eig(C), where each sigma > tol
    DistMatrix<F> U(g);
    HermitianEig( LOWER, C, s, U, tol*tol, frobNorm*frobNorm );
    
    // Sigma := sqrt(Sigma^2)
    {
        const Int localHeight = s.LocalHeight();
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
            s.SetLocal( iLoc, 0, Sqrt(s.GetLocal(iLoc,0)) );
    }

    // (Sigma V) := A^H U
    Gemm( ADJOINT, NORMAL, F(1), A, U, V );

    // Divide each column of (Sigma V) by sigma
    // NOTE: A (potentially better) alternative would be to compute the norm of
    //       each column of V and normalize via it, as it might vary slightly
    //       from the corresponding computed singular value.
    {
        DistMatrix<Real,MR,STAR> s_MR_STAR( g );
        s_MR_STAR.AlignWith( V.DistData() );
        s_MR_STAR = s;
        const Int localWidth = V.LocalWidth();
        const Int localHeight = V.LocalHeight();
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const Real sigma = s_MR_STAR.GetLocal( jLoc, 0 );
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                V.SetLocal( iLoc, jLoc, V.GetLocal(iLoc,jLoc)/sigma );
        }
    }
    A = U;
}

template<typename F>
inline void
Thresholded( Matrix<F>& A, Matrix<BASE(F)>& s, Matrix<F>& V, BASE(F) tol=0 )
{
#ifndef RELEASE
    CallStackEntry entry("svd::Thresholded");
#endif
    if( A.Height() >= A.Width() )
        ThresholdedTall( A, s, V, tol );
    else
        ThresholdedWide( A, s, V, tol );
}

template<typename F>
inline void
Thresholded
( DistMatrix<F>& A, DistMatrix<BASE(F),VR,STAR>& s, DistMatrix<F>& V,
  BASE(F) tol=0 )
{
#ifndef RELEASE
    CallStackEntry entry("svd::Thresholded");
#endif
    EnsurePMRRR();
    if( A.Height() >= A.Width() )
        ThresholdedTall( A, s, V, tol );
    else
        ThresholdedWide( A, s, V, tol );
}

} // namespace svd
} // namespace elem

#endif // ifndef ELEM_LAPACK_SVD_THRESHOLDED_HPP
