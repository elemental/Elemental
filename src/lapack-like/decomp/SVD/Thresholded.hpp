/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_SVD_THRESHOLDED_HPP
#define EL_SVD_THRESHOLDED_HPP

// TODO: Use a relative-truncated HermitianEig for relative thresholding

namespace El {
namespace svd {

template<typename F>
inline void TallAbsoluteThresholded
( Matrix<F>& A, Matrix<Base<F>>& s, Matrix<F>& V, Base<F> tol )
{
    DEBUG_ONLY(
        CallStackEntry cse("svd::TallAbsoluteThresholded");
        if( A.Height() < A.Width() )
            LogicError("A must be at least as tall as it is wide");
        if( tol < 0 )
            LogicError("negative threshold does not make sense");
    )
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    const Real frobNorm = FrobeniusNorm( A );
    if( tol == Real(0) )
    {
        const Real eps = lapack::MachineEpsilon<Real>();
        tol = m*frobNorm*eps;
    }
    if( tol >= frobNorm )
    {
        A.Resize( m, 0 );        
        s.Resize( 0, 1 );
        V.Resize( n, 0 );
        return;
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
    for( Int j=0; j<k; ++j )
    {
        const Real sigma = s.Get( j, 0 );
        for( Int i=0; i<m; ++i )
            A.Set( i, j, A.Get(i,j)/sigma );
    }
}

template<typename F>
inline void TallRelativeThresholded
( Matrix<F>& A, Matrix<Base<F>>& s, Matrix<F>& V, Base<F> relTol )
{
    DEBUG_ONLY(
        CallStackEntry cse("svd::TallRelativeThresholded");
        if( A.Height() < A.Width() )
            LogicError("A must be at least as tall as it is wide");
        if( relTol < 0 )
            LogicError("negative threshold does not make sense");
    )
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();

    // C := A^H A
    Matrix<F> C;
    Herk( LOWER, ADJOINT, F(1), A, C );

    // [V,Sigma^2] := eig(C)
    HermitianEig( LOWER, C, s, V, DESCENDING );
    const Real twoNorm = Sqrt(MaxNorm(s));
    
    // Sigma := sqrt(Sigma^2), where all sigmas > relTol*twoNorm
    for( Int i=0; i<n; ++i )
    {
        const Real sigma = Sqrt(s.Get(i,0));
        if( sigma <= relTol*twoNorm )
        {
            s.Resize( i, 1 );
            V.Resize( n, i );
            break;
        }
        else
            s.Set( i, 0, sigma );
    }
    const int k = s.Height();

    // Y := A V
    Matrix<F> Y;
    Gemm( NORMAL, NORMAL, F(1), A, V, Y );

    // Set each column of A to be the corresponding normalized column of Y
    // NOTE: A (potentially better) alternative would be to compute the norm of
    //       each column of A and normalize via it, as it might vary slightly
    //       from the corresponding computed singular value.
    A = Y;
    for( Int j=0; j<k; ++j )
    {
        const Real sigma = s.Get( j, 0 );
        for( Int i=0; i<m; ++i )
            A.Set( i, j, A.Get(i,j)/sigma );
    }
}

template<typename F>
inline void TallThresholded
( Matrix<F>& A, Matrix<Base<F>>& s, Matrix<F>& V, 
  Base<F> tol, bool relative )
{
    DEBUG_ONLY(CallStackEntry cse("svd::TallThresholded"))
    if( relative )
        TallRelativeThresholded( A, s, V, tol );
    else
        TallAbsoluteThresholded( A, s, V, tol );
}

template<typename F>
inline void
TallAbsoluteThresholded
( DistMatrix<F>& A, DistMatrix<Base<F>,VR,STAR>& s, DistMatrix<F>& V,
  Base<F> tol )
{
    DEBUG_ONLY(
        CallStackEntry cse("svd::TallAbsoluteThresholded");
        if( A.Height() < A.Width() )
            LogicError("A must be at least as tall as it is wide");
        if( tol < 0 )
            LogicError("negative threshold does not make sense");
    )
    typedef Base<F> Real;
    const Grid& g = A.Grid();
    const Int m = A.Height();
    const Int n = A.Width();
    const Real frobNorm = FrobeniusNorm( A );
    if( tol == Real(0) )
    {
        const Real eps = lapack::MachineEpsilon<Real>();
        tol = m*frobNorm*eps;
    }
    if( tol >= frobNorm )
    {
        A.Resize( m, 0 );        
        s.Resize( 0, 1 );
        V.Resize( n, 0 );
        return;
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
TallRelativeThresholded
( DistMatrix<F>& A, DistMatrix<Base<F>,VR,STAR>& s, DistMatrix<F>& V,
  Base<F> relTol )
{
    DEBUG_ONLY(
        CallStackEntry cse("svd::TallRelativeThresholded");
        if( A.Height() < A.Width() )
            LogicError("A must be at least as tall as it is wide");
        if( relTol < 0 )
            LogicError("negative threshold does not make sense");
    )
    typedef Base<F> Real;
    const Grid& g = A.Grid();
    const Int n = A.Width();

    // C := A^H A
    DistMatrix<F> C(g);
    Herk( LOWER, ADJOINT, F(1), A, C );

    // [V,Sigma^2] := eig(C)
    HermitianEig( LOWER, C, s, V, DESCENDING );
    const Real twoNorm = Sqrt(MaxNorm(s));

    // Sigma := sqrt(Sigma^2), where all sigmas > relTol*twoNorm
    DistMatrix<Real,STAR,STAR> s_STAR_STAR( s );
    for( Int i=0; i<n; ++i )
    {
        const Real sigma = Sqrt(s_STAR_STAR.GetLocal(i,0));
        if( sigma <= relTol*twoNorm )
        {
            s_STAR_STAR.Resize( i, 1 );
            V.Resize( n, i );
            break;
        }
        else
            s_STAR_STAR.SetLocal( i, 0, sigma );
    }
    s = s_STAR_STAR;

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
        s_MR_STAR = s_STAR_STAR;
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
inline void TallThresholded
( DistMatrix<F>& A, DistMatrix<Base<F>,VR,STAR>& s, DistMatrix<F>& V,
  Base<F> tol, bool relative )
{
    DEBUG_ONLY(CallStackEntry cse("svd::TallThresholded"))
    if( relative )
        TallRelativeThresholded( A, s, V, tol );
    else
        TallAbsoluteThresholded( A, s, V, tol );
}

template<typename F>
inline void
TallAbsoluteThresholded
( DistMatrix<F,VC,STAR>& A, 
  DistMatrix<Base<F>,STAR,STAR>& s, 
  DistMatrix<F,STAR,STAR>& V,
  Base<F> tol )
{
    DEBUG_ONLY(
        CallStackEntry cse("svd::TallAbsoluteThresholded");
        if( A.Height() < A.Width() )
            LogicError("A must be at least as tall as it is wide");
        if( tol < 0 )
            LogicError("negative threshold does not make sense");
    )
    typedef Base<F> Real;
    const Grid& g = A.Grid();
    const Int m = A.Height();
    const Int n = A.Width();
    const Real frobNorm = FrobeniusNorm( A );
    if( tol == Real(0) )
    {
        const Real eps = lapack::MachineEpsilon<Real>();
        tol = m*frobNorm*eps;
    }
    if( tol >= frobNorm )
    {
        A.Resize( m, 0 );        
        s.Resize( 0, 1 );
        V.Resize( n, 0 );
        return;
    }

    // C := A^H A
    DistMatrix<F,STAR,STAR> C(g);
    Zeros( C, n, n );
    Herk( LOWER, ADJOINT, F(1), A.LockedMatrix(), F(0), C.Matrix() );
    C.SumOver( A.ColComm() );

    // [V,Sigma^2] := eig(C), where each sigma > tol
    HermitianEig( LOWER, C, s, V, tol*tol, frobNorm*frobNorm );
    const int k = s.Height();
    
    // Sigma := sqrt(Sigma^2)
    for( Int i=0; i<k; ++i )
        s.SetLocal( i, 0, Sqrt(s.GetLocal(i,0)) );

    // Y := A V
    DistMatrix<F,VC,STAR> Y(g);
    Y.AlignWith( A );
    Zeros( Y, m, k );
    LocalGemm( NORMAL, NORMAL, F(1), A, V, F(0), Y );

    // Set each column of A to be the corresponding normalized column of Y
    // NOTE: A (potentially better) alternative would be to compute the norm of
    //       each column of A and normalize via it, as it might vary slightly
    //       from the corresponding computed singular value.
    A = Y;
    {
        const Int localHeight = A.LocalHeight();
        for( Int j=0; j<k; ++j )
        {
            const Real sigma = s.GetLocal( j, 0 );
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                A.SetLocal( iLoc, j, A.GetLocal(iLoc,j)/sigma );
        }
    }
}

template<typename F>
inline void
TallRelativeThresholded
( DistMatrix<F,VC,STAR>& A, 
  DistMatrix<Base<F>,STAR,STAR>& s, 
  DistMatrix<F,STAR,STAR>& V,
  Base<F> relTol )
{
    DEBUG_ONLY(
        CallStackEntry cse("svd::TallRelativeThresholded");
        if( A.Height() < A.Width() )
            LogicError("A must be at least as tall as it is wide");
        if( relTol < 0 )
            LogicError("negative threshold does not make sense");
    )
    typedef Base<F> Real;
    const Grid& g = A.Grid();
    const Int m = A.Height();
    const Int n = A.Width();

    // C := A^H A
    DistMatrix<F,STAR,STAR> C(g);
    Zeros( C, n, n );
    Herk( LOWER, ADJOINT, F(1), A.LockedMatrix(), F(0), C.Matrix() );
    C.SumOver( A.ColComm() );

    // [V,Sigma^2] := eig(C)
    HermitianEig( LOWER, C, s, V, DESCENDING );
    const Real twoNorm = Sqrt(MaxNorm(s));
    
    // Sigma := sqrt(Sigma^2), where each sigma > twoNorm*relTol
    for( Int i=0; i<n; ++i )
    {
        const Real sigma = Sqrt(s.GetLocal(i,0));
        if( sigma < relTol*twoNorm )
        {
            s.Resize( i, 1 );
            V.Resize( n, i );
            break;
        }
        else
            s.SetLocal( i, 0, sigma );
    }
    const int k = s.Height();

    // Y := A V
    DistMatrix<F,VC,STAR> Y(g);
    Y.AlignWith( A );
    Zeros( Y, m, k );
    LocalGemm( NORMAL, NORMAL, F(1), A, V, F(0), Y );

    // Set each column of A to be the corresponding normalized column of Y
    // NOTE: A (potentially better) alternative would be to compute the norm of
    //       each column of A and normalize via it, as it might vary slightly
    //       from the corresponding computed singular value.
    A = Y;
    {
        const Int localHeight = A.LocalHeight();
        for( Int j=0; j<k; ++j )
        {
            const Real sigma = s.GetLocal( j, 0 );
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                A.SetLocal( iLoc, j, A.GetLocal(iLoc,j)/sigma );
        }
    }
}

template<typename F>
void TallThresholded
( DistMatrix<F,VC,STAR>& A, 
  DistMatrix<Base<F>,STAR,STAR>& s, 
  DistMatrix<F,STAR,STAR>& V,
  Base<F> tol, bool relative )
{
    DEBUG_ONLY(CallStackEntry cse("svd::TallThresholded"))
    if( relative )
        TallRelativeThresholded( A, s, V, tol );
    else
        TallAbsoluteThresholded( A, s, V, tol );
}

template<typename F>
inline void
WideAbsoluteThresholded
( Matrix<F>& A, Matrix<Base<F>>& s, Matrix<F>& V, Base<F> tol )
{
    DEBUG_ONLY(
        CallStackEntry cse("svd::WideAbsoluteThresholded");
        if( A.Width() < A.Height() )
            LogicError("A must be at least as wide as it is tall");
        if( tol < 0 )
            LogicError("negative threshold does not make sense");
    )
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    const Real frobNorm = FrobeniusNorm( A );
    if( tol == Real(0) )
    {
        const Real eps = lapack::MachineEpsilon<Real>();
        tol = n*frobNorm*eps;
    }
    if( tol >= frobNorm )
    {
        A.Resize( m, 0 );        
        s.Resize( 0, 1 );
        V.Resize( n, 0 );
        return;
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
WideRelativeThresholded
( Matrix<F>& A, Matrix<Base<F>>& s, Matrix<F>& V, Base<F> relTol )
{
    DEBUG_ONLY(
        CallStackEntry cse("svd::WideThresholded");
        if( A.Width() < A.Height() )
            LogicError("A must be at least as wide as it is tall");
        if( relTol < 0 )
            LogicError("negative threshold does not make sense");
    )
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();

    // C := A A^H
    Matrix<F> C;
    Herk( LOWER, NORMAL, F(1), A, C );

    // [U,Sigma^2] := eig(C)
    Matrix<F> U;
    HermitianEig( LOWER, C, s, U, DESCENDING );
    const Real twoNorm = Sqrt(MaxNorm(s));
    
    // Sigma := sqrt(Sigma^2), where each sigma > relTol*twoNorm
    for( Int i=0; i<m; ++i )
    {
        const Real sigma = Sqrt(s.Get(i,0));
        if( sigma <= relTol*twoNorm )
        {
            s.Resize( i, 1 );
            U.Resize( m, i );
            break;
        }
        else
            s.Set( i, 0, sigma );
    }
    const Int k = s.Height();

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
inline void WideThresholded
( Matrix<F>& A, Matrix<Base<F>>& s, Matrix<F>& V, 
  Base<F> tol, bool relative )
{
    DEBUG_ONLY(CallStackEntry cse("svd::WideThresholded"))
    if( relative )
        WideRelativeThresholded( A, s, V, tol );
    else
        WideAbsoluteThresholded( A, s, V, tol );
}

template<typename F>
inline void
WideAbsoluteThresholded
( DistMatrix<F>& A, DistMatrix<Base<F>,VR,STAR>& s, DistMatrix<F>& V,
  Base<F> tol )
{
    DEBUG_ONLY(
        CallStackEntry cse("svd::WideAbsoluteThresholded");
        if( A.Width() < A.Height() )
            LogicError("A must be at least as wide as it is tall");
        if( tol < 0 )
            LogicError("negative threshold does not make sense");
    )
    typedef Base<F> Real;
    const Grid& g = A.Grid();
    const Int m = A.Height();
    const Int n = A.Width();
    const Real frobNorm = FrobeniusNorm( A );
    if( tol == Real(0) )
    {
        const Real eps = lapack::MachineEpsilon<Real>();
        tol = n*frobNorm*eps;
    }
    if( tol >= frobNorm )
    {
        A.Resize( m, 0 );        
        s.Resize( 0, 1 );
        V.Resize( n, 0 );
        return;
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
WideRelativeThresholded
( DistMatrix<F>& A, DistMatrix<Base<F>,VR,STAR>& s, DistMatrix<F>& V,
  Base<F> relTol )
{
    DEBUG_ONLY(
        CallStackEntry cse("svd::WideRelativeThresholded");
        if( A.Width() < A.Height() )
            LogicError("A must be at least as wide as it is tall");
        if( relTol < 0 )
            LogicError("negative threshold does not make sense");
    )
    typedef Base<F> Real;
    const Grid& g = A.Grid();
    const Int m = A.Height();

    // C := A A^H
    DistMatrix<F> C( g );
    Herk( LOWER, NORMAL, F(1), A, C );

    // [U,Sigma^2] := eig(C)
    DistMatrix<F> U(g);
    HermitianEig( LOWER, C, s, U, DESCENDING );
    const Real twoNorm = Sqrt(MaxNorm(s));
    
    // Sigma := sqrt(Sigma^2), where all sigmas > relTol*twoNorm
    DistMatrix<Real,STAR,STAR> s_STAR_STAR( s );
    for( Int i=0; i<m; ++i )
    {
        const Real sigma = Sqrt(s_STAR_STAR.GetLocal(i,0));
        if( sigma <= relTol*twoNorm )
        {
            s_STAR_STAR.Resize( i, 1 );
            U.Resize( m, i );
            break;
        }
        else
            s_STAR_STAR.SetLocal( i, 0, sigma );
    }
    s = s_STAR_STAR;

    // (Sigma V) := A^H U
    Gemm( ADJOINT, NORMAL, F(1), A, U, V );

    // Divide each column of (Sigma V) by sigma
    // NOTE: A (potentially better) alternative would be to compute the norm of
    //       each column of V and normalize via it, as it might vary slightly
    //       from the corresponding computed singular value.
    {
        DistMatrix<Real,MR,STAR> s_MR_STAR( g );
        s_MR_STAR.AlignWith( V.DistData() );
        s_MR_STAR = s_STAR_STAR;
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
inline void WideThresholded
( DistMatrix<F>& A, DistMatrix<Base<F>,VR,STAR>& s, DistMatrix<F>& V,
  Base<F> tol, bool relative )
{
    DEBUG_ONLY(CallStackEntry cse("svd::WideThresholded"))
    if( relative )
        WideRelativeThresholded( A, s, V, tol );
    else
        WideAbsoluteThresholded( A, s, V, tol );
}

// NOTE: [* ,VR] WideThresholded would produce U with different distribution
//       than A. It makes more sense to overwrite A with V'.

template<typename F>
void Thresholded
( Matrix<F>& A, Matrix<Base<F>>& s, Matrix<F>& V, 
  Base<F> tol, bool relative )
{
    DEBUG_ONLY(CallStackEntry cse("svd::Thresholded"))
    if( A.Height() >= A.Width() )
        TallThresholded( A, s, V, tol, relative );
    else
        WideThresholded( A, s, V, tol, relative );
}

template<typename F>
void Thresholded
( DistMatrix<F>& A, DistMatrix<Base<F>,VR,STAR>& s, DistMatrix<F>& V,
  Base<F> tol, bool relative )
{
    DEBUG_ONLY(CallStackEntry cse("svd::Thresholded"))
    if( A.Height() >= A.Width() )
        TallThresholded( A, s, V, tol, relative );
    else
        WideThresholded( A, s, V, tol, relative );
}

} // namespace svd
} // namespace El

#endif // ifndef EL_SVD_THRESHOLDED_HPP
