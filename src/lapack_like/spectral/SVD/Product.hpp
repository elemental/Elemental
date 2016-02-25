/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_SVD_PRODUCT_HPP
#define EL_SVD_PRODUCT_HPP

// TODO: Use a relative-truncated HermitianEig for relative thresholding

namespace El {
namespace svd {

// Compute singular triplets
// =========================

template<typename F>
inline void TallAbsoluteProduct
( const Matrix<F>& A,
        Matrix<F>& U, 
        Matrix<Base<F>>& s,
        Matrix<F>& V,
  Base<F> tol,
  bool avoidU )
{
    DEBUG_ONLY(
      CSE cse("svd::TallAbsoluteProduct");
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
        const Real eps = limits::Epsilon<Real>();
        tol = m*frobNorm*Sqrt(eps);
    }
    if( tol >= frobNorm )
    {
        U.Resize( m, 0 );        
        s.Resize( 0, 1 );
        V.Resize( n, 0 );
        return;
    }

    // C := A^H A
    Matrix<F> C;
    Herk( LOWER, ADJOINT, Real(1), A, C );

    // [V,Sigma^2] := eig(C), where each sigma > tol
    HermitianEigSubset<Real> subset;
    subset.rangeSubset = true;
    subset.lowerBound = tol*tol;
    // NOTE: While the square of the Frobenius norm should be a strict upper
    //       bound with exact computation, it has been observed that it can
    //       be lower than the finite-precision result in practice
    subset.upperBound = 2*frobNorm*frobNorm;
    HermitianEig( LOWER, C, s, V, DESCENDING, subset );
    
    // Sigma := sqrt(Sigma^2)
    const Int k = s.Height();
    for( Int i=0; i<k; ++i )
        s.Set( i, 0, Sqrt(s.Get(i,0)) );

    if( avoidU )
    {
        // Y := A V
        Matrix<F> Y;
        Gemm( NORMAL, NORMAL, F(1), A, V, Y );

        // Set each column of U to be the corresponding normalized column of Y
        U = Y;
        Matrix<Base<F>> colNorms;
        ColumnTwoNorms( U, colNorms );
        DiagonalSolve( RIGHT, NORMAL, colNorms, U );
    }
}

template<typename F>
inline void TallRelativeProduct
( const Matrix<F>& A,
        Matrix<F>& U,
        Matrix<Base<F>>& s,
        Matrix<F>& V,
  Base<F> relTol,
  bool avoidU )
{
    DEBUG_ONLY(
      CSE cse("svd::TallRelativeProduct");
      if( A.Height() < A.Width() )
          LogicError("A must be at least as tall as it is wide");
      if( relTol < 0 )
          LogicError("negative threshold does not make sense");
    )
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    if( relTol == Real(0) )
    {
        const Real eps = limits::Epsilon<Real>();
        relTol = m*Sqrt(eps);
    }

    // C := A^H A
    Matrix<F> C;
    Herk( LOWER, ADJOINT, Real(1), A, C );

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

    if( avoidU )
    {
        // Y := A V
        Matrix<F> Y;
        Gemm( NORMAL, NORMAL, F(1), A, V, Y );

        // Set each column of U to be the corresponding normalized column of Y
        U = Y;
        Matrix<Base<F>> colNorms;
        ColumnTwoNorms( U, colNorms );
        DiagonalSolve( RIGHT, NORMAL, colNorms, U );
    }
}

template<typename F>
inline void TallProduct
( const Matrix<F>& A,
        Matrix<F>& U,
        Matrix<Base<F>>& s,
        Matrix<F>& V, 
  Base<F> tol,
  bool relative,
  bool avoidU )
{
    DEBUG_ONLY(CSE cse("svd::TallProduct"))
    if( relative )
        TallRelativeProduct( A, U, s, V, tol, avoidU );
    else
        TallAbsoluteProduct( A, U, s, V, tol, avoidU );
}

template<typename F>
inline void
TallAbsoluteProduct
( const DistMatrix<F>& A,
        DistMatrix<F>& U,
        ElementalMatrix<Base<F>>& s, 
        DistMatrix<F>& V,
  Base<F> tol,
  bool avoidU )
{
    DEBUG_ONLY(
      CSE cse("svd::TallAbsoluteProduct");
      AssertSameGrids( A, U, s, V );
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
        const Real eps = limits::Epsilon<Real>();
        tol = m*frobNorm*Sqrt(eps);
    }
    if( tol >= frobNorm )
    {
        U.Resize( m, 0 );        
        s.Resize( 0, 1 );
        V.Resize( n, 0 );
        return;
    }

    // C := A^H A
    const Grid& g = A.Grid();
    DistMatrix<F> C(g);
    Herk( LOWER, ADJOINT, Real(1), A, C );

    // [V,Sigma^2] := eig(C), where each sigma > tol
    HermitianEigSubset<Real> subset;
    subset.rangeSubset = true;
    subset.lowerBound = tol*tol;
    // NOTE: While the square of the Frobenius norm should be a strict upper
    //       bound with exact computation, it has been observed that it can
    //       be lower than the finite-precision result in practice
    subset.upperBound = 2*frobNorm*frobNorm;
    HermitianEig( LOWER, C, s, V, DESCENDING, subset );
    
    // Sigma := sqrt(Sigma^2)
    {
        const Int localHeight = s.LocalHeight();
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
            s.SetLocal( iLoc, 0, Sqrt(s.GetLocal(iLoc,0)) );
    }

    if( avoidU )
    {
        // Y := A V
        DistMatrix<F> Y(g);
        Gemm( NORMAL, NORMAL, F(1), A, V, Y );

        // Set each column of U to be the corresponding normalized column of Y
        U = Y;
        DistMatrix<Real,MR,STAR> colNorms(g);
        ColumnTwoNorms( U, colNorms );
        DiagonalSolve( RIGHT, NORMAL, colNorms, U );
    }
}

template<typename F>
inline void
TallAbsoluteProduct
( const ElementalMatrix<F>& APre,
        ElementalMatrix<F>& UPre,
        ElementalMatrix<Base<F>>& s, 
        ElementalMatrix<F>& VPre,
  Base<F> tol,
  bool avoidU )
{
    DEBUG_ONLY(CSE cse("svd::TallAbsoluteProduct"))
    DistMatrixReadProxy<F,F,MC,MR> AProx( APre );
    DistMatrixWriteProxy<F,F,MC,MR> UProx( UPre );
    DistMatrixWriteProxy<F,F,MC,MR> VProx( VPre );
    auto& A = AProx.GetLocked();
    auto& U = UProx.Get();
    auto& V = VProx.Get();
    TallAbsoluteProduct( A, U, s, V, tol, avoidU );
}

template<typename F>
inline void
TallRelativeProduct
( const DistMatrix<F>& A,
        DistMatrix<F>& U,
        ElementalMatrix<Base<F>>& s, 
        DistMatrix<F>& V,
  Base<F> relTol,
  bool avoidU )
{
    DEBUG_ONLY(
      CSE cse("svd::TallRelativeProduct");
      AssertSameGrids( A, U, s, V );
      if( A.Height() < A.Width() )
          LogicError("A must be at least as tall as it is wide");
      if( relTol < 0 )
          LogicError("negative threshold does not make sense");
    )
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    if( relTol == Real(0) )
    {
        const Real eps = limits::Epsilon<Real>();
        relTol = m*Sqrt(eps);
    }

    // C := A^H A
    typedef Base<F> Real;
    const Grid& g = A.Grid();
    DistMatrix<F> C(g);
    Herk( LOWER, ADJOINT, Real(1), A, C );

    // [V,Sigma^2] := eig(C)
    HermitianEig( LOWER, C, s, V, DESCENDING );
    const Real twoNorm = Sqrt(MaxNorm(s));

    // Sigma := sqrt(Sigma^2), where all sigmas > relTol*twoNorm
    DistMatrix<Real,STAR,STAR> s_STAR_STAR( s );
    for( Int i=0; i<n; ++i )
    {
        const Real lambda = s_STAR_STAR.GetLocal(i,0);
        if( lambda <= Real(0) || Sqrt(lambda) <= relTol*twoNorm )
        {
            s_STAR_STAR.Resize( i, 1 );
            V.Resize( n, i );
            break;
        }
        else
            s_STAR_STAR.SetLocal( i, 0, Sqrt(lambda) );
    }
    Copy( s_STAR_STAR, s );

    if( avoidU )
    {
        // Y := A V
        DistMatrix<F> Y(g);
        Gemm( NORMAL, NORMAL, F(1), A, V, Y );

        // Set each column of U to be the corresponding normalized column of Y
        U = Y;
        DistMatrix<Real,MR,STAR> colNorms(g);
        ColumnTwoNorms( U, colNorms );
        DiagonalSolve( RIGHT, NORMAL, colNorms, U );
    }
}

template<typename F>
inline void
TallRelativeProduct
( const ElementalMatrix<F>& APre,
        ElementalMatrix<F>& UPre,
        ElementalMatrix<Base<F>>& s, 
        ElementalMatrix<F>& VPre,
  Base<F> relTol,
  bool avoidU )
{
    DEBUG_ONLY(CSE cse("svd::TallRelativeProduct"))
    DistMatrixReadProxy<F,F,MC,MR> AProx( APre );
    DistMatrixWriteProxy<F,F,MC,MR> UProx( UPre );
    DistMatrixWriteProxy<F,F,MC,MR> VProx( VPre );
    auto& A = AProx.GetLocked();
    auto& U = UProx.Get();
    auto& V = VProx.Get();
    TallRelativeProduct( A, U, s, V, relTol, avoidU );
}

template<typename F>
inline void TallProduct
( const ElementalMatrix<F>& A,
        ElementalMatrix<F>& U,
        ElementalMatrix<Base<F>>& s, 
        ElementalMatrix<F>& V,
  Base<F> tol,
  bool relative,
  bool avoidU )
{
    DEBUG_ONLY(CSE cse("svd::TallProduct"))
    if( relative )
        TallRelativeProduct( A, U, s, V, tol, avoidU );
    else
        TallAbsoluteProduct( A, U, s, V, tol, avoidU );
}

template<typename F>
inline void
TallAbsoluteProduct
( const DistMatrix<F,VC,STAR>& A,
        DistMatrix<F,VC,STAR>& U,
        DistMatrix<Base<F>,STAR,STAR>& s, 
        DistMatrix<F,STAR,STAR>& V,
  Base<F> tol,
  bool avoidU )
{
    DEBUG_ONLY(
      CSE cse("svd::TallAbsoluteProduct");
      AssertSameGrids( A, U, s, V );
      if( A.Height() < A.Width() )
          LogicError("A must be at least as tall as it is wide");
      if( tol < 0 )
          LogicError("negative threshold does not make sense");
    )
    const Int m = A.Height();
    const Int n = A.Width();

    typedef Base<F> Real;
    const Real frobNorm = FrobeniusNorm( A );
    if( tol == Real(0) )
    {
        const Real eps = limits::Epsilon<Real>();
        tol = m*frobNorm*Sqrt(eps);
    }
    if( tol >= frobNorm )
    {
        U.Resize( m, 0 );        
        s.Resize( 0, 1 );
        V.Resize( n, 0 );
        return;
    }

    // C := A^H A
    const Grid& g = A.Grid();
    DistMatrix<F,STAR,STAR> C(g);
    Zeros( C, n, n );
    Herk( LOWER, ADJOINT, Real(1), A.LockedMatrix(), Real(0), C.Matrix() );
    El::AllReduce( C, A.ColComm() );

    // [V,Sigma^2] := eig(C), where each sigma > tol
    HermitianEigSubset<Real> subset;
    subset.rangeSubset = true;
    subset.lowerBound = tol*tol;
    // NOTE: While the square of the Frobenius norm should be a strict upper
    //       bound with exact computation, it has been observed that it can
    //       be lower than the finite-precision result in practice
    subset.upperBound = 2*frobNorm*frobNorm;
    HermitianEig( LOWER, C, s, V, DESCENDING, subset );
    const int k = s.Height();
    
    // Sigma := sqrt(Sigma^2)
    for( Int i=0; i<k; ++i )
        s.SetLocal( i, 0, Sqrt(s.GetLocal(i,0)) );

    if( avoidU )
    {
        // Y := A V
        DistMatrix<F,VC,STAR> Y(g);
        Y.AlignWith( A );
        Zeros( Y, m, k );
        LocalGemm( NORMAL, NORMAL, F(1), A, V, F(0), Y );

        // Set each column of U to be the corresponding normalized column of Y
        U = Y;
        DistMatrix<Real,STAR,STAR> colNorms(g);
        ColumnTwoNorms( U, colNorms );
        DiagonalSolve( RIGHT, NORMAL, colNorms, U );
    }
}

template<typename F>
inline void
TallAbsoluteProduct
( const DistMatrix<F,VC,STAR>& A,
        ElementalMatrix<F>& UPre,
        ElementalMatrix<Base<F>>& sPre, 
        ElementalMatrix<F>& VPre,
  Base<F> tol,
  bool avoidU )
{
    DEBUG_ONLY(CSE cse("svd::TallAbsoluteProduct"))
    DistMatrixWriteProxy<F,F,VC,STAR> UProx( UPre );
    DistMatrixWriteProxy<Base<F>,Base<F>,STAR,STAR> sProx( sPre );
    DistMatrixWriteProxy<F,F,STAR,STAR> VProx( VPre );
    auto& s = sProx.Get();
    auto& U = UProx.Get();
    auto& V = VProx.Get();
    TallAbsoluteProduct( A, U, s, V, tol, avoidU );
}

template<typename F>
inline void
TallRelativeProduct
( const DistMatrix<F,VC,STAR>& A,
        DistMatrix<F,VC,STAR>& U,
        DistMatrix<Base<F>,STAR,STAR>& s, 
        DistMatrix<F,STAR,STAR>& V,
  Base<F> relTol,
  bool avoidU )
{
    DEBUG_ONLY(
      CSE cse("svd::TallRelativeProduct");
      AssertSameGrids( A, U, s, V );
      if( A.Height() < A.Width() )
          LogicError("A must be at least as tall as it is wide");
      if( relTol < 0 )
          LogicError("negative threshold does not make sense");
    )
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    if( relTol == Real(0) )
    {
        const Real eps = limits::Epsilon<Real>();
        relTol = m*Sqrt(eps);
    }

    // C := A^H A
    typedef Base<F> Real;
    const Grid& g = A.Grid();
    DistMatrix<F,STAR,STAR> C(g);
    Zeros( C, n, n );
    Herk( LOWER, ADJOINT, Real(1), A.LockedMatrix(), Real(0), C.Matrix() );
    El::AllReduce( C, A.ColComm() );

    // [V,Sigma^2] := eig(C)
    HermitianEig( LOWER, C, s, V, DESCENDING );
    const Real twoNorm = Sqrt(MaxNorm(s));
    
    // Sigma := sqrt(Sigma^2), where each sigma > twoNorm*relTol
    for( Int i=0; i<n; ++i )
    {
        const Real lambda = s.GetLocal(i,0);
        if( lambda <= Real(0) || Sqrt(lambda) <= relTol*twoNorm )
        {
            s.Resize( i, 1 );
            V.Resize( n, i );
            break;
        }
        else
            s.SetLocal( i, 0, Sqrt(lambda) );
    }
    const int k = s.Height();

    if( avoidU )
    {
        // Y := A V
        DistMatrix<F,VC,STAR> Y(g);
        Y.AlignWith( A );
        Zeros( Y, m, k );
        LocalGemm( NORMAL, NORMAL, F(1), A, V, F(0), Y );

        // Set each column of U to be the corresponding normalized column of Y
        U = Y;
        DistMatrix<Real,STAR,STAR> colNorms(g);
        ColumnTwoNorms( U, colNorms );
        DiagonalSolve( RIGHT, NORMAL, colNorms, U );
    }
}

template<typename F>
inline void
TallRelativeProduct
( const DistMatrix<F,VC,STAR>& A,
        ElementalMatrix<F>& UPre,
        ElementalMatrix<Base<F>>& sPre, 
        ElementalMatrix<F>& VPre,
  Base<F> relTol,
  bool avoidU )
{
    DEBUG_ONLY(CSE cse("svd::TallRelativeProduct"))
    DistMatrixWriteProxy<Base<F>,Base<F>,STAR,STAR> sProx( sPre );
    DistMatrixWriteProxy<F,F,VC,STAR> UProx( UPre );
    DistMatrixWriteProxy<F,F,STAR,STAR> VProx( VPre );
    auto& s = sProx.Get();
    auto& U = UProx.Get();
    auto& V = VProx.Get();
    TallRelativeProduct( A, U, s, V, relTol, avoidU );
}

template<typename F>
void TallProduct
( const DistMatrix<F,VC,STAR>& A,
        ElementalMatrix<F>& U,
        ElementalMatrix<Base<F>>& s, 
        ElementalMatrix<F>& V,
  Base<F> tol,
  bool relative,
  bool avoidU )
{
    DEBUG_ONLY(CSE cse("svd::TallProduct"))
    if( relative )
        TallRelativeProduct( A, U, s, V, tol, avoidU );
    else
        TallAbsoluteProduct( A, U, s, V, tol, avoidU );
}

template<typename F>
inline void
WideAbsoluteProduct
( const Matrix<F>& A,
        Matrix<F>& U,
        Matrix<Base<F>>& s,
        Matrix<F>& V,
  Base<F> tol,
  bool avoidV )
{
    DEBUG_ONLY(
      CSE cse("svd::WideAbsoluteProduct");
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
        const Real eps = limits::Epsilon<Real>();
        tol = n*frobNorm*Sqrt(eps);
    }
    if( tol >= frobNorm )
    {
        U.Resize( m, 0 );        
        s.Resize( 0, 1 );
        V.Resize( n, 0 );
        return;
    }

    // C := A A^H
    Matrix<F> C;
    Herk( LOWER, NORMAL, Real(1), A, C );

    // [U,Sigma^2] := eig(C), where each sigma > tol
    HermitianEigSubset<Real> subset;
    subset.rangeSubset = true;
    subset.lowerBound = tol*tol;
    // NOTE: While the square of the Frobenius norm should be a strict upper
    //       bound with exact computation, it has been observed that it can
    //       be lower than the finite-precision result in practice
    subset.upperBound = 2*frobNorm*frobNorm;
    HermitianEig( LOWER, C, s, U, DESCENDING, subset );
    
    // Sigma := sqrt(Sigma^2)
    const Int k = s.Height();
    for( Int i=0; i<k; ++i )
        s.Set( i, 0, Sqrt(s.Get(i,0)) );

    if( avoidV )
    {
        // (Sigma V) := A^H U
        Gemm( ADJOINT, NORMAL, F(1), A, U, V );

        // Normalize each column of Sigma V
        Matrix<Base<F>> colNorms;
        ColumnTwoNorms( V, colNorms );
        DiagonalSolve( RIGHT, NORMAL, colNorms, V );
    }
}

template<typename F>
inline void
WideRelativeProduct
( const Matrix<F>& A,
        Matrix<F>& U,
        Matrix<Base<F>>& s,
        Matrix<F>& V,
  Base<F> relTol,
  bool avoidV )
{
    DEBUG_ONLY(
      CSE cse("svd::WideProduct");
      if( A.Width() < A.Height() )
          LogicError("A must be at least as wide as it is tall");
      if( relTol < 0 )
          LogicError("negative threshold does not make sense");
    )
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    if( relTol == Real(0) )
    {
        const Real eps = limits::Epsilon<Real>();
        relTol = n*Sqrt(eps);
    }

    // C := A A^H
    Matrix<F> C;
    Herk( LOWER, NORMAL, Real(1), A, C );

    // [U,Sigma^2] := eig(C)
    HermitianEig( LOWER, C, s, U, DESCENDING );
    const Real twoNorm = Sqrt(MaxNorm(s));
    
    // Sigma := sqrt(Sigma^2), where each sigma > relTol*twoNorm
    for( Int i=0; i<m; ++i )
    {
        const Real lambda = s.Get(i,0);
        if( lambda <= Real(0) || Sqrt(lambda) <= relTol*twoNorm )
        {
            s.Resize( i, 1 );
            U.Resize( m, i );
            break;
        }
        else
            s.Set( i, 0, Sqrt(lambda) );
    }

    if( avoidV )
    {
        // (Sigma V) := A^H U
        Gemm( ADJOINT, NORMAL, F(1), A, U, V );

        // Normalize each column of Sigma V
        Matrix<Base<F>> colNorms;
        ColumnTwoNorms( V, colNorms );
        DiagonalSolve( RIGHT, NORMAL, colNorms, V );
    }
}

template<typename F>
inline void WideProduct
( const Matrix<F>& A,
        Matrix<F>& U,
        Matrix<Base<F>>& s,
        Matrix<F>& V, 
  Base<F> tol,
  bool relative,
  bool avoidV )
{
    DEBUG_ONLY(CSE cse("svd::WideProduct"))
    if( relative )
        WideRelativeProduct( A, U, s, V, tol, avoidV );
    else
        WideAbsoluteProduct( A, U, s, V, tol, avoidV );
}

template<typename F>
inline void
WideAbsoluteProduct
( const DistMatrix<F>& A,
        DistMatrix<F>& U,
        ElementalMatrix<Base<F>>& s, 
        DistMatrix<F>& V,
  Base<F> tol,
  bool avoidV )
{
    DEBUG_ONLY(
      CSE cse("svd::WideAbsoluteProduct");
      if( A.Width() < A.Height() )
          LogicError("A must be at least as wide as it is tall");
      if( tol < 0 )
          LogicError("negative threshold does not make sense");
    )
    const Int m = A.Height();
    const Int n = A.Width();

    typedef Base<F> Real;
    const Real frobNorm = FrobeniusNorm( A );
    if( tol == Real(0) )
    {
        const Real eps = limits::Epsilon<Real>();
        tol = n*frobNorm*Sqrt(eps);
    }
    if( tol >= frobNorm )
    {
        U.Resize( m, 0 );        
        s.Resize( 0, 1 );
        V.Resize( n, 0 );
        return;
    }

    // C := A A^H
    const Grid& g = A.Grid();
    DistMatrix<F> C( g );
    Herk( LOWER, NORMAL, Real(1), A, C );

    // [U,Sigma^2] := eig(C), where each sigma > tol
    HermitianEigSubset<Real> subset;
    subset.rangeSubset = true;
    subset.lowerBound = tol*tol;
    // NOTE: While the square of the Frobenius norm should be a strict upper
    //       bound with exact computation, it has been observed that it can
    //       be lower than the finite-precision result in practice
    subset.upperBound = 2*frobNorm*frobNorm;
    HermitianEig( LOWER, C, s, U, DESCENDING, subset );
    
    // Sigma := sqrt(Sigma^2)
    {
        const Int localHeight = s.LocalHeight();
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
            s.SetLocal( iLoc, 0, Sqrt(s.GetLocal(iLoc,0)) );
    }

    if( avoidV )
    {
        // (Sigma V) := A^H U
        Gemm( ADJOINT, NORMAL, F(1), A, U, V );

        // Normalize each column of Sigma V
        DistMatrix<Real,MR,STAR> colNorms(g);
        ColumnTwoNorms( V, colNorms );
        DiagonalSolve( RIGHT, NORMAL, colNorms, V );
    }
}

template<typename F>
inline void
WideAbsoluteProduct
( const ElementalMatrix<F>& APre,
        ElementalMatrix<F>& UPre,
        ElementalMatrix<Base<F>>& s, 
        ElementalMatrix<F>& VPre,
  Base<F> tol,
  bool avoidV )
{
    DEBUG_ONLY(CSE cse("svd::WideAbsoluteProduct"))
    DistMatrixReadProxy<F,F,MC,MR> AProx( APre );
    DistMatrixWriteProxy<F,F,MC,MR> UProx( UPre );
    DistMatrixWriteProxy<F,F,MC,MR> VProx( VPre );
    auto& A = AProx.GetLocked();
    auto& U = UProx.Get();
    auto& V = VProx.Get();
    WideAbsoluteProduct( A, U, s, V, tol, avoidV );
}

template<typename F>
inline void
WideRelativeProduct
( const DistMatrix<F>& A,
        DistMatrix<F>& U,
        ElementalMatrix<Base<F>>& s, 
        DistMatrix<F>& V,
  Base<F> relTol,
  bool avoidV )
{
    DEBUG_ONLY(
      CSE cse("svd::WideRelativeProduct");
      AssertSameGrids( A, U, s, V );
      if( A.Width() < A.Height() )
          LogicError("A must be at least as wide as it is tall");
      if( relTol < 0 )
          LogicError("negative threshold does not make sense");
    )
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    if( relTol == Real(0) )
    {
        const Real eps = limits::Epsilon<Real>();
        relTol = n*Sqrt(eps);
    }

    // C := A A^H
    typedef Base<F> Real;
    const Grid& g = A.Grid();
    DistMatrix<F> C( g );
    Herk( LOWER, NORMAL, Real(1), A, C );

    // [U,Sigma^2] := eig(C)
    HermitianEig( LOWER, C, s, U, DESCENDING );
    const Real twoNorm = Sqrt(MaxNorm(s));
    
    // Sigma := sqrt(Sigma^2), where all sigmas > relTol*twoNorm
    DistMatrix<Real,STAR,STAR> s_STAR_STAR( s );
    for( Int i=0; i<m; ++i )
    {
        const Real lambda = s_STAR_STAR.GetLocal(i,0);
        if( lambda <= Real(0) || Sqrt(lambda) <= relTol*twoNorm )
        {
            s_STAR_STAR.Resize( i, 1 );
            U.Resize( m, i );
            break;
        }
        else
            s_STAR_STAR.SetLocal( i, 0, Sqrt(lambda) );
    }
    Copy( s_STAR_STAR, s );

    if( avoidV )
    {
        // (Sigma V) := A^H U
        Gemm( ADJOINT, NORMAL, F(1), A, U, V );

        // Normalize each column of Sigma V
        DistMatrix<Real,MR,STAR> colNorms(g);
        ColumnTwoNorms( V, colNorms );
        DiagonalSolve( RIGHT, NORMAL, colNorms, V );
    }
}

template<typename F>
inline void
WideRelativeProduct
( const ElementalMatrix<F>& APre,
        ElementalMatrix<F>& UPre,
        ElementalMatrix<Base<F>>& s, 
        ElementalMatrix<F>& VPre,
  Base<F> relTol,
  bool avoidV )
{
    DEBUG_ONLY(CSE cse("svd::WideRelativeProduct"))
    DistMatrixReadProxy<F,F,MC,MR> AProx( APre );
    DistMatrixWriteProxy<F,F,MC,MR> UProx( UPre );
    DistMatrixWriteProxy<F,F,MC,MR> VProx( VPre );
    auto& A = AProx.GetLocked();
    auto& U = UProx.Get();
    auto& V = VProx.Get();
    WideRelativeProduct( A, U, s, V, relTol, avoidV );
}

template<typename F>
inline void WideProduct
( const ElementalMatrix<F>& A,
        ElementalMatrix<F>& U,
        ElementalMatrix<Base<F>>& s, 
        ElementalMatrix<F>& V,
  Base<F> tol,
  bool relative,
  bool avoidV )
{
    DEBUG_ONLY(CSE cse("svd::WideProduct"))
    if( relative )
        WideRelativeProduct( A, U, s, V, tol, avoidV );
    else
        WideAbsoluteProduct( A, U, s, V, tol, avoidV );
}

// NOTE: [* ,VR] WideProduct would produce U with different distribution
//       than A. It makes more sense to overwrite A with V'.
// TODO: Update the above note and the following routines now that A is not
//       overwritten

template<typename F>
void Product
( const Matrix<F>& A,
        Matrix<F>& U, 
        Matrix<Base<F>>& s,
        Matrix<F>& V, 
  Base<F> tol,
  bool relative,
  bool avoidU,
  bool avoidV )
{
    DEBUG_ONLY(CSE cse("svd::Product"))
    // TODO: If m(A) >=~ n(A) but avoidV requested, find a way to make use of it
    if( A.Height() >= A.Width() )
        TallProduct( A, U, s, V, tol, relative, avoidU );
    else
        WideProduct( A, U, s, V, tol, relative, avoidV );
}

template<typename F>
void Product
( const ElementalMatrix<F>& A,
        ElementalMatrix<F>& U,
        ElementalMatrix<Base<F>>& s, 
        ElementalMatrix<F>& V,
  Base<F> tol,
  bool relative,
  bool avoidU,
  bool avoidV )
{
    DEBUG_ONLY(CSE cse("svd::Product"))
    // TODO: If m(A) >=~ n(A) but avoidV requested, find a way to make use of it
    if( A.Height() >= A.Width() )
        TallProduct( A, U, s, V, tol, relative, avoidU );
    else
        WideProduct( A, U, s, V, tol, relative, avoidV );
}

// Compute singular values
// =======================

template<typename F>
inline void TallAbsoluteProduct
( const Matrix<F>& A,
        Matrix<Base<F>>& s,
  Base<F> tol )
{
    DEBUG_ONLY(
      CSE cse("svd::TallAbsoluteProduct");
      if( A.Height() < A.Width() )
          LogicError("A must be at least as tall as it is wide");
      if( tol < 0 )
          LogicError("negative threshold does not make sense");
    )
    typedef Base<F> Real;
    const Int m = A.Height();
    const Real frobNorm = FrobeniusNorm( A );
    if( tol == Real(0) )
    {
        const Real eps = limits::Epsilon<Real>();
        tol = m*frobNorm*Sqrt(eps);
    }
    if( tol >= frobNorm )
    {
        s.Resize( 0, 1 );
        return;
    }

    // C := A^H A
    Matrix<F> C;
    Herk( LOWER, ADJOINT, Real(1), A, C );

    // [Sigma^2] := eig(C), where each sigma > tol
    HermitianEigSubset<Real> subset;
    subset.rangeSubset = true;
    subset.lowerBound = tol*tol;
    // NOTE: While the square of the Frobenius norm should be a strict upper
    //       bound with exact computation, it has been observed that it can
    //       be lower than the finite-precision result in practice
    subset.upperBound = 2*frobNorm*frobNorm;
    HermitianEig( LOWER, C, s, DESCENDING, subset );
    
    // Sigma := sqrt(Sigma^2)
    const Int k = s.Height();
    for( Int i=0; i<k; ++i )
        s.Set( i, 0, Sqrt(s.Get(i,0)) );
}

template<typename F>
inline void TallRelativeProduct
( const Matrix<F>& A,
        Matrix<Base<F>>& s,
  Base<F> relTol )
{
    DEBUG_ONLY(
      CSE cse("svd::TallRelativeProduct");
      if( A.Height() < A.Width() )
          LogicError("A must be at least as tall as it is wide");
      if( relTol < 0 )
          LogicError("negative threshold does not make sense");
    )
    typedef Base<F> Real;
    const Int n = A.Width();

    // C := A^H A
    Matrix<F> C;
    Herk( LOWER, ADJOINT, Real(1), A, C );

    // [Sigma^2] := eig(C)
    HermitianEig( LOWER, C, s, DESCENDING );
    const Real twoNorm = Sqrt(MaxNorm(s));
    
    // Sigma := sqrt(Sigma^2), where all sigmas > relTol*twoNorm
    for( Int i=0; i<n; ++i )
    {
        const Real sigma = Sqrt(s.Get(i,0));
        if( sigma <= relTol*twoNorm )
        {
            s.Resize( i, 1 );
            break;
        }
        else
            s.Set( i, 0, sigma );
    }
}

template<typename F>
inline void TallProduct
( const Matrix<F>& A,
        Matrix<Base<F>>& s,
  Base<F> tol,
  bool relative )
{
    DEBUG_ONLY(CSE cse("svd::TallProduct"))
    if( relative )
        TallRelativeProduct( A, s, tol );
    else
        TallAbsoluteProduct( A, s, tol );
}

template<typename F>
inline void
TallAbsoluteProduct
( const DistMatrix<F>& A,
        ElementalMatrix<Base<F>>& s, 
  Base<F> tol )
{
    DEBUG_ONLY(
      CSE cse("svd::TallAbsoluteProduct");
      AssertSameGrids( A, s );
      if( A.Height() < A.Width() )
          LogicError("A must be at least as tall as it is wide");
      if( tol < 0 )
          LogicError("negative threshold does not make sense");
    )
    typedef Base<F> Real;
    const Int m = A.Height();
    const Real frobNorm = FrobeniusNorm( A );
    if( tol == Real(0) )
    {
        const Real eps = limits::Epsilon<Real>();
        tol = m*frobNorm*Sqrt(eps);
    }
    if( tol >= frobNorm )
    {
        s.Resize( 0, 1 );
        return;
    }

    // C := A^H A
    const Grid& g = A.Grid();
    DistMatrix<F> C(g);
    Herk( LOWER, ADJOINT, Real(1), A, C );

    // [Sigma^2] := eig(C), where each sigma > tol
    HermitianEigSubset<Real> subset;
    subset.rangeSubset = true;
    subset.lowerBound = tol*tol;
    // NOTE: While the square of the Frobenius norm should be a strict upper
    //       bound with exact computation, it has been observed that it can
    //       be lower than the finite-precision result in practice
    subset.upperBound = 2*frobNorm*frobNorm;
    HermitianEig( LOWER, C, s, DESCENDING, subset );
    
    // Sigma := sqrt(Sigma^2)
    const Int localHeight = s.LocalHeight();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        s.SetLocal( iLoc, 0, Sqrt(s.GetLocal(iLoc,0)) );
}

template<typename F>
inline void
TallAbsoluteProduct
( const ElementalMatrix<F>& APre,
        ElementalMatrix<Base<F>>& s, 
  Base<F> tol )
{
    DEBUG_ONLY(CSE cse("svd::TallAbsoluteProduct"))
    DistMatrixReadProxy<F,F,MC,MR> AProx( APre );
    auto& A = AProx.GetLocked();
    TallAbsoluteProduct( A, s, tol );
}

template<typename F>
inline void
TallRelativeProduct
( const DistMatrix<F>& A,
        ElementalMatrix<Base<F>>& s, 
  Base<F> relTol )
{
    DEBUG_ONLY(
      CSE cse("svd::TallRelativeProduct");
      AssertSameGrids( A, s );
      if( A.Height() < A.Width() )
          LogicError("A must be at least as tall as it is wide");
      if( relTol < 0 )
          LogicError("negative threshold does not make sense");
    )
    const Int n = A.Width();

    // C := A^H A
    typedef Base<F> Real;
    const Grid& g = A.Grid();
    DistMatrix<F> C(g);
    Herk( LOWER, ADJOINT, Real(1), A, C );

    // [Sigma^2] := eig(C)
    HermitianEig( LOWER, C, s, DESCENDING );
    const Real twoNorm = Sqrt(MaxNorm(s));

    // Sigma := sqrt(Sigma^2), where all sigmas > relTol*twoNorm
    DistMatrix<Real,STAR,STAR> s_STAR_STAR( s );
    for( Int i=0; i<n; ++i )
    {
        const Real lambda = s_STAR_STAR.GetLocal(i,0);
        if( lambda <= Real(0) || Sqrt(lambda) <= relTol*twoNorm )
        {
            s_STAR_STAR.Resize( i, 1 );
            break;
        }
        else
            s_STAR_STAR.SetLocal( i, 0, Sqrt(lambda) );
    }
    Copy( s_STAR_STAR, s );
}

template<typename F>
inline void
TallRelativeProduct
( const ElementalMatrix<F>& APre,
        ElementalMatrix<Base<F>>& s, 
  Base<F> relTol )
{
    DEBUG_ONLY(CSE cse("svd::TallRelativeProduct"))
    DistMatrixReadProxy<F,F,MC,MR> AProx( APre );
    auto& A = AProx.GetLocked();
    TallRelativeProduct( A, s, relTol );
}

template<typename F>
inline void TallProduct
( const ElementalMatrix<F>& A,
        ElementalMatrix<Base<F>>& s, 
  Base<F> tol,
  bool relative )
{
    DEBUG_ONLY(CSE cse("svd::TallProduct"))
    if( relative )
        TallRelativeProduct( A, s, tol );
    else
        TallAbsoluteProduct( A, s, tol );
}

template<typename F>
inline void
TallAbsoluteProduct
( const DistMatrix<F,VC,STAR>& A,
        DistMatrix<Base<F>,STAR,STAR>& s, 
  Base<F> tol )
{
    DEBUG_ONLY(
      CSE cse("svd::TallAbsoluteProduct");
      AssertSameGrids( A, s );
      if( A.Height() < A.Width() )
          LogicError("A must be at least as tall as it is wide");
      if( tol < 0 )
          LogicError("negative threshold does not make sense");
    )
    const Int m = A.Height();
    const Int n = A.Width();

    typedef Base<F> Real;
    const Real frobNorm = FrobeniusNorm( A );
    if( tol == Real(0) )
    {
        const Real eps = limits::Epsilon<Real>();
        tol = m*frobNorm*Sqrt(eps);
    }
    if( tol >= frobNorm )
    {
        s.Resize( 0, 1 );
        return;
    }

    // C := A^H A
    const Grid& g = A.Grid();
    DistMatrix<F,STAR,STAR> C(g);
    Zeros( C, n, n );
    Herk( LOWER, ADJOINT, Real(1), A.LockedMatrix(), Real(0), C.Matrix() );
    El::AllReduce( C, A.ColComm() );

    // [Sigma^2] := eig(C), where each sigma > tol
    HermitianEigSubset<Real> subset;
    subset.rangeSubset = true;
    subset.lowerBound = tol*tol;
    // NOTE: While the square of the Frobenius norm should be a strict upper
    //       bound with exact computation, it has been observed that it can
    //       be lower than the finite-precision result in practice
    subset.upperBound = 2*frobNorm*frobNorm;
    HermitianEig( LOWER, C, s, DESCENDING, subset );
    const int k = s.Height();
    
    // Sigma := sqrt(Sigma^2)
    for( Int i=0; i<k; ++i )
        s.SetLocal( i, 0, Sqrt(s.GetLocal(i,0)) );
}

template<typename F>
inline void
TallAbsoluteProduct
( const DistMatrix<F,VC,STAR>& A,
        ElementalMatrix<Base<F>>& sPre, 
  Base<F> tol )
{
    DEBUG_ONLY(CSE cse("svd::TallAbsoluteProduct"))
    typedef Base<F> Real;
    DistMatrixWriteProxy<Real,Real,STAR,STAR> sProx( sPre );
    auto& s = sProx.Get();
    TallAbsoluteProduct( A, s, tol );
}

template<typename F>
inline void
TallRelativeProduct
( const DistMatrix<F,VC,STAR>& A,
        DistMatrix<Base<F>,STAR,STAR>& s, 
  Base<F> relTol )
{
    DEBUG_ONLY(
      CSE cse("svd::TallRelativeProduct");
      AssertSameGrids( A, s );
      if( A.Height() < A.Width() )
          LogicError("A must be at least as tall as it is wide");
      if( relTol < 0 )
          LogicError("negative threshold does not make sense");
    )
    const Int m = A.Height();
    const Int n = A.Width();

    // C := A^H A
    typedef Base<F> Real;
    const Grid& g = A.Grid();
    DistMatrix<F,STAR,STAR> C(g);
    Zeros( C, n, n );
    Herk( LOWER, ADJOINT, Real(1), A.LockedMatrix(), Real(0), C.Matrix() );
    El::AllReduce( C, A.ColComm() );

    // [V,Sigma^2] := eig(C)
    HermitianEig( LOWER, C, s, DESCENDING );
    const Real twoNorm = Sqrt(MaxNorm(s));
    
    // Sigma := sqrt(Sigma^2), where each sigma > twoNorm*relTol
    for( Int i=0; i<n; ++i )
    {
        const Real lambda = s.GetLocal(i,0);
        if( lambda <= Real(0) || Sqrt(lambda) <= relTol*twoNorm )
        {
            s.Resize( i, 1 );
            break;
        }
        else
            s.SetLocal( i, 0, Sqrt(lambda) );
    }
    const int k = s.Height();
}

template<typename F>
inline void
TallRelativeProduct
( const DistMatrix<F,VC,STAR>& A,
        ElementalMatrix<Base<F>>& sPre, 
  Base<F> relTol )
{
    DEBUG_ONLY(CSE cse("svd::TallRelativeProduct"))
    typedef Base<F> Real;
    DistMatrixWriteProxy<Real,Real,STAR,STAR> sProx( sPre );
    auto& s = sProx.Get();
    TallRelativeProduct( A, s, relTol );
}

template<typename F>
void TallProduct
( const DistMatrix<F,VC,STAR>& A,
        ElementalMatrix<Base<F>>& s, 
  Base<F> tol,
  bool relative )
{
    DEBUG_ONLY(CSE cse("svd::TallProduct"))
    if( relative )
        TallRelativeProduct( A, s, tol );
    else
        TallAbsoluteProduct( A, s, tol );
}

template<typename F>
inline void
WideAbsoluteProduct
( const Matrix<F>& A,
        Matrix<Base<F>>& s,
  Base<F> tol )
{
    DEBUG_ONLY(
      CSE cse("svd::WideAbsoluteProduct");
      if( A.Width() < A.Height() )
          LogicError("A must be at least as wide as it is tall");
      if( tol < 0 )
          LogicError("negative threshold does not make sense");
    )
    typedef Base<F> Real;
    const Int n = A.Width();
    const Real frobNorm = FrobeniusNorm( A );
    if( tol == Real(0) )
    {
        const Real eps = limits::Epsilon<Real>();
        tol = n*frobNorm*Sqrt(eps);
    }
    if( tol >= frobNorm )
    {
        s.Resize( 0, 1 );
        return;
    }

    // C := A A^H
    Matrix<F> C;
    Herk( LOWER, NORMAL, Real(1), A, C );

    // [Sigma^2] := eig(C), where each sigma > tol
    HermitianEigSubset<Real> subset;
    subset.rangeSubset = true;
    subset.lowerBound = tol*tol;
    // NOTE: While the square of the Frobenius norm should be a strict upper
    //       bound with exact computation, it has been observed that it can
    //       be lower than the finite-precision result in practice
    subset.upperBound = 2*frobNorm*frobNorm;
    HermitianEig( LOWER, C, s, DESCENDING, subset );
    
    // Sigma := sqrt(Sigma^2)
    const Int k = s.Height();
    for( Int i=0; i<k; ++i )
        s.Set( i, 0, Sqrt(s.Get(i,0)) );
}

template<typename F>
inline void
WideRelativeProduct
( const Matrix<F>& A,
        Matrix<Base<F>>& s,
  Base<F> relTol )
{
    DEBUG_ONLY(
      CSE cse("svd::WideProduct");
      if( A.Width() < A.Height() )
          LogicError("A must be at least as wide as it is tall");
      if( relTol < 0 )
          LogicError("negative threshold does not make sense");
    )
    typedef Base<F> Real;
    const Int m = A.Height();

    // C := A A^H
    Matrix<F> C;
    Herk( LOWER, NORMAL, Real(1), A, C );

    // [Sigma^2] := eig(C)
    HermitianEig( LOWER, C, s, DESCENDING );
    const Real twoNorm = Sqrt(MaxNorm(s));
    
    // Sigma := sqrt(Sigma^2), where each sigma > relTol*twoNorm
    for( Int i=0; i<m; ++i )
    {
        const Real lambda = s.Get(i,0);
        if( lambda <= Real(0) || Sqrt(lambda) <= relTol*twoNorm )
        {
            s.Resize( i, 1 );
            break;
        }
        else
            s.Set( i, 0, Sqrt(lambda) );
    }
}

template<typename F>
inline void WideProduct
( const Matrix<F>& A,
        Matrix<Base<F>>& s,
  Base<F> tol,
  bool relative )
{
    DEBUG_ONLY(CSE cse("svd::WideProduct"))
    if( relative )
        WideRelativeProduct( A, s, tol );
    else
        WideAbsoluteProduct( A, s, tol );
}

template<typename F>
inline void
WideAbsoluteProduct
( const DistMatrix<F>& A,
        ElementalMatrix<Base<F>>& s, 
  Base<F> tol )
{
    DEBUG_ONLY(
      CSE cse("svd::WideAbsoluteProduct");
      if( A.Width() < A.Height() )
          LogicError("A must be at least as wide as it is tall");
      if( tol < 0 )
          LogicError("negative threshold does not make sense");
    )
    const Int n = A.Width();

    typedef Base<F> Real;
    const Real frobNorm = FrobeniusNorm( A );
    if( tol == Real(0) )
    {
        const Real eps = limits::Epsilon<Real>();
        tol = n*frobNorm*Sqrt(eps);
    }
    if( tol >= frobNorm )
    {
        s.Resize( 0, 1 );
        return;
    }

    // C := A A^H
    const Grid& g = A.Grid();
    DistMatrix<F> C( g );
    Herk( LOWER, NORMAL, Real(1), A, C );

    // [Sigma^2] := eig(C), where each sigma > tol
    HermitianEigSubset<Real> subset;
    subset.rangeSubset = true;
    subset.lowerBound = tol*tol;
    // NOTE: While the square of the Frobenius norm should be a strict upper
    //       bound with exact computation, it has been observed that it can
    //       be lower than the finite-precision result in practice
    subset.upperBound = 2*frobNorm*frobNorm;
    HermitianEig( LOWER, C, s, DESCENDING, subset );
    
    // Sigma := sqrt(Sigma^2)
    const Int localHeight = s.LocalHeight();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        s.SetLocal( iLoc, 0, Sqrt(s.GetLocal(iLoc,0)) );
}

template<typename F>
inline void
WideAbsoluteProduct
( const ElementalMatrix<F>& APre,
        ElementalMatrix<Base<F>>& s, 
  Base<F> tol )
{
    DEBUG_ONLY(CSE cse("svd::WideAbsoluteProduct"))
    DistMatrixReadProxy<F,F,MC,MR> AProx( APre );
    auto& A = AProx.GetLocked();
    WideAbsoluteProduct( A, s, tol );
}

template<typename F>
inline void
WideRelativeProduct
( const DistMatrix<F>& A,
        ElementalMatrix<Base<F>>& s, 
  Base<F> relTol )
{
    DEBUG_ONLY(
      CSE cse("svd::WideRelativeProduct");
      AssertSameGrids( A, s );
      if( A.Width() < A.Height() )
          LogicError("A must be at least as wide as it is tall");
      if( relTol < 0 )
          LogicError("negative threshold does not make sense");
    )
    const Int m = A.Height();

    // C := A A^H
    typedef Base<F> Real;
    const Grid& g = A.Grid();
    DistMatrix<F> C( g );
    Herk( LOWER, NORMAL, Real(1), A, C );

    // [Sigma^2] := eig(C)
    HermitianEig( LOWER, C, s, DESCENDING );
    const Real twoNorm = Sqrt(MaxNorm(s));
    
    // Sigma := sqrt(Sigma^2), where all sigmas > relTol*twoNorm
    DistMatrix<Real,STAR,STAR> s_STAR_STAR( s );
    for( Int i=0; i<m; ++i )
    {
        const Real lambda = s_STAR_STAR.GetLocal(i,0);
        if( lambda <= Real(0) || Sqrt(lambda) <= relTol*twoNorm )
        {
            s_STAR_STAR.Resize( i, 1 );
            break;
        }
        else
            s_STAR_STAR.SetLocal( i, 0, Sqrt(lambda) );
    }
    Copy( s_STAR_STAR, s );
}

template<typename F>
inline void
WideRelativeProduct
( const ElementalMatrix<F>& APre,
        ElementalMatrix<Base<F>>& s, 
  Base<F> relTol )
{
    DEBUG_ONLY(CSE cse("svd::WideRelativeProduct"))
    DistMatrixReadProxy<F,F,MC,MR> AProx( APre );
    auto& A = AProx.GetLocked();
    WideRelativeProduct( A, s, relTol );
}

template<typename F>
inline void WideProduct
( const ElementalMatrix<F>& A,
        ElementalMatrix<Base<F>>& s, 
  Base<F> tol,
  bool relative )
{
    DEBUG_ONLY(CSE cse("svd::WideProduct"))
    if( relative )
        WideRelativeProduct( A, s, tol );
    else
        WideAbsoluteProduct( A, s, tol );
}

template<typename F>
void Product
( const Matrix<F>& A,
        Matrix<Base<F>>& s,
  Base<F> tol,
  bool relative )
{
    DEBUG_ONLY(CSE cse("svd::Product"))
    if( A.Height() >= A.Width() )
        TallProduct( A, s, tol, relative );
    else
        WideProduct( A, s, tol, relative );
}

template<typename F>
void Product
( const ElementalMatrix<F>& A,
        ElementalMatrix<Base<F>>& s, 
  Base<F> tol,
  bool relative )
{
    DEBUG_ONLY(CSE cse("svd::Product"))
    if( A.Height() >= A.Width() )
        TallProduct( A, s, tol, relative );
    else
        WideProduct( A, s, tol, relative );
}

} // namespace svd
} // namespace El

#endif // ifndef EL_SVD_PRODUCT_HPP
