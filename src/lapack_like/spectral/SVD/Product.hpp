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
SVDInfo TallAbsoluteProduct
( const Matrix<F>& A,
        Matrix<F>& U, 
        Matrix<Base<F>>& s,
        Matrix<F>& V,
  Base<F> tol,
  bool avoidU )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( A.Height() < A.Width() )
          LogicError("A must be at least as tall as it is wide");
      if( tol < 0 )
          LogicError("negative threshold does not make sense");
    )
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    const Real frobNorm = FrobeniusNorm( A );
    SVDInfo info;

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
        return info;
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
    HermitianEigCtrl<F> ctrl;
    ctrl.tridiagEigCtrl.subset = subset;
    ctrl.tridiagEigCtrl.sort = DESCENDING;
    HermitianEig( LOWER, C, s, V, ctrl );
    // TODO(poulson): Incorporate HermitianEigInfo into SVDInfo
    
    // Sigma := sqrt(Sigma^2)
    const Int k = s.Height();
    for( Int i=0; i<k; ++i )
        s(i) = Sqrt(s(i));

    if( !avoidU )
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

    return info;
}

template<typename F>
SVDInfo TallRelativeProduct
( const Matrix<F>& A,
        Matrix<F>& U,
        Matrix<Base<F>>& s,
        Matrix<F>& V,
  Base<F> relTol,
  bool avoidU )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( A.Height() < A.Width() )
          LogicError("A must be at least as tall as it is wide");
      if( relTol < 0 )
          LogicError("negative threshold does not make sense");
    )
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    SVDInfo info;

    if( relTol == Real(0) )
    {
        const Real eps = limits::Epsilon<Real>();
        relTol = m*Sqrt(eps);
    }

    // C := A^H A
    Matrix<F> C;
    Herk( LOWER, ADJOINT, Real(1), A, C );

    // [V,Sigma^2] := eig(C)
    HermitianEigCtrl<F> ctrl;
    ctrl.tridiagEigCtrl.sort = DESCENDING;
    HermitianEig( LOWER, C, s, V, ctrl );
    // TODO(poulson): Incorporate HermitianEigInfo into SVDInfo
    const Real twoNorm = Sqrt(MaxNorm(s));
    
    // Sigma := sqrt(Sigma^2), where all sigmas > relTol*twoNorm
    for( Int i=0; i<n; ++i )
    {
        const Real sigma = Sqrt(s(i));
        if( sigma <= relTol*twoNorm )
        {
            s.Resize( i, 1 );
            V.Resize( n, i );
            break;
        }
        else
            s(i) = sigma;
    }

    if( !avoidU )
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

    return info;
}

template<typename F>
SVDInfo TallProduct
( const Matrix<F>& A,
        Matrix<F>& U,
        Matrix<Base<F>>& s,
        Matrix<F>& V, 
  Base<F> tol,
  bool relative,
  bool avoidU )
{
    DEBUG_CSE
    if( relative )
        return TallRelativeProduct( A, U, s, V, tol, avoidU );
    else
        return TallAbsoluteProduct( A, U, s, V, tol, avoidU );
}

template<typename F>
SVDInfo TallAbsoluteProduct
( const DistMatrix<F>& A,
        DistMatrix<F>& U,
        AbstractDistMatrix<Base<F>>& s, 
        DistMatrix<F>& V,
  Base<F> tol,
  bool avoidU )
{
    DEBUG_CSE
    DEBUG_ONLY(
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
    SVDInfo info;

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
        return info;
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
    HermitianEigCtrl<F> ctrl;
    ctrl.tridiagEigCtrl.subset = subset;
    ctrl.tridiagEigCtrl.sort = DESCENDING;
    HermitianEig( LOWER, C, s, V, ctrl );
    // TODO(poulson): Incorporate HermitianEigInfo into SVDInfo
    
    // Sigma := sqrt(Sigma^2)
    const Int localHeight = s.LocalHeight();
    auto& sLoc = s.Matrix();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        sLoc(iLoc) = Sqrt(sLoc(iLoc));

    if( !avoidU )
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

    return info;
}

template<typename F>
SVDInfo TallAbsoluteProduct
( const AbstractDistMatrix<F>& APre,
        AbstractDistMatrix<F>& UPre,
        AbstractDistMatrix<Base<F>>& s, 
        AbstractDistMatrix<F>& VPre,
  Base<F> tol,
  bool avoidU )
{
    DEBUG_CSE
    DistMatrixReadProxy<F,F,MC,MR> AProx( APre );
    DistMatrixWriteProxy<F,F,MC,MR> UProx( UPre );
    DistMatrixWriteProxy<F,F,MC,MR> VProx( VPre );
    auto& A = AProx.GetLocked();
    auto& U = UProx.Get();
    auto& V = VProx.Get();
    return TallAbsoluteProduct( A, U, s, V, tol, avoidU );
}

template<typename F>
SVDInfo TallRelativeProduct
( const DistMatrix<F>& A,
        DistMatrix<F>& U,
        AbstractDistMatrix<Base<F>>& s, 
        DistMatrix<F>& V,
  Base<F> relTol,
  bool avoidU )
{
    DEBUG_CSE
    DEBUG_ONLY(
      AssertSameGrids( A, U, s, V );
      if( A.Height() < A.Width() )
          LogicError("A must be at least as tall as it is wide");
      if( relTol < 0 )
          LogicError("negative threshold does not make sense");
    )
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    const Grid& g = A.Grid();
    SVDInfo info;

    if( relTol == Real(0) )
    {
        const Real eps = limits::Epsilon<Real>();
        relTol = m*Sqrt(eps);
    }

    // C := A^H A
    DistMatrix<F> C(g);
    Herk( LOWER, ADJOINT, Real(1), A, C );

    // [V,Sigma^2] := eig(C)
    HermitianEigCtrl<F> ctrl;
    ctrl.tridiagEigCtrl.sort = DESCENDING;
    HermitianEig( LOWER, C, s, V, ctrl );
    // TODO(poulson): Incorporate HermitianEigInfo into SVDInfo
    const Real twoNorm = Sqrt(MaxNorm(s));

    // Sigma := sqrt(Sigma^2), where all sigmas > relTol*twoNorm
    DistMatrix<Real,STAR,STAR> s_STAR_STAR( s );
    auto& sLoc = s_STAR_STAR.Matrix();
    for( Int i=0; i<n; ++i )
    {
        const Real lambda = sLoc(i);
        if( lambda <= Real(0) || Sqrt(lambda) <= relTol*twoNorm )
        {
            s_STAR_STAR.Resize( i, 1 );
            V.Resize( n, i );
            break;
        }
        else
            sLoc(i) = Sqrt(lambda);
    }
    Copy( s_STAR_STAR, s );

    if( !avoidU )
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

    return info;
}

template<typename F>
SVDInfo TallRelativeProduct
( const AbstractDistMatrix<F>& APre,
        AbstractDistMatrix<F>& UPre,
        AbstractDistMatrix<Base<F>>& s, 
        AbstractDistMatrix<F>& VPre,
  Base<F> relTol,
  bool avoidU )
{
    DEBUG_CSE
    DistMatrixReadProxy<F,F,MC,MR> AProx( APre );
    DistMatrixWriteProxy<F,F,MC,MR> UProx( UPre );
    DistMatrixWriteProxy<F,F,MC,MR> VProx( VPre );
    auto& A = AProx.GetLocked();
    auto& U = UProx.Get();
    auto& V = VProx.Get();
    return TallRelativeProduct( A, U, s, V, relTol, avoidU );
}

template<typename F>
SVDInfo TallProduct
( const AbstractDistMatrix<F>& A,
        AbstractDistMatrix<F>& U,
        AbstractDistMatrix<Base<F>>& s, 
        AbstractDistMatrix<F>& V,
  Base<F> tol,
  bool relative,
  bool avoidU )
{
    DEBUG_CSE
    if( relative )
        return TallRelativeProduct( A, U, s, V, tol, avoidU );
    else
        return TallAbsoluteProduct( A, U, s, V, tol, avoidU );
}

template<typename F>
SVDInfo TallAbsoluteProduct
( const DistMatrix<F,VC,STAR>& A,
        DistMatrix<F,VC,STAR>& U,
        DistMatrix<Base<F>,STAR,STAR>& s, 
        DistMatrix<F,STAR,STAR>& V,
  Base<F> tol,
  bool avoidU )
{
    DEBUG_CSE
    DEBUG_ONLY(
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
    SVDInfo info;

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
        return info;
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
    HermitianEigCtrl<F> ctrl;
    ctrl.tridiagEigCtrl.subset = subset;
    ctrl.tridiagEigCtrl.sort = DESCENDING;
    HermitianEig( LOWER, C, s, V, ctrl );
    // TODO(poulson): Incorporate HermitianEigInfo into SVDInfo
    const int k = s.Height();
    
    // Sigma := sqrt(Sigma^2)
    auto& sLoc = s.Matrix();
    for( Int i=0; i<k; ++i )
        sLoc(i) = Sqrt(sLoc(i));

    if( !avoidU )
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

    return info;
}

template<typename F>
SVDInfo TallAbsoluteProduct
( const DistMatrix<F,VC,STAR>& A,
        AbstractDistMatrix<F>& UPre,
        AbstractDistMatrix<Base<F>>& sPre, 
        AbstractDistMatrix<F>& VPre,
  Base<F> tol,
  bool avoidU )
{
    DEBUG_CSE
    DistMatrixWriteProxy<F,F,VC,STAR> UProx( UPre );
    DistMatrixWriteProxy<Base<F>,Base<F>,STAR,STAR> sProx( sPre );
    DistMatrixWriteProxy<F,F,STAR,STAR> VProx( VPre );
    auto& s = sProx.Get();
    auto& U = UProx.Get();
    auto& V = VProx.Get();
    return TallAbsoluteProduct( A, U, s, V, tol, avoidU );
}

template<typename F>
SVDInfo TallRelativeProduct
( const DistMatrix<F,VC,STAR>& A,
        DistMatrix<F,VC,STAR>& U,
        DistMatrix<Base<F>,STAR,STAR>& s, 
        DistMatrix<F,STAR,STAR>& V,
  Base<F> relTol,
  bool avoidU )
{
    DEBUG_CSE
    DEBUG_ONLY(
      AssertSameGrids( A, U, s, V );
      if( A.Height() < A.Width() )
          LogicError("A must be at least as tall as it is wide");
      if( relTol < 0 )
          LogicError("negative threshold does not make sense");
    )
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    const Grid& g = A.Grid();
    SVDInfo info;

    if( relTol == Real(0) )
    {
        const Real eps = limits::Epsilon<Real>();
        relTol = m*Sqrt(eps);
    }

    // C := A^H A
    DistMatrix<F,STAR,STAR> C(g);
    Zeros( C, n, n );
    Herk( LOWER, ADJOINT, Real(1), A.LockedMatrix(), Real(0), C.Matrix() );
    El::AllReduce( C, A.ColComm() );

    // [V,Sigma^2] := eig(C)
    HermitianEigCtrl<F> ctrl;
    ctrl.tridiagEigCtrl.sort = DESCENDING;
    HermitianEig( LOWER, C, s, V, ctrl );
    // TODO(poulson): Incorporate HermitianEigInfo into SVDInfo
    const Real twoNorm = Sqrt(MaxNorm(s));
    
    // Sigma := sqrt(Sigma^2), where each sigma > twoNorm*relTol
    auto& sLoc = s.Matrix();
    for( Int i=0; i<n; ++i )
    {
        const Real lambda = sLoc(i);
        if( lambda <= Real(0) || Sqrt(lambda) <= relTol*twoNorm )
        {
            s.Resize( i, 1 );
            V.Resize( n, i );
            break;
        }
        else
            sLoc(i) = Sqrt(lambda);
    }
    const int k = s.Height();

    if( !avoidU )
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

    return info;
}

template<typename F>
SVDInfo TallRelativeProduct
( const DistMatrix<F,VC,STAR>& A,
        AbstractDistMatrix<F>& UPre,
        AbstractDistMatrix<Base<F>>& sPre, 
        AbstractDistMatrix<F>& VPre,
  Base<F> relTol,
  bool avoidU )
{
    DEBUG_CSE
    DistMatrixWriteProxy<Base<F>,Base<F>,STAR,STAR> sProx( sPre );
    DistMatrixWriteProxy<F,F,VC,STAR> UProx( UPre );
    DistMatrixWriteProxy<F,F,STAR,STAR> VProx( VPre );
    auto& s = sProx.Get();
    auto& U = UProx.Get();
    auto& V = VProx.Get();
    return TallRelativeProduct( A, U, s, V, relTol, avoidU );
}

template<typename F>
SVDInfo TallProduct
( const DistMatrix<F,VC,STAR>& A,
        AbstractDistMatrix<F>& U,
        AbstractDistMatrix<Base<F>>& s, 
        AbstractDistMatrix<F>& V,
  Base<F> tol,
  bool relative,
  bool avoidU )
{
    DEBUG_CSE
    if( relative )
        return TallRelativeProduct( A, U, s, V, tol, avoidU );
    else
        return TallAbsoluteProduct( A, U, s, V, tol, avoidU );
}

template<typename F>
SVDInfo WideAbsoluteProduct
( const Matrix<F>& A,
        Matrix<F>& U,
        Matrix<Base<F>>& s,
        Matrix<F>& V,
  Base<F> tol,
  bool avoidV )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( A.Width() < A.Height() )
          LogicError("A must be at least as wide as it is tall");
      if( tol < 0 )
          LogicError("negative threshold does not make sense");
    )
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    const Real frobNorm = FrobeniusNorm( A );
    SVDInfo info;

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
        return info;
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
    HermitianEigCtrl<F> ctrl;
    ctrl.tridiagEigCtrl.subset = subset;
    ctrl.tridiagEigCtrl.sort = DESCENDING;
    HermitianEig( LOWER, C, s, U, ctrl );
    // TODO(poulson): Incorporate HermitianEigInfo into SVDInfo
    
    // Sigma := sqrt(Sigma^2)
    const Int k = s.Height();
    for( Int i=0; i<k; ++i )
        s(i) = Sqrt(s(i));

    if( !avoidV )
    {
        // (Sigma V) := A^H U
        Gemm( ADJOINT, NORMAL, F(1), A, U, V );

        // Normalize each column of Sigma V
        Matrix<Base<F>> colNorms;
        ColumnTwoNorms( V, colNorms );
        DiagonalSolve( RIGHT, NORMAL, colNorms, V );
    }

    return info;
}

template<typename F>
SVDInfo WideRelativeProduct
( const Matrix<F>& A,
        Matrix<F>& U,
        Matrix<Base<F>>& s,
        Matrix<F>& V,
  Base<F> relTol,
  bool avoidV )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( A.Width() < A.Height() )
          LogicError("A must be at least as wide as it is tall");
      if( relTol < 0 )
          LogicError("negative threshold does not make sense");
    )
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    SVDInfo info;

    if( relTol == Real(0) )
    {
        const Real eps = limits::Epsilon<Real>();
        relTol = n*Sqrt(eps);
    }

    // C := A A^H
    Matrix<F> C;
    Herk( LOWER, NORMAL, Real(1), A, C );

    // [U,Sigma^2] := eig(C)
    HermitianEigCtrl<F> ctrl;
    ctrl.tridiagEigCtrl.sort = DESCENDING;
    HermitianEig( LOWER, C, s, U, ctrl );
    // TODO(poulson): Incorporate HermitianEigInfo into SVDInfo
    const Real twoNorm = Sqrt(MaxNorm(s));
    
    // Sigma := sqrt(Sigma^2), where each sigma > relTol*twoNorm
    for( Int i=0; i<m; ++i )
    {
        const Real lambda = s(i);
        if( lambda <= Real(0) || Sqrt(lambda) <= relTol*twoNorm )
        {
            s.Resize( i, 1 );
            U.Resize( m, i );
            break;
        }
        else
            s(i) = Sqrt(lambda);
    }

    if( !avoidV )
    {
        // (Sigma V) := A^H U
        Gemm( ADJOINT, NORMAL, F(1), A, U, V );

        // Normalize each column of Sigma V
        Matrix<Base<F>> colNorms;
        ColumnTwoNorms( V, colNorms );
        DiagonalSolve( RIGHT, NORMAL, colNorms, V );
    }

    return info;
}

template<typename F>
SVDInfo WideProduct
( const Matrix<F>& A,
        Matrix<F>& U,
        Matrix<Base<F>>& s,
        Matrix<F>& V, 
  Base<F> tol,
  bool relative,
  bool avoidV )
{
    DEBUG_CSE
    if( relative )
        return WideRelativeProduct( A, U, s, V, tol, avoidV );
    else
        return WideAbsoluteProduct( A, U, s, V, tol, avoidV );
}

template<typename F>
SVDInfo WideAbsoluteProduct
( const DistMatrix<F>& A,
        DistMatrix<F>& U,
        AbstractDistMatrix<Base<F>>& s, 
        DistMatrix<F>& V,
  Base<F> tol,
  bool avoidV )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( A.Width() < A.Height() )
          LogicError("A must be at least as wide as it is tall");
      if( tol < 0 )
          LogicError("negative threshold does not make sense");
    )
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    SVDInfo info;

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
        return info;
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
    HermitianEigCtrl<F> ctrl;
    ctrl.tridiagEigCtrl.subset = subset;
    ctrl.tridiagEigCtrl.sort = DESCENDING;
    HermitianEig( LOWER, C, s, U, ctrl );
    
    // Sigma := sqrt(Sigma^2)
    const Int localHeight = s.LocalHeight();
    auto& sLoc = s.Matrix();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        sLoc(iLoc) = Sqrt(sLoc(iLoc));

    if( !avoidV )
    {
        // (Sigma V) := A^H U
        Gemm( ADJOINT, NORMAL, F(1), A, U, V );

        // Normalize each column of Sigma V
        DistMatrix<Real,MR,STAR> colNorms(g);
        ColumnTwoNorms( V, colNorms );
        DiagonalSolve( RIGHT, NORMAL, colNorms, V );
    }

    return info;
}

template<typename F>
SVDInfo WideAbsoluteProduct
( const AbstractDistMatrix<F>& APre,
        AbstractDistMatrix<F>& UPre,
        AbstractDistMatrix<Base<F>>& s, 
        AbstractDistMatrix<F>& VPre,
  Base<F> tol,
  bool avoidV )
{
    DEBUG_CSE
    DistMatrixReadProxy<F,F,MC,MR> AProx( APre );
    DistMatrixWriteProxy<F,F,MC,MR> UProx( UPre );
    DistMatrixWriteProxy<F,F,MC,MR> VProx( VPre );
    auto& A = AProx.GetLocked();
    auto& U = UProx.Get();
    auto& V = VProx.Get();
    return WideAbsoluteProduct( A, U, s, V, tol, avoidV );
}

template<typename F>
SVDInfo WideRelativeProduct
( const DistMatrix<F>& A,
        DistMatrix<F>& U,
        AbstractDistMatrix<Base<F>>& s, 
        DistMatrix<F>& V,
  Base<F> relTol,
  bool avoidV )
{
    DEBUG_CSE
    DEBUG_ONLY(
      AssertSameGrids( A, U, s, V );
      if( A.Width() < A.Height() )
          LogicError("A must be at least as wide as it is tall");
      if( relTol < 0 )
          LogicError("negative threshold does not make sense");
    )
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    const Grid& g = A.Grid();
    SVDInfo info;

    if( relTol == Real(0) )
    {
        const Real eps = limits::Epsilon<Real>();
        relTol = n*Sqrt(eps);
    }

    // C := A A^H
    DistMatrix<F> C( g );
    Herk( LOWER, NORMAL, Real(1), A, C );

    // [U,Sigma^2] := eig(C)
    HermitianEigCtrl<F> ctrl;
    ctrl.tridiagEigCtrl.sort = DESCENDING;
    HermitianEig( LOWER, C, s, U, ctrl );
    const Real twoNorm = Sqrt(MaxNorm(s));
    
    // Sigma := sqrt(Sigma^2), where all sigmas > relTol*twoNorm
    DistMatrix<Real,STAR,STAR> s_STAR_STAR( s );
    auto& sLoc = s_STAR_STAR.Matrix();
    for( Int i=0; i<m; ++i )
    {
        const Real lambda = sLoc(i);
        if( lambda <= Real(0) || Sqrt(lambda) <= relTol*twoNorm )
        {
            s_STAR_STAR.Resize( i, 1 );
            U.Resize( m, i );
            break;
        }
        else
            sLoc(i) = Sqrt(lambda);
    }
    Copy( s_STAR_STAR, s );

    if( !avoidV )
    {
        // (Sigma V) := A^H U
        Gemm( ADJOINT, NORMAL, F(1), A, U, V );

        // Normalize each column of Sigma V
        DistMatrix<Real,MR,STAR> colNorms(g);
        ColumnTwoNorms( V, colNorms );
        DiagonalSolve( RIGHT, NORMAL, colNorms, V );
    }

    return info;
}

template<typename F>
SVDInfo WideRelativeProduct
( const AbstractDistMatrix<F>& APre,
        AbstractDistMatrix<F>& UPre,
        AbstractDistMatrix<Base<F>>& s, 
        AbstractDistMatrix<F>& VPre,
  Base<F> relTol,
  bool avoidV )
{
    DEBUG_CSE
    DistMatrixReadProxy<F,F,MC,MR> AProx( APre );
    DistMatrixWriteProxy<F,F,MC,MR> UProx( UPre );
    DistMatrixWriteProxy<F,F,MC,MR> VProx( VPre );
    auto& A = AProx.GetLocked();
    auto& U = UProx.Get();
    auto& V = VProx.Get();
    return WideRelativeProduct( A, U, s, V, relTol, avoidV );
}

template<typename F>
SVDInfo WideProduct
( const AbstractDistMatrix<F>& A,
        AbstractDistMatrix<F>& U,
        AbstractDistMatrix<Base<F>>& s, 
        AbstractDistMatrix<F>& V,
  Base<F> tol,
  bool relative,
  bool avoidV )
{
    DEBUG_CSE
    if( relative )
        return WideRelativeProduct( A, U, s, V, tol, avoidV );
    else
        return WideAbsoluteProduct( A, U, s, V, tol, avoidV );
}

// NOTE: [* ,VR] WideProduct would produce U with different distribution
//       than A. It makes more sense to overwrite A with V'.
// TODO: Update the above note and the following routines now that A is not
//       overwritten

template<typename F>
SVDInfo Product
( const Matrix<F>& A,
        Matrix<F>& U, 
        Matrix<Base<F>>& s,
        Matrix<F>& V, 
  Base<F> tol,
  bool relative,
  bool avoidU,
  bool avoidV )
{
    DEBUG_CSE
    // TODO: If m(A) >=~ n(A) but avoidV requested, find a way to make use of it
    if( A.Height() >= A.Width() )
        return TallProduct( A, U, s, V, tol, relative, avoidU );
    else
        return WideProduct( A, U, s, V, tol, relative, avoidV );
}

template<typename F>
SVDInfo Product
( const AbstractDistMatrix<F>& A,
        AbstractDistMatrix<F>& U,
        AbstractDistMatrix<Base<F>>& s, 
        AbstractDistMatrix<F>& V,
  Base<F> tol,
  bool relative,
  bool avoidU,
  bool avoidV )
{
    DEBUG_CSE
    // TODO: If m(A) >=~ n(A) but avoidV requested, find a way to make use of it
    if( A.Height() >= A.Width() )
        return TallProduct( A, U, s, V, tol, relative, avoidU );
    else
        return WideProduct( A, U, s, V, tol, relative, avoidV );
}

// Compute singular values
// =======================

template<typename F>
SVDInfo TallAbsoluteProduct
( const Matrix<F>& A,
        Matrix<Base<F>>& s,
  Base<F> tol )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( A.Height() < A.Width() )
          LogicError("A must be at least as tall as it is wide");
      if( tol < 0 )
          LogicError("negative threshold does not make sense");
    )
    typedef Base<F> Real;
    const Int m = A.Height();
    const Real frobNorm = FrobeniusNorm( A );
    SVDInfo info;

    if( tol == Real(0) )
    {
        const Real eps = limits::Epsilon<Real>();
        tol = m*frobNorm*Sqrt(eps);
    }
    if( tol >= frobNorm )
    {
        s.Resize( 0, 1 );
        return info;
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
    HermitianEigCtrl<F> ctrl;
    ctrl.tridiagEigCtrl.subset = subset;
    ctrl.tridiagEigCtrl.sort = DESCENDING;
    HermitianEig( LOWER, C, s, ctrl );
    
    // Sigma := sqrt(Sigma^2)
    const Int k = s.Height();
    for( Int i=0; i<k; ++i )
        s(i) = Sqrt(s(i));

    return info;
}

template<typename F>
SVDInfo TallRelativeProduct
( const Matrix<F>& A,
        Matrix<Base<F>>& s,
  Base<F> relTol )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( A.Height() < A.Width() )
          LogicError("A must be at least as tall as it is wide");
      if( relTol < 0 )
          LogicError("negative threshold does not make sense");
    )
    typedef Base<F> Real;
    const Int n = A.Width();
    SVDInfo info;

    // C := A^H A
    Matrix<F> C;
    Herk( LOWER, ADJOINT, Real(1), A, C );

    // [Sigma^2] := eig(C)
    HermitianEigCtrl<F> ctrl;
    ctrl.tridiagEigCtrl.sort = DESCENDING;
    HermitianEig( LOWER, C, s, ctrl );
    const Real twoNorm = Sqrt(MaxNorm(s));
    
    // Sigma := sqrt(Sigma^2), where all sigmas > relTol*twoNorm
    for( Int i=0; i<n; ++i )
    {
        const Real sigma = Sqrt(s(i));
        if( sigma <= relTol*twoNorm )
        {
            s.Resize( i, 1 );
            break;
        }
        else
            s(i) = sigma;
    }

    return info;
}

template<typename F>
SVDInfo TallProduct
( const Matrix<F>& A,
        Matrix<Base<F>>& s,
  Base<F> tol,
  bool relative )
{
    DEBUG_CSE
    if( relative )
        return TallRelativeProduct( A, s, tol );
    else
        return TallAbsoluteProduct( A, s, tol );
}

template<typename F>
SVDInfo TallAbsoluteProduct
( const DistMatrix<F>& A,
        AbstractDistMatrix<Base<F>>& s, 
  Base<F> tol )
{
    DEBUG_CSE
    DEBUG_ONLY(
      AssertSameGrids( A, s );
      if( A.Height() < A.Width() )
          LogicError("A must be at least as tall as it is wide");
      if( tol < 0 )
          LogicError("negative threshold does not make sense");
    )
    typedef Base<F> Real;
    const Int m = A.Height();
    const Grid& g = A.Grid();
    const Real frobNorm = FrobeniusNorm( A );
    SVDInfo info;

    if( tol == Real(0) )
    {
        const Real eps = limits::Epsilon<Real>();
        tol = m*frobNorm*Sqrt(eps);
    }
    if( tol >= frobNorm )
    {
        s.Resize( 0, 1 );
        return info;
    }

    // C := A^H A
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
    HermitianEigCtrl<F> ctrl;
    ctrl.tridiagEigCtrl.sort = DESCENDING;
    ctrl.tridiagEigCtrl.subset = subset;
    HermitianEig( LOWER, C, s, ctrl );
    
    // Sigma := sqrt(Sigma^2)
    const Int localHeight = s.LocalHeight();
    auto& sLoc = s.Matrix();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        sLoc(iLoc) = Sqrt(sLoc(iLoc));

    return info;
}

template<typename F>
SVDInfo TallAbsoluteProduct
( const AbstractDistMatrix<F>& APre,
        AbstractDistMatrix<Base<F>>& s, 
  Base<F> tol )
{
    DEBUG_CSE
    DistMatrixReadProxy<F,F,MC,MR> AProx( APre );
    auto& A = AProx.GetLocked();
    return TallAbsoluteProduct( A, s, tol );
}

template<typename F>
SVDInfo TallRelativeProduct
( const DistMatrix<F>& A,
        AbstractDistMatrix<Base<F>>& s, 
  Base<F> relTol )
{
    DEBUG_CSE
    DEBUG_ONLY(
      AssertSameGrids( A, s );
      if( A.Height() < A.Width() )
          LogicError("A must be at least as tall as it is wide");
      if( relTol < 0 )
          LogicError("negative threshold does not make sense");
    )
    typedef Base<F> Real;
    const Int n = A.Width();
    const Grid& g = A.Grid();
    SVDInfo info;

    // C := A^H A
    DistMatrix<F> C(g);
    Herk( LOWER, ADJOINT, Real(1), A, C );

    // [Sigma^2] := eig(C)
    HermitianEigCtrl<F> ctrl;
    ctrl.tridiagEigCtrl.sort = DESCENDING;
    HermitianEig( LOWER, C, s, ctrl );
    const Real twoNorm = Sqrt(MaxNorm(s));

    // Sigma := sqrt(Sigma^2), where all sigmas > relTol*twoNorm
    DistMatrix<Real,STAR,STAR> s_STAR_STAR( s );
    auto& sLoc = s_STAR_STAR.Matrix();
    for( Int i=0; i<n; ++i )
    {
        const Real lambda = sLoc(i);
        if( lambda <= Real(0) || Sqrt(lambda) <= relTol*twoNorm )
        {
            s_STAR_STAR.Resize( i, 1 );
            break;
        }
        else
            sLoc(i) = Sqrt(lambda);
    }
    Copy( s_STAR_STAR, s );

    return info;
}

template<typename F>
SVDInfo TallRelativeProduct
( const AbstractDistMatrix<F>& APre,
        AbstractDistMatrix<Base<F>>& s, 
  Base<F> relTol )
{
    DEBUG_CSE
    DistMatrixReadProxy<F,F,MC,MR> AProx( APre );
    auto& A = AProx.GetLocked();
    return TallRelativeProduct( A, s, relTol );
}

template<typename F>
SVDInfo TallProduct
( const AbstractDistMatrix<F>& A,
        AbstractDistMatrix<Base<F>>& s, 
  Base<F> tol,
  bool relative )
{
    DEBUG_CSE
    if( relative )
        return TallRelativeProduct( A, s, tol );
    else
        return TallAbsoluteProduct( A, s, tol );
}

template<typename F>
SVDInfo TallAbsoluteProduct
( const DistMatrix<F,VC,STAR>& A,
        DistMatrix<Base<F>,STAR,STAR>& s, 
  Base<F> tol )
{
    DEBUG_CSE
    DEBUG_ONLY(
      AssertSameGrids( A, s );
      if( A.Height() < A.Width() )
          LogicError("A must be at least as tall as it is wide");
      if( tol < 0 )
          LogicError("negative threshold does not make sense");
    )
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    const Real frobNorm = FrobeniusNorm( A );
    SVDInfo info;

    if( tol == Real(0) )
    {
        const Real eps = limits::Epsilon<Real>();
        tol = m*frobNorm*Sqrt(eps);
    }
    if( tol >= frobNorm )
    {
        s.Resize( 0, 1 );
        return info;
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
    HermitianEigCtrl<F> ctrl;
    ctrl.tridiagEigCtrl.subset = subset;
    ctrl.tridiagEigCtrl.sort = DESCENDING;
    HermitianEig( LOWER, C, s, ctrl );
    const int k = s.Height();
    
    // Sigma := sqrt(Sigma^2)
    auto& sLoc = s.Matrix();
    for( Int i=0; i<k; ++i )
        sLoc(i) = Sqrt(sLoc(i));

    return info;
}

template<typename F>
SVDInfo TallAbsoluteProduct
( const DistMatrix<F,VC,STAR>& A,
        AbstractDistMatrix<Base<F>>& sPre, 
  Base<F> tol )
{
    DEBUG_CSE
    typedef Base<F> Real;
    DistMatrixWriteProxy<Real,Real,STAR,STAR> sProx( sPre );
    auto& s = sProx.Get();
    return TallAbsoluteProduct( A, s, tol );
}

template<typename F>
SVDInfo TallRelativeProduct
( const DistMatrix<F,VC,STAR>& A,
        DistMatrix<Base<F>,STAR,STAR>& s, 
  Base<F> relTol )
{
    DEBUG_CSE
    DEBUG_ONLY(
      AssertSameGrids( A, s );
      if( A.Height() < A.Width() )
          LogicError("A must be at least as tall as it is wide");
      if( relTol < 0 )
          LogicError("negative threshold does not make sense");
    )
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    const Grid& g = A.Grid();
    SVDInfo info;

    // C := A^H A
    DistMatrix<F,STAR,STAR> C(g);
    Zeros( C, n, n );
    Herk( LOWER, ADJOINT, Real(1), A.LockedMatrix(), Real(0), C.Matrix() );
    El::AllReduce( C, A.ColComm() );

    // [V,Sigma^2] := eig(C)
    HermitianEigCtrl<F> ctrl;
    ctrl.tridiagEigCtrl.sort = DESCENDING;
    HermitianEig( LOWER, C, s, ctrl );
    const Real twoNorm = Sqrt(MaxNorm(s));
    
    // Sigma := sqrt(Sigma^2), where each sigma > twoNorm*relTol
    auto& sLoc = s.Matrix();
    for( Int i=0; i<n; ++i )
    {
        const Real lambda = sLoc(i);
        if( lambda <= Real(0) || Sqrt(lambda) <= relTol*twoNorm )
        {
            s.Resize( i, 1 );
            break;
        }
        else
            sLoc(i) = Sqrt(lambda);
    }
    const int k = s.Height();

    return info;
}

template<typename F>
SVDInfo TallRelativeProduct
( const DistMatrix<F,VC,STAR>& A,
        AbstractDistMatrix<Base<F>>& sPre, 
  Base<F> relTol )
{
    DEBUG_CSE
    typedef Base<F> Real;
    DistMatrixWriteProxy<Real,Real,STAR,STAR> sProx( sPre );
    auto& s = sProx.Get();
    return TallRelativeProduct( A, s, relTol );
}

template<typename F>
SVDInfo TallProduct
( const DistMatrix<F,VC,STAR>& A,
        AbstractDistMatrix<Base<F>>& s, 
  Base<F> tol,
  bool relative )
{
    DEBUG_CSE
    if( relative )
        return TallRelativeProduct( A, s, tol );
    else
        return TallAbsoluteProduct( A, s, tol );
}

template<typename F>
SVDInfo WideAbsoluteProduct
( const Matrix<F>& A,
        Matrix<Base<F>>& s,
  Base<F> tol )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( A.Width() < A.Height() )
          LogicError("A must be at least as wide as it is tall");
      if( tol < 0 )
          LogicError("negative threshold does not make sense");
    )
    typedef Base<F> Real;
    const Int n = A.Width();
    const Real frobNorm = FrobeniusNorm( A );
    SVDInfo info;

    if( tol == Real(0) )
    {
        const Real eps = limits::Epsilon<Real>();
        tol = n*frobNorm*Sqrt(eps);
    }
    if( tol >= frobNorm )
    {
        s.Resize( 0, 1 );
        return info;
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
    HermitianEigCtrl<F> ctrl;
    ctrl.tridiagEigCtrl.subset = subset;
    ctrl.tridiagEigCtrl.sort = DESCENDING;
    HermitianEig( LOWER, C, s, ctrl );
    
    // Sigma := sqrt(Sigma^2)
    const Int k = s.Height();
    for( Int i=0; i<k; ++i )
        s(i) = Sqrt(s(i));

    return info;
}

template<typename F>
SVDInfo WideRelativeProduct
( const Matrix<F>& A,
        Matrix<Base<F>>& s,
  Base<F> relTol )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( A.Width() < A.Height() )
          LogicError("A must be at least as wide as it is tall");
      if( relTol < 0 )
          LogicError("negative threshold does not make sense");
    )
    typedef Base<F> Real;
    const Int m = A.Height();
    SVDInfo info;

    // C := A A^H
    Matrix<F> C;
    Herk( LOWER, NORMAL, Real(1), A, C );

    // [Sigma^2] := eig(C)
    HermitianEigCtrl<F> ctrl;
    ctrl.tridiagEigCtrl.sort = DESCENDING;
    HermitianEig( LOWER, C, s, ctrl );
    const Real twoNorm = Sqrt(MaxNorm(s));
    
    // Sigma := sqrt(Sigma^2), where each sigma > relTol*twoNorm
    for( Int i=0; i<m; ++i )
    {
        const Real lambda = s(i);
        if( lambda <= Real(0) || Sqrt(lambda) <= relTol*twoNorm )
        {
            s.Resize( i, 1 );
            break;
        }
        else
            s(i) = Sqrt(lambda);
    }

    return info;
}

template<typename F>
SVDInfo WideProduct
( const Matrix<F>& A,
        Matrix<Base<F>>& s,
  Base<F> tol,
  bool relative )
{
    DEBUG_CSE
    if( relative )
        return WideRelativeProduct( A, s, tol );
    else
        return WideAbsoluteProduct( A, s, tol );
}

template<typename F>
SVDInfo WideAbsoluteProduct
( const DistMatrix<F>& A,
        AbstractDistMatrix<Base<F>>& s, 
  Base<F> tol )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( A.Width() < A.Height() )
          LogicError("A must be at least as wide as it is tall");
      if( tol < 0 )
          LogicError("negative threshold does not make sense");
    )
    typedef Base<F> Real;
    const Int n = A.Width();
    const Real frobNorm = FrobeniusNorm( A );
    SVDInfo info;

    if( tol == Real(0) )
    {
        const Real eps = limits::Epsilon<Real>();
        tol = n*frobNorm*Sqrt(eps);
    }
    if( tol >= frobNorm )
    {
        s.Resize( 0, 1 );
        return info;
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
    HermitianEigCtrl<F> ctrl;
    ctrl.tridiagEigCtrl.subset = subset;
    ctrl.tridiagEigCtrl.sort = DESCENDING;
    HermitianEig( LOWER, C, s, ctrl );
    
    // Sigma := sqrt(Sigma^2)
    const Int localHeight = s.LocalHeight();
    auto& sLoc = s.Matrix();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        sLoc(iLoc) = Sqrt(sLoc(iLoc));

    return info;
}

template<typename F>
SVDInfo WideAbsoluteProduct
( const AbstractDistMatrix<F>& APre,
        AbstractDistMatrix<Base<F>>& s, 
  Base<F> tol )
{
    DEBUG_CSE
    DistMatrixReadProxy<F,F,MC,MR> AProx( APre );
    auto& A = AProx.GetLocked();
    return WideAbsoluteProduct( A, s, tol );
}

template<typename F>
SVDInfo WideRelativeProduct
( const DistMatrix<F>& A,
        AbstractDistMatrix<Base<F>>& s, 
  Base<F> relTol )
{
    DEBUG_CSE
    DEBUG_ONLY(
      AssertSameGrids( A, s );
      if( A.Width() < A.Height() )
          LogicError("A must be at least as wide as it is tall");
      if( relTol < 0 )
          LogicError("negative threshold does not make sense");
    )
    const Int m = A.Height();
    typedef Base<F> Real;
    const Grid& g = A.Grid();
    SVDInfo info;

    // C := A A^H
    DistMatrix<F> C( g );
    Herk( LOWER, NORMAL, Real(1), A, C );

    // [Sigma^2] := eig(C)
    HermitianEigCtrl<F> ctrl;
    ctrl.tridiagEigCtrl.sort = DESCENDING;
    HermitianEig( LOWER, C, s, ctrl );
    const Real twoNorm = Sqrt(MaxNorm(s));
    
    // Sigma := sqrt(Sigma^2), where all sigmas > relTol*twoNorm
    DistMatrix<Real,STAR,STAR> s_STAR_STAR( s );
    auto& sLoc = s_STAR_STAR.Matrix();
    for( Int i=0; i<m; ++i )
    {
        const Real lambda = sLoc(i);
        if( lambda <= Real(0) || Sqrt(lambda) <= relTol*twoNorm )
        {
            s_STAR_STAR.Resize( i, 1 );
            break;
        }
        else
            sLoc(i) = Sqrt(lambda);
    }
    Copy( s_STAR_STAR, s );

    return info;
}

template<typename F>
SVDInfo WideRelativeProduct
( const AbstractDistMatrix<F>& APre,
        AbstractDistMatrix<Base<F>>& s, 
  Base<F> relTol )
{
    DEBUG_CSE
    DistMatrixReadProxy<F,F,MC,MR> AProx( APre );
    auto& A = AProx.GetLocked();
    return WideRelativeProduct( A, s, relTol );
}

template<typename F>
SVDInfo WideProduct
( const AbstractDistMatrix<F>& A,
        AbstractDistMatrix<Base<F>>& s, 
  Base<F> tol,
  bool relative )
{
    DEBUG_CSE
    if( relative )
        return WideRelativeProduct( A, s, tol );
    else
        return WideAbsoluteProduct( A, s, tol );
}

template<typename F>
SVDInfo Product
( const Matrix<F>& A,
        Matrix<Base<F>>& s,
  Base<F> tol,
  bool relative )
{
    DEBUG_CSE
    if( A.Height() >= A.Width() )
        return TallProduct( A, s, tol, relative );
    else
        return WideProduct( A, s, tol, relative );
}

template<typename F>
SVDInfo Product
( const AbstractDistMatrix<F>& A,
        AbstractDistMatrix<Base<F>>& s, 
  Base<F> tol,
  bool relative )
{
    DEBUG_CSE
    if( A.Height() >= A.Width() )
        return TallProduct( A, s, tol, relative );
    else
        return WideProduct( A, s, tol, relative );
}

} // namespace svd
} // namespace El

#endif // ifndef EL_SVD_PRODUCT_HPP
