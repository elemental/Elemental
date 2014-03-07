/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_PSEUDOSPECTRUM_HPP
#define ELEM_PSEUDOSPECTRUM_HPP

#include ELEM_FROBENIUSNORM_INC
#include ELEM_MAXNORM_INC
#include ELEM_ONENORM_INC
#include ELEM_TWONORMESTIMATE_INC
#include ELEM_HESSENBERG_INC
#include ELEM_SCHUR_INC

#include "./Pseudospectrum/Power.hpp"
#include "./Pseudospectrum/Lanczos.hpp"
#include "./Pseudospectrum/KrylovSpectral.hpp"
#include "./Pseudospectrum/Analytic.hpp"

namespace elem {

namespace pspec {

template<typename F>
inline bool
NumericallyNormal( const Matrix<F>& U, Base<F> tol )
{
    auto w = U.GetDiagonal();
    const Base<F> diagFrob = FrobeniusNorm( w );
    const Base<F> upperFrob = FrobeniusNorm( U );
    const Base<F> offDiagFrob = Sqrt(upperFrob*upperFrob-diagFrob*diagFrob);
    return offDiagFrob <= tol*diagFrob;
}

template<typename F>
inline bool
NumericallyNormal( const DistMatrix<F>& U, Base<F> tol )
{
    auto w = U.GetDiagonal();
    const Base<F> diagFrob = FrobeniusNorm( w );
    const Base<F> upperFrob = FrobeniusNorm( U );
    const Base<F> offDiagFrob = Sqrt(upperFrob*upperFrob-diagFrob*diagFrob);
    return offDiagFrob <= tol*diagFrob;
}

template<typename T>
inline void
ReshapeIntoGrids
( Int xSize, Int ySize,
  const Matrix<T>& invNorms,   const Matrix<Int>& itCounts,
        Matrix<T>& invNormMap,       Matrix<Int>& itCountMap )
{
#if 0    
    invNormMap.Resize( xSize, ySize );
    itCountMap.Resize( xSize, ySize );
    for( Int j=0; j<xSize; ++j )
    {
        auto normGridSub = View( invNormMap, 0, j, ySize, 1 );
        auto countGridSub = View( itCountMap, 0, j, ySize, 1 );
        auto shiftSub = LockedView( invNorms, j*ySize, 0, ySize, 1 );
        auto countSub = LockedView( itCounts, j*ySize, 0, ySize, 1 );
        normGridSub = shiftSub;
        countGridSub = countSub;
    }
#else
    // The sequential case can be optimized much more heavily than in parallel
    invNormMap.Resize( xSize, ySize, xSize );
    itCountMap.Resize( xSize, ySize, xSize );
    MemCopy( invNormMap.Buffer(), invNorms.LockedBuffer(), xSize*ySize );
    MemCopy( itCountMap.Buffer(), itCounts.LockedBuffer(), xSize*ySize );
#endif
}

template<typename T>
inline void
ReshapeIntoGrids
( Int xSize, Int ySize,
  const DistMatrix<T,VR,STAR>& invNorms, 
  const DistMatrix<Int,VR,STAR>& itCounts,
        DistMatrix<T>& invNormMap, 
        DistMatrix<Int>& itCountMap )
{
    invNormMap.SetGrid( invNorms.Grid() );
    itCountMap.SetGrid( invNorms.Grid() );
    invNormMap.Resize( xSize, ySize );
    itCountMap.Resize( xSize, ySize );
    for( Int j=0; j<xSize; ++j )
    {
        auto normGridSub = View( invNormMap, 0, j, ySize, 1 );
        auto countGridSub = View( itCountMap, 0, j, ySize, 1 );
        auto shiftSub = LockedView( invNorms, j*ySize, 0, ySize, 1 );
        auto countSub = LockedView( itCounts, j*ySize, 0, ySize, 1 );
        normGridSub = shiftSub;
        countGridSub = countSub;
    }
}

} // namespace pspec

template<typename F>
inline Matrix<Int>
TriangularPseudospectrum
( const Matrix<F>& U, const Matrix<Complex<BASE(F)> >& shifts, 
  Matrix<BASE(F)>& invNorms, bool lanczos=true, Int krylovSize=10, 
  bool reorthog=true, bool deflate=true, Int maxIts=1000, BASE(F) tol=1e-6, 
  bool progress=false )
{
    DEBUG_ONLY(CallStackEntry cse("TriangularPseudospectrum"))
    typedef Base<F> Real;
    typedef Complex<Real> C;

    Matrix<C> UCpx;
    if( IsComplex<F>::val )
        UCpx = LockedView( U );
    else
    {
        const Int n = U.Height();
        UCpx.Resize( n, n );
        for( Int j=0; j<n; ++j )
            for( Int i=0; i<n; ++i )
                UCpx.Set( i, j, U.Get(i,j) );
    }

    // Check if the off-diagonal is sufficiently small; if so, compute the 
    // pseudospectrum analytically from the eigenvalues. This also takes care
    // of the case where the matrix is a constant multiple of the identity 
    // matrix, which, after shifting, can lead to the zero matrix, which would 
    // cause problems for the Lanczos convergence criteria.
    Matrix<Int> itCounts;
    if( pspec::NumericallyNormal( UCpx, tol ) )
    {
        if( progress )
            std::cout << "Matrix was numerically normal" << std::endl;
        auto w = UCpx.GetDiagonal();
        pspec::Analytic( w, shifts, invNorms );
        Zeros( itCounts, shifts.Height(), 1 );        
        return itCounts;
    }

    if( lanczos )
    {
        if( krylovSize > 1 )
            itCounts =
               pspec::TriangularKrylovSpectral
               ( UCpx, shifts, invNorms, krylovSize, reorthog, deflate, maxIts, 
                 tol, progress );
        else
            itCounts = 
               pspec::TriangularLanczos
               ( UCpx, shifts, invNorms, deflate, maxIts, tol, progress );
    }
    else
        itCounts =
           pspec::TriangularPower
           ( UCpx, shifts, invNorms, deflate, maxIts, tol, progress );

    return itCounts;
}

template<typename F>
inline Matrix<Int>
HessenbergPseudospectrum
( const Matrix<F>& H, const Matrix<Complex<BASE(F)> >& shifts, 
  Matrix<BASE(F)>& invNorms, bool lanczos=true, Int krylovSize=10, 
  bool reorthog=true, bool deflate=true, Int maxIts=1000, BASE(F) tol=1e-6, 
  bool progress=false )
{
    DEBUG_ONLY(CallStackEntry cse("HessenbergPseudospectrum"))
    typedef Base<F> Real;
    typedef Complex<Real> C;

    Matrix<C> HCpx;
    if( IsComplex<F>::val )
        HCpx = LockedView( H );
    else
    {
        const Int n = H.Height();
        HCpx.Resize( n, n );
        for( Int j=0; j<n; ++j )
            for( Int i=0; i<n; ++i )
                HCpx.Set( i, j, H.Get(i,j) );
    }

    // TODO: Check if the subdiagonal is numerically zero, and, if so, revert to
    //       TriangularPseudospectrum?
    Matrix<Int> itCounts;
    if( lanczos )
    {
        if( krylovSize > 1 )
            itCounts =
               pspec::HessenbergKrylovSpectral
               ( HCpx, shifts, invNorms, krylovSize, reorthog, deflate, maxIts, 
                 tol, progress );
        else
            itCounts = 
               pspec::HessenbergLanczos
               ( HCpx, shifts, invNorms, deflate, maxIts, tol, progress );
    }
    else
        itCounts =
           pspec::HessenbergPower
           ( HCpx, shifts, invNorms, deflate, maxIts, tol, progress );

    return itCounts;
}

template<typename F>
inline DistMatrix<Int,VR,STAR>
TriangularPseudospectrum
( const DistMatrix<F>& U, const DistMatrix<Complex<BASE(F)>,VR,STAR>& shifts,
  DistMatrix<BASE(F),VR,STAR>& invNorms, bool lanczos=true, Int krylovSize=10, 
  bool reorthog=true, bool deflate=true, Int maxIts=1000, BASE(F) tol=1e-6, 
  bool progress=false )
{
    DEBUG_ONLY(CallStackEntry cse("TriangularPseudospectrum"))
    typedef Base<F> Real;
    typedef Complex<Real> C;

    const Grid& g = U.Grid();
    DistMatrix<C> UCpx(g);
    if( IsComplex<F>::val )
        UCpx = LockedView( U );
    else
    {
        UCpx.AlignWith( U );
        const Int n = U.Height();
        UCpx.Resize( n, n );
        const Int mLocal = U.LocalHeight();
        const Int nLocal = U.LocalWidth();
        for( Int jLoc=0; jLoc<nLocal; ++jLoc )
            for( Int iLoc=0; iLoc<mLocal; ++iLoc )
                UCpx.SetLocal( iLoc, jLoc, U.GetLocal(iLoc,jLoc) );
    }

    // Check if the off-diagonal is sufficiently small; if so, compute the 
    // pseudospectrum analytically from the eigenvalues. This also takes care
    // of the case where the matrix is a constant multiple of the identity 
    // matrix, which, after shifting, can lead to the zero matrix, which would 
    // cause problems for the Lanczos convergence criteria.
    DistMatrix<Int,VR,STAR> itCounts(g);
    if( pspec::NumericallyNormal( UCpx, tol ) )
    {
        if( progress && U.Grid().Rank() == 0 )
            std::cout << "Matrix was numerically normal" << std::endl;
        auto w = UCpx.GetDiagonal();
        DistMatrix<C,STAR,STAR> w_STAR_STAR( w );
        pspec::Analytic( w_STAR_STAR, shifts, invNorms );
        itCounts.AlignWith( shifts );
        Zeros( itCounts, shifts.Height(), 1 );
        return itCounts;
    }

    if( lanczos )
    {
        if( krylovSize > 1 )
            itCounts =
               pspec::TriangularKrylovSpectral
               ( UCpx, shifts, invNorms, krylovSize, reorthog, deflate, maxIts, 
                 tol, progress );
        else
            itCounts = 
               pspec::TriangularLanczos
               ( UCpx, shifts, invNorms, deflate, maxIts, tol, progress );
    }
    else
        itCounts =
           pspec::TriangularPower
           ( UCpx, shifts, invNorms, deflate, maxIts, tol, progress );

    return itCounts;
}

template<typename F>
inline DistMatrix<Int,VR,STAR>
HessenbergPseudospectrum
( const DistMatrix<F>& H, const DistMatrix<Complex<BASE(F)>,VR,STAR>& shifts,
  DistMatrix<BASE(F),VR,STAR>& invNorms, bool lanczos=true, Int krylovSize=10, 
  bool reorthog=true, bool deflate=true, Int maxIts=1000, BASE(F) tol=1e-6, 
  bool progress=false )
{
    DEBUG_ONLY(CallStackEntry cse("HessenbergPseudospectrum"))
    typedef Base<F> Real;
    typedef Complex<Real> C;

    const Grid& g = H.Grid();
    DistMatrix<C> HCpx(g);
    if( IsComplex<F>::val )
        HCpx = LockedView( H );
    else
    {
        HCpx.AlignWith( H );
        const Int n = H.Height();
        HCpx.Resize( n, n );
        const Int mLocal = H.LocalHeight();
        const Int nLocal = H.LocalWidth();
        for( Int jLoc=0; jLoc<nLocal; ++jLoc )
            for( Int iLoc=0; iLoc<mLocal; ++iLoc )
                HCpx.SetLocal( iLoc, jLoc, H.GetLocal(iLoc,jLoc) );
    }

    // TODO: Check if the subdiagonal is sufficiently small, and, if so, revert
    //       to TriangularPseudospectrum
    DistMatrix<Int,VC,STAR> itCounts(g);
    if( lanczos )
    {
        if( krylovSize > 1 )
            itCounts =
               pspec::HessenbergKrylovSpectral
               ( HCpx, shifts, invNorms, krylovSize, reorthog, deflate, maxIts, 
                 tol, progress );
        else
            itCounts = 
               pspec::HessenbergLanczos
               ( HCpx, shifts, invNorms, deflate, maxIts, tol, progress );
    }
    else
        itCounts =
           pspec::HessenbergPower
           ( HCpx, shifts, invNorms, deflate, maxIts, tol, progress );

    return itCounts;
}

namespace pspec {

template<typename F>
inline Matrix<Int>
Triangular
( const Matrix<F>& A, const Matrix<Complex<BASE(F)> >& shifts, 
  Matrix<BASE(F)>& invNorms, bool lanczos=true, Int krylovSize=10, 
  bool reorthog=true, bool deflate=true, Int maxIts=1000, BASE(F) tol=1e-6, 
  bool progress=false )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::Triangular"))
    typedef Base<F> Real;
    typedef Complex<Real> C;

    const Int n = A.Height();
    Matrix<C> U( n, n );
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<n; ++i )
            U.Set( i, j, A.Get(i,j) );

    Matrix<C> w;
    schur::QR( U, w );

    return TriangularPseudospectrum
           ( U, shifts, invNorms, lanczos, krylovSize, reorthog, deflate, 
             maxIts, tol, progress );
}

template<typename F>
inline Matrix<Int>
Hessenberg
( const Matrix<F>& A, const Matrix<Complex<BASE(F)> >& shifts, 
  Matrix<BASE(F)>& invNorms, bool lanczos=true, Int krylovSize=10, 
  bool reorthog=true, bool deflate=true, Int maxIts=1000, BASE(F) tol=1e-6, 
  bool progress=false )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::Hessenberg"))
    typedef Base<F> Real;
    typedef Complex<Real> C;

    const Int n = A.Height();
    Matrix<C> H( n, n );
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<n; ++i )
            H.Set( i, j, A.Get(i,j) );

    Hessenberg( UPPER, H );

    return HessenbergPseudospectrum
           ( H, shifts, invNorms, lanczos, krylovSize, reorthog, deflate, 
             maxIts, tol, progress );
}

template<typename F>
inline DistMatrix<Int,VR,STAR>
Triangular
( const DistMatrix<F>& A, const DistMatrix<Complex<BASE(F)>,VR,STAR>& shifts,
  DistMatrix<BASE(F),VR,STAR>& invNorms, bool lanczos=true, Int krylovSize=10, 
  bool reorthog=true, bool deflate=true, Int maxIts=1000, BASE(F) tol=1e-6, 
  bool progress=false )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::Triangular"))
    typedef Base<F> Real;
    typedef Complex<Real> C;

    const Grid& g = A.Grid();
    DistMatrix<C> U(g);
    U.AlignWith( A );
    const Int n = A.Height();
    U.Resize( n, n );
    const Int mLocal = A.LocalHeight();
    const Int nLocal = A.LocalWidth();
    for( Int jLoc=0; jLoc<nLocal; ++jLoc )
        for( Int iLoc=0; iLoc<mLocal; ++iLoc )
            U.SetLocal( iLoc, jLoc, A.GetLocal(iLoc,jLoc) );

    // We don't actually need the Schur vectors, but SDC requires their 
    // computation in order to form the full triangular factor
    DistMatrix<C> X(g);
    DistMatrix<C,VR,STAR> w(g);
    const bool formATR = true;
    // TODO: Expose these as options
    const Int cutoff = 256;
    const Int maxInnerIts = 2;
    const Int maxOuterIts = 10;
    const Base<F> signTol=tol/10;
    const Base<F> relTol=tol/10;
    const Base<F> spreadFactor=1e-6;
    const bool random=true;
    schur::SDC
    ( U, w, X, formATR, cutoff, maxInnerIts, maxOuterIts, signTol, relTol, 
      spreadFactor, random, progress );
    X.Empty();

    return TriangularPseudospectrum
           ( U, shifts, invNorms, lanczos, krylovSize, reorthog, deflate, 
             maxIts, tol, progress );
}

template<typename F>
inline DistMatrix<Int,VR,STAR>
Hessenberg
( const DistMatrix<F>& A, const DistMatrix<Complex<BASE(F)>,VR,STAR>& shifts,
  DistMatrix<BASE(F),VR,STAR>& invNorms, bool lanczos=true, Int krylovSize=10, 
  bool reorthog=true, bool deflate=true, Int maxIts=1000, BASE(F) tol=1e-6, 
  bool progress=false )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::Hessenberg"))
    typedef Base<F> Real;
    typedef Complex<Real> C;

    const Grid& g = A.Grid();
    DistMatrix<C> H(g);
    H.AlignWith( A );
    const Int n = A.Height();
    H.Resize( n, n );
    const Int mLocal = A.LocalHeight();
    const Int nLocal = A.LocalWidth();
    for( Int jLoc=0; jLoc<nLocal; ++jLoc )
        for( Int iLoc=0; iLoc<mLocal; ++iLoc )
            H.SetLocal( iLoc, jLoc, A.GetLocal(iLoc,jLoc) );

    Hessenberg( UPPER, H );

    return HessenbergPseudospectrum
           ( H, shifts, invNorms, lanczos, krylovSize, reorthog, deflate, 
             maxIts, tol, progress );
}

} // namespace pspec

template<typename F>
inline Matrix<Int>
Pseudospectrum
( const Matrix<F>& A, const Matrix<Complex<BASE(F)> >& shifts, 
  Matrix<BASE(F)>& invNorms, bool schur=true, bool lanczos=true, 
  Int krylovSize=10, bool reorthog=true, bool deflate=true, Int maxIts=1000, 
  BASE(F) tol=1e-6, bool progress=false )
{
    DEBUG_ONLY(CallStackEntry cse("Pseudospectrum"))
    if( schur )
        return pspec::Triangular
        ( A, shifts, invNorms, lanczos, krylovSize, reorthog, deflate, maxIts, 
          tol, progress );
    else
        return pspec::Hessenberg
        ( A, shifts, invNorms, lanczos, krylovSize, reorthog, deflate, maxIts, 
          tol, progress );
}

template<typename F>
inline DistMatrix<Int,VR,STAR>
Pseudospectrum
( const DistMatrix<F>& A, const DistMatrix<Complex<BASE(F)>,VR,STAR>& shifts,
  DistMatrix<BASE(F),VR,STAR>& invNorms, 
  bool schur=false, bool lanczos=true, Int krylovSize=10, 
  bool reorthog=true, bool deflate=true, Int maxIts=1000, BASE(F) tol=1e-6, 
  bool progress=false )
{
    DEBUG_ONLY(CallStackEntry cse("Pseudospectrum"))
    if( schur )
        return pspec::Triangular
        ( A, shifts, invNorms, lanczos, krylovSize, reorthog, deflate, maxIts, 
          tol, progress );
    else
        return pspec::Hessenberg
        ( A, shifts, invNorms, lanczos, krylovSize, reorthog, deflate, maxIts, 
          tol, progress );
}

template<typename F>
inline Matrix<Int>
TriangularPseudospectrum
( const Matrix<F>& U, Matrix<BASE(F)>& invNormMap, 
  Complex<BASE(F)> center, BASE(F) xWidth, BASE(F) yWidth,
  Int xSize, Int ySize, bool lanczos=true, Int krylovSize=10, 
  bool reorthog=true, bool deflate=true, Int maxIts=1000, BASE(F) tol=1e-6, 
  bool progress=false )
{
    DEBUG_ONLY(CallStackEntry cse("TriangularPseudospectrum"))
    typedef Base<F> Real;
    typedef Complex<Real> C;

    const Real xStep = xWidth/(xSize-1);
    const Real yStep = yWidth/(ySize-1);
    const C corner = center - C(xWidth/2,yWidth/2);
    Matrix<C> shifts( xSize*ySize, 1, U.Grid() );
    for( Int j=0; j<xSize*ySize; ++j )
    {
        const Int x = j / ySize;
        const Int y = j % ySize;
        shifts.Set( j, 0, corner+C(x*xStep,y*yStep) );
    }

    // Form the vector of invNorms
    Matrix<Real> invNorms;
    auto itCounts = 
        TriangularPseudospectrum
        ( U, shifts, invNorms, lanczos, krylovSize, reorthog, deflate, maxIts, 
          tol, progress );

    // Rearrange the vectors into grids
    Matrix<Int> itCountMap; 
    pspec::ReshapeIntoGrids
    ( xSize, ySize, invNorms, itCounts, invNormMap, itCountMap );
    return itCountMap;
}

template<typename F>
inline Matrix<Int>
HessenbergPseudospectrum
( const Matrix<F>& H, Matrix<BASE(F)>& invNormMap, 
  Complex<BASE(F)> center, BASE(F) xWidth, BASE(F) yWidth,
  Int xSize, Int ySize, bool lanczos=true, Int krylovSize=10, 
  bool reorthog=true, bool deflate=true, Int maxIts=1000, BASE(F) tol=1e-6, 
  bool progress=false )
{
    DEBUG_ONLY(CallStackEntry cse("HessenbergPseudospectrum"))
    typedef Base<F> Real;
    typedef Complex<Real> C;

    const Real xStep = xWidth/(xSize-1);
    const Real yStep = yWidth/(ySize-1);
    const C corner = center - C(xWidth/2,yWidth/2);
    Matrix<C> shifts( xSize*ySize, 1, H.Grid() );
    for( Int j=0; j<xSize*ySize; ++j )
    {
        const Int x = j / ySize;
        const Int y = j % ySize;
        shifts.Set( j, 0, corner+C(x*xStep,y*yStep) );
    }

    // Form the vector of invNorms
    Matrix<Real> invNorms;
    auto itCounts = 
        HessenbergPseudospectrum
        ( H, shifts, invNorms, lanczos, krylovSize, reorthog, deflate, maxIts, 
          tol, progress );

    // Rearrange the vectors into grids
    Matrix<Int> itCountMap; 
    pspec::ReshapeIntoGrids
    ( xSize, ySize, invNorms, itCounts, invNormMap, itCountMap );
    return itCountMap;
}

template<typename F>
inline DistMatrix<Int>
TriangularPseudospectrum
( const DistMatrix<F>& U, DistMatrix<BASE(F)>& invNormMap, 
  Complex<BASE(F)> center, BASE(F) xWidth, BASE(F) yWidth, Int xSize, Int ySize,
  bool lanczos=true, Int krylovSize=10, bool reorthog=true, bool deflate=true, 
  Int maxIts=1000, BASE(F) tol=1e-6, bool progress=false )
{
    DEBUG_ONLY(CallStackEntry cse("TriangularPseudospectrum"))
    typedef Base<F> Real;
    typedef Complex<Real> C;
    const Grid& g = U.Grid();

    const Real xStep = xWidth/(xSize-1);
    const Real yStep = yWidth/(ySize-1);
    const C corner = center - C(xWidth/2,yWidth/2);
    DistMatrix<C,VR,STAR> shifts( xSize*ySize, 1, g );
    const Int numLocShifts = shifts.LocalHeight();
    for( Int jLoc=0; jLoc<numLocShifts; ++jLoc )
    {
        const Int j = shifts.ColShift() + jLoc*shifts.ColStride();
        const Int x = j / ySize;
        const Int y = j % ySize;
        shifts.SetLocal( jLoc, 0, corner+C(x*xStep,y*yStep) );
    }

    // Form the vector of invNorms
    DistMatrix<Real,VR,STAR> invNorms(g);
    auto itCounts = 
        TriangularPseudospectrum
        ( U, shifts, invNorms, lanczos, krylovSize, reorthog, deflate, maxIts, 
          tol, progress );

    // Rearrange the vectors into grids
    DistMatrix<Int> itCountMap(g); 
    pspec::ReshapeIntoGrids
    ( xSize, ySize, invNorms, itCounts, invNormMap, itCountMap );
    return itCountMap;
}

template<typename F>
inline DistMatrix<Int>
HessenbergPseudospectrum
( const DistMatrix<F>& H, DistMatrix<BASE(F)>& invNormMap, 
  Complex<BASE(F)> center, BASE(F) xWidth, BASE(F) yWidth, Int xSize, Int ySize,
  bool lanczos=true, Int krylovSize=10, bool reorthog=true, bool deflate=true, 
  Int maxIts=1000, BASE(F) tol=1e-6, bool progress=false )
{
    DEBUG_ONLY(CallStackEntry cse("HessenbergPseudospectrum"))
    typedef Base<F> Real;
    typedef Complex<Real> C;
    const Grid& g = H.Grid();

    const Real xStep = xWidth/(xSize-1);
    const Real yStep = yWidth/(ySize-1);
    const C corner = center - C(xWidth/2,yWidth/2);
    DistMatrix<C,VR,STAR> shifts( xSize*ySize, 1, g );
    const Int numLocShifts = shifts.LocalHeight();
    for( Int jLoc=0; jLoc<numLocShifts; ++jLoc )
    {
        const Int j = shifts.ColShift() + jLoc*shifts.ColStride();
        const Int x = j / ySize;
        const Int y = j % ySize;
        shifts.SetLocal( jLoc, 0, corner+C(x*xStep,y*yStep) );
    }

    // Form the vector of invNorms
    DistMatrix<Real,VR,STAR> invNorms(g);
    auto itCounts = 
        HessenbergPseudospectrum
        ( H, shifts, invNorms, lanczos, krylovSize, reorthog, deflate, maxIts, 
          tol, progress );

    // Rearrange the vectors into grids
    DistMatrix<Int> itCountMap(g); 
    pspec::ReshapeIntoGrids
    ( xSize, ySize, invNorms, itCounts, invNormMap, itCountMap );
    return itCountMap;
}

template<typename F>
inline Matrix<Int>
Pseudospectrum
( const Matrix<F>& A, Matrix<BASE(F)>& invNormMap, 
  Complex<BASE(F)> center, BASE(F) xWidth, BASE(F) yWidth,
  Int xSize, Int ySize, bool schur=true, bool lanczos=true, Int krylovSize=10, 
  bool reorthog=true, bool deflate=true, Int maxIts=1000, BASE(F) tol=1e-6, 
  bool progress=false )
{
    DEBUG_ONLY(CallStackEntry cse("Pseudospectrum"))
    typedef Base<F> Real;
    typedef Complex<Real> C;

    const Real xStep = xWidth/(xSize-1);
    const Real yStep = yWidth/(ySize-1);
    const C corner = center - C(xWidth/2,yWidth/2);
    Matrix<C> shifts( xSize*ySize, 1, A.Grid() );
    for( Int j=0; j<xSize*ySize; ++j )
    {
        const Int x = j / ySize;
        const Int y = j % ySize;
        shifts.Set( j, 0, corner+C(x*xStep,y*yStep) );
    }

    // Form the vector of invNorms
    Matrix<Real> invNorms;
    auto itCounts = 
        Pseudospectrum
        ( A, shifts, invNorms, schur, lanczos, krylovSize, reorthog, deflate, 
          maxIts, tol, progress );

    // Rearrange the vectors into grids
    Matrix<Int> itCountMap; 
    pspec::ReshapeIntoGrids
    ( xSize, ySize, invNorms, itCounts, invNormMap, itCountMap );
    return itCountMap;
}

template<typename F>
inline DistMatrix<Int>
Pseudospectrum
( const DistMatrix<F>& A, DistMatrix<BASE(F)>& invNormMap, 
  Complex<BASE(F)> center, BASE(F) xWidth, BASE(F) yWidth, Int xSize, Int ySize,
  bool schur=false, bool lanczos=true, Int krylovSize=10, bool reorthog=true, 
  bool deflate=true, Int maxIts=1000, BASE(F) tol=1e-6, bool progress=false )
{
    DEBUG_ONLY(CallStackEntry cse("Pseudospectrum"))
    typedef Base<F> Real;
    typedef Complex<Real> C;
    const Grid& g = A.Grid();

    const Real xStep = xWidth/(xSize-1);
    const Real yStep = yWidth/(ySize-1);
    const C corner = center - C(xWidth/2,yWidth/2);
    DistMatrix<C,VR,STAR> shifts( xSize*ySize, 1, g );
    const Int numLocShifts = shifts.LocalHeight();
    for( Int jLoc=0; jLoc<numLocShifts; ++jLoc )
    {
        const Int j = shifts.ColShift() + jLoc*shifts.ColStride();
        const Int x = j / ySize;
        const Int y = j % ySize;
        shifts.SetLocal( jLoc, 0, corner+C(x*xStep,y*yStep) );
    }

    // Form the vector of invNorms
    DistMatrix<Real,VR,STAR> invNorms(g);
    auto itCounts =
        Pseudospectrum
        ( A, shifts, invNorms, schur, lanczos, krylovSize, reorthog, deflate, 
          maxIts, tol, progress );

    // Rearrange the vectors into grids
    DistMatrix<Int> itCountMap(g); 
    pspec::ReshapeIntoGrids
    ( xSize, ySize, invNorms, itCounts, invNormMap, itCountMap );
    return itCountMap;
}

template<typename F>
inline Matrix<Int>
TriangularPseudospectrum
( const Matrix<F>& U, Matrix<BASE(F)>& invNormMap, 
  Complex<BASE(F)> center,
  Int xSize, Int ySize, bool lanczos=true, Int krylovSize=10, 
  bool reorthog=true, bool deflate=true, Int maxIts=1000, BASE(F) tol=1e-6, 
  bool progress=false )
{
    DEBUG_ONLY(CallStackEntry cse("TriangularPseudospectrum"))

    auto diag = U.GetDiagonal();
    const Base<F> radius = MaxNorm( diag );
    const Base<F> oneNorm = OneNorm( U );

    // Essentially three cases are handled here:
    // 1) The zero matrix (force the pseudospectrum width to 1)
    // 2) Typical matrices (use a small multiple of the spectral radius)
    // 3) Highly non-normal matrices (e.g., triangular with zero main diagonal)
    Base<F> width;
    if( oneNorm == Base<F>(0) )
    {
        width = 1;
        if( progress )
            std::cout << "Setting width to 1 to handle zero matrix" 
                      << std::endl;
    }
    else if( radius >= 0.2*oneNorm )
    {
        width = 2.5*radius;
        if( progress )
            std::cout << "Setting width to " << width 
                      << " based on the spectral radius, " << radius 
                      << std::endl;
    }
    else
    {
        width = 0.8*oneNorm;
        if( progress )
            std::cout << "Setting width to " << width 
                      << " based on the one norm, " << oneNorm << std::endl;
    }

    return TriangularPseudospectrum
           ( U, invNormMap, center, width, width, xSize, ySize, 
             lanczos, krylovSize, reorthog, deflate, maxIts, tol, progress );
}

template<typename F>
inline Matrix<Int>
HessenbergPseudospectrum
( const Matrix<F>& H, Matrix<BASE(F)>& invNormMap, 
  Complex<BASE(F)> center,
  Int xSize, Int ySize, bool lanczos=true, Int krylovSize=10, 
  bool reorthog=true, bool deflate=true, Int maxIts=1000, BASE(F) tol=1e-6, 
  bool progress=false )
{
    DEBUG_ONLY(CallStackEntry cse("HessenbergPseudospectrum"))

    const Base<F> infNorm = InfinityNorm( H );
    const Base<F> oneNorm = OneNorm( H );
    Base<F> width;
    if( oneNorm == Base<F>(0) )
    {
        width = 1;
        if( progress )
            std::cout << "Setting width to 1 to handle zero matrix" 
                      << std::endl;
    }
    else
    {
        width = 0.8*Max(oneNorm,infNorm);
        if( progress )
            std::cout << "Setting width to " << width 
                      << " based on the one norm, " << oneNorm 
                      << ", and infinity norm, " << infNorm << std::endl;
    }

    return HessenbergPseudospectrum
           ( H, invNormMap, center, width, width, xSize, ySize, 
             lanczos, krylovSize, reorthog, deflate, maxIts, tol, progress );
}

template<typename F>
inline DistMatrix<Int>
TriangularPseudospectrum
( const DistMatrix<F>& U, DistMatrix<BASE(F)>& invNormMap, 
  Complex<BASE(F)> center, Int xSize, Int ySize,
  bool lanczos=true, Int krylovSize=10, bool reorthog=true, bool deflate=true, 
  Int maxIts=1000, BASE(F) tol=1e-6, bool progress=false )
{
    DEBUG_ONLY(CallStackEntry cse("TriangularPseudospectrum"))
    auto diag = U.GetDiagonal();
    const Base<F> radius = MaxNorm( diag );
    const Base<F> oneNorm = OneNorm( U );

    // Essentially three cases are handled here:
    // 1) The zero matrix (force the pseudospectrum width to 1)
    // 2) Typical matrices (use a small multiple of the spectral radius)
    // 3) Highly non-normal matrices (e.g., triangular with zero main diagonal)
    Base<F> width;
    if( oneNorm == Base<F>(0) && radius == Base<F>(0) )
    {
        width = 1;
        if( progress && U.Grid().Rank() == 0 )
            std::cout << "Setting width to 1 to handle zero matrix"
                      << std::endl;
    }
    else if( radius >= 0.2*oneNorm )
    {
        width = 2.5*radius;
        if( progress && U.Grid().Rank() == 0 )
            std::cout << "Setting width to " << width 
                      << " based on the spectral radius, " << radius 
                      << std::endl;
    }
    else
    {
        width = 0.8*oneNorm;
        if( progress && U.Grid().Rank() == 0 )
            std::cout << "Setting width to " << width
                      << " based on the one norm, " << oneNorm << std::endl;
    }

    return TriangularPseudospectrum
           ( U, invNormMap, center, width, width, xSize, ySize, 
             lanczos, krylovSize, reorthog, deflate, maxIts, tol, progress );
}

template<typename F>
inline DistMatrix<Int>
HessenbergPseudospectrum
( const DistMatrix<F>& H, DistMatrix<BASE(F)>& invNormMap, 
  Complex<BASE(F)> center, Int xSize, Int ySize,
  bool lanczos=true, Int krylovSize=10, bool reorthog=true, bool deflate=true, 
  Int maxIts=1000, BASE(F) tol=1e-6, bool progress=false )
{
    DEBUG_ONLY(CallStackEntry cse("HessenbergPseudospectrum"))

    const Base<F> oneNorm = OneNorm( H );
    const Base<F> infNorm = InfinityNorm( H );
    Base<F> width;
    if( oneNorm == Base<F>(0) )
    {
        width = 1;
        if( progress && H.Grid().Rank() == 0 )
            std::cout << "Setting width to 1 to handle zero matrix"
                      << std::endl;
    }
    else
    {
        width = 0.8*Max(oneNorm,infNorm);
        if( progress && H.Grid().Rank() == 0 )
            std::cout << "Setting width to " << width 
                      << " based on the one norm, " << oneNorm 
                      << ", and infinity norm, " << infNorm << std::endl;
    }

    return HessenbergPseudospectrum
           ( H, invNormMap, center, width, width, xSize, ySize, 
             lanczos, krylovSize, reorthog, deflate, maxIts, tol, progress );
}

template<typename F>
inline Matrix<Int>
Pseudospectrum
( const Matrix<F>& A, Matrix<BASE(F)>& invNormMap, 
  Complex<BASE(F)> center,
  Int xSize, Int ySize, bool schur=true, bool lanczos=true, Int krylovSize=10, 
  bool reorthog=true, bool deflate=true, Int maxIts=1000, BASE(F) tol=1e-6, 
  bool progress=false )
{
    DEBUG_ONLY(CallStackEntry cse("Pseudospectrum"))
    typedef Base<F> Real;
    typedef Complex<Real> C;

    const Int n = A.Height();
    Matrix<C> B( n, n );
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<n; ++i )
            B.Set( i, j, A.Get(i,j) );

    if( schur )
    {
        Matrix<C> w;
        schur::QR( B, w );
        return TriangularPseudospectrum
               ( B, invNormMap, center, xSize, ySize, 
                 lanczos, krylovSize, reorthog, deflate, maxIts, tol, 
                 progress );
    }
    else
    {
        Hessenberg( UPPER, B );
        return HessenbergPseudospectrum
               ( B, invNormMap, center, xSize, ySize, 
                 lanczos, krylovSize, reorthog, deflate, maxIts, tol, 
                 progress );
    }
}

template<typename F>
inline DistMatrix<Int>
Pseudospectrum
( const DistMatrix<F>& A, DistMatrix<BASE(F)>& invNormMap, 
  Complex<BASE(F)> center, Int xSize, Int ySize,
  bool schur=false, bool lanczos=true, Int krylovSize=10, bool reorthog=true, 
  bool deflate=true, Int maxIts=1000, BASE(F) tol=1e-6, bool progress=false )
{
    DEBUG_ONLY(CallStackEntry cse("Pseudospectrum"))
    typedef Base<F> Real;
    typedef Complex<Real> C;

    const Grid& g = A.Grid();
    DistMatrix<C> B(g);
    B.AlignWith( A );
    const Int n = A.Height();
    B.Resize( n, n );
    const Int mLocal = A.LocalHeight();
    const Int nLocal = A.LocalWidth();
    for( Int jLoc=0; jLoc<nLocal; ++jLoc )
        for( Int iLoc=0; iLoc<mLocal; ++iLoc )
            B.SetLocal( iLoc, jLoc, A.GetLocal(iLoc,jLoc) );

    if( schur )
    {
        // We don't actually need the Schur vectors, but SDC requires their
        // computation in order to form the full triangular factor
        DistMatrix<C> X(g);
        DistMatrix<C,VR,STAR> w(g);
        const bool formATR = true;
        // TODO: Expose these as options
        const Int cutoff = 256;
        const Int maxInnerIts = 2;
        const Int maxOuterIts = 10;
        const Base<F> signTol=tol/10;
        const Base<F> relTol=tol/10;
        const Base<F> spreadFactor=1e-6;
        const bool random=true;
        schur::SDC
        ( B, w, X, formATR, cutoff, maxInnerIts, maxOuterIts, signTol, relTol,
          spreadFactor, random, progress );
        X.Empty();
 
        return TriangularPseudospectrum
               ( B, invNormMap, center, xSize, ySize,
                 lanczos, krylovSize, reorthog, deflate, maxIts, tol, 
                 progress );
    }
    else
    {
        Hessenberg( UPPER, B );
        return HessenbergPseudospectrum
               ( B, invNormMap, center, xSize, ySize, 
                 lanczos, krylovSize, reorthog, deflate, maxIts, tol, 
                 progress );
    }
}

} // namespace elem

#endif // ifndef ELEM_PSEUDOSPECTRUM_HPP
