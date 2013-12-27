/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_PSEUDOSPECTRUM_HPP
#define ELEM_LAPACK_PSEUDOSPECTRUM_HPP

#include "elemental/lapack-like/Schur.hpp"

#include "elemental/lapack-like/Pseudospectrum/ShiftedTrsm.hpp"
#include "elemental/lapack-like/Pseudospectrum/Power.hpp"
#include "elemental/lapack-like/Pseudospectrum/Lanczos.hpp"

namespace elem {

namespace pspec {

template<typename T>
inline void
ReshapeIntoGrids
( Int xSize, Int ySize,
  const Matrix<T>& invNorms,   const Matrix<Int>& itCounts,
        Matrix<T>& invNormMap,       Matrix<Int>& itCountMap )
{
#if 0    
    invNormMap.ResizeTo( xSize, ySize );
    itCountMap.ResizeTo( xSize, ySize );
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
    invNormMap.ResizeTo( xSize, ySize, xSize );
    itCountMap.ResizeTo( xSize, ySize, xSize );
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
    invNormMap.ResizeTo( xSize, ySize );
    itCountMap.ResizeTo( xSize, ySize );
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
  Matrix<BASE(F)>& invNorms, bool lanczos=true, bool deflate=true,
  Int maxIts=1000, BASE(F) tol=1e-6, bool progress=false )
{
    DEBUG_ONLY(CallStackEntry cse("TriangularPseudospectrum"))
    typedef Base<F> Real;
    typedef Complex<Real> C;
    const Int n = U.Height();

    Matrix<C> UCpx;
    if( IsComplex<F>::val )
        UCpx = LockedView( U );
    else
    {
        UCpx.ResizeTo( n, n );
        for( Int j=0; j<n; ++j )
            for( Int i=0; i<n; ++i )
                UCpx.Set( i, j, U.Get(i,j) );
    }

    Matrix<Int> itCounts;
    if( lanczos )
        itCounts = 
           pspec::TriangularLanczos
           ( UCpx, shifts, invNorms, deflate, maxIts, tol, progress );
    else
        itCounts =
           pspec::TriangularPower
           ( UCpx, shifts, invNorms, deflate, maxIts, tol, progress );

    return itCounts;
}

template<typename F>
inline DistMatrix<Int,VR,STAR>
TriangularPseudospectrum
( const DistMatrix<F>& U, const DistMatrix<Complex<BASE(F)>,VR,STAR>& shifts,
  DistMatrix<BASE(F),VR,STAR>& invNorms, bool lanczos=true, bool deflate=true,
  Int maxIts=1000, BASE(F) tol=1e-6, bool progress=false )
{
    DEBUG_ONLY(CallStackEntry cse("TriangularPseudospectrum"))
    typedef Base<F> Real;
    typedef Complex<Real> C;
    const Int n = U.Height();
    const Int mLocal = U.LocalHeight();
    const Int nLocal = U.LocalWidth();
    const Grid& g = U.Grid();

    DistMatrix<C> UCpx(g);
    if( IsComplex<F>::val )
        UCpx = LockedView( U );
    else
    {
        UCpx.AlignWith( U );
        UCpx.ResizeTo( n, n );
        for( Int jLoc=0; jLoc<nLocal; ++jLoc )
            for( Int iLoc=0; iLoc<mLocal; ++iLoc )
                UCpx.SetLocal( iLoc, jLoc, U.GetLocal(iLoc,jLoc) );
    }

    DistMatrix<Int,VR,STAR> itCounts(g);
    if( lanczos )
        itCounts = 
           pspec::TriangularLanczos
           ( UCpx, shifts, invNorms, deflate, maxIts, tol, progress );
    else
        itCounts =
           pspec::TriangularPower
           ( UCpx, shifts, invNorms, deflate, maxIts, tol, progress );

    return itCounts;
}

template<typename F>
inline Matrix<Int>
Pseudospectrum
( const Matrix<F>& A, const Matrix<Complex<BASE(F)> >& shifts, 
  Matrix<BASE(F)>& invNorms, bool lanczos=true, bool deflate=true,
  Int maxIts=1000, BASE(F) tol=1e-6, bool progress=false )
{
    DEBUG_ONLY(CallStackEntry cse("Pseudospectrum"))
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
           ( U, shifts, invNorms, lanczos, deflate, maxIts, tol, progress );
}

template<typename F>
inline DistMatrix<Int,VR,STAR>
Pseudospectrum
( const DistMatrix<F>& A, const DistMatrix<Complex<BASE(F)>,VR,STAR>& shifts,
  DistMatrix<BASE(F),VR,STAR>& invNorms, bool lanczos=true, bool deflate=true,
  Int maxIts=1000, BASE(F) tol=1e-6, bool progress=false )
{
    DEBUG_ONLY(CallStackEntry cse("Pseudospectrum"))
    typedef Base<F> Real;
    typedef Complex<Real> C;
    const Int n = A.Height();
    const Int mLocal = A.LocalHeight();
    const Int nLocal = A.LocalWidth();
    const Grid& g = A.Grid();

    DistMatrix<C> U(g);
    U.AlignWith( A );
    U.ResizeTo( n, n );
    for( Int jLoc=0; jLoc<nLocal; ++jLoc )
        for( Int iLoc=0; iLoc<mLocal; ++iLoc )
            U.SetLocal( iLoc, jLoc, A.GetLocal(iLoc,jLoc) );

    // We don't actually need the Schur vectors, but SDC requires their 
    // computation in order to form the full triangular factor
    DistMatrix<C> X(g);
    DistMatrix<C,VR,STAR> w(g);
    schur::SDC( U, w, X );
    X.Empty();

    return TriangularPseudospectrum
           ( U, shifts, invNorms, lanczos, deflate, maxIts, tol, progress );
}

template<typename F>
inline Matrix<Int>
TriangularPseudospectrum
( const Matrix<F>& U, Matrix<BASE(F)>& invNormMap, 
  Complex<BASE(F)> center, BASE(F) xWidth, BASE(F) yWidth,
  Int xSize, Int ySize, bool lanczos=true, bool deflate=true, 
  Int maxIts=1000, BASE(F) tol=1e-6, bool progress=false )
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
        ( U, shifts, invNorms, lanczos, deflate, maxIts, tol, progress );

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
  bool lanczos=true, bool deflate=true, 
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
        ( U, shifts, invNorms, lanczos, deflate, maxIts, tol, progress );

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
  Int xSize, Int ySize, bool lanczos=true, bool deflate=true, 
  Int maxIts=1000, BASE(F) tol=1e-6, bool progress=false )
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
        ( A, shifts, invNorms, lanczos, deflate, maxIts, tol, progress );

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
  bool lanczos=true, bool deflate=true, 
  Int maxIts=1000, BASE(F) tol=1e-6, bool progress=false )
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
        ( A, shifts, invNorms, lanczos, deflate, maxIts, tol, progress );

    // Rearrange the vectors into grids
    DistMatrix<Int> itCountMap(g); 
    pspec::ReshapeIntoGrids
    ( xSize, ySize, invNorms, itCounts, invNormMap, itCountMap );
    return itCountMap;
}

// TODO: Chunked pseudospectrum driver
// TODO: Spectrally-centered pseudospectrum driver
// TODO: operator= which can convert between datatypes

} // namespace elem

#endif // ifndef ELEM_LAPACK_PSEUDOSPECTRUM_HPP
