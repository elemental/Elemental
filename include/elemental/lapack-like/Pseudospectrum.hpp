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

    Matrix<Int> itCounts;
    if( lanczos )
        itCounts = 
           pspec::TriangularLanczos
           ( U, shifts, invNorms, deflate, maxIts, tol, progress );
    else
        itCounts =
           pspec::TriangularPower
           ( U, shifts, invNorms, deflate, maxIts, tol, progress );

    return itCounts;
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

    DistMatrix<Int,VR,STAR> itCounts(g);
    if( lanczos )
        itCounts = 
           pspec::TriangularLanczos
           ( U, shifts, invNorms, deflate, maxIts, tol, progress );
    else
        itCounts =
           pspec::TriangularPower
           ( U, shifts, invNorms, deflate, maxIts, tol, progress );

    return itCounts;
}

template<typename F>
inline Matrix<Int>
Pseudospectrum
( const Matrix<F>& A, Matrix<BASE(F)>& invNormMap, 
  Complex<BASE(F)> center, BASE(F) halfWidth,
  Int xSize, Int ySize, bool lanczos=true, bool deflate=true, 
  Int maxIts=1000, BASE(F) tol=1e-6, bool progress=false )
{
    DEBUG_ONLY(CallStackEntry cse("Pseudospectrum"))
    typedef Base<F> Real;
    typedef Complex<Real> C;

    if( halfWidth == Real(0) )
        halfWidth = FrobeniusNorm( A );
    const Real xStep = 2*halfWidth/(xSize-1);
    const Real yStep = 2*halfWidth/(ySize-1);
    const C corner = center - C(halfWidth,halfWidth);
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

    // Rearrange the vector into a grid 
    invNormMap.ResizeTo( xSize, ySize );
    Matrix<Int> itCountMap( xSize, ySize );
    for( Int j=0; j<xSize; ++j )
    {
        auto normGridSub = View( invNormMap, 0, j, ySize, 1 );
        auto countGridSub = View( itCountMap, 0, j, ySize, 1 );
        auto shiftSub = LockedView( invNorms, j*ySize, 0, ySize, 1 );
        auto countSub = LockedView( itCounts, j*ySize, 0, ySize, 1 );
        normGridSub = shiftSub;
        countGridSub = countSub;
    }

    return itCountMap;
}

template<typename F>
inline DistMatrix<Int>
Pseudospectrum
( const DistMatrix<F>& A, DistMatrix<BASE(F)>& invNormMap, 
  Complex<BASE(F)> center, BASE(F) halfWidth, Int xSize, Int ySize, 
  bool lanczos=true, bool deflate=true, 
  Int maxIts=1000, BASE(F) tol=1e-6, bool progress=false )
{
    DEBUG_ONLY(CallStackEntry cse("Pseudospectrum"))
    typedef Base<F> Real;
    typedef Complex<Real> C;
    const Grid& g = A.Grid();

    if( halfWidth == Real(0) )
        halfWidth = FrobeniusNorm( A );
    const Real xStep = 2*halfWidth/(xSize-1);
    const Real yStep = 2*halfWidth/(ySize-1);
    const C corner = center - C(halfWidth,halfWidth);
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

    // Rearrange the vector into a grid 
    invNormMap.ResizeTo( xSize, ySize );
    DistMatrix<Int> itCountMap( xSize, ySize, g );
    for( Int j=0; j<xSize; ++j )
    {
        auto normGridSub = View( invNormMap, 0, j, ySize, 1 );
        auto countGridSub = View( itCountMap, 0, j, ySize, 1 );
        auto shiftSub = LockedView( invNorms, j*ySize, 0, ySize, 1 );
        auto countSub = LockedView( itCounts, j*ySize, 0, ySize, 1 );
        normGridSub = shiftSub;
        countGridSub = countSub;
    }

    return itCountMap;
}

template<typename F>
inline Matrix<Int>
TriangularPseudospectrum
( const Matrix<F>& U, Matrix<BASE(F)>& invNormMap, 
  Complex<BASE(F)> center, BASE(F) halfWidth,
  Int xSize, Int ySize, bool lanczos=true, bool deflate=true, 
  Int maxIts=1000, BASE(F) tol=1e-6, bool progress=false )
{
    DEBUG_ONLY(CallStackEntry cse("TriangularPseudospectrum"))
    typedef Base<F> Real;
    typedef Complex<Real> C;

    if( halfWidth == Real(0) )
        halfWidth = FrobeniusNorm( U );
    const Real xStep = 2*halfWidth/(xSize-1);
    const Real yStep = 2*halfWidth/(ySize-1);
    const C corner = center - C(halfWidth,halfWidth);
    Matrix<C> shifts( xSize*ySize, 1, U.Grid() );
    for( Int j=0; j<xSize*ySize; ++j )
    {
        const Int x = j / ySize;
        const Int y = j % ySize;
        shifts.Set( j, 0, corner+C(x*xStep,y*yStep) );
    }

    // Form the vector of invNorms
    Matrix<Real> invNorms;
    Matrix<Int> itCounts;
    if( lanczos )
        itCounts = 
           pspec::TriangularLanczos
           ( U, shifts, invNorms, deflate, maxIts, tol, progress );
    else
        itCounts =
           pspec::TriangularPower
           ( U, shifts, invNorms, deflate, maxIts, tol, progress );

    // Rearrange the vector into a grid 
    invNormMap.ResizeTo( xSize, ySize );
    Matrix<Int> itCountMap( xSize, ySize );
    for( Int j=0; j<xSize; ++j )
    {
        auto normGridSub = View( invNormMap, 0, j, ySize, 1 );
        auto countGridSub = View( itCountMap, 0, j, ySize, 1 );
        auto shiftSub = LockedView( invNorms, j*ySize, 0, ySize, 1 );
        auto countSub = LockedView( itCounts, j*ySize, 0, ySize, 1 );
        normGridSub = shiftSub;
        countGridSub = countSub;
    }

    return itCountMap;
}

template<typename F>
inline DistMatrix<Int>
TriangularPseudospectrum
( const DistMatrix<F>& U, DistMatrix<BASE(F)>& invNormMap, 
  Complex<BASE(F)> center, BASE(F) halfWidth, Int xSize, Int ySize, 
  bool lanczos=true, bool deflate=true, 
  Int maxIts=1000, BASE(F) tol=1e-6, bool progress=false )
{
    DEBUG_ONLY(CallStackEntry cse("TriangularPseudospectrum"))
    typedef Base<F> Real;
    typedef Complex<Real> C;
    const Grid& g = U.Grid();

    if( halfWidth == Real(0) )
        halfWidth = FrobeniusNorm( U );
    const Real xStep = 2*halfWidth/(xSize-1);
    const Real yStep = 2*halfWidth/(ySize-1);
    const C corner = center - C(halfWidth,halfWidth);
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
    DistMatrix<Int,VR,STAR> itCounts(g); 
    if( lanczos )
        itCounts = 
           pspec::TriangularLanczos
           ( U, shifts, invNorms, deflate, maxIts, tol, progress );
    else
        itCounts =
           pspec::TriangularPower
           ( U, shifts, invNorms, deflate, maxIts, tol, progress );

    // Rearrange the vector into a grid 
    invNormMap.ResizeTo( xSize, ySize );
    DistMatrix<Int> itCountMap( xSize, ySize, g );
    for( Int j=0; j<xSize; ++j )
    {
        auto normGridSub = View( invNormMap, 0, j, ySize, 1 );
        auto countGridSub = View( itCountMap, 0, j, ySize, 1 );
        auto shiftSub = LockedView( invNorms, j*ySize, 0, ySize, 1 );
        auto countSub = LockedView( itCounts, j*ySize, 0, ySize, 1 );
        normGridSub = shiftSub;
        countGridSub = countSub;
    }

    return itCountMap;
}

} // namespace elem

#endif // ifndef ELEM_LAPACK_PSEUDOSPECTRUM_HPP
