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

#include ELEM_COPY_INC
#include ELEM_MAXNORM_INC
#include ELEM_ONENORM_INC
#include ELEM_TWONORMESTIMATE_INC
#include ELEM_HESSENBERG_INC
#include ELEM_SCHUR_INC

#include "./Pseudospectrum/Util.hpp"
#include "./Pseudospectrum/Power.hpp"
#include "./Pseudospectrum/Lanczos.hpp"
#include "./Pseudospectrum/IRA.hpp"
#include "./Pseudospectrum/IRL.hpp"
#include "./Pseudospectrum/Analytic.hpp"

namespace elem {

template<typename F>
inline Matrix<Int>
TriangularPseudospectrum
( const Matrix<F>& U, const Matrix<Complex<BASE(F)> >& shifts, 
  Matrix<BASE(F)>& invNorms, 
  PseudospecCtrl<BASE(F)> psCtrl=PseudospecCtrl<BASE(F)>() )
{
    DEBUG_ONLY(CallStackEntry cse("TriangularPseudospectrum"))
    typedef Base<F> Real;
    typedef Complex<Real> C;

    psCtrl.schur = true;

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
    if( pspec::NumericallyNormal( UCpx, psCtrl.tol ) )
    {
        Matrix<Int> itCounts;
        if( psCtrl.progress )
            std::cout << "Matrix was numerically normal" << std::endl;
        auto w = UCpx.GetDiagonal();
        pspec::Analytic( w, shifts, invNorms, psCtrl.snapCtrl );
        Zeros( itCounts, shifts.Height(), 1 );        
        return itCounts;
    }

    if( psCtrl.arnoldi )
    {
        if( psCtrl.basisSize > 1 )
            return pspec::IRA( UCpx, shifts, invNorms, psCtrl );
        else
            return pspec::Lanczos( UCpx, shifts, invNorms, psCtrl );
    }
    else
        return pspec::Power( UCpx, shifts, invNorms, psCtrl );
}

template<typename Real>
inline Matrix<Int>
QuasiTriangularPseudospectrum
( const Matrix<Real>& U, const Matrix<Complex<Real> >& w,
  const Matrix<Complex<Real> >& shifts, 
  Matrix<Real>& invNorms, 
  PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() )
{
    DEBUG_ONLY(CallStackEntry cse("QuasiTriangularPseudospectrum"))

    psCtrl.schur = true;

    // Check if the off-diagonal is sufficiently small; if so, compute the 
    // pseudospectrum analytically from the eigenvalues. This also takes care
    // of the case where the matrix is a constant multiple of the identity 
    // matrix, which, after shifting, can lead to the zero matrix, which would 
    // cause problems for the Lanczos convergence criteria.
    if( pspec::NumericallyNormal( U, w, psCtrl.tol ) )
    {
        Matrix<Int> itCounts;
        if( psCtrl.progress )
            std::cout << "Matrix was numerically normal" << std::endl;
        pspec::Analytic( w, shifts, invNorms, psCtrl.snapCtrl );
        Zeros( itCounts, shifts.Height(), 1 );        
        return itCounts;
    }

    return pspec::IRA( U, shifts, invNorms, psCtrl );
}

template<typename F>
inline Matrix<Int>
HessenbergPseudospectrum
( const Matrix<F>& H, const Matrix<Complex<BASE(F)> >& shifts, 
  Matrix<BASE(F)>& invNorms, 
  PseudospecCtrl<BASE(F)> psCtrl=PseudospecCtrl<BASE(F)>() )
{
    DEBUG_ONLY(CallStackEntry cse("HessenbergPseudospectrum"))
    typedef Base<F> Real;
    typedef Complex<Real> C;

    psCtrl.schur = false;

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
    //       triangular version of Pseudospectrum?
    if( psCtrl.arnoldi )
    {
        if( psCtrl.basisSize > 1 )
            return pspec::IRA( HCpx, shifts, invNorms, psCtrl );
        else
            return pspec::Lanczos( HCpx, shifts, invNorms, psCtrl );
    }
    else
        return pspec::Power( HCpx, shifts, invNorms, psCtrl );
}

template<typename F>
inline DistMatrix<Int,VR,STAR>
TriangularPseudospectrum
( const DistMatrix<F>& U, const DistMatrix<Complex<BASE(F)>,VR,STAR>& shifts,
  DistMatrix<BASE(F),VR,STAR>& invNorms, 
  PseudospecCtrl<BASE(F)> psCtrl=PseudospecCtrl<BASE(F)>() )
{
    DEBUG_ONLY(CallStackEntry cse("TriangularPseudospectrum"))
    typedef Base<F> Real;
    typedef Complex<Real> C;

    psCtrl.schur = true;

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
    if( pspec::NumericallyNormal( UCpx, psCtrl.tol ) )
    {
        DistMatrix<Int,VR,STAR> itCounts(g);
        if( psCtrl.progress && U.Grid().Rank() == 0 )
            std::cout << "Matrix was numerically normal" << std::endl;
        auto w = UCpx.GetDiagonal();
        DistMatrix<C,STAR,STAR> w_STAR_STAR( w );
        pspec::Analytic( w_STAR_STAR, shifts, invNorms, psCtrl.snapCtrl );
        itCounts.AlignWith( shifts );
        Zeros( itCounts, shifts.Height(), 1 );
        return itCounts;
    }

    if( psCtrl.arnoldi )
    {
        if( psCtrl.basisSize > 1 )
            return pspec::IRA( UCpx, shifts, invNorms, psCtrl );
        else
            return pspec::Lanczos( UCpx, shifts, invNorms, psCtrl );
    }
    else
        return pspec::Power( UCpx, shifts, invNorms, psCtrl );
}

template<typename Real>
inline DistMatrix<Int,VR,STAR>
QuasiTriangularPseudospectrum
( const DistMatrix<Real>& U, const DistMatrix<Complex<Real>,VR,STAR>& w,
  const DistMatrix<Complex<Real>,VR,STAR>& shifts,
  DistMatrix<Real,VR,STAR>& invNorms, 
  PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() )
{
    DEBUG_ONLY(CallStackEntry cse("QuasiTriangularPseudospectrum"))
    const Grid& g = U.Grid();

    psCtrl.schur = true;

    // Check if the off-diagonal is sufficiently small; if so, compute the 
    // pseudospectrum analytically from the eigenvalues. This also takes care
    // of the case where the matrix is a constant multiple of the identity 
    // matrix, which, after shifting, can lead to the zero matrix, which would 
    // cause problems for the Lanczos convergence criteria.
    if( pspec::NumericallyNormal( U, w, psCtrl.tol ) )
    {
        DistMatrix<Int,VR,STAR> itCounts(g);
        if( psCtrl.progress && U.Grid().Rank() == 0 )
            std::cout << "Matrix was numerically normal" << std::endl;
        DistMatrix<Complex<Real>,STAR,STAR> w_STAR_STAR( w );
        pspec::Analytic( w_STAR_STAR, shifts, invNorms, psCtrl.snapCtrl );
        itCounts.AlignWith( shifts );
        Zeros( itCounts, shifts.Height(), 1 );
        return itCounts;
    }

    return pspec::IRA( U, shifts, invNorms, psCtrl );
}

template<typename F>
inline DistMatrix<Int,VR,STAR>
HessenbergPseudospectrum
( const DistMatrix<F>& H, const DistMatrix<Complex<BASE(F)>,VR,STAR>& shifts,
  DistMatrix<BASE(F),VR,STAR>& invNorms, 
  PseudospecCtrl<BASE(F)> psCtrl=PseudospecCtrl<BASE(F)>() )
{
    DEBUG_ONLY(CallStackEntry cse("HessenbergPseudospectrum"))
    typedef Base<F> Real;
    typedef Complex<Real> C;

    psCtrl.schur = false;

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
    if( psCtrl.arnoldi )
    {
        if( psCtrl.basisSize > 1 )
            return pspec::IRA( HCpx, shifts, invNorms, psCtrl );
        else
            return pspec::Lanczos( HCpx, shifts, invNorms, psCtrl );
    }
    else
        return pspec::Power( HCpx, shifts, invNorms, psCtrl );
}

template<typename Real>
inline Matrix<Int>
Pseudospectrum
( const Matrix<Real>& A, const Matrix<Complex<Real> >& shifts, 
  Matrix<Real>& invNorms, 
  PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() )
{
    DEBUG_ONLY(CallStackEntry cse("Pseudospectrum"))
    typedef Complex<Real> C;

    if( psCtrl.forceComplexSchur )
    {
        Matrix<C> ACpx;
        Copy( A, ACpx );
        return Pseudospectrum( ACpx, shifts, invNorms, psCtrl );
    }

    if( psCtrl.schur )
    {
        Matrix<Real> U( A );
        Matrix<C> w;
        const bool fullTriangle = true;
        schur::QR( U, w, fullTriangle );
        if( psCtrl.forceComplexPs )
        {
            Matrix<C> UCpx;
            schur::RealToComplex( U, UCpx );
            return TriangularPseudospectrum( UCpx, shifts, invNorms, psCtrl );
        }
        return QuasiTriangularPseudospectrum( U, w, shifts, invNorms, psCtrl );
    }
    else
        LogicError("Real Hessenberg algorithm not yet supported");
}

template<typename Real>
inline Matrix<Int>
Pseudospectrum
( const Matrix<Complex<Real> >& A, const Matrix<Complex<Real> >& shifts, 
  Matrix<Real>& invNorms, PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() )
{
    DEBUG_ONLY(CallStackEntry cse("Pseudospectrum"))
    typedef Complex<Real> C;

    Matrix<C> U( A );
    if( psCtrl.schur )
    {
        Matrix<C> w;
        const bool fullTriangle = true;
        schur::QR( U, w, fullTriangle );
        return TriangularPseudospectrum( U, shifts, invNorms, psCtrl );
    }
    else
    {
        Hessenberg( UPPER, U );
        return HessenbergPseudospectrum( U, shifts, invNorms, psCtrl );
    }
}

template<typename Real>
inline DistMatrix<Int,VR,STAR>
Pseudospectrum
( const DistMatrix<Real>& A, const DistMatrix<Complex<Real>,VR,STAR>& shifts,
  DistMatrix<Real,VR,STAR>& invNorms, 
  PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() )
{
    DEBUG_ONLY(CallStackEntry cse("Pseudospectrum"))
    typedef Complex<Real> C;
    const Grid& g = A.Grid();

    if( psCtrl.forceComplexSchur )
    {
        DistMatrix<C> ACpx(g);
        Copy( A, ACpx );
        return Pseudospectrum( ACpx, shifts, invNorms, psCtrl );
    }

    if( psCtrl.schur )
    {
        DistMatrix<Real> U( A );
        DistMatrix<C,VR,STAR> w(g);
        const bool fullTriangle = true;
#ifdef ELEM_HAVE_SCALAPACK
        schur::QR( U, w, fullTriangle );
#else
        // We don't actually need the Schur vectors, but SDC requires their 
        // computation in order to form the full triangular factor
        DistMatrix<Real> X(g);
        schur::SDC( U, w, X, fullTriangle, psCtrl.sdcCtrl );
        X.Empty();
#endif
        if( psCtrl.forceComplexPs )
        {
            DistMatrix<C> UCpx(g);
            schur::RealToComplex( U, UCpx );
            return TriangularPseudospectrum( UCpx, shifts, invNorms, psCtrl );
        }
        return QuasiTriangularPseudospectrum( U, w, shifts, invNorms, psCtrl );
    }
    else
        LogicError("Real Hessenberg algorithm not yet supported");
}

template<typename Real>
inline DistMatrix<Int,VR,STAR>
Pseudospectrum
( const DistMatrix<Complex<Real> >& A, 
  const DistMatrix<Complex<Real>,VR,STAR>& shifts,
  DistMatrix<Real,VR,STAR>& invNorms, 
  PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() )
{
    DEBUG_ONLY(CallStackEntry cse("Pseudospectrum"))
    typedef Complex<Real> C;

    const Grid& g = A.Grid();
    DistMatrix<C> U( A );

    if( psCtrl.schur )
    {
        DistMatrix<C,VR,STAR> w(g);
        const bool fullTriangle = true;
#ifdef ELEM_HAVE_SCALAPACK
        schur::QR( U, w, fullTriangle );
#else
        // We don't actually need the Schur vectors, but SDC requires their 
        // computation in order to form the full triangular factor
        DistMatrix<C> X(g);
        schur::SDC( U, w, X, fullTriangle, psCtrl.sdcCtrl );
        X.Empty();
#endif
        return TriangularPseudospectrum( U, shifts, invNorms, psCtrl );
    }
    else
    {
        Hessenberg( UPPER, U );
        return HessenbergPseudospectrum( U, shifts, invNorms, psCtrl );
    }
}

// Treat each pixel as being located a cell center and tesselate a box with
// said square cells
template<typename F>
inline Matrix<Int>
TriangularPseudospectrum
( const Matrix<F>& U, Matrix<BASE(F)>& invNormMap, 
  Complex<BASE(F)> center, BASE(F) realWidth, BASE(F) imagWidth,
  Int realSize, Int imagSize, 
  PseudospecCtrl<BASE(F)> psCtrl=PseudospecCtrl<BASE(F)>() )
{
    DEBUG_ONLY(CallStackEntry cse("TriangularPseudospectrum"))
    typedef Base<F> Real;
    typedef Complex<Real> C;

    psCtrl.snapCtrl.realSize = realSize;
    psCtrl.snapCtrl.imagSize = imagSize;

    const Real realStep = realWidth/realSize;
    const Real imagStep = imagWidth/imagSize;
    const C corner = center + C(-realWidth/2,imagWidth/2);
    Matrix<C> shifts( realSize*imagSize, 1, U.Grid() );
    for( Int j=0; j<realSize*imagSize; ++j )
    {
        const Int x = j / imagSize;
        const Int y = j % imagSize;
        shifts.Set( j, 0, corner+C((x+0.5)*realStep,-(y+0.5)*imagStep) );
    }

    // Form the vector of invNorms
    Matrix<Real> invNorms;
    auto itCounts = 
        TriangularPseudospectrum( U, shifts, invNorms, psCtrl );

    // Rearrange the vectors into grids
    Matrix<Int> itCountMap; 
    pspec::ReshapeIntoGrid( realSize, imagSize, invNorms, invNormMap );
    pspec::ReshapeIntoGrid( realSize, imagSize, itCounts, itCountMap );
    return itCountMap;
}

template<typename Real>
inline Matrix<Int>
QuasiTriangularPseudospectrum
( const Matrix<Real>& U, const Matrix<Complex<Real> >& w, 
  Matrix<Real>& invNormMap, 
  Complex<Real> center, Real realWidth, Real imagWidth,
  Int realSize, Int imagSize, 
  PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() )
{
    DEBUG_ONLY(CallStackEntry cse("QuasiTriangularPseudospectrum"))
    typedef Complex<Real> C;

    psCtrl.snapCtrl.realSize = realSize;
    psCtrl.snapCtrl.imagSize = imagSize;

    const Real realStep = realWidth/realSize;
    const Real imagStep = imagWidth/imagSize;
    const C corner = center + C(-realWidth/2,imagWidth/2);
    Matrix<C> shifts( realSize*imagSize, 1, U.Grid() );
    for( Int j=0; j<realSize*imagSize; ++j )
    {
        const Int x = j / imagSize;
        const Int y = j % imagSize;
        shifts.Set( j, 0, corner+C((x+0.5)*realStep,-(y+0.5)*imagStep) );
    }

    // Form the vector of invNorms
    Matrix<Real> invNorms;
    auto itCounts = 
        QuasiTriangularPseudospectrum( U, w, shifts, invNorms, psCtrl );

    // Rearrange the vectors into grids
    Matrix<Int> itCountMap; 
    pspec::ReshapeIntoGrid( realSize, imagSize, invNorms, invNormMap );
    pspec::ReshapeIntoGrid( realSize, imagSize, itCounts, itCountMap );
    return itCountMap;
}

template<typename F>
inline Matrix<Int>
HessenbergPseudospectrum
( const Matrix<F>& H, Matrix<BASE(F)>& invNormMap, 
  Complex<BASE(F)> center, BASE(F) realWidth, BASE(F) imagWidth,
  Int realSize, Int imagSize, 
  PseudospecCtrl<BASE(F)> psCtrl=PseudospecCtrl<BASE(F)>() )
{
    DEBUG_ONLY(CallStackEntry cse("HessenbergPseudospectrum"))
    typedef Base<F> Real;
    typedef Complex<Real> C;

    psCtrl.snapCtrl.realSize = realSize;
    psCtrl.snapCtrl.imagSize = imagSize;

    const Real realStep = realWidth/realSize;
    const Real imagStep = imagWidth/imagSize;
    const C corner = center + C(-realWidth/2,imagWidth/2);
    Matrix<C> shifts( realSize*imagSize, 1, H.Grid() );
    for( Int j=0; j<realSize*imagSize; ++j )
    {
        const Int x = j / imagSize;
        const Int y = j % imagSize;
        shifts.Set( j, 0, corner+C((x+0.5)*realStep,-(y+0.5)*imagStep) );
    }

    // Form the vector of invNorms
    Matrix<Real> invNorms;
    auto itCounts = HessenbergPseudospectrum( H, shifts, invNorms, psCtrl );

    // Rearrange the vectors into grids
    Matrix<Int> itCountMap; 
    pspec::ReshapeIntoGrid( realSize, imagSize, invNorms, invNormMap );
    pspec::ReshapeIntoGrid( realSize, imagSize, itCounts, itCountMap );
    return itCountMap;
}

template<typename F>
inline DistMatrix<Int>
TriangularPseudospectrum
( const DistMatrix<F>& U, DistMatrix<BASE(F)>& invNormMap, 
  Complex<BASE(F)> center, BASE(F) realWidth, BASE(F) imagWidth, 
  Int realSize, Int imagSize, 
  PseudospecCtrl<BASE(F)> psCtrl=PseudospecCtrl<BASE(F)>() )
{
    DEBUG_ONLY(CallStackEntry cse("TriangularPseudospectrum"))
    typedef Base<F> Real;
    typedef Complex<Real> C;
    const Grid& g = U.Grid();

    psCtrl.snapCtrl.realSize = realSize;
    psCtrl.snapCtrl.imagSize = imagSize;

    const Real realStep = realWidth/realSize;
    const Real imagStep = imagWidth/imagSize;
    const C corner = center + C(-realWidth/2,imagWidth/2);
    DistMatrix<C,VR,STAR> shifts( realSize*imagSize, 1, g );
    const Int numLocShifts = shifts.LocalHeight();
    for( Int iLoc=0; iLoc<numLocShifts; ++iLoc )
    {
        const Int i = shifts.GlobalRow(iLoc);
        const Int x = i / imagSize;
        const Int y = i % imagSize;
        shifts.SetLocal
        ( iLoc, 0, corner+C((x+0.5)*realStep,-(y+0.5)*imagStep) );
    }

    // Form the vector of invNorms
    DistMatrix<Real,VR,STAR> invNorms(g);
    auto itCounts = TriangularPseudospectrum( U, shifts, invNorms, psCtrl );

    // Rearrange the vectors into grids
    DistMatrix<Int> itCountMap(g); 
    pspec::ReshapeIntoGrid( realSize, imagSize, invNorms, invNormMap );
    pspec::ReshapeIntoGrid( realSize, imagSize, itCounts, itCountMap );
    return itCountMap;
}

template<typename Real>
inline DistMatrix<Int>
QuasiTriangularPseudospectrum
( const DistMatrix<Real>& U, const DistMatrix<Complex<Real>,VR,STAR>& w,
  DistMatrix<Real>& invNormMap, 
  Complex<Real> center, Real realWidth, Real imagWidth, 
  Int realSize, Int imagSize, 
  PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() )
{
    DEBUG_ONLY(CallStackEntry cse("QuasiTriangularPseudospectrum"))
    typedef Complex<Real> C;
    const Grid& g = U.Grid();

    psCtrl.snapCtrl.realSize = realSize;
    psCtrl.snapCtrl.imagSize = imagSize;

    const Real realStep = realWidth/realSize;
    const Real imagStep = imagWidth/imagSize;
    const C corner = center + C(-realWidth/2,imagWidth/2);
    DistMatrix<C,VR,STAR> shifts( realSize*imagSize, 1, g );
    const Int numLocShifts = shifts.LocalHeight();
    for( Int iLoc=0; iLoc<numLocShifts; ++iLoc )
    {
        const Int i = shifts.GlobalRow(iLoc);
        const Int x = i / imagSize;
        const Int y = i % imagSize;
        shifts.SetLocal
        ( iLoc, 0, corner+C((x+0.5)*realStep,-(y+0.5)*imagStep) );
    }

    // Form the vector of invNorms
    DistMatrix<Real,VR,STAR> invNorms(g);
    auto itCounts = 
        QuasiTriangularPseudospectrum( U, w, shifts, invNorms, psCtrl );

    // Rearrange the vectors into grids
    DistMatrix<Int> itCountMap(g); 
    pspec::ReshapeIntoGrid( realSize, imagSize, invNorms, invNormMap );
    pspec::ReshapeIntoGrid( realSize, imagSize, itCounts, itCountMap );
    return itCountMap;
}

template<typename F>
inline DistMatrix<Int>
HessenbergPseudospectrum
( const DistMatrix<F>& H, DistMatrix<BASE(F)>& invNormMap, 
  Complex<BASE(F)> center, BASE(F) realWidth, BASE(F) imagWidth, 
  Int realSize, Int imagSize, 
  PseudospecCtrl<BASE(F)> psCtrl=PseudospecCtrl<BASE(F)>() )
{
    DEBUG_ONLY(CallStackEntry cse("HessenbergPseudospectrum"))
    typedef Base<F> Real;
    typedef Complex<Real> C;
    const Grid& g = H.Grid();

    psCtrl.snapCtrl.realSize = realSize;
    psCtrl.snapCtrl.imagSize = imagSize;

    const Real realStep = realWidth/realSize;
    const Real imagStep = imagWidth/imagSize;
    const C corner = center + C(-realWidth/2,imagWidth/2);
    DistMatrix<C,VR,STAR> shifts( realSize*imagSize, 1, g );
    const Int numLocShifts = shifts.LocalHeight();
    for( Int iLoc=0; iLoc<numLocShifts; ++iLoc )
    {
        const Int i = shifts.GlobalRow(iLoc);
        const Int x = i / imagSize;
        const Int y = i % imagSize;
        shifts.SetLocal
        ( iLoc, 0, corner+C((x+0.5)*realStep,-(y+0.5)*imagStep) );
    }

    // Form the vector of invNorms
    DistMatrix<Real,VR,STAR> invNorms(g);
    auto itCounts = HessenbergPseudospectrum( H, shifts, invNorms, psCtrl );

    // Rearrange the vectors into grids
    DistMatrix<Int> itCountMap(g); 
    pspec::ReshapeIntoGrid( realSize, imagSize, invNorms, invNormMap );
    pspec::ReshapeIntoGrid( realSize, imagSize, itCounts, itCountMap );
    return itCountMap;
}

template<typename F>
inline Matrix<Int>
Pseudospectrum
( const Matrix<F>& A, Matrix<BASE(F)>& invNormMap, 
  Complex<BASE(F)> center, BASE(F) realWidth, BASE(F) imagWidth,
  Int realSize, Int imagSize, 
  PseudospecCtrl<BASE(F)> psCtrl=PseudospecCtrl<BASE(F)>() )
{
    DEBUG_ONLY(CallStackEntry cse("Pseudospectrum"))
    typedef Base<F> Real;
    typedef Complex<Real> C;

    psCtrl.snapCtrl.realSize = realSize;
    psCtrl.snapCtrl.imagSize = imagSize;

    const Real realStep = realWidth/realSize;
    const Real imagStep = imagWidth/imagSize;
    const C corner = center + C(-realWidth/2,imagWidth/2);
    Matrix<C> shifts( realSize*imagSize, 1, A.Grid() );
    for( Int j=0; j<realSize*imagSize; ++j )
    {
        const Int x = j / imagSize;
        const Int y = j % imagSize;
        shifts.Set( j, 0, corner+C((x+0.5)*realStep,-(y+0.5)*imagStep) );
    }

    // Form the vector of invNorms
    Matrix<Real> invNorms;
    auto itCounts = Pseudospectrum( A, shifts, invNorms, psCtrl );

    // Rearrange the vectors into grids
    Matrix<Int> itCountMap; 
    pspec::ReshapeIntoGrid( realSize, imagSize, invNorms, invNormMap );
    pspec::ReshapeIntoGrid( realSize, imagSize, itCounts, itCountMap );
    return itCountMap;
}

template<typename F>
inline DistMatrix<Int>
Pseudospectrum
( const DistMatrix<F>& A, DistMatrix<BASE(F)>& invNormMap, 
  Complex<BASE(F)> center, BASE(F) realWidth, BASE(F) imagWidth, 
  Int realSize, Int imagSize, 
  PseudospecCtrl<BASE(F)> psCtrl=PseudospecCtrl<BASE(F)>() )
{
    DEBUG_ONLY(CallStackEntry cse("Pseudospectrum"))
    typedef Base<F> Real;
    typedef Complex<Real> C;
    const Grid& g = A.Grid();

    psCtrl.snapCtrl.realSize = realSize;
    psCtrl.snapCtrl.imagSize = imagSize;

    const Real realStep = realWidth/realSize;
    const Real imagStep = imagWidth/imagSize;
    const C corner = center + C(-realWidth/2,imagWidth/2);
    DistMatrix<C,VR,STAR> shifts( realSize*imagSize, 1, g );
    const Int numLocShifts = shifts.LocalHeight();
    for( Int iLoc=0; iLoc<numLocShifts; ++iLoc )
    {
        const Int i = shifts.GlobalRow(iLoc);
        const Int x = i / imagSize;
        const Int y = i % imagSize;
        shifts.SetLocal
        ( iLoc, 0, corner+C((x+0.5)*realStep,-(y+0.5)*imagStep) );
    }

    // Form the vector of invNorms
    DistMatrix<Real,VR,STAR> invNorms(g);
    auto itCounts = Pseudospectrum( A, shifts, invNorms, psCtrl );

    // Rearrange the vectors into grids
    DistMatrix<Int> itCountMap(g); 
    pspec::ReshapeIntoGrid( realSize, imagSize, invNorms, invNormMap );
    pspec::ReshapeIntoGrid( realSize, imagSize, itCounts, itCountMap );
    return itCountMap;
}

template<typename F>
inline Matrix<Int>
TriangularPseudospectrum
( const Matrix<F>& U, Matrix<BASE(F)>& invNormMap, 
  Complex<BASE(F)> center, Int realSize, Int imagSize, 
  PseudospecCtrl<BASE(F)> psCtrl=PseudospecCtrl<BASE(F)>() )
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
        if( psCtrl.progress )
            std::cout << "Setting width to 1 to handle zero matrix" 
                      << std::endl;
    }
    else if( radius >= 0.2*oneNorm )
    {
        width = 2.5*radius;
        if( psCtrl.progress )
            std::cout << "Setting width to " << width 
                      << " based on the spectral radius, " << radius 
                      << std::endl;
    }
    else
    {
        width = 0.8*oneNorm;
        if( psCtrl.progress )
            std::cout << "Setting width to " << width 
                      << " based on the one norm, " << oneNorm << std::endl;
    }

    return TriangularPseudospectrum
           ( U, invNormMap, center, width, width, realSize, imagSize, psCtrl );
}

template<typename Real>
inline Matrix<Int>
QuasiTriangularPseudospectrum
( const Matrix<Real>& U, const Matrix<Complex<Real> >& w,
  Matrix<Real>& invNormMap, 
  Complex<Real> center,
  Int realSize, Int imagSize, 
  PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() )
{
    DEBUG_ONLY(CallStackEntry cse("QuasiTriangularPseudospectrum"))

    const Real radius = MaxNorm( w );
    const Real oneNorm = OneNorm( U );

    // Essentially three cases are handled here:
    // 1) The zero matrix (force the pseudospectrum width to 1)
    // 2) Typical matrices (use a small multiple of the spectral radius)
    // 3) Highly non-normal matrices (e.g., triangular with zero main diagonal)
    Real width;
    if( oneNorm == Real(0) )
    {
        width = 1;
        if( psCtrl.progress )
            std::cout << "Setting width to 1 to handle zero matrix" 
                      << std::endl;
    }
    else if( radius >= 0.2*oneNorm )
    {
        width = 2.5*radius;
        if( psCtrl.progress )
            std::cout << "Setting width to " << width 
                      << " based on the spectral radius, " << radius 
                      << std::endl;
    }
    else
    {
        width = 0.8*oneNorm;
        if( psCtrl.progress )
            std::cout << "Setting width to " << width 
                      << " based on the one norm, " << oneNorm << std::endl;
    }

    return QuasiTriangularPseudospectrum
           ( U, w, invNormMap, center, width, width, realSize, imagSize, 
             psCtrl );
}

template<typename F>
inline Matrix<Int>
HessenbergPseudospectrum
( const Matrix<F>& H, Matrix<BASE(F)>& invNormMap, 
  Complex<BASE(F)> center, Int realSize, Int imagSize, 
  PseudospecCtrl<BASE(F)> psCtrl=PseudospecCtrl<BASE(F)>() )
{
    DEBUG_ONLY(CallStackEntry cse("HessenbergPseudospectrum"))

    const Base<F> infNorm = InfinityNorm( H );
    const Base<F> oneNorm = OneNorm( H );
    Base<F> width;
    if( oneNorm == Base<F>(0) )
    {
        width = 1;
        if( psCtrl.progress )
            std::cout << "Setting width to 1 to handle zero matrix" 
                      << std::endl;
    }
    else
    {
        width = 0.8*Max(oneNorm,infNorm);
        if( psCtrl.progress )
            std::cout << "Setting width to " << width 
                      << " based on the one norm, " << oneNorm 
                      << ", and infinity norm, " << infNorm << std::endl;
    }

    return HessenbergPseudospectrum
           ( H, invNormMap, center, width, width, realSize, imagSize, psCtrl );
}

template<typename F>
inline DistMatrix<Int>
TriangularPseudospectrum
( const DistMatrix<F>& U, DistMatrix<BASE(F)>& invNormMap, 
  Complex<BASE(F)> center, Int realSize, Int imagSize,
  PseudospecCtrl<BASE(F)> psCtrl=PseudospecCtrl<BASE(F)>() )
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
        if( psCtrl.progress && U.Grid().Rank() == 0 )
            std::cout << "Setting width to 1 to handle zero matrix"
                      << std::endl;
    }
    else if( radius >= 0.2*oneNorm )
    {
        width = 2.5*radius;
        if( psCtrl.progress && U.Grid().Rank() == 0 )
            std::cout << "Setting width to " << width 
                      << " based on the spectral radius, " << radius 
                      << std::endl;
    }
    else
    {
        width = 0.8*oneNorm;
        if( psCtrl.progress && U.Grid().Rank() == 0 )
            std::cout << "Setting width to " << width
                      << " based on the one norm, " << oneNorm << std::endl;
    }

    return TriangularPseudospectrum
           ( U, invNormMap, center, width, width, realSize, imagSize, psCtrl );
}

template<typename Real>
inline DistMatrix<Int>
QuasiTriangularPseudospectrum
( const DistMatrix<Real>& U, const DistMatrix<Complex<Real>,VR,STAR>& w,
  DistMatrix<Real>& invNormMap, 
  Complex<Real> center, Int realSize, Int imagSize,
  PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() )
{
    DEBUG_ONLY(CallStackEntry cse("QuasiTriangularPseudospectrum"))

    const Real radius = MaxNorm( w );
    const Real oneNorm = OneNorm( U );

    // Essentially three cases are handled here:
    // 1) The zero matrix (force the pseudospectrum width to 1)
    // 2) Typical matrices (use a small multiple of the spectral radius)
    // 3) Highly non-normal matrices (e.g., triangular with zero main diagonal)
    Real width;
    if( oneNorm == Real(0) )
    {
        width = 1;
        if( psCtrl.progress && U.Grid().Rank() == 0 )
            std::cout << "Setting width to 1 to handle zero matrix"
                      << std::endl;
    }
    else if( radius >= 0.2*oneNorm )
    {
        width = 2.5*radius;
        if( psCtrl.progress && U.Grid().Rank() == 0 )
            std::cout << "Setting width to " << width 
                      << " based on the spectral radius, " << radius 
                      << std::endl;
    }
    else
    {
        width = 0.8*oneNorm;
        if( psCtrl.progress && U.Grid().Rank() == 0 )
            std::cout << "Setting width to " << width
                      << " based on the one norm, " << oneNorm << std::endl;
    }

    return QuasiTriangularPseudospectrum
           ( U, w, invNormMap, center, width, width, realSize, imagSize, 
             psCtrl );
}

template<typename F>
inline DistMatrix<Int>
HessenbergPseudospectrum
( const DistMatrix<F>& H, DistMatrix<BASE(F)>& invNormMap, 
  Complex<BASE(F)> center, Int realSize, Int imagSize,
  PseudospecCtrl<BASE(F)> psCtrl=PseudospecCtrl<BASE(F)>() )
{
    DEBUG_ONLY(CallStackEntry cse("HessenbergPseudospectrum"))

    const Base<F> oneNorm = OneNorm( H );
    const Base<F> infNorm = InfinityNorm( H );
    Base<F> width;
    if( oneNorm == Base<F>(0) )
    {
        width = 1;
        if( psCtrl.progress && H.Grid().Rank() == 0 )
            std::cout << "Setting width to 1 to handle zero matrix"
                      << std::endl;
    }
    else
    {
        width = 0.8*Max(oneNorm,infNorm);
        if( psCtrl.progress && H.Grid().Rank() == 0 )
            std::cout << "Setting width to " << width 
                      << " based on the one norm, " << oneNorm 
                      << ", and infinity norm, " << infNorm << std::endl;
    }

    return HessenbergPseudospectrum
           ( H, invNormMap, center, width, width, realSize, imagSize, psCtrl );
}

template<typename Real>
inline Matrix<Int>
Pseudospectrum
( const Matrix<Complex<Real> >& A, Matrix<Real>& invNormMap, 
  Complex<Real> center, Int realSize, Int imagSize, 
  PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() )
{
    DEBUG_ONLY(CallStackEntry cse("Pseudospectrum"))
    typedef Complex<Real> C;

    Matrix<C> B( A );
    if( psCtrl.schur )
    {
        Matrix<C> w;
        const bool fullTriangle = true;
        schur::QR( B, w, fullTriangle );
        return TriangularPseudospectrum
               ( B, invNormMap, center, realSize, imagSize, psCtrl );
    }
    else
    {
        Hessenberg( UPPER, B );
        return HessenbergPseudospectrum
               ( B, invNormMap, center, realSize, imagSize, psCtrl );
    }
}

template<typename Real>
inline Matrix<Int>
Pseudospectrum
( const Matrix<Real>& A, Matrix<Real>& invNormMap, 
  Complex<Real> center, Int realSize, Int imagSize, 
  PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() )
{
    DEBUG_ONLY(CallStackEntry cse("Pseudospectrum"))
    typedef Complex<Real> C;

    if( psCtrl.forceComplexSchur )
    {
        Matrix<C> ACpx;
        Copy( A, ACpx );
        return Pseudospectrum
               ( ACpx, invNormMap, center, realSize, imagSize, psCtrl );
    }

    if( psCtrl.schur )
    {
        Matrix<Real> B( A );
        Matrix<C> w;
        const bool fullTriangle = true;
        schur::QR( B, w, fullTriangle );
        if( psCtrl.forceComplexPs )
        {
            Matrix<C> BCpx;
            schur::RealToComplex( B, BCpx );
            return TriangularPseudospectrum
                   ( BCpx, invNormMap, center, realSize, imagSize, psCtrl );
        }
        return QuasiTriangularPseudospectrum
               ( B, w, invNormMap, center, realSize, imagSize, psCtrl );
    }
    else
        LogicError("Real Hessenberg algorithm not yet supported");
}

template<typename Real>
inline DistMatrix<Int>
Pseudospectrum
( const DistMatrix<Complex<Real> >& A, DistMatrix<Real>& invNormMap, 
  Complex<Real> center, Int realSize, Int imagSize,
  PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() )
{
    DEBUG_ONLY(CallStackEntry cse("Pseudospectrum"))
    typedef Complex<Real> C;

    const Grid& g = A.Grid();
    DistMatrix<C> B( A );

    if( psCtrl.schur )
    {
        DistMatrix<C,VR,STAR> w(g);
        const bool fullTriangle = true;
#ifdef ELEM_HAVE_SCALAPACK
        schur::QR( B, w, fullTriangle );
#else
        // We don't actually need the Schur vectors, but SDC requires their
        // computation in order to form the full triangular factor
        DistMatrix<C> X(g);
        schur::SDC( B, w, X, fullTriangle psCtrl.sdcCtrl );
        X.Empty();
#endif
 
        return TriangularPseudospectrum
               ( B, invNormMap, center, realSize, imagSize, psCtrl );
    }
    else
    {
        Hessenberg( UPPER, B );
        return HessenbergPseudospectrum
               ( B, invNormMap, center, realSize, imagSize, psCtrl );
    }
}

template<typename Real>
inline DistMatrix<Int>
Pseudospectrum
( const DistMatrix<Real>& A, DistMatrix<Real>& invNormMap, 
  Complex<Real> center, Int realSize, Int imagSize,
  PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() )
{
    DEBUG_ONLY(CallStackEntry cse("Pseudospectrum"))
    typedef Complex<Real> C;
    const Grid& g = A.Grid();

    if( psCtrl.forceComplexSchur )
    {
        DistMatrix<C> ACpx(g);
        Copy( A, ACpx );
        return Pseudospectrum
               ( ACpx, invNormMap, center, realSize, imagSize, psCtrl );
    }

    if( psCtrl.schur )
    {
        DistMatrix<Real> B( A );
        DistMatrix<C,VR,STAR> w(g);
        const bool fullTriangle = true;
#ifdef ELEM_HAVE_SCALAPACK
        schur::QR( B, w, fullTriangle );
#else
        // We don't actually need the Schur vectors, but SDC requires their
        // computation in order to form the full triangular factor
        DistMatrix<Real> X(g);
        schur::SDC( B, w, X, fullTriangle, psCtrl.sdcCtrl );
        X.Empty();
#endif
        if( psCtrl.forceComplexPs ) 
        {
            DistMatrix<C> BCpx(g);
            schur::RealToComplex( B, BCpx );
            return TriangularPseudospectrum
                   ( BCpx, invNormMap, center, realSize, imagSize, psCtrl );
        }
        return QuasiTriangularPseudospectrum
               ( B, w, invNormMap, center, realSize, imagSize, psCtrl );
    }
    else
        LogicError("Real Hessenberg algorithm not yet supported");
}

} // namespace elem

#endif // ifndef ELEM_PSEUDOSPECTRUM_HPP
