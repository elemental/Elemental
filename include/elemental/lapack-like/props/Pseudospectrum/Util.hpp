/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_PSEUDOSPECTRUM_UTIL_HPP
#define ELEM_PSEUDOSPECTRUM_UTIL_HPP

#include ELEM_ENTRYWISEMAP_INC
#include ELEM_MULTISHIFTQUASITRSM_INC
#include ELEM_MULTISHIFTTRSM_INC
#include ELEM_MULTISHIFTHESSSOLVE_INC
#include ELEM_FROBENIUSNORM_INC
#include ELEM_ZERONORM_INC
#include ELEM_ONES_INC

namespace elem {

// Configurations for how often and what format numerical (num) and image (img)
// snapshots of the pseudospectral estimates should be saved
struct SnapshotCtrl
{
    Int realSize, imagSize;
    Int imgFreq, numFreq;
    Int imgSaveCount, numSaveCount;
    std::string imgBase, numBase;
    FileFormat imgFormat, numFormat;

    SnapshotCtrl()
    : realSize(0), imagSize(0),
      imgFreq(0), numFreq(0), imgSaveCount(0), numSaveCount(0),
      imgBase("ps"), numBase("ps"), imgFormat(PNG), numFormat(ASCII_MATLAB)
    { }
};

template<typename Real>
struct PseudospecCtrl
{
    // Preprocessing configuration
    bool schur; // begin with reduction to Schur form?
    bool forceComplexSchur;
    bool forceComplexPs;
    SdcCtrl<Real> sdcCtrl;

    // Convergence and deflation criteria
    Int maxIts;
    Real tol;
    bool deflate; 

    // (Implicitly Restarted) Arnoldi/Lanczos. If basisSize > 1, then
    // there is implicit restarting
    bool arnoldi;
    Int basisSize;
    bool reorthog; // only matters for IRL, which isn't currently used

    // Whether or not to print progress information at each iteration
    bool progress;

    SnapshotCtrl snapCtrl;

    PseudospecCtrl()
    : schur(true), forceComplexSchur(false), forceComplexPs(false), sdcCtrl(),
      maxIts(200), tol(1e-6), deflate(true),
      arnoldi(true), basisSize(10), reorthog(true),
      progress(false), snapCtrl()
    { }
};

namespace pspec {

template<typename F>
inline bool
NumericallyNormal( const Matrix<F>& U, BASE(F) tol )
{
    auto w = U.GetDiagonal();
    const Base<F> diagFrob = FrobeniusNorm( w );
    const Base<F> upperFrob = FrobeniusNorm( U );
    const Base<F> offDiagFrob = Sqrt(upperFrob*upperFrob-diagFrob*diagFrob);
    return offDiagFrob <= tol*diagFrob;
}

template<typename F>
inline bool
NumericallyNormal
( const Matrix<F>& U, const Matrix<Complex<BASE(F)> >& w, BASE(F) tol )
{
    const Base<F> eigFrob = FrobeniusNorm( w );
    const Base<F> upperFrob = FrobeniusNorm( U );
    const Base<F> strictlyUpperFrob = Sqrt(upperFrob*upperFrob-eigFrob*eigFrob);
    return strictlyUpperFrob <= tol*eigFrob;
}

template<typename F>
inline bool
NumericallyNormal( const DistMatrix<F>& U, BASE(F) tol )
{
    auto w = U.GetDiagonal();
    const Base<F> diagFrob = FrobeniusNorm( w );
    const Base<F> upperFrob = FrobeniusNorm( U );
    const Base<F> offDiagFrob = Sqrt(upperFrob*upperFrob-diagFrob*diagFrob);
    return offDiagFrob <= tol*diagFrob;
}

template<typename F>
inline bool
NumericallyNormal
( const DistMatrix<F>& U, const DistMatrix<Complex<BASE(F)>,VR,STAR>& w, 
  BASE(F) tol )
{
    const Base<F> eigFrob = FrobeniusNorm( w );
    const Base<F> upperFrob = FrobeniusNorm( U );
    const Base<F> strictlyUpperFrob = Sqrt(upperFrob*upperFrob-eigFrob*eigFrob);
    return strictlyUpperFrob <= tol*eigFrob;
}

template<typename T>
inline void
ReshapeIntoGrid( Int realSize, Int imagSize, const Matrix<T>& x, Matrix<T>& X )
{
#if 0    
    X.Resize( realSize, imagSize );
    for( Int j=0; j<realSize; ++j )
    {
        auto XSub = View( X, 0, j, imagSize, 1 );
        auto xSub = LockedView( x, j*imagSize, 0, imagSize, 1 );
        XSub = xSub;
    }
#else
    // The sequential case can be optimized much more heavily than in parallel
    X.Resize( realSize, imagSize, realSize );
    MemCopy( X.Buffer(), x.LockedBuffer(), realSize*imagSize );
#endif
}

template<typename T>
inline void
ReshapeIntoGrid
( Int realSize, Int imagSize, const DistMatrix<T,VR,STAR>& x, DistMatrix<T>& X )
{
    X.SetGrid( x.Grid() );
    X.Resize( realSize, imagSize );
    for( Int j=0; j<realSize; ++j )
    {
        auto XSub = View( X, 0, j, imagSize, 1 );
        auto xSub = LockedView( x, j*imagSize, 0, imagSize, 1 );
        XSub = xSub;
    }
}

template<typename T>
inline void
RestoreOrdering
( const Matrix<Int>& preimage, Matrix<T>& x )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::RestoreOrdering"))
    auto xCopy = x;
    const Int numShifts = preimage.Height();
    for( Int j=0; j<numShifts; ++j )
    {
        const Int dest = preimage.Get(j,0);
        x.Set( dest, 0, xCopy.Get(j,0) );
    }
}

template<typename T1,typename T2>
inline void
RestoreOrdering
( const Matrix<Int>& preimage, Matrix<T1>& x, Matrix<T2>& y )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::RestoreOrdering"))
    auto xCopy = x;
    auto yCopy = y;
    const Int numShifts = preimage.Height();
    for( Int j=0; j<numShifts; ++j )
    {
        const Int dest = preimage.Get(j,0);
        x.Set( dest, 0, xCopy.Get(j,0) );
        y.Set( dest, 0, yCopy.Get(j,0) );
    }
}

template<typename T>
inline void
RestoreOrdering
( const DistMatrix<Int,VR,STAR>& preimage,
        DistMatrix<T,  VR,STAR>& x )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::RestoreOrdering"))
    DistMatrix<Int,STAR,STAR> preimageCopy( preimage );
    DistMatrix<T,STAR,STAR> xCopy( x );
    const Int numShifts = preimage.Height();
    for( Int j=0; j<numShifts; ++j )
    {
        const Int dest = preimageCopy.Get(j,0);
        x.Set( dest, 0, xCopy.Get(j,0) );
    }
}

template<typename T1,typename T2>
inline void
RestoreOrdering
( const DistMatrix<Int,VR,STAR>& preimage,
        DistMatrix<T1, VR,STAR>& x,
        DistMatrix<T2, VR,STAR>& y )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::RestoreOrdering"))
    DistMatrix<Int,STAR,STAR> preimageCopy( preimage );
    DistMatrix<T1, STAR,STAR> xCopy( x );
    DistMatrix<T2, STAR,STAR> yCopy( y );
    const Int numShifts = preimage.Height();
    for( Int j=0; j<numShifts; ++j )
    {
        const Int dest = preimageCopy.Get(j,0);
        x.Set( dest, 0, xCopy.Get(j,0) );
        y.Set( dest, 0, yCopy.Get(j,0) );
    }
}

template<typename F>
inline Base<F> NormCap()
{ return Base<F>(1)/lapack::MachineEpsilon<Base<F>>(); }

template<typename Real>
inline bool HasNan( const std::vector<Real>& x )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::HasNan"))
    bool hasNan = false;
    const Int n = x.size();
    for( Int j=0; j<n; ++j )
        if( std::isnan( x[j] ) )
            hasNan = true;
    return hasNan;
}

template<typename F>
inline bool HasNan( const Matrix<F>& H )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::HasNan"))
    bool hasNan = false;
    const Int m = H.Height();
    const Int n = H.Width();
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<m; ++i )
            if( std::isnan(H.GetRealPart(i,j)) ||
                std::isnan(H.GetImagPart(i,j)) )
                hasNan = true;
    return hasNan;
}

template<typename T1,typename T2>
inline void
ExtractList
( const std::vector<std::vector<T1>>& vecList, std::vector<T2>& list, Int i )
{
    DEBUG_ONLY(
        CallStackEntry cse("pspec::ExtractList");
        if( vecList.size() != 0 && vecList[0].size() <= i )
            LogicError("Invalid index");
    )
    const Int numVecs = vecList.size();
    list.resize( numVecs );
    for( Int k=0; k<numVecs; ++k )
        list[k] = vecList[k][i];
}

template<typename T1,typename T2>
inline void
ExtractList
( const std::vector<Matrix<T1>>& matList, std::vector<T2>& list, Int i, Int j )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::ExtractList"))
    const Int numMats = matList.size();
    list.resize( numMats );
    for( Int k=0; k<numMats; ++k )
        list[k] = matList[k].Get( i, j );
}

template<typename T1,typename T2>
inline void
PlaceList
( std::vector<std::vector<T1>>& vecList, const std::vector<T2>& list, Int i )
{
    DEBUG_ONLY(
        CallStackEntry cse("pspec::PlaceList");
        if( vecList.size() != 0 && vecList[0].size() <= i )
            LogicError("Invalid index");
        if( vecList.size() != list.size() )
            LogicError("List sizes do not match");
    )
    const Int numVecs = vecList.size();
    for( Int k=0; k<numVecs; ++k )
        vecList[k][i] = list[k];
}

template<typename T1,typename T2>
inline void
PlaceList
( std::vector<Matrix<T1>>& matList, const std::vector<T2>& list, Int i, Int j )
{
    DEBUG_ONLY(
        CallStackEntry cse("pspec::PlaceList");
        if( matList.size() != list.size() )
            LogicError("List sizes do not match");
    )
    const Int numMats = matList.size();
    for( Int k=0; k<numMats; ++k )
        matList[k].Set( i, j, list[k] );
}

template<typename T1,typename T2>
inline void
UpdateList
( std::vector<Matrix<T1>>& matList, const std::vector<T2>& list, Int i, Int j )
{
    DEBUG_ONLY(
        CallStackEntry cse("pspec::UpdateList");
        if( matList.size() != list.size() )
            LogicError("List sizes do not match");
    )
    const Int numMats = matList.size();
    for( Int k=0; k<numMats; ++k )
        matList[k].Update( i, j, list[k] );
}

template<typename T1,typename T2>
inline void
PushBackList
( std::vector<std::vector<T1>>& vecList, const std::vector<T2>& list )
{
    DEBUG_ONLY(
        CallStackEntry cse("pspec::PushBackList"); 
        if( vecList.size() != list.size() )
            LogicError("List sizes do not match");
    )
    const Int numVecs = vecList.size();
    for( Int k=0; k<numVecs; ++k )
        vecList[k].push_back( list[k] );
}

template<typename F,typename FComp>
inline void
ColumnSubtractions
( const std::vector<FComp>& components,
  const Matrix<F>& X, Matrix<F>& Y )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::ColumnSubtractions"))
    const Int numShifts = Y.Width();
    if( numShifts == 0 )
        return;
    const Int m = Y.Height();
    for( Int j=0; j<numShifts; ++j )
    {
        const F gamma = components[j];
        blas::Axpy( m, -gamma, X.LockedBuffer(0,j), 1, Y.Buffer(0,j), 1 );
    }
}

template<typename Real>
inline void
ColumnSubtractions
( const std::vector<Complex<Real> >& components,
  const Matrix<Real>& XReal, const Matrix<Real>& XImag,
        Matrix<Real>& YReal,       Matrix<Real>& YImag )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::ColumnSubtractions"))
    const Int numShifts = YReal.Width();
    if( numShifts == 0 )
        return;
    const Int m = YReal.Height();
    for( Int j=0; j<numShifts; ++j )
    {
        const Complex<Real> gamma = components[j];
        blas::Axpy
        ( m, -gamma.real(), XReal.LockedBuffer(0,j), 1, YReal.Buffer(0,j), 1 );
        blas::Axpy
        ( m,  gamma.imag(), XImag.LockedBuffer(0,j), 1, YReal.Buffer(0,j), 1 );
        blas::Axpy
        ( m, -gamma.real(), XImag.LockedBuffer(0,j), 1, YImag.Buffer(0,j), 1 );
        blas::Axpy
        ( m, -gamma.imag(), XReal.LockedBuffer(0,j), 1, YImag.Buffer(0,j), 1 );
    }
}

template<typename F,typename FComp>
inline void
ColumnSubtractions
( const std::vector<FComp>& components,
  const DistMatrix<F>& X, DistMatrix<F>& Y )
{
    DEBUG_ONLY(
        CallStackEntry cse("pspec::ColumnSubtractions");
        if( X.ColAlign() != Y.ColAlign() || X.RowAlign() != Y.RowAlign() )
            LogicError("X and Y should have been aligned");
    )
    ColumnSubtractions( components, X.LockedMatrix(), Y.Matrix() );
}

template<typename Real>
inline void
ColumnSubtractions
( const std::vector<Complex<Real> >& components,
  const DistMatrix<Real>& XReal, const DistMatrix<Real>& XImag,
        DistMatrix<Real>& YReal,       DistMatrix<Real>& YImag )
{
    DEBUG_ONLY(
        CallStackEntry cse("pspec::ColumnSubtractions");
        if( XReal.ColAlign() != YReal.ColAlign() || 
            XReal.RowAlign() != YReal.RowAlign() )
            LogicError("X and Y should have been aligned");
    )
    ColumnSubtractions
    ( components, XReal.LockedMatrix(), XImag.LockedMatrix(), 
                  YReal.Matrix(),       YImag.Matrix() );
}

template<typename F>
inline void
ColumnNorms( const Matrix<F>& X, Matrix<BASE(F)>& norms )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::ColumnNorms"))
    typedef Base<F> Real;
    const Int m = X.Height();
    const Int n = X.Width();
    norms.Resize( n, 1 );
    for( Int j=0; j<n; ++j )
    {
        const Real alpha = blas::Nrm2( m, X.LockedBuffer(0,j), 1 );
        norms.Set( j, 0, alpha );
    }
}

template<typename Real>
inline void
ColumnNorms
( const Matrix<Real>& XReal, const Matrix<Real>& XImag, Matrix<Real>& norms )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::ColumnNorms"))
    const Int m = XReal.Height();
    const Int n = XReal.Width();
    norms.Resize( n, 1 );
    for( Int j=0; j<n; ++j )
    {
        const Real alpha = blas::Nrm2( m, XReal.LockedBuffer(0,j), 1 );
        const Real beta  = blas::Nrm2( m, XImag.LockedBuffer(0,j), 1 );
        norms.Set( j, 0, lapack::SafeNorm(alpha,beta) );
    }
}

template<typename F,Dist U,Dist V>
inline void
ColumnNorms( const DistMatrix<F,U,V>& X, DistMatrix<BASE(F),V,STAR>& norms )
{
    DEBUG_ONLY(
        CallStackEntry cse("pspec::ColumnNorms");
        if( X.RowAlign() != norms.ColAlign() )
            LogicError("Invalid norms alignment");
    )
    typedef Base<F> Real;
    const Int n = X.Width();
    const Int mLocal = X.LocalHeight();
    const Int nLocal = X.LocalWidth();

    // TODO: Switch to more stable parallel norm computation using scaling
    norms.Resize( n, 1 );
    for( Int jLoc=0; jLoc<nLocal; ++jLoc )
    {
        const Base<F> localNorm = blas::Nrm2(mLocal,X.LockedBuffer(0,jLoc),1);
        norms.SetLocal( jLoc, 0, localNorm*localNorm );
    }

    mpi::AllReduce( norms.Buffer(), nLocal, mpi::SUM, X.ColComm() );
    for( Int jLoc=0; jLoc<nLocal; ++jLoc )
    {
        const Real alpha = norms.GetLocal(jLoc,0);
        norms.SetLocal( jLoc, 0, Sqrt(alpha) );
    }
}

template<typename Real,Dist U,Dist V>
inline void
ColumnNorms
( const DistMatrix<Real,U,V>& XReal, 
  const DistMatrix<Real,U,V>& XImag, DistMatrix<Real,V,STAR>& norms )
{
    DEBUG_ONLY(
        CallStackEntry cse("pspec::ColumnNorms");
        if( XReal.RowAlign() != norms.ColAlign() )
            LogicError("Invalid norms alignment");
    )
    const Int n = XReal.Width();
    const Int mLocal = XReal.LocalHeight();
    const Int nLocal = XReal.LocalWidth();

    // TODO: Switch to more stable parallel norm computation using scaling
    norms.Resize( n, 1 );
    for( Int jLoc=0; jLoc<nLocal; ++jLoc )
    {
        const Real alpha = blas::Nrm2(mLocal,XReal.LockedBuffer(0,jLoc),1);
        const Real beta = blas::Nrm2(mLocal,XImag.LockedBuffer(0,jLoc),1);
        const Real gamma = lapack::SafeNorm(alpha,beta);
        norms.SetLocal( jLoc, 0, gamma*gamma );
    }

    mpi::AllReduce( norms.Buffer(), nLocal, mpi::SUM, XReal.ColComm() );
    for( Int jLoc=0; jLoc<nLocal; ++jLoc )
    {
        const Real alpha = norms.GetLocal(jLoc,0);
        norms.SetLocal( jLoc, 0, Sqrt(alpha) );
    }
}

template<typename F>
inline void
ColumnNorms( const Matrix<F>& X, std::vector<BASE(F)>& norms )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::ColumnNorms"))
    typedef Base<F> Real;
    Matrix<Real> normCol;
    ColumnNorms( X, normCol );

    const Int numShifts = X.Width();
    norms.resize( numShifts );
    for( Int j=0; j<numShifts; ++j )
        norms[j] = normCol.Get(j,0);
}

template<typename Real>
inline void
ColumnNorms
( const Matrix<Real>& XReal, 
  const Matrix<Real>& XImag, std::vector<Real>& norms )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::ColumnNorms"))
    Matrix<Real> normCol;
    ColumnNorms( XReal, XImag, normCol );

    const Int numShifts = XReal.Width();
    norms.resize( numShifts );
    for( Int j=0; j<numShifts; ++j )
        norms[j] = normCol.Get(j,0);
}

template<typename F>
inline void
ColumnNorms( const DistMatrix<F>& X, std::vector<BASE(F)>& norms )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::ColumnNorms"))
    typedef Base<F> Real;
    DistMatrix<Real,MR,STAR> normCol( X.Grid() );
    ColumnNorms( X, normCol );

    const Int numLocShifts = X.LocalWidth();
    norms.resize( numLocShifts );
    for( Int jLoc=0; jLoc<numLocShifts; ++jLoc )
        norms[jLoc] = normCol.GetLocal(jLoc,0);
}

template<typename Real>
inline void
ColumnNorms
( const DistMatrix<Real>& XReal, 
  const DistMatrix<Real>& XImag, std::vector<Real>& norms )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::ColumnNorms"))
    DistMatrix<Real,MR,STAR> normCol( XReal.Grid() );
    ColumnNorms( XReal, XImag, normCol );

    const Int numLocShifts = XReal.LocalWidth();
    norms.resize( numLocShifts );
    for( Int jLoc=0; jLoc<numLocShifts; ++jLoc )
        norms[jLoc] = normCol.GetLocal(jLoc,0);
}

template<typename F>
inline void
InnerProducts
( const Matrix<F>& X, const Matrix<F>& Y, std::vector<BASE(F)>& innerProds )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::InnerProducts"))
    typedef Base<F> Real;
    const Int numShifts = X.Width();
    const Int m = X.Height();
    innerProds.resize( numShifts );
    for( Int j=0; j<numShifts; ++j )
    {
        const Real alpha =
            RealPart(blas::Dot( m, X.LockedBuffer(0,j), 1,
                                   Y.LockedBuffer(0,j), 1 ));
        innerProds[j] = alpha;
    }
}

template<typename Real>
inline void
InnerProducts
( const Matrix<Real>& XReal, const Matrix<Real>& XImag,
  const Matrix<Real>& YReal, const Matrix<Real>& YImag, 
        std::vector<Real>& innerProds )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::InnerProducts"))
    const Int numShifts = XReal.Width();
    const Int m = XReal.Height();
    innerProds.resize( numShifts );
    for( Int j=0; j<numShifts; ++j )
    {
        const Real alpha =
            blas::Dot( m, XReal.LockedBuffer(0,j), 1,
                          YReal.LockedBuffer(0,j), 1 );
        const Real beta = 
            blas::Dot( m, XImag.LockedBuffer(0,j), 1,
                          YImag.LockedBuffer(0,j), 1 );
        innerProds[j] = alpha + beta;
    }
}

template<typename F>
inline void
InnerProducts
( const Matrix<F>& X, const Matrix<F>& Y, std::vector<F>& innerProds )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::InnerProducts"))
    const Int numShifts = X.Width();
    const Int m = X.Height();
    innerProds.resize( numShifts );
    for( Int j=0; j<numShifts; ++j )
    {
        const F alpha =
            blas::Dot( m, X.LockedBuffer(0,j), 1,
                          Y.LockedBuffer(0,j), 1 );
        innerProds[j] = alpha;
    }
}

template<typename Real>
inline void
InnerProducts
( const Matrix<Real>& XReal, const Matrix<Real>& XImag,
  const Matrix<Real>& YReal, const Matrix<Real>& YImag, 
        std::vector<Complex<Real> >& innerProds )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::InnerProducts"))
    const Int numShifts = XReal.Width();
    const Int m = XReal.Height();
    innerProds.resize( numShifts );
    for( Int j=0; j<numShifts; ++j )
    {
        const Real alpha =
            blas::Dot( m, XReal.LockedBuffer(0,j), 1,
                          YReal.LockedBuffer(0,j), 1 );
        const Real beta = 
            blas::Dot( m, XImag.LockedBuffer(0,j), 1,
                          YImag.LockedBuffer(0,j), 1 );
        const Real delta = 
            blas::Dot( m, XReal.LockedBuffer(0,j), 1,
                          YImag.LockedBuffer(0,j), 1 );
        const Real gamma =
            blas::Dot( m, XImag.LockedBuffer(0,j), 1,
                          YReal.LockedBuffer(0,j), 1 );
        // Keep in mind that XImag should be conjugated
        innerProds[j] = Complex<Real>(alpha+beta,delta-gamma);
    }
}

template<typename F>
inline void
InnerProducts
( const DistMatrix<F>& X, const DistMatrix<F>& Y, 
  std::vector<BASE(F)>& innerProds )
{
    DEBUG_ONLY(
        CallStackEntry cse("pspec::InnerProducts");
        if( X.ColAlign() != Y.ColAlign() || X.RowAlign() != Y.RowAlign() )
            LogicError("X and Y should have been aligned");
    )
    InnerProducts( X.LockedMatrix(), Y.LockedMatrix(), innerProds );
    const Int numLocShifts = X.LocalWidth();
    mpi::AllReduce( innerProds.data(), numLocShifts, mpi::SUM, X.ColComm() );
}

template<typename Real>
inline void
InnerProducts
( const DistMatrix<Real>& XReal, const DistMatrix<Real>& XImag,
  const DistMatrix<Real>& YReal, const DistMatrix<Real>& YImag,
  std::vector<Real>& innerProds )
{
    DEBUG_ONLY(
        CallStackEntry cse("pspec::InnerProducts");
        if( XReal.ColAlign() != YReal.ColAlign() || 
            XReal.RowAlign() != YReal.RowAlign() )
            LogicError("X and Y should have been aligned");
    )
    InnerProducts
    ( XReal.LockedMatrix(), XImag.LockedMatrix(), 
      YReal.LockedMatrix(), YImag.LockedMatrix(), innerProds );
    const Int numLocShifts = XReal.LocalWidth();
    mpi::AllReduce
    ( innerProds.data(), numLocShifts, mpi::SUM, XReal.ColComm() );
}

template<typename F>
inline void
InnerProducts
( const DistMatrix<F>& X, const DistMatrix<F>& Y, std::vector<F>& innerProds )
{
    DEBUG_ONLY(
        CallStackEntry cse("pspec::InnerProducts");
        if( X.ColAlign() != Y.ColAlign() || X.RowAlign() != Y.RowAlign() )
            LogicError("X and Y should have been aligned");
    )
    InnerProducts( X.LockedMatrix(), Y.LockedMatrix(), innerProds );
    const Int numLocShifts = X.LocalWidth();
    mpi::AllReduce( innerProds.data(), numLocShifts, mpi::SUM, X.ColComm() );
}

template<typename Real>
inline void
InnerProducts
( const DistMatrix<Real>& XReal, const DistMatrix<Real>& XImag,
  const DistMatrix<Real>& YReal, const DistMatrix<Real>& YImag,
        std::vector<Complex<Real> >& innerProds )
{
    DEBUG_ONLY(
        CallStackEntry cse("pspec::InnerProducts");
        if( XReal.ColAlign() != YReal.ColAlign() || 
            XReal.RowAlign() != YReal.RowAlign() )
            LogicError("X and Y should have been aligned");
    )
    InnerProducts
    ( XReal.LockedMatrix(), XImag.LockedMatrix(), 
      YReal.LockedMatrix(), YImag.LockedMatrix(), innerProds );
    const Int numLocShifts = XReal.LocalWidth();
    mpi::AllReduce
    ( innerProds.data(), numLocShifts, mpi::SUM, XReal.ColComm() );
}

template<typename F>
inline void
InvBetaScale( const std::vector<BASE(F)>& scales, Matrix<F>& Y )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::InvBetaScale"))
    const Int numShifts = Y.Width();
    if( numShifts == 0 )
        return;
    const Int m = Y.Height();
    for( Int j=0; j<numShifts; ++j )
        blas::Scal( m, F(1)/scales[j], Y.Buffer(0,j), 1 );
}

template<typename F>
inline void
InvBetaScale( const std::vector<BASE(F)>& scales, DistMatrix<F>& Y )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::InvBetaScale"))
    InvBetaScale( scales, Y.Matrix() );
}

template<typename F>
inline void
FixColumns( Matrix<F>& X )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::FixColumns"))
    typedef Base<F> Real;
    Matrix<Real> norms;
    ColumnNorms( X, norms );
    const Int m = X.Height();
    const Int n = X.Width();
    for( Int j=0; j<n; ++j )
    {
        auto x = View( X, 0, j, m, 1 );
        Real norm = norms.Get(j,0);
        if( norm == Real(0) )
        {
            MakeGaussian( x );
            norm = FrobeniusNorm( x );
        }
        Scale( Real(1)/norm, x );
    }
}

template<typename F,Dist U,Dist V>
inline void
FixColumns( DistMatrix<F,U,V>& X )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::FixColumns"))
    typedef Base<F> Real;
    DistMatrix<Real,V,STAR> norms( X.Grid() );
    ColumnNorms( X, norms );
    const Int m = X.Height();
    const Int nLocal = X.LocalWidth();
    for( Int jLoc=0; jLoc<nLocal; ++jLoc )
    {
        const Int j = X.GlobalCol(jLoc);
        auto x = View( X, 0, j, m, 1 );
        Real norm = norms.GetLocal(jLoc,0);
        if( norm == Real(0) )
        {
            MakeGaussian( x );
            norm = FrobeniusNorm( x );
        }
        Scale( Real(1)/norm, x );
    }
}

template<typename Real>
inline void CapEstimates( Matrix<Real>& activeEsts )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::CapEstimates"))
    const Real normCap = NormCap<Real>();
    const Int n = activeEsts.Height();
    for( Int j=0; j<n; ++j )
    {
        Real alpha = activeEsts.Get(j,0);
        if( std::isnan(alpha) || alpha >= normCap )
            alpha = normCap;
        activeEsts.Set( j, 0, alpha );
    }
}

template<typename Real>
inline void CapEstimates( DistMatrix<Real,MR,STAR>& activeEsts )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::CapEstimates"))
    CapEstimates( activeEsts.Matrix() );
}

template<typename Real>
inline Matrix<Int>
FindConverged
( const Matrix<Real>& lastActiveEsts,
  const Matrix<Real>& activeEsts,
        Matrix<Int >& activeItCounts,
        Real maxDiff )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::FindConverged"))
    const Real normCap = NormCap<Real>();

    const Int numActiveShifts=activeEsts.Height();
    Matrix<Int> activeConverged;
    Zeros( activeConverged, numActiveShifts, 1 );

    for( Int j=0; j<numActiveShifts; ++j )
    {
        const Real lastEst = lastActiveEsts.Get(j,0);
        const Real currEst = activeEsts.Get(j,0);
        bool converged = false;
        if( currEst >= normCap )
            converged = true;
        else if( Abs(currEst) > 0 )
            converged = (Abs(lastEst-currEst)/Abs(currEst) <= maxDiff);

        if( converged )
            activeConverged.Set( j, 0, 1 );
        else
            activeItCounts.Update( j, 0, 1 );
    }
    return activeConverged;
}

template<typename Real>
inline DistMatrix<Int,MR,STAR>
FindConverged
( const DistMatrix<Real,MR,STAR>& lastActiveEsts,
  const DistMatrix<Real,MR,STAR>& activeEsts,
        DistMatrix<Int, VR,STAR>& activeItCounts,
        Real maxDiff )
{
    DEBUG_ONLY(
        CallStackEntry cse("pspec::FindConverged");
        if( activeItCounts.ColAlign()%activeEsts.ColStride() !=
            activeEsts.ColAlign() )
            LogicError("Invalid column alignment");
    )
    const Real normCap = NormCap<Real>();

    DistMatrix<Int,MR,STAR> activeConverged( activeEsts.Grid() );
    activeConverged.AlignWith( activeEsts );
    Zeros( activeConverged, activeEsts.Height(), 1 );

    const Int numLocShifts=activeEsts.LocalHeight();
    for( Int iLoc=0; iLoc<numLocShifts; ++iLoc )
    {
        const Real lastEst = lastActiveEsts.GetLocal(iLoc,0);
        const Real currEst = activeEsts.GetLocal(iLoc,0);
        bool converged = false;
        if( currEst >= normCap )
            converged = true;
        else if( Abs(currEst) > 0 )
            converged = (Abs(lastEst-currEst)/Abs(currEst) <= maxDiff);

        if( converged )
            activeConverged.SetLocal( iLoc, 0, 1 );
        else
        {
            const Int i = activeEsts.GlobalRow(iLoc);
            activeItCounts.Update( i, 0, 1 );
        }
    }

    return activeConverged;
}

template<typename Real>
inline void
Snapshot
( const Matrix<Real>& estimates, const Matrix<Int>& preimage, 
  Int numIts, bool deflate, SnapshotCtrl& snapCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::Snapshot"));
    if( snapCtrl.realSize != 0 && snapCtrl.imagSize != 0 )
    {
        const bool numSave = 
            ( snapCtrl.numFreq > 0 && 
              snapCtrl.numSaveCount >= snapCtrl.numFreq );
        const bool imgSave = 
            ( snapCtrl.imgFreq > 0 && 
              snapCtrl.imgSaveCount >= snapCtrl.imgFreq );
        Matrix<Real> invNorms;
        Matrix<Real> estMap; 
        if( numSave || imgSave )
        {
            invNorms = estimates;
            if( deflate )
                RestoreOrdering( preimage, invNorms );
            pspec::ReshapeIntoGrid
            ( snapCtrl.realSize, snapCtrl.imagSize, invNorms, estMap );
        }
        if( numSave )
        {
            std::ostringstream os;
            os << snapCtrl.numBase << "-" << numIts;
            Write( estMap, os.str(), snapCtrl.numFormat );
            snapCtrl.numSaveCount = 0;
        }
        if( imgSave )
        {
            EntrywiseMap( estMap, []( Real alpha ) { return Log(alpha); } );
            std::ostringstream os;
            os << snapCtrl.imgBase << "-" << numIts;
            Write( estMap, os.str(), snapCtrl.imgFormat );
            auto colorMap = GetColorMap();
            SetColorMap( GRAYSCALE_DISCRETE );
            Write( estMap, os.str()+"-discrete", snapCtrl.imgFormat );
            SetColorMap( colorMap );
            snapCtrl.imgSaveCount = 0;
        }
    }
}

template<typename Real>
inline void
Snapshot
( const DistMatrix<Real,MR,STAR>& estimates, 
  const DistMatrix<Int,    VR,STAR>& preimage, 
  Int numIts, bool deflate, SnapshotCtrl& snapCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::Snapshot"));
    if( snapCtrl.realSize != 0 && snapCtrl.imagSize != 0 )
    {
        const bool numSave = 
            ( snapCtrl.numFreq > 0 && 
              snapCtrl.numSaveCount >= snapCtrl.numFreq );
        const bool imgSave = 
            ( snapCtrl.imgFreq > 0 && 
              snapCtrl.imgSaveCount >= snapCtrl.imgFreq );
        DistMatrix<Real,VR,STAR> invNorms(estimates.Grid());
        DistMatrix<Real> estMap(estimates.Grid()); 
        if( numSave || imgSave )
        {
            invNorms = estimates;
            if( deflate )
                RestoreOrdering( preimage, invNorms );
            pspec::ReshapeIntoGrid
            ( snapCtrl.realSize, snapCtrl.imagSize, invNorms, estMap );
        }
        if( numSave )
        {
            std::ostringstream os;
            os << snapCtrl.numBase << "-" << numIts;
            Write( estMap, os.str(), snapCtrl.numFormat );
            snapCtrl.numSaveCount = 0;
        }
        if( imgSave )
        {
            EntrywiseMap( estMap, []( Real alpha ) { return Log(alpha); } );
            std::ostringstream os;
            os << snapCtrl.imgBase << "-" << numIts;
            Write( estMap, os.str(), snapCtrl.imgFormat );
            auto colorMap = GetColorMap();
            SetColorMap( GRAYSCALE_DISCRETE );
            Write( estMap, os.str()+"-discrete", snapCtrl.imgFormat );
            SetColorMap( colorMap );
            snapCtrl.imgSaveCount = 0;
        }
    }
}

} // namespace pspec
} // namespace elem

#endif // ifndef ELEM_PSEUDOSPECTRUM_UTIL_HPP
