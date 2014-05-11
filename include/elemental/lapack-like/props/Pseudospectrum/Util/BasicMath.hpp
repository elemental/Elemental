/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_PSEUDOSPECTRUM_UTIL_BASICMATH_HPP
#define ELEM_PSEUDOSPECTRUM_UTIL_BASICMATH_HPP

#include ELEM_FROBENIUSNORM_INC

namespace elem {
namespace pspec {

template<typename F>
inline bool
TriangIsNormal( const Matrix<F>& U, Base<F> tol )
{
    auto w = U.GetDiagonal();
    const Base<F> diagFrob = FrobeniusNorm( w );
    const Base<F> upperFrob = FrobeniusNorm( U );
    const Base<F> offDiagFrob = Sqrt(upperFrob*upperFrob-diagFrob*diagFrob);
    return offDiagFrob <= tol*diagFrob;
}

template<typename F>
inline bool
TriangIsNormal( const DistMatrix<F>& U, Base<F> tol )
{
    auto w = U.GetDiagonal();
    const Base<F> diagFrob = FrobeniusNorm( w );
    const Base<F> upperFrob = FrobeniusNorm( U );
    const Base<F> offDiagFrob = Sqrt(upperFrob*upperFrob-diagFrob*diagFrob);
    return offDiagFrob <= tol*diagFrob;
}

template<typename F>
inline bool
QuasiTriangIsNormal( const Matrix<F>& U, Base<F> tol )
{
    const auto w = schur::QuasiTriangEig( U );
    const Base<F> eigFrob = FrobeniusNorm( w );
    const Base<F> upperFrob = FrobeniusNorm( U );
    const Base<F> strictlyUpperFrob = Sqrt(upperFrob*upperFrob-eigFrob*eigFrob);
    return strictlyUpperFrob <= tol*eigFrob;
}

template<typename F>
inline bool
QuasiTriangIsNormal( const DistMatrix<F>& U, Base<F> tol )
{
    const auto w = schur::QuasiTriangEig( U );
    const Base<F> eigFrob = FrobeniusNorm( w );
    const Base<F> upperFrob = FrobeniusNorm( U );
    const Base<F> strictlyUpperFrob = Sqrt(upperFrob*upperFrob-eigFrob*eigFrob);
    return strictlyUpperFrob <= tol*eigFrob;
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
( const std::vector<Complex<Real>>& components,
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
( const std::vector<Complex<Real>>& components,
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
ColumnNorms( const Matrix<F>& X, Matrix<Base<F>>& norms )
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
ColumnNorms( const DistMatrix<F,U,V>& X, DistMatrix<Base<F>,V,STAR>& norms )
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
ColumnNorms( const Matrix<F>& X, std::vector<Base<F>>& norms )
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
ColumnNorms( const DistMatrix<F>& X, std::vector<Base<F>>& norms )
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
( const Matrix<F>& X, const Matrix<F>& Y, std::vector<Base<F>>& innerProds )
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
        std::vector<Complex<Real>>& innerProds )
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
  std::vector<Base<F>>& innerProds )
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
        std::vector<Complex<Real>>& innerProds )
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
InvBetaScale( const std::vector<Base<F>>& scales, Matrix<F>& Y )
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
InvBetaScale( const std::vector<Base<F>>& scales, DistMatrix<F>& Y )
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

} // namespace pspec
} // namespace elem

#endif // ifndef ELEM_PSEUDOSPECTRUM_UTIL_BASICMATH_HPP
