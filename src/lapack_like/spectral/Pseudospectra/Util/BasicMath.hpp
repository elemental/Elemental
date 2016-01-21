/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_PSEUDOSPECTRA_UTIL_BASICMATH_HPP
#define EL_PSEUDOSPECTRA_UTIL_BASICMATH_HPP

namespace El {
namespace pspec {

template<typename F>
inline bool
TriangIsNormal( const Matrix<F>& U, Base<F> tol )
{
    const Base<F> diagFrob = FrobeniusNorm(GetDiagonal(U));
    const Base<F> upperFrob = FrobeniusNorm( U );
    const Base<F> offDiagFrob = Sqrt(upperFrob*upperFrob-diagFrob*diagFrob);
    return offDiagFrob <= tol*diagFrob;
}

template<typename F>
inline bool
TriangIsNormal( const ElementalMatrix<F>& UPre, Base<F> tol )
{
    DistMatrixReadProxy<F,F,MC,MR> UProx( UPre );
    auto& U = UProx.GetLocked();

    const Base<F> diagFrob = FrobeniusNorm(GetDiagonal(U));
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
QuasiTriangIsNormal( const ElementalMatrix<F>& U, Base<F> tol )
{
    const auto w = schur::QuasiTriangEig( U );
    const Base<F> eigFrob = FrobeniusNorm( w );
    const Base<F> upperFrob = FrobeniusNorm( U );
    const Base<F> strictlyUpperFrob = Sqrt(upperFrob*upperFrob-eigFrob*eigFrob);
    return strictlyUpperFrob <= tol*eigFrob;
}

template<typename F>
inline Base<F> NormCap()
{ return Base<F>(1)/limits::Epsilon<Base<F>>(); }

template<typename F>
inline bool HasNan( const Matrix<F>& H )
{
    DEBUG_ONLY(CSE cse("pspec::HasNan"))
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
( const Matrix<FComp>& components,
  const Matrix<F>& X, Matrix<F>& Y )
{
    DEBUG_ONLY(
      CSE cse("pspec::ColumnSubtractions");
      if( components.Width() != 1 )
          LogicError("components assumed to be a column vector");
    )
    const Int numShifts = Y.Width();
    const Int m = Y.Height();
    for( Int j=0; j<numShifts; ++j )
    {
        const F gamma = components.Get(j,0);
        blas::Axpy( m, -gamma, X.LockedBuffer(0,j), 1, Y.Buffer(0,j), 1 );
    }
}

template<typename Real>
inline void
ColumnSubtractions
( const Matrix<Complex<Real>>& components,
  const Matrix<Real>& XReal, const Matrix<Real>& XImag,
        Matrix<Real>& YReal,       Matrix<Real>& YImag )
{
    DEBUG_ONLY(
      CSE cse("pspec::ColumnSubtractions");
      if( components.Width() != 1 )
          LogicError("components assumed to be a column vector");
    )
    const Int numShifts = YReal.Width();
    const Int m = YReal.Height();
    for( Int j=0; j<numShifts; ++j )
    {
        const Complex<Real> gamma = components.Get(j,0);
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
( const Matrix<FComp>& components,
  const DistMatrix<F>& X, DistMatrix<F>& Y )
{
    DEBUG_ONLY(
      CSE cse("pspec::ColumnSubtractions");
      if( X.ColAlign() != Y.ColAlign() || X.RowAlign() != Y.RowAlign() )
          LogicError("X and Y should have been aligned");
    )
    ColumnSubtractions( components, X.LockedMatrix(), Y.Matrix() );
}

template<typename Real>
inline void
ColumnSubtractions
( const Matrix<Complex<Real>>& components,
  const DistMatrix<Real>& XReal, const DistMatrix<Real>& XImag,
        DistMatrix<Real>& YReal,       DistMatrix<Real>& YImag )
{
    DEBUG_ONLY(
      CSE cse("pspec::ColumnSubtractions");
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
InnerProducts
( const Matrix<F>& X, const Matrix<F>& Y, Matrix<Base<F>>& innerProds )
{
    DEBUG_ONLY(CSE cse("pspec::InnerProducts"))
    typedef Base<F> Real;
    const Int numShifts = X.Width();
    const Int m = X.Height();
    innerProds.Resize( numShifts, 1 );
    for( Int j=0; j<numShifts; ++j )
    {
        const Real alpha =
            RealPart(blas::Dot( m, X.LockedBuffer(0,j), 1,
                                   Y.LockedBuffer(0,j), 1 ));
        innerProds.Set( j, 0, alpha );
    }
}

template<typename Real>
inline void
InnerProducts
( const Matrix<Real>& XReal, const Matrix<Real>& XImag,
  const Matrix<Real>& YReal, const Matrix<Real>& YImag, 
        Matrix<Real>& innerProds )
{
    DEBUG_ONLY(CSE cse("pspec::InnerProducts"))
    const Int numShifts = XReal.Width();
    const Int m = XReal.Height();
    innerProds.Resize( numShifts, 1 );
    for( Int j=0; j<numShifts; ++j )
    {
        const Real alpha =
            blas::Dot( m, XReal.LockedBuffer(0,j), 1,
                          YReal.LockedBuffer(0,j), 1 );
        const Real beta = 
            blas::Dot( m, XImag.LockedBuffer(0,j), 1,
                          YImag.LockedBuffer(0,j), 1 );
        innerProds.Set( j, 0, alpha+beta );
    }
}

template<typename F>
inline void
InnerProducts
( const Matrix<F>& X, const Matrix<F>& Y, Matrix<F>& innerProds )
{
    DEBUG_ONLY(CSE cse("pspec::InnerProducts"))
    const Int numShifts = X.Width();
    const Int m = X.Height();
    innerProds.Resize( numShifts, 1 );
    for( Int j=0; j<numShifts; ++j )
    {
        const F alpha =
            blas::Dot( m, X.LockedBuffer(0,j), 1,
                          Y.LockedBuffer(0,j), 1 );
        innerProds.Set( j, 0, alpha );
    }
}

template<typename Real>
inline void
InnerProducts
( const Matrix<Real>& XReal, const Matrix<Real>& XImag,
  const Matrix<Real>& YReal, const Matrix<Real>& YImag, 
        Matrix<Complex<Real>>& innerProds )
{
    DEBUG_ONLY(CSE cse("pspec::InnerProducts"))
    const Int numShifts = XReal.Width();
    const Int m = XReal.Height();
    innerProds.Resize( numShifts, 1 );
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
        innerProds.Set( j, 0, Complex<Real>(alpha+beta,delta-gamma) );
    }
}

// TODO: Use the appropriate distribution for 'innerProds'
template<typename F>
inline void
InnerProducts
( const DistMatrix<F>& X, const DistMatrix<F>& Y, Matrix<Base<F>>& innerProds )
{
    DEBUG_ONLY(
      CSE cse("pspec::InnerProducts");
      if( X.ColAlign() != Y.ColAlign() || X.RowAlign() != Y.RowAlign() )
          LogicError("X and Y should have been aligned");
    )
    InnerProducts( X.LockedMatrix(), Y.LockedMatrix(), innerProds );
    const Int numLocShifts = X.LocalWidth();
    mpi::AllReduce( innerProds.Buffer(), numLocShifts, mpi::SUM, X.ColComm() );
}

template<typename Real>
inline void
InnerProducts
( const DistMatrix<Real>& XReal, const DistMatrix<Real>& XImag,
  const DistMatrix<Real>& YReal, const DistMatrix<Real>& YImag,
  Matrix<Real>& innerProds )
{
    DEBUG_ONLY(
      CSE cse("pspec::InnerProducts");
      if( XReal.ColAlign() != YReal.ColAlign() || 
          XReal.RowAlign() != YReal.RowAlign() )
          LogicError("X and Y should have been aligned");
    )
    InnerProducts
    ( XReal.LockedMatrix(), XImag.LockedMatrix(), 
      YReal.LockedMatrix(), YImag.LockedMatrix(), innerProds );
    const Int numLocShifts = XReal.LocalWidth();
    mpi::AllReduce
    ( innerProds.Buffer(), numLocShifts, mpi::SUM, XReal.ColComm() );
}

template<typename F>
inline void
InnerProducts
( const DistMatrix<F>& X, const DistMatrix<F>& Y, Matrix<F>& innerProds )
{
    DEBUG_ONLY(
      CSE cse("pspec::InnerProducts");
      if( X.ColAlign() != Y.ColAlign() || X.RowAlign() != Y.RowAlign() )
          LogicError("X and Y should have been aligned");
    )
    InnerProducts( X.LockedMatrix(), Y.LockedMatrix(), innerProds );
    const Int numLocShifts = X.LocalWidth();
    mpi::AllReduce( innerProds.Buffer(), numLocShifts, mpi::SUM, X.ColComm() );
}

template<typename Real>
inline void
InnerProducts
( const DistMatrix<Real>& XReal, const DistMatrix<Real>& XImag,
  const DistMatrix<Real>& YReal, const DistMatrix<Real>& YImag,
        Matrix<Complex<Real>>& innerProds )
{
    DEBUG_ONLY(
      CSE cse("pspec::InnerProducts");
      if( XReal.ColAlign() != YReal.ColAlign() || 
          XReal.RowAlign() != YReal.RowAlign() )
          LogicError("X and Y should have been aligned");
    )
    InnerProducts
    ( XReal.LockedMatrix(), XImag.LockedMatrix(), 
      YReal.LockedMatrix(), YImag.LockedMatrix(), innerProds );
    const Int numLocShifts = XReal.LocalWidth();
    mpi::AllReduce
    ( innerProds.Buffer(), numLocShifts, mpi::SUM, XReal.ColComm() );
}

template<typename F>
inline void
FixColumns( Matrix<F>& X )
{
    DEBUG_ONLY(CSE cse("pspec::FixColumns"))
    typedef Base<F> Real;
    Matrix<Real> norms;
    ColumnTwoNorms( X, norms );
    const Int n = X.Width();
    for( Int j=0; j<n; ++j )
    {
        auto x = X( ALL, IR(j) );
        Real norm = norms.Get(j,0);
        if( norm == Real(0) )
        {
            MakeGaussian( x );
            norm = FrobeniusNorm( x );
        }
        x *= 1/norm;
    }
}

template<typename F,Dist U,Dist V>
inline void
FixColumns( DistMatrix<F,U,V>& X )
{
    DEBUG_ONLY(CSE cse("pspec::FixColumns"))
    typedef Base<F> Real;
    DistMatrix<Real,V,STAR> norms( X.Grid() );
    ColumnTwoNorms( X, norms );
    const Int nLocal = X.LocalWidth();
    for( Int jLoc=0; jLoc<nLocal; ++jLoc )
    {
        const Int j = X.GlobalCol(jLoc);
        auto x = X( ALL, IR(j) );
        Real norm = norms.GetLocal(jLoc,0);
        if( norm == Real(0) )
        {
            MakeGaussian( x );
            norm = FrobeniusNorm( x );
        }
        x *= 1/norm;
    }
}

template<typename Real>
inline void CapEstimates( Matrix<Real>& activeEsts )
{
    DEBUG_ONLY(CSE cse("pspec::CapEstimates"))
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
    DEBUG_ONLY(CSE cse("pspec::CapEstimates"))
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
    DEBUG_ONLY(CSE cse("pspec::FindConverged"))
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
      CSE cse("pspec::FindConverged");
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
} // namespace El

#endif // ifndef EL_PSEUDOSPECTRA_UTIL_BASICMATH_HPP
