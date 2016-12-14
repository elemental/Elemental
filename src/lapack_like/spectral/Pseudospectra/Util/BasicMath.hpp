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

template<typename Field>
bool TriangIsNormal( const Matrix<Field>& U, Base<Field> tol )
{
    const Base<Field> diagFrob = FrobeniusNorm(GetDiagonal(U));
    const Base<Field> upperFrob = FrobeniusNorm( U );
    const Base<Field> offDiagFrob = Sqrt(upperFrob*upperFrob-diagFrob*diagFrob);
    return offDiagFrob <= tol*diagFrob;
}

template<typename Field>
bool TriangIsNormal( const AbstractDistMatrix<Field>& UPre, Base<Field> tol )
{
    DistMatrixReadProxy<Field,Field,MC,MR> UProx( UPre );
    auto& U = UProx.GetLocked();

    const Base<Field> diagFrob = FrobeniusNorm(GetDiagonal(U));
    const Base<Field> upperFrob = FrobeniusNorm( U );
    const Base<Field> offDiagFrob = Sqrt(upperFrob*upperFrob-diagFrob*diagFrob);
    return offDiagFrob <= tol*diagFrob;
}

template<typename Field>
bool QuasiTriangIsNormal( const Matrix<Field>& U, Base<Field> tol )
{
    const auto w = schur::QuasiTriangEig( U );
    const Base<Field> eigFrob = FrobeniusNorm( w );
    const Base<Field> upperFrob = FrobeniusNorm( U );
    const Base<Field> strictlyUpperFrob =
      Sqrt(upperFrob*upperFrob-eigFrob*eigFrob);
    return strictlyUpperFrob <= tol*eigFrob;
}

template<typename Field>
bool QuasiTriangIsNormal( const AbstractDistMatrix<Field>& U, Base<Field> tol )
{
    const auto w = schur::QuasiTriangEig( U );
    const Base<Field> eigFrob = FrobeniusNorm( w );
    const Base<Field> upperFrob = FrobeniusNorm( U );
    const Base<Field> strictlyUpperFrob =
      Sqrt(upperFrob*upperFrob-eigFrob*eigFrob);
    return strictlyUpperFrob <= tol*eigFrob;
}

template<typename Field>
Base<Field> NormCap() { return Base<Field>(1)/limits::Epsilon<Base<Field>>(); }

template<typename Field>
bool HasNan( const Matrix<Field>& H )
{
    EL_DEBUG_CSE
    bool hasNan = false;
    const Int m = H.Height();
    const Int n = H.Width();
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<m; ++i )
            if( !limits::IsFinite(RealPart(H(i,j))) ||
                !limits::IsFinite(ImagPart(H(i,j))) )
                hasNan = true;
    return hasNan;
}

template<typename Field,typename FieldComp>
void ColumnSubtractions
( const Matrix<FieldComp>& components,
  const Matrix<Field>& X, Matrix<Field>& Y )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( components.Width() != 1 )
          LogicError("components assumed to be a column vector");
    )
    const Int numShifts = Y.Width();
    const Int m = Y.Height();
    for( Int j=0; j<numShifts; ++j )
    {
        const Field gamma = components(j);
        blas::Axpy( m, -gamma, X.LockedBuffer(0,j), 1, Y.Buffer(0,j), 1 );
    }
}

template<typename Real>
void ColumnSubtractions
( const Matrix<Complex<Real>>& components,
  const Matrix<Real>& XReal, const Matrix<Real>& XImag,
        Matrix<Real>& YReal,       Matrix<Real>& YImag )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( components.Width() != 1 )
          LogicError("components assumed to be a column vector");
    )
    const Int numShifts = YReal.Width();
    const Int m = YReal.Height();
    for( Int j=0; j<numShifts; ++j )
    {
        const Complex<Real> gamma = components(j);
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

template<typename Field,typename FieldComp>
void ColumnSubtractions
( const Matrix<FieldComp>& components,
  const DistMatrix<Field>& X, DistMatrix<Field>& Y )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( X.ColAlign() != Y.ColAlign() || X.RowAlign() != Y.RowAlign() )
          LogicError("X and Y should have been aligned");
    )
    ColumnSubtractions( components, X.LockedMatrix(), Y.Matrix() );
}

template<typename Real>
void ColumnSubtractions
( const Matrix<Complex<Real>>& components,
  const DistMatrix<Real>& XReal, const DistMatrix<Real>& XImag,
        DistMatrix<Real>& YReal,       DistMatrix<Real>& YImag )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( XReal.ColAlign() != YReal.ColAlign() ||
          XReal.RowAlign() != YReal.RowAlign() )
          LogicError("X and Y should have been aligned");
    )
    ColumnSubtractions
    ( components, XReal.LockedMatrix(), XImag.LockedMatrix(),
                  YReal.Matrix(),       YImag.Matrix() );
}

template<typename Field>
void InnerProducts
( const Matrix<Field>& X,
  const Matrix<Field>& Y,
        Matrix<Base<Field>>& innerProds )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    const Int numShifts = X.Width();
    const Int m = X.Height();
    innerProds.Resize( numShifts, 1 );
    for( Int j=0; j<numShifts; ++j )
    {
        const Real alpha =
            RealPart(blas::Dot( m, X.LockedBuffer(0,j), 1,
                                   Y.LockedBuffer(0,j), 1 ));
        innerProds(j) = alpha;
    }
}

template<typename Real>
void InnerProducts
( const Matrix<Real>& XReal, const Matrix<Real>& XImag,
  const Matrix<Real>& YReal, const Matrix<Real>& YImag,
        Matrix<Real>& innerProds )
{
    EL_DEBUG_CSE
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
        innerProds(j) = alpha+beta;
    }
}

template<typename Field>
void InnerProducts
( const Matrix<Field>& X, const Matrix<Field>& Y, Matrix<Field>& innerProds )
{
    EL_DEBUG_CSE
    const Int numShifts = X.Width();
    const Int m = X.Height();
    innerProds.Resize( numShifts, 1 );
    for( Int j=0; j<numShifts; ++j )
    {
        const Field alpha =
            blas::Dot( m, X.LockedBuffer(0,j), 1,
                          Y.LockedBuffer(0,j), 1 );
        innerProds(j) = alpha;
    }
}

template<typename Real>
void InnerProducts
( const Matrix<Real>& XReal, const Matrix<Real>& XImag,
  const Matrix<Real>& YReal, const Matrix<Real>& YImag,
        Matrix<Complex<Real>>& innerProds )
{
    EL_DEBUG_CSE
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
        innerProds(j) = Complex<Real>(alpha+beta,delta-gamma);
    }
}

// TODO(poulson): Use the appropriate distribution for 'innerProds'
template<typename Field>
void InnerProducts
( const DistMatrix<Field>& X,
  const DistMatrix<Field>& Y,
        Matrix<Base<Field>>& innerProds )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( X.ColAlign() != Y.ColAlign() || X.RowAlign() != Y.RowAlign() )
          LogicError("X and Y should have been aligned");
    )
    InnerProducts( X.LockedMatrix(), Y.LockedMatrix(), innerProds );
    const Int numLocShifts = X.LocalWidth();
    mpi::AllReduce( innerProds.Buffer(), numLocShifts, mpi::SUM, X.ColComm() );
}

template<typename Real>
void InnerProducts
( const DistMatrix<Real>& XReal, const DistMatrix<Real>& XImag,
  const DistMatrix<Real>& YReal, const DistMatrix<Real>& YImag,
  Matrix<Real>& innerProds )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
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

template<typename Field>
void InnerProducts
( const DistMatrix<Field>& X,
  const DistMatrix<Field>& Y,
        Matrix<Field>& innerProds )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( X.ColAlign() != Y.ColAlign() || X.RowAlign() != Y.RowAlign() )
          LogicError("X and Y should have been aligned");
    )
    InnerProducts( X.LockedMatrix(), Y.LockedMatrix(), innerProds );
    const Int numLocShifts = X.LocalWidth();
    mpi::AllReduce( innerProds.Buffer(), numLocShifts, mpi::SUM, X.ColComm() );
}

template<typename Real>
void InnerProducts
( const DistMatrix<Real>& XReal, const DistMatrix<Real>& XImag,
  const DistMatrix<Real>& YReal, const DistMatrix<Real>& YImag,
        Matrix<Complex<Real>>& innerProds )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
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

template<typename Field>
void FixColumns( Matrix<Field>& X )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    Matrix<Real> norms;
    ColumnTwoNorms( X, norms );
    const Int n = X.Width();
    for( Int j=0; j<n; ++j )
    {
        auto x = X( ALL, IR(j) );
        Real norm = norms(j);
        if( norm == Real(0) )
        {
            MakeGaussian( x );
            norm = FrobeniusNorm( x );
        }
        x *= 1/norm;
    }
}

template<typename Field,Dist U,Dist V>
void FixColumns( DistMatrix<Field,U,V>& X )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    DistMatrix<Real,V,STAR> norms( X.Grid() );
    ColumnTwoNorms( X, norms );
    const Int nLocal = X.LocalWidth();
    auto& normsLoc = norms.Matrix();
    for( Int jLoc=0; jLoc<nLocal; ++jLoc )
    {
        const Int j = X.GlobalCol(jLoc);
        auto x = X( ALL, IR(j) );
        Real norm = normsLoc(jLoc);
        if( norm == Real(0) )
        {
            MakeGaussian( x );
            norm = FrobeniusNorm( x );
        }
        x *= 1/norm;
    }
}

template<typename Real>
void CapEstimates( Matrix<Real>& activeEsts )
{
    EL_DEBUG_CSE
    const Real normCap = NormCap<Real>();
    const Int n = activeEsts.Height();
    for( Int j=0; j<n; ++j )
    {
        Real alpha = activeEsts(j);
        if( !limits::IsFinite(alpha) || alpha >= normCap )
            alpha = normCap;
        activeEsts(j) = alpha;
    }
}

template<typename Real>
void CapEstimates( DistMatrix<Real,MR,STAR>& activeEsts )
{
    EL_DEBUG_CSE
    CapEstimates( activeEsts.Matrix() );
}

template<typename Real>
Matrix<Int>
FindConverged
( const Matrix<Real>& lastActiveEsts,
  const Matrix<Real>& activeEsts,
        Matrix<Int >& activeItCounts,
        Real maxDiff )
{
    EL_DEBUG_CSE
    const Real normCap = NormCap<Real>();

    const Int numActiveShifts=activeEsts.Height();
    Matrix<Int> activeConverged;
    Zeros( activeConverged, numActiveShifts, 1 );

    for( Int j=0; j<numActiveShifts; ++j )
    {
        const Real lastEst = lastActiveEsts(j);
        const Real currEst = activeEsts(j);
        bool converged = false;
        if( currEst >= normCap )
            converged = true;
        else if( Abs(currEst) > 0 )
            converged = (Abs(lastEst-currEst)/Abs(currEst) <= maxDiff);

        if( converged )
            activeConverged(j) = 1;
        else
            ++activeItCounts(j);
    }
    return activeConverged;
}

template<typename Real>
DistMatrix<Int,MR,STAR>
FindConverged
( const DistMatrix<Real,MR,STAR>& lastActiveEsts,
  const DistMatrix<Real,MR,STAR>& activeEsts,
        DistMatrix<Int, VR,STAR>& activeItCounts,
        Real maxDiff )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( activeItCounts.ColAlign()%activeEsts.ColStride() !=
          activeEsts.ColAlign() )
          LogicError("Invalid column alignment");
    )
    const Real normCap = NormCap<Real>();

    DistMatrix<Int,MR,STAR> activeConverged( activeEsts.Grid() );
    activeConverged.AlignWith( activeEsts );
    Zeros( activeConverged, activeEsts.Height(), 1 );

    auto& activeEstsLoc = activeEsts.LockedMatrix();
    auto& lastActiveEstsLoc = lastActiveEsts.LockedMatrix();
    auto& activeConvergedLoc = activeConverged.Matrix();

    const Int numLocShifts=activeEsts.LocalHeight();
    for( Int iLoc=0; iLoc<numLocShifts; ++iLoc )
    {
        const Real lastEst = lastActiveEstsLoc(iLoc);
        const Real currEst = activeEstsLoc(iLoc);
        bool converged = false;
        if( currEst >= normCap )
            converged = true;
        else if( Abs(currEst) > 0 )
            converged = (Abs(lastEst-currEst)/Abs(currEst) <= maxDiff);

        if( converged )
            activeConvergedLoc(iLoc) = 1;
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
