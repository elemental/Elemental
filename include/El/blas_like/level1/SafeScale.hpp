/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_BLAS_SAFE_SCALE_HPP
#define EL_BLAS_SAFE_SCALE_HPP

namespace El {

template<typename Real,
         typename=EnableIf<IsReal<Real>>>
bool SafeScaleStep
(       Real& numerator,
        Real& denominator,
        Real& scaleStep,
  const Real& zero,
  const Real& smallNum,
  const Real& bigNum )
{
    EL_DEBUG_CSE

    Real shrunkDenominator = denominator * smallNum;
    if( Abs(shrunkDenominator) > Abs(numerator) && numerator != zero )
    {
        scaleStep = smallNum;
        denominator = shrunkDenominator;
        return false;
    }

    Real shrunkNumerator = numerator / bigNum;
    if( Abs(shrunkNumerator) > Abs(denominator) )
    {
        scaleStep = bigNum;
        numerator = shrunkNumerator;
        return false;
    }

    scaleStep = numerator / denominator;
    return true;
}

template<typename Field>
void SafeScale( Base<Field> numerator, Base<Field> denominator, Field& alpha )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    const Real zero(0);
    const Real smallNum = limits::SafeMin<Real>();
    const Real bigNum = Real(1) / smallNum;

    bool done = false;
    Real scaleStep;
    while( !done )
    {
        done =
          SafeScaleStep
          ( numerator, denominator, scaleStep, zero, smallNum, bigNum );
        alpha *= scaleStep;
    }
}

template<typename Field>
void SafeScale
( Base<Field> numerator, Base<Field> denominator, Matrix<Field>& A )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    const Real zero(0);
    const Real smallNum = limits::SafeMin<Real>();
    const Real bigNum = Real(1) / smallNum;

    bool done = false;
    Real scaleStep;
    while( !done )
    {
        done =
          SafeScaleStep
          ( numerator, denominator, scaleStep, zero, smallNum, bigNum );
        A *= scaleStep;
    }
}

template<typename Field>
void SafeScaleTrapezoid
( Base<Field> numerator, Base<Field> denominator,
  UpperOrLower uplo, Matrix<Field>& A, Int offset )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    const Real zero(0);
    const Real smallNum = limits::SafeMin<Real>();
    const Real bigNum = Real(1) / smallNum;

    bool done = false;
    Real scaleStep;
    while( !done )
    {
        done =
          SafeScaleStep
          ( numerator, denominator, scaleStep, zero, smallNum, bigNum );
        ScaleTrapezoid( scaleStep, uplo, A, offset );
    }
}

template<typename Field>
void SafeScale
( Base<Field> numerator, Base<Field> denominator, AbstractDistMatrix<Field>& A )
{
    EL_DEBUG_CSE
    SafeScale( numerator, denominator, A.Matrix() );
}

template<typename Field>
void SafeScaleTrapezoid
( Base<Field> numerator, Base<Field> denominator,
  UpperOrLower uplo, AbstractDistMatrix<Field>& A, Int offset )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    const Real zero(0);
    const Real smallNum = limits::SafeMin<Real>();
    const Real bigNum = Real(1) / smallNum;

    bool done = false;
    Real scaleStep;
    while( !done )
    {
        done =
          SafeScaleStep
          ( numerator, denominator, scaleStep, zero, smallNum, bigNum );
        ScaleTrapezoid( scaleStep, uplo, A, offset );
    }
}

template<typename Field>
void SafeScale
( Base<Field> numerator, Base<Field> denominator, SparseMatrix<Field>& A )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    const Real zero(0);
    const Real smallNum = limits::SafeMin<Real>();
    const Real bigNum = Real(1) / smallNum;

    bool done = false;
    Real scaleStep;
    while( !done )
    {
        done =
          SafeScaleStep
          ( numerator, denominator, scaleStep, zero, smallNum, bigNum );
        A *= scaleStep;
    }
}

template<typename Field>
void SafeScaleTrapezoid
( Base<Field> numerator, Base<Field> denominator,
  UpperOrLower uplo, SparseMatrix<Field>& A, Int offset )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    const Real zero(0);
    const Real smallNum = limits::SafeMin<Real>();
    const Real bigNum = Real(1) / smallNum;

    bool done = false;
    Real scaleStep;
    while( !done )
    {
        done =
          SafeScaleStep
          ( numerator, denominator, scaleStep, zero, smallNum, bigNum );
        ScaleTrapezoid( scaleStep, uplo, A, offset );
    }
}

template<typename Field>
void SafeScale
( Base<Field> numerator, Base<Field> denominator, DistSparseMatrix<Field>& A )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    const Real zero(0);
    const Real smallNum = limits::SafeMin<Real>();
    const Real bigNum = Real(1) / smallNum;

    bool done = false;
    Real scaleStep;
    while( !done )
    {
        done =
          SafeScaleStep
          ( numerator, denominator, scaleStep, zero, smallNum, bigNum );
        A *= scaleStep;
    }
}

template<typename Field>
void SafeScaleTrapezoid
( Base<Field> numerator, Base<Field> denominator,
  UpperOrLower uplo, DistSparseMatrix<Field>& A, Int offset )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    const Real zero(0);
    const Real smallNum = limits::SafeMin<Real>();
    const Real bigNum = Real(1) / smallNum;

    bool done = false;
    Real scaleStep;
    while( !done )
    {
        done =
          SafeScaleStep
          ( numerator, denominator, scaleStep, zero, smallNum, bigNum );
        ScaleTrapezoid( scaleStep, uplo, A, offset );
    }
}

template<typename Field>
void SafeScale
( Base<Field> numerator, Base<Field> denominator, DistMultiVec<Field>& A )
{
    EL_DEBUG_CSE
    SafeScale( numerator, denominator, A.Matrix() );
}

template<typename Field>
void SafeScaleHermitianTridiag
( Base<Field> numerator, Base<Field> denominator,
  Matrix<Base<Field>>& d, Matrix<Field>& e )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    const Real zero(0);
    const Real smallNum = limits::SafeMin<Real>();
    const Real bigNum = Real(1) / smallNum;

    bool done = false;
    Real scaleStep;
    while( !done )
    {
        done =
          SafeScaleStep
          ( numerator, denominator, scaleStep, zero, smallNum, bigNum );
        d *= scaleStep;
        e *= scaleStep;
    }
}

} // namespace El

#endif // ifndef EL_BLAS_SAFE_SCALE_HPP
