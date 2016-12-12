/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {

template<typename Field>
Base<Field>
TwoNormEstimate( const Matrix<Field>& A, Base<Field> tol, Int maxIts )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    const Int m = A.Height();
    const Int n = A.Width();

    Matrix<Field> x, y;
    Gaussian( y, n, 1 );

    Int numIts=0;
    Real estimate=0, lastEst;
    do
    {
        lastEst = estimate;
        Gemv( NORMAL, Field(1), A, y, x );
        Real xNorm = FrobeniusNorm( x );
        if( xNorm == Real(0) )
        {
            Gaussian( x, m, 1 );
            xNorm = FrobeniusNorm( x );
        }
        x *= Real(1)/xNorm;
        Gemv( ADJOINT, Field(1), A, x, y );
        estimate = FrobeniusNorm( y );
    } while( ++numIts < maxIts && Abs(estimate-lastEst) > tol*Max(m,n) );

    if( Abs(estimate-lastEst) > tol*Max(m,n) )
        RuntimeError("Two-norm estimate did not converge in time");

    return estimate;
}

template<typename Field>
Base<Field> TwoNormEstimate
( const AbstractDistMatrix<Field>& APre, Base<Field> tol, Int maxIts )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;

    DistMatrixReadProxy<Field,Field,MC,MR> AProx( APre );
    auto& A = AProx.GetLocked();

    const Int m = A.Height();
    const Int n = A.Width();

    const Grid& g = APre.Grid();
    DistMatrix<Field> x(g), y(g);
    Gaussian( y, n, 1 );

    Int numIts=0;
    Real estimate=0, lastEst;
    do
    {
        lastEst = estimate;
        Gemv( NORMAL, Field(1), A, y, x );
        Real xNorm = FrobeniusNorm( x );
        if( xNorm == Real(0) )
        {
            Gaussian( x, m, 1 );
            xNorm = FrobeniusNorm( x );
        }
        x *= Real(1)/xNorm;
        Gemv( ADJOINT, Field(1), A, x, y );
        estimate = FrobeniusNorm( y );
    } while( ++numIts < maxIts && Abs(estimate-lastEst) > tol*Max(m,n) );

    if( Abs(estimate-lastEst) > tol*Max(m,n) )
        RuntimeError("Two-norm estimate did not converge in time");

    return estimate;
}

template<typename Field>
Base<Field> TwoNormEstimate( const SparseMatrix<Field>& A, Int basisSize )
{
    auto extremal = ExtremalSingValEst( A, basisSize );
    return extremal.second;
}

template<typename Field>
Base<Field> TwoNormEstimate( const DistSparseMatrix<Field>& A, Int basisSize )
{
    auto extremal = ExtremalSingValEst( A, basisSize );
    return extremal.second;
}

template<typename Field>
Base<Field>
HermitianTwoNormEstimate( const SparseMatrix<Field>& A, Int basisSize )
{
    auto extremal = HermitianExtremalSingValEst( A, basisSize );
    return extremal.second;
}

template<typename Field>
Base<Field>
HermitianTwoNormEstimate( const DistSparseMatrix<Field>& A, Int basisSize )
{
    auto extremal = HermitianExtremalSingValEst( A, basisSize );
    return extremal.second;
}

template<typename Field>
Base<Field> HermitianTwoNormEstimate
( UpperOrLower uplo,
  const Matrix<Field>& A,
  Base<Field> tol,
  Int maxIts )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    const Int n = A.Height();

    Matrix<Field> x, y;
    Zeros( x, n, 1 );
    Gaussian( y, n, 1 );

    Int numIts=0;
    Real estimate=0, lastEst;
    do
    {
        lastEst = estimate;
        Hemv( uplo, Field(1), A, y, Field(0), x );
        Real xNorm = FrobeniusNorm( x );
        if( xNorm == Real(0) )
        {
            Gaussian( x, n, 1 );
            xNorm = FrobeniusNorm( x );
        }
        x *= Real(1)/xNorm;
        Hemv( uplo, Field(1), A, x, Field(0), y );
        estimate = FrobeniusNorm( y );
    } while( ++numIts < maxIts && Abs(estimate-lastEst) > tol*n );

    if( Abs(estimate-lastEst) > tol*n )
        RuntimeError("Two-norm estimate did not converge in time");

    return estimate;
}

template<typename Field>
Base<Field> HermitianTwoNormEstimate
( UpperOrLower uplo,
  const AbstractDistMatrix<Field>& APre,
  Base<Field> tol,
  Int maxIts )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;

    DistMatrixReadProxy<Field,Field,MC,MR> AProx( APre );
    auto& A = AProx.GetLocked();

    const Int n = A.Height();

    const Grid& g = APre.Grid();
    DistMatrix<Field> x(g), y(g);
    Zeros( x, n, 1 );
    Gaussian( y, n, 1 );

    Int numIts=0;
    Real estimate=0, lastEst;
    do
    {
        lastEst = estimate;
        Hemv( uplo, Field(1), A, y, Field(0), x );
        Real xNorm = FrobeniusNorm( x );
        if( xNorm == Real(0) )
        {
            Gaussian( x, n, 1 );
            xNorm = FrobeniusNorm( x );
        }
        x *= Real(1)/xNorm;
        Hemv( uplo, Field(1), A, x, Field(0), y );
        estimate = FrobeniusNorm( y );
    } while( ++numIts < maxIts && Abs(estimate-lastEst) > tol*n );

    if( Abs(estimate-lastEst) > tol*n )
        RuntimeError("Two-norm estimate did not converge in time");

    return estimate;
}

template<typename Field>
Base<Field> SymmetricTwoNormEstimate
( UpperOrLower uplo, const Matrix<Field>& A, Base<Field> tol, Int maxIts )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    const Int n = A.Height();

    Matrix<Field> x, y;
    Zeros( x, n, 1 );
    Gaussian( y, n, 1 );

    Int numIts=0;
    Real estimate=0, lastEst;
    do
    {
        lastEst = estimate;
        Symv( uplo, Field(1), A, y, Field(0), x );
        Real xNorm = FrobeniusNorm( x );
        if( xNorm == Real(0) )
        {
            Gaussian( x, n, 1 );
            xNorm = FrobeniusNorm( x );
        }
        x *= Real(1)/xNorm;
        Conjugate( x );
        Symv( uplo, Field(1), A, x, Field(0), y );
        Conjugate( y );
        estimate = FrobeniusNorm( y );
    } while( ++numIts < maxIts && Abs(estimate-lastEst) > tol*n );

    if( Abs(estimate-lastEst) > tol*n )
        RuntimeError("Two-norm estimate did not converge in time");

    return estimate;
}

template<typename Field>
Base<Field> SymmetricTwoNormEstimate
( UpperOrLower uplo,
  const AbstractDistMatrix<Field>& APre,
  Base<Field> tol,
  Int maxIts )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;

    DistMatrixReadProxy<Field,Field,MC,MR> AProx( APre );
    auto& A = AProx.GetLocked();

    const Int n = A.Height();

    const Grid& g = APre.Grid();
    DistMatrix<Field> x(g), y(g);
    Zeros( x, n, 1 );
    Gaussian( y, n, 1 );

    Int numIts=0;
    Real estimate=0, lastEst;
    do
    {
        lastEst = estimate;
        Symv( uplo, Field(1), A, y, Field(0), x );
        Real xNorm = FrobeniusNorm( x );
        if( xNorm == Real(0) )
        {
            Gaussian( x, n, 1 );
            xNorm = FrobeniusNorm( x );
        }
        x *= Real(1)/xNorm;
        Conjugate( x );
        Symv( uplo, Field(1), A, x, Field(0), y );
        Conjugate( y );
        estimate = FrobeniusNorm( y );
    } while( ++numIts < maxIts && Abs(estimate-lastEst) > tol*n );

    if( Abs(estimate-lastEst) > tol*n )
        RuntimeError("Two-norm estimate did not converge in time");

    return estimate;
}

#define PROTO(Field) \
  template Base<Field> TwoNormEstimate \
  ( const Matrix<Field>& A, Base<Field> tol, Int maxIts ); \
  template Base<Field> TwoNormEstimate \
  ( const AbstractDistMatrix<Field>& A, Base<Field> tol, Int maxIts ); \
  template Base<Field> TwoNormEstimate \
  ( const SparseMatrix<Field>& A, Int basisSize ); \
  template Base<Field> TwoNormEstimate \
  ( const DistSparseMatrix<Field>& A, Int basisSize ); \
  template Base<Field> HermitianTwoNormEstimate \
  ( const SparseMatrix<Field>& A, Int basisSize ); \
  template Base<Field> HermitianTwoNormEstimate \
  ( const DistSparseMatrix<Field>& A, Int basisSize ); \
  template Base<Field> HermitianTwoNormEstimate \
  ( UpperOrLower uplo, const Matrix<Field>& A, Base<Field> tol, Int maxIts ); \
  template Base<Field> HermitianTwoNormEstimate \
  ( UpperOrLower uplo, const AbstractDistMatrix<Field>& A, Base<Field> tol, \
    Int maxIts ); \
  template Base<Field> SymmetricTwoNormEstimate \
  ( UpperOrLower uplo, const Matrix<Field>& A, Base<Field> tol, Int maxIts ); \
  template Base<Field> SymmetricTwoNormEstimate \
  ( UpperOrLower uplo, const AbstractDistMatrix<Field>& A, Base<Field> tol, \
    Int maxIts );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
