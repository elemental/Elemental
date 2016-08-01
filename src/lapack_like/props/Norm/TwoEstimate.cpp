/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {

template<typename F>
Base<F> TwoNormEstimate( const Matrix<F>& A, Base<F> tol, Int maxIts )
{
    DEBUG_CSE
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();

    Matrix<F> x, y;
    Gaussian( y, n, 1 );
    
    Int numIts=0;
    Real estimate=0, lastEst;
    do
    {
        lastEst = estimate;
        Gemv( NORMAL, F(1), A, y, x );
        Real xNorm = FrobeniusNorm( x );
        if( xNorm == Real(0) )
        {
            Gaussian( x, m, 1 );    
            xNorm = FrobeniusNorm( x );
        }
        x *= Real(1)/xNorm;
        Gemv( ADJOINT, F(1), A, x, y );
        estimate = FrobeniusNorm( y );
    } while( ++numIts < maxIts && Abs(estimate-lastEst) > tol*Max(m,n) );

    if( Abs(estimate-lastEst) > tol*Max(m,n) )
        RuntimeError("Two-norm estimate did not converge in time");

    return estimate;
}

template<typename F>
Base<F> TwoNormEstimate
( const ElementalMatrix<F>& APre, Base<F> tol, Int maxIts )
{
    DEBUG_CSE
    typedef Base<F> Real;

    DistMatrixReadProxy<F,F,MC,MR> AProx( APre );
    auto& A = AProx.GetLocked();

    const Int m = A.Height();
    const Int n = A.Width();

    const Grid& g = APre.Grid();
    DistMatrix<F> x(g), y(g);
    Gaussian( y, n, 1 );
    
    Int numIts=0;
    Real estimate=0, lastEst;
    do
    {
        lastEst = estimate;
        Gemv( NORMAL, F(1), A, y, x );
        Real xNorm = FrobeniusNorm( x );
        if( xNorm == Real(0) )
        {
            Gaussian( x, m, 1 );    
            xNorm = FrobeniusNorm( x );
        }
        x *= Real(1)/xNorm;
        Gemv( ADJOINT, F(1), A, x, y );
        estimate = FrobeniusNorm( y );
    } while( ++numIts < maxIts && Abs(estimate-lastEst) > tol*Max(m,n) );

    if( Abs(estimate-lastEst) > tol*Max(m,n) )
        RuntimeError("Two-norm estimate did not converge in time");

    return estimate;
}

template<typename F>
Base<F> TwoNormEstimate( const SparseMatrix<F>& A, Int basisSize )
{
    auto extremal = ExtremalSingValEst( A, basisSize );
    return extremal.second;
}

template<typename F>
Base<F> TwoNormEstimate( const DistSparseMatrix<F>& A, Int basisSize )
{
    auto extremal = ExtremalSingValEst( A, basisSize );
    return extremal.second;
}

template<typename F>
Base<F> HermitianTwoNormEstimate( const SparseMatrix<F>& A, Int basisSize )
{
    auto extremal = HermitianExtremalSingValEst( A, basisSize );
    return extremal.second;
}

template<typename F>
Base<F> HermitianTwoNormEstimate( const DistSparseMatrix<F>& A, Int basisSize )
{
    auto extremal = HermitianExtremalSingValEst( A, basisSize );
    return extremal.second;
}

template<typename F>
Base<F> HermitianTwoNormEstimate
( UpperOrLower uplo,
  const Matrix<F>& A,
  Base<F> tol,
  Int maxIts )
{
    DEBUG_CSE
    typedef Base<F> Real;
    const Int n = A.Height();

    Matrix<F> x, y;
    Zeros( x, n, 1 );
    Gaussian( y, n, 1 );
    
    Int numIts=0;
    Real estimate=0, lastEst;
    do
    {
        lastEst = estimate;
        Hemv( uplo, F(1), A, y, F(0), x );
        Real xNorm = FrobeniusNorm( x );
        if( xNorm == Real(0) )
        {
            Gaussian( x, n, 1 );    
            xNorm = FrobeniusNorm( x );
        }
        x *= Real(1)/xNorm;
        Hemv( uplo, F(1), A, x, F(0), y );
        estimate = FrobeniusNorm( y );
    } while( ++numIts < maxIts && Abs(estimate-lastEst) > tol*n );

    if( Abs(estimate-lastEst) > tol*n )
        RuntimeError("Two-norm estimate did not converge in time");

    return estimate;
}

template<typename F>
Base<F> HermitianTwoNormEstimate
( UpperOrLower uplo,
  const ElementalMatrix<F>& APre,
  Base<F> tol, 
  Int maxIts )
{
    DEBUG_CSE
    typedef Base<F> Real;

    DistMatrixReadProxy<F,F,MC,MR> AProx( APre );
    auto& A = AProx.GetLocked();

    const Int n = A.Height();

    const Grid& g = APre.Grid();
    DistMatrix<F> x(g), y(g);
    Zeros( x, n, 1 );
    Gaussian( y, n, 1 );
    
    Int numIts=0;
    Real estimate=0, lastEst;
    do
    {
        lastEst = estimate;
        Hemv( uplo, F(1), A, y, F(0), x );
        Real xNorm = FrobeniusNorm( x );
        if( xNorm == Real(0) )
        {
            Gaussian( x, n, 1 );    
            xNorm = FrobeniusNorm( x );
        }
        x *= Real(1)/xNorm;
        Hemv( uplo, F(1), A, x, F(0), y );
        estimate = FrobeniusNorm( y );
    } while( ++numIts < maxIts && Abs(estimate-lastEst) > tol*n );

    if( Abs(estimate-lastEst) > tol*n )
        RuntimeError("Two-norm estimate did not converge in time");

    return estimate;
}

template<typename F>
Base<F> SymmetricTwoNormEstimate
( UpperOrLower uplo, const Matrix<F>& A, Base<F> tol, Int maxIts )
{
    DEBUG_CSE
    typedef Base<F> Real;
    const Int n = A.Height();

    Matrix<F> x, y;
    Zeros( x, n, 1 );
    Gaussian( y, n, 1 );
    
    Int numIts=0;
    Real estimate=0, lastEst;
    do
    {
        lastEst = estimate;
        Symv( uplo, F(1), A, y, F(0), x );
        Real xNorm = FrobeniusNorm( x );
        if( xNorm == Real(0) )
        {
            Gaussian( x, n, 1 );    
            xNorm = FrobeniusNorm( x );
        }
        x *= Real(1)/xNorm;
        Conjugate( x );
        Symv( uplo, F(1), A, x, F(0), y );
        Conjugate( y );
        estimate = FrobeniusNorm( y );
    } while( ++numIts < maxIts && Abs(estimate-lastEst) > tol*n );

    if( Abs(estimate-lastEst) > tol*n )
        RuntimeError("Two-norm estimate did not converge in time");

    return estimate;
}

template<typename F>
Base<F> SymmetricTwoNormEstimate
( UpperOrLower uplo,
  const ElementalMatrix<F>& APre,
  Base<F> tol, 
  Int maxIts )
{
    DEBUG_CSE
    typedef Base<F> Real;

    DistMatrixReadProxy<F,F,MC,MR> AProx( APre );
    auto& A = AProx.GetLocked();

    const Int n = A.Height();

    const Grid& g = APre.Grid();
    DistMatrix<F> x(g), y(g);
    Zeros( x, n, 1 );
    Gaussian( y, n, 1 );
    
    Int numIts=0;
    Real estimate=0, lastEst;
    do
    {
        lastEst = estimate;
        Symv( uplo, F(1), A, y, F(0), x );
        Real xNorm = FrobeniusNorm( x );
        if( xNorm == Real(0) )
        {
            Gaussian( x, n, 1 );    
            xNorm = FrobeniusNorm( x );
        }
        x *= Real(1)/xNorm;
        Conjugate( x );
        Symv( uplo, F(1), A, x, F(0), y );
        Conjugate( y );
        estimate = FrobeniusNorm( y );
    } while( ++numIts < maxIts && Abs(estimate-lastEst) > tol*n );

    if( Abs(estimate-lastEst) > tol*n )
        RuntimeError("Two-norm estimate did not converge in time");

    return estimate;
}

#define PROTO(F) \
  template Base<F> TwoNormEstimate \
  ( const Matrix<F>& A, Base<F> tol, Int maxIts ); \
  template Base<F> TwoNormEstimate \
  ( const ElementalMatrix<F>& A, Base<F> tol, Int maxIts ); \
  template Base<F> TwoNormEstimate \
  ( const SparseMatrix<F>& A, Int basisSize ); \
  template Base<F> TwoNormEstimate \
  ( const DistSparseMatrix<F>& A, Int basisSize ); \
  template Base<F> HermitianTwoNormEstimate \
  ( const SparseMatrix<F>& A, Int basisSize ); \
  template Base<F> HermitianTwoNormEstimate \
  ( const DistSparseMatrix<F>& A, Int basisSize ); \
  template Base<F> HermitianTwoNormEstimate \
  ( UpperOrLower uplo, const Matrix<F>& A, Base<F> tol, Int maxIts ); \
  template Base<F> HermitianTwoNormEstimate \
  ( UpperOrLower uplo, const ElementalMatrix<F>& A, Base<F> tol, \
    Int maxIts ); \
  template Base<F> SymmetricTwoNormEstimate \
  ( UpperOrLower uplo, const Matrix<F>& A, Base<F> tol, Int maxIts ); \
  template Base<F> SymmetricTwoNormEstimate \
  ( UpperOrLower uplo, const ElementalMatrix<F>& A, Base<F> tol, \
    Int maxIts );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
