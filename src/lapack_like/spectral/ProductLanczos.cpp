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
void ProductLanczos
( const SparseMatrix<Field>& A,
        Matrix<Base<Field>>& T,
        Int basisSize )
{
    EL_DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();

    Matrix<Field> S;

    // Cache the adjoint
    // -----------------
    SparseMatrix<Field> AAdj;
    Adjoint( A, AAdj );

    if( m >= n )
    {
        auto applyA =
          [&]( const Matrix<Field>& X, Matrix<Field>& Y )
          {
              Zeros( S, m, X.Width() );
              Multiply( NORMAL, Field(1), A,    X, Field(0), S );
              Zeros( Y, n, X.Width() );
              Multiply( NORMAL, Field(1), AAdj, S, Field(0), Y );
          };
        Lanczos<Field>( n, applyA, T, basisSize );
    }
    else
    {
        auto applyA =
          [&]( const Matrix<Field>& X, Matrix<Field>& Y )
          {
              Zeros( S, n, X.Width() );
              Multiply( NORMAL, Field(1), AAdj, X, Field(0), S );
              Zeros( Y, m, X.Width() );
              Multiply( NORMAL, Field(1), A,    S, Field(0), Y );
          };
        Lanczos<Field>( m, applyA, T, basisSize );
    }
}

template<typename Field>
Base<Field> ProductLanczosDecomp
( const SparseMatrix<Field>& A,
        Matrix<Field>& V,
        Matrix<Base<Field>>& T,
        Matrix<Field>& v,
        Int basisSize )
{
    EL_DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();

    Matrix<Field> S;

    // Cache the adjoint
    // -----------------
    SparseMatrix<Field> AAdj;
    Adjoint( A, AAdj );

    if( m >= n )
    {
        auto applyA =
          [&]( const Matrix<Field>& X, Matrix<Field>& Y )
          {
              Zeros( S, m, X.Width() );
              Multiply( NORMAL, Field(1), A,    X, Field(0), S );
              Zeros( Y, n, X.Width() );
              Multiply( NORMAL, Field(1), AAdj, S, Field(0), Y );
          };
        return LanczosDecomp( n, applyA, V, T, v, basisSize );
    }
    else
    {
        auto applyA =
          [&]( const Matrix<Field>& X, Matrix<Field>& Y )
          {
              Zeros( S, n, X.Width() );
              Multiply( NORMAL, Field(1), AAdj, X, Field(0), S );
              Zeros( Y, m, X.Width() );
              Multiply( NORMAL, Field(1), A,    S, Field(0), Y );
          };
        return LanczosDecomp( m, applyA, V, T, v, basisSize );
    }
}

template<typename Field>
void ProductLanczos
( const DistSparseMatrix<Field>& A,
        AbstractDistMatrix<Base<Field>>& T,
        Int basisSize )
{
    EL_DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();
    const Grid& grid = A.Grid();

    // Cache the adjoint
    // -----------------
    DistSparseMatrix<Field> AAdj(grid);
    Adjoint( A, AAdj );

    DistMultiVec<Field> S(grid);

    if( m >= n )
    {
        auto applyA =
          [&]( const DistMultiVec<Field>& X, DistMultiVec<Field>& Y )
          {
              Zeros( S, m, X.Width() );
              Multiply( NORMAL, Field(1), A,    X, Field(0), S );
              Zeros( Y, n, X.Width() );
              Multiply( NORMAL, Field(1), AAdj, S, Field(0), Y );
          };
        Lanczos<Field>( n, applyA, T, basisSize );
    }
    else
    {
        auto applyA =
          [&]( const DistMultiVec<Field>& X, DistMultiVec<Field>& Y )
          {
              Zeros( S, n, X.Width() );
              Multiply( NORMAL, Field(1), AAdj, X, Field(0), S );
              Zeros( Y, m, X.Width() );
              Multiply( NORMAL, Field(1), A,    S, Field(0), Y );
          };
        Lanczos<Field>( m, applyA, T, basisSize );
    }
}

template<typename Field>
Base<Field> ProductLanczosDecomp
( const DistSparseMatrix<Field>& A,
        DistMultiVec<Field>& V,
        AbstractDistMatrix<Base<Field>>& T,
        DistMultiVec<Field>& v,
        Int basisSize )
{
    EL_DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();
    const Grid& grid = A.Grid();

    DistMultiVec<Field> S(grid);

    // Cache the adjoint
    // -----------------
    DistSparseMatrix<Field> AAdj(grid);
    Adjoint( A, AAdj );

    if( m >= n )
    {
        auto applyA =
          [&]( const DistMultiVec<Field>& X, DistMultiVec<Field>& Y )
          {
              Zeros( S, m, X.Width() );
              Multiply( NORMAL, Field(1), A,    X, Field(0), S );
              Zeros( Y, n, X.Width() );
              Multiply( NORMAL, Field(1), AAdj, S, Field(0), Y );
          };
        return LanczosDecomp( n, applyA, V, T, v, basisSize );
    }
    else
    {
        auto applyA =
          [&]( const DistMultiVec<Field>& X, DistMultiVec<Field>& Y )
          {
              Zeros( S, n, X.Width() );
              Multiply( NORMAL, Field(1), AAdj, X, Field(0), S );
              Zeros( Y, m, X.Width() );
              Multiply( NORMAL, Field(1), A,    S, Field(0), Y );
          };
        return LanczosDecomp( m, applyA, V, T, v, basisSize );
    }
}

#define PROTO(Field) \
  template void ProductLanczos \
  ( const SparseMatrix<Field>& A, \
          Matrix<Base<Field>>& T, \
          Int basisSize ); \
  template void ProductLanczos \
  ( const DistSparseMatrix<Field>& A, \
          AbstractDistMatrix<Base<Field>>& T, \
          Int basisSize ); \
  template Base<Field> ProductLanczosDecomp \
  ( const SparseMatrix<Field>& A, \
          Matrix<Field>& V, \
          Matrix<Base<Field>>& T, \
          Matrix<Field>& v, \
          Int basisSize ); \
  template Base<Field> ProductLanczosDecomp \
  ( const DistSparseMatrix<Field>& A, \
          DistMultiVec<Field>& V, \
          AbstractDistMatrix<Base<Field>>& T, \
          DistMultiVec<Field>& v, \
          Int basisSize );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
