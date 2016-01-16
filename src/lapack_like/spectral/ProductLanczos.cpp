/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename F>
void ProductLanczos
( const SparseMatrix<F>& A,
        Matrix<Base<F>>& T,
        Int basisSize )
{
    DEBUG_ONLY(CSE cse("ProductLanczos"))
    const Int m = A.Height();
    const Int n = A.Width();

    Matrix<F> S;

    // Cache the adjoint
    // -----------------
    SparseMatrix<F> AAdj;
    Adjoint( A, AAdj );

    if( m >= n )
    {
        auto applyA =
          [&]( const Matrix<F>& X, Matrix<F>& Y )
          {
              Zeros( S, m, X.Width() );
              Multiply( NORMAL, F(1), A,    X, F(0), S );
              Zeros( Y, n, X.Width() );
              Multiply( NORMAL, F(1), AAdj, S, F(0), Y );
          };
        Lanczos<F>( n, applyA, T, basisSize );
    }
    else
    {
        auto applyA =
          [&]( const Matrix<F>& X, Matrix<F>& Y )
          {
              Zeros( S, n, X.Width() );
              Multiply( NORMAL, F(1), AAdj, X, F(0), S );
              Zeros( Y, m, X.Width() );
              Multiply( NORMAL, F(1), A,    S, F(0), Y );
          };
        Lanczos<F>( m, applyA, T, basisSize );
    }
}

template<typename F>
Base<F> ProductLanczosDecomp
( const SparseMatrix<F>& A,
        Matrix<F>& V, 
        Matrix<Base<F>>& T,
        Matrix<F>& v,
        Int basisSize )
{
    DEBUG_ONLY(CSE cse("ProductLanczosDecomp"))
    const Int m = A.Height();
    const Int n = A.Width();

    Matrix<F> S;
    
    // Cache the adjoint
    // -----------------
    SparseMatrix<F> AAdj;
    Adjoint( A, AAdj );
    
    if( m >= n )
    {   
        auto applyA =
          [&]( const Matrix<F>& X, Matrix<F>& Y )
          {   
              Zeros( S, m, X.Width() );
              Multiply( NORMAL, F(1), A,    X, F(0), S );
              Zeros( Y, n, X.Width() );
              Multiply( NORMAL, F(1), AAdj, S, F(0), Y );
          };
        return LanczosDecomp( n, applyA, V, T, v, basisSize );
    }
    else
    {
        auto applyA =
          [&]( const Matrix<F>& X, Matrix<F>& Y )
          {   
              Zeros( S, n, X.Width() );
              Multiply( NORMAL, F(1), AAdj, X, F(0), S );
              Zeros( Y, m, X.Width() );
              Multiply( NORMAL, F(1), A,    S, F(0), Y );
          };
        return LanczosDecomp( m, applyA, V, T, v, basisSize );
    }
}

template<typename F>
void ProductLanczos
( const DistSparseMatrix<F>& A,
        ElementalMatrix<Base<F>>& T,
        Int basisSize )
{
    DEBUG_ONLY(CSE cse("ProductLanczos"))
    const Int m = A.Height();
    const Int n = A.Width();
    mpi::Comm comm = A.Comm();

    // Cache the adjoint
    // -----------------
    DistSparseMatrix<F> AAdj(comm);
    Adjoint( A, AAdj );
    
    DistMultiVec<F> S(comm);

    if( m >= n )
    {
        auto applyA =
          [&]( const DistMultiVec<F>& X, DistMultiVec<F>& Y )
          {
              Zeros( S, m, X.Width() );
              Multiply( NORMAL, F(1), A,    X, F(0), S );
              Zeros( Y, n, X.Width() );
              Multiply( NORMAL, F(1), AAdj, S, F(0), Y );
          };
        Lanczos<F>( n, applyA, T, basisSize );
    }
    else
    {
        auto applyA =
          [&]( const DistMultiVec<F>& X, DistMultiVec<F>& Y )
          {
              Zeros( S, n, X.Width() );
              Multiply( NORMAL, F(1), AAdj, X, F(0), S );
              Zeros( Y, m, X.Width() );
              Multiply( NORMAL, F(1), A,    S, F(0), Y );
          };
        Lanczos<F>( m, applyA, T, basisSize );
    }
}

template<typename F>
Base<F> ProductLanczosDecomp
( const DistSparseMatrix<F>& A,
        DistMultiVec<F>& V, 
        ElementalMatrix<Base<F>>& T,
        DistMultiVec<F>& v,
        Int basisSize )
{
    DEBUG_ONLY(CSE cse("ProductLanczosDecomp"))
    const Int m = A.Height();
    const Int n = A.Width();
    mpi::Comm comm = A.Comm();

    DistMultiVec<F> S(comm);
    
    // Cache the adjoint
    // -----------------
    DistSparseMatrix<F> AAdj(comm);
    Adjoint( A, AAdj );
    
    if( m >= n )
    {   
        auto applyA =
          [&]( const DistMultiVec<F>& X, DistMultiVec<F>& Y )
          {   
              Zeros( S, m, X.Width() );
              Multiply( NORMAL, F(1), A,    X, F(0), S );
              Zeros( Y, n, X.Width() );
              Multiply( NORMAL, F(1), AAdj, S, F(0), Y );
          };
        return LanczosDecomp( n, applyA, V, T, v, basisSize );
    }
    else
    {
        auto applyA =
          [&]( const DistMultiVec<F>& X, DistMultiVec<F>& Y )
          {   
              Zeros( S, n, X.Width() );
              Multiply( NORMAL, F(1), AAdj, X, F(0), S );
              Zeros( Y, m, X.Width() );
              Multiply( NORMAL, F(1), A,    S, F(0), Y );
          };
        return LanczosDecomp( m, applyA, V, T, v, basisSize );
    }
}

#define PROTO(F) \
  template void ProductLanczos \
  ( const SparseMatrix<F>& A, \
          Matrix<Base<F>>& T, \
          Int basisSize ); \
  template void ProductLanczos \
  ( const DistSparseMatrix<F>& A, \
          ElementalMatrix<Base<F>>& T, \
          Int basisSize ); \
  template Base<F> ProductLanczosDecomp \
  ( const SparseMatrix<F>& A, \
          Matrix<F>& V, \
          Matrix<Base<F>>& T, \
          Matrix<F>& v, \
          Int basisSize ); \
  template Base<F> ProductLanczosDecomp \
  ( const DistSparseMatrix<F>& A, \
          DistMultiVec<F>& V, \
          ElementalMatrix<Base<F>>& T, \
          DistMultiVec<F>& v, \
          Int basisSize );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
