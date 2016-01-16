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
void Lanczos
( const SparseMatrix<F>& A,
        Matrix<Base<F>>& T,
        Int basisSize )
{
    DEBUG_ONLY(CSE cse("Lanczos"))
    const Int n = A.Height();
    if( n != A.Width() )
        LogicError("A was not square");
    
    auto applyA =
      [&]( const Matrix<F>& X, Matrix<F>& Y )
      {
          Zeros( Y, n, X.Width() );
          Multiply( NORMAL, F(1), A, X, F(0), Y );
      };
    Lanczos<F>( n, applyA, T, basisSize );
}

template<typename F>
Base<F> LanczosDecomp
( const SparseMatrix<F>& A,
        Matrix<F>& V, 
        Matrix<Base<F>>& T,
        Matrix<F>& v,
        Int basisSize )
{
    DEBUG_ONLY(CSE cse("LanczosDecomp"))
    const Int n = A.Height();
    if( n != A.Width() )
        LogicError("A was not square");

    auto applyA =
      [&]( const Matrix<F>& X, Matrix<F>& Y )
      {
          Zeros( Y, n, X.Width() );
          Multiply( NORMAL, F(1), A, X, F(0), Y );
      };
    return LanczosDecomp( n, applyA, V, T, v, basisSize );
}

template<typename F>
void Lanczos
( const DistSparseMatrix<F>& A,
        ElementalMatrix<Base<F>>& T,
        Int basisSize )
{
    DEBUG_ONLY(CSE cse("Lanczos"))
    const Int n = A.Height();
    if( n != A.Width() )
        LogicError("A was not square");

    auto applyA =
      [&]( const DistMultiVec<F>& X, DistMultiVec<F>& Y )
      {
          Zeros( Y, n, X.Width() );
          Multiply( NORMAL, F(1), A, X, F(0), Y );
      };
    Lanczos<F>( n, applyA, T, basisSize );
}

template<typename F>
Base<F> LanczosDecomp
( const DistSparseMatrix<F>& A,
        DistMultiVec<F>& V, 
        ElementalMatrix<Base<F>>& T,
        DistMultiVec<F>& v,
        Int basisSize )
{
    DEBUG_ONLY(CSE cse("LanczosDecomp"))
    const Int n = A.Height();
    if( n != A.Width() )
        LogicError("A was not square");

    auto applyA =
      [&]( const DistMultiVec<F>& X, DistMultiVec<F>& Y )
      {
          Zeros( Y, n, X.Width() );
          Multiply( NORMAL, F(1), A, X, F(0), Y );
      };
    return LanczosDecomp( n, applyA, V, T, v, basisSize );
}

#define PROTO(F) \
  template void Lanczos \
  ( const SparseMatrix<F>& A, \
          Matrix<Base<F>>& T, \
          Int basisSize ); \
  template void Lanczos \
  ( const DistSparseMatrix<F>& A, \
          ElementalMatrix<Base<F>>& T, \
          Int basisSize ); \
  template Base<F> LanczosDecomp \
  ( const SparseMatrix<F>& A, \
          Matrix<F>& V, \
          Matrix<Base<F>>& T, \
          Matrix<F>& v, \
          Int basisSize ); \
  template Base<F> LanczosDecomp \
  ( const DistSparseMatrix<F>& A, \
          DistMultiVec<F>& V, \
          ElementalMatrix<Base<F>>& T, \
          DistMultiVec<F>& v, \
          Int basisSize );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
