/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

#include "./Gemm/NN.hpp"
#include "./Gemm/NT.hpp"
#include "./Gemm/TN.hpp"
#include "./Gemm/TT.hpp"

namespace El {

template<typename T>
void Gemm
( Orientation orientA, Orientation orientB,
  T alpha, const Matrix<T>& A, const Matrix<T>& B, T beta, Matrix<T>& C )
{
    DEBUG_ONLY(
      CallStackEntry cse("Gemm");
      if( orientA == NORMAL && orientB == NORMAL )
      {
          if( A.Height() != C.Height() ||
              B.Width()  != C.Width()  ||
              A.Width()  != B.Height() )
              LogicError("Nonconformal GemmNN");
      }
      else if( orientA == NORMAL )
      {
          if( A.Height() != C.Height() ||
              B.Height() != C.Width()  ||
              A.Width()  != B.Width() )
              LogicError("Nonconformal GemmN(T/C)");
      }
      else if( orientB == NORMAL )
      {
          if( A.Width()  != C.Height() ||
              B.Width()  != C.Width()  ||
              A.Height() != B.Height() )
              LogicError("Nonconformal Gemm(T/C)N");
      }
      else
      {
          if( A.Width()  != C.Height() ||
              B.Height() != C.Width()  ||
              A.Height() != B.Width() )
              LogicError("Nonconformal Gemm(T/C)(T/C)");
      }
    )
    const char transA = OrientationToChar( orientA );
    const char transB = OrientationToChar( orientB );
    const Int m = C.Height();
    const Int n = C.Width();
    const Int k = ( orientA == NORMAL ? A.Width() : A.Height() );
    if( k != 0 )
    {
        blas::Gemm
        ( transA, transB, m, n, k,
          alpha, A.LockedBuffer(), A.LDim(), B.LockedBuffer(), B.LDim(),
          beta,  C.Buffer(),       C.LDim() );
    }
    else
    {
        Scale( beta, C );
    }
}

template<typename T>
void Gemm
( Orientation orientA, Orientation orientB,
  T alpha, const Matrix<T>& A, const Matrix<T>& B, Matrix<T>& C )
{
    DEBUG_ONLY(CallStackEntry cse("Gemm"))
    const Int m = ( orientA==NORMAL ? A.Height() : A.Width() );
    const Int n = ( orientB==NORMAL ? B.Width() : B.Height() );
    Zeros( C, m, n );
    Gemm( orientA, orientB, alpha, A, B, T(0), C );
}

template<typename T>
void Gemm
( Orientation orientA, Orientation orientB,
  T alpha, const AbstractDistMatrix<T>& A, const AbstractDistMatrix<T>& B,
  T beta,        AbstractDistMatrix<T>& C, GemmAlgorithm alg )
{
    DEBUG_ONLY(CallStackEntry cse("Gemm"))
    if( orientA == NORMAL && orientB == NORMAL )
    {
        if( alg == GEMM_CANNON )
            gemm::Cannon_NN( alpha, A, B, beta, C );
        else 
            gemm::SUMMA_NN( alpha, A, B, beta, C, alg );
    }
    else if( orientA == NORMAL )
    {
        gemm::SUMMA_NT( orientB, alpha, A, B, beta, C, alg );
    }
    else if( orientB == NORMAL )
    {
        gemm::SUMMA_TN( orientA, alpha, A, B, beta, C, alg );
    }
    else
    {
        gemm::SUMMA_TT( orientA, orientB, alpha, A, B, beta, C, alg );
    }
}

template<typename T>
void Gemm
( Orientation orientA, Orientation orientB,
  T alpha, const AbstractDistMatrix<T>& A, const AbstractDistMatrix<T>& B,
                 AbstractDistMatrix<T>& C, GemmAlgorithm alg )
{
    DEBUG_ONLY(CallStackEntry cse("Gemm"))
    const Int m = ( orientA==NORMAL ? A.Height() : A.Width() );
    const Int n = ( orientB==NORMAL ? B.Width() : B.Height() );
    Zeros( C, m, n );
    Gemm( orientA, orientB, alpha, A, B, T(0), C, alg );
}

#define PROTO(T) \
  template void Gemm \
  ( Orientation orientA, Orientation orientB, \
    T alpha, const Matrix<T>& A, const Matrix<T>& B, \
    T beta,        Matrix<T>& C ); \
  template void Gemm \
  ( Orientation orientA, Orientation orientB, \
    T alpha, const Matrix<T>& A, const Matrix<T>& B, \
                   Matrix<T>& C ); \
  template void Gemm \
  ( Orientation orientA, Orientation orientB, \
    T alpha, const AbstractDistMatrix<T>& A, const AbstractDistMatrix<T>& B, \
    T beta,        AbstractDistMatrix<T>& C, GemmAlgorithm alg ); \
  template void Gemm \
  ( Orientation orientA, Orientation orientB, \
    T alpha, const AbstractDistMatrix<T>& A, const AbstractDistMatrix<T>& B, \
                   AbstractDistMatrix<T>& C, GemmAlgorithm alg );

#define EL_ENABLE_QUAD
#include "El/macros/Instantiate.h"

} // namespace El
