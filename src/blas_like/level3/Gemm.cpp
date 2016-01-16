/*
   Copyright (c) 2009-2016, Jack Poulson
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
  T alpha, const Matrix<T>& A,
           const Matrix<T>& B, 
  T beta,        Matrix<T>& C )
{
    DEBUG_ONLY(CSE cse("Gemm"))
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
    const char transA = OrientationToChar( orientA );
    const char transB = OrientationToChar( orientB );
    const Int m = C.Height();
    const Int n = C.Width();
    const Int k = ( orientA == NORMAL ? A.Width() : A.Height() );
    if( k != 0 )
    {
        blas::Gemm
        ( transA, transB, m, n, k,
          alpha, A.LockedBuffer(), A.LDim(),
                 B.LockedBuffer(), B.LDim(),
          beta,  C.Buffer(),       C.LDim() );
    }
    else
    {
        C *= beta;
    }
}

template<typename T>
void Gemm
( Orientation orientA, Orientation orientB,
  T alpha, const Matrix<T>& A,
           const Matrix<T>& B, 
                 Matrix<T>& C )
{
    DEBUG_ONLY(CSE cse("Gemm"))
    const Int m = ( orientA==NORMAL ? A.Height() : A.Width() );
    const Int n = ( orientB==NORMAL ? B.Width() : B.Height() );
    Zeros( C, m, n );
    Gemm( orientA, orientB, alpha, A, B, T(0), C );
}

template<typename T>
void Gemm
( Orientation orientA, Orientation orientB,
  T alpha, const ElementalMatrix<T>& A,
           const ElementalMatrix<T>& B,
  T beta,        ElementalMatrix<T>& C, 
  GemmAlgorithm alg )
{
    DEBUG_ONLY(CSE cse("Gemm"))
    C *= beta;
    if( orientA == NORMAL && orientB == NORMAL )
    {
        if( alg == GEMM_CANNON )
            gemm::Cannon_NN( alpha, A, B, C );
        else 
            gemm::SUMMA_NN( alpha, A, B, C, alg );
    }
    else if( orientA == NORMAL )
    {
        gemm::SUMMA_NT( orientB, alpha, A, B, C, alg );
    }
    else if( orientB == NORMAL )
    {
        gemm::SUMMA_TN( orientA, alpha, A, B, C, alg );
    }
    else
    {
        gemm::SUMMA_TT( orientA, orientB, alpha, A, B, C, alg );
    }
}

template<typename T>
void Gemm
( Orientation orientA, Orientation orientB,
  T alpha, const ElementalMatrix<T>& A,
           const ElementalMatrix<T>& B,
                 ElementalMatrix<T>& C, GemmAlgorithm alg )
{
    DEBUG_ONLY(CSE cse("Gemm"))
    const Int m = ( orientA==NORMAL ? A.Height() : A.Width() );
    const Int n = ( orientB==NORMAL ? B.Width() : B.Height() );
    Zeros( C, m, n );
    Gemm( orientA, orientB, alpha, A, B, T(0), C, alg );
}

template<typename T>
void LocalGemm
( Orientation orientA, Orientation orientB,
  T alpha, const ElementalMatrix<T>& A,
           const ElementalMatrix<T>& B,
  T beta,        ElementalMatrix<T>& C )
{
    DEBUG_ONLY(
      CSE cse("LocalGemm");
      if( orientA == NORMAL && orientB == NORMAL )
      {
          if( A.ColDist() != C.ColDist() ||
              A.RowDist() != B.ColDist() ||
              B.RowDist() != C.RowDist() )
              LogicError
              ("Tried to form C[",C.ColDist(),",",C.RowDist(),"] := "
               "A[",A.ColDist(),",",A.RowDist(),"] "
               "B[",B.ColDist(),",",B.RowDist(),"]");
          if( A.ColAlign() != C.ColAlign() )
              LogicError("A's cols must align with C's rows");
          if( A.RowAlign() != B.ColAlign() )
              LogicError("A's rows must align with B's cols");
          if( B.RowAlign() != C.RowAlign() )
              LogicError("B's rows must align with C's rows");
          if( A.Height() != C.Height() ||
              A.Width() != B.Height() ||
              B.Width() != C.Width() )
              LogicError
              ("Nonconformal LocalGemmNN:\n",
               DimsString(A,"A"),"\n",
               DimsString(B,"B"),"\n",
               DimsString(C,"C"));
      }
      else if( orientA == NORMAL )
      {
          if( A.ColDist() != C.ColDist() ||
              A.RowDist() != B.RowDist() ||
              B.ColDist() != C.RowDist() )
              LogicError
              ("Tried to form C[",C.ColDist(),",",C.RowDist(),"] := "
               "A[",A.ColDist(),",",A.RowDist(),"] "
               "B[",B.ColDist(),",",B.RowDist(),"]'");
          if( A.ColAlign() != C.ColAlign() )
              LogicError("A's cols must align with C's rows");
          if( A.RowAlign() != B.RowAlign() )
              LogicError("A's rows must align with B's rows");
          if( B.ColAlign() != C.RowAlign() )
              LogicError("B's cols must align with C's rows");
          if( A.Height() != C.Height() ||
              A.Width() != B.Width() ||
              B.Height() != C.Width() )
              LogicError
              ("Nonconformal LocalGemmNT:\n",
               DimsString(A,"A"),"\n",
               DimsString(B,"B"),"\n",
               DimsString(C,"C"));
      }
      else if( orientB == NORMAL )
      {
          if( A.RowDist() != C.ColDist() ||
              A.ColDist() != B.ColDist() ||
              B.RowDist() != C.RowDist() )
              LogicError
              ("Tried to form C[",C.ColDist(),",",C.RowDist(),"] := "
               "A[",A.ColDist(),",",A.RowDist(),"]' "
               "B[",B.ColDist(),",",B.RowDist(),"]");
          if( A.RowAlign() != C.ColAlign() )
              LogicError("A's rows must align with C's cols");
          if( A.ColAlign() != B.ColAlign() )
              LogicError("A's cols must align with B's cols");
          if( B.RowAlign() != C.RowAlign() )
              LogicError("B's rows must align with C's rows");
          if( A.Width() != C.Height() ||
              A.Height() != B.Height() ||
              B.Width() != C.Width() )
              LogicError
              ("Nonconformal LocalGemmTN:\n",
               DimsString(A,"A"),"\n",
               DimsString(B,"B"),"\n",
               DimsString(C,"C"));
      }
      else
      {
          if( A.RowDist() != C.ColDist() ||
              A.ColDist() != B.RowDist() ||
              B.ColDist() != C.RowDist() )
              LogicError
              ("Tried to form C[",C.ColDist(),",",C.RowDist(),"] := "
               "A[",A.ColDist(),",",A.RowDist(),"]' "
               "B[",B.ColDist(),",",B.RowDist(),"]'");
          if( A.RowAlign() != C.ColAlign() )
              LogicError("A's rows must align with C's cols");
          if( A.ColAlign() != B.RowAlign() )
              LogicError("A's cols must align with B's rows");
          if( B.ColAlign() != C.RowAlign() )
              LogicError("B's cols must align with C's rows");
          if( A.Width() != C.Height() ||
              A.Height() != B.Width() ||
              B.Height() != C.Width() )
              LogicError
              ("Nonconformal LocalGemmTT:\n",
               DimsString(A,"A"),"\n",
               DimsString(B,"B"),"\n",
               DimsString(C,"C"));
      }
    )
    Gemm
    ( orientA , orientB,
      alpha, A.LockedMatrix(), B.LockedMatrix(), beta, C.Matrix() );
}

template<typename T>
void LocalGemm
( Orientation orientA, Orientation orientB,
  T alpha, const ElementalMatrix<T>& A,
           const ElementalMatrix<T>& B,
                 ElementalMatrix<T>& C )
{
    DEBUG_ONLY(CSE cse("LocalGemm"))
    const Int m = ( orientA==NORMAL ? A.Height() : A.Width() );
    const Int n = ( orientB==NORMAL ? B.Width() : B.Height() );
    Zeros( C, m, n );
    LocalGemm( orientA, orientB, alpha, A, B, T(0), C );
}

#define PROTO(T) \
  template void Gemm \
  ( Orientation orientA, Orientation orientB, \
    T alpha, const Matrix<T>& A, \
             const Matrix<T>& B, \
    T beta,        Matrix<T>& C ); \
  template void Gemm \
  ( Orientation orientA, Orientation orientB, \
    T alpha, const Matrix<T>& A, \
             const Matrix<T>& B, \
                   Matrix<T>& C ); \
  template void Gemm \
  ( Orientation orientA, Orientation orientB, \
    T alpha, const ElementalMatrix<T>& A, \
             const ElementalMatrix<T>& B, \
    T beta,        ElementalMatrix<T>& C, GemmAlgorithm alg ); \
  template void Gemm \
  ( Orientation orientA, Orientation orientB, \
    T alpha, const ElementalMatrix<T>& A, \
             const ElementalMatrix<T>& B, \
                   ElementalMatrix<T>& C, GemmAlgorithm alg ); \
  template void LocalGemm \
  ( Orientation orientA, Orientation orientB, \
    T alpha, const ElementalMatrix<T>& A, \
             const ElementalMatrix<T>& B, \
    T beta,        ElementalMatrix<T>& C ); \
  template void LocalGemm \
  ( Orientation orientA, Orientation orientB, \
    T alpha, const ElementalMatrix<T>& A, \
             const ElementalMatrix<T>& B, \
                   ElementalMatrix<T>& C );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

} // namespace El
