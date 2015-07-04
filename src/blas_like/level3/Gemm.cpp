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

template<typename Ring>
void Gemm
( Ring alpha, const OrientedMatrix<Ring>& A, 
              const OrientedMatrix<Ring>& B, 
  Ring beta,                Matrix<Ring>& C )
{
    DEBUG_ONLY(CSE cse("Gemm"))
    if( A.orient == NORMAL && B.orient == NORMAL )
    {
        if( A.matrix.Height() != C.Height() ||
            B.matrix.Width()  != C.Width()  ||
            A.matrix.Width()  != B.matrix.Height() )
            LogicError("Nonconformal GemmNN");
    }
    else if( A.orient == NORMAL )
    {
        if( A.matrix.Height() != C.Height() ||
            B.matrix.Height() != C.Width()  ||
            A.matrix.Width()  != B.matrix.Width() )
            LogicError("Nonconformal GemmN(T/C)");
    }
    else if( B.orient == NORMAL )
    {
        if( A.matrix.Width()  != C.Height() ||
            B.matrix.Width()  != C.Width()  ||
            A.matrix.Height() != B.matrix.Height() )
            LogicError("Nonconformal Gemm(T/C)N");
    }
    else
    {
        if( A.matrix.Width()  != C.Height() ||
            B.matrix.Height() != C.Width()  ||
            A.matrix.Height() != B.matrix.Width() )
            LogicError("Nonconformal Gemm(T/C)(T/C)");
    }
    const char transA = OrientationToChar( A.orient );
    const char transB = OrientationToChar( B.orient );
    const Int m = C.Height();
    const Int n = C.Width();
    const Int k = ( A.orient == NORMAL ? A.matrix.Width() : A.matrix.Height() );
    if( k != 0 )
    {
        blas::Gemm
        ( transA, transB, m, n, k,
          alpha, A.matrix.LockedBuffer(), A.matrix.LDim(), 
                 B.matrix.LockedBuffer(), B.matrix.LDim(),
          beta,  C.Buffer(),       C.LDim() );
    }
    else
    {
        C *= beta;
    }
}

template<typename Ring>
void Gemm
( Ring alpha, const OrientedMatrix<Ring>& A, 
              const OrientedMatrix<Ring>& B, 
                            Matrix<Ring>& C )
{
    DEBUG_ONLY(CSE cse("Gemm"))
    const Int m = ( A.orient==NORMAL ? A.matrix.Height() : A.matrix.Width() );
    const Int n = ( B.orient==NORMAL ? B.matrix.Width() : B.matrix.Height() );
    Zeros( C, m, n );
    Gemm( alpha, A, B, Ring(0), C );
}

template<typename Ring>
void Gemm
( Ring alpha, const OrientedAbstractDistMatrix<Ring>& A, 
              const OrientedAbstractDistMatrix<Ring>& B,
  Ring beta,                AbstractDistMatrix<Ring>& C, 
  GemmAlgorithm alg )
{
    DEBUG_ONLY(CSE cse("Gemm"))
    C *= beta;
    if( A.orient == NORMAL && B.orient == NORMAL )
    {
        if( alg == GEMM_CANNON )
            gemm::Cannon_NN( alpha, A.matrix, B.matrix, C );
        else 
            gemm::SUMMA_NN( alpha, A.matrix, B.matrix, C, alg );
    }
    else if( A.orient == NORMAL )
    {
        gemm::SUMMA_NT( B.orient, alpha, A.matrix, B.matrix, C, alg );
    }
    else if( B.orient == NORMAL )
    {
        gemm::SUMMA_TN( A.orient, alpha, A.matrix, B.matrix, C, alg );
    }
    else
    {
        gemm::SUMMA_TT( A.orient, B.orient, alpha, A.matrix, B.matrix, C, alg );
    }
}

template<typename Ring>
void Gemm
( Ring alpha, const OrientedAbstractDistMatrix<Ring>& A, 
              const OrientedAbstractDistMatrix<Ring>& B,
                            AbstractDistMatrix<Ring>& C, 
  GemmAlgorithm alg )
{
    DEBUG_ONLY(CSE cse("Gemm"))
    const Int m = ( A.orient==NORMAL ? A.matrix.Height() : A.matrix.Width() );
    const Int n = ( B.orient==NORMAL ? B.matrix.Width() : B.matrix.Height() );
    Zeros( C, m, n );
    Gemm( alpha, A, B, Ring(0), C, alg );
}

template<typename Ring>
void LocalGemm
( Ring alpha, const OrientedAbstractDistMatrix<Ring>& AOr,
              const OrientedAbstractDistMatrix<Ring>& BOr,
  Ring beta,                AbstractDistMatrix<Ring>& C )
{
    DEBUG_ONLY(CSE cse("LocalGemm"))
    const auto orientA = AOr.orient;
    const auto orientB = BOr.orient;
    const auto& A = AOr.matrix;
    const auto& B = BOr.matrix;
    DEBUG_ONLY(
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
    ( alpha, A.LockedMatrix().Orient(orientA),
             B.LockedMatrix().Orient(orientB), 
      beta,  C.Matrix() );
}

template<typename Ring>
void LocalGemm
( Ring alpha, const OrientedAbstractDistMatrix<Ring>& A,
              const OrientedAbstractDistMatrix<Ring>& B,
                            AbstractDistMatrix<Ring>& C )
{
    DEBUG_ONLY(CSE cse("LocalGemm"))
    const Int m = ( A.orient==NORMAL ? A.matrix.Height() : A.matrix.Width() );
    const Int n = ( B.orient==NORMAL ? B.matrix.Width() : B.matrix.Height() );
    Zeros( C, m, n );
    LocalGemm( alpha, A, B, Ring(0), C );
}

// Deprecated
// ==========
template<typename Ring>
void Gemm
( Orientation orientA, Orientation orientB,
  Ring alpha, const AbstractDistMatrix<Ring>& A, 
              const AbstractDistMatrix<Ring>& B,
  Ring beta,        AbstractDistMatrix<Ring>& C, 
  GemmAlgorithm alg )
{
    DEBUG_ONLY(CSE cse("Gemm"))
    Gemm( alpha, A.Orient(orientA), B.Orient(orientB), beta, C, alg );
}

template<typename Ring>
void Gemm
( Orientation orientA, Orientation orientB,
  Ring alpha, const Matrix<Ring>& A, 
              const Matrix<Ring>& B, 
                    Matrix<Ring>& C )
{
    DEBUG_ONLY(CSE cse("Gemm"))
    Gemm( alpha, A.Orient(orientA), B.Orient(orientB), C );
}

template<typename Ring>
void Gemm
( Orientation orientA, Orientation orientB,
  Ring alpha, const Matrix<Ring>& A, const Matrix<Ring>& B, 
  Ring beta,        Matrix<Ring>& C )
{
    DEBUG_ONLY(CSE cse("Gemm"))
    Gemm( alpha, A.Orient(orientA), B.Orient(orientB), beta, C );
}


template<typename Ring>
void Gemm
( Orientation orientA, Orientation orientB,
  Ring alpha, const AbstractDistMatrix<Ring>& A, 
              const AbstractDistMatrix<Ring>& B,
                    AbstractDistMatrix<Ring>& C, 
  GemmAlgorithm alg )
{
    DEBUG_ONLY(CSE cse("Gemm"))
    Gemm( alpha, A.Orient(orientA), B.Orient(orientB), C, alg );
}

template<typename Ring>
void LocalGemm
( Orientation orientA, Orientation orientB,
  Ring alpha, const AbstractDistMatrix<Ring>& A,
              const AbstractDistMatrix<Ring>& B,
  Ring beta,        AbstractDistMatrix<Ring>& C )
{
    DEBUG_ONLY(CSE cse("LocalGemm"))
    LocalGemm( alpha, A.Orient(orientA), B.Orient(orientB), beta, C );
}


template<typename Ring>
void LocalGemm
( Orientation orientA, Orientation orientB,
  Ring alpha, const AbstractDistMatrix<Ring>& A,
              const AbstractDistMatrix<Ring>& B,
                    AbstractDistMatrix<Ring>& C )
{
    DEBUG_ONLY(CSE cse("LocalGemm"))
    LocalGemm( alpha, A.Orient(orientA), B.Orient(orientB), C );
}

#define PROTO(Ring) \
  template void Gemm \
  ( Ring alpha, const OrientedMatrix<Ring>& A, \
                const OrientedMatrix<Ring>& B, \
    Ring beta,                Matrix<Ring>& C ); \
  template void Gemm \
  ( Ring alpha, const OrientedMatrix<Ring>& A, \
                const OrientedMatrix<Ring>& B, \
                              Matrix<Ring>& C ); \
  template void Gemm \
  ( Ring alpha, const OrientedAbstractDistMatrix<Ring>& A, \
                const OrientedAbstractDistMatrix<Ring>& B, \
    Ring beta,                AbstractDistMatrix<Ring>& C, \
    GemmAlgorithm alg ); \
  template void Gemm \
  ( Ring alpha, const OrientedAbstractDistMatrix<Ring>& A, \
                const OrientedAbstractDistMatrix<Ring>& B, \
                              AbstractDistMatrix<Ring>& C, \
    GemmAlgorithm alg ); \
  template void LocalGemm \
  ( Ring alpha, const OrientedAbstractDistMatrix<Ring>& A, \
                const OrientedAbstractDistMatrix<Ring>& B, \
    Ring beta,                AbstractDistMatrix<Ring>& C ); \
  template void LocalGemm \
  ( Ring alpha, const OrientedAbstractDistMatrix<Ring>& A, \
                const OrientedAbstractDistMatrix<Ring>& B, \
                              AbstractDistMatrix<Ring>& C ); \
  /* Deprecated */ \
  template void Gemm \
  ( Orientation orientA, Orientation orientB, \
    Ring alpha, const Matrix<Ring>& A, \
                const Matrix<Ring>& B, \
    Ring beta,        Matrix<Ring>& C ); \
  template void Gemm \
  ( Orientation orientA, Orientation orientB, \
    Ring alpha, const Matrix<Ring>& A, \
                const Matrix<Ring>& B, \
                      Matrix<Ring>& C ); \
  template void Gemm \
  ( Orientation orientA, Orientation orientB, \
    Ring alpha, const AbstractDistMatrix<Ring>& A, \
                const AbstractDistMatrix<Ring>& B, \
    Ring beta,        AbstractDistMatrix<Ring>& C, \
    GemmAlgorithm alg ); \
  template void Gemm \
  ( Orientation orientA, Orientation orientB, \
    Ring alpha, const AbstractDistMatrix<Ring>& A, \
                const AbstractDistMatrix<Ring>& B, \
                      AbstractDistMatrix<Ring>& C, \
    GemmAlgorithm alg ); \
  template void LocalGemm \
  ( Orientation orientA, Orientation orientB, \
    Ring alpha, const AbstractDistMatrix<Ring>& A, \
                const AbstractDistMatrix<Ring>& B, \
    Ring beta,        AbstractDistMatrix<Ring>& C ); \
  template void LocalGemm \
  ( Orientation orientA, Orientation orientB, \
    Ring alpha, const AbstractDistMatrix<Ring>& A, \
                const AbstractDistMatrix<Ring>& B, \
                      AbstractDistMatrix<Ring>& C );

#define EL_ENABLE_QUAD
#include "El/macros/Instantiate.h"

} // namespace El
