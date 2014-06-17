/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"



#include "./Gemm/NN.hpp"
#include "./Gemm/NT.hpp"
#include "./Gemm/TN.hpp"
#include "./Gemm/TT.hpp"

namespace El {

template<typename T>
void Gemm
( Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const Matrix<T>& A, const Matrix<T>& B, T beta, Matrix<T>& C )
{
    DEBUG_ONLY(
        CallStackEntry cse("Gemm");
        if( orientationOfA == NORMAL && orientationOfB == NORMAL )
        {
            if( A.Height() != C.Height() ||
                B.Width()  != C.Width()  ||
                A.Width()  != B.Height() )
                LogicError("Nonconformal GemmNN");
        }
        else if( orientationOfA == NORMAL )
        {
            if( A.Height() != C.Height() ||
                B.Height() != C.Width()  ||
                A.Width()  != B.Width() )
                LogicError("Nonconformal GemmN(T/C)");
        }
        else if( orientationOfB == NORMAL )
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
    const char transA = OrientationToChar( orientationOfA );
    const char transB = OrientationToChar( orientationOfB );
    const Int m = C.Height();
    const Int n = C.Width();
    const Int k = ( orientationOfA == NORMAL ? A.Width() : A.Height() );
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
( Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const Matrix<T>& A, const Matrix<T>& B, Matrix<T>& C )
{
    DEBUG_ONLY(CallStackEntry cse("Gemm"))
    const Int m = ( orientationOfA==NORMAL ? A.Height() : A.Width() );
    const Int n = ( orientationOfB==NORMAL ? B.Width() : B.Height() );
    Zeros( C, m, n );
    Gemm( orientationOfA, orientationOfB, alpha, A, B, T(0), C );
}

template<typename T>
void Gemm
( Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const DistMatrix<T>& A, const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C, GemmAlgorithm alg )
{
    DEBUG_ONLY(CallStackEntry cse("Gemm"))
    if( orientationOfA == NORMAL && orientationOfB == NORMAL )
    {
        if( alg == GEMM_CANNON )
            gemm::Cannon_NN( alpha, A, B, beta, C );
        else 
            gemm::SUMMA_NN( alpha, A, B, beta, C, alg );
    }
    else if( orientationOfA == NORMAL )
    {
        gemm::SUMMA_NT( orientationOfB, alpha, A, B, beta, C, alg );
    }
    else if( orientationOfB == NORMAL )
    {
        gemm::SUMMA_TN( orientationOfA, alpha, A, B, beta, C, alg );
    }
    else
    {
        gemm::SUMMA_TT
        ( orientationOfA, orientationOfB, alpha, A, B, beta, C, alg );
    }
}

template<typename T>
void Gemm
( Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const DistMatrix<T>& A, const DistMatrix<T>& B,
                 DistMatrix<T>& C, GemmAlgorithm alg )
{
    DEBUG_ONLY(CallStackEntry cse("Gemm"))
    const Int m = ( orientationOfA==NORMAL ? A.Height() : A.Width() );
    const Int n = ( orientationOfB==NORMAL ? B.Width() : B.Height() );
    Zeros( C, m, n );
    Gemm( orientationOfA, orientationOfB, alpha, A, B, T(0), C, alg );
}

#define PROTO(T) \
  template void Gemm \
  ( Orientation orientationA, Orientation orientationB, \
    T alpha, const Matrix<T>& A, const Matrix<T>& B, \
    T beta,        Matrix<T>& C ); \
  template void Gemm \
  ( Orientation orientationA, Orientation orientationB, \
    T alpha, const Matrix<T>& A, const Matrix<T>& B, \
                   Matrix<T>& C ); \
  template void Gemm \
  ( Orientation orientationA, Orientation orientationB, \
    T alpha, const DistMatrix<T>& A, const DistMatrix<T>& B, \
    T beta,        DistMatrix<T>& C, GemmAlgorithm alg ); \
  template void Gemm \
  ( Orientation orientationA, Orientation orientationB, \
    T alpha, const DistMatrix<T>& A, const DistMatrix<T>& B, \
                   DistMatrix<T>& C, GemmAlgorithm alg );

PROTO(Int)
PROTO(float)
PROTO(double)
PROTO(Complex<float>)
PROTO(Complex<double>)

} // namespace El
