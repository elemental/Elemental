/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

#include "./Gemm/NN.hpp"
#include "./Gemm/NT.hpp"
#include "./Gemm/TN.hpp"
#include "./Gemm/TT.hpp"

namespace elem {

template<typename T>
inline void
Gemm
( Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const Matrix<T>& A, const Matrix<T>& B, T beta, Matrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("Gemm");
    if( orientationOfA == NORMAL && orientationOfB == NORMAL )
    {
        if( A.Height() != C.Height() ||
            B.Width()  != C.Width()  ||
            A.Width()  != B.Height() )
            throw std::logic_error("Nonconformal GemmNN");
    }
    else if( orientationOfA == NORMAL )
    {
        if( A.Height() != C.Height() ||
            B.Height() != C.Width()  ||
            A.Width()  != B.Width() )
            throw std::logic_error("Nonconformal GemmN(T/C)");
    }
    else if( orientationOfB == NORMAL )
    {
        if( A.Width()  != C.Height() ||
            B.Width()  != C.Width()  ||
            A.Height() != B.Height() )
            throw std::logic_error("Nonconformal Gemm(T/C)N");
    }
    else
    {
        if( A.Width()  != C.Height() ||
            B.Height() != C.Width()  ||
            A.Height() != B.Width() )
            throw std::logic_error("Nonconformal Gemm(T/C)(T/C)");
    }
#endif
    const char transA = OrientationToChar( orientationOfA );
    const char transB = OrientationToChar( orientationOfB );
    const int m = C.Height();
    const int n = C.Width();
    const int k = ( orientationOfA == NORMAL ? A.Width() : A.Height() );
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
#ifndef RELEASE
    PopCallStack();
#endif
}

namespace internal {

template<typename T>
inline void
GemmA
( Orientation orientationOfA, 
  Orientation orientationOfB,
  T alpha, const DistMatrix<T>& A,
           const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("internal::GemmA");
#endif
    if( orientationOfA == NORMAL && orientationOfB == NORMAL )
    {
        GemmNNA( alpha, A, B, beta, C );
    }
    else if( orientationOfA == NORMAL )
    {
        GemmNTA( orientationOfB, alpha, A, B, beta, C );
    }
    else if( orientationOfB == NORMAL )
    {
        GemmTNA( orientationOfA, alpha, A, B, beta, C );
    }
    else
    {
        GemmTTA( orientationOfA, orientationOfB, alpha, A, B, beta, C );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
GemmB
( Orientation orientationOfA, 
  Orientation orientationOfB,
  T alpha, const DistMatrix<T>& A,
           const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("internal::GemmB");
#endif
    if( orientationOfA == NORMAL && orientationOfB == NORMAL )
    {
        GemmNNB( alpha, A, B, beta, C );
    }
    else if( orientationOfA == NORMAL )
    {
        GemmNTB( orientationOfB, alpha, A, B, beta, C );
    }
    else if( orientationOfB == NORMAL )
    {
        GemmTNB( orientationOfA, alpha, A, B, beta, C );
    }
    else
    {
        GemmTTB( orientationOfA, orientationOfB, alpha, A, B, beta, C );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
GemmC
( Orientation orientationOfA, 
  Orientation orientationOfB,
  T alpha, const DistMatrix<T>& A,
           const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("internal::GemmC");
#endif
    if( orientationOfA == NORMAL && orientationOfB == NORMAL )
    {
        GemmNNC( alpha, A, B, beta, C );
    }
    else if( orientationOfA == NORMAL )
    {
        GemmNTC( orientationOfB, alpha, A, B, beta, C );
    }
    else if( orientationOfB == NORMAL )
    {
        GemmTNC( orientationOfA, alpha, A, B, beta, C );
    }
    else
    {
        GemmTTC( orientationOfA, orientationOfB, alpha, A, B, beta, C );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
GemmDot
( Orientation orientationOfA, 
  Orientation orientationOfB,
  T alpha, const DistMatrix<T>& A,
           const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("internal::GemmDot");
#endif
    if( orientationOfA == NORMAL && orientationOfB == NORMAL )
        GemmNNDot( alpha, A, B, beta, C );
    else
        throw std::logic_error("GemmDot only implemented for NN case");
    // This code will be enabled when the routines are implemented
    /*
    else if( orientationOfA == NORMAL )
    {
        GemmNTDot( orientationOfB, alpha, A, B, beta, C );
    }
    else if( orientationOfB == NORMAL )
    {
        GemmTNDot( orientationOfA, alpha, A, B, beta, C );
    }
    else
    {
        GemmTTDot( orientationOfA, orientationOfB,
                                   alpha, A, B, beta, C );
    }
    */
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace internal

template<typename T>
inline void
Gemm
( Orientation orientationOfA, 
  Orientation orientationOfB,
  T alpha, const DistMatrix<T>& A,
           const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("Gemm");
#endif
    if( orientationOfA == NORMAL && orientationOfB == NORMAL )
    {
        internal::GemmNN( alpha, A, B, beta, C );
    }
    else if( orientationOfA == NORMAL )
    {
        internal::GemmNT( orientationOfB, alpha, A, B, beta, C );
    }
    else if( orientationOfB == NORMAL )
    {
        internal::GemmTN( orientationOfA, alpha, A, B, beta, C );
    }
    else
    {
        internal::GemmTT
        ( orientationOfA, orientationOfB, alpha, A, B, beta, C );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
