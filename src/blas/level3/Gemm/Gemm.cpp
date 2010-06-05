/*
   This file is part of elemental, a library for distributed-memory dense 
   linear algebra.

   Copyright (C) 2009-2010 Jack Poulson <jack.poulson@gmail.com>

   This program is released under the terms of the license contained in the 
   file LICENSE.
*/
#include "elemental/blas_internal.hpp"
using namespace std;
using namespace elemental;

template<typename T>
void
elemental::blas::Gemm
( Orientation orientationOfA, 
  Orientation orientationOfB,
  T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("blas::Gemm");
#endif
    if( orientationOfA == Normal && orientationOfB == Normal )
    {
        blas::internal::GemmNN( alpha, A, B, beta, C );
    }
    else if( orientationOfA == Normal )
    {
        blas::internal::GemmNT( orientationOfB, alpha, A, B, beta, C );
    }
    else if( orientationOfB == Normal )
    {
        blas::internal::GemmTN
        ( orientationOfA, alpha, A, B, beta, C );
    }
    else
    {
        blas::internal::GemmTT
        ( orientationOfA, orientationOfB, alpha, A, B, beta, C );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::blas::internal::GemmA
( Orientation orientationOfA, 
  Orientation orientationOfB,
  T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("blas::internal::GemmA");
#endif
    if( orientationOfA == Normal && orientationOfB == Normal )
    {
        blas::internal::GemmNNA( alpha, A, B, beta, C );
    }
    else if( orientationOfA == Normal )
    {
        blas::internal::GemmNTA( orientationOfB, alpha, A, B, beta, C );
    }
    else if( orientationOfB == Normal )
    {
        blas::internal::GemmTNA( orientationOfA, alpha, A, B, beta, C );
    }
    else
    {
        blas::internal::GemmTTA
        ( orientationOfA, orientationOfB, alpha, A, B, beta, C );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::blas::internal::GemmB
( Orientation orientationOfA, 
  Orientation orientationOfB,
  T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("blas::internal::GemmB");
#endif
    if( orientationOfA == Normal && orientationOfB == Normal )
    {
        blas::internal::GemmNNB( alpha, A, B, beta, C );
    }
    else if( orientationOfA == Normal )
    {
        blas::internal::GemmNTB( orientationOfB, alpha, A, B, beta, C );
    }
    else if( orientationOfB == Normal )
    {
        blas::internal::GemmTNB( orientationOfA, alpha, A, B, beta, C );
    }
    else
    {
        blas::internal::GemmTTB
        ( orientationOfA, orientationOfB, alpha, A, B, beta, C );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::blas::internal::GemmC
( Orientation orientationOfA, 
  Orientation orientationOfB,
  T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("blas::internal::GemmC");
#endif
    if( orientationOfA == Normal && orientationOfB == Normal )
    {
        blas::internal::GemmNNC( alpha, A, B, beta, C );
    }
    else if( orientationOfA == Normal )
    {
        blas::internal::GemmNTC( orientationOfB, alpha, A, B, beta, C );
    }
    else if( orientationOfB == Normal )
    {
        blas::internal::GemmTNC( orientationOfA, alpha, A, B, beta, C );
    }
    else
    {
        blas::internal::GemmTTC
        ( orientationOfA, orientationOfB, alpha, A, B, beta, C );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::blas::internal::GemmDot
( Orientation orientationOfA, 
  Orientation orientationOfB,
  T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("blas::internal::GemmDot");
#endif
    if( orientationOfA == Normal && orientationOfB == Normal )
        blas::internal::GemmNNDot( alpha, A, B, beta, C );
    else
        throw "GemmDot currently only implemented for NN case.";
    // This code will be enabled when the routines are implemented
    /*
    else if( orientationOfA == Normal )
    {
        blas::internal::GemmNTDot( orientationOfB, alpha, A, B, beta, C );
    }
    else if( orientationOfB == Normal )
    {
        blas::internal::GemmTNDot( orientationOfA, alpha, A, B, beta, C );
    }
    else
    {
        blas::internal::GemmTTDot( orientationOfA, orientationOfB,
                                   alpha, A, B, beta, C );
    }
    */
#ifndef RELEASE
    PopCallStack();
#endif
}

template void elemental::blas::Gemm
( Orientation orientationOfA, 
  Orientation orientationOfB,
  float alpha, const DistMatrix<float,MC,MR>& A,
               const DistMatrix<float,MC,MR>& B,
  float beta,        DistMatrix<float,MC,MR>& C );

template void elemental::blas::internal::GemmA
( Orientation orientationOfA, 
  Orientation orientationOfB,
  float alpha, const DistMatrix<float,MC,MR>& A,
               const DistMatrix<float,MC,MR>& B,
  float beta,        DistMatrix<float,MC,MR>& C );

template void elemental::blas::internal::GemmB
( Orientation orientationOfA, 
  Orientation orientationOfB,
  float alpha, const DistMatrix<float,MC,MR>& A,
               const DistMatrix<float,MC,MR>& B,
  float beta,        DistMatrix<float,MC,MR>& C );

template void elemental::blas::internal::GemmC
( Orientation orientationOfA, 
  Orientation orientationOfB,
  float alpha, const DistMatrix<float,MC,MR>& A,
               const DistMatrix<float,MC,MR>& B,
  float beta,        DistMatrix<float,MC,MR>& C );

template void elemental::blas::internal::GemmDot
( Orientation orientationOfA, 
  Orientation orientationOfB,
  float alpha, const DistMatrix<float,MC,MR>& A,
               const DistMatrix<float,MC,MR>& B,
  float beta,        DistMatrix<float,MC,MR>& C );

template void elemental::blas::Gemm
( Orientation orientationOfA, 
  Orientation orientationOfB,
  double alpha, const DistMatrix<double,MC,MR>& A,
                const DistMatrix<double,MC,MR>& B,
  double beta,        DistMatrix<double,MC,MR>& C );

template void elemental::blas::internal::GemmA
( Orientation orientationOfA, 
  Orientation orientationOfB,
  double alpha, const DistMatrix<double,MC,MR>& A,
                const DistMatrix<double,MC,MR>& B,
  double beta,        DistMatrix<double,MC,MR>& C );

template void elemental::blas::internal::GemmB
( Orientation orientationOfA, 
  Orientation orientationOfB,
  double alpha, const DistMatrix<double,MC,MR>& A,
                const DistMatrix<double,MC,MR>& B,
  double beta,        DistMatrix<double,MC,MR>& C );

template void elemental::blas::internal::GemmC
( Orientation orientationOfA, 
  Orientation orientationOfB,
  double alpha, const DistMatrix<double,MC,MR>& A,
                const DistMatrix<double,MC,MR>& B,
  double beta,        DistMatrix<double,MC,MR>& C );

template void elemental::blas::internal::GemmDot
( Orientation orientationOfA, 
  Orientation orientationOfB,
  double alpha, const DistMatrix<double,MC,MR>& A,
                const DistMatrix<double,MC,MR>& B,
  double beta,        DistMatrix<double,MC,MR>& C );

#ifndef WITHOUT_COMPLEX
template void elemental::blas::Gemm
( Orientation orientationOfA, 
  Orientation orientationOfB,
  scomplex alpha, const DistMatrix<scomplex,MC,MR>& A,
                  const DistMatrix<scomplex,MC,MR>& B,
  scomplex beta,        DistMatrix<scomplex,MC,MR>& C );

template void elemental::blas::internal::GemmA
( Orientation orientationOfA, 
  Orientation orientationOfB,
  scomplex alpha, const DistMatrix<scomplex,MC,MR>& A,
                  const DistMatrix<scomplex,MC,MR>& B,
  scomplex beta,        DistMatrix<scomplex,MC,MR>& C );

template void elemental::blas::internal::GemmB
( Orientation orientationOfA, 
  Orientation orientationOfB,
  scomplex alpha, const DistMatrix<scomplex,MC,MR>& A,
                  const DistMatrix<scomplex,MC,MR>& B,
  scomplex beta,        DistMatrix<scomplex,MC,MR>& C );

template void elemental::blas::internal::GemmC
( Orientation orientationOfA, 
  Orientation orientationOfB,
  scomplex alpha, const DistMatrix<scomplex,MC,MR>& A,
                  const DistMatrix<scomplex,MC,MR>& B,
  scomplex beta,        DistMatrix<scomplex,MC,MR>& C );

template void elemental::blas::internal::GemmDot
( Orientation orientationOfA,
  Orientation orientationOfB,
  scomplex alpha, const DistMatrix<scomplex,MC,MR>& A,
                  const DistMatrix<scomplex,MC,MR>& B,
  scomplex beta,        DistMatrix<scomplex,MC,MR>& C );

template void elemental::blas::Gemm
( Orientation orientationOfA, 
  Orientation orientationOfB,
  dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& A,
                  const DistMatrix<dcomplex,MC,MR>& B,
  dcomplex beta,        DistMatrix<dcomplex,MC,MR>& C );

template void elemental::blas::internal::GemmA
( Orientation orientationOfA, 
  Orientation orientationOfB,
  dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& A,
                  const DistMatrix<dcomplex,MC,MR>& B,
  dcomplex beta,        DistMatrix<dcomplex,MC,MR>& C );

template void elemental::blas::internal::GemmB
( Orientation orientationOfA, 
  Orientation orientationOfB,
  dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& A,
                  const DistMatrix<dcomplex,MC,MR>& B,
  dcomplex beta,        DistMatrix<dcomplex,MC,MR>& C );

template void elemental::blas::internal::GemmC
( Orientation orientationOfA, 
  Orientation orientationOfB,
  dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& A,
                  const DistMatrix<dcomplex,MC,MR>& B,
  dcomplex beta,        DistMatrix<dcomplex,MC,MR>& C );

template void elemental::blas::internal::GemmDot
( Orientation orientationOfA,
  Orientation orientationOfB,
  dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& A,
                  const DistMatrix<dcomplex,MC,MR>& B,
  dcomplex beta,        DistMatrix<dcomplex,MC,MR>& C );
#endif
