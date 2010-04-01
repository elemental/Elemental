/*
   Copyright 2009-2010 Jack Poulson

   This file is part of Elemental.

   Elemental is free software: you can redistribute it and/or modify it under
   the terms of the GNU Lesser General Public License as published by the
   Free Software Foundation; either version 3 of the License, or 
   (at your option) any later version.

   Elemental is distributed in the hope that it will be useful, but 
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with Elemental. If not, see <http://www.gnu.org/licenses/>.
*/
#include "ElementalBLASInternal.h"
using namespace std;
using namespace Elemental;

template<typename T>
void
Elemental::BLAS::Gemm
( const Orientation orientationOfA, 
  const Orientation orientationOfB,
  const T alpha, const DistMatrix<T,MC,MR>& A,
                 const DistMatrix<T,MC,MR>& B,
  const T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("BLAS::Gemm");
#endif
    if( orientationOfA == Normal && orientationOfB == Normal )
    {
        BLAS::Internal::GemmNN( alpha, A, B, beta, C );
    }
    else if( orientationOfA == Normal )
    {
        BLAS::Internal::GemmNT( orientationOfB, alpha, A, B, beta, C );
    }
    else if( orientationOfB == Normal )
    {
        BLAS::Internal::GemmTN( orientationOfA, alpha, A, B, beta, C );
    }
    else
    {
        BLAS::Internal::GemmTT( orientationOfA, orientationOfB,
                                alpha, A, B, beta, C );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::BLAS::Internal::GemmA
( const Orientation orientationOfA, 
  const Orientation orientationOfB,
  const T alpha, const DistMatrix<T,MC,MR>& A,
                 const DistMatrix<T,MC,MR>& B,
  const T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("BLAS::Internal::GemmA");
#endif
    if( orientationOfA == Normal && orientationOfB == Normal )
    {
        BLAS::Internal::GemmNNA( alpha, A, B, beta, C );
    }
    else if( orientationOfA == Normal )
    {
        BLAS::Internal::GemmNTA( orientationOfB, alpha, A, B, beta, C );
    }
    else if( orientationOfB == Normal )
    {
        BLAS::Internal::GemmTNA( orientationOfA, alpha, A, B, beta, C );
    }
    else
    {
        BLAS::Internal::GemmTTA( orientationOfA, orientationOfB,
                                 alpha, A, B, beta, C );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::BLAS::Internal::GemmB
( const Orientation orientationOfA, 
  const Orientation orientationOfB,
  const T alpha, const DistMatrix<T,MC,MR>& A,
                 const DistMatrix<T,MC,MR>& B,
  const T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("BLAS::Internal::GemmB");
#endif
    if( orientationOfA == Normal && orientationOfB == Normal )
    {
        BLAS::Internal::GemmNNB( alpha, A, B, beta, C );
    }
    else if( orientationOfA == Normal )
    {
        BLAS::Internal::GemmNTB( orientationOfB, alpha, A, B, beta, C );
    }
    else if( orientationOfB == Normal )
    {
        BLAS::Internal::GemmTNB( orientationOfA, alpha, A, B, beta, C );
    }
    else
    {
        BLAS::Internal::GemmTTB( orientationOfA, orientationOfB, 
                                 alpha, A, B, beta, C );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::BLAS::Internal::GemmC
( const Orientation orientationOfA, 
  const Orientation orientationOfB,
  const T alpha, const DistMatrix<T,MC,MR>& A,
                 const DistMatrix<T,MC,MR>& B,
  const T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("BLAS::Internal::GemmC");
#endif
    if( orientationOfA == Normal && orientationOfB == Normal )
    {
        BLAS::Internal::GemmNNC( alpha, A, B, beta, C );
    }
    else if( orientationOfA == Normal )
    {
        BLAS::Internal::GemmNTC( orientationOfB, alpha, A, B, beta, C );
    }
    else if( orientationOfB == Normal )
    {
        BLAS::Internal::GemmTNC( orientationOfA, alpha, A, B, beta, C );
    }
    else
    {
        BLAS::Internal::GemmTTC( orientationOfA, orientationOfB,
                                 alpha, A, B, beta, C );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::BLAS::Internal::GemmDot
( const Orientation orientationOfA, 
  const Orientation orientationOfB,
  const T alpha, const DistMatrix<T,MC,MR>& A,
                 const DistMatrix<T,MC,MR>& B,
  const T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("BLAS::Internal::GemmDot");
#endif
    if( orientationOfA == Normal && orientationOfB == Normal )
        BLAS::Internal::GemmNNDot( alpha, A, B, beta, C );
    else
        throw "GemmDot currently only implemented for NN case.";
    // This code will be enabled when the routines are implemented
    /*
    else if( orientationOfA == Normal )
    {
        BLAS::Internal::GemmNTDot( orientationOfB, alpha, A, B, beta, C );
    }
    else if( orientationOfB == Normal )
    {
        BLAS::Internal::GemmTNDot( orientationOfA, alpha, A, B, beta, C );
    }
    else
    {
        BLAS::Internal::GemmTTDot( orientationOfA, orientationOfB,
                                   alpha, A, B, beta, C );
    }
    */
#ifndef RELEASE
    PopCallStack();
#endif
}

template void Elemental::BLAS::Gemm
( const Orientation orientationOfA, 
  const Orientation orientationOfB,
  const float alpha, const DistMatrix<float,MC,MR>& A,
                     const DistMatrix<float,MC,MR>& B,
  const float beta,        DistMatrix<float,MC,MR>& C );

template void Elemental::BLAS::Internal::GemmA
( const Orientation orientationOfA, 
  const Orientation orientationOfB,
  const float alpha, const DistMatrix<float,MC,MR>& A,
                     const DistMatrix<float,MC,MR>& B,
  const float beta,        DistMatrix<float,MC,MR>& C );

template void Elemental::BLAS::Internal::GemmB
( const Orientation orientationOfA, 
  const Orientation orientationOfB,
  const float alpha, const DistMatrix<float,MC,MR>& A,
                     const DistMatrix<float,MC,MR>& B,
  const float beta,        DistMatrix<float,MC,MR>& C );

template void Elemental::BLAS::Internal::GemmC
( const Orientation orientationOfA, 
  const Orientation orientationOfB,
  const float alpha, const DistMatrix<float,MC,MR>& A,
                     const DistMatrix<float,MC,MR>& B,
  const float beta,        DistMatrix<float,MC,MR>& C );

template void Elemental::BLAS::Internal::GemmDot
( const Orientation orientationOfA, 
  const Orientation orientationOfB,
  const float alpha, const DistMatrix<float,MC,MR>& A,
                     const DistMatrix<float,MC,MR>& B,
  const float beta,        DistMatrix<float,MC,MR>& C );

template void Elemental::BLAS::Gemm
( const Orientation orientationOfA, 
  const Orientation orientationOfB,
  const double alpha, const DistMatrix<double,MC,MR>& A,
                      const DistMatrix<double,MC,MR>& B,
  const double beta,        DistMatrix<double,MC,MR>& C );

template void Elemental::BLAS::Internal::GemmA
( const Orientation orientationOfA, 
  const Orientation orientationOfB,
  const double alpha, const DistMatrix<double,MC,MR>& A,
                      const DistMatrix<double,MC,MR>& B,
  const double beta,        DistMatrix<double,MC,MR>& C );

template void Elemental::BLAS::Internal::GemmB
( const Orientation orientationOfA, 
  const Orientation orientationOfB,
  const double alpha, const DistMatrix<double,MC,MR>& A,
                      const DistMatrix<double,MC,MR>& B,
  const double beta,        DistMatrix<double,MC,MR>& C );

template void Elemental::BLAS::Internal::GemmC
( const Orientation orientationOfA, 
  const Orientation orientationOfB,
  const double alpha, const DistMatrix<double,MC,MR>& A,
                      const DistMatrix<double,MC,MR>& B,
  const double beta,        DistMatrix<double,MC,MR>& C );

template void Elemental::BLAS::Internal::GemmDot
( const Orientation orientationOfA, 
  const Orientation orientationOfB,
  const double alpha, const DistMatrix<double,MC,MR>& A,
                      const DistMatrix<double,MC,MR>& B,
  const double beta,        DistMatrix<double,MC,MR>& C );

#ifndef WITHOUT_COMPLEX
template void Elemental::BLAS::Gemm
( const Orientation orientationOfA, 
  const Orientation orientationOfB,
  const scomplex alpha, const DistMatrix<scomplex,MC,MR>& A,
                        const DistMatrix<scomplex,MC,MR>& B,
  const scomplex beta,        DistMatrix<scomplex,MC,MR>& C );

template void Elemental::BLAS::Internal::GemmA
( const Orientation orientationOfA, 
  const Orientation orientationOfB,
  const scomplex alpha, const DistMatrix<scomplex,MC,MR>& A,
                        const DistMatrix<scomplex,MC,MR>& B,
  const scomplex beta,        DistMatrix<scomplex,MC,MR>& C );

template void Elemental::BLAS::Internal::GemmB
( const Orientation orientationOfA, 
  const Orientation orientationOfB,
  const scomplex alpha, const DistMatrix<scomplex,MC,MR>& A,
                        const DistMatrix<scomplex,MC,MR>& B,
  const scomplex beta,        DistMatrix<scomplex,MC,MR>& C );

template void Elemental::BLAS::Internal::GemmC
( const Orientation orientationOfA, 
  const Orientation orientationOfB,
  const scomplex alpha, const DistMatrix<scomplex,MC,MR>& A,
                        const DistMatrix<scomplex,MC,MR>& B,
  const scomplex beta,        DistMatrix<scomplex,MC,MR>& C );

template void Elemental::BLAS::Internal::GemmDot
( const Orientation orientationOfA,
  const Orientation orientationOfB,
  const scomplex alpha, const DistMatrix<scomplex,MC,MR>& A,
                        const DistMatrix<scomplex,MC,MR>& B,
  const scomplex beta,        DistMatrix<scomplex,MC,MR>& C );

template void Elemental::BLAS::Gemm
( const Orientation orientationOfA, 
  const Orientation orientationOfB,
  const dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& A,
                        const DistMatrix<dcomplex,MC,MR>& B,
  const dcomplex beta,        DistMatrix<dcomplex,MC,MR>& C );

template void Elemental::BLAS::Internal::GemmA
( const Orientation orientationOfA, 
  const Orientation orientationOfB,
  const dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& A,
                        const DistMatrix<dcomplex,MC,MR>& B,
  const dcomplex beta,        DistMatrix<dcomplex,MC,MR>& C );

template void Elemental::BLAS::Internal::GemmB
( const Orientation orientationOfA, 
  const Orientation orientationOfB,
  const dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& A,
                        const DistMatrix<dcomplex,MC,MR>& B,
  const dcomplex beta,        DistMatrix<dcomplex,MC,MR>& C );

template void Elemental::BLAS::Internal::GemmC
( const Orientation orientationOfA, 
  const Orientation orientationOfB,
  const dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& A,
                        const DistMatrix<dcomplex,MC,MR>& B,
  const dcomplex beta,        DistMatrix<dcomplex,MC,MR>& C );

template void Elemental::BLAS::Internal::GemmDot
( const Orientation orientationOfA,
  const Orientation orientationOfB,
  const dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& A,
                        const DistMatrix<dcomplex,MC,MR>& B,
  const dcomplex beta,        DistMatrix<dcomplex,MC,MR>& C );
#endif
