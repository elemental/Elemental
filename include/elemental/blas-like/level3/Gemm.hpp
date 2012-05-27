/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

    - Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

    - Neither the name of the owner nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/

#include "./Gemm/NN.hpp"
#include "./Gemm/NT.hpp"
#include "./Gemm/TN.hpp"
#include "./Gemm/TT.hpp"

namespace elem {

template<typename T>
inline void
Gemm
( Orientation orientationOfA, 
  Orientation orientationOfB,
  T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
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

template<typename T>
inline void
internal::GemmA
( Orientation orientationOfA, 
  Orientation orientationOfB,
  T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("internal::GemmA");
#endif
    if( orientationOfA == NORMAL && orientationOfB == NORMAL )
    {
        internal::GemmNNA( alpha, A, B, beta, C );
    }
    else if( orientationOfA == NORMAL )
    {
        internal::GemmNTA( orientationOfB, alpha, A, B, beta, C );
    }
    else if( orientationOfB == NORMAL )
    {
        internal::GemmTNA( orientationOfA, alpha, A, B, beta, C );
    }
    else
    {
        internal::GemmTTA
        ( orientationOfA, orientationOfB, alpha, A, B, beta, C );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
internal::GemmB
( Orientation orientationOfA, 
  Orientation orientationOfB,
  T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("internal::GemmB");
#endif
    if( orientationOfA == NORMAL && orientationOfB == NORMAL )
    {
        internal::GemmNNB( alpha, A, B, beta, C );
    }
    else if( orientationOfA == NORMAL )
    {
        internal::GemmNTB( orientationOfB, alpha, A, B, beta, C );
    }
    else if( orientationOfB == NORMAL )
    {
        internal::GemmTNB( orientationOfA, alpha, A, B, beta, C );
    }
    else
    {
        internal::GemmTTB
        ( orientationOfA, orientationOfB, alpha, A, B, beta, C );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
internal::GemmC
( Orientation orientationOfA, 
  Orientation orientationOfB,
  T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("internal::GemmC");
#endif
    if( orientationOfA == NORMAL && orientationOfB == NORMAL )
    {
        internal::GemmNNC( alpha, A, B, beta, C );
    }
    else if( orientationOfA == NORMAL )
    {
        internal::GemmNTC( orientationOfB, alpha, A, B, beta, C );
    }
    else if( orientationOfB == NORMAL )
    {
        internal::GemmTNC( orientationOfA, alpha, A, B, beta, C );
    }
    else
    {
        internal::GemmTTC
        ( orientationOfA, orientationOfB, alpha, A, B, beta, C );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
internal::GemmDot
( Orientation orientationOfA, 
  Orientation orientationOfB,
  T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("internal::GemmDot");
#endif
    if( orientationOfA == NORMAL && orientationOfB == NORMAL )
        internal::GemmNNDot( alpha, A, B, beta, C );
    else
        throw std::logic_error("GemmDot only implemented for NN case");
    // This code will be enabled when the routines are implemented
    /*
    else if( orientationOfA == NORMAL )
    {
        internal::GemmNTDot( orientationOfB, alpha, A, B, beta, C );
    }
    else if( orientationOfB == NORMAL )
    {
        internal::GemmTNDot( orientationOfA, alpha, A, B, beta, C );
    }
    else
    {
        internal::GemmTTDot( orientationOfA, orientationOfB,
                                   alpha, A, B, beta, C );
    }
    */
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
