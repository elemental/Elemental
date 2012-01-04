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

#include "./Trrk/LocalTrrk.hpp"
#include "./Trrk/TrrkNN.hpp"
#include "./Trrk/TrrkNT.hpp"
#include "./Trrk/TrrkTN.hpp"
#include "./Trrk/TrrkTT.hpp"

namespace elemental {

template<typename T>
inline void
Trrk
( UpperOrLower uplo, 
  Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const Matrix<T>& A, const Matrix<T>& B,
  T beta,        Matrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("Trrk");
#endif
    if( orientationOfA==NORMAL && orientationOfB==NORMAL )
        internal::TrrkNN( uplo, alpha, A, B, beta, C );
    else if( orientationOfA==NORMAL )
        internal::TrrkNT( uplo, orientationOfB, alpha, A, B, beta, C );
    else if( orientationOfB==NORMAL )
        internal::TrrkTN( uplo, orientationOfA, alpha, A, B, beta, C );
    else
        internal::TrrkTT
        ( uplo, orientationOfA, orientationOfB, alpha, A, B, beta, C );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Trrk
( UpperOrLower uplo, 
  Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("Trrk");
#endif
    if( orientationOfA==NORMAL && orientationOfB==NORMAL )
        internal::TrrkNN( uplo, alpha, A, B, beta, C );
    else if( orientationOfA==NORMAL )
        internal::TrrkNT( uplo, orientationOfB, alpha, A, B, beta, C );
    else if( orientationOfB==NORMAL )
        internal::TrrkTN( uplo, orientationOfA, alpha, A, B, beta, C );
    else
        internal::TrrkTT
        ( uplo, orientationOfA, orientationOfB, alpha, A, B, beta, C );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elemental
