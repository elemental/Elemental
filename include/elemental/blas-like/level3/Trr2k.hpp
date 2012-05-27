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

#include "./Trr2k/Local.hpp"
#include "./Trr2k/NNNN.hpp"
#include "./Trr2k/NNNT.hpp"
#include "./Trr2k/NNTN.hpp"
#include "./Trr2k/NNTT.hpp"
#include "./Trr2k/NTNN.hpp"
#include "./Trr2k/NTNT.hpp"
#include "./Trr2k/NTTN.hpp"
#include "./Trr2k/NTTT.hpp"
#include "./Trr2k/TNNN.hpp"
#include "./Trr2k/TNNT.hpp"
#include "./Trr2k/TNTN.hpp"
#include "./Trr2k/TNTT.hpp"
#include "./Trr2k/TTNN.hpp"
#include "./Trr2k/TTNT.hpp"
#include "./Trr2k/TTTN.hpp"
#include "./Trr2k/TTTT.hpp"

namespace elem {

// This will be enabled as soon as the underlying routines are written
/*
template<typename T>
inline void
Trr2k
( UpperOrLower uplo, 
  Orientation orientationOfA, Orientation orientationOfB,
  Orientation orientationOfC, Orientation orientationOfD,
  T alpha, const Matrix<T>& A, const Matrix<T>& B,
           const Matrix<T>& C, const Matrix<T>& D,
  T beta,        Matrix<T>& E )
{
#ifndef RELEASE
    PushCallStack("Trr2k");
#endif
    const bool normalA = orientationOfA == NORMAL;
    const bool normalB = orientationOfB == NORMAL;
    const bool normalC = orientationOfC == NORMAL;
    const bool normalD = orientationOfD == NORMAL;
    int subcase = 8*normalA + 4*normalB + 2*normalC + normalD;
    switch( subcase )
    {
    case 0: 
        internal::Trr2kNNNN( uplo, alpha, A, B, C, D, beta, E );
        break;
    case 1:
        internal::Trr2kNNNT
        ( uplo, orientationOfD, alpha, A, B, C, D, beta, E );
        break;
    case 2:
        internal::Trr2kNNTN
        ( uplo, orientationOfC, alpha, A, B, C, D, beta, E );
        break;
    case 3:
        internal::Trr2kNNTT
        ( uplo, orientationOfC, orientationOfD, alpha, A, B, C, D, beta, E );
        break;
    case 4:
        internal::Trr2kNTNN
        ( uplo, orientationOfB, alpha, A, B, C, D, beta, E );
        break;
    case 5:
        internal::Trr2kNTNT
        ( uplo, orientationOfB, orientationOfD, alpha, A, B, C, D, beta, E );
        break;
    case 6:
        internal::Trr2kNTTN
        ( uplo, orientationOfB, orientationOfC, alpha, A, B, C, D, beta, E );
        break;
    case 7:
        internal::Trr2kNTTT
        ( uplo, orientationOfB, orientationOfC, orientationOfD, 
          alpha, A, B, C, D, beta, E );
        break;
    case 8:
        internal::Trr2kTNNN
        ( uplo, orientationOfA, alpha, A, B, C, D, beta, E );
        break;
    case 9:
        internal::Trr2kTNNT
        ( uplo, orientationOfA, orientationOfD, alpha, A, B, C, D, beta, E );
        break;
    case 10:
        internal::Trr2kTNTN
        ( uplo, orientationOfA, orientationOfC, alpha, A, B, C, D, beta, E );
        break;
    case 11:
        internal::Trr2kTNTT
        ( uplo, orientationOfA, orientationOfC, orientationOfD,
          alpha, A, B, C, D, beta, E );
        break;
    case 12:
        internal::Trr2kTTNN
        ( uplo, orientationOfA, orientationOfB, alpha, A, B, C, D, beta, E );
        break;
    case 13:
        internal::Trr2kTTNT
        ( uplo, orientationOfA, orientationOfB, orientationOfD,
          alpha, A, B, C, D, beta, E );
        break;
    case 14:
        internal::Trr2kTTTN
        ( uplo, orientationOfA, orientationOfB, orientationOfC,
          alpha, A, B, C, D, beta, E );
        break;
    case 15:
        internal::Trr2kTTTN
        ( uplo, orientationOfA, orientationOfB, orientationOfC, orientationOfD,
          alpha, A, B, C, D, beta, E );
        break;
    default:
        throw std::logic_error("Impossible subcase");
    }
#ifndef RELEASE
    PopCallStack();
#endif
}
*/

template<typename T>
inline void
Trr2k
( UpperOrLower uplo, 
  Orientation orientationOfA, Orientation orientationOfB,
  Orientation orientationOfC, Orientation orientationOfD,
  T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
           const DistMatrix<T,MC,MR>& C,
           const DistMatrix<T,MC,MR>& D,
  T beta,        DistMatrix<T,MC,MR>& E )
{
#ifndef RELEASE
    PushCallStack("Trr2k");
#endif
    const bool normalA = orientationOfA == NORMAL;
    const bool normalB = orientationOfB == NORMAL;
    const bool normalC = orientationOfC == NORMAL;
    const bool normalD = orientationOfD == NORMAL;
    int subcase = 8*normalA + 4*normalB + 2*normalC + normalD;
    switch( subcase )
    {
    case 0: 
        internal::Trr2kNNNN( uplo, alpha, A, B, C, D, beta, E );
        break;
    case 1:
        internal::Trr2kNNNT
        ( uplo, orientationOfD, alpha, A, B, C, D, beta, E );
        break;
    case 2:
        internal::Trr2kNNTN
        ( uplo, orientationOfC, alpha, A, B, C, D, beta, E );
        break;
    case 3:
        internal::Trr2kNNTT
        ( uplo, orientationOfC, orientationOfD, alpha, A, B, C, D, beta, E );
        break;
    case 4:
        internal::Trr2kNTNN
        ( uplo, orientationOfB, alpha, A, B, C, D, beta, E );
        break;
    case 5:
        internal::Trr2kNTNT
        ( uplo, orientationOfB, orientationOfD, alpha, A, B, C, D, beta, E );
        break;
    case 6:
        internal::Trr2kNTTN
        ( uplo, orientationOfB, orientationOfC, alpha, A, B, C, D, beta, E );
        break;
    case 7:
        internal::Trr2kNTTT
        ( uplo, orientationOfB, orientationOfC, orientationOfD, 
          alpha, A, B, C, D, beta, E );
        break;
    case 8:
        internal::Trr2kTNNN
        ( uplo, orientationOfA, alpha, A, B, C, D, beta, E );
        break;
    case 9:
        internal::Trr2kTNNT
        ( uplo, orientationOfA, orientationOfD, alpha, A, B, C, D, beta, E );
        break;
    case 10:
        internal::Trr2kTNTN
        ( uplo, orientationOfA, orientationOfC, alpha, A, B, C, D, beta, E );
        break;
    case 11:
        internal::Trr2kTNTT
        ( uplo, orientationOfA, orientationOfC, orientationOfD,
          alpha, A, B, C, D, beta, E );
        break;
    case 12:
        internal::Trr2kTTNN
        ( uplo, orientationOfA, orientationOfB, alpha, A, B, C, D, beta, E );
        break;
    case 13:
        internal::Trr2kTTNT
        ( uplo, orientationOfA, orientationOfB, orientationOfD,
          alpha, A, B, C, D, beta, E );
        break;
    case 14:
        internal::Trr2kTTTN
        ( uplo, orientationOfA, orientationOfB, orientationOfC,
          alpha, A, B, C, D, beta, E );
        break;
    case 15:
        internal::Trr2kTTTN
        ( uplo, orientationOfA, orientationOfB, orientationOfC, orientationOfD,
          alpha, A, B, C, D, beta, E );
        break;
    default:
        throw std::logic_error("Impossible subcase");
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
