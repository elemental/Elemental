/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace elem {
namespace internal {

// Distributed E := alpha (A^{T/H} B + C D) + beta E
template<typename T>
inline void
Trr2kTNNN
( UpperOrLower uplo,
  Orientation orientationOfA,
  T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B,
           const DistMatrix<T,MC,MR>& C, const DistMatrix<T,MC,MR>& D,
  T beta,        DistMatrix<T,MC,MR>& E )
{
#ifndef RELEASE
    PushCallStack("internal::Trr2kTNNN");
#endif
    Trr2kNNTN( uplo, orientationOfA, alpha, C, D, A, B, beta, E );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace internal
} // namespace elem
