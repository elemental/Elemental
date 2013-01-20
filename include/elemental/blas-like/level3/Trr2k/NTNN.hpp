/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef BLAS_TRR2K_NTNN_HPP
#define BLAS_TRR2K_NTNN_HPP

#include "elemental/blas-like/level3/Trr2k/NNNT.hpp"

namespace elem {
namespace internal {

// Distributed E := alpha (A B^{T/H} + C D) + beta E
template<typename T>
inline void
Trr2kNTNN
( UpperOrLower uplo,
  Orientation orientationOfB,
  T alpha, const DistMatrix<T>& A, const DistMatrix<T>& B,
           const DistMatrix<T>& C, const DistMatrix<T>& D,
  T beta,        DistMatrix<T>& E )
{
#ifndef RELEASE
    PushCallStack("internal::Trr2kNTNN");
#endif
    Trr2kNNNT( uplo, orientationOfB, alpha, C, D, A, B, beta, E );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace internal
} // namespace elem

#endif // ifndef BLAS_TRR2K_NTNN_HPP
