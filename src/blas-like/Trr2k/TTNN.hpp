/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef ELEM_TRR2K_TTNN_HPP
#define ELEM_TRR2K_TTNN_HPP

#include "./NNTT.hpp"

namespace elem {
namespace internal {

// Distributed E := alpha (A^{T/H} B^{T/H} + C D) + beta E
template<typename T>
void Trr2kTTNN
( UpperOrLower uplo,
  Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const DistMatrix<T>& A, const DistMatrix<T>& B,
           const DistMatrix<T>& C, const DistMatrix<T>& D,
  T beta,        DistMatrix<T>& E )
{
    DEBUG_ONLY(CallStackEntry cse("internal::Trr2kTTNN"))
    Trr2kNNTT
    ( uplo, orientationOfA, orientationOfB, alpha, C, D, A, B, beta, E );
}

} // namespace internal
} // namespace elem

#endif // ifndef ELEM_TRR2K_TTNN_HPP
