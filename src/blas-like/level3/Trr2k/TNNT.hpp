/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_TRR2K_TNNT_HPP
#define EL_TRR2K_TNNT_HPP

#include "./NTTN.hpp"

namespace El {
namespace trr2k {

// Distributed E := alpha (A^{T/H} B + C D^{T/H}) + beta E
template<typename T>
void Trr2kTNNT
( UpperOrLower uplo,
  Orientation orientationOfA, Orientation orientationOfD,
  T alpha, const DistMatrix<T>& A, const DistMatrix<T>& B,
           const DistMatrix<T>& C, const DistMatrix<T>& D,
  T beta,        DistMatrix<T>& E )
{
    DEBUG_ONLY(CallStackEntry cse("trr2k::Trr2kTNNT"))
    Trr2kNTTN
    ( uplo, orientationOfD, orientationOfA, alpha, C, D, A, B, beta, E );
}

} // namespace trr2k
} // namespace El

#endif // ifndef EL_TRR2K_TNNT_HPP
