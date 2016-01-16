/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

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

namespace El {

template<typename T>
void Trr2k
( UpperOrLower uplo, 
  Orientation orientA, Orientation orientB,
  Orientation orientC, Orientation orientD,
  T alpha, const ElementalMatrix<T>& A, const ElementalMatrix<T>& B,
  T beta,  const ElementalMatrix<T>& C, const ElementalMatrix<T>& D,
  T gamma,       ElementalMatrix<T>& E )
{
    DEBUG_ONLY(CSE cse("Trr2k"))
    const bool normalA = orientA == NORMAL;
    const bool normalB = orientB == NORMAL;
    const bool normalC = orientC == NORMAL;
    const bool normalD = orientD == NORMAL;
    Int subcase = 8*normalA + 4*normalB + 2*normalC + normalD;
    ScaleTrapezoid( gamma, uplo, E );
    switch( subcase )
    {
    case 0: 
        trr2k::Trr2kNNNN( uplo, alpha, A, B, beta, C, D, E );
        break;
    case 1:
        trr2k::Trr2kNNNT
        ( uplo, orientD, alpha, A, B, beta, C, D, E );
        break;
    case 2:
        trr2k::Trr2kNNTN
        ( uplo, orientC, alpha, A, B, beta, C, D, E );
        break;
    case 3:
        trr2k::Trr2kNNTT
        ( uplo, orientC, orientD, 
          alpha, A, B, beta, C, D, E );
        break;
    case 4:
        trr2k::Trr2kNTNN
        ( uplo, orientB, alpha, A, B, beta, C, D, E );
        break;
    case 5:
        trr2k::Trr2kNTNT
        ( uplo, orientB, orientD, 
          alpha, A, B, beta, C, D, E );
        break;
    case 6:
        trr2k::Trr2kNTTN
        ( uplo, orientB, orientC, 
          alpha, A, B, beta, C, D, E );
        break;
    case 7:
        trr2k::Trr2kNTTT
        ( uplo, orientB, orientC, orientD, 
          alpha, A, B, beta, C, D, E );
        break;
    case 8:
        trr2k::Trr2kTNNN
        ( uplo, orientA, alpha, A, B, beta, C, D, E );
        break;
    case 9:
        trr2k::Trr2kTNNT
        ( uplo, orientA, orientD, 
          alpha, A, B, beta, C, D, E );
        break;
    case 10:
        trr2k::Trr2kTNTN
        ( uplo, orientA, orientC, 
          alpha, A, B, beta, C, D, E );
        break;
    case 11:
        trr2k::Trr2kTNTT
        ( uplo, orientA, orientC, orientD,
          alpha, A, B, beta, C, D, E );
        break;
    case 12:
        trr2k::Trr2kTTNN
        ( uplo, orientA, orientB, 
          alpha, A, B, beta, C, D, E );
        break;
    case 13:
        trr2k::Trr2kTTNT
        ( uplo, orientA, orientB, orientD,
          alpha, A, B, beta, C, D, E );
        break;
    case 14:
        trr2k::Trr2kTTTN
        ( uplo, orientA, orientB, orientC,
          alpha, A, B, beta, C, D, E );
        break;
    case 15:
        trr2k::Trr2kTTTT
        ( uplo, orientA, orientB, orientC, orientD,
          alpha, A, B, beta, C, D, E );
        break;
    default:
        LogicError("Impossible subcase");
    }
}

#define PROTO(T) \
  template void Trr2k \
  ( UpperOrLower uplo, \
    Orientation orientA, Orientation orientB, \
    Orientation orientC, Orientation orientD, \
    T alpha, const ElementalMatrix<T>& A, const ElementalMatrix<T>& B, \
    T beta,  const ElementalMatrix<T>& C, const ElementalMatrix<T>& D, \
    T gamma,       ElementalMatrix<T>& E ); \
  template void LocalTrr2k \
  ( UpperOrLower uplo, \
    Orientation orientA, Orientation orientB, \
    Orientation orientC, Orientation orientD, \
    T alpha, const ElementalMatrix<T>& A, const ElementalMatrix<T>& B, \
    T beta,  const ElementalMatrix<T>& C, const ElementalMatrix<T>& D, \
    T gamma,       ElementalMatrix<T>& E );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

} // namespace El
