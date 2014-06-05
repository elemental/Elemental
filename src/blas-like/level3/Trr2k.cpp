/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

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

// This will be enabled as soon as the underlying routines are written
/*
template<typename T>
void Trr2k
( UpperOrLower uplo, 
  Orientation orientationOfA, Orientation orientationOfB,
  Orientation orientationOfC, Orientation orientationOfD,
  T alpha, const Matrix<T>& A, const Matrix<T>& B,
           const Matrix<T>& C, const Matrix<T>& D,
  T beta,        Matrix<T>& E )
{
    DEBUG_ONLY(CallStackEntry cse("Trr2k"))
    const bool normalA = orientationOfA == NORMAL;
    const bool normalB = orientationOfB == NORMAL;
    const bool normalC = orientationOfC == NORMAL;
    const bool normalD = orientationOfD == NORMAL;
    Int subcase = 8*normalA + 4*normalB + 2*normalC + normalD;
    switch( subcase )
    {
    case 0: 
        trr2k::Trr2kNNNN( uplo, alpha, A, B, C, D, beta, E );
        break;
    case 1:
        trr2k::Trr2kNNNT
        ( uplo, orientationOfD, alpha, A, B, C, D, beta, E );
        break;
    case 2:
        trr2k::Trr2kNNTN
        ( uplo, orientationOfC, alpha, A, B, C, D, beta, E );
        break;
    case 3:
        trr2k::Trr2kNNTT
        ( uplo, orientationOfC, orientationOfD, alpha, A, B, C, D, beta, E );
        break;
    case 4:
        trr2k::Trr2kNTNN
        ( uplo, orientationOfB, alpha, A, B, C, D, beta, E );
        break;
    case 5:
        trr2k::Trr2kNTNT
        ( uplo, orientationOfB, orientationOfD, alpha, A, B, C, D, beta, E );
        break;
    case 6:
        trr2k::Trr2kNTTN
        ( uplo, orientationOfB, orientationOfC, alpha, A, B, C, D, beta, E );
        break;
    case 7:
        trr2k::Trr2kNTTT
        ( uplo, orientationOfB, orientationOfC, orientationOfD, 
          alpha, A, B, C, D, beta, E );
        break;
    case 8:
        trr2k::Trr2kTNNN
        ( uplo, orientationOfA, alpha, A, B, C, D, beta, E );
        break;
    case 9:
        trr2k::Trr2kTNNT
        ( uplo, orientationOfA, orientationOfD, alpha, A, B, C, D, beta, E );
        break;
    case 10:
        trr2k::Trr2kTNTN
        ( uplo, orientationOfA, orientationOfC, alpha, A, B, C, D, beta, E );
        break;
    case 11:
        trr2k::Trr2kTNTT
        ( uplo, orientationOfA, orientationOfC, orientationOfD,
          alpha, A, B, C, D, beta, E );
        break;
    case 12:
        trr2k::Trr2kTTNN
        ( uplo, orientationOfA, orientationOfB, alpha, A, B, C, D, beta, E );
        break;
    case 13:
        trr2k::Trr2kTTNT
        ( uplo, orientationOfA, orientationOfB, orientationOfD,
          alpha, A, B, C, D, beta, E );
        break;
    case 14:
        trr2k::Trr2kTTTN
        ( uplo, orientationOfA, orientationOfB, orientationOfC,
          alpha, A, B, C, D, beta, E );
        break;
    case 15:
        trr2k::Trr2kTTTN
        ( uplo, orientationOfA, orientationOfB, orientationOfC, orientationOfD,
          alpha, A, B, C, D, beta, E );
        break;
    default:
        LogicError("Impossible subcase");
    }
}
*/

template<typename T>
void Trr2k
( UpperOrLower uplo, 
  Orientation orientationOfA, Orientation orientationOfB,
  Orientation orientationOfC, Orientation orientationOfD,
  T alpha, const DistMatrix<T>& A, const DistMatrix<T>& B,
           const DistMatrix<T>& C, const DistMatrix<T>& D,
  T beta,        DistMatrix<T>& E )
{
    DEBUG_ONLY(CallStackEntry cse("Trr2k"))
    const bool normalA = orientationOfA == NORMAL;
    const bool normalB = orientationOfB == NORMAL;
    const bool normalC = orientationOfC == NORMAL;
    const bool normalD = orientationOfD == NORMAL;
    Int subcase = 8*normalA + 4*normalB + 2*normalC + normalD;
    switch( subcase )
    {
    case 0: 
        trr2k::Trr2kNNNN( uplo, alpha, A, B, C, D, beta, E );
        break;
    case 1:
        trr2k::Trr2kNNNT
        ( uplo, orientationOfD, alpha, A, B, C, D, beta, E );
        break;
    case 2:
        trr2k::Trr2kNNTN
        ( uplo, orientationOfC, alpha, A, B, C, D, beta, E );
        break;
    case 3:
        trr2k::Trr2kNNTT
        ( uplo, orientationOfC, orientationOfD, alpha, A, B, C, D, beta, E );
        break;
    case 4:
        trr2k::Trr2kNTNN
        ( uplo, orientationOfB, alpha, A, B, C, D, beta, E );
        break;
    case 5:
        trr2k::Trr2kNTNT
        ( uplo, orientationOfB, orientationOfD, alpha, A, B, C, D, beta, E );
        break;
    case 6:
        trr2k::Trr2kNTTN
        ( uplo, orientationOfB, orientationOfC, alpha, A, B, C, D, beta, E );
        break;
    case 7:
        trr2k::Trr2kNTTT
        ( uplo, orientationOfB, orientationOfC, orientationOfD, 
          alpha, A, B, C, D, beta, E );
        break;
    case 8:
        trr2k::Trr2kTNNN
        ( uplo, orientationOfA, alpha, A, B, C, D, beta, E );
        break;
    case 9:
        trr2k::Trr2kTNNT
        ( uplo, orientationOfA, orientationOfD, alpha, A, B, C, D, beta, E );
        break;
    case 10:
        trr2k::Trr2kTNTN
        ( uplo, orientationOfA, orientationOfC, alpha, A, B, C, D, beta, E );
        break;
    case 11:
        trr2k::Trr2kTNTT
        ( uplo, orientationOfA, orientationOfC, orientationOfD,
          alpha, A, B, C, D, beta, E );
        break;
    case 12:
        trr2k::Trr2kTTNN
        ( uplo, orientationOfA, orientationOfB, alpha, A, B, C, D, beta, E );
        break;
    case 13:
        trr2k::Trr2kTTNT
        ( uplo, orientationOfA, orientationOfB, orientationOfD,
          alpha, A, B, C, D, beta, E );
        break;
    case 14:
        trr2k::Trr2kTTTN
        ( uplo, orientationOfA, orientationOfB, orientationOfC,
          alpha, A, B, C, D, beta, E );
        break;
    case 15:
        trr2k::Trr2kTTTT
        ( uplo, orientationOfA, orientationOfB, orientationOfC, orientationOfD,
          alpha, A, B, C, D, beta, E );
        break;
    default:
        LogicError("Impossible subcase");
    }
}

#define PROTO(T) \
  template void Trr2k \
  ( UpperOrLower uplo, \
    Orientation orientationOfA, Orientation orientationOfB, \
    Orientation orientationOfC, Orientation orientationOfD, \
    T alpha, const DistMatrix<T>& A, const DistMatrix<T>& B, \
             const DistMatrix<T>& C, const DistMatrix<T>& D, \
    T beta,        DistMatrix<T>& E ); \
  template void LocalTrr2k \
  ( UpperOrLower uplo, \
    T alpha, const DistMatrix<T,MC,  STAR>& A, \
             const DistMatrix<T,STAR,MR  >& B, \
             const DistMatrix<T,MC,  STAR>& C, \
             const DistMatrix<T,STAR,MR  >& D, \
    T beta,        DistMatrix<T>& E ); \
  template void LocalTrr2k \
  ( UpperOrLower uplo, Orientation orientationOfD, \
    T alpha, const DistMatrix<T,MC,  STAR>& A,  \
             const DistMatrix<T,STAR,MR  >& B, \
             const DistMatrix<T,MC,  STAR>& C, \
             const DistMatrix<T,MR,  STAR>& D, \
    T beta,        DistMatrix<T>& E ); \
  template void LocalTrr2k \
  ( UpperOrLower uplo, Orientation orientationOfC, \
    T alpha, const DistMatrix<T,MC,  STAR>& A,  \
             const DistMatrix<T,STAR,MR  >& B, \
             const DistMatrix<T,STAR,MC  >& C, \
             const DistMatrix<T,STAR,MR  >& D, \
    T beta,        DistMatrix<T>& E ); \
  template void LocalTrr2k \
  ( UpperOrLower uplo, Orientation orientationOfC, Orientation orientationOfD, \
    T alpha, const DistMatrix<T,MC,  STAR>& A, \
             const DistMatrix<T,STAR,MR  >& B, \
             const DistMatrix<T,STAR,MC  >& C, \
             const DistMatrix<T,MR,  STAR>& D, \
    T beta,        DistMatrix<T>& E ); \
  template void LocalTrr2k \
  ( UpperOrLower uplo, Orientation orientationOfB, \
    T alpha, const DistMatrix<T,MC,  STAR>& A, \
             const DistMatrix<T,MR,  STAR>& B, \
             const DistMatrix<T,MC,  STAR>& C, \
             const DistMatrix<T,STAR,MR  >& D, \
    T beta,        DistMatrix<T>& E ); \
  template void LocalTrr2k \
  ( UpperOrLower uplo, Orientation orientationOfB, Orientation orientationOfD, \
    T alpha, const DistMatrix<T,MC,STAR>& A, \
             const DistMatrix<T,MR,STAR>& B, \
             const DistMatrix<T,MC,STAR>& C, \
             const DistMatrix<T,MR,STAR>& D, \
    T beta,        DistMatrix<T>& E ); \
  template void LocalTrr2k \
  ( UpperOrLower uplo, Orientation orientationOfB, Orientation orientationOfC, \
    T alpha, const DistMatrix<T,MC,  STAR>& A, \
             const DistMatrix<T,MR,  STAR>& B, \
             const DistMatrix<T,STAR,MC  >& C, \
             const DistMatrix<T,STAR,MR  >& D, \
    T beta,        DistMatrix<T>& E ); \
  template void LocalTrr2k \
  ( UpperOrLower uplo,          Orientation orientationOfB, \
    Orientation orientationOfC, Orientation orientationOfD, \
    T alpha, const DistMatrix<T,MC,  STAR>& A, \
             const DistMatrix<T,MR,  STAR>& B, \
             const DistMatrix<T,STAR,MC  >& C, \
             const DistMatrix<T,MR,  STAR>& D, \
    T beta,        DistMatrix<T>& E ); \
  template void LocalTrr2k \
  ( UpperOrLower uplo, Orientation orientationOfA, \
    T alpha, const DistMatrix<T,STAR,MC  >& A, \
             const DistMatrix<T,STAR,MR  >& B, \
             const DistMatrix<T,MC,  STAR>& C, \
             const DistMatrix<T,STAR,MR  >& D, \
    T beta,        DistMatrix<T>& E ); \
  template void LocalTrr2k \
  ( UpperOrLower uplo, Orientation orientationOfA, Orientation orientationOfD, \
    T alpha, const DistMatrix<T,STAR,MC  >& A, \
             const DistMatrix<T,STAR,MR>& B, \
             const DistMatrix<T,MC,  STAR>& C, \
             const DistMatrix<T,MR,STAR>& D, \
    T beta,        DistMatrix<T>& E ); \
  template void LocalTrr2k \
  ( UpperOrLower uplo, Orientation orientationOfA, Orientation orientationOfC, \
    T alpha, const DistMatrix<T,STAR,MC>& A, \
             const DistMatrix<T,STAR,MR>& B, \
             const DistMatrix<T,STAR,MC>& C, \
             const DistMatrix<T,STAR,MR>& D, \
    T beta,        DistMatrix<T>& E ); \
  template void LocalTrr2k \
  ( UpperOrLower uplo,          Orientation orientationOfA, \
    Orientation orientationOfC, Orientation orientationOfD, \
    T alpha, const DistMatrix<T,STAR,MC>& A, \
             const DistMatrix<T,STAR,MR>& B, \
             const DistMatrix<T,STAR,MC>& C, \
             const DistMatrix<T,MR,STAR>& D, \
    T beta,        DistMatrix<T>& E ); \
  template void LocalTrr2k \
  ( UpperOrLower uplo, Orientation orientationOfA, Orientation orientationOfB, \
    T alpha, const DistMatrix<T,STAR,MC  >& A, \
             const DistMatrix<T,MR,  STAR>& B, \
             const DistMatrix<T,MC,  STAR>& C, \
             const DistMatrix<T,STAR,MR  >& D, \
    T beta,        DistMatrix<T>& E ); \
  template void LocalTrr2k \
  ( UpperOrLower uplo,          Orientation orientationOfA, \
    Orientation orientationOfB, Orientation orientationOfD, \
    T alpha, const DistMatrix<T,STAR,MC  >& A, \
             const DistMatrix<T,MR,  STAR>& B, \
             const DistMatrix<T,MC,  STAR>& C, \
             const DistMatrix<T,MR,  STAR>& D, \
    T beta,        DistMatrix<T>& E ); \
  template void LocalTrr2k \
  ( UpperOrLower uplo,          Orientation orientationOfA, \
    Orientation orientationOfB, Orientation orientationOfC, \
    T alpha, const DistMatrix<T,STAR,MC  >& A, \
             const DistMatrix<T,MR,  STAR>& B, \
             const DistMatrix<T,STAR,MC  >& C, \
             const DistMatrix<T,STAR,MR  >& D, \
    T beta,        DistMatrix<T>& E ); \
  template void LocalTrr2k \
  ( UpperOrLower uplo, \
    Orientation orientationOfA, Orientation orientationOfB, \
    Orientation orientationOfC, Orientation orientationOfD, \
    T alpha, const DistMatrix<T,STAR,MC  >& A, \
             const DistMatrix<T,MR,  STAR>& B, \
             const DistMatrix<T,STAR,MC  >& C, \
             const DistMatrix<T,MR,  STAR>& D, \
    T beta,        DistMatrix<T>& E  ); \

#ifndef EL_DISABLE_FLOAT
PROTO(float);
#ifndef EL_DISABLE_COMPLEX
PROTO(Complex<float>);
#endif // ifndef EL_DISABLE_COMPLEX
#endif // ifndef EL_DISABLE_FLOAT

PROTO(double);
#ifndef EL_DISABLE_COMPLEX
PROTO(Complex<double>);
#endif // ifndef EL_DISABLE_COMPLEX

} // namespace El
