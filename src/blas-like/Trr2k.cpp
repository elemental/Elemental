/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "elemental-lite.hpp"

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
void Trr2k
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
void Trr2k
( UpperOrLower uplo, 
  Orientation orientationOfA, Orientation orientationOfB,
  Orientation orientationOfC, Orientation orientationOfD,
  T alpha, const DistMatrix<T>& A, const DistMatrix<T>& B,
           const DistMatrix<T>& C, const DistMatrix<T>& D,
  T beta,        DistMatrix<T>& E )
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
        internal::Trr2kTTTT
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

#ifndef DISABLE_FLOAT
template void Trr2k
( UpperOrLower uplo,
  Orientation orientationOfA, Orientation orientationOfB,
  Orientation orientationOfC, Orientation orientationOfD,
  float alpha, const DistMatrix<float>& A, const DistMatrix<float>& B,
               const DistMatrix<float>& C, const DistMatrix<float>& D,
  float beta,        DistMatrix<float>& E );
template void LocalTrr2k
( UpperOrLower uplo,
  float alpha, const DistMatrix<float,MC,  STAR>& A, 
               const DistMatrix<float,STAR,MR>& B,
               const DistMatrix<float,MC,  STAR>& C, 
               const DistMatrix<float,STAR,MR>& D,
  float beta,        DistMatrix<float,MC,  MR  >& E );
template void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfD,
  float alpha, const DistMatrix<float,MC,  STAR>& A, 
               const DistMatrix<float,STAR,MR>& B,
               const DistMatrix<float,MC,  STAR>& C, 
               const DistMatrix<float,MR,STAR>& D,
  float beta,        DistMatrix<float,MC,  MR  >& E  );
template void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfC,
  float alpha, const DistMatrix<float,MC,  STAR>& A, 
               const DistMatrix<float,STAR,MR>& B,
               const DistMatrix<float,STAR,MC  >& C, 
               const DistMatrix<float,STAR,MR>& D,
  float beta,        DistMatrix<float,MC,  MR  >& E  );
template void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfC,
  Orientation orientationOfD,
  float alpha, const DistMatrix<float,MC,  STAR>& A, 
               const DistMatrix<float,STAR,MR>& B,
               const DistMatrix<float,STAR,MC  >& C, 
               const DistMatrix<float,MR,STAR>& D,
  float beta,        DistMatrix<float,MC,  MR  >& E  );
template void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfB,
  float alpha, const DistMatrix<float,MC,  STAR>& A, 
               const DistMatrix<float,MR,STAR>& B,
               const DistMatrix<float,MC,  STAR>& C, 
               const DistMatrix<float,STAR,MR>& D,
  float beta,        DistMatrix<float,MC,  MR  >& E  );
template void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfB,
  Orientation orientationOfD,
  float alpha, const DistMatrix<float,MC,STAR>& A, 
               const DistMatrix<float,MR,STAR>& B,
               const DistMatrix<float,MC,STAR>& C, 
               const DistMatrix<float,MR,STAR>& D,
  float beta,        DistMatrix<float,MC,MR  >& E  );
template void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfB,
  Orientation orientationOfC,
  float alpha, const DistMatrix<float,MC,  STAR>& A, 
               const DistMatrix<float,MR,STAR>& B,
               const DistMatrix<float,STAR,MC  >& C, 
               const DistMatrix<float,STAR,MR>& D,
  float beta,        DistMatrix<float,MC,  MR  >& E  );
template void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfB,
  Orientation orientationOfC,
  Orientation orientationOfD,
  float alpha, const DistMatrix<float,MC,  STAR>& A, 
               const DistMatrix<float,MR,STAR>& B,
               const DistMatrix<float,STAR,MC  >& C, 
               const DistMatrix<float,MR,STAR>& D,
  float beta,        DistMatrix<float,MC,  MR  >& E  );
template void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfA,
  float alpha, const DistMatrix<float,STAR,MC  >& A, 
               const DistMatrix<float,STAR,MR>& B,
               const DistMatrix<float,MC,  STAR>& C, 
               const DistMatrix<float,STAR,MR>& D,
  float beta,        DistMatrix<float,MC,  MR  >& E  );
template void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfD,
  float alpha, const DistMatrix<float,STAR,MC  >& A, 
               const DistMatrix<float,STAR,MR>& B,
               const DistMatrix<float,MC,  STAR>& C, 
               const DistMatrix<float,MR,STAR>& D,
  float beta,        DistMatrix<float,MC,  MR  >& E  );
template void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfC,
  float alpha, const DistMatrix<float,STAR,MC>& A, 
               const DistMatrix<float,STAR,MR>& B,
               const DistMatrix<float,STAR,MC>& C, 
               const DistMatrix<float,STAR,MR>& D,
  float beta,        DistMatrix<float,MC,  MR>& E  );
template void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfC,
  Orientation orientationOfD,
  float alpha, const DistMatrix<float,STAR,MC>& A, 
               const DistMatrix<float,STAR,MR>& B,
               const DistMatrix<float,STAR,MC>& C, 
               const DistMatrix<float,MR,STAR>& D,
  float beta,        DistMatrix<float,MC,  MR>& E  );
template void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfB,
  float alpha, const DistMatrix<float,STAR,MC  >& A, 
               const DistMatrix<float,MR,STAR>& B,
               const DistMatrix<float,MC,  STAR>& C, 
               const DistMatrix<float,STAR,MR>& D,
  float beta,        DistMatrix<float,MC,  MR  >& E  );
template void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfB,
  Orientation orientationOfD,
  float alpha, const DistMatrix<float,STAR,MC  >& A, 
               const DistMatrix<float,MR,STAR>& B,
               const DistMatrix<float,MC,  STAR>& C, 
               const DistMatrix<float,MR,STAR>& D,
  float beta,        DistMatrix<float,MC,  MR  >& E  );
template void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfB,
  Orientation orientationOfC,
  float alpha, const DistMatrix<float,STAR,MC>& A, 
               const DistMatrix<float,MR,STAR>& B,
               const DistMatrix<float,STAR,MC>& C, 
               const DistMatrix<float,STAR,MR>& D,
  float beta,        DistMatrix<float,MC,  MR>& E  );
template void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfB,
  Orientation orientationOfC,
  Orientation orientationOfD,
  float alpha, const DistMatrix<float,STAR,MC>& A, 
               const DistMatrix<float,MR,STAR>& B,
               const DistMatrix<float,STAR,MC>& C, 
               const DistMatrix<float,MR,STAR>& D,
  float beta,        DistMatrix<float,MC,  MR>& E  );
#ifndef DISABLE_COMPLEX
template void Trr2k
( UpperOrLower uplo,
  Orientation orientationOfA, Orientation orientationOfB,
  Orientation orientationOfC, Orientation orientationOfD,
  Complex<float> alpha, const DistMatrix<Complex<float> >& A, 
                        const DistMatrix<Complex<float> >& B,
                        const DistMatrix<Complex<float> >& C, 
                        const DistMatrix<Complex<float> >& D,
  Complex<float> beta,        DistMatrix<Complex<float> >& E );
template void LocalTrr2k
( UpperOrLower uplo,
  Complex<float> alpha, const DistMatrix<Complex<float>,MC,  STAR>& A, 
                        const DistMatrix<Complex<float>,STAR,MR>& B,
                        const DistMatrix<Complex<float>,MC,  STAR>& C, 
                        const DistMatrix<Complex<float>,STAR,MR>& D,
  Complex<float> beta,        DistMatrix<Complex<float>,MC,  MR  >& E );
template void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfD,
  Complex<float> alpha, const DistMatrix<Complex<float>,MC,  STAR>& A, 
                        const DistMatrix<Complex<float>,STAR,MR>& B,
                        const DistMatrix<Complex<float>,MC,  STAR>& C, 
                        const DistMatrix<Complex<float>,MR,STAR>& D,
  Complex<float> beta,        DistMatrix<Complex<float>,MC,  MR  >& E  );
template void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfC,
  Complex<float> alpha, const DistMatrix<Complex<float>,MC,  STAR>& A, 
                        const DistMatrix<Complex<float>,STAR,MR>& B,
                        const DistMatrix<Complex<float>,STAR,MC  >& C, 
                        const DistMatrix<Complex<float>,STAR,MR>& D,
  Complex<float> beta,        DistMatrix<Complex<float>,MC,  MR  >& E  );
template void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfC,
  Orientation orientationOfD,
  Complex<float> alpha, const DistMatrix<Complex<float>,MC,  STAR>& A, 
                        const DistMatrix<Complex<float>,STAR,MR>& B,
                        const DistMatrix<Complex<float>,STAR,MC  >& C, 
                        const DistMatrix<Complex<float>,MR,STAR>& D,
  Complex<float> beta,        DistMatrix<Complex<float>,MC,  MR  >& E  );
template void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfB,
  Complex<float> alpha, const DistMatrix<Complex<float>,MC,  STAR>& A, 
                        const DistMatrix<Complex<float>,MR,STAR>& B,
                        const DistMatrix<Complex<float>,MC,  STAR>& C, 
                        const DistMatrix<Complex<float>,STAR,MR>& D,
  Complex<float> beta,        DistMatrix<Complex<float>,MC,  MR  >& E  );
template void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfB,
  Orientation orientationOfD,
  Complex<float> alpha, const DistMatrix<Complex<float>,MC,STAR>& A, 
                        const DistMatrix<Complex<float>,MR,STAR>& B,
                        const DistMatrix<Complex<float>,MC,STAR>& C, 
                        const DistMatrix<Complex<float>,MR,STAR>& D,
  Complex<float> beta,        DistMatrix<Complex<float>,MC,MR  >& E  );
template void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfB,
  Orientation orientationOfC,
  Complex<float> alpha, const DistMatrix<Complex<float>,MC,  STAR>& A, 
                        const DistMatrix<Complex<float>,MR,STAR>& B,
                        const DistMatrix<Complex<float>,STAR,MC  >& C, 
                        const DistMatrix<Complex<float>,STAR,MR>& D,
  Complex<float> beta,        DistMatrix<Complex<float>,MC,  MR  >& E  );
template void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfB,
  Orientation orientationOfC,
  Orientation orientationOfD,
  Complex<float> alpha, const DistMatrix<Complex<float>,MC,  STAR>& A, 
                        const DistMatrix<Complex<float>,MR,STAR>& B,
                        const DistMatrix<Complex<float>,STAR,MC  >& C, 
                        const DistMatrix<Complex<float>,MR,STAR>& D,
  Complex<float> beta,        DistMatrix<Complex<float>,MC,  MR  >& E  );
template void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfA,
  Complex<float> alpha, const DistMatrix<Complex<float>,STAR,MC  >& A, 
                        const DistMatrix<Complex<float>,STAR,MR>& B,
                        const DistMatrix<Complex<float>,MC,  STAR>& C, 
                        const DistMatrix<Complex<float>,STAR,MR>& D,
  Complex<float> beta,        DistMatrix<Complex<float>,MC,  MR  >& E  );
template void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfD,
  Complex<float> alpha, const DistMatrix<Complex<float>,STAR,MC  >& A, 
                        const DistMatrix<Complex<float>,STAR,MR>& B,
                        const DistMatrix<Complex<float>,MC,  STAR>& C, 
                        const DistMatrix<Complex<float>,MR,STAR>& D,
  Complex<float> beta,        DistMatrix<Complex<float>,MC,  MR  >& E  );
template void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfC,
  Complex<float> alpha, const DistMatrix<Complex<float>,STAR,MC>& A, 
                        const DistMatrix<Complex<float>,STAR,MR>& B,
                        const DistMatrix<Complex<float>,STAR,MC>& C, 
                        const DistMatrix<Complex<float>,STAR,MR>& D,
  Complex<float> beta,        DistMatrix<Complex<float>,MC,  MR>& E  );
template void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfC,
  Orientation orientationOfD,
  Complex<float> alpha, const DistMatrix<Complex<float>,STAR,MC>& A, 
                        const DistMatrix<Complex<float>,STAR,MR>& B,
                        const DistMatrix<Complex<float>,STAR,MC>& C, 
                        const DistMatrix<Complex<float>,MR,STAR>& D,
  Complex<float> beta,        DistMatrix<Complex<float>,MC,  MR>& E  );
template void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfB,
  Complex<float> alpha, const DistMatrix<Complex<float>,STAR,MC  >& A, 
                        const DistMatrix<Complex<float>,MR,STAR>& B,
                        const DistMatrix<Complex<float>,MC,  STAR>& C, 
                        const DistMatrix<Complex<float>,STAR,MR>& D,
  Complex<float> beta,        DistMatrix<Complex<float>,MC,  MR  >& E  );
template void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfB,
  Orientation orientationOfD,
  Complex<float> alpha, const DistMatrix<Complex<float>,STAR,MC  >& A, 
                        const DistMatrix<Complex<float>,MR,STAR>& B,
                        const DistMatrix<Complex<float>,MC,  STAR>& C, 
                        const DistMatrix<Complex<float>,MR,STAR>& D,
  Complex<float> beta,        DistMatrix<Complex<float>,MC,  MR  >& E  );
template void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfB,
  Orientation orientationOfC,
  Complex<float> alpha, const DistMatrix<Complex<float>,STAR,MC>& A, 
                        const DistMatrix<Complex<float>,MR,STAR>& B,
                        const DistMatrix<Complex<float>,STAR,MC>& C, 
                        const DistMatrix<Complex<float>,STAR,MR>& D,
  Complex<float> beta,        DistMatrix<Complex<float>,MC,  MR>& E  );
template void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfB,
  Orientation orientationOfC,
  Orientation orientationOfD,
  Complex<float> alpha, const DistMatrix<Complex<float>,STAR,MC>& A, 
                        const DistMatrix<Complex<float>,MR,STAR>& B,
                        const DistMatrix<Complex<float>,STAR,MC>& C, 
                        const DistMatrix<Complex<float>,MR,STAR>& D,
  Complex<float> beta,        DistMatrix<Complex<float>,MC,  MR>& E  );
#endif // ifndef DISABLE_COMPLEX
#endif // ifndef DISABLE_FLOAT

template void Trr2k
( UpperOrLower uplo,
  Orientation orientationOfA, Orientation orientationOfB,
  Orientation orientationOfC, Orientation orientationOfD,
  double alpha, const DistMatrix<double>& A, const DistMatrix<double>& B,
               const DistMatrix<double>& C, const DistMatrix<double>& D,
  double beta,        DistMatrix<double>& E );
template void LocalTrr2k
( UpperOrLower uplo,
  double alpha, const DistMatrix<double,MC,  STAR>& A, 
                const DistMatrix<double,STAR,MR>& B,
                const DistMatrix<double,MC,  STAR>& C, 
                const DistMatrix<double,STAR,MR>& D,
  double beta,        DistMatrix<double,MC,  MR  >& E );
template void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfD,
  double alpha, const DistMatrix<double,MC,  STAR>& A, 
                const DistMatrix<double,STAR,MR>& B,
                const DistMatrix<double,MC,  STAR>& C, 
                const DistMatrix<double,MR,STAR>& D,
  double beta,        DistMatrix<double,MC,  MR  >& E  );
template void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfC,
  double alpha, const DistMatrix<double,MC,  STAR>& A, 
                const DistMatrix<double,STAR,MR>& B,
                const DistMatrix<double,STAR,MC  >& C, 
                const DistMatrix<double,STAR,MR>& D,
  double beta,        DistMatrix<double,MC,  MR  >& E  );
template void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfC,
  Orientation orientationOfD,
  double alpha, const DistMatrix<double,MC,  STAR>& A, 
                const DistMatrix<double,STAR,MR>& B,
                const DistMatrix<double,STAR,MC  >& C, 
                const DistMatrix<double,MR,STAR>& D,
  double beta,        DistMatrix<double,MC,  MR  >& E  );
template void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfB,
  double alpha, const DistMatrix<double,MC,  STAR>& A, 
               const DistMatrix<double,MR,STAR>& B,
               const DistMatrix<double,MC,  STAR>& C, 
               const DistMatrix<double,STAR,MR>& D,
  double beta,        DistMatrix<double,MC,  MR  >& E  );
template void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfB,
  Orientation orientationOfD,
  double alpha, const DistMatrix<double,MC,STAR>& A, 
                const DistMatrix<double,MR,STAR>& B,
                const DistMatrix<double,MC,STAR>& C, 
                const DistMatrix<double,MR,STAR>& D,
  double beta,        DistMatrix<double,MC,MR  >& E  );
template void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfB,
  Orientation orientationOfC,
  double alpha, const DistMatrix<double,MC,  STAR>& A, 
                const DistMatrix<double,MR,STAR>& B,
                const DistMatrix<double,STAR,MC  >& C, 
                const DistMatrix<double,STAR,MR>& D,
  double beta,        DistMatrix<double,MC,  MR  >& E  );
template void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfB,
  Orientation orientationOfC,
  Orientation orientationOfD,
  double alpha, const DistMatrix<double,MC,  STAR>& A, 
                const DistMatrix<double,MR,STAR>& B,
                const DistMatrix<double,STAR,MC  >& C, 
                const DistMatrix<double,MR,STAR>& D,
  double beta,        DistMatrix<double,MC,  MR  >& E  );
template void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfA,
  double alpha, const DistMatrix<double,STAR,MC  >& A, 
                const DistMatrix<double,STAR,MR>& B,
                const DistMatrix<double,MC,  STAR>& C, 
                const DistMatrix<double,STAR,MR>& D,
  double beta,        DistMatrix<double,MC,  MR  >& E  );
template void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfD,
  double alpha, const DistMatrix<double,STAR,MC  >& A, 
                const DistMatrix<double,STAR,MR>& B,
                const DistMatrix<double,MC,  STAR>& C, 
                const DistMatrix<double,MR,STAR>& D,
  double beta,        DistMatrix<double,MC,  MR  >& E  );
template void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfC,
  double alpha, const DistMatrix<double,STAR,MC>& A, 
                const DistMatrix<double,STAR,MR>& B,
                const DistMatrix<double,STAR,MC>& C, 
                const DistMatrix<double,STAR,MR>& D,
  double beta,        DistMatrix<double,MC,  MR>& E  );
template void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfC,
  Orientation orientationOfD,
  double alpha, const DistMatrix<double,STAR,MC>& A, 
                const DistMatrix<double,STAR,MR>& B,
                const DistMatrix<double,STAR,MC>& C, 
                const DistMatrix<double,MR,STAR>& D,
  double beta,        DistMatrix<double,MC,  MR>& E  );
template void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfB,
  double alpha, const DistMatrix<double,STAR,MC  >& A, 
                const DistMatrix<double,MR,STAR>& B,
                const DistMatrix<double,MC,  STAR>& C, 
                const DistMatrix<double,STAR,MR>& D,
  double beta,        DistMatrix<double,MC,  MR  >& E  );
template void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfB,
  Orientation orientationOfD,
  double alpha, const DistMatrix<double,STAR,MC  >& A, 
                const DistMatrix<double,MR,STAR>& B,
                const DistMatrix<double,MC,  STAR>& C, 
                const DistMatrix<double,MR,STAR>& D,
  double beta,        DistMatrix<double,MC,  MR  >& E  );
template void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfB,
  Orientation orientationOfC,
  double alpha, const DistMatrix<double,STAR,MC>& A, 
                const DistMatrix<double,MR,STAR>& B,
                const DistMatrix<double,STAR,MC>& C, 
                const DistMatrix<double,STAR,MR>& D,
  double beta,        DistMatrix<double,MC,  MR>& E  );
template void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfB,
  Orientation orientationOfC,
  Orientation orientationOfD,
  double alpha, const DistMatrix<double,STAR,MC>& A, 
                const DistMatrix<double,MR,STAR>& B,
                const DistMatrix<double,STAR,MC>& C, 
                const DistMatrix<double,MR,STAR>& D,
  double beta,        DistMatrix<double,MC,  MR>& E  );
#ifndef DISABLE_COMPLEX
template void Trr2k
( UpperOrLower uplo,
  Orientation orientationOfA, Orientation orientationOfB,
  Orientation orientationOfC, Orientation orientationOfD,
  Complex<double> alpha, const DistMatrix<Complex<double> >& A, 
                         const DistMatrix<Complex<double> >& B,
                         const DistMatrix<Complex<double> >& C, 
                         const DistMatrix<Complex<double> >& D,
  Complex<double> beta,        DistMatrix<Complex<double> >& E );
template void LocalTrr2k
( UpperOrLower uplo,
  Complex<double> alpha, const DistMatrix<Complex<double>,MC,  STAR>& A, 
                         const DistMatrix<Complex<double>,STAR,MR>& B,
                         const DistMatrix<Complex<double>,MC,  STAR>& C, 
                         const DistMatrix<Complex<double>,STAR,MR>& D,
  Complex<double> beta,        DistMatrix<Complex<double>,MC,  MR  >& E );
template void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfD,
  Complex<double> alpha, const DistMatrix<Complex<double>,MC,  STAR>& A, 
                         const DistMatrix<Complex<double>,STAR,MR>& B,
                         const DistMatrix<Complex<double>,MC,  STAR>& C, 
                         const DistMatrix<Complex<double>,MR,STAR>& D,
  Complex<double> beta,        DistMatrix<Complex<double>,MC,  MR  >& E  );
template void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfC,
  Complex<double> alpha, const DistMatrix<Complex<double>,MC,  STAR>& A, 
                         const DistMatrix<Complex<double>,STAR,MR>& B,
                         const DistMatrix<Complex<double>,STAR,MC  >& C, 
                         const DistMatrix<Complex<double>,STAR,MR>& D,
  Complex<double> beta,        DistMatrix<Complex<double>,MC,  MR  >& E  );
template void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfC,
  Orientation orientationOfD,
  Complex<double> alpha, const DistMatrix<Complex<double>,MC,  STAR>& A, 
                         const DistMatrix<Complex<double>,STAR,MR>& B,
                         const DistMatrix<Complex<double>,STAR,MC  >& C, 
                         const DistMatrix<Complex<double>,MR,STAR>& D,
  Complex<double> beta,        DistMatrix<Complex<double>,MC,  MR  >& E  );
template void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfB,
  Complex<double> alpha, const DistMatrix<Complex<double>,MC,  STAR>& A, 
                         const DistMatrix<Complex<double>,MR,STAR>& B,
                         const DistMatrix<Complex<double>,MC,  STAR>& C, 
                         const DistMatrix<Complex<double>,STAR,MR>& D,
  Complex<double> beta,        DistMatrix<Complex<double>,MC,  MR  >& E  );
template void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfB,
  Orientation orientationOfD,
  Complex<double> alpha, const DistMatrix<Complex<double>,MC,STAR>& A, 
                         const DistMatrix<Complex<double>,MR,STAR>& B,
                         const DistMatrix<Complex<double>,MC,STAR>& C, 
                         const DistMatrix<Complex<double>,MR,STAR>& D,
  Complex<double> beta,        DistMatrix<Complex<double>,MC,MR  >& E  );
template void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfB,
  Orientation orientationOfC,
  Complex<double> alpha, const DistMatrix<Complex<double>,MC,  STAR>& A, 
                         const DistMatrix<Complex<double>,MR,STAR>& B,
                         const DistMatrix<Complex<double>,STAR,MC  >& C, 
                         const DistMatrix<Complex<double>,STAR,MR>& D,
  Complex<double> beta,        DistMatrix<Complex<double>,MC,  MR  >& E  );
template void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfB,
  Orientation orientationOfC,
  Orientation orientationOfD,
  Complex<double> alpha, const DistMatrix<Complex<double>,MC,  STAR>& A, 
                         const DistMatrix<Complex<double>,MR,STAR>& B,
                         const DistMatrix<Complex<double>,STAR,MC  >& C, 
                         const DistMatrix<Complex<double>,MR,STAR>& D,
  Complex<double> beta,        DistMatrix<Complex<double>,MC,  MR  >& E  );
template void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfA,
  Complex<double> alpha, const DistMatrix<Complex<double>,STAR,MC  >& A, 
                         const DistMatrix<Complex<double>,STAR,MR>& B,
                         const DistMatrix<Complex<double>,MC,  STAR>& C, 
                         const DistMatrix<Complex<double>,STAR,MR>& D,
  Complex<double> beta,        DistMatrix<Complex<double>,MC,  MR  >& E  );
template void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfD,
  Complex<double> alpha, const DistMatrix<Complex<double>,STAR,MC  >& A, 
                         const DistMatrix<Complex<double>,STAR,MR>& B,
                         const DistMatrix<Complex<double>,MC,  STAR>& C, 
                         const DistMatrix<Complex<double>,MR,STAR>& D,
  Complex<double> beta,        DistMatrix<Complex<double>,MC,  MR  >& E  );
template void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfC,
  Complex<double> alpha, const DistMatrix<Complex<double>,STAR,MC>& A, 
                         const DistMatrix<Complex<double>,STAR,MR>& B,
                         const DistMatrix<Complex<double>,STAR,MC>& C, 
                         const DistMatrix<Complex<double>,STAR,MR>& D,
  Complex<double> beta,        DistMatrix<Complex<double>,MC,  MR>& E  );
template void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfC,
  Orientation orientationOfD,
  Complex<double> alpha, const DistMatrix<Complex<double>,STAR,MC>& A, 
                         const DistMatrix<Complex<double>,STAR,MR>& B,
                         const DistMatrix<Complex<double>,STAR,MC>& C, 
                         const DistMatrix<Complex<double>,MR,STAR>& D,
  Complex<double> beta,        DistMatrix<Complex<double>,MC,  MR>& E  );
template void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfB,
  Complex<double> alpha, const DistMatrix<Complex<double>,STAR,MC  >& A, 
                         const DistMatrix<Complex<double>,MR,STAR>& B,
                         const DistMatrix<Complex<double>,MC,  STAR>& C, 
                         const DistMatrix<Complex<double>,STAR,MR>& D,
  Complex<double> beta,        DistMatrix<Complex<double>,MC,  MR  >& E  );
template void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfB,
  Orientation orientationOfD,
  Complex<double> alpha, const DistMatrix<Complex<double>,STAR,MC  >& A, 
                         const DistMatrix<Complex<double>,MR,STAR>& B,
                         const DistMatrix<Complex<double>,MC,  STAR>& C, 
                         const DistMatrix<Complex<double>,MR,STAR>& D,
  Complex<double> beta,        DistMatrix<Complex<double>,MC,  MR  >& E  );
template void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfB,
  Orientation orientationOfC,
  Complex<double> alpha, const DistMatrix<Complex<double>,STAR,MC>& A, 
                         const DistMatrix<Complex<double>,MR,STAR>& B,
                         const DistMatrix<Complex<double>,STAR,MC>& C, 
                         const DistMatrix<Complex<double>,STAR,MR>& D,
  Complex<double> beta,        DistMatrix<Complex<double>,MC,  MR>& E  );
template void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfB,
  Orientation orientationOfC,
  Orientation orientationOfD,
  Complex<double> alpha, const DistMatrix<Complex<double>,STAR,MC>& A, 
                         const DistMatrix<Complex<double>,MR,STAR>& B,
                         const DistMatrix<Complex<double>,STAR,MC>& C, 
                         const DistMatrix<Complex<double>,MR,STAR>& D,
  Complex<double> beta,        DistMatrix<Complex<double>,MC,  MR>& E  );
#endif // ifndef DISABLE_COMPLEX

} // namespace elem
