/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "elemental-lite.hpp"

#include "./Trrk/Local.hpp"
#include "./Trrk/NN.hpp"
#include "./Trrk/NT.hpp"
#include "./Trrk/TN.hpp"
#include "./Trrk/TT.hpp"

namespace elem {

template<typename T>
void Trrk
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
void Trrk
( UpperOrLower uplo, Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const DistMatrix<T>& A, const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C )
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

#ifndef DISABLE_FLOAT
template void Trrk
( UpperOrLower uplo, 
  Orientation orientationOfA, Orientation orientationOfB,
  float alpha, const Matrix<float>& A, const Matrix<float>& B,
  float beta,        Matrix<float>& C );
template void Trrk
( UpperOrLower uplo, 
  Orientation orientationOfA, Orientation orientationOfB,
  float alpha, const DistMatrix<float>& A, const DistMatrix<float>& B,
  float beta,        DistMatrix<float>& C );
template void LocalTrrk
( UpperOrLower uplo,
  float alpha, const DistMatrix<float,MC,  STAR>& A,
               const DistMatrix<float,STAR,MR  >& B,
  float beta,        DistMatrix<float,MC,  MR  >& C );
template void LocalTrrk
( UpperOrLower uplo,
  Orientation orientationOfB,
  float alpha, const DistMatrix<float,MC,STAR>& A,
               const DistMatrix<float,MR,STAR>& B,
  float beta,        DistMatrix<float>& C );
template void LocalTrrk
( UpperOrLower uplo,
  Orientation orientationOfA,
  float alpha, const DistMatrix<float,STAR,MC>& A,
               const DistMatrix<float,STAR,MR>& B,
  float beta,        DistMatrix<float,MC,  MR>& C );
template void LocalTrrk
( UpperOrLower uplo,
  Orientation orientationOfA, Orientation orientationOfB,
  float alpha, const DistMatrix<float,STAR,MC  >& A,
               const DistMatrix<float,MR,  STAR>& B,
  float beta,        DistMatrix<float,MC,  MR  >& E );
namespace internal {
template void TrrkNN
( UpperOrLower uplo,
  float alpha, const Matrix<float>& A, const Matrix<float>& B,
  float beta,        Matrix<float>& C );
template void TrrkNT
( UpperOrLower uplo,
  Orientation orientationOfB,
  float alpha, const Matrix<float>& A, const Matrix<float>& B,
  float beta,        Matrix<float>& C );
template void TrrkTN
( UpperOrLower uplo,
  Orientation orientationOfA,
  float alpha, const Matrix<float>& A, const Matrix<float>& B,
  float beta,        Matrix<float>& C );
template void TrrkTT
( UpperOrLower uplo,
  Orientation orientationOfA, Orientation orientationOfB,
  float alpha, const Matrix<float>& A, const Matrix<float>& B,
  float beta,        Matrix<float>& C );
} // namespace internal
#ifndef DISABLE_COMPLEX
template void Trrk
( UpperOrLower uplo, 
  Orientation orientationOfA, Orientation orientationOfB,
  Complex<float> alpha, const Matrix<Complex<float> >& A, 
                        const Matrix<Complex<float> >& B,
  Complex<float> beta,        Matrix<Complex<float> >& C );
template void Trrk
( UpperOrLower uplo, 
  Orientation orientationOfA, Orientation orientationOfB,
  Complex<float> alpha, const DistMatrix<Complex<float> >& A, 
                        const DistMatrix<Complex<float> >& B,
  Complex<float> beta,        DistMatrix<Complex<float> >& C );
template void LocalTrrk
( UpperOrLower uplo,
  Complex<float> alpha, const DistMatrix<Complex<float>,MC,  STAR>& A,
                        const DistMatrix<Complex<float>,STAR,MR  >& B,
  Complex<float> beta,        DistMatrix<Complex<float>,MC,  MR  >& C );
template void LocalTrrk
( UpperOrLower uplo,
  Orientation orientationOfB,
  Complex<float> alpha, const DistMatrix<Complex<float>,MC,STAR>& A,
                        const DistMatrix<Complex<float>,MR,STAR>& B,
  Complex<float> beta,        DistMatrix<Complex<float> >& C );
template void LocalTrrk
( UpperOrLower uplo,
  Orientation orientationOfA,
  Complex<float> alpha, const DistMatrix<Complex<float>,STAR,MC>& A,
                        const DistMatrix<Complex<float>,STAR,MR>& B,
  Complex<float> beta,        DistMatrix<Complex<float>,MC,  MR>& C );
template void LocalTrrk
( UpperOrLower uplo,
  Orientation orientationOfA, Orientation orientationOfB,
  Complex<float> alpha, const DistMatrix<Complex<float>,STAR,MC  >& A,
                        const DistMatrix<Complex<float>,MR,  STAR>& B,
  Complex<float> beta,        DistMatrix<Complex<float>,MC,  MR  >& E );
namespace internal {
template void TrrkNN
( UpperOrLower uplo,
  Complex<float> alpha, const Matrix<Complex<float> >& A, 
                        const Matrix<Complex<float> >& B,
  Complex<float> beta,        Matrix<Complex<float> >& C );
template void TrrkNT
( UpperOrLower uplo,
  Orientation orientationOfB,
  Complex<float> alpha, const Matrix<Complex<float> >& A, 
                        const Matrix<Complex<float> >& B,
  Complex<float> beta,        Matrix<Complex<float> >& C );
template void TrrkTN
( UpperOrLower uplo,
  Orientation orientationOfA,
  Complex<float> alpha, const Matrix<Complex<float> >& A, 
                        const Matrix<Complex<float> >& B,
  Complex<float> beta,        Matrix<Complex<float> >& C );
template void TrrkTT
( UpperOrLower uplo,
  Orientation orientationOfA, Orientation orientationOfB,
  Complex<float> alpha, const Matrix<Complex<float> >& A, 
                        const Matrix<Complex<float> >& B,
  Complex<float> beta,        Matrix<Complex<float> >& C );
} // namespace internal
#endif // ifndef DISABLE_COMPLEX
#endif // ifndef DISABLE_FLOAT

template void Trrk
( UpperOrLower uplo, 
  Orientation orientationOfA, Orientation orientationOfB,
  double alpha, const Matrix<double>& A, const Matrix<double>& B,
  double beta,        Matrix<double>& C );
template void Trrk
( UpperOrLower uplo, 
  Orientation orientationOfA, Orientation orientationOfB,
  double alpha, const DistMatrix<double>& A, const DistMatrix<double>& B,
  double beta,        DistMatrix<double>& C );
template void LocalTrrk
( UpperOrLower uplo,
  double alpha, const DistMatrix<double,MC,  STAR>& A,
                const DistMatrix<double,STAR,MR  >& B,
  double beta,        DistMatrix<double,MC,  MR  >& C );
template void LocalTrrk
( UpperOrLower uplo,
  Orientation orientationOfB,
  double alpha, const DistMatrix<double,MC,STAR>& A,
                const DistMatrix<double,MR,STAR>& B,
  double beta,        DistMatrix<double>& C );
template void LocalTrrk
( UpperOrLower uplo,
  Orientation orientationOfA,
  double alpha, const DistMatrix<double,STAR,MC>& A,
                const DistMatrix<double,STAR,MR>& B,
  double beta,        DistMatrix<double,MC,  MR>& C );
template void LocalTrrk
( UpperOrLower uplo,
  Orientation orientationOfA, Orientation orientationOfB,
  double alpha, const DistMatrix<double,STAR,MC  >& A,
                const DistMatrix<double,MR,  STAR>& B,
  double beta,        DistMatrix<double,MC,  MR  >& E );
namespace internal {
template void TrrkNN
( UpperOrLower uplo,
  double alpha, const Matrix<double>& A, const Matrix<double>& B,
  double beta,        Matrix<double>& C );
template void TrrkNT
( UpperOrLower uplo,
  Orientation orientationOfB,
  double alpha, const Matrix<double>& A, const Matrix<double>& B,
  double beta,        Matrix<double>& C );
template void TrrkTN
( UpperOrLower uplo,
  Orientation orientationOfA,
  double alpha, const Matrix<double>& A, const Matrix<double>& B,
  double beta,        Matrix<double>& C );
template void TrrkTT
( UpperOrLower uplo,
  Orientation orientationOfA, Orientation orientationOfB,
  double alpha, const Matrix<double>& A, const Matrix<double>& B,
  double beta,        Matrix<double>& C );
} // namespace internal
#ifndef DISABLE_COMPLEX
template void Trrk
( UpperOrLower uplo, 
  Orientation orientationOfA, Orientation orientationOfB,
  Complex<double> alpha, const Matrix<Complex<double> >& A, 
                         const Matrix<Complex<double> >& B,
  Complex<double> beta,        Matrix<Complex<double> >& C );
template void Trrk
( UpperOrLower uplo, 
  Orientation orientationOfA, Orientation orientationOfB,
  Complex<double> alpha, const DistMatrix<Complex<double> >& A, 
                         const DistMatrix<Complex<double> >& B,
  Complex<double> beta,        DistMatrix<Complex<double> >& C );
template void LocalTrrk
( UpperOrLower uplo,
  Complex<double> alpha, const DistMatrix<Complex<double>,MC,  STAR>& A,
                         const DistMatrix<Complex<double>,STAR,MR  >& B,
  Complex<double> beta,        DistMatrix<Complex<double>,MC,  MR  >& C );
template void LocalTrrk
( UpperOrLower uplo,
  Orientation orientationOfB,
  Complex<double> alpha, const DistMatrix<Complex<double>,MC,STAR>& A,
                         const DistMatrix<Complex<double>,MR,STAR>& B,
  Complex<double> beta,        DistMatrix<Complex<double> >& C );
template void LocalTrrk
( UpperOrLower uplo,
  Orientation orientationOfA,
  Complex<double> alpha, const DistMatrix<Complex<double>,STAR,MC>& A,
                         const DistMatrix<Complex<double>,STAR,MR>& B,
  Complex<double> beta,        DistMatrix<Complex<double>,MC,  MR>& C );
template void LocalTrrk
( UpperOrLower uplo,
  Orientation orientationOfA, Orientation orientationOfB,
  Complex<double> alpha, const DistMatrix<Complex<double>,STAR,MC  >& A,
                         const DistMatrix<Complex<double>,MR,  STAR>& B,
  Complex<double> beta,        DistMatrix<Complex<double>,MC,  MR  >& E );
namespace internal {
template void TrrkNN
( UpperOrLower uplo,
  Complex<double> alpha, const Matrix<Complex<double> >& A, 
                         const Matrix<Complex<double> >& B,
  Complex<double> beta,        Matrix<Complex<double> >& C );
template void TrrkNT
( UpperOrLower uplo,
  Orientation orientationOfB,
  Complex<double> alpha, const Matrix<Complex<double> >& A, 
                         const Matrix<Complex<double> >& B,
  Complex<double> beta,        Matrix<Complex<double> >& C );
template void TrrkTN
( UpperOrLower uplo,
  Orientation orientationOfA,
  Complex<double> alpha, const Matrix<Complex<double> >& A, 
                         const Matrix<Complex<double> >& B,
  Complex<double> beta,        Matrix<Complex<double> >& C );
template void TrrkTT
( UpperOrLower uplo,
  Orientation orientationOfA, Orientation orientationOfB,
  Complex<double> alpha, const Matrix<Complex<double> >& A, 
                         const Matrix<Complex<double> >& B,
  Complex<double> beta,        Matrix<Complex<double> >& C );
} // namespace internal
#endif // ifndef DISABLE_COMPLEX

} // namespace elem
