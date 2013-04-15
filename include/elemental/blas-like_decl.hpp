/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef BLAS_DECL_HPP
#define BLAS_DECL_HPP

namespace elem {

template<typename T>
void Trrk
( UpperOrLower uplo, 
  Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const Matrix<T>& A, const Matrix<T>& B,
  T beta,        Matrix<T>& C );
template<typename T>
void Trrk
( UpperOrLower uplo, 
  Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const DistMatrix<T>& A, const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C );
template<typename T>
void LocalTrrk
( UpperOrLower uplo,
  T alpha, const DistMatrix<T,MC,  STAR>& A,
           const DistMatrix<T,STAR,MR  >& B,
  T beta,        DistMatrix<T,MC,  MR  >& C );
template<typename T>
void LocalTrrk
( UpperOrLower uplo,
  Orientation orientationOfB,
  T alpha, const DistMatrix<T,MC,STAR>& A,
           const DistMatrix<T,MR,STAR>& B,
  T beta,        DistMatrix<T>& C );
template<typename T>
void LocalTrrk
( UpperOrLower uplo,
  Orientation orientationOfA,
  T alpha, const DistMatrix<T,STAR,MC>& A,
           const DistMatrix<T,STAR,MR>& B,
  T beta,        DistMatrix<T,MC,  MR>& C );
template<typename T>
void LocalTrrk
( UpperOrLower uplo,
  Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const DistMatrix<T,STAR,MC  >& A,
           const DistMatrix<T,MR,  STAR>& B,
  T beta,        DistMatrix<T,MC,  MR  >& C );

template<typename T>
void Trr2k
( UpperOrLower uplo, 
  Orientation orientationOfA, Orientation orientationOfB,
  Orientation orientationOfC, Orientation orientationOfD,
  T alpha, const Matrix<T>& A, const Matrix<T>& B,
           const Matrix<T>& C, const Matrix<T>& D,
  T beta,        Matrix<T>& E );
template<typename T>
void Trr2k
( UpperOrLower uplo, 
  Orientation orientationOfA, Orientation orientationOfB,
  Orientation orientationOfC, Orientation orientationOfD,
  T alpha, const DistMatrix<T>& A, const DistMatrix<T>& B,
           const DistMatrix<T>& C, const DistMatrix<T>& D,
  T beta,        DistMatrix<T>& E );

template<typename T>
void LocalTrr2k
( UpperOrLower uplo,
  T alpha, const DistMatrix<T,MC,  STAR>& A, const DistMatrix<T,STAR,MR>& B,
           const DistMatrix<T,MC,  STAR>& C, const DistMatrix<T,STAR,MR>& D,
  T beta,        DistMatrix<T,MC,  MR  >& E );
template<typename T>
void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfD,
  T alpha, const DistMatrix<T,MC,  STAR>& A, const DistMatrix<T,STAR,MR>& B,
           const DistMatrix<T,MC,  STAR>& C, const DistMatrix<T,MR,STAR>& D,
  T beta,        DistMatrix<T,MC,  MR  >& E  );
template<typename T>
void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfC,
  T alpha, const DistMatrix<T,MC,  STAR>& A, const DistMatrix<T,STAR,MR>& B,
           const DistMatrix<T,STAR,MC  >& C, const DistMatrix<T,STAR,MR>& D,
  T beta,        DistMatrix<T,MC,  MR  >& E  );
template<typename T>
void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfC,
  Orientation orientationOfD,
  T alpha, const DistMatrix<T,MC,  STAR>& A, const DistMatrix<T,STAR,MR>& B,
           const DistMatrix<T,STAR,MC  >& C, const DistMatrix<T,MR,STAR>& D,
  T beta,        DistMatrix<T,MC,  MR  >& E  );
template<typename T>
void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfB,
  T alpha, const DistMatrix<T,MC,  STAR>& A, const DistMatrix<T,MR,STAR>& B,
           const DistMatrix<T,MC,  STAR>& C, const DistMatrix<T,STAR,MR>& D,
  T beta,        DistMatrix<T,MC,  MR  >& E  );
template<typename T>
void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfB,
  Orientation orientationOfD,
  T alpha, const DistMatrix<T,MC,STAR>& A, const DistMatrix<T,MR,STAR>& B,
           const DistMatrix<T,MC,STAR>& C, const DistMatrix<T,MR,STAR>& D,
  T beta,        DistMatrix<T,MC,MR  >& E  );
template<typename T>
void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfB,
  Orientation orientationOfC,
  T alpha, const DistMatrix<T,MC,  STAR>& A, const DistMatrix<T,MR,STAR>& B,
           const DistMatrix<T,STAR,MC  >& C, const DistMatrix<T,STAR,MR>& D,
  T beta,        DistMatrix<T,MC,  MR  >& E  );
template<typename T>
void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfB,
  Orientation orientationOfC,
  Orientation orientationOfD,
  T alpha, const DistMatrix<T,MC,  STAR>& A, const DistMatrix<T,MR,STAR>& B,
           const DistMatrix<T,STAR,MC  >& C, const DistMatrix<T,MR,STAR>& D,
  T beta,        DistMatrix<T,MC,  MR  >& E  );
template<typename T>
void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfA,
  T alpha, const DistMatrix<T,STAR,MC  >& A, const DistMatrix<T,STAR,MR>& B,
           const DistMatrix<T,MC,  STAR>& C, const DistMatrix<T,STAR,MR>& D,
  T beta,        DistMatrix<T,MC,  MR  >& E  );
template<typename T>
void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfD,
  T alpha, const DistMatrix<T,STAR,MC  >& A, const DistMatrix<T,STAR,MR>& B,
           const DistMatrix<T,MC,  STAR>& C, const DistMatrix<T,MR,STAR>& D,
  T beta,        DistMatrix<T,MC,  MR  >& E  );
template<typename T>
void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfC,
  T alpha, const DistMatrix<T,STAR,MC>& A, const DistMatrix<T,STAR,MR>& B,
           const DistMatrix<T,STAR,MC>& C, const DistMatrix<T,STAR,MR>& D,
  T beta,        DistMatrix<T,MC,  MR>& E  );
template<typename T>
void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfC,
  Orientation orientationOfD,
  T alpha, const DistMatrix<T,STAR,MC>& A, const DistMatrix<T,STAR,MR>& B,
           const DistMatrix<T,STAR,MC>& C, const DistMatrix<T,MR,STAR>& D,
  T beta,        DistMatrix<T,MC,  MR>& E  );
template<typename T>
void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfB,
  T alpha, const DistMatrix<T,STAR,MC  >& A, const DistMatrix<T,MR,STAR>& B,
           const DistMatrix<T,MC,  STAR>& C, const DistMatrix<T,STAR,MR>& D,
  T beta,        DistMatrix<T,MC,  MR  >& E  );
template<typename T>
void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfB,
  Orientation orientationOfD,
  T alpha, const DistMatrix<T,STAR,MC  >& A, const DistMatrix<T,MR,STAR>& B,
           const DistMatrix<T,MC,  STAR>& C, const DistMatrix<T,MR,STAR>& D,
  T beta,        DistMatrix<T,MC,  MR  >& E  );
template<typename T>
void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfB,
  Orientation orientationOfC,
  T alpha, const DistMatrix<T,STAR,MC>& A, const DistMatrix<T,MR,STAR>& B,
           const DistMatrix<T,STAR,MC>& C, const DistMatrix<T,STAR,MR>& D,
  T beta,        DistMatrix<T,MC,  MR>& E  );
template<typename T>
void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfB,
  Orientation orientationOfC,
  Orientation orientationOfD,
  T alpha, const DistMatrix<T,STAR,MC>& A, const DistMatrix<T,MR,STAR>& B,
           const DistMatrix<T,STAR,MC>& C, const DistMatrix<T,MR,STAR>& D,
  T beta,        DistMatrix<T,MC,  MR>& E  );


namespace internal {

template<typename T>
void TrrkNN
( UpperOrLower uplo,
  T alpha, const Matrix<T>& A, const Matrix<T>& B,
  T beta,        Matrix<T>& C );
template<typename T>
void TrrkNT
( UpperOrLower uplo,
  Orientation orientationOfB,
  T alpha, const Matrix<T>& A, const Matrix<T>& B,
  T beta,        Matrix<T>& C );
template<typename T>
void TrrkTN
( UpperOrLower uplo,
  Orientation orientationOfA,
  T alpha, const Matrix<T>& A, const Matrix<T>& B,
  T beta,        Matrix<T>& C );
template<typename T>
void TrrkTT
( UpperOrLower uplo,
  Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const Matrix<T>& A, const Matrix<T>& B,
  T beta,        Matrix<T>& C );

} // namespace internal

//----------------------------------------------------------------------------//
// Tuning parameters                                                          //
//----------------------------------------------------------------------------//

template<typename T> void SetLocalSymvBlocksize( int blocksize );
template<> void SetLocalSymvBlocksize<float>( int blocksize );
template<> void SetLocalSymvBlocksize<double>( int blocksize );
template<> void SetLocalSymvBlocksize<Complex<float> >( int blocksize );
template<> void SetLocalSymvBlocksize<Complex<double> >( int blocksize );

template<typename T> void SetLocalTrrkBlocksize( int blocksize );
template<> void SetLocalTrrkBlocksize<float>( int blocksize );
template<> void SetLocalTrrkBlocksize<double>( int blocksize );
template<> void 
SetLocalTrrkBlocksize<Complex<float> >( int blocksize );
template<> void 
SetLocalTrrkBlocksize<Complex<double> >( int blocksize );

template<typename T> void SetLocalTrr2kBlocksize( int blocksize );
template<> void SetLocalTrr2kBlocksize<float>( int blocksize );
template<> void SetLocalTrr2kBlocksize<double>( int blocksize );
template<> void SetLocalTrr2kBlocksize<Complex<float> >( int blocksize );
template<> void SetLocalTrr2kBlocksize<Complex<double> >( int blocksize );

template<typename T> int LocalSymvBlocksize();
template<> int LocalSymvBlocksize<float>();
template<> int LocalSymvBlocksize<double>();
template<> int LocalSymvBlocksize<scomplex>();
template<> int LocalSymvBlocksize<dcomplex>();

template<typename T> int LocalTrrkBlocksize();
template<> int LocalTrrkBlocksize<float>();
template<> int LocalTrrkBlocksize<double>();
template<> int LocalTrrkBlocksize<scomplex>();
template<> int LocalTrrkBlocksize<dcomplex>();

template<typename T> int LocalTrr2kBlocksize();
template<> int LocalTrr2kBlocksize<float>();
template<> int LocalTrr2kBlocksize<double>();
template<> int LocalTrr2kBlocksize<scomplex>();
template<> int LocalTrr2kBlocksize<dcomplex>();

} // namespace elem

#endif // ifndef BLAS_DECL_HPP
