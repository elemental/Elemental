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
#ifndef ELEMENTAL_BASIC_HPP
#define ELEMENTAL_BASIC_HPP 1

#include "elemental/core/partitioning.hpp"

namespace elem {

//----------------------------------------------------------------------------//
// Tuning parameters                                                          //
//----------------------------------------------------------------------------//

template<typename T> void SetLocalHemvBlocksize( int blocksize );
template<> void SetLocalHemvBlocksize<float>( int blocksize );
template<> void SetLocalHemvBlocksize<double>( int blocksize );
template<> void SetLocalHemvBlocksize<Complex<float> >( int blocksize );
template<> void SetLocalHemvBlocksize<Complex<double> >( int blocksize );

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

template<typename T> int LocalHemvBlocksize();
template<> int LocalHemvBlocksize<float>();
template<> int LocalHemvBlocksize<double>();
template<> int LocalHemvBlocksize<scomplex>();
template<> int LocalHemvBlocksize<dcomplex>();

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

//----------------------------------------------------------------------------//
// Level 1 BLAS-like functionality                                            //
//----------------------------------------------------------------------------//

//
// Axpy (Alpha X Plus Y):
//
// Y := alpha X + Y
//

// Serial version
template<typename T>
void
Axpy( T alpha, const Matrix<T>& X, Matrix<T>& Y );

// Parallel version
template<typename T, Distribution U, Distribution V>
void
Axpy( T alpha, const DistMatrix<T,U,V>& X, DistMatrix<T,U,V>& Y );

//
// Copy:
//
// Y := X
//

// Serial version
template<typename T>
void
Copy( const Matrix<T>& X, Matrix<T>& Y );

// Parallel version
template<typename T, 
         Distribution U, Distribution V,
         Distribution W, Distribution Z >
void
Copy( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B );

//
// DiagonalScale
//
// Performs either X := op(D) X or
//                 X := X op(D)

// Serial version
template<typename T>
void DiagonalScale
( Side side, Orientation orientation, 
  const Matrix<T>& d, Matrix<T>& X );

// Parallel version
template<typename T,
         Distribution U,Distribution V,
         Distribution W,Distribution Z>
void DiagonalScale
( Side side, Orientation orientation, 
  const DistMatrix<T,U,V>& d, DistMatrix<T,W,Z>& X );

//
// DiagonalSolve
//
// Performs either X := op(D)^-1 X or 
//                 X := X op(D)^-1

// Serial version
template<typename F>
void DiagonalSolve
( Side side, Orientation orientation, 
  const Matrix<F>& d, Matrix<F>& X, bool checkIfSingular=false );

// Parallel version
template<typename F,
         Distribution U,Distribution V,
         Distribution W,Distribution Z>
void DiagonalSolve
( Side side, Orientation orientation, 
  const DistMatrix<F,U,V>& d, DistMatrix<F,W,Z>& X, 
  bool checkIfSingular=false );

//
// Dot: 
// 
// Returns (x,y) = x^H y.
//
// Though the standard BLAS interface only defines DOT for real 
// datatypes, it is naturally generalized to an inner product over the
// complex field. Recall that the conjugate symmetry of inner products 
// requires that (x,y) = conj(y,x), so that (x,x) = conj( (x,x) ) => 
// (x,x) is real. This requires that we choose (x,x) = conj(x)^T * x.
//

// Serial version
template<typename T>
T
Dot( const Matrix<T>& x, const Matrix<T>& y );

// Parallel version
template<typename T, 
         Distribution U, Distribution V,
         Distribution W, Distribution Z >
T
Dot( const DistMatrix<T,U,V>& x, const DistMatrix<T,W,Z>& y );

//
// Dotc:
//
// Returns (x,y) = x^H y.
//
// This is the sister routine to 'Dot'; while 'Dot' is originally defined 
// only over the reals, 'Dotc' was defined only over the complex field. 
// They both have been extended to the same function, so from our point of 
// view they are identical.
//

// Serial version
template<typename T>
T
Dotc( const Matrix<T>& x, const Matrix<T>& y );

// Parallel version
template<typename T, 
         Distribution U, Distribution V,
         Distribution W, Distribution Z >
T
Dotc( const DistMatrix<T,U,V>& x, const DistMatrix<T,W,Z>& y );

//
// Dotu: 
//
// Returns x^T y.
//
// Note: in the complex case, this is NOT an inner product.
//

// Serial version
template<typename T>
T
Dotu( const Matrix<T>& x, const Matrix<T>& y );

// Parallel version
template<typename T, Distribution U, Distribution V,
                     Distribution W, Distribution Z >
T
Dotu( const DistMatrix<T,U,V>& x, const DistMatrix<T,W,Z>& y );

//
// Nrm2 (2-norm):
//
// || x ||_2 = sqrt( x^H x ).
//

// Serial version for real datatypes
template<typename R>
R Nrm2( const Matrix<R>& x ); 

// Serial version for complex datatypes
template<typename R>
R Nrm2( const Matrix<Complex<R> >& x );

// Parallel version for real datatypes
template<typename R>
R Nrm2( const DistMatrix<R,MC,MR>& x );

// Parallel version for complex datatypes
template<typename R>
R Nrm2( const DistMatrix<Complex<R>, MC, MR >& x );

// 
// Scal:
//
// X := alpha X
//

// Serial version
template<typename T>
void Scal( T alpha, Matrix<T>& X );
    
// Parallel version
template<typename T, Distribution U, Distribution V>
void Scal( T alpha, DistMatrix<T,U,V>& X );
    
//----------------------------------------------------------------------------//
// Level 1 BLAS-like extensions                                               //
//----------------------------------------------------------------------------//

//
// Conjugate: 
//
// Conjugates a matrix. The in-place version performs A := Conjugate(A), 
// while the out-of-place sets B := Conjugate(A).
//

// In-place serial version for real datatypes. 
// Note: this is a no-op.
template<typename Z>
void Conjugate( Matrix<Z>& A );

// In-place serial version for complex datatypes.
template<typename Z>
void Conjugate( Matrix<Complex<Z> >& A );

// In-place parallel version
template<typename T, Distribution U, Distribution V>
void Conjugate( DistMatrix<T,U,V>& A );

// Out-of-place serial version.
template<typename T>
void Conjugate( const Matrix<T>& A, Matrix<T>& B );

// Out-of-place parallel version.
template<typename T, 
         Distribution U, Distribution V,
         Distribution W, Distribution Z>
void Conjugate( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B );

//
// Adjoint:
//
// B := A^H
//

// Serial version
template<typename T>
void Adjoint( const Matrix<T>& A, Matrix<T>& B );

// Parallel version
template<typename T, 
         Distribution U, Distribution V,
         Distribution W, Distribution Z>
void Adjoint( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B );

//
// Transpose:
//
// B := A^T
//

// Serial version
template<typename T>
void Transpose( const Matrix<T>& A, Matrix<T>& B );

// Parallel version
template<typename T, 
         Distribution U, Distribution V,
         Distribution W, Distribution Z>
void Transpose( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B );

//----------------------------------------------------------------------------//
// Level 2 BLAS-like routines                                                 //
//----------------------------------------------------------------------------//

//
// Gemv (GEneral Matrix-Vector multiply):
//
// y := alpha orientation( A ) x + beta y,
// where orientation( A ) is determined by 'orientation'.
//

// Serial version
template<typename T>
void Gemv
( Orientation orientation, 
  T alpha, const Matrix<T>& A, const Matrix<T>& x, T beta, Matrix<T>& y );

// Parallel version
template<typename T>
void Gemv
( Orientation orientationOfA,
  T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& x,
  T beta,        DistMatrix<T,MC,MR>& y );

//
// Ger (GEneral Rank-one update):
//
// A := alpha x y^H + A
//

// Serial version
template<typename T>
void Ger( T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A );

// Parallel version
template<typename T>
void Ger
( T alpha, const DistMatrix<T,MC,MR>& x, const DistMatrix<T,MC,MR>& y,
                 DistMatrix<T,MC,MR>& A );

//
// Gerc (GEneral Rank-one Conjugated update):
//
// A := alpha x y^H + A
//
// This is identical to Ger because both have been extended to work for both
// real and complex datatypes.
//

// Serial version
template<typename T>
void Gerc( T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A );

// Parallel version
template<typename T>
void Gerc
( T alpha, const DistMatrix<T,MC,MR>& x, const DistMatrix<T,MC,MR>& y,
                 DistMatrix<T,MC,MR>& A );

//
// Geru (GEneral Rank-one Unconjugated update):
//
// A := alpha x y^T + A
//

// Serial version
template<typename T>
void Geru( T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A );

// Parallel version
template<typename T>
void Geru
( T alpha, const DistMatrix<T,MC,MR>& x, const DistMatrix<T,MC,MR>& y,
                 DistMatrix<T,MC,MR>& A );

//
// Hemv (HErmitian Matrix-Vector multiply):
//
// Implicitly performs
//   y := alpha A x + beta y,
// where only the triangle specified by 'uplo' is referenced and the other
// triangle is implied by the Hermitian assumption.

// Serial version
template<typename T>
void Hemv
( UpperOrLower uplo,
  T alpha, const Matrix<T>& A, const Matrix<T>& x, T beta, Matrix<T>& y );

// Parallel version
template<typename T>
void Hemv
( UpperOrLower uplo,
  T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& x,
  T beta,        DistMatrix<T,MC,MR>& y );

//
// Her (HErmitian Rank-one update):
//
// Implicitly performs
//   A := alpha x x^H + A,
// where only the triangle specified by 'uplo' is updated.
//

// Serial version
template<typename T>
void Her( UpperOrLower uplo, T alpha, const Matrix<T>& x, Matrix<T>& A );

// Parallel version
template<typename T>
void Her
( UpperOrLower uplo, 
  T alpha, const DistMatrix<T,MC,MR>& x, DistMatrix<T,MC,MR>& A );

//
// Her2 (HErmitian Rank-2 update):
//
// Implicitly performs
//   A := alpha ( x y^H + y x^H ) + A,
// where only the triangle specified by 'uplo' is updated.
//

// Serial version
template<typename T>
void Her2
( UpperOrLower uplo, 
  T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A );

// Parallel version
template<typename T>
void Her2
( UpperOrLower uplo,
  T alpha, const DistMatrix<T,MC,MR>& x, const DistMatrix<T,MC,MR>& y,
                 DistMatrix<T,MC,MR>& A );

//
// Symv (SYmmetric Matrix-Vector multiply):
//
// Implicitly performs
//   y := alpha A x + beta y,
// where only the triangle specified by 'uplo' is referenced and the other
// triangle is implied by the symmetry assumption.
//

// Serial version
template<typename T>
void Symv
( UpperOrLower uplo,
  T alpha, const Matrix<T>& A, const Matrix<T>& x, T beta, Matrix<T>& y );

// Parallel version
template<typename T>
void Symv
( UpperOrLower uplo,
  T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& x,
  T beta,        DistMatrix<T,MC,MR>& y );

//
// Syr (SYmmetric Rank-one update):
//
// Implicitly performs the update
//   A := alpha x x^T + A,
// where only the triangle specified by 'uplo' is updated.
//

// Serial version
template<typename T>
void Syr( UpperOrLower uplo, T alpha, const Matrix<T>& x, Matrix<T>& A );

// Parallel version
template<typename T>
void Syr
( UpperOrLower uplo,
  T alpha, const DistMatrix<T,MC,MR>& x, DistMatrix<T,MC,MR>& A );

//
// Syr2 (SYmmetric Rank-2 update):
//
// Implicitly perform the update
//   A := alpha ( x y^T + y x^T ) + A
// where only the triangle specified by 'uplo' is updated.
//

// Serial version
template<typename T>
void Syr2
( UpperOrLower uplo, 
  T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A );

// Parallel version
template<typename T>
void Syr2
( UpperOrLower uplo,
  T alpha, const DistMatrix<T,MC,MR>& x, const DistMatrix<T,MC,MR>& y,
                 DistMatrix<T,MC,MR>& A );

//
// Trmv (TRiangular Matrix-Vector multiply):
//
// Performs the update
//   x := orientation( A ) x,
// where 'uplo' determines whether or not A is to be implicitly treated as 
// lower or upper triangular, and 'diagonal' specifies whether it has an 
// implicit unit diagonal.
//

// Serial version
template<typename T>
void Trmv
( UpperOrLower uplo, Orientation orientation, Diagonal diagonal,
  const Matrix<T>& A, Matrix<T>& x );

// Parallel version
template<typename T>
void Trmv
( UpperOrLower uplo, Orientation orientation, Diagonal diagonal,
  const DistMatrix<T,MC,MR>& A, DistMatrix<T,MC,MR>& x );

//
// Trsv (TRiangular Solve with a Vector):
//
// Performs the update
//   x := orientation( A )^-1 x,
// where 'uplo' determines whether or not A is to be implicitly treated as 
// lower or upper triangular, and 'diagonal' specifies whether it has an 
// implicit unit diagonal.
//

// Serial version
template<typename F>
void Trsv
( UpperOrLower uplo, Orientation orientation, Diagonal diagonal,
  const Matrix<F>& A, Matrix<F>& x );

// Parallel version
template<typename F>
void Trsv
( UpperOrLower uplo, Orientation orientation, Diagonal diagonal,
  const DistMatrix<F,MC,MR>& A, DistMatrix<F,MC,MR>& x );

//----------------------------------------------------------------------------//
// Level 3 BLAS-like routines                                                 //
//----------------------------------------------------------------------------//

//
// Gemm (GEneral Matrix-Matrix multiplication):
//
// C := alpha orientationOfA( A ) orientationOfB( B ) + beta C
//

// Serial version
template<typename T>
void Gemm
( Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const Matrix<T>& A, const Matrix<T>& B, T beta, Matrix<T>& C );

// Parallel version
template<typename T>
void Gemm
( Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C );

//
// Hemm (HErmitian Matrix-Matrix multiplication):
//
// Performs the update
//   C := alpha A B + beta C,  { side = LEFT }
// or
//   C := alpha B A + beta C,  { side = RIGHT }
// where only the triangle of 'A' specified by 'uplo' is referenced, and the
// other triangle is implied by the Hermitian assumption.
//

// Serial version
template<typename T>
void Hemm
( Side side, UpperOrLower uplo,
  T alpha, const Matrix<T>& A, const Matrix<T>& B, T beta, Matrix<T>& C );

// Parallel version
template<typename T>
void Hemm
( Side side, UpperOrLower uplo,
  T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C );

//
// Her2k (HErmitian Rank-2K update):
//
// Performs the update
//   C := alpha ( A B^H + B A^H ) + beta C, { orientation = NORMAL }
// or
//   C := alpha ( A^H B + B^H A ) + beta C, { orientation = ADJOINT }
// where only the triangle of C specified by 'uplo' is updated.
//

// Serial version
template<typename T>
void Her2k
( UpperOrLower uplo, Orientation orientation,
  T alpha, const Matrix<T>& A, const Matrix<T>& B, T beta, Matrix<T>& C );

// Parallel version
template<typename T>
void Her2k
( UpperOrLower uplo, Orientation orientation,
  T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C );

//
// Herk (HErmitian Rank-K update):
//
// Performs the update
//   C := alpha A B^H + beta C,  { orientation = NORMAL }
// or
//   C := alpha A^H B + beta C,  { orientation = ADJOINT }
// where only the triangle of C specified by 'uplo' is updated.
//

// Serial version
template<typename T>
void Herk
( UpperOrLower uplo, Orientation orientation,
  T alpha, const Matrix<T>& A, T beta, Matrix<T>& C );

// Parallel version
template<typename T>
void Herk
( UpperOrLower uplo, Orientation orientation,
  T alpha, const DistMatrix<T,MC,MR>& A, T beta, DistMatrix<T,MC,MR>& C );

//
// Hetrmm (HErmitian TRiangular Matrix-Matrix multiply):
//
// Either L := tril(L' L) or U := triu(U U')
//
// NOTE: This is not a standard BLAS routine and is in fact the same operation
//       as the LAPACK routine ?LAUUM. Since it is similar in spirit to many
//       other BLAS routines, I have instead placed it here.
//

// Serial version
template<typename T>
void Hetrmm( UpperOrLower uplo, Matrix<T>& A );

// Parallel version
template<typename T>
void Hetrmm( UpperOrLower uplo, DistMatrix<T,MC,MR>& A );

//
// Symm (SYmmetric Matrix-Matrix multiplication):
//
// Performs the update
//   C := alpha A B + beta C,  { side = LEFT }
// or
//   C := alpha B A + beta C,  { side = RIGHT }
// where only the triangle of A specified by 'uplo' is referenced, and the 
// other triangle is implied by the symmetry assumption.
//

// Serial version
template<typename T>
void Symm
( Side side, UpperOrLower uplo,
  T alpha, const Matrix<T>& A, const Matrix<T>& B, T beta, Matrix<T>& C ); 

// Parallel version
template<typename T>
void Symm
( Side side, UpperOrLower uplo,
  T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C );

//
// Syr2k (SYmmetric Rank-2K update):
//
// Performs the update
//   C := alpha ( A B^H + B A^H ) + beta C,  { orientation = NORMAL }
// or
//   C := alpha ( A^H B + B^H A ) + beta C,  { orientation = TRANSPOSE }
// where only the triangle of C specified by 'uplo' is updated.
//

// Serial version
template<typename T>
void Syr2k
( UpperOrLower uplo, Orientation orientation,
  T alpha, const Matrix<T>& A, const Matrix<T>& B, T beta, Matrix<T>& C );

// Parallel version
template<typename T>
void Syr2k
( UpperOrLower uplo, Orientation orientation,
  T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C );

//
// Syrk (SYmmetric Rank-K update):
//
// Performs the update
//   C := alpha A B^H + beta C,  { orientation = NORMAL }
// or
//   C := alpha A^H B + beta C,  { orientation = TRANSPOSE }
// where only the triangle of C specified by 'uplo' is updated.
//

// Serial version
template<typename T>
void Syrk
( UpperOrLower uplo, Orientation orientation,
  T alpha, const Matrix<T>& A, T beta, Matrix<T>& C );

// Parallel version
template<typename T>
void Syrk
( UpperOrLower uplo, Orientation orientation,
  T alpha, const DistMatrix<T,MC,MR>& A, T beta, DistMatrix<T,MC,MR>& C );

//
// Trmm (TRiangular Matrix-Matrix multiplication):
//
// Performs the update
//   B := alpha orientation( A ) B,  { side = LEFT }
// or
//   B := alpha B orientation( A ),  { side = RIGHT }
// where 'uplo' determines whether A is assumed to be upper or lower 
// triangular and 'diagonal' determines whether A has an implicit unit
// diagonal.
//

// Serial version
template<typename T>
void Trmm
( Side side, UpperOrLower uplo, Orientation orientation, Diagonal diagonal,
  T alpha, const Matrix<T>& A, Matrix<T>& B );

// Parallel version
template<typename T>
void Trmm
( Side side, UpperOrLower uplo, Orientation orientation, Diagonal diagonal,
  T alpha, const DistMatrix<T,MC,MR>& A, DistMatrix<T,MC,MR>& B );

//
// Trr2k (TRiangular Rank-2K update):
//
// Performs the update:
//
//   E := alpha (op(A) op(B) + op(C) op(D)) + beta E
//
// where op(X) is determined by 'orientationOfX', and only the triangle of 
// X specified by 'uplo' is updated.
//

// TODO: Serial version

// Parallel version 
template<typename T>
void Trr2k
( UpperOrLower uplo, 
  Orientation orientationOfA, Orientation orientationOfB,
  Orientation orientationOfC, Orientation orientationOfD,
  T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B, 
           const DistMatrix<T,MC,MR>& C, const DistMatrix<T,MC,MR>& D,
  T beta,        DistMatrix<T,MC,MR>& E );

//
// Trrk (TRiangular Rank-K update):
//
// Performs the update:
//
//   C := alpha op(A) op(B) + beta C
//
// where op(A) and op(B) are respectively determined by 'orientationOfA'
// and 'orientationOfB'. Only the triangle of C specified by 'uplo' is updated.
//

// Serial version
template<typename T>
void Trrk
( UpperOrLower uplo, Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const Matrix<T>& A, const Matrix<T>& B, 
  T beta,        Matrix<T>& C );

// Parallel version 
template<typename T>
void Trrk
( UpperOrLower uplo, Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B, 
  T beta,        DistMatrix<T,MC,MR>& C );

//
// Trsm (TRiangular Solve with Multiple right-hand sides):
//
// Performs the update
//   B := alpha orientation( A )^-1 B,  { side = LEFT }
// or
//   B := alpha B orientation( A )^-1,  { side = RIGHT }
// where 'uplo' determines whether A is assumed to be upper or lower
// triangular and 'diagonal' determines whether A has an implicit unit
// diagonal.
//

// Serial version
template<typename F>
void Trsm
( Side side, UpperOrLower uplo, Orientation orientation, Diagonal diagonal,
  F alpha, const Matrix<F>& A, Matrix<F>& B, 
  bool checkIfSingular=false ); 
        
// Parallel version
template<typename F>
void Trsm
( Side side, UpperOrLower uplo, Orientation orientation, Diagonal diagonal,
  F alpha, const DistMatrix<F,MC,MR>& A, DistMatrix<F,MC,MR>& B,
  bool checkIfSingular=false );

} // namespace elem

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

#include "./basic/internal.hpp"
#include "./basic/level1.hpp"
#include "./basic/level2.hpp"
#include "./basic/level3.hpp"

namespace elem {

//----------------------------------------------------------------------------//
// Local BLAS-like routines: Level 1                                          //
//----------------------------------------------------------------------------//

template<typename T>
inline void
Axpy( T alpha, const Matrix<T>& X, Matrix<T>& Y )
{
#ifndef RELEASE
    PushCallStack("Axpy");
#endif
    // If X and Y are vectors, we can allow one to be a column and the other
    // to be a row. Otherwise we force X and Y to be the same dimension.
    if( (X.Height()==1 || X.Width()==1) && (Y.Height()==1 || Y.Width()==1) )
    {
        const unsigned XLength = ( X.Width()==1 ? X.Height() : X.Width() );
#ifndef RELEASE
        const unsigned YLength = ( Y.Width()==1 ? Y.Height() : Y.Width() );
        if( XLength != YLength )
            throw std::logic_error("Nonconformal Axpy");
#endif
        if( X.Width()==1 && Y.Width()==1 )
        {
            blas::Axpy
            ( XLength, alpha, X.LockedBuffer(0,0), 1, Y.Buffer(0,0), 1 );
        }
        else if( X.Width()==1 )
        {
            blas::Axpy
            ( XLength, alpha, X.LockedBuffer(0,0), 1, Y.Buffer(0,0), Y.LDim() );
        }
        else if( Y.Width()==1 )
        {
            blas::Axpy
            ( XLength, alpha, X.LockedBuffer(0,0), X.LDim(), Y.Buffer(0,0), 1 );
        }
        else
        {
            blas::Axpy
            ( XLength, alpha, 
              X.LockedBuffer(0,0), X.LDim(), Y.Buffer(0,0), Y.LDim() );
        }
    }
    else 
    {
#ifndef RELEASE
        if( X.Height() != Y.Height() || X.Width() != Y.Width() )
            throw std::logic_error("Nonconformal Axpy");
#endif
        if( X.Width() <= X.Height() )
        {
            for( int j=0; j<X.Width(); ++j )
            {
                blas::Axpy
                ( X.Height(), alpha, X.LockedBuffer(0,j), 1, Y.Buffer(0,j), 1 );
            }
        }
        else
        {
            for( int i=0; i<X.Height(); ++i )
            {
                blas::Axpy
                ( X.Width(), alpha, X.LockedBuffer(i,0), X.LDim(),
                                    Y.Buffer(i,0),       Y.LDim() );
            }
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Copy( const Matrix<T>& A, Matrix<T>& B )
{
#ifndef RELEASE
    PushCallStack("Copy");
#endif
    B = A;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void 
DiagonalScale
( Side side, Orientation orientation, const Matrix<T>& d, Matrix<T>& X )
{
#ifndef RELEASE
    PushCallStack("DiagonalScale");
#endif
    const int m = X.Height();
    const int n = X.Width();
    const int ldim = X.LDim();
    if( side == LEFT )
    {
        for( int i=0; i<m; ++i )
        {
            const T delta = d.Get(i,0);
            T* XBuffer = X.Buffer(i,0);
            if( orientation == ADJOINT )
                for( int j=0; j<n; ++j )
                    XBuffer[j*ldim] *= Conj(delta);
            else
                for( int j=0; j<n; ++j )
                    XBuffer[j*ldim] *= delta;
        }
    }
    else
    {
        for( int j=0; j<n; ++j )
        {
            const T delta = d.Get(j,0);
            T* XBuffer = X.Buffer(0,j);
            if( orientation == ADJOINT )
                for( int i=0; i<m; ++i )
                    XBuffer[i] *= Conj(delta);
            else
                for( int i=0; i<m; ++i )
                    XBuffer[i] *= delta;
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void 
DiagonalSolve
( Side side, Orientation orientation, const Matrix<F>& d, Matrix<F>& X,
  bool checkIfSingular )
{
#ifndef RELEASE
    PushCallStack("DiagonalSolve");
#endif
    const int m = X.Height();
    const int n = X.Width();
    const int ldim = X.LDim();
    if( side == LEFT )
    {
        for( int i=0; i<m; ++i )
        {
            const F delta = d.Get(i,0);
            if( checkIfSingular && delta == (F)0 )
                throw SingularMatrixException();
            const F deltaInv = static_cast<F>(1)/delta;
            F* XBuffer = X.Buffer(i,0);
            if( orientation == ADJOINT )
                for( int j=0; j<n; ++j )
                    XBuffer[j*ldim] *= Conj(deltaInv);
            else
                for( int j=0; j<n; ++j )
                    XBuffer[j*ldim] *= deltaInv;
        }
    }
    else
    {
        for( int j=0; j<n; ++j )
        {
            const F delta = d.Get(j,0);
            if( checkIfSingular && delta == (F)0 )
                throw SingularMatrixException();
            const F deltaInv = static_cast<F>(1)/delta;
            F* XBuffer = X.Buffer(0,j);
            if( orientation == ADJOINT )
                for( int i=0; i<m; ++i )
                    XBuffer[i] *= Conj(deltaInv);
            else
                for( int i=0; i<m; ++i )
                    XBuffer[i] *= deltaInv;
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline T
Dot( const Matrix<T>& x, const Matrix<T>& y )
{
#ifndef RELEASE
    PushCallStack("Dot");
    if( (x.Height() != 1 && x.Width() != 1) ||
        (y.Height() != 1 && y.Width() != 1) )
        throw std::logic_error("Expected vector inputs");
    int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
    int yLength = ( y.Width() == 1 ? y.Height() : y.Width() );
    if( xLength != yLength )
        throw std::logic_error("x and y must be the same length");
#endif
    T dotProduct;
    if( x.Width() == 1 && y.Width() == 1 )
    {
        dotProduct = blas::Dot
                     ( x.Height(), x.LockedBuffer(), 1,
                                   y.LockedBuffer(), 1 );
    }
    else if( x.Width() == 1 )
    {
        dotProduct = blas::Dot
                     ( x.Height(), x.LockedBuffer(), 1,
                                   y.LockedBuffer(), y.LDim() );
    }
    else if( y.Width() == 1 )
    {
        dotProduct = blas::Dot
                     ( x.Width(), x.LockedBuffer(), x.LDim(),
                                  y.LockedBuffer(), 1        );
    }
    else
    {
        dotProduct = blas::Dot
                     ( x.Width(), x.LockedBuffer(), x.LDim(),
                                  y.LockedBuffer(), y.LDim() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return dotProduct;
}

template<typename T>
inline T
Dotc( const Matrix<T>& x, const Matrix<T>& y )
{
#ifndef RELEASE
    PushCallStack("Dotc");
    if( (x.Height() != 1 && x.Width() != 1) ||
        (y.Height() != 1 && y.Width() != 1) )
        throw std::logic_error("Expected vector inputs");
    int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
    int yLength = ( y.Width() == 1 ? y.Height() : y.Width() );
    if( xLength != yLength )
        throw std::logic_error("x and y must be the same length");
#endif
    T dotProduct;
    if( x.Width() == 1 && y.Width() == 1 )
    {
        dotProduct = blas::Dotc
                     ( x.Height(), x.LockedBuffer(), 1,
                                   y.LockedBuffer(), 1 );
    }
    else if( x.Width() == 1 )
    {
        dotProduct = blas::Dotc
                     ( x.Height(), x.LockedBuffer(), 1,
                                   y.LockedBuffer(), y.LDim() );
    }
    else if( y.Width() == 1 )
    {
        dotProduct = blas::Dotc
                     ( x.Width(), x.LockedBuffer(), x.LDim(),
                                  y.LockedBuffer(), 1        );
    }
    else
    {
        dotProduct = blas::Dotc
                     ( x.Width(), x.LockedBuffer(), x.LDim(),
                                  y.LockedBuffer(), y.LDim() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return dotProduct;
}

template<typename T>
inline T
Dotu( const Matrix<T>& x, const Matrix<T>& y )
{
#ifndef RELEASE
    PushCallStack("Dotu");
    if( (x.Height() != 1 && x.Width() != 1) ||
        (y.Height() != 1 && y.Width() != 1) )
        throw std::logic_error("Expected vector inputs");
    int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
    int yLength = ( y.Width() == 1 ? y.Height() : y.Width() );
    if( xLength != yLength )
        throw std::logic_error("x and y must be the same length");
#endif
    T dotProduct;
    if( x.Width() == 1 && y.Width() == 1 )
    {
        dotProduct = blas::Dotu
                     ( x.Height(), x.LockedBuffer(), 1,
                                   y.LockedBuffer(), 1 );
    }
    else if( x.Width() == 1 )
    {
        dotProduct = blas::Dotu
                     ( x.Height(), x.LockedBuffer(), 1,
                                   y.LockedBuffer(), y.LDim() );
    }
    else if( y.Width() == 1 )
    {
        dotProduct = blas::Dotu
                     ( x.Width(), x.LockedBuffer(), x.LDim(),
                                  y.LockedBuffer(), 1        );
    }
    else
    {
        dotProduct = blas::Dotu
                     ( x.Width(), x.LockedBuffer(), x.LDim(),
                                  y.LockedBuffer(), y.LDim() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return dotProduct;
}

template<typename R>
inline R
Nrm2( const Matrix<R>& x )
{
#ifndef RELEASE
    PushCallStack("Nrm2");
    if( x.Height() != 1 && x.Width() != 1 )
        throw std::logic_error("Expected vector input");
#endif
    R norm;
    if( x.Width() == 1 )
        norm = blas::Nrm2( x.Height(), x.LockedBuffer(), 1 );
    else
        norm = blas::Nrm2( x.Width(), x.LockedBuffer(), x.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

template<typename R>
inline R
Nrm2( const Matrix<Complex<R> >& x )
{
#ifndef RELEASE
    PushCallStack("Nrm2");
    if( x.Height() != 1 && x.Width() != 1 )
        throw std::logic_error("Expected vector input");
#endif
    R norm;
    if( x.Width() == 1 )
        norm = blas::Nrm2( x.Height(), x.LockedBuffer(), 1 );
    else
        norm = blas::Nrm2( x.Width(), x.LockedBuffer(), x.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

template<typename T>
inline void
Scal( T alpha, Matrix<T>& X )
{
#ifndef RELEASE
    PushCallStack("Scal");
#endif
    if( alpha != (T)1 )
    {
        if( alpha == (T)0 )
            for( int j=0; j<X.Width(); ++j )
                for( int i=0; i<X.Height(); ++i )
                    X.Set(i,j,0);
        else
            for( int j=0; j<X.Width(); ++j )
                blas::Scal( X.Height(), alpha, X.Buffer(0,j), 1 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

//----------------------------------------------------------------------------//
// Local BLAS-like routines: Level 1 (extensions)                             //
//----------------------------------------------------------------------------//

// Default case is for real datatypes
template<typename Z>
inline void
Conjugate( Matrix<Z>& A )
{ }

// Specialization is to complex datatypes
template<typename Z>
inline void
Conjugate( Matrix<Complex<Z> >& A )
{
#ifndef RELEASE
    PushCallStack("Conjugate (in-place)");
#endif
    const int m = A.Height();
    const int n = A.Width();
    for( int j=0; j<n; ++j )
        for( int i=0; i<m; ++i )
            A.Set(i,j,Conj(A.Get(i,j)));
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Conjugate( const Matrix<T>& A, Matrix<T>& B )
{
#ifndef RELEASE
    PushCallStack("Conjugate");
#endif
    const int m = A.Height();
    const int n = A.Width();
    B.ResizeTo( m, n );
    for( int j=0; j<n; ++j )
        for( int i=0; i<m; ++i )
            B.Set(i,j,Conj(A.Get(i,j)));
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Adjoint( const Matrix<T>& A, Matrix<T>& B )
{
#ifndef RELEASE
    PushCallStack("Adjoint");
#endif
    const int m = A.Height();
    const int n = A.Width();
    if( !B.Viewing() )
        B.ResizeTo( n, m );
    else if( B.Height() != n || B.Width() != m )
        throw std::logic_error
        ("If Adjoint'ing into a view, it must be the right size");
    for( int j=0; j<n; ++j )
        for( int i=0; i<m; ++i )
            B.Set(j,i,Conj(A.Get(i,j)));
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Transpose( const Matrix<T>& A, Matrix<T>& B )
{
#ifndef RELEASE
    PushCallStack("Transpose");
#endif
    const int m = A.Height();
    const int n = A.Width();
    if( !B.Viewing() )
        B.ResizeTo( n, m );
    else if( B.Height() != n || B.Width() != m )
        throw std::logic_error
        ("If Transposing into a view, it must be the right size");
    for( int j=0; j<n; ++j )
        for( int i=0; i<m; ++i )
            B.Set(j,i,A.Get(i,j));
#ifndef RELEASE
    PopCallStack();
#endif
}

//----------------------------------------------------------------------------//
// Local BLAS-like routines: Level 2                                          //
//----------------------------------------------------------------------------//

template<typename T>
inline void
Gemv
( Orientation orientation,
  T alpha, const Matrix<T>& A, const Matrix<T>& x, T beta, Matrix<T>& y )
{
#ifndef RELEASE
    PushCallStack("Gemv");
    if( ( x.Height() != 1 && x.Width() != 1 ) ||
        ( y.Height() != 1 && y.Width() != 1 ) )
    {
        std::ostringstream msg;
        msg << "x and y must be vectors:\n"
            << "  x ~ " << x.Height() << " x " << x.Width() << "\n"
            << "  y ~ " << y.Height() << " x " << y.Width();
        throw std::logic_error( msg.str().c_str() );
    }
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    const int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
    if( orientation == NORMAL )
    {
        if( A.Height() != yLength || A.Width() != xLength )
        {
            std::ostringstream msg;
            msg << "A must conform with x and y:\n"
                << "  A ~ " << A.Height() << " x " << A.Width() << "\n"
                << "  x ~ " << x.Height() << " x " << x.Width() << "\n"
                << "  y ~ " << y.Height() << " x " << y.Width();
            throw std::logic_error( msg.str().c_str() );
        }
    }
    else
    {
        if( A.Width() != yLength || A.Height() != xLength )
        {
            std::ostringstream msg;
            msg << "A must conform with x and y:\n"
                << "  A ~ " << A.Height() << " x " << A.Width() << "\n"
                << "  x ~ " << x.Height() << " x " << x.Width() << "\n"
                << "  y ~ " << y.Height() << " x " << y.Width();
            throw std::logic_error( msg.str().c_str() );
        }
    }
#endif
    const char transChar = OrientationToChar( orientation );
    const int m = A.Height();
    const int n = A.Width();
    const int k = ( transChar == 'N' ? n : m );
    const int incx = ( x.Width()==1 ? 1 : x.LDim() );
    const int incy = ( y.Width()==1 ? 1 : y.LDim() );
    if( k != 0 )
    {
        blas::Gemv
        ( transChar, m, n, 
          alpha, A.LockedBuffer(), A.LDim(), x.LockedBuffer(), incx, 
          beta,  y.Buffer(), incy );
    }
    else
    {
        Scal( beta, y );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Ger( T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("Ger");
    if( ( x.Height() != 1 && x.Width() != 1 ) ||
        ( y.Height() != 1 && y.Width() != 1 ) )
        throw std::logic_error("x and y must be vectors");
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    const int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
    if( xLength != A.Height() || yLength != A.Width() )
    {
        std::ostringstream msg;
        msg << "Nonconformal Ger:\n"
            << "  x ~ " << x.Height() << " x " << x.Width() << "\n"
            << "  y ~ " << y.Height() << " x " << y.Width() << "\n"
            << "  A ~ " << A.Height() << " x " << A.Width();
        throw std::logic_error( msg.str().c_str() );
    }
#endif
    const int m = A.Height(); 
    const int n = A.Width();
    const int incx = ( x.Width()==1 ? 1 : x.LDim() );
    const int incy = ( y.Width()==1 ? 1 : y.LDim() );
    blas::Ger
    ( m, n, alpha, x.LockedBuffer(), incx, y.LockedBuffer(), incy, 
                   A.Buffer(), A.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Gerc( T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("Gerc");
    if( ( x.Height() != 1 && x.Width() != 1 ) ||
        ( y.Height() != 1 && y.Width() != 1 ) )
        throw std::logic_error("x and y must be vectors");
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    const int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
    if( xLength != A.Height() || yLength != A.Width() )
        throw std::logic_error("Nonconformal Gerc");
#endif
    const int m = A.Height(); 
    const int n = A.Width();
    const int incx = ( x.Width()==1 ? 1 : x.LDim() );
    const int incy = ( y.Width()==1 ? 1 : y.LDim() );
    blas::Gerc
    ( m, n, alpha, x.LockedBuffer(), incx, y.LockedBuffer(), incy, 
                   A.Buffer(), A.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Geru( T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("Geru");
    if( ( x.Height() != 1 && x.Width() != 1 ) ||
        ( y.Height() != 1 && y.Width() != 1 ) )
        throw std::logic_error("x and y must be vectors");
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    const int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
    if( xLength != A.Height() || yLength != A.Width() )
        throw std::logic_error("Nonconformal Geru");
#endif
    const int m = A.Height(); 
    const int n = A.Width();
    const int incx = ( x.Width()==1 ? 1 : x.LDim() );
    const int incy = ( y.Width()==1 ? 1 : y.LDim() );
    blas::Geru
    ( m, n, alpha, x.LockedBuffer(), incx, y.LockedBuffer(), incy, 
                   A.Buffer(), A.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Hemv
( UpperOrLower uplo,
  T alpha, const Matrix<T>& A, const Matrix<T>& x, T beta, Matrix<T>& y )
{
#ifndef RELEASE
    PushCallStack("Hemv");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    if( ( x.Height() != 1 && x.Width() != 1 ) ||
        ( y.Height() != 1 && y.Width() != 1 ) )
        throw std::logic_error("x and y must be vectors");
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    const int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
    if( A.Height() != xLength || A.Height() != yLength )
        throw std::logic_error("A must conform with x and y");
#endif
    const char uploChar = UpperOrLowerToChar( uplo );
    const int m = A.Height();
    const int incx = ( x.Width()==1 ? 1 : x.LDim() );
    const int incy = ( y.Width()==1 ? 1 : y.LDim() );
    blas::Hemv
    ( uploChar, m,
      alpha, A.LockedBuffer(), A.LDim(), x.LockedBuffer(), incx,
      beta,  y.Buffer(), incy );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Her( UpperOrLower uplo, T alpha, const Matrix<T>& x, Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("Her");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    if( x.Width() != 1 && x.Height() != 1 )
        throw std::logic_error("x must be a vector");
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    if( xLength != A.Height() )
        throw std::logic_error("x must conform with A");
#endif
    const char uploChar = UpperOrLowerToChar( uplo );
    const int m = A.Height();
    const int incx = ( x.Width()==1 ? 1 : x.LDim() );
    blas::Her
    ( uploChar, m, alpha, x.LockedBuffer(), incx, A.Buffer(), A.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Her2
( UpperOrLower uplo,
  T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("Her2");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    if( (x.Width() != 1 && x.Height() != 1) || 
        (y.Width() != 1 && y.Height() != 1) )
        throw std::logic_error("x and y must be vectors");
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    const int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
    if( xLength != A.Height() || yLength != A.Height() )
        throw std::logic_error("x and y must conform with A");
#endif
    const char uploChar = UpperOrLowerToChar( uplo );
    const int m = A.Height();
    const int incx = ( x.Width()==1 ? 1 : x.LDim() );
    const int incy = ( y.Width()==1 ? 1 : y.LDim() );
    blas::Her2
    ( uploChar, m,
      alpha, x.LockedBuffer(), incx, y.LockedBuffer(), incy,
             A.Buffer(), A.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Symv
( UpperOrLower uplo,
  T alpha, const Matrix<T>& A, const Matrix<T>& x, T beta, Matrix<T>& y )
{
#ifndef RELEASE
    PushCallStack("Symv");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    if( ( x.Height() != 1 && x.Width() != 1 ) ||
        ( y.Height() != 1 && y.Width() != 1 ) )
        throw std::logic_error("x and y must be vectors");
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    const int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
    if( A.Height() != xLength || A.Height() != yLength )
        throw std::logic_error("A must conform with x and y");
#endif
    const char uploChar = UpperOrLowerToChar( uplo );
    const int m = A.Height();
    const int incx = ( x.Width()==1 ? 1 : x.LDim() );
    const int incy = ( y.Width()==1 ? 1 : y.LDim() );
    blas::Symv
    ( uploChar, m, 
      alpha, A.LockedBuffer(), A.LDim(), x.LockedBuffer(), incx, 
      beta,  y.Buffer(), incy );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Syr( UpperOrLower uplo, T alpha, const Matrix<T>& x, Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("Syr");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    if( x.Width() != 1 && x.Height() != 1 )
        throw std::logic_error("x must be a vector");
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    if( xLength != A.Height() )
        throw std::logic_error("x must conform with A");
#endif
    const char uploChar = UpperOrLowerToChar( uplo );
    const int m = A.Height();
    const int incx = ( x.Width()==1 ? 1 : x.LDim() );
    blas::Syr
    ( uploChar, m, alpha, x.LockedBuffer(), incx, A.Buffer(), A.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Syr2
( UpperOrLower uplo, 
  T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("Syr2");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    if( (x.Width() != 1 && x.Height() != 1) || 
        (y.Width() != 1 && y.Height() != 1) )
        throw std::logic_error("x and y must be vectors");
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    const int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
    if( xLength != A.Height() || yLength != A.Height() )
        throw std::logic_error("x and y must conform with A");
#endif
    const char uploChar = UpperOrLowerToChar( uplo );
    const int m = A.Height();
    const int incx = ( x.Width()==1 ? 1 : x.LDim() );
    const int incy = ( y.Width()==1 ? 1 : y.LDim() );
    blas::Syr2
    ( uploChar, m,
      alpha, x.LockedBuffer(), incx, y.LockedBuffer(), incy,
             A.Buffer(), A.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Trmv
( UpperOrLower uplo, Orientation orientation, Diagonal diagonal,
  const Matrix<T>& A, Matrix<T>& x )
{
#ifndef RELEASE
    PushCallStack("Trmv");
    if( x.Height() != 1 && x.Width() != 1 )
        throw std::logic_error("x must be a vector");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    if( xLength != A.Height() )
        throw std::logic_error("x must conform with A");
#endif
    const char uploChar = UpperOrLowerToChar( uplo );
    const char transChar = OrientationToChar( orientation );
    const char diagChar = DiagonalToChar( diagonal );
    const int m = A.Height();
    const int incx = ( x.Width()==1 ? 1 : x.LDim() );
    blas::Trmv
    ( uploChar, transChar, diagChar, m,
      A.LockedBuffer(), A.LDim(), x.Buffer(), incx );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
Trsv
( UpperOrLower uplo, Orientation orientation, Diagonal diagonal,
  const Matrix<F>& A, Matrix<F>& x )
{
#ifndef RELEASE
    PushCallStack("Trsv");
    if( x.Height() != 1 && x.Width() != 1 )
        throw std::logic_error("x must be a vector");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    if( xLength != A.Height() )
        throw std::logic_error("x must conform with A");
#endif
    const char uploChar = UpperOrLowerToChar( uplo );
    const char transChar = OrientationToChar( orientation );
    const char diagChar = DiagonalToChar( diagonal );
    const int m = A.Height();
    const int incx = ( x.Width()==1 ? 1 : x.LDim() );
    blas::Trsv
    ( uploChar, transChar, diagChar, m,
      A.LockedBuffer(), A.LDim(), x.Buffer(), incx );
#ifndef RELEASE
    PopCallStack();
#endif
}

//----------------------------------------------------------------------------//
// Local BLAS-like routines: Level 3                                          //
//----------------------------------------------------------------------------//

template<typename T>
inline void
Gemm
( Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const Matrix<T>& A, const Matrix<T>& B, T beta, Matrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("Gemm");
    if( orientationOfA == NORMAL && orientationOfB == NORMAL )
    {
        if( A.Height() != C.Height() ||
            B.Width()  != C.Width()  ||
            A.Width()  != B.Height() )
            throw std::logic_error("Nonconformal GemmNN");
    }
    else if( orientationOfA == NORMAL )
    {
        if( A.Height() != C.Height() ||
            B.Height() != C.Width()  ||
            A.Width()  != B.Width() )
            throw std::logic_error("Nonconformal GemmN(T/C)");
    }
    else if( orientationOfB == NORMAL )
    {
        if( A.Width()  != C.Height() ||
            B.Width()  != C.Width()  ||
            A.Height() != B.Height() )
            throw std::logic_error("Nonconformal Gemm(T/C)N");
    }
    else
    {
        if( A.Width()  != C.Height() ||
            B.Height() != C.Width()  ||
            A.Height() != B.Width() )
            throw std::logic_error("Nonconformal Gemm(T/C)(T/C)");
    }
#endif
    const char transA = OrientationToChar( orientationOfA );
    const char transB = OrientationToChar( orientationOfB );
    const int m = C.Height();
    const int n = C.Width();
    const int k = ( orientationOfA == NORMAL ? A.Width() : A.Height() );
    if( k != 0 )
    {
        blas::Gemm
        ( transA, transB, m, n, k, 
          alpha, A.LockedBuffer(), A.LDim(), B.LockedBuffer(), B.LDim(),
          beta,  C.Buffer(),       C.LDim() );
    }
    else
    {
        Scal( beta, C );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Hemm
( Side side, UpperOrLower uplo,
  T alpha, const Matrix<T>& A, const Matrix<T>& B, T beta, Matrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("Hemm");
#endif
    const char sideChar = SideToChar( side );
    const char uploChar = UpperOrLowerToChar( uplo );
    blas::Hemm
    ( sideChar, uploChar, C.Height(), C.Width(),
      alpha, A.LockedBuffer(), A.LDim(), 
             B.LockedBuffer(), B.LDim(),
      beta,  C.Buffer(),       C.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Her2k
( UpperOrLower uplo, Orientation orientation,
  T alpha, const Matrix<T>& A, const Matrix<T>& B, T beta, Matrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("Her2k");
    if( orientation == NORMAL )
    {
        if( A.Height() != C.Height() || A.Height() != C.Width() ||
            B.Height() != C.Height() ||B.Height() != C.Width() )
            throw std::logic_error("Nonconformal Her2k");
    }
    else if( orientation == ADJOINT )
    {
        if( A.Width() != C.Height() || A.Width() != C.Width() ||
            B.Width() != C.Height() || B.Width() != C.Width() )
            throw std::logic_error("Nonconformal Her2k");
    }
    else
        throw std::logic_error
        ("Her2k only accepts NORMAL and ADJOINT options");
#endif
    const char uploChar = UpperOrLowerToChar( uplo );
    const char transChar = OrientationToChar( orientation );
    const int k = ( orientation == NORMAL ? A.Width() : A.Height() );
    blas::Her2k
    ( uploChar, transChar, C.Height(), k, 
      alpha, A.LockedBuffer(), A.LDim(), 
             B.LockedBuffer(), B.LDim(),
      beta,  C.Buffer(),       C.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Herk
( UpperOrLower uplo, Orientation orientation,
  T alpha, const Matrix<T>& A, T beta, Matrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("Herk");
    if( orientation == NORMAL )
    {
        if( A.Height() != C.Height() || A.Height() != C.Width() )
            throw std::logic_error("Nonconformal Herk");
    }
    else if( orientation == ADJOINT )
    {
        if( A.Width() != C.Height() || A.Width() != C.Width() )
            throw std::logic_error("Nonconformal Herk");
    }
    else
        throw std::logic_error("Herk only accepts NORMAL and ADJOINT options.");
#endif
    const char uploChar = UpperOrLowerToChar( uplo );
    const char transChar = OrientationToChar( orientation );
    const int k = ( orientation == NORMAL ? A.Width() : A.Height() );
    blas::Herk
    ( uploChar, transChar, C.Height(), k, 
      alpha, A.LockedBuffer(), A.LDim(), 
      beta,  C.Buffer(),       C.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Hetrmm( UpperOrLower uplo, Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("Hetrmm");
#endif
    const char uploChar = UpperOrLowerToChar( uplo );
    blas::Hetrmm( uploChar, A.Height(), A.Buffer(), A.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Symm
( Side side, UpperOrLower uplo,
  T alpha, const Matrix<T>& A, const Matrix<T>& B, T beta, Matrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("Symm");
#endif
    const char sideChar = SideToChar( side );
    const char uploChar = UpperOrLowerToChar( uplo );
    blas::Symm
    ( sideChar, uploChar, C.Height(), C.Width(),
      alpha, A.LockedBuffer(), A.LDim(), 
             B.LockedBuffer(), B.LDim(),
      beta,  C.Buffer(),       C.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Syr2k
( UpperOrLower uplo, Orientation orientation,
  T alpha, const Matrix<T>& A, const Matrix<T>& B, T beta, Matrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("Syr2k");
    if( orientation == NORMAL )
    {
        if( A.Height() != C.Height() || A.Height() != C.Width() ||
            B.Height() != C.Height() ||B.Height() != C.Width()    )
            throw std::logic_error("Nonconformal Syr2k");
    }
    else if( orientation == TRANSPOSE )
    {
        if( A.Width() != C.Height() || A.Width() != C.Width() ||
            B.Width() != C.Height() || B.Width() != C.Width()   )
            throw std::logic_error("Nonconformal Syr2k");
    }
    else
        throw std::logic_error
        ("Syr2k only accepts NORMAL and TRANSPOSE options");
#endif
    const char uploChar = UpperOrLowerToChar( uplo );
    const char transChar = OrientationToChar( orientation );
    const int k = ( orientation == NORMAL ? A.Width() : A.Height() );
    blas::Syr2k
    ( uploChar, transChar, C.Height(), k, 
      alpha, A.LockedBuffer(), A.LDim(), 
             B.LockedBuffer(), B.LDim(),
      beta,  C.Buffer(),       C.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Syrk
( UpperOrLower uplo, Orientation orientation,
  T alpha, const Matrix<T>& A, T beta, Matrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("Syrk");
    if( orientation == NORMAL )
    {
        if( A.Height() != C.Height() || A.Height() != C.Width() )
            throw std::logic_error("Nonconformal Syrk");
    }
    else if( orientation == TRANSPOSE )
    {
        if( A.Width() != C.Height() || A.Width() != C.Width() )
            throw std::logic_error("Nonconformal Syrk");
    }
    else
        throw std::logic_error
        ("Syrk only accepts NORMAL and TRANSPOSE options");
#endif
    const char uploChar = UpperOrLowerToChar( uplo );
    const char transChar = OrientationToChar( orientation );
    const int k = ( orientation == NORMAL ? A.Width() : A.Height() );
    blas::Syrk
    ( uploChar, transChar, C.Height(), k, 
      alpha, A.LockedBuffer(), A.LDim(), 
      beta,  C.Buffer(),       C.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Trmm
( Side side, UpperOrLower uplo, 
  Orientation orientation, Diagonal diagonal,
  T alpha, const Matrix<T>& A, Matrix<T>& B )
{
#ifndef RELEASE
    PushCallStack("Trmm");
    if( A.Height() != A.Width() )
        throw std::logic_error("Triangular matrix must be square");
    if( side == LEFT )
    {
        if( A.Height() != B.Height() )
            throw std::logic_error("Nonconformal Trmm");
    }
    else
    {
        if( A.Height() != B.Width() )
            throw std::logic_error("Nonconformal Trmm");
    }
#endif
    const char sideChar = SideToChar( side );
    const char uploChar = UpperOrLowerToChar( uplo );
    const char transChar = OrientationToChar( orientation );
    const char diagChar = DiagonalToChar( diagonal );
    blas::Trmm
    ( sideChar, uploChar, transChar, diagChar, B.Height(), B.Width(),
      alpha, A.LockedBuffer(), A.LDim(), B.Buffer(), B.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
Trsm
( Side side, UpperOrLower uplo,
  Orientation orientation,Diagonal diagonal,
  F alpha, const Matrix<F>& A, Matrix<F>& B,
  bool checkIfSingular )
{
#ifndef RELEASE
    PushCallStack("Trsm");
    if( A.Height() != A.Width() )
        throw std::logic_error("Triangular matrix must be square");
    if( side == LEFT )
    {
        if( A.Height() != B.Height() )
            throw std::logic_error("Nonconformal Trsm");
    }
    else
    {
        if( A.Height() != B.Width() )
            throw std::logic_error("Nonconformal Trsm");
    }
#endif
    const char sideChar = SideToChar( side );
    const char uploChar = UpperOrLowerToChar( uplo );
    const char transChar = OrientationToChar( orientation );
    const char diagChar = DiagonalToChar( diagonal );
    if( checkIfSingular && diagonal != UNIT )
    {
        const int n = A.Height();
        for( int j=0; j<n; ++j )
            if( A.Get(j,j) == (F)0 )
                throw SingularMatrixException();
    }
    blas::Trsm
    ( sideChar, uploChar, transChar, diagChar, B.Height(), B.Width(),
      alpha, A.LockedBuffer(), A.LDim(), B.Buffer(), B.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

//----------------------------------------------------------------------------//
// Distributed BLAS-like routines: Level 1                                    //
//----------------------------------------------------------------------------//

template<typename T,Distribution U,Distribution V>
inline void
Axpy( T alpha, const DistMatrix<T,U,V>& X, DistMatrix<T,U,V>& Y )
{
#ifndef RELEASE
    PushCallStack("Axpy");
    if( X.Grid() != Y.Grid() )
        throw std::logic_error
        ("X and Y must be distributed over the same grid");
    if( X.ColAlignment() != Y.ColAlignment() ||
        X.RowAlignment() != Y.RowAlignment() )
        throw std::logic_error("Axpy requires X and Y be aligned");
#endif
    Axpy( alpha, X.LockedLocalMatrix(), Y.LocalMatrix() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V,
                    Distribution W,Distribution Z>
inline void
Copy( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B )
{
#ifndef RELEASE
    PushCallStack("Copy");
#endif
    B = A;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V,
                    Distribution W,Distribution Z>
inline void
DiagonalScale
( Side side, Orientation orientation, 
  const DistMatrix<T,U,V>& d, DistMatrix<T,W,Z>& X )
{
#ifndef RELEASE
    PushCallStack("DiagonalScale");
#endif
    if( side == LEFT )
    {
        if( U == W && V == STAR && d.ColAlignment() == X.ColAlignment() )
        {
            DiagonalScale
            ( LEFT, orientation, d.LockedLocalMatrix(), X.LocalMatrix() );
        }
        else
        {
            DistMatrix<T,W,STAR> d_W_STAR( X.Grid() );
            d_W_STAR = d;
            DiagonalScale
            ( LEFT, orientation, 
              d_W_STAR.LockedLocalMatrix(), X.LocalMatrix() );
        }
    }
    else
    {
        if( U == Z && V == STAR && d.ColAlignment() == X.RowAlignment() )
        {
            DiagonalScale
            ( RIGHT, orientation, d.LockedLocalMatrix(), X.LocalMatrix() );
        }
        else
        {
            DistMatrix<T,Z,STAR> d_Z_STAR( X.Grid() );
            d_Z_STAR = d;
            DiagonalScale
            ( RIGHT, orientation, 
              d_Z_STAR.LockedLocalMatrix(), X.LocalMatrix() );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F,Distribution U,Distribution V,
                    Distribution W,Distribution Z>
inline void
DiagonalSolve
( Side side, Orientation orientation, 
  const DistMatrix<F,U,V>& d, DistMatrix<F,W,Z>& X,
  bool checkIfSingular )
{
#ifndef RELEASE
    PushCallStack("DiagonalSolve");
#endif
    if( side == LEFT )
    {
        if( U == W && V == STAR && d.ColAlignment() == X.ColAlignment() )
        {
            DiagonalSolve
            ( LEFT, orientation, d.LockedLocalMatrix(), X.LocalMatrix(),
              checkIfSingular );
        }
        else
        {
            DistMatrix<F,W,STAR> d_W_STAR( X.Grid() );
            d_W_STAR = d;
            DiagonalSolve
            ( LEFT, orientation, 
              d_W_STAR.LockedLocalMatrix(), X.LocalMatrix(), checkIfSingular );
        }
    }
    else
    {
        if( U == Z && V == STAR && d.ColAlignment() == X.RowAlignment() )
        {
            DiagonalSolve
            ( RIGHT, orientation, d.LockedLocalMatrix(), X.LocalMatrix(),
              checkIfSingular );
        }
        else
        {
            DistMatrix<F,Z,STAR> d_Z_STAR( X.Grid() );
            d_Z_STAR = d;
            DiagonalSolve
            ( RIGHT, orientation, 
              d_Z_STAR.LockedLocalMatrix(), X.LocalMatrix(), checkIfSingular );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

// Our extended Dotc is equivalent to our extended Dot, 
// but we are burdened with consistency
template<typename T,Distribution U,Distribution V,
                    Distribution W,Distribution Z>
inline T
Dotc( const DistMatrix<T,U,V>& x, const DistMatrix<T,W,Z>& y )
{
#ifndef RELEASE
    PushCallStack("Dotc");
#endif
    T globalDot = Dot( x, y );
#ifndef RELEASE
    PopCallStack();
#endif
    return globalDot;
}

template<typename T,Distribution U,Distribution V>
inline void
Scal( T alpha, DistMatrix<T,U,V>& A )
{
#ifndef RELEASE
    PushCallStack("Scal");
#endif
    Scal( alpha, A.LocalMatrix() );
#ifndef RELEASE
    PopCallStack();
#endif
}

//----------------------------------------------------------------------------//
// Distributed BLAS-like routines: Level 1 (extensions)                       //
//----------------------------------------------------------------------------//

template<typename T,Distribution U,Distribution V>
inline void
Conjugate( DistMatrix<T,U,V>& A )
{
#ifndef RELEASE
    PushCallStack("Conjugate (in-place)");
#endif
    Conjugate( A.LocalMatrix() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V,
                    Distribution W,Distribution Z>
inline void
Conjugate( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B )
{
#ifndef RELEASE
    PushCallStack("Conjugate");
#endif
    B = A;
    Conjugate( B ); 
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V,
                    Distribution W,Distribution Z>
inline void
Adjoint( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B )
{
#ifndef RELEASE
    PushCallStack("Adjoint");
#endif
    if( U == Z && V == W && 
        A.ColAlignment() == B.RowAlignment() && 
        A.RowAlignment() == B.ColAlignment() )
    {
        Adjoint( A.LockedLocalMatrix(), B.LocalMatrix() );
    }
    else
    {
        DistMatrix<T,Z,W> C( B.Grid() );
        if( B.ConstrainedColAlignment() )
            C.AlignRowsWith( B );
        if( B.ConstrainedRowAlignment() )
            C.AlignColsWith( B );
        C = A;

        if( !B.Viewing() )
        {
            if( !B.ConstrainedColAlignment() )
                B.AlignColsWith( C );
            if( !B.ConstrainedRowAlignment() )
                B.AlignRowsWith( C );
            B.ResizeTo( A.Width(), A.Height() );
        }
        else if( B.Height() != A.Width() || B.Width() != A.Height() )
        {
            throw std::logic_error
            ("If Adjoint'ing into a view, it must be the right size");
        }
        Adjoint( C.LockedLocalMatrix(), B.LocalMatrix() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V,
                    Distribution W,Distribution Z>
inline void
Transpose( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B )
{
#ifndef RELEASE
    PushCallStack("Transpose");
#endif
    if( U == Z && V == W && 
        A.ColAlignment() == B.RowAlignment() && 
        A.RowAlignment() == B.ColAlignment() )
    {
        Transpose( A.LockedLocalMatrix(), B.LocalMatrix() );
    }
    else
    {
        DistMatrix<T,Z,W> C( B.Grid() );
        if( B.ConstrainedColAlignment() )
            C.AlignRowsWith( B );
        if( B.ConstrainedRowAlignment() )
            C.AlignColsWith( B );
        C = A;

        if( !B.Viewing() )
        {
            if( !B.ConstrainedColAlignment() )
                B.AlignColsWith( C );
            if( !B.ConstrainedRowAlignment() )
                B.AlignRowsWith( C );
            B.ResizeTo( A.Width(), A.Height() );
        }
        else if( B.Height() != A.Width() || B.Width() != A.Height() )
        {
            throw std::logic_error
            ("If Transposing into a view, it must be the right size");
        }
        Transpose( C.LockedLocalMatrix(), B.LocalMatrix() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

//----------------------------------------------------------------------------//
// Distributed BLAS-like routines: Level 2                                    //
//----------------------------------------------------------------------------//

// Our extended Ger and Gerc are equivalent
template<typename T>
inline void
Gerc
( T alpha, const DistMatrix<T,MC,MR>& x,
           const DistMatrix<T,MC,MR>& y,
                 DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("Gerc");
#endif
    Ger( alpha, x, y, A );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem

#endif /* ELEMENTAL_BASIC_HPP */

