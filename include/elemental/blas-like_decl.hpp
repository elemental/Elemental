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
void Axpy( T alpha, const Matrix<T>& X, Matrix<T>& Y );
template<typename T>
void Axpy( typename Base<T>::type alpha, const Matrix<T>& X, Matrix<T>& Y );

// Parallel version
template<typename T, Distribution U, Distribution V>
void Axpy( T alpha, const DistMatrix<T,U,V>& X, DistMatrix<T,U,V>& Y );
template<typename T, Distribution U, Distribution V>
void Axpy
( typename Base<T>::type alpha, 
  const DistMatrix<T,U,V>& X, DistMatrix<T,U,V>& Y );

//
// AxpyTriangle (Alpha X Plus Y):
//
// Y := alpha X + Y,
//
// where X and Y are triangular.
//

// Serial version
template<typename T>
void AxpyTriangle
( UpperOrLower uplo, T alpha, const Matrix<T>& X, Matrix<T>& Y );
template<typename T>
void AxpyTriangle
( UpperOrLower uplo, typename Base<T>::type alpha, 
  const Matrix<T>& X, Matrix<T>& Y );

// Parallel version
template<typename T, Distribution U, Distribution V>
void AxpyTriangle
( UpperOrLower uplo, T alpha, 
  const DistMatrix<T,U,V>& X, DistMatrix<T,U,V>& Y );
template<typename T, Distribution U, Distribution V>
void AxpyTriangle
( UpperOrLower uplo, typename Base<T>::type alpha, 
  const DistMatrix<T,U,V>& X, DistMatrix<T,U,V>& Y );

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
( LeftOrRight side, Orientation orientation, 
  const Matrix<T>& d, Matrix<T>& X );
template<typename T>
void DiagonalScale
( LeftOrRight side, Orientation orientation,
  const Matrix<typename Base<T>::type>& d, Matrix<T>& X );

// Parallel version
template<typename T,
         Distribution U,Distribution V,
         Distribution W,Distribution Z>
void DiagonalScale
( LeftOrRight side, Orientation orientation, 
  const DistMatrix<T,U,V>& d, DistMatrix<T,W,Z>& X );
template<typename T,
         Distribution U,Distribution V,
         Distribution W,Distribution Z>
void DiagonalScale
( LeftOrRight side, Orientation orientation, 
  const DistMatrix<typename Base<T>::type,U,V>& d, DistMatrix<T,W,Z>& X );

//
// DiagonalSolve
//
// Performs either X := op(D)^-1 X or 
//                 X := X op(D)^-1

// Serial version
template<typename F>
void DiagonalSolve
( LeftOrRight side, Orientation orientation, 
  const Matrix<F>& d, Matrix<F>& X, bool checkIfSingular=false );
template<typename F>
void DiagonalSolve
( LeftOrRight side, Orientation orientation, 
  const Matrix<typename Base<F>::type>& d, Matrix<F>& X, 
  bool checkIfSingular=false );

// Parallel version
template<typename F,
         Distribution U,Distribution V,
         Distribution W,Distribution Z>
void DiagonalSolve
( LeftOrRight side, Orientation orientation, 
  const DistMatrix<F,U,V>& d, DistMatrix<F,W,Z>& X, 
  bool checkIfSingular=false );
template<typename F,
         Distribution U,Distribution V,
         Distribution W,Distribution Z>
void DiagonalSolve
( LeftOrRight side, Orientation orientation, 
  const DistMatrix<typename Base<F>::type,U,V>& d, DistMatrix<F,W,Z>& X, 
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

template<typename F>
typename Base<F>::type Nrm2( const Matrix<F>& x );
template<typename F>
typename Base<F>::type Nrm2( const DistMatrix<F>& x );

// 
// Scale (also Scal for legacy reasons):
//
// X := alpha X
//

// Serial version
template<typename T>
void Scal( T alpha, Matrix<T>& X );
template<typename T>
void Scale( T alpha, Matrix<T>& X );
template<typename T>
void Scal( typename Base<T>::type alpha, Matrix<T>& X );
template<typename T>
void Scale( typename Base<T>::type alpha, Matrix<T>& X );
   
// Parallel version
template<typename T,Distribution U,Distribution V>
void Scal( T alpha, DistMatrix<T,U,V>& X );
template<typename T,Distribution U,Distribution V>
void Scale( T alpha, DistMatrix<T,U,V>& X );
template<typename T,Distribution U,Distribution V>
void Scal( typename Base<T>::type alpha, DistMatrix<T,U,V>& X );
template<typename T,Distribution U,Distribution V>
void Scale( typename Base<T>::type alpha, DistMatrix<T,U,V>& X );
    
//----------------------------------------------------------------------------//
// Level 1 BLAS-like extensions                                               //
//----------------------------------------------------------------------------//

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
         Distribution U,Distribution V,
         Distribution W,Distribution Z>
void Adjoint( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B );


//
// Conjugate: 
//
// Conjugates a matrix. The in-place version performs A := Conjugate(A), 
// while the out-of-place sets B := Conjugate(A).
//

// In-place versions
template<typename Z>
void Conjugate( Matrix<Z>& A );
template<typename Z>
void Conjugate( Matrix<Complex<Z> >& A );
template<typename T,Distribution U,Distribution V>
void Conjugate( DistMatrix<T,U,V>& A );

// Out-of-place versions
template<typename T>
void Conjugate( const Matrix<T>& A, Matrix<T>& B );
template<typename T, 
         Distribution U,Distribution V,
         Distribution W,Distribution Z>
void Conjugate( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B );

//
// MakeTriangular:
//
// Force the matrix to be either lower or upper triangular
//
template<typename T>
void MakeTriangular( UpperOrLower uplo, Matrix<T>& A );
template<typename T,Distribution U,Distribution V>
void MakeTriangular( UpperOrLower uplo, DistMatrix<T,U,V>& A );

// MakeUnitTriangular?

//
// MakeTrapezoidal:
//
// Force the matrix to be lower or upper trapezoidal, with the diagonal
// defined relative to either the top-left or bottom-right corner of the matrix
// (based on the 'side' parameter). The 'offset' parameter determines where the
// last nonzero diagonal is, with '0' meaning the main diagonal, '1' meaning 
// the superdiagonal, '-1' meaning the subdiagonal, and likewise for all other
// integer values.
//
template<typename T>
void MakeTrapezoidal
( LeftOrRight side, UpperOrLower uplo, int offset, Matrix<T>& A );
template<typename T,Distribution U,Distribution V>
void MakeTrapezoidal
( LeftOrRight side, UpperOrLower uplo, int offset, DistMatrix<T,U,V>& A );

// TODO: MakeTriangular

//
// MakeHermitian:
//
// Turn an implicitly Hermitian matrix into an explicitly Hermitian one by 
// forcing the diagonal to be real and conjugate-transposing the specified 
// strictly triangular portion of the matrix into the other.
//

template<typename T>
void MakeHermitian( UpperOrLower uplo, Matrix<T>& A );
template<typename T>
void MakeHermitian( UpperOrLower uplo, DistMatrix<T>& A );

//
// MakeReal:
//
// Modify a matrix to ensure that all of its imaginary components are zero.
//

template<typename T>
void MakeReal( Matrix<T>& A );
template<typename T,Distribution U,Distribution V>
void MakeReal( DistMatrix<T,U,V>& A );

//
// MakeSymmetric:
//
// Turn an implicitly symmetric matrix into an explicitly symmetric one by 
// forcing the diagonal to be real and conjugate-transposing the specified 
// strictly triangular portion of the matrix into the other.
//

template<typename T>
void MakeSymmetric( UpperOrLower uplo, Matrix<T>& A );
template<typename T>
void MakeSymmetric( UpperOrLower uplo, DistMatrix<T>& A );

//
// ScaleTrapezoid:
//
// Scale only a trapezoidal portion of a matrix, using the parameter convention
// described above.
//
template<typename T>
void ScaleTrapezoid
( T alpha, LeftOrRight side, UpperOrLower uplo, int offset, 
  Matrix<T>& A );
template<typename T,Distribution U,Distribution V>
void ScaleTrapezoid
( T alpha, LeftOrRight side, UpperOrLower uplo, int offset, 
  DistMatrix<T,U,V>& A );

// TODO: ScaleTriangle

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
         Distribution U,Distribution V,
         Distribution W,Distribution Z>
void Transpose( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B );

//
// Zero the contents of a matrix
//
template<typename T>
void Zero( Matrix<T>& A );
template<typename T,Distribution U,Distribution V>
void Zero( DistMatrix<T,U,V>& A );

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
  T alpha, const DistMatrix<T>& A, const DistMatrix<T>& x,
  T beta,                                DistMatrix<T>& y );
template<typename T>
void Gemv
( Orientation orientationOfA,
  T alpha, const DistMatrix<T>& A, const DistMatrix<T,VC,STAR>& x,
  T beta,                                DistMatrix<T,VC,STAR>& y );

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
( T alpha, const DistMatrix<T>& x, const DistMatrix<T>& y,
                 DistMatrix<T>& A );

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
( T alpha, const DistMatrix<T>& x, const DistMatrix<T>& y,
                 DistMatrix<T>& A );

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
( T alpha, const DistMatrix<T>& x, const DistMatrix<T>& y,
                 DistMatrix<T>& A );

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
  T alpha, const DistMatrix<T>& A, const DistMatrix<T>& x,
  T beta,        DistMatrix<T>& y );

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
( UpperOrLower uplo, T alpha, const DistMatrix<T>& x, DistMatrix<T>& A );

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
  T alpha, const DistMatrix<T>& x, const DistMatrix<T>& y,
                 DistMatrix<T>& A );

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
  T alpha, const DistMatrix<T>& A, const DistMatrix<T>& x,
  T beta,        DistMatrix<T>& y );

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
( UpperOrLower uplo, T alpha, const DistMatrix<T>& x, DistMatrix<T>& A );

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
  T alpha, const DistMatrix<T>& x, const DistMatrix<T>& y,
                 DistMatrix<T>& A );

//
// Trmv (TRiangular Matrix-Vector multiply):
//
// Performs the update
//   x := orientation( A ) x,
// where 'uplo' determines whether or not A is to be implicitly treated as 
// lower or upper triangular, and 'diag' specifies whether it has an 
// implicit unit diagonal.
//

// Serial version
template<typename T>
void Trmv
( UpperOrLower uplo, Orientation orientation, UnitOrNonUnit diag,
  const Matrix<T>& A, Matrix<T>& x );

// Parallel version
template<typename T>
void Trmv
( UpperOrLower uplo, Orientation orientation, UnitOrNonUnit diag,
  const DistMatrix<T>& A, DistMatrix<T>& x );

//
// Trsv (TRiangular Solve with a Vector):
//
// Performs the update
//   x := orientation( A )^-1 x,
// where 'uplo' determines whether or not A is to be implicitly treated as 
// lower or upper triangular, and 'diag' specifies whether it has an 
// implicit unit diagonal.
//

// Serial version
template<typename F>
void Trsv
( UpperOrLower uplo, Orientation orientation, UnitOrNonUnit diag,
  const Matrix<F>& A, Matrix<F>& x );

// Parallel version
template<typename F>
void Trsv
( UpperOrLower uplo, Orientation orientation, UnitOrNonUnit diag,
  const DistMatrix<F>& A, DistMatrix<F>& x );

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
  T alpha, const DistMatrix<T>& A, const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C );

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
( LeftOrRight side, UpperOrLower uplo,
  T alpha, const Matrix<T>& A, const Matrix<T>& B, T beta, Matrix<T>& C );

// Parallel version
template<typename T>
void Hemm
( LeftOrRight side, UpperOrLower uplo,
  T alpha, const DistMatrix<T>& A, const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C );

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
  T alpha, const DistMatrix<T>& A, const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C );

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
  T alpha, const DistMatrix<T>& A, T beta, DistMatrix<T>& C );

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
( LeftOrRight side, UpperOrLower uplo,
  T alpha, const Matrix<T>& A, const Matrix<T>& B, T beta, Matrix<T>& C ); 

// Parallel version
template<typename T>
void Symm
( LeftOrRight side, UpperOrLower uplo,
  T alpha, const DistMatrix<T>& A, const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C );

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
  T alpha, const DistMatrix<T>& A, const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C );

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
  T alpha, const DistMatrix<T>& A, T beta, DistMatrix<T>& C );

//
// Trmm (TRiangular Matrix-Matrix multiplication):
//
// Performs the update
//   B := alpha orientation( A ) B,  { side = LEFT }
// or
//   B := alpha B orientation( A ),  { side = RIGHT }
// where 'uplo' determines whether A is assumed to be upper or lower 
// triangular and 'diag' determines whether A has an implicit unit
// diagonal.
//

// Serial version
template<typename T>
void Trmm
( LeftOrRight side, UpperOrLower uplo, 
  Orientation orientation, UnitOrNonUnit diag,
  T alpha, const Matrix<T>& A, Matrix<T>& B );

// Parallel version
template<typename T>
void Trmm
( LeftOrRight side, UpperOrLower uplo, 
  Orientation orientation, UnitOrNonUnit diag,
  T alpha, const DistMatrix<T>& A, DistMatrix<T>& B );

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
  T alpha, const DistMatrix<T>& A, const DistMatrix<T>& B, 
           const DistMatrix<T>& C, const DistMatrix<T>& D,
  T beta,        DistMatrix<T>& E );

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
  T alpha, const DistMatrix<T>& A, const DistMatrix<T>& B, 
  T beta,        DistMatrix<T>& C );

//
// Trsm (TRiangular Solve with Multiple right-hand sides):
//
// Performs the update
//   B := alpha orientation( A )^-1 B,  { side = LEFT }
// or
//   B := alpha B orientation( A )^-1,  { side = RIGHT }
// where 'uplo' determines whether A is assumed to be upper or lower
// triangular and 'diag' determines whether A has an implicit unit
// diagonal.
//

// Serial version
template<typename F>
void Trsm
( LeftOrRight side, UpperOrLower uplo, 
  Orientation orientation, UnitOrNonUnit diag,
  F alpha, const Matrix<F>& A, Matrix<F>& B, 
  bool checkIfSingular=false ); 
        
// Parallel version
template<typename F>
void Trsm
( LeftOrRight side, UpperOrLower uplo, 
  Orientation orientation, UnitOrNonUnit diag,
  F alpha, const DistMatrix<F>& A, DistMatrix<F>& B,
  bool checkIfSingular=false );

//
// Trtrsm (TRiangular TRiangular Solve with Multiple right-hand sides):
//
// Performs the update
//   B := alpha orientation( A )^-1 B,  { side = LEFT }
// or
//   B := alpha B orientation( A )^-1,  { side = RIGHT }
// where 'uplo' determines whether A and B are assumed to be upper or lower
// triangular and 'diag' determines whether A has an implicit unit
// diagonal.
//

// Serial version
template<typename F>
void Trtrsm
( LeftOrRight side, UpperOrLower uplo, 
  Orientation orientation, UnitOrNonUnit diag,
  F alpha, const Matrix<F>& A, Matrix<F>& B, 
  bool checkIfSingular=false ); 
        
// Parallel version
template<typename F>
void Trtrsm
( LeftOrRight side, UpperOrLower uplo, 
  Orientation orientation, UnitOrNonUnit diag,
  F alpha, const DistMatrix<F>& A, DistMatrix<F>& B,
  bool checkIfSingular=false );

// Trtrmm (TRiangular TRiangular Matrix-Matrix multiply):
//
// Either L := tril(L^[T/H] L) or U := triu(U U^[T/H])
//
// NOTE: This is not a standard BLAS routine and is similar to the LAPACK
//       routine ?LAUUM. 
//

// Serial version
template<typename T>
void Trtrmm( Orientation orientation, UpperOrLower uplo, Matrix<T>& A );

// Parallel version
template<typename T>
void Trtrmm( Orientation orientation, UpperOrLower uplo, DistMatrix<T>& A );

// Trdtrmm (TRiangular Diagonal TRiangular Matrix-Matrix multiply):
//
// Either L := tril(L^T inv(D) L) or U := triu(U inv(D) U^T)
//
// NOTE: This is not a standard BLAS routine and is similar to the LAPACK
//       routine ?LAUUM. 
//

// Serial version
template<typename F>
void Trdtrmm( Orientation orientation, UpperOrLower uplo, Matrix<F>& A );

// Parallel version
template<typename F>
void Trdtrmm( Orientation orientation, UpperOrLower uplo, DistMatrix<F>& A );

// TwoSidedTrmm
//
// Either A := L^H A L or A := U A U^H
//
// NOTE: This is not a standard BLAS routine and is similar to the LAPACK 
//       routines ?sygst and ?hegst.

// Serial version
template<typename F>
void TwoSidedTrmm
( UpperOrLower uplo, UnitOrNonUnit diag, Matrix<F>& A, const Matrix<F>& B );

// Parallel version
template<typename F>
void TwoSidedTrmm
( UpperOrLower uplo, UnitOrNonUnit diag, 
  DistMatrix<F>& A, const DistMatrix<F>& B );

// TwoSidedTrsm
//
// Either A := inv(L) A inv(L)^H or A := inv(U)^H A inv(U)
//
// NOTE: This is not a standard BLAS routine and is similar to the LAPACK 
//       routines ?sygst and ?hegst.

// Serial version
template<typename F>
void TwoSidedTrsm
( UpperOrLower uplo, UnitOrNonUnit diag, Matrix<F>& A, const Matrix<F>& B );

// Parallel version
template<typename F>
void TwoSidedTrsm
( UpperOrLower uplo, UnitOrNonUnit diag, 
  DistMatrix<F>& A, const DistMatrix<F>& B );

} // namespace elem
