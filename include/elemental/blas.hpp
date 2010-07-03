/*
   This file is part of Elemental, a library for distributed-memory dense 
   linear algebra.

   Copyright (c) 2009-2010 Jack Poulson <jack.poulson@gmail.com>.
   All rights reserved.

   This file is released under the terms of the license contained in the file
   LICENSE-PURE.
*/
#ifndef ELEMENTAL_BLAS_HPP
#define ELEMENTAL_BLAS_HPP 1

#include "elemental/partitioning.hpp"

namespace elemental {
namespace blas {

//----------------------------------------------------------------------------//
// Level 1 BLAS                                                               //
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
R
Nrm2( const Matrix<R>& x ); 

#ifndef WITHOUT_COMPLEX
// Serial version for complex datatypes
template<typename R>
R
Nrm2( const Matrix< std::complex<R> >& x );
#endif

// Parallel version for real datatypes
template<typename R>
R
Nrm2( const DistMatrix<R,MC,MR>& x );

#ifndef WITHOUT_COMPLEX
// Parallel version for complex datatypes
template<typename R>
R
Nrm2( const DistMatrix< std::complex<R>, MC, MR >& x );
#endif

// 
// Scal:
//
// X := alpha X
//

// Serial version
template<typename T>
void
Scal( T alpha, Matrix<T>& X );
    
// Parallel version
template<typename T, Distribution U, Distribution V>
void
Scal
( T alpha, DistMatrix<T,U,V>& A );
    
//----------------------------------------------------------------------------//
// Level 1 BLAS extensions                                                    //
//----------------------------------------------------------------------------//

//
// Conj: 
//
// Conjugates a matrix. The in-place version performs A := Conj(A), while the 
// out-of-place sets B := Conj(A).
//

// In-place serial version for real datatypes. 
// Note: this is a no-op.
template<typename R>
void
Conj( Matrix<R>& A );

#ifndef WITHOUT_COMPLEX
// In-place serial version for complex datatypes.
template<typename R>
void
Conj( Matrix< std::complex<R> >& A );
#endif

// In-place parallel version
template<typename T, Distribution U, Distribution V>
void
Conj( DistMatrix<T,U,V>& A );

// Out-of-place serial version.
template<typename T>
void
Conj( const Matrix<T>& A, Matrix<T>& B );


// Out-of-place parallel version.
template<typename T, 
         Distribution U, Distribution V,
         Distribution W, Distribution Z>
void
Conj( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B );

//
// ConjTrans:
//
// B := A^H
//

// Serial version
template<typename T>
void
ConjTrans( const Matrix<T>& A, Matrix<T>& B );

// Parallel version
template<typename T, 
         Distribution U, Distribution V,
         Distribution W, Distribution Z>
void
ConjTrans( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B );

//
// Trans:
//
// B := A^T
//

// Serial version
template<typename T>
void
Trans( const Matrix<T>& A, Matrix<T>& B );

// Parallel version
template<typename T, 
         Distribution U, Distribution V,
         Distribution W, Distribution Z>
void
Trans( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B );

//----------------------------------------------------------------------------//
// Level 2 BLAS                                                               //
//----------------------------------------------------------------------------//

//
// Gemv (GEneral Matrix-Vector multiply):
//
// y := alpha orientation( A ) x + beta y,
// where orientation( A ) is determined by 'orientation'.
//

// Serial version
template<typename T>
void
Gemv
( Orientation orientation, 
  T alpha, const Matrix<T>& A, const Matrix<T>& x, T beta, Matrix<T>& y );

// Parallel version
template<typename T>
void
Gemv
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
void
Ger( T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A );

// Parallel version
template<typename T>
void
Ger
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
void
Gerc( T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A );

// Parallel version
template<typename T>
void
Gerc
( T alpha, const DistMatrix<T,MC,MR>& x, const DistMatrix<T,MC,MR>& y,
                 DistMatrix<T,MC,MR>& A );

//
// Geru (GEneral Rank-one Unconjugated update):
//
// A := alpha x y^T + A
//

// Serial version
template<typename T>
void
Geru( T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A );

// Parallel version
template<typename T>
void
Geru
( T alpha, const DistMatrix<T,MC,MR>& x, const DistMatrix<T,MC,MR>& y,
                 DistMatrix<T,MC,MR>& A );

//
// Hemv (HErmitian Matrix-Vector multiply):
//
// Implicitly performs
//   y := alpha A x + beta y,
// where only the triangle specified by 'shape' is referenced and the other
// triangle is implied by the Hermitian assumption.

// Serial version
template<typename T>
void
Hemv( Shape shape,
      T alpha, const Matrix<T>& A, const Matrix<T>& x, T beta, Matrix<T>& y );

// Parallel version
template<typename T>
void
Hemv
( Shape shape,
  T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& x,
  T beta,        DistMatrix<T,MC,MR>& y );

//
// Her (HErmitian Rank-one update):
//
// Implicitly performs
//   A := alpha x x^H + A,
// where only the triangle specified by 'shape' is updated.
//

// Serial version
template<typename T>
void
Her( Shape shape, T alpha, const Matrix<T>& x, Matrix<T>& A );

// Parallel version
template<typename T>
void
Her
( Shape shape, T alpha, const DistMatrix<T,MC,MR>& x, DistMatrix<T,MC,MR>& A );

//
// Her2 (HErmitian Rank-2 update):
//
// Implicitly performs
//   A := alpha ( x y^H + y x^H ) + A,
// where only the triangle specified by 'shape' is updated.
//

// Serial version
template<typename T>
void
Her2
( Shape shape, T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A );

// Parallel version
template<typename T>
void
Her2
( Shape shape,
  T alpha, const DistMatrix<T,MC,MR>& x, const DistMatrix<T,MC,MR>& y,
                 DistMatrix<T,MC,MR>& A );

//
// Symv (SYmmetric Matrix-Vector multiply):
//
// Implicitly performs
//   y := alpha A x + beta y,
// where only the triangle specified by 'shape' is referenced and the other
// triangle is implied by the symmetry assumption.
//

// Serial version
template<typename T>
void
Symv
( Shape shape,
  T alpha, const Matrix<T>& A, const Matrix<T>& x, T beta, Matrix<T>& y );

// Parallel version
template<typename T>
void
Symv
( Shape shape,
  T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& x,
  T beta,        DistMatrix<T,MC,MR>& y );

//
// Syr (SYmmetric Rank-one update):
//
// Implicitly performs the update
//   A := alpha x x^T + A,
// where only the triangle specified by 'shape' is updated.
//

// Serial version
template<typename T>
void
Syr( Shape shape, T alpha, const Matrix<T>& x, Matrix<T>& A );

// Parallel version
template<typename T>
void
Syr
( Shape shape,
  T alpha, const DistMatrix<T,MC,MR>& x, DistMatrix<T,MC,MR>& A );

//
// Syr2 (SYmmetric Rank-2 update):
//
// Implicitly perform the update
//   A := alpha ( x y^T + y x^T ) + A
// where only the triangle specified by 'shape' is updated.
//

// Serial version
template<typename T>
void
Syr2
( Shape shape, T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A );

// Parallel version
template<typename T>
void
Syr2
( Shape shape,
  T alpha, const DistMatrix<T,MC,MR>& x, const DistMatrix<T,MC,MR>& y,
                 DistMatrix<T,MC,MR>& A );

//
// Trmv (TRiangular Matrix-Vector multiply):
//
// Performs the update
//   x := orientation( A ) x,
// where 'shape' determines whether or not A is to be implicitly treated as 
// lower or upper triangular, and 'diagonal' specifies whether it has an 
// implicit unit diagonal.
//

// Serial version
template<typename T>
void
Trmv
( Shape shape, Orientation orientation, Diagonal diagonal,
  const Matrix<T>& A, Matrix<T>& x );

// Parallel version
template<typename T>
void
Trmv
( Shape shape, Orientation orientation, Diagonal diagonal,
  const DistMatrix<T,MC,MR>& A, DistMatrix<T,MC,MR>& x );

//
// Trsv (TRiangular Solve with a Vector):
//
// Performs the update
//   x := orientation( A )^-1 x,
// where 'shape' determines whether or not A is to be implicitly treated as 
// lower or upper triangular, and 'diagonal' specifies whether it has an 
// implicit unit diagonal.
//

// Serial version
template<typename T>
void
Trsv
( Shape shape, Orientation orientation, Diagonal diagonal,
  const Matrix<T>& A, Matrix<T>& x );

// Parallel version
template<typename T>
void
Trsv
( Shape shape, Orientation orientation, Diagonal diagonal,
  const DistMatrix<T,MC,MR>& A, DistMatrix<T,MC,MR>& x );

//----------------------------------------------------------------------------//
// Level 3 BLAS                                                               //
//----------------------------------------------------------------------------//

//
// Gemm (GEneral Matrix-Matrix multiplication):
//
// C := alpha orientationOfA( A ) orientationOfB( B ) + beta C
//

// Serial version
template<typename T>
void
Gemm
( Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const Matrix<T>& A, const Matrix<T>& B, T beta, Matrix<T>& C );

// Parallel version
template<typename T>
void
Gemm
( Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C );

//
// Hemm (HErmitian Matrix-Matrix multiplication):
//
// Performs the update
//   C := alpha A B + beta C,  { side = Left }
// or
//   C := alpha B A + beta C,  { side = Right }
// where only the triangle of 'A' specified by 'shape' is referenced, and the
// other triangle is implied by the Hermitian assumption.
//

// Serial version
template<typename T>
void
Hemm
( Side side, Shape shape,
  T alpha, const Matrix<T>& A, const Matrix<T>& B, T beta, Matrix<T>& C );

// Parallel version
template<typename T>
void
Hemm
( Side side, Shape shape,
  T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C );

//
// Her2k (HErmitian Rank-2K update):
//
// Performs the update
//   C := alpha ( A B^H + B A^H ) + beta C, { orientation = Normal }
// or
//   C := alpha ( A^H B + B^H A ) + beta C, { orientation = ConjugateTranspose }
// where only the triangle of C specified by 'shape' is updated.
//

// Serial version
template<typename T>
void
Her2k
( Shape shape, Orientation orientation,
  T alpha, const Matrix<T>& A, const Matrix<T>& B, T beta, Matrix<T>& C );

// Parallel version
template<typename T>
void
Her2k
( Shape shape, Orientation orientation,
  T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C );

//
// Herk (HErmitian Rank-K update):
//
// Performs the update
//   C := alpha A B^H + beta C,  { orientation = Normal }
// or
//   C := alpha A^H B + beta C,  { orientation = ConjugateTranspose }
// where only the triangle of C specified by 'shape' is updated.
//

// Serial version
template<typename T>
void
Herk
( Shape shape, Orientation orientation,
  T alpha, const Matrix<T>& A, T beta, Matrix<T>& C );

// Parallel version
template<typename T>
void
Herk
( Shape shape, Orientation orientation,
  T alpha, const DistMatrix<T,MC,MR>& A, T beta, DistMatrix<T,MC,MR>& C );

//
// Symm (SYmmetric Matrix-Matrix multiplication):
//
// Performs the update
//   C := alpha A B + beta C,  { side = Left }
// or
//   C := alpha B A + beta C,  { side = Right }
// where only the triangle of A specified by 'shape' is referenced, and the 
// other triangle is implied by the symmetry assumption.
//

// Serial version
template<typename T>
void
Symm
( Side side, Shape shape,
  T alpha, const Matrix<T>& A, const Matrix<T>& B, T beta, Matrix<T>& C ); 

// Parallel version
template<typename T>
void
Symm
( Side side, Shape shape,
  T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C );

//
// Syr2k (SYmmetric Rank-2K update):
//
// Performs the update
//   C := alpha ( A B^H + B A^H ) + beta C,  { orientation = Normal }
// or
//   C := alpha ( A^H B + B^H A ) + beta C,  { orientation = Transpose }
// where only the triangle of C specified by 'shape' is updated.
//

// Serial version
template<typename T>
void
Syr2k
( Shape shape, Orientation orientation,
  T alpha, const Matrix<T>& A, const Matrix<T>& B, T beta, Matrix<T>& C );

// Parallel version
template<typename T>
void
Syr2k
( Shape shape, Orientation orientation,
  T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C );

//
// Syrk (SYmmetric Rank-K update):
//
// Performs the update
//   C := alpha A B^H + beta C,  { orientation = Normal }
// or
//   C := alpha A^H B + beta C,  { orientation = Transpose }
// where only the triangle of C specified by 'shape' is updated.
//

// Serial version
template<typename T>
void
Syrk
( Shape shape, Orientation orientation,
  T alpha, const Matrix<T>& A, T beta, Matrix<T>& C );

// Parallel version
template<typename T>
void
Syrk
( Shape shape, Orientation orientation,
  T alpha, const DistMatrix<T,MC,MR>& A, T beta, DistMatrix<T,MC,MR>& C );

//
// Trmm (TRiangular Matrix-Matrix multiplication):
//
// Performs the update
//   B := alpha orientation( A ) B,  { side = Left }
// or
//   B := alpha B orientation( A ),  { side = Right }
// where 'shape' determines whether A is assumed to be upper or lower 
// triangular and 'diagonal' determines whether A has an implicit unit
// diagonal.
//

// Serial version
template<typename T>
void
Trmm
( Side side, Shape shape, Orientation orientation, Diagonal diagonal,
  T alpha, const Matrix<T>& A, Matrix<T>& B );

// Parallel version
template<typename T>
void
Trmm
( Side side, Shape shape, Orientation orientation, Diagonal diagonal,
  T alpha, const DistMatrix<T,MC,MR>& A, DistMatrix<T,MC,MR>& B );

//
// Trsm (TRiangular Solve with Multiple right-hand sides):
//
// Performs the update
//   B := alpha orientation( A )^-1 B,  { side = Left }
// or
//   B := alpha B orientation( A )^-1,  { side = Right }
// where 'shape' determines whether A is assumed to be upper or lower
// triangular and 'diagonal' determines whether A has an implicit unit
// diagonal.
//

// Serial version
template<typename T>
void
Trsm
( Side side, Shape shape, Orientation orientation, Diagonal diagonal,
  T alpha, const Matrix<T>& A, Matrix<T>& B ); 
        
// Parallel version
template<typename T>
void
Trsm
( Side side, Shape shape, Orientation orientation, Diagonal diagonal,
  T alpha, const DistMatrix<T,MC,MR>& A, DistMatrix<T,MC,MR>& B );

} // blas
} // elemental

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
// Local BLAS: Level 1                                                        //
//----------------------------------------------------------------------------//

template<typename T>
inline void
elemental::blas::Axpy
( T alpha, const Matrix<T>& X, Matrix<T>& Y )
{
#ifndef RELEASE
    PushCallStack("blas::Axpy");
#endif
    // If X and Y are vectors, we can allow one to be a column and the other
    // to be a row. Otherwise we force X and Y to be the same dimension.
    if( (X.Height()==1 || X.Width()==1) && (Y.Height()==1 || Y.Width()==1) )
    {
        const unsigned XLength = ( X.Width()==1 ? X.Height() : X.Width() );
#ifndef RELEASE
        const unsigned YLength = ( Y.Width()==1 ? Y.Height() : Y.Width() );
        if( XLength != YLength )
            throw std::logic_error( "Nonconformal Axpy." );
#endif
        if( X.Width()==1 && Y.Width()==1 )
        {
            elemental::wrappers::blas::Axpy
            ( XLength, alpha, X.LockedBuffer(0,0), 1, Y.Buffer(0,0), 1 );
        }
        else if( X.Width()==1 )
        {
            elemental::wrappers::blas::Axpy
            ( XLength, alpha, X.LockedBuffer(0,0), 1, Y.Buffer(0,0), Y.LDim() );
        }
        else if( Y.Width()==1 )
        {
            elemental::wrappers::blas::Axpy
            ( XLength, alpha, X.LockedBuffer(0,0), X.LDim(), Y.Buffer(0,0), 1 );
        }
        else
        {
            elemental::wrappers::blas::Axpy
            ( XLength, alpha, 
              X.LockedBuffer(0,0), X.LDim(), Y.Buffer(0,0), Y.LDim() );
        }
    }
    else 
    {
#ifndef RELEASE
        if( X.Height() != Y.Height() || X.Width() != Y.Width() )
            throw std::logic_error( "Nonconformal Axpy." );
#endif
        if( X.Width() <= X.Height() )
        {
            for( int j=0; j<X.Width(); ++j )
            {
                elemental::wrappers::blas::Axpy
                ( X.Height(), alpha, X.LockedBuffer(0,j), 1, Y.Buffer(0,j), 1 );
            }
        }
        else
        {
            for( int i=0; i<X.Height(); ++i )
            {
                elemental::wrappers::blas::Axpy
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
elemental::blas::Copy
( const Matrix<T>& A, Matrix<T>& B )
{
#ifndef RELEASE
    PushCallStack("blas::Copy");
#endif
    B = A;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline T
elemental::blas::Dot
( const Matrix<T>& x, const Matrix<T>& y )
{
#ifndef RELEASE
    PushCallStack("blas::Dot");
    if( (x.Height() != 1 && x.Width() != 1) ||
        (y.Height() != 1 && y.Width() != 1) )
        throw std::logic_error( "Expected vector inputs." );
    int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
    int yLength = ( y.Width() == 1 ? y.Height() : y.Width() );
    if( xLength != yLength )
        throw std::logic_error( "x and y must be the same length." );
#endif
    T dotProduct;
    if( x.Width() == 1 && y.Width() == 1 )
    {
        dotProduct = wrappers::blas::Dot
                     ( x.Height(), x.LockedBuffer(), 1,
                                   y.LockedBuffer(), 1 );
    }
    else if( x.Width() == 1 )
    {
        dotProduct = wrappers::blas::Dot
                     ( x.Height(), x.LockedBuffer(), 1,
                                   y.LockedBuffer(), y.LDim() );
    }
    else if( y.Width() == 1 )
    {
        dotProduct = wrappers::blas::Dot
                     ( x.Width(), x.LockedBuffer(), x.LDim(),
                                  y.LockedBuffer(), 1        );
    }
    else
    {
        dotProduct = wrappers::blas::Dot
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
elemental::blas::Dotc
( const Matrix<T>& x, const Matrix<T>& y )
{
#ifndef RELEASE
    PushCallStack("blas::Dotc");
    if( (x.Height() != 1 && x.Width() != 1) ||
        (y.Height() != 1 && y.Width() != 1) )
        throw std::logic_error( "Expected vector inputs." );
    int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
    int yLength = ( y.Width() == 1 ? y.Height() : y.Width() );
    if( xLength != yLength )
        throw std::logic_error( "x and y must be the same length." );
#endif
    T dotProduct;
    if( x.Width() == 1 && y.Width() == 1 )
    {
        dotProduct = wrappers::blas::Dotc
                     ( x.Height(), x.LockedBuffer(), 1,
                                   y.LockedBuffer(), 1 );
    }
    else if( x.Width() == 1 )
    {
        dotProduct = wrappers::blas::Dotc
                     ( x.Height(), x.LockedBuffer(), 1,
                                   y.LockedBuffer(), y.LDim() );
    }
    else if( y.Width() == 1 )
    {
        dotProduct = wrappers::blas::Dotc
                     ( x.Width(), x.LockedBuffer(), x.LDim(),
                                  y.LockedBuffer(), 1        );
    }
    else
    {
        dotProduct = wrappers::blas::Dotc
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
elemental::blas::Dotu
( const Matrix<T>& x, const Matrix<T>& y )
{
#ifndef RELEASE
    PushCallStack("blas::Dotu");
    if( (x.Height() != 1 && x.Width() != 1) ||
        (y.Height() != 1 && y.Width() != 1) )
        throw std::logic_error( "Expected vector inputs." );
    int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
    int yLength = ( y.Width() == 1 ? y.Height() : y.Width() );
    if( xLength != yLength )
        throw std::logic_error( "x and y must be the same length." );
#endif
    T dotProduct;
    if( x.Width() == 1 && y.Width() == 1 )
    {
        dotProduct = wrappers::blas::Dotu
                     ( x.Height(), x.LockedBuffer(), 1,
                                   y.LockedBuffer(), 1 );
    }
    else if( x.Width() == 1 )
    {
        dotProduct = wrappers::blas::Dotu
                     ( x.Height(), x.LockedBuffer(), 1,
                                   y.LockedBuffer(), y.LDim() );
    }
    else if( y.Width() == 1 )
    {
        dotProduct = wrappers::blas::Dotu
                     ( x.Width(), x.LockedBuffer(), x.LDim(),
                                  y.LockedBuffer(), 1        );
    }
    else
    {
        dotProduct = wrappers::blas::Dotu
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
elemental::blas::Nrm2
( const Matrix<R>& x )
{
#ifndef RELEASE
    PushCallStack("blas::Nrm2");
    if( x.Height() != 1 && x.Width() != 1 )
        throw std::logic_error( "Expected vector input." );
#endif
    R norm;
    if( x.Width() == 1 )
    {
        norm = wrappers::blas::Nrm2
               ( x.Height(), x.LockedBuffer(), 1 );
    }
    else
    {
        norm = wrappers::blas::Nrm2
               ( x.Width(), x.LockedBuffer(), x.LDim() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

#ifndef WITHOUT_COMPLEX
template<typename R>
inline R
elemental::blas::Nrm2
( const Matrix< std::complex<R> >& x )
{
#ifndef RELEASE
    PushCallStack("blas::Nrm2");
    if( x.Height() != 1 && x.Width() != 1 )
        throw std::logic_error( "Expected vector input." );
#endif
    R norm;
    if( x.Width() == 1 )
    {
        norm = wrappers::blas::Nrm2
               ( x.Height(), x.LockedBuffer(), 1 );
    }
    else
    {
        norm = wrappers::blas::Nrm2
               ( x.Width(), x.LockedBuffer(), x.LDim() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}
#endif

template<typename T>
inline void
elemental::blas::Scal
( T alpha, Matrix<T>& X )
{
#ifndef RELEASE
    PushCallStack("blas::Scal");
#endif
    if( alpha != (T)1 )
    {
        if( alpha == (T)0 )
        {
            for( int j=0; j<X.Width(); ++j )
                for( int i=0; i<X.Height(); ++i )
                    X(i,j) = 0;
        }
        else
        {
            for( int j=0; j<X.Width(); ++j )
            {
                wrappers::blas::Scal
                ( X.Height(), alpha, X.Buffer(0,j), 1 );
            }
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

//----------------------------------------------------------------------------//
// Local BLAS: Level 1 (extensions)                                           //
//----------------------------------------------------------------------------//

// Default case is for real datatypes
template<typename R>
inline void
elemental::blas::Conj
( Matrix<R>& A )
{ }

#ifndef WITHOUT_COMPLEX
// Specialization is to complex datatypes
template<typename R>
inline void
elemental::blas::Conj
( Matrix< std::complex<R> >& A )
{
#ifndef RELEASE
    PushCallStack("blas::Conj (in-place)");
#endif
    const int m = A.Height();
    const int n = A.Width();
    for( int j=0; j<n; ++j )
        for( int i=0; i<m; ++i )
            A(i,j) = elemental::Conj( A(i,j) );
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif

template<typename T>
inline void
elemental::blas::Conj
( const Matrix<T>& A, Matrix<T>& B )
{
#ifndef RELEASE
    PushCallStack("blas::Conj");
#endif
    const int m = A.Height();
    const int n = A.Width();
    B.ResizeTo( m, n );
    for( int j=0; j<n; ++j )
        for( int i=0; i<m; ++i )
            B(i,j) = elemental::Conj( A(i,j) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::blas::ConjTrans
( const Matrix<T>& A, Matrix<T>& B )
{
#ifndef RELEASE
    PushCallStack("blas::ConjTrans");
#endif
    const int m = A.Height();
    const int n = A.Width();
    B.ResizeTo( n, m );
    for( int j=0; j<n; ++j )
        for( int i=0; i<m; ++i )
            B(j,i) = elemental::Conj( A(i,j) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::blas::Trans
( const Matrix<T>& A, Matrix<T>& B )
{
#ifndef RELEASE
    PushCallStack("blas::Trans");
#endif
    const int m = A.Height();
    const int n = A.Width();
    B.ResizeTo( n, m );
    for( int j=0; j<n; ++j )
        for( int i=0; i<m; ++i )
            B(j,i) = A(i,j);
#ifndef RELEASE
    PopCallStack();
#endif
}

//----------------------------------------------------------------------------//
// Local BLAS: Level 2                                                        //
//----------------------------------------------------------------------------//

template<typename T>
inline void
elemental::blas::Gemv
( Orientation orientation,
  T alpha, const Matrix<T>& A, const Matrix<T>& x, T beta, Matrix<T>& y )
{
#ifndef RELEASE
    PushCallStack("blas::Gemv");
    if( ( x.Height() != 1 && x.Width() != 1 ) ||
        ( y.Height() != 1 && y.Width() != 1 ) )
    {
        std::ostringstream msg;
        msg << "x and y must be vectors: " << std::endl
            << "  x ~ " << x.Height() << " x " << x.Width() << std::endl
            << "  y ~ " << y.Height() << " x " << y.Width() << std::endl;
        throw std::logic_error( msg.str() );
    }
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    const int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
    if( orientation == Normal )
    {
        if( A.Height() != yLength || A.Width() != xLength )
        {
            std::ostringstream msg;
            msg << "A must conform with x and y:" << std::endl
                << "  A ~ " << A.Height() << " x " << A.Width() << std::endl
                << "  x ~ " << x.Height() << " x " << x.Width() << std::endl
                << "  y ~ " << y.Height() << " x " << y.Width() << std::endl;
            throw std::logic_error( msg.str() );
        }
    }
    else
    {
        if( A.Width() != yLength || A.Height() != xLength )
        {
            std::ostringstream msg;
            msg << "A must conform with x and y:" << std::endl
                << "  A ~ " << A.Height() << " x " << A.Width() << std::endl
                << "  x ~ " << x.Height() << " x " << x.Width() << std::endl
                << "  y ~ " << y.Height() << " x " << y.Width() << std::endl;
            throw std::logic_error( msg.str() );
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
        wrappers::blas::Gemv
        ( transChar, m, n, 
          alpha, A.LockedBuffer(), A.LDim(), x.LockedBuffer(), incx, 
          beta,  y.Buffer(), incy );
    }
    else
    {
        blas::Scal( beta, y );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::blas::Ger
( T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("blas::Ger");
    if( ( x.Height() != 1 && x.Width() != 1 ) ||
        ( y.Height() != 1 && y.Width() != 1 ) )
        throw std::logic_error( "x and y must be vectors." );
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    const int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
    if( xLength != A.Height() || yLength != A.Width() )
    {
        std::ostringstream msg;
        msg << "Nonconformal Ger: " << std::endl
            << "  x ~ " << x.Height() << " x " << x.Width() << std::endl
            << "  y ~ " << y.Height() << " x " << y.Width() << std::endl
            << "  A ~ " << A.Height() << " x " << A.Width() << std::endl;
        throw std::logic_error( msg.str() );
    }
#endif
    const int m = A.Height(); 
    const int n = A.Width();
    const int incx = ( x.Width()==1 ? 1 : x.LDim() );
    const int incy = ( y.Width()==1 ? 1 : y.LDim() );
    wrappers::blas::Ger
    ( m, n, alpha, x.LockedBuffer(), incx, y.LockedBuffer(), incy, 
                   A.Buffer(), A.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::blas::Gerc
( T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("blas::Gerc");
    if( ( x.Height() != 1 && x.Width() != 1 ) ||
        ( y.Height() != 1 && y.Width() != 1 ) )
        throw std::logic_error( "x and y must be vectors." );
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    const int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
    if( xLength != A.Height() || yLength != A.Width() )
        throw std::logic_error( "Nonconformal Gerc." );
#endif
    const int m = A.Height(); 
    const int n = A.Width();
    const int incx = ( x.Width()==1 ? 1 : x.LDim() );
    const int incy = ( y.Width()==1 ? 1 : y.LDim() );
    wrappers::blas::Gerc
    ( m, n, alpha, x.LockedBuffer(), incx, y.LockedBuffer(), incy, 
                   A.Buffer(), A.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::blas::Geru
( T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("blas::Geru");
    if( ( x.Height() != 1 && x.Width() != 1 ) ||
        ( y.Height() != 1 && y.Width() != 1 ) )
        throw std::logic_error( "x and y must be vectors." );
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    const int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
    if( xLength != A.Height() || yLength != A.Width() )
        throw std::logic_error( "Nonconformal Geru." );
#endif
    const int m = A.Height(); 
    const int n = A.Width();
    const int incx = ( x.Width()==1 ? 1 : x.LDim() );
    const int incy = ( y.Width()==1 ? 1 : y.LDim() );
    wrappers::blas::Geru
    ( m, n, alpha, x.LockedBuffer(), incx, y.LockedBuffer(), incy, 
                   A.Buffer(), A.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::blas::Hemv
( Shape shape,
  T alpha, const Matrix<T>& A, const Matrix<T>& x, T beta, Matrix<T>& y )
{
#ifndef RELEASE
    PushCallStack("blas::Hemv");
    if( A.Height() != A.Width() )
        throw std::logic_error( "A must be square." );
    if( ( x.Height() != 1 && x.Width() != 1 ) ||
        ( y.Height() != 1 && y.Width() != 1 ) )
        throw std::logic_error( "x and y must be vectors." );
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    const int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
    if( A.Height() != xLength || A.Height() != yLength )
        throw std::logic_error( "A must conform with x and y." );
#endif
    const char uploChar = ShapeToChar( shape );
    const int m = A.Height();
    const int incx = ( x.Width()==1 ? 1 : x.LDim() );
    const int incy = ( y.Width()==1 ? 1 : y.LDim() );
    wrappers::blas::Hemv
    ( uploChar, m,
      alpha, A.LockedBuffer(), A.LDim(), x.LockedBuffer(), incx,
      beta,  y.Buffer(), incy );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::blas::Her
( Shape shape, T alpha, const Matrix<T>& x, Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("blas::Her");
    if( A.Height() != A.Width() )
        throw std::logic_error( "A must be square." );
    if( x.Width() != 1 && x.Height() != 1 )
        throw std::logic_error( "x must be a vector." );
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    if( xLength != A.Height() )
        throw std::logic_error( "x must conform with A." );
#endif
    const char uploChar = ShapeToChar( shape );
    const int m = A.Height();
    const int incx = ( x.Width()==1 ? 1 : x.LDim() );
    wrappers::blas::Her
    ( uploChar, m,
      alpha, x.LockedBuffer(), incx, A.Buffer(), A.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::blas::Her2
( Shape shape,
  T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("blas::Her2");
    if( A.Height() != A.Width() )
        throw std::logic_error( "A must be square." );
    if( (x.Width() != 1 && x.Height() != 1) || 
        (y.Width() != 1 && y.Height() != 1) )
        throw std::logic_error( "x and y must be vectors." );
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    const int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
    if( xLength != A.Height() || yLength != A.Height() )
        throw std::logic_error( "x and y must conform with A." );
#endif
    const char uploChar = ShapeToChar( shape );
    const int m = A.Height();
    const int incx = ( x.Width()==1 ? 1 : x.LDim() );
    const int incy = ( y.Width()==1 ? 1 : y.LDim() );
    wrappers::blas::Her2
    ( uploChar, m,
      alpha, x.LockedBuffer(), incx, y.LockedBuffer(), incy,
             A.Buffer(), A.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::blas::Symv
( Shape shape,
  T alpha, const Matrix<T>& A, const Matrix<T>& x, T beta, Matrix<T>& y )
{
#ifndef RELEASE
    PushCallStack("blas::Symv");
    if( A.Height() != A.Width() )
        throw std::logic_error( "A must be square." );
    if( ( x.Height() != 1 && x.Width() != 1 ) ||
        ( y.Height() != 1 && y.Width() != 1 ) )
        throw std::logic_error( "x and y must be vectors." );
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    const int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
    if( A.Height() != xLength || A.Height() != yLength )
        throw std::logic_error( "A must conform with x and y." );
#endif
    const char uploChar = ShapeToChar( shape );
    const int m = A.Height();
    const int incx = ( x.Width()==1 ? 1 : x.LDim() );
    const int incy = ( y.Width()==1 ? 1 : y.LDim() );
    wrappers::blas::Symv
    ( uploChar, m, 
      alpha, A.LockedBuffer(), A.LDim(), x.LockedBuffer(), incx, 
      beta,  y.Buffer(), incy );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::blas::Syr
( Shape shape, T alpha, const Matrix<T>& x, Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("blas::Syr");
    if( A.Height() != A.Width() )
        throw std::logic_error( "A must be square." );
    if( x.Width() != 1 && x.Height() != 1 )
        throw std::logic_error( "x must be a vector." );
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    if( xLength != A.Height() )
        throw std::logic_error( "x must conform with A." );
#endif
    const char uploChar = ShapeToChar( shape );
    const int m = A.Height();
    const int incx = ( x.Width()==1 ? 1 : x.LDim() );
    wrappers::blas::Syr
    ( uploChar, m,
      alpha, x.LockedBuffer(), incx, A.Buffer(), A.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::blas::Syr2
( Shape shape, T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("blas::Syr2");
    if( A.Height() != A.Width() )
        throw std::logic_error( "A must be square." );
    if( (x.Width() != 1 && x.Height() != 1) || 
        (y.Width() != 1 && y.Height() != 1) )
        throw std::logic_error( "x and y must be vectors." );
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    const int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
    if( xLength != A.Height() || yLength != A.Height() )
        throw std::logic_error( "x and y must conform with A." );
#endif
    const char uploChar = ShapeToChar( shape );
    const int m = A.Height();
    const int incx = ( x.Width()==1 ? 1 : x.LDim() );
    const int incy = ( y.Width()==1 ? 1 : y.LDim() );
    wrappers::blas::Syr2
    ( uploChar, m,
      alpha, x.LockedBuffer(), incx, y.LockedBuffer(), incy,
             A.Buffer(), A.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::blas::Trmv
( Shape shape, Orientation orientation, Diagonal diagonal,
  const Matrix<T>& A, Matrix<T>& x )
{
#ifndef RELEASE
    PushCallStack("blas::Trmv");
    if( x.Height() != 1 && x.Width() != 1 )
        throw std::logic_error( "x must be a vector." );
    if( A.Height() != A.Width() )
        throw std::logic_error( "A must be square." );
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    if( xLength != A.Height() )
        throw std::logic_error( "x must conform with A." );
#endif
    const char uploChar = ShapeToChar( shape );
    const char transChar = OrientationToChar( orientation );
    const char diagChar = DiagonalToChar( diagonal );
    const int m = A.Height();
    const int incx = ( x.Width()==1 ? 1 : x.LDim() );
    wrappers::blas::Trmv
    ( uploChar, transChar, diagChar, m,
      A.LockedBuffer(), A.LDim(), x.Buffer(), incx );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::blas::Trsv
( Shape shape, Orientation orientation, Diagonal diagonal,
  const Matrix<T>& A, Matrix<T>& x )
{
#ifndef RELEASE
    PushCallStack("blas::Trsv");
    if( x.Height() != 1 && x.Width() != 1 )
        throw std::logic_error( "x must be a vector." );
    if( A.Height() != A.Width() )
        throw std::logic_error( "A must be square." );
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    if( xLength != A.Height() )
        throw std::logic_error( "x must conform with A." );
#endif
    const char uploChar = ShapeToChar( shape );
    const char transChar = OrientationToChar( orientation );
    const char diagChar = DiagonalToChar( diagonal );
    const int m = A.Height();
    const int incx = ( x.Width()==1 ? 1 : x.LDim() );
    wrappers::blas::Trsv
    ( uploChar, transChar, diagChar, m,
      A.LockedBuffer(), A.LDim(), x.Buffer(), incx );
#ifndef RELEASE
    PopCallStack();
#endif
}

//----------------------------------------------------------------------------//
// Local BLAS: Level 3                                                        //
//----------------------------------------------------------------------------//

template<typename T>
inline void
elemental::blas::Gemm
( Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const Matrix<T>& A, const Matrix<T>& B, T beta, Matrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("blas::Gemm");
    if( orientationOfA == Normal && orientationOfB == Normal )
    {
        if( A.Height() != C.Height() ||
            B.Width()  != C.Width()  ||
            A.Width()  != B.Height() )
            throw std::logic_error( "Nonconformal GemmNN." );
    }
    else if( orientationOfA == Normal )
    {
        if( A.Height() != C.Height() ||
            B.Height() != C.Width()  ||
            A.Width()  != B.Width() )
            throw std::logic_error( "Nonconformal GemmN(T/C)." );
    }
    else if( orientationOfB == Normal )
    {
        if( A.Width()  != C.Height() ||
            B.Width()  != C.Width()  ||
            A.Height() != B.Height() )
            throw std::logic_error( "Nonconformal Gemm(T/C)N." );
    }
    else
    {
        if( A.Width()  != C.Height() ||
            B.Height() != C.Width()  ||
            A.Height() != B.Width() )
            throw std::logic_error( "Nonconformal Gemm(T/C)(T/C)." );
    }
#endif
    const char transA = OrientationToChar( orientationOfA );
    const char transB = OrientationToChar( orientationOfB );
    const int m = C.Height();
    const int n = C.Width();
    const int k = ( orientationOfA == Normal ? A.Width() : A.Height() );
    if( k != 0 )
    {
        wrappers::blas::Gemm
        ( transA, transB, m, n, k, 
          alpha, A.LockedBuffer(), A.LDim(), B.LockedBuffer(), B.LDim(),
          beta,  C.Buffer(),       C.LDim() );
    }
    else
    {
        blas::Scal( beta, C );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::blas::Hemm
( Side side, Shape shape,
  T alpha, const Matrix<T>& A, const Matrix<T>& B, T beta, Matrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("blas::Hemm");
#endif
    const char sideChar = SideToChar( side );
    const char shapeChar = ShapeToChar( shape );
    wrappers::blas::Hemm
    ( sideChar, shapeChar, C.Height(), C.Width(),
      alpha, A.LockedBuffer(), A.LDim(), 
             B.LockedBuffer(), B.LDim(),
      beta,  C.Buffer(),       C.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::blas::Her2k
( Shape shape, Orientation orientation,
  T alpha, const Matrix<T>& A, const Matrix<T>& B, T beta, Matrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("blas::Her2k");
    if( orientation == Normal )
    {
        if( A.Height() != C.Height() || A.Height() != C.Width() ||
            B.Height() != C.Height() ||B.Height() != C.Width() )
            throw std::logic_error( "Nonconformal Her2k." );
    }
    else if( orientation == ConjugateTranspose )
    {
        if( A.Width() != C.Height() || A.Width() != C.Width() ||
            B.Width() != C.Height() || B.Width() != C.Width() )
            throw std::logic_error( "Nonconformal Her2k." );
    }
    else
        throw std::logic_error
        ( "Her2k only accepts Normal and ConjugateTranspose options." );
#endif
    const char uplo = ShapeToChar( shape );
    const char trans = OrientationToChar( orientation );
    const int k = ( orientation == Normal ? A.Width() : A.Height() );
    wrappers::blas::Her2k
    ( uplo, trans, C.Height(), k, 
      alpha, A.LockedBuffer(), A.LDim(), 
             B.LockedBuffer(), B.LDim(),
      beta,  C.Buffer(),       C.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::blas::Herk
( Shape shape, Orientation orientation,
  T alpha, const Matrix<T>& A, T beta, Matrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("blas::Herk");
    if( orientation == Normal )
    {
        if( A.Height() != C.Height() || A.Height() != C.Width() )
            throw std::logic_error( "Nonconformal Herk." );
    }
    else if( orientation == ConjugateTranspose )
    {
        if( A.Width() != C.Height() || A.Width() != C.Width() )
            throw std::logic_error( "Nonconformal Herk." );
    }
    else
        throw std::logic_error
        ( "Herk only accepts Normal and ConjugateTranpose options." );
#endif
    const char uplo = ShapeToChar( shape );
    const char trans = OrientationToChar( orientation );
    const int k = ( orientation == Normal ? A.Width() : A.Height() );
    wrappers::blas::Herk
    ( uplo, trans, C.Height(), k, 
      alpha, A.LockedBuffer(), A.LDim(), 
      beta,  C.Buffer(),       C.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::blas::Symm
( Side side, Shape shape,
  T alpha, const Matrix<T>& A, const Matrix<T>& B, T beta, Matrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("blas::Symm");
#endif
    const char sideChar = SideToChar( side );
    const char shapeChar = ShapeToChar( shape );
    wrappers::blas::Symm
    ( sideChar, shapeChar, C.Height(), C.Width(),
      alpha, A.LockedBuffer(), A.LDim(), 
             B.LockedBuffer(), B.LDim(),
      beta,  C.Buffer(),       C.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::blas::Syr2k
( Shape shape, Orientation orientation,
  T alpha, const Matrix<T>& A, const Matrix<T>& B, T beta, Matrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("blas::Syr2k");
    if( orientation == Normal )
    {
        if( A.Height() != C.Height() || A.Height() != C.Width() ||
            B.Height() != C.Height() ||B.Height() != C.Width()    )
            throw std::logic_error( "Nonconformal Syr2k." );
    }
    else if( orientation == Transpose )
    {
        if( A.Width() != C.Height() || A.Width() != C.Width() ||
            B.Width() != C.Height() || B.Width() != C.Width()   )
            throw std::logic_error( "Nonconformal Syr2k." );
    }
    else
        throw std::logic_error
        ( "Syr2k only accepts Normal and Tranpose options." );
#endif
    const char uplo = ShapeToChar( shape );
    const char trans = OrientationToChar( orientation );
    const int k = ( orientation == Normal ? A.Width() : A.Height() );
    wrappers::blas::Syr2k
    ( uplo, trans, C.Height(), k, 
      alpha, A.LockedBuffer(), A.LDim(), 
             B.LockedBuffer(), B.LDim(),
      beta,  C.Buffer(),       C.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::blas::Syrk
( Shape shape, Orientation orientation,
  T alpha, const Matrix<T>& A, T beta, Matrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("blas::Syrk");
    if( orientation == Normal )
    {
        if( A.Height() != C.Height() || A.Height() != C.Width() )
            throw std::logic_error( "Nonconformal Syrk." );
    }
    else if( orientation == Transpose )
    {
        if( A.Width() != C.Height() || A.Width() != C.Width() )
            throw std::logic_error( "Nonconformal Syrk." );
    }
    else
        throw std::logic_error
        ( "Syrk only accepts Normal and Tranpose options." );
#endif
    const char uplo = ShapeToChar( shape );
    const char trans = OrientationToChar( orientation );
    const int k = ( orientation == Normal ? A.Width() : A.Height() );
    wrappers::blas::Syrk
    ( uplo, trans, C.Height(), k, 
      alpha, A.LockedBuffer(), A.LDim(), 
      beta,  C.Buffer(),       C.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::blas::Trmm
( Side side, Shape shape, 
  Orientation orientation, Diagonal diagonal,
  T alpha, const Matrix<T>& A, Matrix<T>& B )
{
#ifndef RELEASE
    PushCallStack("blas::Trmm");
    if( A.Height() != A.Width() )
        throw std::logic_error( "Triangular matrix must be square." );
    if( side == Left )
    {
        if( A.Height() != B.Height() )
            throw std::logic_error( "Nonconformal Trmm." );
    }
    else
    {
        if( A.Height() != B.Width() )
            throw std::logic_error( "Nonconformal Trmm." );
    }
#endif
    const char sideChar = SideToChar( side );
    const char uploChar = ShapeToChar( shape );
    const char transChar = OrientationToChar( orientation );
    const char diagChar = DiagonalToChar( diagonal );
    wrappers::blas::Trmm
    ( sideChar, uploChar, transChar, diagChar, B.Height(), B.Width(),
      alpha, A.LockedBuffer(), A.LDim(), B.Buffer(), B.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::blas::Trsm
( Side side, Shape shape,
  Orientation orientation,Diagonal diagonal,
  T alpha, const Matrix<T>& A, Matrix<T>& B )
{
#ifndef RELEASE
    PushCallStack("blas::Trsm");
    if( A.Height() != A.Width() )
        throw std::logic_error( "Triangular matrix must be square." );
    if( side == Left )
    {
        if( A.Height() != B.Height() )
            throw std::logic_error( "Nonconformal Trsm." );
    }
    else
    {
        if( A.Height() != B.Width() )
            throw std::logic_error( "Nonconformal Trsm." );
    }
#endif
    const char sideChar = SideToChar( side );
    const char uploChar = ShapeToChar( shape );
    const char transChar = OrientationToChar( orientation );
    const char diagChar = DiagonalToChar( diagonal );
    wrappers::blas::Trsm
    ( sideChar, uploChar, transChar, diagChar, B.Height(), B.Width(),
      alpha, A.LockedBuffer(), A.LDim(), B.Buffer(), B.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

//----------------------------------------------------------------------------//
// Distributed BLAS: Level 1                                                  //
//----------------------------------------------------------------------------//

template<typename T, elemental::Distribution U, elemental::Distribution V>
inline void
elemental::blas::Axpy
( T alpha, const DistMatrix<T,U,V>& X, DistMatrix<T,U,V>& Y )
{
#ifndef RELEASE
    PushCallStack("blas::Axpy");
    if( X.GetGrid() != Y.GetGrid() )
        throw std::logic_error
        ( "X and Y must be distributed over the same grid." );
    if( X.ColAlignment() != Y.ColAlignment() ||
        X.RowAlignment() != Y.RowAlignment() )
        throw std::logic_error( "Axpy requires X and Y be aligned." );
#endif
    blas::Axpy( alpha, X.LockedLocalMatrix(), Y.LocalMatrix() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, elemental::Distribution U, elemental::Distribution V,
                     elemental::Distribution W, elemental::Distribution Z >
inline void
elemental::blas::Copy
( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B )
{
#ifndef RELEASE
    PushCallStack("blas::Copy");
#endif
    B = A;
#ifndef RELEASE
    PopCallStack();
#endif
}

// Our extended Dotc is equivalent to our extended Dot, 
// but we are burdened with consistency
template<typename T, elemental::Distribution U, elemental::Distribution V,
                     elemental::Distribution W, elemental::Distribution Z >
inline T
elemental::blas::Dotc
( const DistMatrix<T,U,V>& x, const DistMatrix<T,W,Z>& y )
{
#ifndef RELEASE
    PushCallStack("blas::Dotc");
#endif
    T globalDot = blas::Dot( x, y );
#ifndef RELEASE
    PopCallStack();
#endif
    return globalDot;
}

template<typename T, elemental::Distribution U, elemental::Distribution V>
inline void
elemental::blas::Scal
( T alpha, DistMatrix<T,U,V>& A )
{
#ifndef RELEASE
    PushCallStack("blas::Scal");
#endif
    blas::Scal( alpha, A.LocalMatrix() );
#ifndef RELEASE
    PopCallStack();
#endif
}

//----------------------------------------------------------------------------//
// Distributed BLAS: Level 1 (extensions)                                     //
//----------------------------------------------------------------------------//

template<typename T, elemental::Distribution U, elemental::Distribution V>
inline void
elemental::blas::Conj
( DistMatrix<T,U,V>& A )
{
#ifndef RELEASE
    PushCallStack("blas::Conj (in-place)");
#endif
    blas::Conj( A.LocalMatrix() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, elemental::Distribution U, elemental::Distribution V,
                     elemental::Distribution W, elemental::Distribution Z >
inline void
elemental::blas::Conj
( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B )
{
#ifndef RELEASE
    PushCallStack("blas::Conj");
#endif
    B = A;
    blas::Conj( B ); 
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, elemental::Distribution U, elemental::Distribution V,
                     elemental::Distribution W, elemental::Distribution Z >
inline void
elemental::blas::ConjTrans
( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B )
{
#ifndef RELEASE
    PushCallStack("blas::ConjTrans");
#endif
    DistMatrix<T,Z,W> C( B.GetGrid() );
    if( B.ConstrainedColAlignment() )
        C.AlignRowsWith( B );
    if( B.ConstrainedRowAlignment() )
        C.AlignColsWith( B );

    C = A;

    if( !B.ConstrainedColAlignment() )
        B.AlignColsWith( C );
    if( !B.ConstrainedRowAlignment() )
        B.AlignRowsWith( C );

    B.ResizeTo( A.Width(), A.Height() );
    blas::ConjTrans( C.LockedLocalMatrix(), B.LocalMatrix() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, elemental::Distribution U, elemental::Distribution V,
                     elemental::Distribution W, elemental::Distribution Z >
inline void
elemental::blas::Trans
( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B )
{
#ifndef RELEASE
    PushCallStack("blas::Trans");
#endif
    DistMatrix<T,Z,W> C( B.GetGrid() );
    if( B.ConstrainedColAlignment() )
        C.AlignRowsWith( B );
    if( B.ConstrainedRowAlignment() )
        C.AlignColsWith( B );

    C = A;

    if( !B.ConstrainedColAlignment() )
        B.AlignColsWith( C );
    if( !B.ConstrainedRowAlignment() )
        B.AlignRowsWith( C );

    B.ResizeTo( A.Width(), A.Height() );
    blas::Trans( C.LockedLocalMatrix(), B.LocalMatrix() );
#ifndef RELEASE
    PopCallStack();
#endif
}

//----------------------------------------------------------------------------//
// Distributed BLAS: Level 2                                                  //
//----------------------------------------------------------------------------//

// Our extended Ger and Gerc are equivalent
template<typename T>
inline void
elemental::blas::Gerc
( T alpha, const DistMatrix<T,MC,MR>& x,
           const DistMatrix<T,MC,MR>& y,
                 DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("blas::Gerc");
#endif
    blas::Ger( alpha, x, y, A );
#ifndef RELEASE
    PopCallStack();
#endif
}

#endif /* ELEMENTAL_BLAS_HPP */

