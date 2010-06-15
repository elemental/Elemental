/*
   This file is part of elemental, a library for distributed-memory dense 
   linear algebra.

   Copyright (C) 2009-2010 Jack Poulson <jack.poulson@gmail.com>

   This program is released under the terms of the license contained in the 
   file LICENSE.
*/
#ifndef ELEMENTAL_LAPACK_HPP
#define ELEMENTAL_LAPACK_HPP 1

#include "elemental/blas.hpp"
#include "elemental/wrappers/lapack.hpp"

namespace elemental {
namespace lapack {

//----------------------------------------------------------------------------//
// Chol:                                                                      //
//                                                                            //
// Overwrite a triangle of A with the Cholesky factor of A. 'shape'           //
// determines whether it is the upper or lower triangle.                      //
//----------------------------------------------------------------------------//

// Serial version
template<typename T>
void
Chol
( Shape shape, Matrix<T>& A );

// Parallel version
template<typename T>
void
Chol
( Shape shape, DistMatrix<T,MC,MR>& A );

//----------------------------------------------------------------------------//
// GaussElim (Gaussian Elimination):                                          //
//                                                                            //
// Uses an LU factorization with partial pivoting to overwrite B := A^-1 B    //
//----------------------------------------------------------------------------//

// TODO: Add a serial version

// Parallel version
template<typename T>
void
GaussElim
( DistMatrix<T,MC,MR>& A, DistMatrix<T,MC,MR>& B );

//----------------------------------------------------------------------------//
// Hegst (HErmitian GEneralized to STandard eigenvalue problem):              //
//                                                                            //
// If bothOnLeft,                                                             //
//   reduce the problem A B X = X Lamda to A X = X Lambda                     //
// If ~bothOnLeft,                                                            //
//   reduce the problem A X = B X Lambda to A X = X Lambda                    //
//                                                                            //
// D contains the Cholesky factor of B in the triangle corresponding to the   //
// parameter 'shape'.                                                         //
//----------------------------------------------------------------------------//

// Serial version
template<typename T>
void
Hegst
( bool bothOnLeft, Shape shape, 
  Matrix<T>& A, const Matrix<T>& D );

// Parallel version
template<typename T>
void
Hegst
( bool bothOnLeft, Shape shape, 
  DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B );

//----------------------------------------------------------------------------//
// LU (LU factorization with partial pivoting):                               //
//                                                                            //
// Overwrite A with its LU factorization after partial pivoting: P A = L U.   //
// P is compressed into the vector p by storing the location of the nonzero   //
// element of each row.                                                       //
//----------------------------------------------------------------------------//

// Serial version
template<typename T>
void
LU
( Matrix<T>& A, Matrix<int>& p );

// Parallel version
template<typename T>
void
LU
( DistMatrix<T,MC,MR>& A, DistMatrix<int,VC,Star>& p );

//----------------------------------------------------------------------------//
// Tridiag (Householder tridiagonalization):                                  //
//                                                                            //
// The diagonal and sub/super-diagonal of A are overwritten with a similar    //
// tridiagonal matrix that is found by successively applying Householder      //
// reflections to zero the matrix outside of the tridiagonal band.            //
//                                                                            //
// 'shape' decided which triangle of A specifies the Hermitian matrix, and on //
// exit the transforms are stored above the super/sub-diagonal and are        //
// implicitly one on the super/sub-diagonal.                                  //
//----------------------------------------------------------------------------//

// Serial version for real datatypes
template<typename R>
void
Tridiag
( Shape shape, Matrix<R>& A, Matrix<R>& d, Matrix<R>& e, Matrix<R>& t );

#ifndef WITHOUT_COMPLEX
// Serial version for complex datatypes
template<typename R>
void
Tridiag
( Shape shape, 
  Matrix< std::complex<R> >& A,
  Matrix< R               >& d,
  Matrix< R               >& e,
  Matrix< std::complex<R> >& t );
#endif

// Parallel version for real datatypes
template<typename R>
void
Tridiag
( Shape shape, 
  DistMatrix<R,MC,MR  >& A,
  DistMatrix<R,MD,Star>& d,
  DistMatrix<R,MD,Star>& e,
  DistMatrix<R,MD,Star>& t );

#ifndef WITHOUT_COMPLEX
// Parallel version for complex datatypes
template<typename R>
void
Tridiag
( Shape shape,
  DistMatrix<std::complex<R>,MC,MR  >& A,
  DistMatrix<R,              MD,Star>& d,
  DistMatrix<R,              MD,Star>& e,
  DistMatrix<std::complex<R>,MD,Star>& t );
#endif

//----------------------------------------------------------------------------//
// Trinv (TRiangular INVersion):                                              //
//                                                                            //
// Inverts a triangular matrix. 'shape' determines whether A is assumed to be //
// upper or lower triangular, and 'diagonal' determines whether or not A is   //
// to be treated as having a unit diagonal.                                   //
//----------------------------------------------------------------------------//

// Serial version
template<typename T>
void
Trinv
( Shape shape, Diagonal diagonal, Matrix<T>& A );

// Parallel version
template<typename T>
void
Trinv
( Shape shape, Diagonal diagonal, DistMatrix<T,MC,MR>& A  );

//----------------------------------------------------------------------------//
// UT (UT transform):                                                         //
//                                                                            //
// Applies the accumulated Householder transforms that are stored in the      //
// triangle of H specified by 'shape' to the matrix A.                        //
//                                                                            //
// If 'shape' is set to 'Lower', then offset determines the diagonal that the //
// transforms are stored above (they are implicitly one on that diagonal).    //
// Due to the conventions of the LAPACK routines 'chetrd' and 'zhetrd', the   //
// transforms are assumed to be accumulated left-to-right.                    //
//                                                                            //
// If 'shape' is set to 'Upper', then offset determines the diagonal that the //
// transforms are stored below (they are implicitly one on that diagonal).    //
// Due to the conventions of the LAPACK routines 'chetrd' and 'zhetrd', the   //
// transforms are assumed to be accumulated right-to-left.                    //
//----------------------------------------------------------------------------//

// TODO: Add serial versions

// Parallel version
template<typename T>
void
UT
( Shape shape, 
  Orientation orientation,
  int offset,
  const DistMatrix<T,MC,MR>& H,
        DistMatrix<T,MC,MR>& A );

} // lapack
} // elemental

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

template<typename T>
inline void
elemental::lapack::Chol
( Shape shape, Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("lapack::Chol");
    if( A.Height() != A.Width() )
        throw std::logic_error( "A must be square." );
#endif
    const char uplo = ShapeToChar( shape );
    wrappers::lapack::Chol( uplo, A.Height(), A.Buffer(), A.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::lapack::Hegst
( bool bothOnLeft, Shape shape, Matrix<T>& A, const Matrix<T>& B )
{
#ifndef RELEASE
    PushCallStack("lapack::Hegst");
    if( A.Height() != A.Width() )
        throw std::logic_error( "A must be square." );
    if( B.Height() != B.Width() )
        throw std::logic_error( "B must be square." );
    if( A.Height() != B.Height() )
        throw std::logic_error( "A and B must be the same size." );
#endif
    const int itype = ( bothOnLeft ? 2 : 1 );
    const char uplo = ShapeToChar( shape );
    wrappers::lapack::Hegst
    ( itype, uplo, A.Height(), 
      A.Buffer(), A.LDim(), B.LockedBuffer(), B.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::lapack::LU
( Matrix<T>& A, Matrix<int>& p )
{
#ifndef RELEASE
    PushCallStack("lapack::LU");
    if( p.Height() != A.Height() )
        throw std::logic_error( "A and p must be the same height." );
#endif
    wrappers::lapack::LU
    ( A.Height(), A.Width(), A.Buffer(), A.LDim(), p.Buffer() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline void
elemental::lapack::Tridiag
( Shape shape, Matrix<R>& A, Matrix<R>& d, Matrix<R>& e, Matrix<R>& t )
{
#ifndef RELEASE
    PushCallStack("lapack::Tridiag");
    if( A.Height() != A.Width() )
        throw std::logic_error( "A must be square." );
    if( d.Height() != A.Height() || d.Width() != 1 )
        throw std::logic_error( "d must be a column vector of length n." );
    if( e.Height() != A.Height()-1 || e.Width() != 1 )
        throw std::logic_error( "e must be a column vector of length n-1." );
    if( t.Height() != A.Height()-1 || t.Width() != 1 )
        throw std::logic_error( "t must be a column vector of length n-1." );
#endif
    const char uplo = ShapeToChar( shape );
    wrappers::lapack::Tridiag
    ( uplo, A.Height(), A.Buffer(), A.LDim(),
      d.Buffer(), e.Buffer(), t.Buffer() );
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
template<typename R>
inline void
elemental::lapack::Tridiag
( Shape shape, 
  Matrix< std::complex<R> >& A, 
  Matrix< R               >& d, 
  Matrix< R               >& e, 
  Matrix< std::complex<R> >& t )
{
#ifndef RELEASE
    PushCallStack("lapack::Tridiag");
    if( A.Height() != A.Width() )
        throw std::logic_error( "A must be square." );
    if( d.Height() != A.Height() || d.Width() != 1 )
        throw std::logic_error( "d must be a column vector of length n." );
    if( e.Height() != A.Height()-1 || e.Width() != 1 )
        throw std::logic_error( "e must be a column vector of length n-1." );
    if( t.Height() != A.Height()-1 || t.Width() != 1 )
        throw std::logic_error( "t must be a column vector of length n-1." );
#endif
    const char uplo = ShapeToChar( shape );
    wrappers::lapack::Tridiag
    ( uplo, A.Height(), A.Buffer(), A.LDim(),
      d.Buffer(), e.Buffer(), t.Buffer() );
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX

template<typename T>
inline void
elemental::lapack::Trinv
( Shape shape, Diagonal diagonal, Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("lapack::Trinv");
    if( A.Height() != A.Width() )
        throw std::logic_error( "A must be square." );
#endif
    const char uplo = ShapeToChar( shape );
    const char diag = DiagonalToChar( diagonal );
    wrappers::lapack::Trinv
    ( uplo, diag, A.Height(), A.Buffer(), A.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

#endif /* ELEMENTAL_LAPACK_HPP */

