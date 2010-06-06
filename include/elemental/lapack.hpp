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

//--------------------------------------------------------------------//
// Local LAPACK                                                       //
//--------------------------------------------------------------------//
template<typename T>
void
Chol
( Shape shape, Matrix<T>& A );

template<typename T>
void
LU
( Matrix<T>& A, Matrix<int>& p );

template<typename R>
void
Tridiag
( Shape shape, Matrix<R>& A, Matrix<R>& d, Matrix<R>& e, Matrix<R>& t );

#ifndef WITHOUT_COMPLEX
template<typename R>
void
Tridiag
( Shape shape, 
  Matrix< std::complex<R> >& A,
  Matrix< R               >& d,
  Matrix< R               >& e,
  Matrix< std::complex<R> >& t );
#endif

template<typename T>
void
Trinv
( Shape shape, Diagonal diagonal, Matrix<T>& A );

//--------------------------------------------------------------------//
// Distributed LAPACK                                                 //
//--------------------------------------------------------------------//
template<typename T>
void
Chol
( Shape shape, DistMatrix<T,MC,MR>& A );

template<typename T>
void
GaussElim
( DistMatrix<T,MC,MR>& A, DistMatrix<T,MC,MR>& B );

template<typename T>
void
LU
( DistMatrix<T,MC,MR>& A, DistMatrix<int,VC,Star>& p );

template<typename R>
void
Tridiag
( Shape shape, 
  DistMatrix<R,MC,MR  >& A,
  DistMatrix<R,MD,Star>& d,
  DistMatrix<R,MD,Star>& e,
  DistMatrix<R,MD,Star>& t );

#ifndef WITHOUT_COMPLEX
template<typename R>
void
Tridiag
( Shape shape,
  DistMatrix<std::complex<R>,MC,MR  >& A,
  DistMatrix<R,              MD,Star>& d,
  DistMatrix<R,              MD,Star>& e,
  DistMatrix<std::complex<R>,MD,Star>& t );
#endif

template<typename T>
void
Trinv
( Shape shape, Diagonal diagonal, DistMatrix<T,MC,MR>& A  );

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

