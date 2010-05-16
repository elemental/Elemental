/*
   Copyright 2009-2010 Jack Poulson

   This file is part of Elemental.

   Elemental is free software: you can redistribute it and/or modify it under
   the terms of the GNU Lesser General Public License as published by the
   Free Software Foundation; either version 3 of the License, or 
   (at your option) any later version.

   Elemental is distributed in the hope that it will be useful, but 
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with Elemental. If not, see <http://www.gnu.org/licenses/>.
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

template<typename T>
void
Tridiag
( Shape shape, Matrix<T>& A, Matrix<T>& d, Matrix<T>& e, Matrix<T>& tau );

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

template<typename T>
void
Tridiag
( Shape shape, 
  DistMatrix<T,MC,MR  >& A,
  DistMatrix<T,MD,Star>& d,
  DistMatrix<T,MD,Star>& e,
  DistMatrix<T,MD,Star>& t );

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
        throw "A must be square.";
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
        throw "A and p must be the same height.";
#endif
    wrappers::lapack::LU
    ( A.Height(), A.Width(), A.Buffer(), A.LDim(), p.Buffer() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::lapack::Tridiag
( Shape shape, Matrix<T>& A, Matrix<T>& d, Matrix<T>& e, Matrix<T>& t )
{
#ifndef RELEASE
    PushCallStack("lapack::Tridiag");
    if( A.Height() != A.Width() )
        throw "A must be square.";
    if( d.Height() != A.Height() || d.Width() != 1 )
        throw "d must be a column vector of length n.";
    if( e.Height() != A.Height()-1 || e.Width() != 1 )
        throw "e must be a column vector of length n-1.";
    if( t.Height() != A.Height()-1 || t.Width() != 1 )
        throw "t must be a column vector of length n-1.";
#endif
    const char uplo = ShapeToChar( shape );
    wrappers::lapack::Tridiag
    ( uplo, A.Height(), A.Buffer(), A.LDim(),
      d.Buffer(), e.Buffer(), t.Buffer() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::lapack::Trinv
( Shape shape, Diagonal diagonal, Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("lapack::Trinv");
    if( A.Height() != A.Width() )
        throw "A must be square.";
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

