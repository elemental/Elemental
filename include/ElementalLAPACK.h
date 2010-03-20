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
#ifndef ELEMENTAL_LAPACK_H
#define ELEMENTAL_LAPACK_H 1

#include "ElementalBLAS.h"
#include "wrappers/LAPACK.h"

namespace Elemental
{
    namespace LAPACK
    {
        //--------------------------------------------------------------------//
        // Local LAPACK                                                       //
        //--------------------------------------------------------------------//
        template<typename T>
        void
        Chol
        ( const Shape shape, Matrix<T>& A );

        template<typename T>
        void
        LU
        ( Matrix<T>& A, Matrix<int>& p );

        template<typename T>
        void
        Tridiag
        ( const Shape shape, Matrix<T>& A, 
          Matrix<T>& d, Matrix<T>& e, Matrix<T>& tau );

        template<typename T>
        void
        Trinv
        ( const Shape shape, const Diagonal diagonal, Matrix<T>& A );

        //--------------------------------------------------------------------//
        // Distributed LAPACK                                                 //
        //--------------------------------------------------------------------//
        template<typename T>
        void
        Chol
        ( const Shape shape, DistMatrix<T,MC,MR>& A );

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
        ( const Shape shape, 
          DistMatrix<T,MC,MR  >& A,
          DistMatrix<T,MD,Star>& d,
          DistMatrix<T,MD,Star>& e,
          DistMatrix<T,MD,Star>& t );

        template<typename T>
        void
        Trinv
        ( const Shape shape, 
          const Diagonal diagonal,
          DistMatrix<T,MC,MR>& A  );
    }
}

//----------------------------------------------------------------------------//

template<typename T>
inline void
Elemental::LAPACK::Chol
( const Shape shape, Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("LAPACK::Chol");
    if( A.Height() != A.Width() )
    {
        std::cerr << "A must be square." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
#endif
    const char uplo = ShapeToChar( shape );
    Elemental::wrappers::LAPACK::Chol( uplo, A.Height(), A.Buffer(), A.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Elemental::LAPACK::LU
( Matrix<T>& A, Matrix<int>& p )
{
#ifndef RELEASE
    PushCallStack("LAPACK::LU");
    if( p.Height() != A.Height() )
    {
        std::cerr << "A and p must be the same height." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
#endif
    Elemental::wrappers::LAPACK::LU
    ( A.Height(), A.Width(), A.Buffer(), A.LDim(), p.Buffer() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Elemental::LAPACK::Tridiag
( const Shape shape, Matrix<T>& A,
  Matrix<T>& d, Matrix<T>& e, Matrix<T>& t )
{
#ifndef RELEASE
    PushCallStack("LAPACK::Tridiag");
    if( A.Height() != A.Width() )
    {
        std::cerr << "A must be square." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
    if( d.Height() != A.Height() || d.Width() != 1 )
    {
        std::cerr << "d must be a column vector of length n." << std::endl;
        DumpCallStack(); 
        throw std::exception();
    }
    if( e.Height() != A.Height()-1 || e.Width() != 1 )
    {
        std::cerr << "e must be a column vector of length n-1." << std::endl;
        DumpCallStack(); 
        throw std::exception();
    }
    if( t.Height() != A.Height()-1 || t.Width() != 1 )
    {
        std::cerr << "t must be a column vector of length n-1." << std::endl;
        DumpCallStack(); 
        throw std::exception();
    }
#endif
    const char uplo = ShapeToChar( shape );
    Elemental::wrappers::LAPACK::Tridiag
    ( uplo, A.Height(), A.Buffer(), A.LDim(),
      d.Buffer(), e.Buffer(), t.Buffer()     );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Elemental::LAPACK::Trinv
( const Shape shape, const Diagonal diagonal, Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("LAPACK::Trinv");
    if( A.Height() != A.Width() )
    {
        std::cerr << "A must be square." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
#endif
    const char uplo = ShapeToChar( shape );
    const char diag = DiagonalToChar( diagonal );
    Elemental::wrappers::LAPACK::Trinv
    ( uplo, diag, A.Height(), A.Buffer(), A.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

#endif /* ELEMENTAL_LAPACK_H */

