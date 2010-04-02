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
#include "Elemental/BLASInternal.hpp"
using namespace Elemental;

template<typename T>
void
Elemental::BLAS::Hemm
( const Side side, const Shape shape,
  const T alpha, const DistMatrix<T,MC,MR>& A,
                 const DistMatrix<T,MC,MR>& B,
  const T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("BLAS::Hemm");
#endif
    if( side == Left && shape == Lower )
    {
        BLAS::Internal::HemmLL( alpha, A, B, beta, C );
    }
    else if( side == Left )
    {
        BLAS::Internal::HemmLU( alpha, A, B, beta, C );
    }
    else if( shape == Lower )
    {
        BLAS::Internal::HemmRL( alpha, A, B, beta, C );
    }
    else
    {
        BLAS::Internal::HemmRU( alpha, A, B, beta, C );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void Elemental::BLAS::Hemm
( const Side side, const Shape shape,
  const float alpha, const DistMatrix<float,MC,MR>& A,
                     const DistMatrix<float,MC,MR>& B,
  const float beta,        DistMatrix<float,MC,MR>& C );

template void Elemental::BLAS::Hemm
( const Side side, const Shape shape,
  const double alpha, const DistMatrix<double,MC,MR>& A,
                      const DistMatrix<double,MC,MR>& B,
  const double beta,        DistMatrix<double,MC,MR>& C );

#ifndef WITHOUT_COMPLEX
template void Elemental::BLAS::Hemm
( const Side side, const Shape shape,
  const scomplex alpha, const DistMatrix<scomplex,MC,MR>& A,
                        const DistMatrix<scomplex,MC,MR>& B,
  const scomplex beta,        DistMatrix<scomplex,MC,MR>& C );

template void Elemental::BLAS::Hemm
( const Side side, const Shape shape,
  const dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& A,
                        const DistMatrix<dcomplex,MC,MR>& B,
  const dcomplex beta,        DistMatrix<dcomplex,MC,MR>& C );
#endif
