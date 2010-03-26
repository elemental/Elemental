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
#include "ElementalBLASInternal.h"
using namespace Elemental;

template<typename T>
void
Elemental::BLAS::Trsm
( const Side side, 
  const Shape shape, 
  const Orientation orientation, 
  const Diagonal diagonal,
  const T alpha, 
  const DistMatrix<T,MC,MR>& A,
        DistMatrix<T,MC,MR>& X   )
{
#ifndef RELEASE
    PushCallStack("BLAS::Trsm");
#endif
    if( side == Left && shape == Lower )
    {
        if( orientation == Normal )
            BLAS::Internal::TrsmLLN( diagonal, alpha, A, X );
        else
            BLAS::Internal::TrsmLLT( orientation, diagonal, alpha, A, X );
    }
    else if( side == Left && shape == Upper )
    {
        if( orientation == Normal )
            BLAS::Internal::TrsmLUN( diagonal, alpha, A, X );
        else
            BLAS::Internal::TrsmLUT( orientation, diagonal, alpha, A, X );
    }
    else if( side == Right && shape == Lower )
    {
        if( orientation == Normal )
            BLAS::Internal::TrsmRLN( diagonal, alpha, A, X );
        else
            BLAS::Internal::TrsmRLT( orientation, diagonal, alpha, A, X );
    }
    else if( side == Right && shape == Upper )
    {
        if( orientation == Normal )
            BLAS::Internal::TrsmRUN( diagonal, alpha, A, X );
        else
            BLAS::Internal::TrsmRUT( orientation, diagonal, alpha, A, X );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void Elemental::BLAS::Trsm
( const Side side, 
  const Shape shape, 
  const Orientation orientation, 
  const Diagonal diagonal,
  const float alpha, 
  const DistMatrix<float,MC,MR>& A,
        DistMatrix<float,MC,MR>& X );

template void Elemental::BLAS::Trsm
( const Side side, 
  const Shape shape, 
  const Orientation orientation, 
  const Diagonal diagonal,
  const double alpha, 
  const DistMatrix<double,MC,MR>& A,
        DistMatrix<double,MC,MR>& X );

#ifndef WITHOUT_COMPLEX
template void Elemental::BLAS::Trsm
( const Side side, 
  const Shape shape, 
  const Orientation orientation, 
  const Diagonal diagonal,
  const scomplex alpha, 
  const DistMatrix<scomplex,MC,MR>& A,
        DistMatrix<scomplex,MC,MR>& X );

template void Elemental::BLAS::Trsm
( const Side side, 
  const Shape shape, 
  const Orientation orientation, 
  const Diagonal diagonal,
  const dcomplex alpha, 
  const DistMatrix<dcomplex,MC,MR>& A,
        DistMatrix<dcomplex,MC,MR>& X );
#endif

