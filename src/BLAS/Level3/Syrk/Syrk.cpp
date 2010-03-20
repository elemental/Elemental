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
#include "ElementalBLAS_Internal.h"
using namespace std;
using namespace Elemental;

template<typename T>
void
Elemental::BLAS::Syrk
( const Shape shape, 
  const Orientation orientation,
  const T alpha, const DistMatrix<T,MC,MR>& A,
  const T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("BLAS::Syrk");
    if( orientation == ConjugateTranspose )
    {
        if( A.GetGrid().VCRank() == 0 )
            cerr << "Syrk accepts Normal and Transpose options." << endl;
        DumpCallStack();
        throw exception();
    }
#endif
    if( shape == Lower && orientation == Normal )
    {
        BLAS::Internal::SyrkLN( alpha, A, beta, C );
    }
    else if( shape == Lower )
    {
        BLAS::Internal::SyrkLT( alpha, A, beta, C );
    }
    else if( shape == Upper && orientation == Normal )
    {
        BLAS::Internal::SyrkUN( alpha, A, beta, C );
    }
    else
    {
        BLAS::Internal::SyrkUT( alpha, A, beta, C );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void Elemental::BLAS::Syrk
( const Shape shape, const Orientation orientation,
  const float alpha, const DistMatrix<float,MC,MR>& A,
  const float beta,        DistMatrix<float,MC,MR>& C );

template void Elemental::BLAS::Syrk
( const Shape shape, const Orientation orientation,
  const double alpha, const DistMatrix<double,MC,MR>& A,
  const double beta,        DistMatrix<double,MC,MR>& C );

#ifndef WITHOUT_COMPLEX
template void Elemental::BLAS::Syrk
( const Shape shape, const Orientation orientation,
  const scomplex alpha, const DistMatrix<scomplex,MC,MR>& A,
  const scomplex beta,        DistMatrix<scomplex,MC,MR>& C );

template void Elemental::BLAS::Syrk
( const Shape shape, const Orientation orientation,
  const dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& A,
  const dcomplex beta,        DistMatrix<dcomplex,MC,MR>& C );
#endif

