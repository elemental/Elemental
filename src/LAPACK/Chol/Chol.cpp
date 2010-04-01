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
#include "ElementalLAPACKInternal.h"
using namespace Elemental;

// The mainline Cholesky wraps the variant 2 algorithm
template<typename T>
void
Elemental::LAPACK::Chol
( const Shape shape, DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("LAPACK::Chol");
#endif
    LAPACK::Internal::CholVar2( shape, A );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::LAPACK::Internal::CholVar2
( const Shape shape, DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("LAPACK::Internal::CholVar2");
#endif
    if( shape == Lower )
        LAPACK::Internal::CholLVar2( A );
    else
        LAPACK::Internal::CholUVar2( A );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::LAPACK::Internal::CholVar3
( const Shape shape, DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("LAPACK::Internal::CholVar3");
#endif
    if( shape == Lower )
        LAPACK::Internal::CholLVar3( A );
    else
        LAPACK::Internal::CholUVar3( A );
#ifndef RELEASE
    PopCallStack();
#endif
}

template void Elemental::LAPACK::Chol
( const Shape shape, DistMatrix<float,MC,MR>& A );

template void Elemental::LAPACK::Internal::CholVar2
( const Shape shape, DistMatrix<float,MC,MR>& A );

template void Elemental::LAPACK::Internal::CholVar3
( const Shape shape, DistMatrix<float,MC,MR>& A );

template void Elemental::LAPACK::Chol
( const Shape shape, DistMatrix<double,MC,MR>& A );

template void Elemental::LAPACK::Internal::CholVar2
( const Shape shape, DistMatrix<double,MC,MR>& A );

template void Elemental::LAPACK::Internal::CholVar3
( const Shape shape, DistMatrix<double,MC,MR>& A );

#ifndef WITHOUT_COMPLEX
template void Elemental::LAPACK::Chol
( const Shape shape, DistMatrix<scomplex,MC,MR>& A );

template void Elemental::LAPACK::Internal::CholVar2
( const Shape shape, DistMatrix<scomplex,MC,MR>& A );

template void Elemental::LAPACK::Internal::CholVar3
( const Shape shape, DistMatrix<scomplex,MC,MR>& A );

template void Elemental::LAPACK::Chol
( const Shape shape, DistMatrix<dcomplex,MC,MR>& A );

template void Elemental::LAPACK::Internal::CholVar2
( const Shape shape, DistMatrix<dcomplex,MC,MR>& A );

template void Elemental::LAPACK::Internal::CholVar3
( const Shape shape, DistMatrix<dcomplex,MC,MR>& A );
#endif

