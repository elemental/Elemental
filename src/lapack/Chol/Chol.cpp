/*
   Copyright (C) 2009-2010 Jack Poulson <jack.poulson@gmail.com>

   This file is part of Elemental.

   Elemental is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   Elemental is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with Elemental.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "elemental/lapack_internal.hpp"
using namespace elemental;

template<typename T>
void
elemental::lapack::Chol
( Shape shape, DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("lapack::Chol");
#endif
    lapack::internal::CholVar3( shape, A );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::lapack::internal::CholVar2
( Shape shape, DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::CholVar2");
#endif
    if( shape == Lower )
        lapack::internal::CholLVar2( A );
    else
        lapack::internal::CholUVar2( A );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::lapack::internal::CholVar2Naive
( Shape shape, DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::CholVar2Naive");
#endif
    if( shape == Lower )
        lapack::internal::CholLVar2Naive( A );
    else
        lapack::internal::CholUVar2Naive( A );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::lapack::internal::CholVar3
( Shape shape, DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::CholVar3");
#endif
    if( shape == Lower )
        lapack::internal::CholLVar3( A );
    else
        lapack::internal::CholUVar3( A );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::lapack::internal::CholVar3Naive
( Shape shape, DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::CholVar3Naive");
#endif
    if( shape == Lower )
        lapack::internal::CholLVar3Naive( A );
    else
        lapack::internal::CholUVar3Naive( A );
#ifndef RELEASE
    PopCallStack();
#endif
}

template void elemental::lapack::Chol
( Shape shape, DistMatrix<float,MC,MR>& A );

template void elemental::lapack::internal::CholVar2
( Shape shape, DistMatrix<float,MC,MR>& A );

template void elemental::lapack::internal::CholVar2Naive
( Shape shape, DistMatrix<float,MC,MR>& A );

template void elemental::lapack::internal::CholVar3
( Shape shape, DistMatrix<float,MC,MR>& A );

template void elemental::lapack::internal::CholVar3Naive
( Shape shape, DistMatrix<float,MC,MR>& A );

template void elemental::lapack::Chol
( Shape shape, DistMatrix<double,MC,MR>& A );

template void elemental::lapack::internal::CholVar2
( Shape shape, DistMatrix<double,MC,MR>& A );

template void elemental::lapack::internal::CholVar2Naive
( Shape shape, DistMatrix<double,MC,MR>& A );

template void elemental::lapack::internal::CholVar3
( Shape shape, DistMatrix<double,MC,MR>& A );

template void elemental::lapack::internal::CholVar3Naive
( Shape shape, DistMatrix<double,MC,MR>& A );

#ifndef WITHOUT_COMPLEX
template void elemental::lapack::Chol
( Shape shape, DistMatrix<scomplex,MC,MR>& A );

template void elemental::lapack::internal::CholVar2
( Shape shape, DistMatrix<scomplex,MC,MR>& A );

template void elemental::lapack::internal::CholVar2Naive
( Shape shape, DistMatrix<scomplex,MC,MR>& A );

template void elemental::lapack::internal::CholVar3
( Shape shape, DistMatrix<scomplex,MC,MR>& A );

template void elemental::lapack::internal::CholVar3Naive
( Shape shape, DistMatrix<scomplex,MC,MR>& A );

template void elemental::lapack::Chol
( Shape shape, DistMatrix<dcomplex,MC,MR>& A );

template void elemental::lapack::internal::CholVar2
( Shape shape, DistMatrix<dcomplex,MC,MR>& A );

template void elemental::lapack::internal::CholVar2Naive
( Shape shape, DistMatrix<dcomplex,MC,MR>& A );

template void elemental::lapack::internal::CholVar3
( Shape shape, DistMatrix<dcomplex,MC,MR>& A );

template void elemental::lapack::internal::CholVar3Naive
( Shape shape, DistMatrix<dcomplex,MC,MR>& A );
#endif

