/*
   This file is part of elemental, a library for distributed-memory dense 
   linear algebra.

   Copyright (C) 2009-2010 Jack Poulson <jack.poulson@gmail.com>

   This program is released under the terms of the license contained in the 
   file LICENSE.
*/
#include "elemental/lapack_internal.hpp"
using namespace elemental;

// The mainline Cholesky wraps the variant 2 algorithm
template<typename T>
void
elemental::lapack::Chol
( Shape shape, DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("lapack::Chol");
#endif
    lapack::internal::CholVar2( shape, A );
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

template void elemental::lapack::Chol
( Shape shape, DistMatrix<float,MC,MR>& A );

template void elemental::lapack::internal::CholVar2
( Shape shape, DistMatrix<float,MC,MR>& A );

template void elemental::lapack::internal::CholVar3
( Shape shape, DistMatrix<float,MC,MR>& A );

template void elemental::lapack::Chol
( Shape shape, DistMatrix<double,MC,MR>& A );

template void elemental::lapack::internal::CholVar2
( Shape shape, DistMatrix<double,MC,MR>& A );

template void elemental::lapack::internal::CholVar3
( Shape shape, DistMatrix<double,MC,MR>& A );

#ifndef WITHOUT_COMPLEX
template void elemental::lapack::Chol
( Shape shape, DistMatrix<scomplex,MC,MR>& A );

template void elemental::lapack::internal::CholVar2
( Shape shape, DistMatrix<scomplex,MC,MR>& A );

template void elemental::lapack::internal::CholVar3
( Shape shape, DistMatrix<scomplex,MC,MR>& A );

template void elemental::lapack::Chol
( Shape shape, DistMatrix<dcomplex,MC,MR>& A );

template void elemental::lapack::internal::CholVar2
( Shape shape, DistMatrix<dcomplex,MC,MR>& A );

template void elemental::lapack::internal::CholVar3
( Shape shape, DistMatrix<dcomplex,MC,MR>& A );
#endif

