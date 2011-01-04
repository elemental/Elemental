/*
   Copyright (c) 2009-2011, Jack Poulson
   All rights reserved.

   This file is part of Elemental.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

    - Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

    - Neither the name of the owner nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/
#include "elemental/lapack_internal.hpp"
using namespace elemental;

template<typename F> // F represents a real or complex field
void
elemental::lapack::Chol
( Shape shape, DistMatrix<F,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("lapack::Chol");
#endif
    lapack::internal::CholVar3( shape, A );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F> // F represents a real or complex field
void
elemental::lapack::internal::CholVar2
( Shape shape, DistMatrix<F,MC,MR>& A )
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

template<typename F> // F represents a real or complex field
void
elemental::lapack::internal::CholVar2Naive
( Shape shape, DistMatrix<F,MC,MR>& A )
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

template<typename F> // F represents a real or complex field
void
elemental::lapack::internal::CholVar3
( Shape shape, DistMatrix<F,MC,MR>& A )
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

template<typename F> // F represents a real or complex field
void
elemental::lapack::internal::CholVar3Naive
( Shape shape, DistMatrix<F,MC,MR>& A )
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

