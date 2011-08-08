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
#include "elemental/advanced_internal.hpp"
using namespace elemental;
using namespace elemental::imports;

template<typename R> // representation of a real number
R
advanced::HermitianNorm
( Shape shape, const Matrix<R>& A, NormType type )
{
#ifndef RELEASE
    PushCallStack("advanced::HermitianNorm");
#endif
    R norm = 0;
    switch( type )
    {
    case InfinityNorm:
        norm = advanced::internal::HermitianInfinityNorm( shape, A );
        break;
    case FrobeniusNorm: 
        norm = advanced::internal::HermitianFrobeniusNorm( shape, A );
        break;
    case MaxNorm:
        norm = advanced::internal::HermitianMaxNorm( shape, A );
        break;
    case OneNorm:
        norm = advanced::internal::HermitianOneNorm( shape, A );
        break;
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

#ifndef WITHOUT_COMPLEX
template<typename R> // representation of a real number
R
advanced::HermitianNorm
( Shape shape, const Matrix<std::complex<R> >& A, NormType type )
{
#ifndef RELEASE
    PushCallStack("advanced::HermitianNorm");
#endif
    R norm = 0;
    switch( type )
    {
    case InfinityNorm:
        norm = advanced::internal::HermitianInfinityNorm( shape, A );
        break;
    case FrobeniusNorm: 
        norm = advanced::internal::HermitianFrobeniusNorm( shape, A );
        break;
    case MaxNorm:
        norm = advanced::internal::HermitianMaxNorm( shape, A );
        break;
    case OneNorm:
        norm = advanced::internal::HermitianOneNorm( shape, A );
        break;
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}
#endif // WITHOUT_COMPLEX

template<typename R> // representation of a real number
R
advanced::HermitianNorm
( Shape shape, const DistMatrix<R,MC,MR>& A, NormType type )
{
#ifndef RELEASE
    PushCallStack("advanced::HermitianNorm");
#endif
    R norm = 0;
    switch( type )
    {
    case InfinityNorm:
        norm = advanced::internal::HermitianInfinityNorm( shape, A );
        break;
    case FrobeniusNorm: 
        norm = advanced::internal::HermitianFrobeniusNorm( shape, A );
        break;
    case MaxNorm:
        norm = advanced::internal::HermitianMaxNorm( shape, A );
        break;
    case OneNorm:
        norm = advanced::internal::HermitianOneNorm( shape, A );
        break;
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

#ifndef WITHOUT_COMPLEX
template<typename R> // representation of a real number
R
advanced::HermitianNorm
( Shape shape, const DistMatrix<std::complex<R>,MC,MR>& A, NormType type )
{
#ifndef RELEASE
    PushCallStack("advanced::HermitianNorm");
#endif
    R norm = 0;
    switch( type )
    {
    case InfinityNorm:
        norm = advanced::internal::HermitianInfinityNorm( shape, A );
        break;
    case FrobeniusNorm: 
        norm = advanced::internal::HermitianFrobeniusNorm( shape, A );
        break;
    case MaxNorm:
        norm = advanced::internal::HermitianMaxNorm( shape, A );
        break;
    case OneNorm:
        norm = advanced::internal::HermitianOneNorm( shape, A );
        break;
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}
#endif // WITHOUT_COMPLEX

template float elemental::advanced::HermitianNorm
( Shape shape, const Matrix<float>& A, NormType type );
template double elemental::advanced::HermitianNorm
( Shape shape, const Matrix<double>& A, NormType type );
#ifndef WITHOUT_COMPLEX
template float elemental::advanced::HermitianNorm
( Shape shape, const Matrix<std::complex<float> >& A, NormType type );
template double elemental::advanced::HermitianNorm
( Shape shape, const Matrix<std::complex<double> >& A, NormType type );
#endif

template float elemental::advanced::HermitianNorm
( Shape shape, const DistMatrix<float,MC,MR>& A, NormType type );
template double elemental::advanced::HermitianNorm
( Shape shape, const DistMatrix<double,MC,MR>& A, NormType type );
#ifndef WITHOUT_COMPLEX
template float elemental::advanced::HermitianNorm
( Shape shape, const DistMatrix<std::complex<float>,MC,MR>& A, NormType type );
template double elemental::advanced::HermitianNorm
( Shape shape, const DistMatrix<std::complex<double>,MC,MR>& A, NormType type );
#endif
