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
using namespace std;
using namespace elemental;

template<typename F> // F represents a real or complex field
void
elemental::lapack::Hegst
( Side side, Shape shape, 
  DistMatrix<F,MC,MR>& A, const DistMatrix<F,MC,MR>& B )
{
#ifndef RELEASE
    PushCallStack("lapack::Hegst");
#endif
    if( side==Left )
    {
        if( shape == Lower )
            lapack::internal::HegstLLVar4( A, B );
        else
            lapack::internal::HegstLUVar4( A, B );
    }
    else
    {
        if( shape == Lower )
            lapack::internal::HegstRLVar4( A, B );
        else
            lapack::internal::HegstRUVar4( A, B );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void elemental::lapack::Hegst
( Side side, Shape shape, 
  DistMatrix<float,MC,MR>& A, const DistMatrix<float,MC,MR>& B );
template void elemental::lapack::Hegst
( Side side, Shape shape,
  DistMatrix<double,MC,MR>& A, const DistMatrix<double,MC,MR>& B );

#ifndef WITHOUT_COMPLEX
template void elemental::lapack::Hegst
( Side side, Shape shape,
  DistMatrix<scomplex,MC,MR>& A, const DistMatrix<scomplex,MC,MR>& B );
template void elemental::lapack::Hegst
( Side side, Shape shape,
  DistMatrix<dcomplex,MC,MR>& A, const DistMatrix<dcomplex,MC,MR>& B );
#endif

