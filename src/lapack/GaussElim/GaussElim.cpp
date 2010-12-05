/*
   Copyright (c) 2009-2010, Jack Poulson
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
using namespace elemental::lapack::internal;

template<typename T>
void
elemental::lapack::GaussElim
( DistMatrix<T,MC,MR>& A, DistMatrix<T,MC,MR>& B )
{
#ifndef RELEASE
    PushCallStack("lapack::GaussElim");
    if( A.Grid() != B.Grid() )
        throw logic_error( "A and B must be distributed over the same grid." );
    if( A.Height() != A.Width() )
        throw logic_error( "A must be square." );
    if( A.Height() != B.Height() )
        throw logic_error( "A and B must be the same height." );
#endif
    lapack::internal::ReduceToRowEchelon( A, B );
    if( B.Width() == 1 )
        blas::Trsv( Upper, Normal, NonUnit, A, B );
    else
        blas::Trsm( Left, Upper, Normal, NonUnit, (T)1, A, B );
#ifndef RELEASE
    PopCallStack();
#endif
}

template void
elemental::lapack::GaussElim
( DistMatrix<float,MC,MR>& A, DistMatrix<float,MC,MR>& B );

template void
elemental::lapack::GaussElim
( DistMatrix<double,MC,MR>& A, DistMatrix<double,MC,MR>& B );

#ifndef WITHOUT_COMPLEX
template void
elemental::lapack::GaussElim
( DistMatrix<scomplex,MC,MR>& A, DistMatrix<scomplex,MC,MR>& B );

template void
elemental::lapack::GaussElim
( DistMatrix<dcomplex,MC,MR>& A, DistMatrix<dcomplex,MC,MR>& B );
#endif

