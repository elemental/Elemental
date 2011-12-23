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


#include "./Bidiag/LocalBidiag.hpp"
#include "./Bidiag/PanelBidiagL.hpp"
#include "./Bidiag/PanelBidiagU.hpp"
#include "./Bidiag/BidiagL.hpp"
#include "./Bidiag/BidiagU.hpp"

template<typename R> 
inline void elemental::advanced::Bidiag( DistMatrix<R,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("advanced::Bidiag");
#endif
    if( IsComplex<R>::val )
        throw std::logic_error("Called real routine with complex datatype");
    if( A.Height() >= A.Width() )
        advanced::internal::BidiagU( A );
    else
        advanced::internal::BidiagL( A );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R> 
inline void elemental::advanced::Bidiag
( DistMatrix<std::complex<R>,MC,  MR  >& A,
  DistMatrix<std::complex<R>,STAR,STAR>& tP,
  DistMatrix<std::complex<R>,STAR,STAR>& tQ )
{
#ifndef RELEASE
    PushCallStack("advanced::Bidiag");
#endif
    if( A.Height() >= A.Width() )
        advanced::internal::BidiagU( A, tP, tQ );
    else
        advanced::internal::BidiagL( A, tP, tQ );
#ifndef RELEASE
    PopCallStack();
#endif
}

