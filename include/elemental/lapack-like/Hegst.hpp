/*
   Copyright (c) 2009-2012, Jack Poulson
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

#include "./Hegst/LLVar1.hpp"
#include "./Hegst/LLVar2.hpp"
#include "./Hegst/LLVar4.hpp"
#include "./Hegst/LLVar5.hpp"
#include "./Hegst/LUVar1.hpp"
#include "./Hegst/LUVar2.hpp"
#include "./Hegst/LUVar4.hpp"
#include "./Hegst/LUVar5.hpp"
#include "./Hegst/RLVar1.hpp"
#include "./Hegst/RLVar2.hpp"
#include "./Hegst/RLVar3.hpp"
#include "./Hegst/RLVar4.hpp"
#include "./Hegst/RLVar5.hpp"
#include "./Hegst/RUVar1.hpp"
#include "./Hegst/RUVar2.hpp"
#include "./Hegst/RUVar3.hpp"
#include "./Hegst/RUVar4.hpp"
#include "./Hegst/RUVar5.hpp"

namespace elem {

template<typename F> 
inline void
Hegst
( LeftOrRight side, UpperOrLower uplo, 
  DistMatrix<F>& A, const DistMatrix<F>& B )
{
#ifndef RELEASE
    PushCallStack("Hegst");
#endif
    if( side == LEFT )
    {
        if( uplo == LOWER )
            internal::HegstLLVar4( A, B );
        else
            internal::HegstLUVar4( A, B );
    }
    else
    {
        if( uplo == LOWER )
            internal::HegstRLVar4( A, B );
        else
            internal::HegstRUVar4( A, B );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
