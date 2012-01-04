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

#include "./Trmm/TrmmUtil.hpp"
#include "./Trmm/TrmmLLN.hpp"
#include "./Trmm/TrmmLLT.hpp"
#include "./Trmm/TrmmLUN.hpp"
#include "./Trmm/TrmmLUT.hpp"
#include "./Trmm/TrmmRLN.hpp"
#include "./Trmm/TrmmRLT.hpp"
#include "./Trmm/TrmmRUN.hpp"
#include "./Trmm/TrmmRUT.hpp"

namespace elemental {

template<typename T>
inline void
Trmm
( Side side, UpperOrLower uplo, Orientation orientation, Diagonal diagonal,
  T alpha, const DistMatrix<T,MC,MR>& A, DistMatrix<T,MC,MR>& X )
{
#ifndef RELEASE
    PushCallStack("Trmm");
#endif
    if( side == LEFT && uplo == LOWER )
    {
        if( orientation == NORMAL )
            internal::TrmmLLN( diagonal, alpha, A, X );
        else
            internal::TrmmLLT( orientation, diagonal, alpha, A, X );
    }
    else if( side == LEFT )
    {
        if( orientation == NORMAL )
            internal::TrmmLUN( diagonal, alpha, A, X );
        else
            internal::TrmmLUT( orientation, diagonal, alpha, A, X );
    }
    else if( uplo == LOWER )
    {
        if( orientation == NORMAL )
            internal::TrmmRLN( diagonal, alpha, A, X );
        else
            internal::TrmmRLT( orientation, diagonal, alpha, A, X );
    }
    else
    {
        if( orientation == NORMAL )
            internal::TrmmRUN( diagonal, alpha, A, X );
        else
            internal::TrmmRUT( orientation, diagonal, alpha, A, X );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elemental
