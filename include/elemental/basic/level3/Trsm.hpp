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

#include "./Trsm/TrsmLLN.hpp"
#include "./Trsm/TrsmLLT.hpp"
#include "./Trsm/TrsmLUN.hpp"
#include "./Trsm/TrsmLUT.hpp"
#include "./Trsm/TrsmRLN.hpp"
#include "./Trsm/TrsmRLT.hpp"
#include "./Trsm/TrsmRUN.hpp"
#include "./Trsm/TrsmRUT.hpp"

namespace elemental {

template<typename F>
inline void
Trsm
( Side side, 
  UpperOrLower uplo, 
  Orientation orientation, 
  Diagonal diagonal,
  F alpha, 
  const DistMatrix<F,MC,MR>& A,
        DistMatrix<F,MC,MR>& X,
  bool checkIfSingular )
{
#ifndef RELEASE
    PushCallStack("Trsm");
#endif
    const int p = X.Grid().Size();
    if( side == LEFT && uplo == LOWER )
    {
        if( orientation == NORMAL )
        {
            if( X.Width() > 5*p )
                internal::TrsmLLNLarge
                ( diagonal, alpha, A, X, checkIfSingular );
            else
                internal::TrsmLLNMedium
                ( diagonal, alpha, A, X, checkIfSingular );
        }
        else
        {
            if( X.Width() > 5*p )
                internal::TrsmLLTLarge
                ( orientation, diagonal, alpha, A, X, checkIfSingular );
            else
                internal::TrsmLLTMedium
                ( orientation, diagonal, alpha, A, X, checkIfSingular );
        }
    }
    else if( side == LEFT && uplo == UPPER )
    {
        if( orientation == NORMAL )
        {
            if( X.Width() > 5*p )
                internal::TrsmLUNLarge
                ( diagonal, alpha, A, X, checkIfSingular );
            else
                internal::TrsmLUNMedium
                ( diagonal, alpha, A, X, checkIfSingular );
        }
        else
        {
            if( X.Width() > 5*p )
                internal::TrsmLUTLarge
                ( orientation, diagonal, alpha, A, X, checkIfSingular );
            else
                internal::TrsmLUTMedium
                ( orientation, diagonal, alpha, A, X, checkIfSingular );
        }
    }
    else if( side == RIGHT && uplo == LOWER )
    {
        if( orientation == NORMAL )
            internal::TrsmRLN( diagonal, alpha, A, X, checkIfSingular );
        else
            internal::TrsmRLT
            ( orientation, diagonal, alpha, A, X, checkIfSingular );
    }
    else if( side == RIGHT && uplo == UPPER )
    {
        if( orientation == NORMAL )
            internal::TrsmRUN( diagonal, alpha, A, X, checkIfSingular );
        else
            internal::TrsmRUT
            ( orientation, diagonal, alpha, A, X, checkIfSingular );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elemental
