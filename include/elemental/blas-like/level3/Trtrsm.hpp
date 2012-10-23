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

#include "./Trtrsm/LLN.hpp"

namespace elem {

template<typename F>
inline void
Trtrsm
( LeftOrRight side, UpperOrLower uplo,
  Orientation orientation, UnitOrNonUnit diag,
  F alpha, const Matrix<F>& A, Matrix<F>& X,
  bool checkIfSingular )
{
#ifndef RELEASE
    PushCallStack("Trtrsm");
    if( A.Height() != A.Width() || X.Height() != X.Width() )
        throw std::logic_error("Triangular matrices must be square");
    if( A.Height() != X.Height() )
        throw std::logic_error("Nonconformal Trtrsm");
#endif
    if( side == LEFT && uplo == LOWER )
    {
        if( orientation == NORMAL )
            internal::TrtrsmLLN( diag, alpha, A, X, checkIfSingular );
        else
            throw std::logic_error("This option not yet implemented");
        /*
            internal::TrtrsmLLT
            ( orientation, diag, alpha, A, X, checkIfSingular );
        */
    }
    else
        throw std::logic_error("This option not yet implemented");
    /*
    else if( side == LEFT && uplo == UPPER )
    {
        if( orientation == NORMAL )
            internal::TrtrsmLUN( diag, alpha, A, X, checkIfSingular );
        else
            internal::TrtrsmLUT
            ( orientation, diag, alpha, A, X, checkIfSingular );
    }
    else if( side == RIGHT && uplo == LOWER )
    {
        if( orientation == NORMAL )
            internal::TrtrsmRLN( diag, alpha, A, X, checkIfSingular );
        else
            internal::TrtrsmRLT
            ( orientation, diag, alpha, A, X, checkIfSingular );
    }
    else if( side == RIGHT && uplo == UPPER )
    {
        if( orientation == NORMAL )
            internal::TrtrsmRUN( diag, alpha, A, X, checkIfSingular );
        else
            internal::TrtrsmRUT
            ( orientation, diag, alpha, A, X, checkIfSingular );
    }
    */
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
Trtrsm
( LeftOrRight side, UpperOrLower uplo, 
  Orientation orientation, UnitOrNonUnit diag,
  F alpha, const DistMatrix<F>& A, DistMatrix<F>& X,
  bool checkIfSingular )
{
#ifndef RELEASE
    PushCallStack("Trtrsm");
#endif
    if( side == LEFT && uplo == LOWER )
    {
        if( orientation == NORMAL )
            internal::TrtrsmLLN( diag, alpha, A, X, checkIfSingular );
        else
            throw std::logic_error("This option not yet implemented");
        /*
            internal::TrtrsmLLT
            ( orientation, diag, alpha, A, X, checkIfSingular );
        */
    }
    else
        throw std::logic_error("This option not yet implemented");
    /*
    else if( side == LEFT && uplo == UPPER )
    {
        if( orientation == NORMAL )
            internal::TrtrsmLUN( diag, alpha, A, X, checkIfSingular );
        else
            internal::TrtrsmLUT
            ( orientation, diag, alpha, A, X, checkIfSingular );
    }
    else if( side == RIGHT && uplo == LOWER )
    {
        if( orientation == NORMAL )
            internal::TrtrsmRLN( diag, alpha, A, X, checkIfSingular );
        else
            internal::TrtrsmRLT
            ( orientation, diag, alpha, A, X, checkIfSingular );
    }
    else if( side == RIGHT && uplo == UPPER )
    {
        if( orientation == NORMAL )
            internal::TrtrsmRUN( diag, alpha, A, X, checkIfSingular );
        else
            internal::TrtrsmRUT
            ( orientation, diag, alpha, A, X, checkIfSingular );
    }
    */
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
