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

#include "./ApplyPackedReflectors/UTUtil.hpp"
#include "./ApplyPackedReflectors/ApplyPackedReflectorsLLHB.hpp"
#include "./ApplyPackedReflectors/ApplyPackedReflectorsLLHF.hpp"
#include "./ApplyPackedReflectors/ApplyPackedReflectorsLLVB.hpp"
#include "./ApplyPackedReflectors/ApplyPackedReflectorsLLVF.hpp"
#include "./ApplyPackedReflectors/ApplyPackedReflectorsLUHB.hpp"
#include "./ApplyPackedReflectors/ApplyPackedReflectorsLUHF.hpp"
#include "./ApplyPackedReflectors/ApplyPackedReflectorsLUVB.hpp"
#include "./ApplyPackedReflectors/ApplyPackedReflectorsLUVF.hpp"
#include "./ApplyPackedReflectors/ApplyPackedReflectorsRLHB.hpp"
#include "./ApplyPackedReflectors/ApplyPackedReflectorsRLHF.hpp"
#include "./ApplyPackedReflectors/ApplyPackedReflectorsRLVB.hpp"
#include "./ApplyPackedReflectors/ApplyPackedReflectorsRLVF.hpp"
#include "./ApplyPackedReflectors/ApplyPackedReflectorsRUHB.hpp"
#include "./ApplyPackedReflectors/ApplyPackedReflectorsRUHF.hpp"
#include "./ApplyPackedReflectors/ApplyPackedReflectorsRUVB.hpp"
#include "./ApplyPackedReflectors/ApplyPackedReflectorsRUVF.hpp"

template<typename R> // representation of a real number
inline void
elemental::advanced::ApplyPackedReflectors
( Side side, UpperOrLower uplo, 
  VectorDirection direction, ForwardOrBackward order,
  int offset,
  const DistMatrix<R,MC,MR>& H, 
        DistMatrix<R,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("advanced::ApplyPackedReflectors");
#endif
    // Since the complex version does not have the same argument list, there is
    // currently no good way to ensure that this version is not called with 
    // complex datatypes. Until C++11 compilers are commonplace, we cannot
    // use static_assert either.
    if( IsComplex<R>::val )
        throw std::logic_error("Called real routine with complex datatype");

    if( side == LEFT )
    {
        if( uplo == LOWER )
        {
            if( direction == VERTICAL && order == FORWARD )
                advanced::internal::ApplyPackedReflectorsLLVF( offset, H, A );
            else if( direction == VERTICAL )
                advanced::internal::ApplyPackedReflectorsLLVB( offset, H, A );
            else if( order == FORWARD )
                advanced::internal::ApplyPackedReflectorsLLHF( offset, H, A );
            else
                advanced::internal::ApplyPackedReflectorsLLHB( offset, H, A );
        }
        else
        {
            if( direction == VERTICAL && order == FORWARD )
                advanced::internal::ApplyPackedReflectorsLUVF( offset, H, A );
            else if( direction == VERTICAL )
                advanced::internal::ApplyPackedReflectorsLUVB( offset, H, A );
            else if( order == FORWARD )
                advanced::internal::ApplyPackedReflectorsLUHF( offset, H, A );
            else
                advanced::internal::ApplyPackedReflectorsLUHB( offset, H, A );
        }
    }
    else
    {
        if( uplo == LOWER )
        {
            if( direction == VERTICAL && order == FORWARD )
                advanced::internal::ApplyPackedReflectorsRLVF( offset, H, A );
            else if( direction == VERTICAL )
                advanced::internal::ApplyPackedReflectorsRLVB( offset, H, A );
            else if( order == FORWARD )
                advanced::internal::ApplyPackedReflectorsRLHF( offset, H, A );
            else
                advanced::internal::ApplyPackedReflectorsRLHB( offset, H, A );
        }
        else
        {
            if( direction == VERTICAL && order == FORWARD )
                advanced::internal::ApplyPackedReflectorsRUVF( offset, H, A );
            else if( direction == VERTICAL )
                advanced::internal::ApplyPackedReflectorsRUVB( offset, H, A );
            else if( order == FORWARD )
                advanced::internal::ApplyPackedReflectorsRUHF( offset, H, A );
            else
                advanced::internal::ApplyPackedReflectorsRUHB( offset, H, A );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R> // representation of a real number
inline void
elemental::advanced::ApplyPackedReflectors
( Side side, UpperOrLower uplo, 
  VectorDirection direction, ForwardOrBackward order, 
  Conjugation conjugation,
  int offset,
  const DistMatrix<std::complex<R>,MC,MR  >& H, 
  const DistMatrix<std::complex<R>,MD,STAR>& t,
        DistMatrix<std::complex<R>,MC,MR  >& A )
{
#ifndef RELEASE
    PushCallStack("advanced::ApplyPackedReflectors");
#endif
    if( side == LEFT )
    {
        if( uplo == LOWER )
        {
            if( direction == VERTICAL && order == FORWARD )
                advanced::internal::ApplyPackedReflectorsLLVF
                ( conjugation, offset, H, t, A );
            else if( direction == VERTICAL )
                advanced::internal::ApplyPackedReflectorsLLVB
                ( conjugation, offset, H, t, A );
            else if( order == FORWARD )
                advanced::internal::ApplyPackedReflectorsLLHF
                ( conjugation, offset, H, t, A );
            else
                advanced::internal::ApplyPackedReflectorsLLHB
                ( conjugation, offset, H, t, A );
        }
        else
        {
            if( direction == VERTICAL && order == FORWARD )
                advanced::internal::ApplyPackedReflectorsLUVF
                ( conjugation, offset, H, t, A );
            else if( direction == VERTICAL )
                advanced::internal::ApplyPackedReflectorsLUVB
                ( conjugation, offset, H, t, A );
            else if( order == FORWARD )
                advanced::internal::ApplyPackedReflectorsLUHF
                ( conjugation, offset, H, t, A );
            else
                advanced::internal::ApplyPackedReflectorsLUHB
                ( conjugation, offset, H, t, A );
        }
    }
    else
    {
        if( uplo == LOWER )
        {
            if( direction == VERTICAL && order == FORWARD )
                advanced::internal::ApplyPackedReflectorsRLVF
                ( conjugation, offset, H, t, A );
            else if( direction == VERTICAL )
                advanced::internal::ApplyPackedReflectorsRLVB
                ( conjugation, offset, H, t, A );
            else if( order == FORWARD )
                advanced::internal::ApplyPackedReflectorsRLHF
                ( conjugation, offset, H, t, A );
            else
                advanced::internal::ApplyPackedReflectorsRLHB
                ( conjugation, offset, H, t, A );
        }
        else
        {
            if( direction == VERTICAL && order == FORWARD )
                advanced::internal::ApplyPackedReflectorsRUVF
                ( conjugation, offset, H, t, A );
            else if( direction == VERTICAL )
                advanced::internal::ApplyPackedReflectorsRUVB
                ( conjugation, offset, H, t, A );
            else if( order == FORWARD )
                advanced::internal::ApplyPackedReflectorsRUHF
                ( conjugation, offset, H, t, A );
            else
                advanced::internal::ApplyPackedReflectorsRUHB
                ( conjugation, offset, H, t, A );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R> // representation of a real number
inline void
elemental::advanced::ApplyPackedReflectors
( Side side, UpperOrLower uplo, 
  VectorDirection direction, ForwardOrBackward order,
  Conjugation conjugation,
  int offset,
  const DistMatrix<std::complex<R>,MC,  MR  >& H, 
  const DistMatrix<std::complex<R>,STAR,STAR>& t,
        DistMatrix<std::complex<R>,MC,  MR  >& A )
{
#ifndef RELEASE
    PushCallStack("advanced::ApplyPackedReflectors");
#endif
    DistMatrix<std::complex<R>,MD,STAR> tDiag(A.Grid());
    tDiag.AlignWithDiagonal( A, offset );
    tDiag = t;
    advanced::ApplyPackedReflectors
    ( side, uplo, direction, order, conjugation, offset, H, tDiag, A );
#ifndef RELEASE
    PopCallStack();
#endif
}
