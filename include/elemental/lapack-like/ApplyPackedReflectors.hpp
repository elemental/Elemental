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

#include "./ApplyPackedReflectors/Util.hpp"
#include "./ApplyPackedReflectors/LLHB.hpp"
#include "./ApplyPackedReflectors/LLHF.hpp"
#include "./ApplyPackedReflectors/LLVB.hpp"
#include "./ApplyPackedReflectors/LLVF.hpp"
#include "./ApplyPackedReflectors/LUHB.hpp"
#include "./ApplyPackedReflectors/LUHF.hpp"
#include "./ApplyPackedReflectors/LUVB.hpp"
#include "./ApplyPackedReflectors/LUVF.hpp"
#include "./ApplyPackedReflectors/RLHB.hpp"
#include "./ApplyPackedReflectors/RLHF.hpp"
#include "./ApplyPackedReflectors/RLVB.hpp"
#include "./ApplyPackedReflectors/RLVF.hpp"
#include "./ApplyPackedReflectors/RUHB.hpp"
#include "./ApplyPackedReflectors/RUHF.hpp"
#include "./ApplyPackedReflectors/RUVB.hpp"
#include "./ApplyPackedReflectors/RUVF.hpp"

namespace elem {

template<typename R> 
inline void
ApplyPackedReflectors
( LeftOrRight side, UpperOrLower uplo, 
  VerticalOrHorizontal dir, ForwardOrBackward order,
  int offset, const Matrix<R>& H, Matrix<R>& A )
{
#ifndef RELEASE
    PushCallStack("ApplyPackedReflectors");
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
            if( dir == VERTICAL && order == FORWARD )
                internal::ApplyPackedReflectorsLLVF( offset, H, A );
            else if( dir == VERTICAL )
                internal::ApplyPackedReflectorsLLVB( offset, H, A );
            else if( order == FORWARD )
                internal::ApplyPackedReflectorsLLHF( offset, H, A );
            else
                internal::ApplyPackedReflectorsLLHB( offset, H, A );
        }
        else
        {
            if( dir == VERTICAL && order == FORWARD )
                internal::ApplyPackedReflectorsLUVF( offset, H, A );
            else if( dir == VERTICAL )
                internal::ApplyPackedReflectorsLUVB( offset, H, A );
            else if( order == FORWARD )
                internal::ApplyPackedReflectorsLUHF( offset, H, A );
            else
                internal::ApplyPackedReflectorsLUHB( offset, H, A );
        }
    }
    else
    {
        if( uplo == LOWER )
        {
            if( dir == VERTICAL && order == FORWARD )
                internal::ApplyPackedReflectorsRLVF( offset, H, A );
            else if( dir == VERTICAL )
                internal::ApplyPackedReflectorsRLVB( offset, H, A );
            else if( order == FORWARD )
                internal::ApplyPackedReflectorsRLHF( offset, H, A );
            else
                internal::ApplyPackedReflectorsRLHB( offset, H, A );
        }
        else
        {
            if( dir == VERTICAL && order == FORWARD )
                internal::ApplyPackedReflectorsRUVF( offset, H, A );
            else if( dir == VERTICAL )
                internal::ApplyPackedReflectorsRUVB( offset, H, A );
            else if( order == FORWARD )
                internal::ApplyPackedReflectorsRUHF( offset, H, A );
            else
                internal::ApplyPackedReflectorsRUHB( offset, H, A );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R> 
inline void
ApplyPackedReflectors
( LeftOrRight side, UpperOrLower uplo, 
  VerticalOrHorizontal dir, ForwardOrBackward order,
  int offset,
  const DistMatrix<R,MC,MR>& H, 
        DistMatrix<R,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("ApplyPackedReflectors");
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
            if( dir == VERTICAL && order == FORWARD )
                internal::ApplyPackedReflectorsLLVF( offset, H, A );
            else if( dir == VERTICAL )
                internal::ApplyPackedReflectorsLLVB( offset, H, A );
            else if( order == FORWARD )
                internal::ApplyPackedReflectorsLLHF( offset, H, A );
            else
                internal::ApplyPackedReflectorsLLHB( offset, H, A );
        }
        else
        {
            if( dir == VERTICAL && order == FORWARD )
                internal::ApplyPackedReflectorsLUVF( offset, H, A );
            else if( dir == VERTICAL )
                internal::ApplyPackedReflectorsLUVB( offset, H, A );
            else if( order == FORWARD )
                internal::ApplyPackedReflectorsLUHF( offset, H, A );
            else
                internal::ApplyPackedReflectorsLUHB( offset, H, A );
        }
    }
    else
    {
        if( uplo == LOWER )
        {
            if( dir == VERTICAL && order == FORWARD )
                internal::ApplyPackedReflectorsRLVF( offset, H, A );
            else if( dir == VERTICAL )
                internal::ApplyPackedReflectorsRLVB( offset, H, A );
            else if( order == FORWARD )
                internal::ApplyPackedReflectorsRLHF( offset, H, A );
            else
                internal::ApplyPackedReflectorsRLHB( offset, H, A );
        }
        else
        {
            if( dir == VERTICAL && order == FORWARD )
                internal::ApplyPackedReflectorsRUVF( offset, H, A );
            else if( dir == VERTICAL )
                internal::ApplyPackedReflectorsRUVB( offset, H, A );
            else if( order == FORWARD )
                internal::ApplyPackedReflectorsRUHF( offset, H, A );
            else
                internal::ApplyPackedReflectorsRUHB( offset, H, A );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R> 
inline void
ApplyPackedReflectors
( LeftOrRight side, UpperOrLower uplo, 
  VerticalOrHorizontal dir, ForwardOrBackward order, 
  Conjugation conjugation,
  int offset,
  const Matrix<Complex<R> >& H, 
  const Matrix<Complex<R> >& t,
        Matrix<Complex<R> >& A )
{
#ifndef RELEASE
    PushCallStack("ApplyPackedReflectors");
#endif
    if( side == LEFT )
    {
        if( uplo == LOWER )
        {
            if( dir == VERTICAL && order == FORWARD )
                internal::ApplyPackedReflectorsLLVF
                ( conjugation, offset, H, t, A );
            else if( dir == VERTICAL )
                internal::ApplyPackedReflectorsLLVB
                ( conjugation, offset, H, t, A );
            else if( order == FORWARD )
                internal::ApplyPackedReflectorsLLHF
                ( conjugation, offset, H, t, A );
            else
                internal::ApplyPackedReflectorsLLHB
                ( conjugation, offset, H, t, A );
        }
        else
        {
            if( dir == VERTICAL && order == FORWARD )
                internal::ApplyPackedReflectorsLUVF
                ( conjugation, offset, H, t, A );
            else if( dir == VERTICAL )
                internal::ApplyPackedReflectorsLUVB
                ( conjugation, offset, H, t, A );
            else if( order == FORWARD )
                internal::ApplyPackedReflectorsLUHF
                ( conjugation, offset, H, t, A );
            else
                internal::ApplyPackedReflectorsLUHB
                ( conjugation, offset, H, t, A );
        }
    }
    else
    {
        if( uplo == LOWER )
        {
            if( dir == VERTICAL && order == FORWARD )
                internal::ApplyPackedReflectorsRLVF
                ( conjugation, offset, H, t, A );
            else if( dir == VERTICAL )
                internal::ApplyPackedReflectorsRLVB
                ( conjugation, offset, H, t, A );
            else if( order == FORWARD )
                internal::ApplyPackedReflectorsRLHF
                ( conjugation, offset, H, t, A );
            else
                internal::ApplyPackedReflectorsRLHB
                ( conjugation, offset, H, t, A );
        }
        else
        {
            if( dir == VERTICAL && order == FORWARD )
                internal::ApplyPackedReflectorsRUVF
                ( conjugation, offset, H, t, A );
            else if( dir == VERTICAL )
                internal::ApplyPackedReflectorsRUVB
                ( conjugation, offset, H, t, A );
            else if( order == FORWARD )
                internal::ApplyPackedReflectorsRUHF
                ( conjugation, offset, H, t, A );
            else
                internal::ApplyPackedReflectorsRUHB
                ( conjugation, offset, H, t, A );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R> 
inline void
ApplyPackedReflectors
( LeftOrRight side, UpperOrLower uplo, 
  VerticalOrHorizontal dir, ForwardOrBackward order, 
  Conjugation conjugation,
  int offset,
  const DistMatrix<Complex<R>,MC,MR  >& H, 
  const DistMatrix<Complex<R>,MD,STAR>& t,
        DistMatrix<Complex<R>,MC,MR  >& A )
{
#ifndef RELEASE
    PushCallStack("ApplyPackedReflectors");
#endif
    if( side == LEFT )
    {
        if( uplo == LOWER )
        {
            if( dir == VERTICAL && order == FORWARD )
                internal::ApplyPackedReflectorsLLVF
                ( conjugation, offset, H, t, A );
            else if( dir == VERTICAL )
                internal::ApplyPackedReflectorsLLVB
                ( conjugation, offset, H, t, A );
            else if( order == FORWARD )
                internal::ApplyPackedReflectorsLLHF
                ( conjugation, offset, H, t, A );
            else
                internal::ApplyPackedReflectorsLLHB
                ( conjugation, offset, H, t, A );
        }
        else
        {
            if( dir == VERTICAL && order == FORWARD )
                internal::ApplyPackedReflectorsLUVF
                ( conjugation, offset, H, t, A );
            else if( dir == VERTICAL )
                internal::ApplyPackedReflectorsLUVB
                ( conjugation, offset, H, t, A );
            else if( order == FORWARD )
                internal::ApplyPackedReflectorsLUHF
                ( conjugation, offset, H, t, A );
            else
                internal::ApplyPackedReflectorsLUHB
                ( conjugation, offset, H, t, A );
        }
    }
    else
    {
        if( uplo == LOWER )
        {
            if( dir == VERTICAL && order == FORWARD )
                internal::ApplyPackedReflectorsRLVF
                ( conjugation, offset, H, t, A );
            else if( dir == VERTICAL )
                internal::ApplyPackedReflectorsRLVB
                ( conjugation, offset, H, t, A );
            else if( order == FORWARD )
                internal::ApplyPackedReflectorsRLHF
                ( conjugation, offset, H, t, A );
            else
                internal::ApplyPackedReflectorsRLHB
                ( conjugation, offset, H, t, A );
        }
        else
        {
            if( dir == VERTICAL && order == FORWARD )
                internal::ApplyPackedReflectorsRUVF
                ( conjugation, offset, H, t, A );
            else if( dir == VERTICAL )
                internal::ApplyPackedReflectorsRUVB
                ( conjugation, offset, H, t, A );
            else if( order == FORWARD )
                internal::ApplyPackedReflectorsRUHF
                ( conjugation, offset, H, t, A );
            else
                internal::ApplyPackedReflectorsRUHB
                ( conjugation, offset, H, t, A );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R> 
inline void
ApplyPackedReflectors
( LeftOrRight side, UpperOrLower uplo, 
  VerticalOrHorizontal dir, ForwardOrBackward order,
  Conjugation conjugation,
  int offset,
  const DistMatrix<Complex<R>,MC,  MR  >& H, 
  const DistMatrix<Complex<R>,STAR,STAR>& t,
        DistMatrix<Complex<R>,MC,  MR  >& A )
{
#ifndef RELEASE
    PushCallStack("ApplyPackedReflectors");
#endif
    DistMatrix<Complex<R>,MD,STAR> tDiag(A.Grid());
    tDiag.AlignWithDiagonal( A, offset );
    tDiag = t;
    ApplyPackedReflectors
    ( side, uplo, dir, order, conjugation, offset, H, tDiag, A );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
