/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
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
  const DistMatrix<R>& H, 
        DistMatrix<R>& A )
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
  const DistMatrix<Complex<R> >& H, 
  const DistMatrix<Complex<R>,MD,STAR>& t,
        DistMatrix<Complex<R> >& A )
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
  const DistMatrix<Complex<R> >& H, 
  const DistMatrix<Complex<R>,STAR,STAR>& t,
        DistMatrix<Complex<R> >& A )
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
