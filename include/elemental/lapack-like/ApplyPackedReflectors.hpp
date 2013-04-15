/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_APPLYPACKEDREFLECTORS_HPP
#define LAPACK_APPLYPACKEDREFLECTORS_HPP

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
                apply_packed_reflectors::LLVF( offset, H, A );
            else if( dir == VERTICAL )
                apply_packed_reflectors::LLVB( offset, H, A );
            else if( order == FORWARD )
                apply_packed_reflectors::LLHF( offset, H, A );
            else
                apply_packed_reflectors::LLHB( offset, H, A );
        }
        else
        {
            if( dir == VERTICAL && order == FORWARD )
                apply_packed_reflectors::LUVF( offset, H, A );
            else if( dir == VERTICAL )
                apply_packed_reflectors::LUVB( offset, H, A );
            else if( order == FORWARD )
                apply_packed_reflectors::LUHF( offset, H, A );
            else
                apply_packed_reflectors::LUHB( offset, H, A );
        }
    }
    else
    {
        if( uplo == LOWER )
        {
            if( dir == VERTICAL && order == FORWARD )
                apply_packed_reflectors::RLVF( offset, H, A );
            else if( dir == VERTICAL )
                apply_packed_reflectors::RLVB( offset, H, A );
            else if( order == FORWARD )
                apply_packed_reflectors::RLHF( offset, H, A );
            else
                apply_packed_reflectors::RLHB( offset, H, A );
        }
        else
        {
            if( dir == VERTICAL && order == FORWARD )
                apply_packed_reflectors::RUVF( offset, H, A );
            else if( dir == VERTICAL )
                apply_packed_reflectors::RUVB( offset, H, A );
            else if( order == FORWARD )
                apply_packed_reflectors::RUHF( offset, H, A );
            else
                apply_packed_reflectors::RUHB( offset, H, A );
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
                apply_packed_reflectors::LLVF( offset, H, A );
            else if( dir == VERTICAL )
                apply_packed_reflectors::LLVB( offset, H, A );
            else if( order == FORWARD )
                apply_packed_reflectors::LLHF( offset, H, A );
            else
                apply_packed_reflectors::LLHB( offset, H, A );
        }
        else
        {
            if( dir == VERTICAL && order == FORWARD )
                apply_packed_reflectors::LUVF( offset, H, A );
            else if( dir == VERTICAL )
                apply_packed_reflectors::LUVB( offset, H, A );
            else if( order == FORWARD )
                apply_packed_reflectors::LUHF( offset, H, A );
            else
                apply_packed_reflectors::LUHB( offset, H, A );
        }
    }
    else
    {
        if( uplo == LOWER )
        {
            if( dir == VERTICAL && order == FORWARD )
                apply_packed_reflectors::RLVF( offset, H, A );
            else if( dir == VERTICAL )
                apply_packed_reflectors::RLVB( offset, H, A );
            else if( order == FORWARD )
                apply_packed_reflectors::RLHF( offset, H, A );
            else
                apply_packed_reflectors::RLHB( offset, H, A );
        }
        else
        {
            if( dir == VERTICAL && order == FORWARD )
                apply_packed_reflectors::RUVF( offset, H, A );
            else if( dir == VERTICAL )
                apply_packed_reflectors::RUVB( offset, H, A );
            else if( order == FORWARD )
                apply_packed_reflectors::RUHF( offset, H, A );
            else
                apply_packed_reflectors::RUHB( offset, H, A );
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
                apply_packed_reflectors::LLVF( conjugation, offset, H, t, A );
            else if( dir == VERTICAL )
                apply_packed_reflectors::LLVB( conjugation, offset, H, t, A );
            else if( order == FORWARD )
                apply_packed_reflectors::LLHF( conjugation, offset, H, t, A );
            else
                apply_packed_reflectors::LLHB( conjugation, offset, H, t, A );
        }
        else
        {
            if( dir == VERTICAL && order == FORWARD )
                apply_packed_reflectors::LUVF( conjugation, offset, H, t, A );
            else if( dir == VERTICAL )
                apply_packed_reflectors::LUVB( conjugation, offset, H, t, A );
            else if( order == FORWARD )
                apply_packed_reflectors::LUHF( conjugation, offset, H, t, A );
            else
                apply_packed_reflectors::LUHB( conjugation, offset, H, t, A );
        }
    }
    else
    {
        if( uplo == LOWER )
        {
            if( dir == VERTICAL && order == FORWARD )
                apply_packed_reflectors::RLVF( conjugation, offset, H, t, A );
            else if( dir == VERTICAL )
                apply_packed_reflectors::RLVB( conjugation, offset, H, t, A );
            else if( order == FORWARD )
                apply_packed_reflectors::RLHF( conjugation, offset, H, t, A );
            else
                apply_packed_reflectors::RLHB( conjugation, offset, H, t, A );
        }
        else
        {
            if( dir == VERTICAL && order == FORWARD )
                apply_packed_reflectors::RUVF( conjugation, offset, H, t, A );
            else if( dir == VERTICAL )
                apply_packed_reflectors::RUVB( conjugation, offset, H, t, A );
            else if( order == FORWARD )
                apply_packed_reflectors::RUHF( conjugation, offset, H, t, A );
            else
                apply_packed_reflectors::RUHB( conjugation, offset, H, t, A );
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
                apply_packed_reflectors::LLVF( conjugation, offset, H, t, A );
            else if( dir == VERTICAL )
                apply_packed_reflectors::LLVB( conjugation, offset, H, t, A );
            else if( order == FORWARD )
                apply_packed_reflectors::LLHF( conjugation, offset, H, t, A );
            else
                apply_packed_reflectors::LLHB( conjugation, offset, H, t, A );
        }
        else
        {
            if( dir == VERTICAL && order == FORWARD )
                apply_packed_reflectors::LUVF( conjugation, offset, H, t, A );
            else if( dir == VERTICAL )
                apply_packed_reflectors::LUVB( conjugation, offset, H, t, A );
            else if( order == FORWARD )
                apply_packed_reflectors::LUHF( conjugation, offset, H, t, A );
            else
                apply_packed_reflectors::LUHB( conjugation, offset, H, t, A );
        }
    }
    else
    {
        if( uplo == LOWER )
        {
            if( dir == VERTICAL && order == FORWARD )
                apply_packed_reflectors::RLVF( conjugation, offset, H, t, A );
            else if( dir == VERTICAL )
                apply_packed_reflectors::RLVB( conjugation, offset, H, t, A );
            else if( order == FORWARD )
                apply_packed_reflectors::RLHF( conjugation, offset, H, t, A );
            else
                apply_packed_reflectors::RLHB( conjugation, offset, H, t, A );
        }
        else
        {
            if( dir == VERTICAL && order == FORWARD )
                apply_packed_reflectors::RUVF( conjugation, offset, H, t, A );
            else if( dir == VERTICAL )
                apply_packed_reflectors::RUVB( conjugation, offset, H, t, A );
            else if( order == FORWARD )
                apply_packed_reflectors::RUHF( conjugation, offset, H, t, A );
            else
                apply_packed_reflectors::RUHB( conjugation, offset, H, t, A );
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

#endif // ifndef LAPACK_APPLYPACKEDREFLECTORS_HPP
