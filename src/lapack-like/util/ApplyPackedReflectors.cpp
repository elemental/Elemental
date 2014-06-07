/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

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

namespace El {

template<typename F> 
void ApplyPackedReflectors
( LeftOrRight side, UpperOrLower uplo, 
  VerticalOrHorizontal dir, ForwardOrBackward order, 
  Conjugation conjugation,
  Int offset, const Matrix<F>& H, const Matrix<F>& t, Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("ApplyPackedReflectors"))
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
}

template<typename F> 
void ApplyPackedReflectors
( LeftOrRight side, UpperOrLower uplo, 
  VerticalOrHorizontal dir, ForwardOrBackward order, 
  Conjugation conjugation,
  Int offset,
  const DistMatrix<F>& H, const DistMatrix<F,MD,STAR>& t, DistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("ApplyPackedReflectors"))
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
}

template<typename F> 
void ApplyPackedReflectors
( LeftOrRight side, UpperOrLower uplo, 
  VerticalOrHorizontal dir, ForwardOrBackward order,
  Conjugation conjugation,
  Int offset,
  const DistMatrix<F>& H, const DistMatrix<F,STAR,STAR>& t, DistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("ApplyPackedReflectors"))
    DistMatrix<F,MD,STAR> tDiag(A.Grid());
    tDiag.SetRoot( A.DiagonalRoot(offset) );
    tDiag.AlignCols( A.DiagonalAlign(offset) );
    tDiag = t;
    ApplyPackedReflectors
    ( side, uplo, dir, order, conjugation, offset, H, tDiag, A );
}

#define PROTO(F) \
  template void ApplyPackedReflectors \
  ( LeftOrRight side, UpperOrLower uplo, \
    VerticalOrHorizontal dir, ForwardOrBackward order, \
    Conjugation conjugation, Int offset, \
    const Matrix<F>& H, const Matrix<F>& t, \
          Matrix<F>& A ); \
  template void ApplyPackedReflectors \
  ( LeftOrRight side, UpperOrLower uplo, \
    VerticalOrHorizontal dir, ForwardOrBackward order, \
    Conjugation conjugation, Int offset, \
    const DistMatrix<F>& H, const DistMatrix<F,MD,STAR>& t, \
          DistMatrix<F>& A ); \
  template void ApplyPackedReflectors \
  ( LeftOrRight side, UpperOrLower uplo, \
    VerticalOrHorizontal dir, ForwardOrBackward order, \
    Conjugation conjugation, Int offset, \
    const DistMatrix<F>& H, const DistMatrix<F,STAR,STAR>& t, \
          DistMatrix<F>& A ); 

PROTO(float)
PROTO(double)
PROTO(Complex<float>)
PROTO(Complex<double>)

} // namespace El
