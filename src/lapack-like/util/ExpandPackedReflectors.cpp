/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

#include "./ApplyPackedReflectors/Util.hpp"
#include "./ExpandPackedReflectors/LV.hpp"

namespace El {

template<typename F> 
void ExpandPackedReflectors
( UpperOrLower uplo, VerticalOrHorizontal dir, Conjugation conjugation,
  Int offset, Matrix<F>& H, const Matrix<F>& t )
{
    DEBUG_ONLY(CallStackEntry cse("ExpandPackedReflectors"))
    if( uplo == LOWER && dir == VERTICAL )
        expand_packed_reflectors::LV( conjugation, offset, H, t );
    else
        LogicError("This option is not yet supported");
}

template<typename F> 
void ExpandPackedReflectors
( UpperOrLower uplo, VerticalOrHorizontal dir, Conjugation conjugation,
  Int offset, DistMatrix<F>& H, const DistMatrix<F,MD,STAR>& t )
{
    DEBUG_ONLY(CallStackEntry cse("ExpandPackedReflectors"))
    if( uplo == LOWER && dir == VERTICAL )
        expand_packed_reflectors::LV( conjugation, offset, H, t );
    else
        LogicError("This option is not yet supported");
}

template<typename F> 
void ExpandPackedReflectors
( UpperOrLower uplo, VerticalOrHorizontal dir, Conjugation conjugation,
  Int offset, DistMatrix<F>& H, const DistMatrix<F,STAR,STAR>& t )
{
    DEBUG_ONLY(CallStackEntry cse("ExpandPackedReflectors"))
    DistMatrix<F,MD,STAR> tDiag(H.Grid());
    tDiag.SetRoot( H.DiagonalRoot(offset) );
    tDiag.AlignCols( H.DiagonalAlign(offset) );
    tDiag = t;
    ExpandPackedReflectors( uplo, dir, conjugation, offset, H, tDiag );
}

#define PROTO(F) \
  template void ExpandPackedReflectors \
  ( UpperOrLower uplo, VerticalOrHorizontal dir, Conjugation conjugation, \
    Int offset, Matrix<F>& H, const Matrix<F>& t ); \
  template void ExpandPackedReflectors \
  ( UpperOrLower uplo, VerticalOrHorizontal dir, Conjugation conjugation, \
    Int offset, DistMatrix<F>& H, const DistMatrix<F,MD,STAR>& t ); \
  template void ExpandPackedReflectors \
  ( UpperOrLower uplo, VerticalOrHorizontal dir, Conjugation conjugation, \
    Int offset, DistMatrix<F>& H, const DistMatrix<F,STAR,STAR>& t ); 

PROTO(float)
PROTO(double)
PROTO(Complex<float>)
PROTO(Complex<double>)

} // namespace El
