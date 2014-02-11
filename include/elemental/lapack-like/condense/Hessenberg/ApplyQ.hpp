/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_HESSENBERG_APPLYQ_HPP
#define ELEM_HESSENBERG_APPLYQ_HPP

#include ELEM_APPLYPACKEDREFLECTORS_INC

namespace elem {
namespace hessenberg {

template<typename F>
inline void
ApplyQ
( UpperOrLower uplo, LeftOrRight side, Orientation orientation, 
  const Matrix<F>& A, const Matrix<F>& t, Matrix<F>& H )
{
    DEBUG_ONLY(CallStackEntry cse("hessenberg::ApplyQ"))
    const bool normal = (orientation==NORMAL);
    const bool onLeft = (side==LEFT);
    const ForwardOrBackward direction = ( normal==onLeft ? BACKWARD : FORWARD );
    if( uplo == LOWER )
    {
        const Conjugation conjugation = ( normal ? UNCONJUGATED : CONJUGATED );
        ApplyPackedReflectors
        ( side, UPPER, HORIZONTAL, direction, conjugation, 1, A, t, H );
    }
    else
    {
        const Conjugation conjugation = ( normal ? CONJUGATED : UNCONJUGATED );
        ApplyPackedReflectors
        ( side, LOWER, VERTICAL, direction, conjugation, -1, A, t, H );
    }
}

template<typename F>
inline void
ApplyQ
( UpperOrLower uplo, LeftOrRight side, Orientation orientation, 
  const DistMatrix<F>& A, const DistMatrix<F,MD,STAR>& t, DistMatrix<F>& H )
{
    DEBUG_ONLY(CallStackEntry cse("hessenberg::ApplyQ"))
    const bool normal = (orientation==NORMAL);
    const bool onLeft = (side==LEFT);
    const ForwardOrBackward direction = ( normal==onLeft ? BACKWARD : FORWARD );
    if( uplo == LOWER )
    {
        const Conjugation conjugation = ( normal ? UNCONJUGATED : CONJUGATED );
        ApplyPackedReflectors
        ( side, UPPER, HORIZONTAL, direction, conjugation, 1, A, t, H );
    }
    else
    {
        const Conjugation conjugation = ( normal ? CONJUGATED : UNCONJUGATED );
        ApplyPackedReflectors
        ( side, LOWER, VERTICAL, direction, conjugation, -1, A, t, H );
    }
}

template<typename F>
inline void
ApplyQ
( UpperOrLower uplo, LeftOrRight side, Orientation orientation, 
  const DistMatrix<F>& A, const DistMatrix<F,STAR,STAR>& t, DistMatrix<F>& H )
{
    DEBUG_ONLY(CallStackEntry cse("hessenberg::ApplyQ"))
    const Int offset = ( uplo==LOWER ? 1 : -1 );
    DistMatrix<F,MD,STAR> tDiag(A.Grid());
    tDiag.SetRoot( A.DiagonalRoot(offset) );
    tDiag.AlignCols( A.DiagonalAlign(offset) );
    tDiag = t;
    ApplyQ( uplo, side, orientation, A, tDiag, H );
}

} // namespace hessenberg
} // namespace elem

#endif // ifndef ELEM_HESSENBERG_APPLYQ_HPP
