/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_HERMITIANTRIDIAG_APPLYQ_HPP
#define ELEM_HERMITIANTRIDIAG_APPLYQ_HPP

#include ELEM_APPLYPACKEDREFLECTORS_INC

namespace elem {
namespace herm_tridiag {

template<typename F>
inline void
ApplyQ
( LeftOrRight side, UpperOrLower uplo, Orientation orientation, 
  const Matrix<F>& A, const Matrix<F>& t, Matrix<F>& B )
{
    DEBUG_ONLY(CallStackEntry cse("herm_tridiag::ApplyQ"))
    const bool normal = (orientation==NORMAL);
    const bool onLeft = (side==LEFT);
    const ForwardOrBackward direction = 
        ( (normal==onLeft) ^ (uplo==UPPER) ? BACKWARD : FORWARD );
    const Conjugation conjugation = ( normal ? CONJUGATED : UNCONJUGATED );
    const Int offset = ( uplo==UPPER ? 1 : -1 );
    ApplyPackedReflectors
    ( side, uplo, VERTICAL, direction, conjugation, offset, A, t, B );
}

template<typename F>
inline void
ApplyQ
( LeftOrRight side, UpperOrLower uplo, Orientation orientation, 
  const DistMatrix<F>& A, const DistMatrix<F,MD,STAR>& t, DistMatrix<F>& B )
{
    DEBUG_ONLY(CallStackEntry cse("hermi_tridiag::ApplyQ"))
    const bool normal = (orientation==NORMAL);
    const bool onLeft = (side==LEFT);
    const ForwardOrBackward direction = 
        ( (normal==onLeft) ^ (uplo==UPPER) ? BACKWARD : FORWARD );
    const Conjugation conjugation = ( normal ? CONJUGATED : UNCONJUGATED );
    const Int offset = ( uplo==UPPER ? 1 : -1 );
    ApplyPackedReflectors
    ( side, uplo, VERTICAL, direction, conjugation, offset, A, t, B );
}

template<typename F>
inline void
ApplyQ
( LeftOrRight side, UpperOrLower uplo, Orientation orientation, 
  const DistMatrix<F>& A, const DistMatrix<F,STAR,STAR>& t, DistMatrix<F>& B )
{
    DEBUG_ONLY(CallStackEntry cse("herm_tridiag::ApplyQ"))
    const Int offset = ( uplo==UPPER ? 1 : -1 );
    DistMatrix<F,MD,STAR> tDiag(A.Grid());
    tDiag.SetRoot( A.DiagonalRoot(offset) );
    tDiag.AlignCols( A.DiagonalAlign(offset) );
    tDiag = t;
    ApplyQ( side, uplo, orientation, A, tDiag, B );
}

} // namespace herm_tridiag
} // namespace elem

#endif // ifndef ELEM_HERMITIANTRIDIAG_APPLYQ_HPP
