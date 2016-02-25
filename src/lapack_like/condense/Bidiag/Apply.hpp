/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_BIDIAG_APPLY_HPP
#define EL_BIDIAG_APPLY_HPP

namespace El {
namespace bidiag {

template<typename F>
void ApplyQ
( LeftOrRight side, Orientation orientation, 
  const Matrix<F>& A, const Matrix<F>& t, Matrix<F>& B )
{
    DEBUG_ONLY(CSE cse("bidiag::ApplyQ"))
    const bool normal = (orientation==NORMAL);
    const bool onLeft = (side==LEFT);
    const ForwardOrBackward direction = ( normal==onLeft ? BACKWARD : FORWARD );
    const Conjugation conjugation = ( normal ? CONJUGATED : UNCONJUGATED );
    const Int offset = ( A.Height()>=A.Width() ? 0 : -1 );
    ApplyPackedReflectors
    ( side, LOWER, VERTICAL, direction, conjugation, offset, A, t, B );
}

template<typename F>
void ApplyP
( LeftOrRight side, Orientation orientation, 
  const Matrix<F>& A, const Matrix<F>& t, Matrix<F>& B )
{
    DEBUG_ONLY(CSE cse("bidiag::ApplyP"))
    const bool normal = (orientation==NORMAL);
    const bool onLeft = (side==LEFT);
    const ForwardOrBackward direction = ( normal==onLeft ? BACKWARD : FORWARD );
    const Conjugation conjugation = ( normal ? UNCONJUGATED : CONJUGATED );
    const Int offset = ( A.Height()>=A.Width() ? 1 : 0 );
    ApplyPackedReflectors
    ( side, UPPER, HORIZONTAL, direction, conjugation, offset, A, t, B );
}

template<typename F>
void ApplyQ
( LeftOrRight side, Orientation orientation, 
  const ElementalMatrix<F>& A, const ElementalMatrix<F>& t, 
        ElementalMatrix<F>& B )
{
    DEBUG_ONLY(CSE cse("bidiag::ApplyQ"))
    const bool normal = (orientation==NORMAL);
    const bool onLeft = (side==LEFT);
    const ForwardOrBackward direction = ( normal==onLeft ? BACKWARD : FORWARD );
    const Conjugation conjugation = ( normal ? CONJUGATED : UNCONJUGATED );
    const Int offset = ( A.Height()>=A.Width() ? 0 : -1 );
    ApplyPackedReflectors
    ( side, LOWER, VERTICAL, direction, conjugation, offset, A, t, B );
}

template<typename F>
void ApplyP
( LeftOrRight side, Orientation orientation, 
  const ElementalMatrix<F>& A, const ElementalMatrix<F>& t, 
        ElementalMatrix<F>& B )
{
    DEBUG_ONLY(CSE cse("bidiag::ApplyP"))
    const bool normal = (orientation==NORMAL);
    const bool onLeft = (side==LEFT);
    const ForwardOrBackward direction = ( normal==onLeft ? BACKWARD : FORWARD );
    const Conjugation conjugation = ( normal ? UNCONJUGATED : CONJUGATED );
    const Int offset = ( A.Height()>=A.Width() ? 1 : 0 );
    ApplyPackedReflectors
    ( side, UPPER, HORIZONTAL, direction, conjugation, offset, A, t, B );
}

} // namespace bidiag
} // namespace El

#endif // ifndef EL_BIDIAG_APPLY_HPP
