/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_HERMITIANTRIDIAG_APPLYQ_HPP
#define EL_HERMITIANTRIDIAG_APPLYQ_HPP

namespace El {
namespace herm_tridiag {

template<typename F>
void ApplyQ
( LeftOrRight side, UpperOrLower uplo, Orientation orientation, 
  const Matrix<F>& A, const Matrix<F>& t, Matrix<F>& B )
{
    DEBUG_ONLY(CSE cse("herm_tridiag::ApplyQ"))
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
void ApplyQ
( LeftOrRight side, UpperOrLower uplo, Orientation orientation, 
  const ElementalMatrix<F>& A, const ElementalMatrix<F>& t, 
        ElementalMatrix<F>& B )
{
    DEBUG_ONLY(CSE cse("hermi_tridiag::ApplyQ"))
    const bool normal = (orientation==NORMAL);
    const bool onLeft = (side==LEFT);
    const ForwardOrBackward direction = 
        ( (normal==onLeft) ^ (uplo==UPPER) ? BACKWARD : FORWARD );
    const Conjugation conjugation = ( normal ? CONJUGATED : UNCONJUGATED );
    const Int offset = ( uplo==UPPER ? 1 : -1 );
    ApplyPackedReflectors
    ( side, uplo, VERTICAL, direction, conjugation, offset, A, t, B );
}

} // namespace herm_tridiag
} // namespace El

#endif // ifndef EL_HERMITIANTRIDIAG_APPLYQ_HPP
