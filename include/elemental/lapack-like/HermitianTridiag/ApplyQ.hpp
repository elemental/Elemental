/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_HERMITIANTRIDIAG_APPLYQ_HPP
#define LAPACK_HERMITIANTRIDIAG_APPLYQ_HPP

#include "elemental/lapack-like/ApplyPackedReflectors.hpp"

namespace elem {
namespace hermitian_tridiag {

template<typename F>
inline void
ApplyQ
( LeftOrRight side, UpperOrLower uplo, Orientation orientation, 
  const Matrix<F>& A, const Matrix<F>& t, Matrix<F>& B )
{
#ifndef RELEASE
    CallStackEntry cse("hermitian_tridiag::ApplyQ");
#endif
    const bool normal = (orientation==NORMAL);
    const bool onLeft = (side==LEFT);
    const ForwardOrBackward direction = ( normal==onLeft ? BACKWARD : FORWARD );
    const Conjugation conjugation = ( normal ? UNCONJUGATED : CONJUGATED );
    const int offset = ( uplo==UPPER ? 1 : -1 );
    ApplyPackedReflectors
    ( side, uplo, VERTICAL, direction, conjugation, offset, A, t, B );
}

template<typename F>
inline void
ApplyQ
( LeftOrRight side, UpperOrLower uplo, Orientation orientation, 
  const DistMatrix<F>& A, const DistMatrix<F,MD,STAR>& t, DistMatrix<F>& B )
{
#ifndef RELEASE
    CallStackEntry cse("hermitian_tridiag::ApplyQ");
#endif
    const bool normal = (orientation==NORMAL);
    const bool onLeft = (side==LEFT);
    const ForwardOrBackward direction = ( normal==onLeft ? BACKWARD : FORWARD );
    const Conjugation conjugation = ( normal ? UNCONJUGATED : CONJUGATED );
    const int offset = ( uplo==UPPER ? 1 : -1 );
    ApplyPackedReflectors
    ( side, uplo, VERTICAL, direction, conjugation, offset, A, t, B );
}

template<typename F>
inline void
ApplyQ
( LeftOrRight side, UpperOrLower uplo, Orientation orientation, 
  const DistMatrix<F>& A, const DistMatrix<F,STAR,STAR>& t, DistMatrix<F>& B )
{
#ifndef RELEASE
    CallStackEntry cse("hermitian_tridiag::ApplyQ");
#endif
    const int offset = ( uplo==UPPER ? 1 : -1 );
    DistMatrix<F,MD,STAR> tDiag(A.Grid());
    tDiag.AlignWithDiagonal( A, offset );
    tDiag = t;
    ApplyQ( side, uplo, orientation, A, tDiag, B );
}

} // namespace hermitian_tridiag
} // namespace elem

#endif // ifndef LAPACK_HERMITIANTRIDIAG_APPLYQ_HPP
