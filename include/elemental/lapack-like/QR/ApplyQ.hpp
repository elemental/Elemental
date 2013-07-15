/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_QR_APPLYQ_HPP
#define LAPACK_QR_APPLYQ_HPP

#include "elemental/lapack-like/ApplyPackedReflectors.hpp"

namespace elem {
namespace qr {

template<typename F>
inline void
ApplyQ
( LeftOrRight side, Orientation orientation, 
  const Matrix<F>& A, const Matrix<F>& t, Matrix<F>& B )
{
#ifndef RELEASE
    CallStackEntry cse("qr::ApplyQ");
#endif
    const bool normal = (orientation==NORMAL);
    const bool onLeft = (side==LEFT);
    const ForwardOrBackward direction = ( normal==onLeft ? BACKWARD : FORWARD );
    const Conjugation conjugation =  ( normal ? UNCONJUGATED : CONJUGATED );
    ApplyPackedReflectors
    ( side, LOWER, VERTICAL, direction, conjugation, 0, A, t, B );
}

template<typename F>
inline void
ApplyQ
( LeftOrRight side, Orientation orientation, 
  const DistMatrix<F>& A, const DistMatrix<F,MD,STAR>& t, DistMatrix<F>& B )
{
#ifndef RELEASE
    CallStackEntry cse("qr::ApplyQ");
#endif
    const bool normal = (orientation==NORMAL);
    const bool onLeft = (side==LEFT);
    const ForwardOrBackward direction = ( normal==onLeft ? BACKWARD : FORWARD );
    const Conjugation conjugation =  ( normal ? UNCONJUGATED : CONJUGATED );
    ApplyPackedReflectors
    ( side, LOWER, VERTICAL, direction, conjugation, 0, A, t, B );
}

template<typename F>
inline void
ApplyQ
( LeftOrRight side, Orientation orientation, 
  const DistMatrix<F>& A, const DistMatrix<F,STAR,STAR>& t, DistMatrix<F>& B )
{
#ifndef RELEASE
    CallStackEntry cse("qr::ApplyQ");
#endif
    DistMatrix<F,MD,STAR> tDiag(A.Grid());
    tDiag.AlignWithDiagonal( A );
    tDiag = t;
    ApplyQ( side, orientation, A, tDiag, B );
}

} // namespace qr
} // namespace elem

#endif // ifndef LAPACK_QR_APPLYQ_HPP
