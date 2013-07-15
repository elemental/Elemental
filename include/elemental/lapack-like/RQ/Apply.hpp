/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_RQ_APPLY_HPP
#define LAPACK_RQ_APPLY_HPP

#include "elemental/lapack-like/ApplyPackedReflectors.hpp"

namespace elem {
namespace rq {

template<typename F>
inline void
Apply
( LeftOrRight side, Orientation orientation, 
  const Matrix<F>& A, const Matrix<F>& t, Matrix<F>& B )
{
#ifndef RELEASE
    CallStackEntry cse("rq::Apply");
#endif
    const ForwardOrBackward direction = 
        ( orientation==NORMAL ? FORWARD : BACKWARD );
    const Conjugation conjugation = 
        ( orientation==NORMAL ? UNCONJUGATED : CONJUGATED );
    const int offset = A.Width()-A.Height();
    ApplyPackedReflectors
    ( side, LOWER, HORIZONTAL, direction, conjugation, offset, A, t, B );
}

template<typename F>
inline void
Apply
( LeftOrRight side, Orientation orientation, 
  const DistMatrix<F>& A, const DistMatrix<F,MD,STAR>& t, DistMatrix<F>& B )
{
#ifndef RELEASE
    CallStackEntry cse("rq::Apply");
#endif
    const ForwardOrBackward direction = 
        ( orientation==NORMAL ? FORWARD : BACKWARD );
    const Conjugation conjugation =
        ( orientation==NORMAL ? UNCONJUGATED : CONJUGATED );
    const int offset = A.Width()-A.Height();
    ApplyPackedReflectors
    ( side, LOWER, HORIZONTAL, direction, conjugation, offset, A, t, B );
}

} // namespace rq
} // namespace elem

#endif // ifndef LAPACK_RQ_APPLY_HPP
