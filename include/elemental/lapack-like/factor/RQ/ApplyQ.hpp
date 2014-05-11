/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_RQ_APPLYQ_HPP
#define ELEM_RQ_APPLYQ_HPP

#include ELEM_DIAGONALSCALE_INC
#include ELEM_APPLYPACKEDREFLECTORS_INC

namespace elem {
namespace rq {

template<typename F>
inline void
ApplyQ
( LeftOrRight side, Orientation orientation, 
  const Matrix<F>& A, const Matrix<F>& t, const Matrix<Base<F>>& d,
  Matrix<F>& B )
{
    DEBUG_ONLY(CallStackEntry cse("rq::ApplyQ"))
    const bool normal = (orientation==NORMAL);
    const bool onLeft = (side==LEFT);

    const bool applyDFirst = normal!=onLeft;
    if( applyDFirst )
    {
        const Int minDim = d.Height();
        if( onLeft )
        {
            auto BBot = View( B, B.Height()-minDim, 0, minDim, B.Width() );
            DiagonalScale( side, orientation, d, BBot );
        }
        else
        {
            auto BRight = View( B, 0, B.Width()-minDim, B.Height(), minDim );
            DiagonalScale( side, orientation, d, BRight );
        }
    }

    const ForwardOrBackward direction = ( normal==onLeft ? BACKWARD : FORWARD );
    const Conjugation conjugation = ( normal ? CONJUGATED : UNCONJUGATED );
    const Int offset = A.Width()-A.Height();
    ApplyPackedReflectors
    ( side, LOWER, HORIZONTAL, direction, conjugation, offset, A, t, B );

    if( !applyDFirst )
    {
        const Int minDim = d.Height();
        if( onLeft )
        {
            auto BBot = View( B, B.Height()-minDim, 0, minDim, B.Width() );
            DiagonalScale( side, orientation, d, BBot );
        }
        else
        {
            auto BRight = View( B, 0, B.Width()-minDim, B.Height(), minDim );
            DiagonalScale( side, orientation, d, BRight );
        }
    }
}

template<typename F,Dist Ut,Dist Vt,Dist Ud,Dist Vd>
inline void
ApplyQ
( LeftOrRight side, Orientation orientation, 
  const DistMatrix<F>& A, const DistMatrix<F,Ut,Vt>& t, 
  const DistMatrix<Base<F>,Ud,Vd>& d, DistMatrix<F>& B )
{
    DEBUG_ONLY(CallStackEntry cse("rq::ApplyQ"))
    const bool normal = (orientation==NORMAL);
    const bool onLeft = (side==LEFT);

    const bool applyDFirst = normal!=onLeft;
    if( applyDFirst )
    {
        const Int minDim = d.Height();
        if( onLeft )
        {
            auto BBot = View( B, B.Height()-minDim, 0, minDim, B.Width() );
            DiagonalScale( side, orientation, d, BBot );
        }
        else
        {
            auto BRight = View( B, 0, B.Width()-minDim, B.Height(), minDim );
            DiagonalScale( side, orientation, d, BRight );
        }
    }

    const ForwardOrBackward direction = ( normal==onLeft ? BACKWARD : FORWARD );
    const Conjugation conjugation = ( normal ? CONJUGATED : UNCONJUGATED );
    const Int offset = A.Width()-A.Height();

    DistMatrix<F,MD,STAR> tDiag(A.Grid());
    tDiag.SetRoot( A.DiagonalRoot(offset) );
    tDiag.AlignCols( A.DiagonalAlign(offset) );
    tDiag = t;
    ApplyPackedReflectors
    ( side, LOWER, HORIZONTAL, direction, conjugation, offset, A, t, B );

    if( !applyDFirst ) 
    {
        const Int minDim = d.Height();
        if( onLeft )
        {
            auto BBot = View( B, B.Height()-minDim, 0, minDim, B.Width() );
            DiagonalScale( side, orientation, d, BBot );
        }
        else
        {
            auto BRight = View( B, 0, B.Width()-minDim, B.Height(), minDim );
            DiagonalScale( side, orientation, d, BRight );
        }
    }
}

} // namespace rq
} // namespace elem

#endif // ifndef ELEM_RQ_APPLYQ_HPP
