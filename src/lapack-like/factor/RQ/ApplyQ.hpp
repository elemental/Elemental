/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_RQ_APPLYQ_HPP
#define EL_RQ_APPLYQ_HPP

namespace El {
namespace rq {

template<typename F>
void ApplyQ
( LeftOrRight side, Orientation orientation, 
  const Matrix<F>& A, const Matrix<F>& t, 
  const Matrix<Base<F>>& d, Matrix<F>& B )
{
    DEBUG_ONLY(CallStackEntry cse("rq::ApplyQ"))
    const bool normal = (orientation==NORMAL);
    const bool onLeft = (side==LEFT);
    const bool applyDFirst = normal!=onLeft;

    const ForwardOrBackward direction = ( normal==onLeft ? BACKWARD : FORWARD );
    const Conjugation conjugation = ( normal ? CONJUGATED : UNCONJUGATED );
    const Int offset = A.Width()-A.Height();
    const Int minDim = Min(A.Height(),A.Width());

    const Int m = B.Height();
    const Int n = B.Width();

    if( applyDFirst )
    {
        if( onLeft )
        {
            auto BBot = View( B, IndexRange(m-minDim,m), IndexRange(0,n) );
            DiagonalScale( side, orientation, d, BBot );
        }
        else
        {
            auto BRight = View( B, IndexRange(0,m), IndexRange(n-minDim,n) );
            DiagonalScale( side, orientation, d, BRight );
        }
    }

    ApplyPackedReflectors
    ( side, LOWER, HORIZONTAL, direction, conjugation, offset, A, t, B );

    if( !applyDFirst )
    {
        if( onLeft )
        {
            auto BBot = View( B, IndexRange(m-minDim,m), IndexRange(0,n) );
            DiagonalScale( side, orientation, d, BBot );
        }
        else
        {
            auto BRight = View( B, IndexRange(0,m), IndexRange(n-minDim,n) );
            DiagonalScale( side, orientation, d, BRight );
        }
    }
}

template<typename F>
void ApplyQ
( LeftOrRight side, Orientation orientation, 
  const AbstractDistMatrix<F>& APre, const AbstractDistMatrix<F>& tPre, 
  const AbstractDistMatrix<Base<F>>& d, AbstractDistMatrix<F>& BPre )
{
    DEBUG_ONLY(CallStackEntry cse("rq::ApplyQ"))
    const bool normal = (orientation==NORMAL);
    const bool onLeft = (side==LEFT);
    const bool applyDFirst = normal!=onLeft;

    const ForwardOrBackward direction = ( normal==onLeft ? BACKWARD : FORWARD );
    const Conjugation conjugation = ( normal ? CONJUGATED : UNCONJUGATED );
    const Int offset = APre.Width()-APre.Height();
    const Int minDim = Min(APre.Height(),APre.Width());

    const Grid& g = APre.Grid();
    DistMatrix<F> A(g), B(g); 
    DistMatrix<F,MD,STAR> t(g);
    Copy( APre, A, READ_PROXY );
    t.SetRoot( A.DiagonalRoot(offset) );
    t.AlignCols( A.DiagonalAlign(offset) );
    Copy( tPre, t, READ_PROXY );
    Copy( BPre, B, READ_WRITE_PROXY );

    const Int m = B.Height();
    const Int n = B.Width();

    if( applyDFirst )
    {
        if( onLeft )
        {
            auto BBot = View( B, IndexRange(m-minDim,m), IndexRange(0,n) );
            DiagonalScale( side, orientation, d, BBot );
        }
        else
        {
            auto BRight = View( B, IndexRange(0,m), IndexRange(n-minDim,n) );
            DiagonalScale( side, orientation, d, BRight );
        }
    }

    ApplyPackedReflectors
    ( side, LOWER, HORIZONTAL, direction, conjugation, offset, A, t, B );

    if( !applyDFirst ) 
    {
        if( onLeft )
        {
            auto BBot = View( B, IndexRange(m-minDim,m), IndexRange(0,n) );
            DiagonalScale( side, orientation, d, BBot );
        }
        else
        {
            auto BRight = View( B, IndexRange(0,m), IndexRange(n-minDim,n) );
            DiagonalScale( side, orientation, d, BRight );
        }
    }
    Copy( B, BPre, RESTORE_READ_WRITE_PROXY );
}

} // namespace rq
} // namespace El

#endif // ifndef EL_RQ_APPLYQ_HPP
