/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_RQ_APPLYQ_HPP
#define EL_RQ_APPLYQ_HPP

namespace El {
namespace rq {

template<typename F>
void ApplyQ
( LeftOrRight side,
  Orientation orientation, 
  const Matrix<F>& A,
  const Matrix<F>& t, 
  const Matrix<Base<F>>& d,
        Matrix<F>& B )
{
    DEBUG_ONLY(CSE cse("rq::ApplyQ"))
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
            auto BBot = B( IR(m-minDim,m), IR(0,n) );
            DiagonalScale( side, orientation, d, BBot );
        }
        else
        {
            auto BRight = B( IR(0,m), IR(n-minDim,n) );
            DiagonalScale( side, orientation, d, BRight );
        }
    }

    ApplyPackedReflectors
    ( side, LOWER, HORIZONTAL, direction, conjugation, offset, A, t, B );

    if( !applyDFirst )
    {
        if( onLeft )
        {
            auto BBot = B( IR(m-minDim,m), IR(0,n) );
            DiagonalScale( side, orientation, d, BBot );
        }
        else
        {
            auto BRight = B( IR(0,m), IR(n-minDim,n) );
            DiagonalScale( side, orientation, d, BRight );
        }
    }
}

template<typename F>
void ApplyQ
( LeftOrRight side,
  Orientation orientation, 
  const ElementalMatrix<F>& APre,
  const ElementalMatrix<F>& tPre, 
  const ElementalMatrix<Base<F>>& d,
        ElementalMatrix<F>& BPre )
{
    DEBUG_ONLY(CSE cse("rq::ApplyQ"))
    const bool normal = (orientation==NORMAL);
    const bool onLeft = (side==LEFT);
    const bool applyDFirst = normal!=onLeft;

    const ForwardOrBackward direction = ( normal==onLeft ? BACKWARD : FORWARD );
    const Conjugation conjugation = ( normal ? CONJUGATED : UNCONJUGATED );
    const Int offset = APre.Width()-APre.Height();
    const Int minDim = Min(APre.Height(),APre.Width());

    DistMatrixReadProxy<F,F,MC,MR> AProx( APre );
    DistMatrixReadWriteProxy<F,F,MC,MR> BProx( BPre );
    auto& A = AProx.GetLocked();
    auto& B = BProx.Get();

    ElementalProxyCtrl tCtrl;
    tCtrl.rootConstrain = true;
    tCtrl.colConstrain = true;
    tCtrl.root = A.DiagonalRoot(offset);
    tCtrl.colAlign = A.DiagonalAlign(offset);

    DistMatrixReadProxy<F,F,MD,STAR> tProx( tPre, tCtrl );
    auto& t = tProx.GetLocked();

    const Int m = B.Height();
    const Int n = B.Width();

    if( applyDFirst )
    {
        if( onLeft )
        {
            auto BBot = B( IR(m-minDim,m), IR(0,n) );
            DiagonalScale( side, orientation, d, BBot );
        }
        else
        {
            auto BRight = B( IR(0,m), IR(n-minDim,n) );
            DiagonalScale( side, orientation, d, BRight );
        }
    }

    ApplyPackedReflectors
    ( side, LOWER, HORIZONTAL, direction, conjugation, offset, A, t, B );

    if( !applyDFirst ) 
    {
        if( onLeft )
        {
            auto BBot = B( IR(m-minDim,m), IR(0,n) );
            DiagonalScale( side, orientation, d, BBot );
        }
        else
        {
            auto BRight = B( IR(0,m), IR(n-minDim,n) );
            DiagonalScale( side, orientation, d, BRight );
        }
    }
}

} // namespace rq
} // namespace El

#endif // ifndef EL_RQ_APPLYQ_HPP
