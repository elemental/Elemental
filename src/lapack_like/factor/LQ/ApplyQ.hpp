/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_LQ_APPLYQ_HPP
#define EL_LQ_APPLYQ_HPP

namespace El {
namespace lq {

template<typename F>
void ApplyQ
( LeftOrRight side,
  Orientation orientation,
  const Matrix<F>& A,
  const Matrix<F>& householderScalars,
  const Matrix<Base<F>>& signature,
        Matrix<F>& B )
{
    EL_DEBUG_CSE
    const bool normal = (orientation==NORMAL);
    const bool onLeft = (side==LEFT);
    const bool applyDFirst = normal!=onLeft;
    const Int minDim = Min(A.Height(),A.Width());

    const ForwardOrBackward direction = ( normal==onLeft ? FORWARD : BACKWARD );
    const Conjugation conjugation = ( normal ? CONJUGATED : UNCONJUGATED );

    const Int m = B.Height();
    const Int n = B.Width();

    if( applyDFirst )
    {
        if( onLeft )
        {
            auto BTop = B( IR(0,minDim), IR(0,n) );
            DiagonalScale( side, orientation, signature, BTop );
        }
        else
        {
            auto BLeft = B( IR(0,m), IR(0,minDim) );
            DiagonalScale( side, orientation, signature, BLeft );
        }
    }

    ApplyPackedReflectors
    ( side, UPPER, HORIZONTAL, direction, conjugation, 0,
      A, householderScalars, B );

    if( !applyDFirst )
    {
        if( onLeft )
        {
            auto BTop = B( IR(0,minDim), IR(0,n) );
            DiagonalScale( side, orientation, signature, BTop );
        }
        else
        {
            auto BLeft = B( IR(0,m), IR(0,minDim) );
            DiagonalScale( side, orientation, signature, BLeft );
        }
    }
}

template<typename F>
void ApplyQ
( LeftOrRight side,
  Orientation orientation,
  const AbstractDistMatrix<F>& APre,
  const AbstractDistMatrix<F>& householderScalarsPre,
  const AbstractDistMatrix<Base<F>>& signature,
        AbstractDistMatrix<F>& BPre )
{
    EL_DEBUG_CSE
    const bool normal = (orientation==NORMAL);
    const bool onLeft = (side==LEFT);
    const bool applyDFirst = normal!=onLeft;
    const Int minDim = Min(APre.Height(),APre.Width());

    const ForwardOrBackward direction = ( normal==onLeft ? FORWARD : BACKWARD );
    const Conjugation conjugation = ( normal ? CONJUGATED : UNCONJUGATED );

    DistMatrixReadProxy<F,F,MC,MR> AProx( APre );
    DistMatrixReadWriteProxy<F,F,MC,MR> BProx( BPre );
    auto& A = AProx.GetLocked();
    auto& B = BProx.Get();

    ElementalProxyCtrl householderScalarsCtrl;
    householderScalarsCtrl.rootConstrain = true;
    householderScalarsCtrl.colConstrain = true;
    householderScalarsCtrl.root = A.DiagonalRoot();
    householderScalarsCtrl.colAlign = A.DiagonalAlign();

    DistMatrixReadProxy<F,F,MD,STAR>
      householderScalarsProx( householderScalarsPre, householderScalarsCtrl );
    auto& householderScalars = householderScalarsProx.GetLocked();

    const Int m = B.Height();
    const Int n = B.Width();

    if( applyDFirst )
    {
        if( onLeft )
        {
            auto BTop = B( IR(0,minDim), IR(0,n) );
            DiagonalScale( side, orientation, signature, BTop );
        }
        else
        {
            auto BLeft = B( IR(0,m), IR(0,minDim) );
            DiagonalScale( side, orientation, signature, BLeft );
        }
    }

    ApplyPackedReflectors
    ( side, UPPER, HORIZONTAL, direction, conjugation, 0,
      A, householderScalars, B );

    if( !applyDFirst )
    {
        if( onLeft )
        {
            auto BTop = B( IR(0,minDim), IR(0,n) );
            DiagonalScale( side, orientation, signature, BTop );
        }
        else
        {
            auto BLeft = B( IR(0,m), IR(0,minDim) );
            DiagonalScale( side, orientation, signature, BLeft );
        }
    }
}

} // namespace lq
} // namespace El

#endif // ifndef EL_LQ_APPLY_HPP
