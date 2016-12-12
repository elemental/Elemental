/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_LDL_PIVOTED_BLOCKED_HPP
#define EL_LDL_PIVOTED_BLOCKED_HPP

namespace El {
namespace ldl {
namespace pivot {

template<typename F>
void
Blocked
( Matrix<F>& A,
  Matrix<F>& dSub,
  Permutation& P,
  bool conjugate=false,
  LDLPivotType pivotType=BUNCH_KAUFMAN_A,
  Base<F> gamma=0 )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( A.Height() != A.Width() )
          LogicError("A must be square");
    )
    const Int n = A.Height();

    P.MakeIdentity( n );
    P.ReserveSwaps( n );

    if( n == 0 )
    {
        dSub.Resize( 0, 1 );
        return;
    }
    Zeros( dSub, n-1, 1 );

    Matrix<F> XB1, YB1;
    const Int bsize = Blocksize();
    Int k=0;
    while( k < n )
    {
        const Int nbProp = Min(bsize,n-k);
        const Range<Int> indB( k, n ), indBSub( k, n-1 );
        auto dSubB = dSub( indBSub, ALL );
        Panel( A, dSubB, P, XB1, YB1, nbProp, k, conjugate, pivotType, gamma );
        const Int nb = XB1.Width();

        // Update the bottom-right panel
        const Range<Int> ind2( k+nb, n ),
                         ind1Pan( 0,  nb  ),
                         ind2Pan( nb, n-k );
        auto A22 = A( ind2, ind2 );
        auto X21 = XB1( ind2Pan, ind1Pan );
        auto Y21 = YB1( ind2Pan, ind1Pan );
        Trrk( LOWER, NORMAL, TRANSPOSE, F(-1), X21, Y21, F(1), A22 );

        k += nb;
    }
}

template<typename F>
void
Blocked
( AbstractDistMatrix<F>& APre,
  AbstractDistMatrix<F>& dSubPre,
  DistPermutation& P,
  bool conjugate=false,
  LDLPivotType pivotType=BUNCH_KAUFMAN_A,
  Base<F> gamma=0 )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      AssertSameGrids( APre, dSubPre );
      if( APre.Height() != APre.Width() )
          LogicError("A must be square");
    )
    const Int n = APre.Height();

    P.MakeIdentity( n );
    P.ReserveSwaps( n );

    if( n == 0 )
    {
        dSubPre.Resize( 0, 1 );
        return;
    }
    dSubPre.Resize( n-1, 1 );

    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre );
    DistMatrixWriteProxy<F,F,MC,STAR> dSubProx( dSubPre );
    auto& A = AProx.Get();
    auto& dSub = dSubProx.Get();

    Zero( dSub );

    const Grid& g = APre.Grid();
    DistMatrix<F,MC,STAR> XB1(g);
    DistMatrix<F,MR,STAR> YB1(g);
    const Int bsize = Blocksize();
    Int k=0;
    while( k < n )
    {
        const Int nbProp = Min(bsize,n-k);
        const Range<Int> indB( k, n ), indBSub( k, n-1 );
        auto dSubB = dSub( indBSub, ALL );
        Panel( A, dSubB, P, XB1, YB1, nbProp, k, conjugate, pivotType, gamma );
        const Int nb = XB1.Width();

        // Update the bottom-right panel
        const Range<Int> ind2( k+nb, n ),
                         ind1Pan( 0,  nb  ),
                         ind2Pan( nb, n-k );
        auto A22 = A( ind2, ind2 );
        auto X21 = XB1( ind2Pan, ind1Pan );
        auto Y21 = YB1( ind2Pan, ind1Pan );
        LocalTrrk( LOWER, TRANSPOSE, F(-1), X21, Y21, F(1), A22 );

        k += nb;
    }
}

} // namespace pivot
} // namespace ldl
} // namespace El

#endif // ifndef EL_LDL_PIVOTED_BLOCKED_HPP
