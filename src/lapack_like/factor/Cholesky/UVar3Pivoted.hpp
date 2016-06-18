/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_CHOLESKY_UVAR3PIVOTED_HPP
#define EL_CHOLESKY_UVAR3PIVOTED_HPP

#include "./LVar3Pivoted.hpp"

namespace El {
namespace cholesky {

template<typename F>
void UUnblockedPivoted( Matrix<F>& A, Permutation& P )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( A.Height() != A.Width() )
          LogicError("A must be square");
    )
    const Int n = A.Height();

    P.MakeIdentity( n );
    P.ReserveSwaps( n );

    for( Int k=0; k<n; ++k )
    {
        const Range<Int> ind1( k,   k+1 ),
                         ind2( k+1, n   ),
                         indB( k,   n   ),
                         indR( k,   n   );

        auto a12 = A( ind1, ind2 );
        auto A22 = A( ind2, ind2 );
        auto ABR = A( indB, indR );

        // Determine the pivot
        const LDLPivot pivot = pivot::Full( ABR );

        // Apply the pivot
        const Int from = k + pivot.from[0];
        HermitianSwap( UPPER, A, k, from );
        P.RowSwap( k, from );

        // a12 := a12 / sqrt(alpha11)
        const Base<F> delta11 = Sqrt(ABR.GetRealPart(0,0));
        const Base<F> delta11Inv = Base<F>(1)/delta11;
        ABR.Set(0,0,delta11);
        a12 *= delta11Inv;

        // A22 := A22 - a12' a12
        // NOTE: This is silly, but currently necessary
        Conjugate( a12 );
        Her( UPPER, -F(1), a12, A22 );
        Conjugate( a12 );
    }
}

template<typename F>
void UUnblockedPivoted
( AbstractDistMatrix<F>& APre,
  DistPermutation& P )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( APre.Height() != APre.Width() )
          LogicError("A must be square");
    )

    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre );
    auto& A = AProx.Get();

    const Int n = A.Height();

    P.MakeIdentity( n );
    P.ReserveSwaps( n );

    for( Int k=0; k<n; ++k )
    {
        const Range<Int> ind1( k,   k+1 ),
                         ind2( k+1, n   ),
                         indB( k,   n   ),
                         indR( k,   n   );

        auto a12 = A( ind1, ind2 );
        auto A22 = A( ind2, ind2 );
        auto ABR = A( indB, indR );

        // Determine the pivot 
        const LDLPivot pivot = pivot::Full( ABR );

        // Apply the pivot
        const Int from = k + pivot.from[0];
        HermitianSwap( UPPER, A, k, from );
        P.RowSwap( k, from );

        // a12 := a12 / sqrt(alpha11)
        const Base<F> delta11 = Sqrt(ABR.GetRealPart(0,0));
        const Base<F> delta11Inv = Base<F>(1)/delta11;
        ABR.Set(0,0,delta11);
        a12 *= delta11Inv;

        // A22 := A22 - a12' a12
        // NOTE: This is silly, but currently necessary
        Conjugate( a12 );
        Her( UPPER, -F(1), a12, A22 );
        Conjugate( a12 );
    }
}

// We must use a lazy algorithm so that the symmetric pivoting does not move
// data from a fully-updated to partially-updated region (and vice-versa)
template<typename F>
void UPanelPivoted
( Matrix<F>& AFull,
  Permutation& PFull, 
  Matrix<F>& X,
  Matrix<F>& Y,
  Int bsize,
  Int off )
{
    DEBUG_CSE
    auto A = AFull( IR(off,END), IR(off,END) );
    const Int n = A.Height();
    DEBUG_ONLY(
      if( A.Width() != n )
          LogicError("A must be square");
    )
    Zeros( X, n, bsize );
    Zeros( Y, n, bsize );

    auto d = GetDiagonal(A);

    for( Int k=0; k<bsize; ++k )
    {
        const Range<Int> ind0( 0,   k   ),
                         ind1( k,   k+1 ),
                         ind2( k+1, n   ),
                         indB( k,   n   ),
                         indR( k,   n   );

        auto a12 = A( ind1, ind2 );
        auto a1R = A( ind1, indR );
        auto ABR = A( indB, indR );
        auto dB = d( indB, ALL );

        auto x10 = X( ind1, ind0 );
        auto x21 = X( ind2, ind1 );
        auto XB0 = X( indB, ind0 );

        auto y21 = Y( ind2, ind1 );
        auto YB0 = Y( indB, ind0 );

        // Determine the pivot 
        const auto pivot = pivot::PanelFull( ABR, dB, XB0, YB0 );
        const Int from = k + pivot.from[0];

        // Apply the pivot
        HermitianSwap( UPPER, AFull, k+off, from+off );
        PFull.RowSwap( k+off, from+off );
        RowSwap( dB, 0, pivot.from[0] );
        RowSwap( XB0, 0, pivot.from[0] );
        RowSwap( YB0, 0, pivot.from[0] );

        // Update A(k,k:end) -= X(k,0:k-1) Y(k:end,0:k-1)^T
        Gemv( NORMAL, F(-1), YB0, x10, F(1), a1R );
        a1R.MakeReal(0,0);

        // a12 := a12 / sqrt(alpha11)
        const Base<F> delta11 = Sqrt(A.GetRealPart(k,k));
        const Base<F> delta11Inv = Base<F>(1)/delta11;
        A.Set(k,k,delta11);
        a12 *= delta11Inv;

        // Store x21 := a12' and y21 := a12^T
        Adjoint( a12, x21 );
        Transpose( a12, y21 );
    }
}

template<typename F>
void
UPanelPivoted
( DistMatrix<F>& AFull,
  DistPermutation& PFull,
  DistMatrix<F,MC,STAR>& X,
  DistMatrix<F,MR,STAR>& Y,
  Int bsize,
  Int off )
{
    DEBUG_CSE
    auto A = AFull( IR(off,END), IR(off,END) );
    const Int n = A.Height();
    DEBUG_ONLY(
      if( A.Width() != n )
          LogicError("A must be square");
    )
    X.AlignWith( A );
    Y.AlignWith( A );
    Zeros( X, n, bsize );
    Zeros( Y, n, bsize );

    auto d = GetDiagonal(A);

    for( Int k=0; k<bsize; ++k )
    {
        const Range<Int> ind0( 0,   k   ),
                         ind1( k,   k+1 ),
                         ind2( k+1, n   ),
                         indB( k,   n   ),
                         indR( k,   n   );

        auto a12 = A( ind1, ind2 );
        auto a1R = A( ind1, indR );
        auto ABR = A( indB, indR );
        auto dB = d( indB, ALL );

        auto x10 = X( ind1, ind0 );
        auto x21 = X( ind2, ind1 );
        auto XB0 = X( indB, ind0 );

        auto y21 = Y( ind2, ind1 );
        auto YB0 = Y( indB, ind0 );

        // Determine the pivot
        const auto pivot = pivot::PanelFull( ABR, dB, XB0, YB0 );
        const Int from = k + pivot.from[0];

        // Apply the pivot
        HermitianSwap( UPPER, AFull, k+off, from+off );
        PFull.RowSwap( k+off, from+off );
        RowSwap( dB, 0, pivot.from[0] );
        RowSwap( XB0, 0, pivot.from[0] );
        RowSwap( YB0, 0, pivot.from[0] );

        // Update A(k,k:end) -= X(k,0:k-1) Y(k:n-off-1,0:k-1)^T
        if( a1R.ColAlign() == a1R.ColRank() )
            LocalGemv( NORMAL, F(-1), YB0, x10, F(1), a1R );
        a1R.MakeReal(0,0);

        // a12 := a12 / sqrt(alpha11)
        const Base<F> delta11 = Sqrt(A.GetRealPart(k,k));
        const Base<F> delta11Inv = Base<F>(1)/delta11;
        A.Set(k,k,delta11);
        a12 *= delta11Inv;

        // Store x21 := a12' and y21 := a12^T
        Adjoint( a12, x21 );
        Transpose( a12, y21 );
    }
}

template<typename F>
void UVar3( Matrix<F>& A, Permutation& P )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( A.Height() != A.Width() )
          LogicError("A must be square");
    )
    const Int n = A.Height();

    P.MakeIdentity( n );
    P.ReserveSwaps( n );

    Matrix<F> XB1, YB1;
    const Int bsize = Blocksize();
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);
        UPanelPivoted( A, P, XB1, YB1, nb, k );

        // Update the bottom-right panel
        const Range<Int> ind2( k+nb, n ),
                         ind1Pan( 0,  nb  ),
                         ind2Pan( nb, n-k );
        auto A22 = A( ind2, ind2 );
        auto X21 = XB1( ind2Pan, ind1Pan );
        auto Y21 = YB1( ind2Pan, ind1Pan );
        Trrk( UPPER, NORMAL, TRANSPOSE, F(-1), X21, Y21, F(1), A22 );
    }
}

template<typename F>
void UVar3
( AbstractDistMatrix<F>& APre,
  DistPermutation& P )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( APre.Height() != APre.Width() )
          LogicError("A must be square");
    )

    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre );
    auto& A = AProx.Get();

    const Int n = A.Height();

    P.MakeIdentity( n );
    P.ReserveSwaps( n );

    const Grid& g = A.Grid();
    DistMatrix<F,MC,STAR> XB1(g);
    DistMatrix<F,MR,STAR> YB1(g);
    const Int bsize = Blocksize();
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);
        UPanelPivoted( A, P, XB1, YB1, nb, k );

        // Update the bottom-right panel
        const Range<Int> ind2( k+nb, n ),
                         ind1Pan( 0,  nb  ),
                         ind2Pan( nb, n-k );
        auto A22 = A( ind2, ind2 );
        auto X21 = XB1( ind2Pan, ind1Pan );
        auto Y21 = YB1( ind2Pan, ind1Pan );
        LocalTrrk( UPPER, TRANSPOSE, F(-1), X21, Y21, F(1), A22 );
    }
}

} // namespace cholesky
} // namespace El

#endif // ifndef EL_CHOLESKY_UVAR3PIVOTED_HPP
