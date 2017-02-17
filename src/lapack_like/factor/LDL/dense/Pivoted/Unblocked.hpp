/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_LDL_PIVOTED_UNBLOCKED_HPP
#define EL_LDL_PIVOTED_UNBLOCKED_HPP

namespace El {
namespace ldl {
namespace pivot {

template<typename F>
LDLPivot
Select( const Matrix<F>& A, LDLPivotType pivotType, Base<F> gamma )
{
    EL_DEBUG_CSE
    LDLPivot pivot;
    switch( pivotType )
    {
    case BUNCH_KAUFMAN_A:
    case BUNCH_KAUFMAN_C: pivot = BunchKaufmanA( A, gamma ); break;
    case BUNCH_KAUFMAN_D: pivot = BunchKaufmanD( A, gamma ); break;
    case BUNCH_PARLETT:   pivot = BunchParlett( A, gamma ); break;
    default: LogicError("This pivot type not yet supported");
    }
    return pivot;
}

template<typename F>
LDLPivot
Select( const DistMatrix<F>& A, LDLPivotType pivotType, Base<F> gamma )
{
    EL_DEBUG_CSE
    LDLPivot pivot;
    switch( pivotType )
    {
    case BUNCH_KAUFMAN_A:
    case BUNCH_KAUFMAN_C: pivot = BunchKaufmanA( A, gamma ); break;
    case BUNCH_KAUFMAN_D: pivot = BunchKaufmanD( A, gamma ); break;
    case BUNCH_PARLETT:   pivot = BunchParlett( A, gamma ); break;
    default: LogicError("This pivot type not yet supported");
    }
    return pivot;
}

// Unblocked sequential pivoted LDL
template<typename F>
void
Unblocked
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

    Matrix<F> Y21;

    Int k=0;
    while( k < n )
    {
        const Range<Int> indB( k, n ), indR( k, n );

        // Determine the pivot (block)
        auto ABR = A( indB, indR );
        if( pivotType == BUNCH_KAUFMAN_C )
        {
            LogicError("Have not yet generalized pivot storage");
            const auto diagMax = VectorMaxAbsLoc( GetDiagonal(ABR) );
            SymmetricSwap( LOWER, A, k, k+diagMax.index, conjugate );
        }
        const LDLPivot pivot = Select( ABR, pivotType, gamma );

        for( Int l=0; l<pivot.nb; ++l )
        {
            const Int from = k + pivot.from[l];
            SymmetricSwap( LOWER, A, k+l, from, conjugate );
            P.Swap( k+l, from );
        }

        // Update trailing submatrix and store pivots
        const Range<Int> ind1( k,          k+pivot.nb ),
                         ind2( k+pivot.nb, n          );
        if( pivot.nb == 1 )
        {
            // Rank-one update: A22 -= a21 inv(delta11) a21'
            const F delta11Inv = F(1)/ABR(0,0);
            auto a21 = A( ind2, ind1 );
            auto A22 = A( ind2, ind2 );
            Syr( LOWER, -delta11Inv, a21, A22, conjugate );
            a21 *= delta11Inv;
        }
        else
        {
            // Rank-two update: A22 -= A21 inv(D11) A21'
            auto D11 = A( ind1, ind1 );
            auto A21 = A( ind2, ind1 );
            auto A22 = A( ind2, ind2 );
            Y21 = A21;

            auto D11Inv = D11;
            Symmetric2x2Inv( LOWER, D11Inv, conjugate );
            MakeSymmetric( LOWER, D11Inv, conjugate );
            Transform2x2Cols( D11Inv, A21, 0, 1 );
            Trr2( LOWER, F(-1), A21, Y21, A22, conjugate );

            // Only leave the main diagonal of D in A, so that routines like
            // Trsm can still be used. Thus, return the subdiagonal.
            dSub(k) = D11(1,0);
            D11(1,0) = 0;
        }
        k += pivot.nb;
    }
}

template<typename F>
void
Unblocked
( AbstractDistMatrix<F>& APre,
  AbstractDistMatrix<F>& dSub,
  DistPermutation& P,
  bool conjugate=false,
  LDLPivotType pivotType=BUNCH_KAUFMAN_A,
  Base<F> gamma=0 )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( APre.Height() != APre.Width() )
          LogicError("A must be square");
      AssertSameGrids( APre, dSub );
    )
    const Int n = APre.Height();
    const Grid& g = APre.Grid();

    P.MakeIdentity( n );
    P.ReserveSwaps( n );

    Zeros( dSub, n-1, 1 );

    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre );
    auto& A = AProx.Get();

    DistMatrix<F> Y21(g);
    DistMatrix<F,STAR,STAR> D11_STAR_STAR(g), D11Inv_STAR_STAR(g);

    Int k=0;
    while( k < n )
    {
        const Range<Int> indB( k, n ), indR( k, n );

        // Determine the pivot (block)
        auto ABR = A( indB, indR );
        if( pivotType == BUNCH_KAUFMAN_C )
        {
            LogicError("Have not yet generalized pivot storage");
            const auto diagMax = VectorMaxAbsLoc( GetDiagonal(ABR) );
            SymmetricSwap( LOWER, A, k, k+diagMax.index, conjugate );
        }
        const LDLPivot pivot = Select( ABR, pivotType, gamma );

        for( Int l=0; l<pivot.nb; ++l )
        {
            const Int from = k + pivot.from[l];
            SymmetricSwap( LOWER, A, k+l, from, conjugate );
            P.Swap( k+l, from );
        }


        // Update trailing submatrix and store pivots
        const Range<Int> ind1( k,          k+pivot.nb ),
                         ind2( k+pivot.nb, n          );
        if( pivot.nb == 1 )
        {
            // Rank-one update: A22 -= a21 inv(delta11) a21'
            const F delta11Inv = F(1)/ABR.Get(0,0);
            auto a21 = A( ind2, ind1 );
            auto A22 = A( ind2, ind2 );
            Syr( LOWER, -delta11Inv, a21, A22, conjugate );
            a21 *= delta11Inv;
        }
        else
        {
            // Rank-two update: A22 -= A21 inv(D11) A21'
            auto D11 = A( ind1, ind1 );
            auto A21 = A( ind2, ind1 );
            auto A22 = A( ind2, ind2 );
            Y21 = A21;
            D11_STAR_STAR = D11;

            D11Inv_STAR_STAR = D11_STAR_STAR;
            Symmetric2x2Inv( LOWER, D11Inv_STAR_STAR.Matrix(), conjugate );
            MakeSymmetric( LOWER, D11Inv_STAR_STAR.Matrix(), conjugate );
            Transform2x2Cols( D11Inv_STAR_STAR, A21, 0, 1 );

            Trr2( LOWER, F(-1), A21, Y21, A22, conjugate );

            // Only leave the main diagonal of D in A, so that routines like
            // Trsm can still be used. Thus, return the subdiagonal.
            dSub.Set( k, 0, D11_STAR_STAR.GetLocal(1,0) );
            D11.Set( 1, 0, 0 );
        }
        k += pivot.nb;
    }
}

} // namespace pivot
} // namespace ldl
} // namespace El

#endif // ifndef EL_LDL_PIVOTED_UNBLOCKED_HPP
