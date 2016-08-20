/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_LDL_PIVOTED_PANEL_HPP
#define EL_LDL_PIVOTED_PANEL_HPP

namespace El {
namespace ldl {
namespace pivot {

template<typename F>
inline LDLPivot
SelectFromPanel
( const Matrix<F>& A,
  const Matrix<F>& X,
  const Matrix<F>& Y, 
  LDLPivotType pivotType,
  Base<F> gamma )
{
    DEBUG_CSE
    LDLPivot pivot;
    switch( pivotType )
    {
    case BUNCH_KAUFMAN_A: 
    case BUNCH_KAUFMAN_C: pivot = PanelBunchKaufmanA( A, X, Y, gamma ); break;
    case BUNCH_KAUFMAN_D: pivot = PanelBunchKaufmanD( A, X, Y, gamma ); break;
    default: LogicError("This pivot type not yet supported");
    }
    return pivot;
}

template<typename F>
inline LDLPivot
SelectFromPanel
( const DistMatrix<F>& A, 
  const DistMatrix<F,MC,STAR>& X, 
  const DistMatrix<F,MR,STAR>& Y, 
  LDLPivotType pivotType, Base<F> gamma )
{
    DEBUG_CSE
    LDLPivot pivot;
    switch( pivotType )
    {
    case BUNCH_KAUFMAN_A: 
    case BUNCH_KAUFMAN_C: pivot = PanelBunchKaufmanA( A, X, Y, gamma ); break;
    case BUNCH_KAUFMAN_D: pivot = PanelBunchKaufmanD( A, X, Y, gamma ); break;
    default: LogicError("This pivot type not yet supported");
    }
    return pivot;
}

// We must use a lazy algorithm so that the symmetric pivoting does not move
// data from a fully-updated to partially-updated region (and vice-versa)
template<typename F>
inline void
Panel
( Matrix<F>& AFull,
  Matrix<F>& dSub,
  Permutation& PFull, 
  Matrix<F>& X,
  Matrix<F>& Y,
  Int bsize,
  Int off=0,
  bool conjugate=false,
  LDLPivotType pivotType=BUNCH_KAUFMAN_A, 
  Base<F> gamma=0 )
{
    DEBUG_CSE
    const Int nFull = AFull.Height();
    auto A = AFull( IR(off,nFull), IR(off,nFull) );
    const Int n = A.Height();
    Zeros( X, n, bsize );
    Zeros( Y, n, bsize );
    if( n == 0 )
        return;
    DEBUG_ONLY(
      if( A.Width() != n )
          LogicError("A must be square");
      if( dSub.Height() != n-1 || dSub.Width() != 1 )
          LogicError("dSub is the wrong size" );
    )

    Int k=0;
    while( k < bsize )
    {
        const IR ind0( 0, k ), indB( k, END ), indR( k, END );

        // Determine the pivot (block)
        auto X0 = X( ALL, ind0 );
        auto Y0 = Y( ALL, ind0 );
        if( pivotType == BUNCH_KAUFMAN_C )
        {
            auto ABR = A( indB, indR );
            // TODO: Form updated diagonal
            LogicError("Have not yet added diagonal update");
            const auto diagMax = VectorMaxAbsLoc( GetDiagonal(ABR) );
            SymmetricSwap
            ( LOWER, AFull, off+k, off+k+diagMax.index, conjugate );
            PFull.Swap( off+k, off+k+diagMax.index );
            RowSwap( X0, k, k+diagMax.index );
            RowSwap( Y0, k, k+diagMax.index );
        }
        const auto pivot = SelectFromPanel( A, X0, Y0, pivotType, gamma );
        const Int from = pivot.from[pivot.nb-1];
        const Int to = k + (pivot.nb-1);
        if( k+pivot.nb > bsize )
        {
            X.Resize( n, bsize-1 );
            Y.Resize( n, bsize-1 );
            break;
        }

        // Apply the symmetric pivot
        SymmetricSwap( LOWER, AFull, off+to, off+from, conjugate );
        PFull.Swap( off+to, off+from );
        RowSwap( X0, to, from );
        RowSwap( Y0, to, from );

        // Update the active columns and then store the new update factors
        const IR ind1( k, k+pivot.nb ), ind2( k+pivot.nb, n );
        if( pivot.nb == 1 ) 
        {
            // Update A(k:end,k) -= X(k:n-1,0:k-1) Y(k,0:k-1)^T
            auto XB0 = X( indB, ind0 ); 
            auto y10 = Y( ind1, ind0 ); 
            auto aB1 = A( indB, ind1 );
            Gemv( NORMAL, F(-1), XB0, y10, F(1), aB1 );
            if( conjugate )
                aB1.MakeReal(0,0);

            // Store x21 := a21/delta11 and y21 := a21
            const F delta11Inv = F(1)/A(k,k);
            auto a21 = A( ind2, ind1 );
            auto x21 = X( ind2, ind1 ); 
            auto y21 = Y( ind2, ind1 ); 
            if( conjugate )
                Conjugate( a21, y21 );
            else
                y21 = a21;
            a21 *= delta11Inv;
            x21 = a21;
        }
        else
        {
            // Update A(k:end,k:k+1) -= X(k:n-1,0:k-1) Y(k:k+1,0:k-1)^T
            // NOTE: top-right entry of AB1 is above-diagonal
            auto XB0 = X( indB, ind0 );
            auto Y10 = Y( ind1, ind0 );
            auto AB1 = A( indB, ind1 );
            const F psi = AB1(0,1);
            Gemm( NORMAL, TRANSPOSE, F(-1), XB0, Y10, F(1), AB1 );
            AB1.Set(0,1,psi);
            if( conjugate )
            {
                AB1.MakeReal(0,0);
                AB1.MakeReal(1,1);
            }

            // Store X21 := A21/D11 and Y21 := A21 or Y21 := Conj(A21)
            auto D11 = A( ind1, ind1 );
            auto A21 = A( ind2, ind1 );
            auto X21 = X( ind2, ind1 );
            auto Y21 = Y( ind2, ind1 );
            if( conjugate )
                Conjugate( A21, Y21 );
            else
                Y21 = A21;
            
            auto D11Inv = D11;
            Symmetric2x2Inv( LOWER, D11Inv, conjugate );
            MakeSymmetric( LOWER, D11Inv, conjugate );
            Transform2x2Cols( D11Inv, A21, 0, 1 );

            X21 = A21;

            // Only leave the main diagonal of D in A, so that routines like
            // Trsm can still be used. Thus, return the subdiagonal.
            dSub.Set( k, 0, D11(1,0) );
            D11.Set( 1, 0, F(0) );
        }
        k += pivot.nb;
    }
}

template<typename F>
inline void
Panel
( DistMatrix<F>& AFull, 
  ElementalMatrix<F>& dSub, 
  DistPermutation& PFull, 
  DistMatrix<F,MC,STAR>& X,
  DistMatrix<F,MR,STAR>& Y,
  Int bsize,
  Int off=0,
  bool conjugate=false,
  LDLPivotType pivotType=BUNCH_KAUFMAN_A,
  Base<F> gamma=0 )
{
    DEBUG_CSE
    const Int nFull = AFull.Height();
    auto A = AFull( IR(off,nFull), IR(off,nFull) );
    const Int n = A.Height();
    X.AlignWith( A );
    Y.AlignWith( A );
    Zeros( X, n, bsize );
    Zeros( Y, n, bsize );

    if( n == 0 )
        return;
    DEBUG_ONLY(
      if( A.Width() != n )
          LogicError("A must be square");
      if( dSub.Height() != n-1 || dSub.Width() != 1 )
          LogicError("dSub is the wrong size" );
    )

    DistMatrix<F,STAR,STAR> D11_STAR_STAR( A.Grid() ),
                            D11Inv_STAR_STAR( A.Grid() );

    Int k=0;
    while( k < bsize )
    {
        const IR ind0( 0, k ), indB( k, END ), indR( k, END );

        // Determine the pivot (block)
        auto X0 = X( ALL, ind0 );
        auto Y0 = Y( ALL, ind0 );
        if( pivotType == BUNCH_KAUFMAN_C )
        {
            auto ABR = A( indB, indR );
            // TODO: Form updated diagonal
            LogicError("Have not yet added diagonal update");
            const auto diagMax = VectorMaxAbsLoc( GetDiagonal(ABR) );
            SymmetricSwap
            ( LOWER, AFull, off+k, off+k+diagMax.index, conjugate );
            PFull.Swap( off+k, off+k+diagMax.index );
            RowSwap( X0, k, k+diagMax.index );
            RowSwap( Y0, k, k+diagMax.index );
        }
        const auto pivot = SelectFromPanel( A, X0, Y0, pivotType, gamma );
        const Int from = pivot.from[pivot.nb-1];
        const Int to = k + (pivot.nb-1);
        if( k+pivot.nb > bsize )
        {
            X.Resize( n, bsize-1 );
            Y.Resize( n, bsize-1 );
            break;
        }

        // Apply the symmetric pivot
        SymmetricSwap( LOWER, AFull, off+to, off+from, conjugate );
        PFull.Swap( off+to, off+from );
        RowSwap( X0, to, from );
        RowSwap( Y0, to, from );

        // Update the active columns and then store the new update factors
        const IR ind1( k, k+pivot.nb ), ind2( k+pivot.nb, n );
        if( pivot.nb == 1 ) 
        {
            // Update A(k:end,k) -= X(k:n-1,0:k-1) Y(k,0:k-1)^T
            auto aB1 = A( indB, ind1 );
            if( aB1.RowAlign() == aB1.RowRank() )
            {
                auto XB0 = X( indB, ind0 );
                auto y10 = Y( ind1, ind0 );
                LocalGemv( NORMAL, F(-1), XB0, y10, F(1), aB1 );
            }
            if( conjugate )
                aB1.MakeReal(0,0);

            // Store x21 := a21/delta11 and y21 := a21
            const F delta11Inv = F(1)/A.Get(k,k);
            auto a21 = A( ind2, ind1 );
            auto x21 = X( ind2, ind1 );
            auto y21 = Y( ind2, ind1 ); 
            if( conjugate )
                Conjugate( a21, y21 );
            else
                y21 = a21;
            a21 *= delta11Inv;
            x21 = a21;
        }
        else
        {
            // Update A(k:end,k:k+1) -= X(k:end,0:k-1) Y(k:k+1,0:k-1)^T
            // NOTE: top-right entry of AB1 is above-diagonal
            auto XB0 = X( indB, ind0 ); 
            auto Y10 = Y( ind1, ind0 ); 
            auto AB1 = A( indB, ind1 );
            // TODO: Make Get and Set local
            const F psi = AB1.Get(0,1);
            LocalGemm( NORMAL, TRANSPOSE, F(-1), XB0, Y10, F(1), AB1 );
            AB1.Set(0,1,psi);
            if( conjugate )
            {
                AB1.MakeReal(0,0);
                AB1.MakeReal(1,1);
            }

            // Store X21 := A21/D11 and Y21 := A21 or Y21 := Conj(A21)
            auto D11 = A( ind1, ind1 );
            auto A21 = A( ind2, ind1 );
            auto X21 = X( ind2, ind1 );
            auto Y21 = Y( ind2, ind1 );
            if( conjugate )
                Conjugate( A21, Y21 );
            else
                Y21 = A21;
            D11_STAR_STAR = D11;

            D11Inv_STAR_STAR = D11_STAR_STAR;
            Symmetric2x2Inv( LOWER, D11Inv_STAR_STAR.Matrix(), conjugate );
            MakeSymmetric( LOWER, D11Inv_STAR_STAR.Matrix(), conjugate );
            Transform2x2Cols( D11Inv_STAR_STAR, A21, 0, 1 );

            X21 = A21;

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

#endif // ifndef EL_LDL_PIVOTED_PANEL_HPP
