/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_CHOLESKY_LVAR3PIVOTED_HPP
#define EL_CHOLESKY_LVAR3PIVOTED_HPP

namespace El {
namespace cholesky {

namespace pivot {

template<typename F>
inline LDLPivot
Full( const Matrix<F>& A )
{
    DEBUG_ONLY(CSE cse("cholesky::pivot::Full"))
    const auto diagMax = VectorMaxAbsLoc( GetDiagonal(A) );
    LDLPivot pivot;
    pivot.nb = 1;
    pivot.from[0] = diagMax.index;
    return pivot;
}

template<typename F>
inline LDLPivot
Full( const DistMatrix<F>& A )
{
    DEBUG_ONLY(CSE cse("cholesky::pivot::Full"))
    const auto diagMax = VectorMaxAbsLoc( GetDiagonal(A) );
    LDLPivot pivot;
    pivot.nb = 1;
    pivot.from[0] = diagMax.index;
    return pivot;
}

template<typename F>
inline LDLPivot
PanelFull
( const Matrix<F>& A,
        Matrix<F>& d,
  const Matrix<F>& X,
  const Matrix<F>& Y )
{
    DEBUG_ONLY(CSE cse("cholesky::pivot::PanelFull"))

    // Form updated diagonal
    const Int height = d.Height();
    const Int k = X.Width();
          F* dBuf = d.Buffer();
    const F* XBuf = X.LockedBuffer();
    const F* YBuf = Y.LockedBuffer();
    const Int XLDim = X.LDim();
    const Int YLDim = Y.LDim();
    for( Int i=0; i<height; ++i )
    {
        dBuf[i] -= XBuf[i+(k-1)*XLDim]*YBuf[i+(k-1)*YLDim];
        dBuf[i] = RealPart(dBuf[i]);
    }

    // Return maximum from it
    auto diagMax = VectorMaxAbsLoc( d );
    LDLPivot pivot;
    pivot.nb = 1;
    pivot.from[0] = diagMax.index;
    return pivot;
}

template<typename F>
inline LDLPivot
PanelFull
( const DistMatrix<F>& A, 
        DistMatrix<F,MD,STAR>& d,
  const DistMatrix<F,MC,STAR>& X,
  const DistMatrix<F,MR,STAR>& Y )
{
    DEBUG_ONLY(
      CSE cse("cholesky::pivot::PanelFull");
      if( A.ColAlign() != X.ColAlign() || A.RowAlign() != Y.ColAlign() )
          LogicError("A, X, and Y are not properly aligned");
    )

    // Form updated diagonal
    if( d.Participating() )
    {
        const Int dLocalHeight = d.LocalHeight();
        const Int k = X.Width();
        F* dBuf = d.Buffer();
        const F* XBuf = X.LockedBuffer();
        const F* YBuf = Y.LockedBuffer();
        const Int XLDim = X.LDim();
        const Int YLDim = Y.LDim();

        const Int dColShift = d.ColShift();
        const Int dColStride = d.ColStride();
        for( Int iLoc=0; iLoc<dLocalHeight; ++iLoc )
        {
            const Int i = dColShift + iLoc*dColStride;
            // TODO: Avoid repeated calls to LocalRow and use strides
            const Int iLocX = X.LocalRow(i);
            const Int iLocY = Y.LocalRow(i);
            dBuf[iLoc] -= XBuf[iLocX+(k-1)*XLDim]*YBuf[iLocY+(k-1)*YLDim];
            dBuf[iLoc] = RealPart(dBuf[iLoc]);
        }
    }

    // Return maximum from it
    auto diagMax = VectorMaxAbsLoc( d );
    LDLPivot pivot;
    pivot.nb = 1;
    pivot.from[0] = diagMax.index;
    
    return pivot;
}

} // namespace pivot

template<typename F>
inline void
LUnblockedPivoted( Matrix<F>& A, Permutation& P )
{
    DEBUG_ONLY(
      CSE cse("cholesky::LUnblockedPivoted");
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

        auto a21 = A( ind2, ind1 );
        auto A22 = A( ind2, ind2 );
        auto ABR = A( indB, indR );

        // Determine the pivot
        const LDLPivot pivot = pivot::Full( ABR );

        // Apply the pivot
        const Int from = k + pivot.from[0];
        HermitianSwap( LOWER, A, k, from );
        P.RowSwap( k, from );

        // a21 := a21 / sqrt(alpha11)
        const Base<F> delta11 = Sqrt(ABR.GetRealPart(0,0));
        const Base<F> delta11Inv = Base<F>(1)/delta11;
        ABR.Set(0,0,delta11);
        a21 *= delta11Inv;

        // A22 -= a21 a21'
        Her( LOWER, F(-1), a21, A22 );
    }
}

template<typename F>
inline void
LUnblockedPivoted
( AbstractDistMatrix<F>& APre,
  DistPermutation& P )
{
    DEBUG_ONLY(
      CSE cse("cholesky::LUnblockedPivoted");
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

        auto a21 = A( ind2, ind1 );
        auto A22 = A( ind2, ind2 );
        auto ABR = A( indB, indR );

        // Determine the pivot
        const LDLPivot pivot = pivot::Full( ABR );

        // Apply the pivot
        const Int from = k + pivot.from[0];
        HermitianSwap( LOWER, A, k, from );
        P.RowSwap( k, from );

        // a21 := a21 / sqrt(alpha11)
        const Base<F> delta11 = Sqrt(ABR.GetRealPart(0,0));
        const Base<F> delta11Inv = Base<F>(1)/delta11;
        ABR.Set(0,0,delta11);
        a21 *= delta11Inv;

        // A22 -= a21 a21'
        Her( LOWER, F(-1), a21, A22 );
    }
}

// We must use a lazy algorithm so that the symmetric pivoting does not move
// data from a fully-updated to partially-updated region (and vice-versa)
template<typename F>
inline void
LPanelPivoted
( Matrix<F>& AFull,
  Permutation& PFull, 
  Matrix<F>& X,
  Matrix<F>& Y,
  Int bsize,
  Int off )
{
    DEBUG_ONLY(CSE cse("cholesky::LPanelPivoted"))
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

        auto a21 = A( ind2, ind1 );
        auto aB1 = A( indB, ind1 );
        auto ABR = A( indB, indR );
        auto dB = d( indB, ALL );

        auto x21 = X( ind2, ind1 );
        auto XB0 = X( indB, ind0 );

        auto y10 = Y( ind1, ind0 );
        auto y21 = Y( ind2, ind1 );
        auto YB0 = Y( indB, ind0 );

        // Determine the pivot 
        const auto pivot = pivot::PanelFull( ABR, dB, XB0, YB0 );
        const Int from = k + pivot.from[0];

        // Apply the pivot
        HermitianSwap( LOWER, AFull, k+off, from+off );
        PFull.RowSwap( k+off, from+off );
        RowSwap( dB, 0, pivot.from[0] );
        RowSwap( XB0, 0, pivot.from[0] );
        RowSwap( YB0, 0, pivot.from[0] );

        // Update A(k:end,k) -= X(k:end,0:k-1) Y(k,0:k-1)^T
        Gemv( NORMAL, F(-1), XB0, y10, F(1), aB1 );
        aB1.MakeReal(0,0);

        // a21 := a21 / sqrt(alpha11)
        const Base<F> delta11 = Sqrt(A.GetRealPart(k,k));
        const Base<F> delta11Inv = Base<F>(1)/delta11;
        A.SetRealPart(k,k,delta11);
        a21 *= delta11Inv;

        // Store x21 := a21 and y21 := conj(a21)
        Conjugate( a21, y21 );
        x21 = a21;
    }
}

template<typename F>
inline void
LPanelPivoted
( DistMatrix<F>& AFull,
  DistPermutation& PFull,
  DistMatrix<F,MC,STAR>& X,
  DistMatrix<F,MR,STAR>& Y,
  Int bsize,
  Int off )
{
    DEBUG_ONLY(CSE cse("cholesky::LPanelPivoted"))
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

        auto a21 = A( ind2, ind1 );
        auto aB1 = A( indB, ind1 );
        auto ABR = A( indB, indR );
        auto dB = d( indB, ALL );

        auto x21 = X( ind2, ind1 );
        auto XB0 = X( indB, ind0 );

        auto y10 = Y( ind1, ind0 );
        auto y21 = Y( ind2, ind1 );
        auto YB0 = Y( indB, ind0 );

        // Determine the pivot
        const auto pivot = pivot::PanelFull( ABR, dB, XB0, YB0 );
        const Int from = k + pivot.from[0];

        // Apply the pivot
        HermitianSwap( LOWER, AFull, k+off, from+off );
        PFull.RowSwap( k+off, from+off );
        RowSwap( dB, 0, pivot.from[0] );
        RowSwap( XB0, 0, pivot.from[0] );
        RowSwap( YB0, 0, pivot.from[0] );

        // A(k:end,k) -= X(k:end,0:k-1) Y(k,0:k-1)^T
        if( aB1.RowAlign() == aB1.RowRank() )
            LocalGemv( NORMAL, F(-1), XB0, y10, F(1), aB1 );
        aB1.MakeReal(0,0);

        // a21 := a21 / sqrt(alpha11)
        const Base<F> delta11 = Sqrt(A.GetRealPart(k,k));
        const Base<F> delta11Inv = Base<F>(1)/delta11;
        A.SetRealPart(k,k,delta11);
        a21 *= delta11Inv;

        // Store x21 := a21 and y21 := conj(a21)
        Conjugate( a21, y21 );
        x21 = a21;
    }
}

template<typename F>
inline void
LVar3( Matrix<F>& A, Permutation& P )
{
    DEBUG_ONLY(
      CSE cse("cholesky::LVar3");
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
        LPanelPivoted( A, P, XB1, YB1, nb, k );

        // Update the bottom-right panel
        const Range<Int> ind2( k+nb, n ),
                         ind1Pan( 0,  nb  ),
                         ind2Pan( nb, n-k );
        auto A22 = A( ind2, ind2 );
        auto X21 = XB1( ind2Pan, ind1Pan );
        auto Y21 = YB1( ind2Pan, ind1Pan );
        Trrk( LOWER, NORMAL, TRANSPOSE, F(-1), X21, Y21, F(1), A22 );
    }
}

template<typename F>
inline void
LVar3
( AbstractDistMatrix<F>& APre,
  DistPermutation& P )
{
    DEBUG_ONLY(
      CSE cse("cholesky::LVar3");
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
        LPanelPivoted( A, P, XB1, YB1, nb, k );

        // Update the bottom-right panel
        const Range<Int> ind2( k+nb, n ),
                         ind1Pan( 0,  nb  ),
                         ind2Pan( nb, n-k );
        auto A22 = A( ind2, ind2 );
        auto X21 = XB1( ind2Pan, ind1Pan );
        auto Y21 = YB1( ind2Pan, ind1Pan );
        LocalTrrk( LOWER, TRANSPOSE, F(-1), X21, Y21, F(1), A22 );
    }
}

} // namespace cholesky
} // namespace El

#endif // ifndef EL_CHOLESKY_LVAR3PIVOTED_HPP
