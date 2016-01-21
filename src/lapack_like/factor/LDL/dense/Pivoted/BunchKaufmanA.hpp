/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_LDL_PIVOTED_BUNCHKAUFMANA_HPP
#define EL_LDL_PIVOTED_BUNCHKAUFMANA_HPP

// TODO: Reference LAPACK's dsytf2 and zhetf2

namespace El {
namespace ldl {
namespace pivot {

template<typename F>
inline LDLPivot
BunchKaufmanA( const Matrix<F>& A, Base<F> gamma )
{
    DEBUG_ONLY(CSE cse("ldl::pivot::BunchKaufmanA"))
    typedef Base<F> Real;
    const Int n = A.Height();
    if( gamma == Real(0) )
        gamma = LDLPivotConstant<Real>( BUNCH_KAUFMAN_A );

    const Real alpha11Abs = Abs(A.Get(0,0));
    const Range<Int> ind1( 0, 1 ),
                     ind2( 1, n );
    const auto a21Max = VectorMaxAbsLoc( A(ind2,ind1) );
    if( a21Max.value == Real(0) && alpha11Abs == Real(0) )
        throw SingularMatrixException();

    LDLPivot pivot;
    if( alpha11Abs >= gamma*a21Max.value )
    {
        pivot.nb = 1;
        pivot.from[0] = 0;
        return pivot;
    }

    // Find maximum off-diag value in row r (exploit symmetry)
    const Int r = a21Max.index + 1;
    const Range<Int> indr0( 0,   r   ),
                     indr1( r,   r+1 ),
                     indr2( r+1, n   );
    const auto leftMax   = VectorMaxAbsLoc( A(indr1,indr0) );
    const auto bottomMax = VectorMaxAbsLoc( A(indr2,indr1) );
    const Real rowMaxVal = Max(leftMax.value,bottomMax.value);

    if( alpha11Abs >= gamma*a21Max.value*(a21Max.value/rowMaxVal) )
    {
        pivot.nb = 1;
        pivot.from[0] = 0;
        return pivot;
    }

    if( Abs(A.Get(r,r)) >= gamma*rowMaxVal )
    { 
        pivot.nb = 1;
        pivot.from[0] = r;
        return pivot;
    }

    // Default to a 2x2 pivot with 0 and r
    pivot.nb = 2;
    pivot.from[0] = 0;
    pivot.from[1] = r;
    return pivot;
}

template<typename F>
inline LDLPivot
BunchKaufmanA( const DistMatrix<F>& A, Base<F> gamma )
{
    DEBUG_ONLY(CSE cse("ldl::pivot::BunchKaufmanA"))
    typedef Base<F> Real;
    const Int n = A.Height();
    if( gamma == Real(0) )
        gamma = LDLPivotConstant<Real>( BUNCH_KAUFMAN_A );

    const Real alpha11Abs = Abs(A.Get(0,0));
    const Range<Int> ind1( 0, 1 ),
                     ind2( 1, n );
    const auto a21Max = VectorMaxAbsLoc( A(ind2,ind1) );
    if( a21Max.value == Real(0) && alpha11Abs == Real(0) )
        throw SingularMatrixException();

    LDLPivot pivot;
    if( alpha11Abs >= gamma*a21Max.value )
    {
        pivot.nb = 1;
        pivot.from[0] = 0;
        return pivot;
    }

    // Find maximum off-diag value in row r (exploit symmetry)
    const Int r = a21Max.index + 1;
    const Range<Int> indr0( 0,   r   ),
                     indr1( r,   r+1 ),
                     indr2( r+1, n   );
    const auto leftMax   = VectorMaxAbsLoc( A(indr1,indr0) );
    const auto bottomMax = VectorMaxAbsLoc( A(indr2,indr1) );
    const Real rowMaxVal = Max(leftMax.value,bottomMax.value);

    if( alpha11Abs >= gamma*a21Max.value*(a21Max.value/rowMaxVal) )
    {
        pivot.nb = 1;
        pivot.from[0] = 0;
        return pivot;
    }

    if( Abs(A.Get(r,r)) >= gamma*rowMaxVal )
    {
        pivot.nb = 1;
        pivot.from[0] = r;
        return pivot;
    }

    // Default to a 2x2 pivot with 0 and r
    pivot.nb = 2;
    pivot.from[0] = 0;
    pivot.from[1] = r;
    return pivot;
}

// TODO: Switch to the simpler panel update scheme used for Cholesky

template<typename F>
inline LDLPivot
PanelBunchKaufmanA
( const Matrix<F>& A, const Matrix<F>& X, const Matrix<F>& Y, Base<F> gamma )
{
    DEBUG_ONLY(CSE cse("ldl::pivot::PanelBunchKaufmanA"))
    typedef Base<F> Real;
    const Int n = A.Height();
    const Int k = X.Width();
    if( gamma == Real(0) )
        gamma = LDLPivotConstant<Real>( BUNCH_KAUFMAN_A );

    const Range<Int> ind0( 0,   k   ),
                     ind1( k,   k+1 ),  ind1Off( 0, 1   ),
                     ind2( k+1, n   ),  ind2Off( 1, n-k ),
                     indB( k,   n   );

    auto aB1 = A( indB, ind1 );
    auto zB1( aB1 );
    // A(k:n-1,k) -= X(k:n-1,0:k-1) Y(k,0:k-1)^T
    {
        auto XB0 = X( indB, ind0 );
        auto y10 = Y( ind1, ind0 );
        Gemv( NORMAL, F(-1), XB0, y10, F(1), zB1 );
    } 

    const Real alpha11Abs = Abs(zB1.Get(0,0));
    const auto a21Max = VectorMaxAbsLoc( zB1(ind2Off,ind1Off) );
    if( a21Max.value == Real(0) && alpha11Abs == Real(0) )
        throw SingularMatrixException();

    LDLPivot pivot;
    if( alpha11Abs >= gamma*a21Max.value )
    {
        pivot.nb = 1;
        pivot.from[0] = k;
        return pivot;
    }

    // Find maximum off-diag value in row r (exploit symmetry)
    const Int r = a21Max.index + (k+1);
    const Range<Int> indrM( k, r   ),
                     indr1( r, r+1 ), indr1Off( 0, 1   ),
                                      indr2Off( 1, n-r ),
                     indrB( r, n   );
    auto aLeft   = A( indr1, indrM );
    auto aBottom = A( indrB, indr1 );
        
    auto zLeft( aLeft );
    auto zBottom( aBottom );
    auto zStrictBottom = zBottom( indr2Off, indr1Off );

    // Update necessary components out-of-place
    // ----------------------------------------

    // A(r,k:r-1) -= X(r,0:k-1) Y(k:r-1,0:k-1)^T
    {
        auto xr10 = X( indr1, ind0 );
        auto YrM0 = Y( indrM, ind0 );
        Gemv( NORMAL, F(-1), YrM0, xr10, F(1), zLeft );
    }

    // A(r:n-1,r) -= X(r:n-1,0:k-1) Y(r,0:k-1)^T
    {
        auto XrB0 = X( indrB, ind0 );
        auto yr10 = Y( indr1, ind0 );
        Gemv( NORMAL, F(-1), XrB0, yr10, F(1), zBottom );
    } 

    const auto leftMax   = VectorMaxAbsLoc( zLeft );
    const auto bottomMax = VectorMaxAbsLoc( zStrictBottom );
    const Real rowMaxVal = Max(leftMax.value,bottomMax.value);

    if( alpha11Abs >= gamma*a21Max.value*(a21Max.value/rowMaxVal) )
    {
        pivot.nb = 1;
        pivot.from[0] = k;
        return pivot;
    }

    if( Abs(zBottom.Get(0,0)) >= gamma*rowMaxVal )
    {
        pivot.nb = 1;
        pivot.from[0] = r;
        return pivot;
    }

    // Default to a 2x2 pivot with k and r
    pivot.nb = 2;
    pivot.from[0] = k;
    pivot.from[1] = r;
    return pivot;
}

template<typename F>
inline LDLPivot
PanelBunchKaufmanA
( const DistMatrix<F>& A, 
  const DistMatrix<F,MC,STAR>& X, const DistMatrix<F,MR,STAR>& Y, 
  Base<F> gamma )
{
    DEBUG_ONLY(CSE cse("ldl::pivot::PanelBunchKaufmanA"))
    typedef Base<F> Real;
    const Int n = A.Height();
    const Int k = X.Width();
    if( A.ColAlign() != X.ColAlign() || A.RowAlign() != Y.ColAlign() )
        LogicError("X and Y were not properly aligned with A");
    if( gamma == Real(0) )
        gamma = LDLPivotConstant<Real>( BUNCH_KAUFMAN_A );

    const Range<Int> ind0( 0,   k   ),
                     ind1( k,   k+1 ),  ind1Off( 0, 1   ),
                     ind2( k+1, n   ),  ind2Off( 1, n-k ),
                     indB( k,   n   );

    auto aB1 = A( indB, ind1 );
    auto zB1( aB1 );
    // A(k:n-1,k) -= X(k:n-1,0:k-1) Y(k,0:k-1)^T
    if( aB1.RowAlign() == aB1.RowRank() )
    {
        auto XB0 = X( indB, ind0 );
        auto y10 = Y( ind1, ind0 );
        LocalGemv( NORMAL, F(-1), XB0, y10, F(1), zB1 );
    } 

    const Real alpha11Abs = Abs(zB1.Get(0,0));
    const auto a21Max = VectorMaxAbsLoc( zB1(ind2Off,ind1Off) );
    if( a21Max.value == Real(0) && alpha11Abs == Real(0) )
        throw SingularMatrixException();

    LDLPivot pivot;
    if( alpha11Abs >= gamma*a21Max.value )
    {
        pivot.nb = 1;
        pivot.from[0] = k;
        return pivot;
    }

    // Find maximum off-diag value in row r (exploit symmetry)
    const Int r = a21Max.index + (k+1);
    const Range<Int> indrM( k, r   ),
                     indr1( r, r+1 ), indr1Off( 0, 1   ),
                                      indr2Off( 1, n-r ),
                     indrB( r, n   );
    auto aLeft   = A( indr1, indrM );
    auto aBottom = A( indrB, indr1 );
        
    auto zLeft( aLeft );
    auto zBottom( aBottom );
    auto zStrictBottom = zBottom( indr2Off, indr1Off );

    // Update necessary components out-of-place
    // ----------------------------------------

    // A(r,k:r-1) -= X(r,0:k-1) Y(k:r-1,0:k-1)^T
    if( aLeft.ColAlign() == aLeft.ColRank() )
    {
        auto xr10 = X( indr1, ind0 );
        auto YrM0 = Y( indrM, ind0 );
        LocalGemv( NORMAL, F(-1), YrM0, xr10, F(1), zLeft );
    }

    // A(r:n-1,r) -= X(r:n-1,0:k-1) Y(r,0:k-1)^T
    if( aBottom.RowAlign() == aBottom.RowRank() )
    {
        auto XrB0 = X( indrB, ind0 );
        auto yr10 = Y( indr1, ind0 );
        LocalGemv( NORMAL, F(-1), XrB0, yr10, F(1), zBottom );
    } 

    const auto leftMax   = VectorMaxAbsLoc( zLeft );
    const auto bottomMax = VectorMaxAbsLoc( zStrictBottom );
    const Real rowMaxVal = Max(leftMax.value,bottomMax.value);

    if( alpha11Abs >= gamma*a21Max.value*(a21Max.value/rowMaxVal) )
    {
        pivot.nb = 1;
        pivot.from[0] = k;
        return pivot;
    }

    if( Abs(zBottom.Get(0,0)) >= gamma*rowMaxVal )
    {
        pivot.nb = 1;
        pivot.from[0] = r;
        return pivot;
    }

    // Default to a 2x2 pivot with k and r
    pivot.nb = 2;
    pivot.from[0] = k;
    pivot.from[1] = r; 
    return pivot;
}

} // namespace pivot
} // namespace ldl
} // namespace El

#endif // ifndef EL_LDL_PIVOTED_BUNCHKAUFMANA_HPP
