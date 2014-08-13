/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_LDL_PIVOTED_HPP
#define EL_LDL_PIVOTED_HPP

// TODO: Reference LAPACK's dsytf2 and zhetf2

namespace El {
namespace ldl {

namespace pivot {

// TODO: BunchKaufmanC (pivot maximum diagonal entry, then run A)

template<typename F>
inline LDLPivot
BunchKaufmanA( const Matrix<F>& A, Base<F> gamma )
{
    DEBUG_ONLY(CallStackEntry cse("ldl::pivot::BunchKaufmanA"))
    typedef Base<F> Real;
    const Int n = A.Height();
    if( gamma == Real(0) )
        gamma = (1+Sqrt(Real(17)))/8;

    const Real alpha11Abs = Abs(A.Get(0,0));
    const IndexRange ind1( 0, 1 ),
                     ind2( 1, n );
    const auto a21Max = VectorMaxAbs( LockedView(A,ind2,ind1) );
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
    const IndexRange indr0( 0,   r   ),
                     indr1( r,   r+1 ),
                     indr2( r+1, n   );
    const auto leftMax   = VectorMaxAbs( LockedView(A,indr1,indr0) );
    const auto bottomMax = VectorMaxAbs( LockedView(A,indr2,indr1) );
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
    DEBUG_ONLY(CallStackEntry cse("ldl::pivot::BunchKaufmanA"))
    typedef Base<F> Real;
    const Int n = A.Height();
    if( gamma == Real(0) )
        gamma = (1+Sqrt(Real(17)))/8;

    const Real alpha11Abs = Abs(A.Get(0,0));
    const IndexRange ind1( 0, 1 ),
                     ind2( 1, n );
    const auto a21Max = VectorMaxAbs( LockedView(A,ind2,ind1) );
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
    const IndexRange indr0( 0,   r   ),
                     indr1( r,   r+1 ),
                     indr2( r+1, n   );
    const auto leftMax   = VectorMaxAbs( LockedView(A,indr1,indr0) );
    const auto bottomMax = VectorMaxAbs( LockedView(A,indr2,indr1) );
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
BunchKaufmanD( const Matrix<F>& A, Base<F> gamma )
{
    DEBUG_ONLY(CallStackEntry cse("ldl::pivot::BunchKaufmanD"))
    typedef Base<F> Real;
    const Int n = A.Height();
    if( gamma == Real(0) )
        gamma = Real(525)/1000;

    const Real alpha11Abs = Abs(A.Get(0,0));
    const IndexRange ind1( 0, 1 ),
                     ind2( 1, n );
    const auto a21Max = VectorMaxAbs( LockedView(A,ind2,ind1) );
    if( a21Max.value == Real(0) && alpha11Abs == Real(0) )
        throw SingularMatrixException();

    LDLPivot pivot;
    if( alpha11Abs >= gamma*a21Max.value )
    {
        pivot.nb = 1;
        pivot.from[0] = 0;
        return pivot;
    }

    // Find maximum value in row r (exploit symmetry)
    const Int r = a21Max.index + 1;
    const IndexRange indr0( 0,   r   ),
                     indr1( r,   r+1 ),
                     indrB( r,   n   );
    const auto leftMax   = VectorMaxAbs( LockedView(A,indr1,indr0) );
    const auto bottomMax = VectorMaxAbs( LockedView(A,indrB,indr1) );
    const Real rowMaxVal = Max(leftMax.value,bottomMax.value);

    if( alpha11Abs >= gamma*a21Max.value*(a21Max.value/rowMaxVal) )
    {
        pivot.nb = 1;
        pivot.from[0] = 0;
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
BunchKaufmanD( const DistMatrix<F>& A, Base<F> gamma )
{
    DEBUG_ONLY(CallStackEntry cse("ldl::pivot::BunchKaufmanD"))
    typedef Base<F> Real;
    const Int n = A.Height();
    if( gamma == Real(0) )
        gamma = Real(525)/1000;

    const Real alpha11Abs = Abs(A.Get(0,0));
    const IndexRange ind1( 0, 1 ),
                     ind2( 1, n );
    const auto a21Max = VectorMaxAbs( LockedView(A,ind2,ind1) );
    if( a21Max.value == Real(0) && alpha11Abs == Real(0) )
        throw SingularMatrixException();

    LDLPivot pivot;
    if( alpha11Abs >= gamma*a21Max.value )
    {
        pivot.nb = 1;
        pivot.from[0] = 0;
        return pivot;
    }

    // Find maximum value in row r (exploit symmetry)
    const Int r = a21Max.index + 1;
    const IndexRange indr0( 0,   r   ),
                     indr1( r,   r+1 ),
                     indrB( r,   n   );
    const auto leftMax   = VectorMaxAbs( LockedView(A,indr1,indr0) );
    const auto bottomMax = VectorMaxAbs( LockedView(A,indrB,indr1) );
    const Real rowMaxVal = Max(leftMax.value,bottomMax.value);

    if( alpha11Abs >= gamma*a21Max.value*(a21Max.value/rowMaxVal) )
    {
        pivot.nb = 1;
        pivot.from[0] = 0;
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
BunchParlett( const Matrix<F>& A, Base<F> gamma )
{
    DEBUG_ONLY(CallStackEntry cse("ldl::pivot::BunchParlett"))
    typedef Base<F> Real;
    if( gamma == Real(0) )
        gamma = (1+Sqrt(Real(17)))/8;

    const ValueInt<Real> diagMax = VectorMaxAbs( A.GetDiagonal() );
    const ValueIntPair<Real> offDiagMax = SymmetricMaxAbs( LOWER, A );

    LDLPivot pivot;
    if( diagMax.value >= gamma*offDiagMax.value )
    {
        pivot.nb = 1;
        pivot.from[0] = diagMax.index;
        return pivot;
    }
    else
    {
        pivot.nb = 2;
        pivot.from[0] = offDiagMax.indices[0];
        pivot.from[1] = offDiagMax.indices[1];
        return pivot;
    }
}

template<typename F>
inline LDLPivot
BunchParlett( const DistMatrix<F>& A, Base<F> gamma )
{
    DEBUG_ONLY(CallStackEntry cse("ldl::pivot::BunchParlett"))
    typedef Base<F> Real;
    if( gamma == Real(0) )
        gamma = (1+Sqrt(Real(17)))/8;

    const ValueInt<Real> diagMax = VectorMaxAbs( A.GetDiagonal() );
    const ValueIntPair<Real> offDiagMax = SymmetricMaxAbs( LOWER, A );

    LDLPivot pivot;
    if( diagMax.value >= gamma*offDiagMax.value )
    {
        pivot.nb = 1;
        pivot.from[0] = diagMax.index;
        return pivot;
    }
    else
    {
        pivot.nb = 2; 
        pivot.from[0] = offDiagMax.indices[0];
        pivot.from[1] = offDiagMax.indices[1];
        return pivot;
    }
}

// TODO: Switch to the simpler panel update scheme used for Cholesky

template<typename F>
inline LDLPivot
PanelBunchKaufmanA
( const Matrix<F>& A, const Matrix<F>& X, const Matrix<F>& Y, Base<F> gamma )
{
    DEBUG_ONLY(CallStackEntry cse("ldl::pivot::PanelBunchKaufmanA"))
    typedef Base<F> Real;
    const Int n = A.Height();
    const Int k = X.Width();
    if( gamma == Real(0) )
        gamma = (1+Sqrt(Real(17)))/8;

    const IndexRange ind0( 0,   k   ),
                     ind1( k,   k+1 ),  ind1Off( 0, 1   ),
                     ind2( k+1, n   ),  ind2Off( 1, n-k ),
                     indB( k,   n   );

    auto aB1 = LockedView( A, indB, ind1 );
    auto zB1( aB1 );
    // A(k:n-1,k) -= X(k:n-1,0:k-1) Y(k,0:k-1)^T
    {
        auto XB0 = LockedView( X, indB, ind0 );
        auto y10 = LockedView( Y, ind1, ind0 );
        Gemv( NORMAL, F(-1), XB0, y10, F(1), zB1 );
    } 

    const Real alpha11Abs = Abs(zB1.Get(0,0));
    const auto a21Max = VectorMaxAbs( LockedView(zB1,ind2Off,ind1Off) );
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
    const IndexRange indrM( k, r   ),
                     indr1( r, r+1 ), indr1Off( 0, 1   ),
                                      indr2Off( 1, n-r ),
                     indrB( r, n   );
    auto aLeft   = LockedView( A, indr1, indrM );
    auto aBottom = LockedView( A, indrB, indr1 );
        
    auto zLeft( aLeft );
    auto zBottom( aBottom );
    auto zStrictBottom = View( zBottom, indr2Off, indr1Off );

    // Update necessary components out-of-place
    // ----------------------------------------

    // A(r,k:r-1) -= X(r,0:k-1) Y(k:r-1,0:k-1)^T
    {
        auto xr10 = LockedView( X, indr1, ind0 );
        auto YrM0 = LockedView( Y, indrM, ind0 );
        Gemv( NORMAL, F(-1), YrM0, xr10, F(1), zLeft );
    }

    // A(r:n-1,r) -= X(r:n-1,0:k-1) Y(r,0:k-1)^T
    {
        auto XrB0 = LockedView( X, indrB, ind0 );
        auto yr10 = LockedView( Y, indr1, ind0 );
        Gemv( NORMAL, F(-1), XrB0, yr10, F(1), zBottom );
    } 

    const auto leftMax   = VectorMaxAbs( zLeft );
    const auto bottomMax = VectorMaxAbs( zStrictBottom );
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
    DEBUG_ONLY(CallStackEntry cse("ldl::pivot::PanelBunchKaufmanA"))
    typedef Base<F> Real;
    const Int n = A.Height();
    const Int k = X.Width();
    if( A.ColAlign() != X.ColAlign() || A.RowAlign() != Y.ColAlign() )
        LogicError("X and Y were not properly aligned with A");
    if( gamma == Real(0) )
        gamma = (1+Sqrt(Real(17)))/8;

    const IndexRange ind0( 0,   k   ),
                     ind1( k,   k+1 ),  ind1Off( 0, 1   ),
                     ind2( k+1, n   ),  ind2Off( 1, n-k ),
                     indB( k,   n   );

    auto aB1 = LockedView( A, indB, ind1 );
    auto zB1( aB1 );
    // A(k:n-1,k) -= X(k:n-1,0:k-1) Y(k,0:k-1)^T
    if( aB1.RowAlign() == aB1.RowRank() )
    {
        auto XB0 = LockedView( X, indB, ind0 );
        auto y10 = LockedView( Y, ind1, ind0 );
        LocalGemv( NORMAL, F(-1), XB0, y10, F(1), zB1 );
    } 

    const Real alpha11Abs = Abs(zB1.Get(0,0));
    const auto a21Max = VectorMaxAbs( LockedView(zB1,ind2Off,ind1Off) );
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
    const IndexRange indrM( k, r   ),
                     indr1( r, r+1 ), indr1Off( 0, 1   ),
                                      indr2Off( 1, n-r ),
                     indrB( r, n   );
    auto aLeft   = LockedView( A, indr1, indrM );
    auto aBottom = LockedView( A, indrB, indr1 );
        
    auto zLeft( aLeft );
    auto zBottom( aBottom );
    auto zStrictBottom = View( zBottom, indr2Off, indr1Off );

    // Update necessary components out-of-place
    // ----------------------------------------

    // A(r,k:r-1) -= X(r,0:k-1) Y(k:r-1,0:k-1)^T
    if( aLeft.ColAlign() == aLeft.ColRank() )
    {
        auto xr10 = LockedView( X, indr1, ind0 );
        auto YrM0 = LockedView( Y, indrM, ind0 );
        LocalGemv( NORMAL, F(-1), YrM0, xr10, F(1), zLeft );
    }

    // A(r:n-1,r) -= X(r:n-1,0:k-1) Y(r,0:k-1)^T
    if( aBottom.RowAlign() == aBottom.RowRank() )
    {
        auto XrB0 = LockedView( X, indrB, ind0 );
        auto yr10 = LockedView( Y, indr1, ind0 );
        LocalGemv( NORMAL, F(-1), XrB0, yr10, F(1), zBottom );
    } 

    const auto leftMax   = VectorMaxAbs( zLeft );
    const auto bottomMax = VectorMaxAbs( zStrictBottom );
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
PanelBunchKaufmanD
( const Matrix<F>& A, const Matrix<F>& X, const Matrix<F>& Y, Base<F> gamma )
{
    DEBUG_ONLY(CallStackEntry cse("ldl::pivot::PanelBunchKaufmanD"))
    typedef Base<F> Real;
    const Int n = A.Height();
    const Int k = X.Width();
    if( gamma == Real(0) )
        gamma = Real(525)/1000;

    const IndexRange ind0( 0, k   ),
                     ind1( k, k+1 ), ind1Off( 0, 1   ),
                                     ind2Off( 1, n-k ), 
                     indB( k, n   );

    auto aB1 = LockedView( A, indB, ind1 );
    auto zB1( aB1 );
    // A(k:n-1,k) -= X(k:n-1,0:k-1) Y(k,0:k-1)^T
    {
        auto XB0  = LockedView( X, indB, ind0 );
        auto y10 = LockedView( Y, ind1, ind0 );
        Gemv( NORMAL, F(-1), XB0, y10, F(1), zB1 );
    } 

    const Real alpha11Abs = Abs(zB1.Get(0,0));
    const auto a21Max = VectorMaxAbs( LockedView(zB1,ind2Off,ind1Off) );
    if( a21Max.value == Real(0) && alpha11Abs == Real(0) )
        throw SingularMatrixException();

    LDLPivot pivot;
    if( alpha11Abs >= gamma*a21Max.value )
    {
        pivot.nb = 1;
        pivot.from[0] = k;
        return pivot;
    }

    // Find maximum value in row r (exploit symmetry)
    const Int r = a21Max.index + (k+1);
    const IndexRange indrM( k, r   ),
                     indr1( r, r+1 ),
                     indrB( r, n   );
    auto aLeft   = LockedView( A, indr1, indrM );
    auto aBottom = LockedView( A, indrB, indr1 );
        
    auto zLeft( aLeft );
    auto zBottom( aBottom );

    // Update necessary components out-of-place
    // ----------------------------------------

    // A(r,k:r-1) -= X(r,0:k-1) Y(k:r-1,0:k-1)^T
    {
        auto xr10 = LockedView( X, indr1, ind0 );
        auto YrM0 = LockedView( Y, indrM, ind0 );
        Gemv( NORMAL, F(-1), YrM0, xr10, F(1), zLeft );
    }

    // A(r:n-1,r) -= X(r:n-1,0:k-1) Y(r,0:k-1)^T
    {
        auto XrB0 = LockedView( X, indrB, ind0 );
        auto yr10 = LockedView( Y, indr1, ind0 );
        Gemv( NORMAL, F(-1), XrB0, yr10, F(1), zBottom );
    } 

    const auto leftMax   = VectorMaxAbs( zLeft );
    const auto bottomMax = VectorMaxAbs( zBottom );
    const Real rowMaxVal = Max(leftMax.value,bottomMax.value);

    if( alpha11Abs >= gamma*a21Max.value*(a21Max.value/rowMaxVal) )
    {
        pivot.nb = 1;
        pivot.from[0] = k;
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
PanelBunchKaufmanD
( const DistMatrix<F>& A, 
  const DistMatrix<F,MC,STAR>& X, const DistMatrix<F,MR,STAR>& Y, 
  Base<F> gamma )
{
    DEBUG_ONLY(CallStackEntry cse("ldl::pivot::PanelBunchKaufmanD"))
    typedef Base<F> Real;
    const Int n = A.Height();
    const Int k = X.Width();
    if( A.ColAlign() != X.ColAlign() || A.RowAlign() != Y.ColAlign() )
        LogicError("X and Y were not properly aligned with A");
    if( gamma == Real(0) )
        gamma = Real(525)/1000;

    const IndexRange ind0( 0, k   ),
                     ind1( k, k+1 ), ind1Off( 0, 1   ),
                                     ind2Off( 1, n-k ), 
                     indB( k, n   );

    auto aB1 = LockedView( A, indB, ind1 );
    auto zB1( aB1 );
    // A(k:n-1,k) -= X(k:n-1,0:k-1) Y(k,0:k-1)^T
    if( aB1.RowAlign() == aB1.RowRank() )
    {
        auto XB0  = LockedView( X, indB, ind0 );
        auto y10 = LockedView( Y, ind1, ind0 );
        LocalGemv( NORMAL, F(-1), XB0, y10, F(1), zB1 );
    } 

    const Real alpha11Abs = Abs(zB1.Get(0,0));
    const auto a21Max = VectorMaxAbs( LockedView(zB1,ind2Off,ind1Off) );
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
    const IndexRange indrM( k, r   ),
                     indr1( r, r+1 ),
                     indrB( r, n   );
    auto aLeft   = LockedView( A, indr1, indrM );
    auto aBottom = LockedView( A, indrB, indr1 );
        
    auto zLeft( aLeft );
    auto zBottom( aBottom );

    // Update necessary components out-of-place
    // ----------------------------------------

    // A(r,k:r-1) -= X(r,0:k-1) Y(k:r-1,0:k-1)^T
    if( aLeft.ColAlign() == aLeft.ColRank() )
    {
        auto xr10 = LockedView( X, indr1, ind0 );
        auto YrM0 = LockedView( Y, indrM, ind0 );
        LocalGemv( NORMAL, F(-1), YrM0, xr10, F(1), zLeft );
    }

    // A(r:n-1,r) -= X(r:n-1,0:k-1) Y(r,0:k-1)^T
    if( aBottom.RowAlign() == aBottom.RowRank() )
    {
        auto XrB0 = LockedView( X, indrB, ind0 );
        auto yr10 = LockedView( Y, indr1, ind0 );
        LocalGemv( NORMAL, F(-1), XrB0, yr10, F(1), zBottom );
    } 

    const auto leftMax   = VectorMaxAbs( zLeft );
    const auto bottomMax = VectorMaxAbs( zBottom );
    const Real rowMaxVal = Max(leftMax.value,bottomMax.value);

    if( alpha11Abs >= gamma*a21Max.value*(a21Max.value/rowMaxVal) )
    {
        pivot.nb = 1;
        pivot.from[0] = k;
        return pivot;
    }

    // Default to a 2x2 pivot with k and r
    pivot.nb = 2;
    pivot.from[0] = k;
    pivot.from[1] = r; 
    return pivot;
}

} // namespace pivot

template<typename F>
inline LDLPivot
ChoosePivot( const Matrix<F>& A, LDLPivotType pivotType, Base<F> gamma )
{
    DEBUG_ONLY(CallStackEntry cse("ldl::ChoosePivot"))
    LDLPivot pivot;
    switch( pivotType )
    {
    case BUNCH_KAUFMAN_A: 
    case BUNCH_KAUFMAN_C:
        pivot = pivot::BunchKaufmanA( A, gamma ); break;
    case BUNCH_KAUFMAN_D: pivot = pivot::BunchKaufmanD( A, gamma ); break;
    case BUNCH_PARLETT:   pivot = pivot::BunchParlett( A, gamma ); break;
    default: LogicError("This pivot type not yet supported");
    }
    return pivot;
}

template<typename F>
inline LDLPivot
ChoosePivot( const DistMatrix<F>& A, LDLPivotType pivotType, Base<F> gamma )
{
    DEBUG_ONLY(CallStackEntry cse("ldl::ChoosePivot"))
    LDLPivot pivot;
    switch( pivotType )
    {
    case BUNCH_KAUFMAN_A: 
    case BUNCH_KAUFMAN_C:
        pivot = pivot::BunchKaufmanA( A, gamma ); break;
    case BUNCH_KAUFMAN_D: pivot = pivot::BunchKaufmanD( A, gamma ); break;
    case BUNCH_PARLETT:   pivot = pivot::BunchParlett( A, gamma ); break;
    default: LogicError("This pivot type not yet supported");
    }
    return pivot;
}

template<typename F>
inline LDLPivot
ChoosePanelPivot
( const Matrix<F>& A, const Matrix<F>& X, const Matrix<F>& Y, 
  LDLPivotType pivotType, Base<F> gamma )
{
    DEBUG_ONLY(CallStackEntry cse("ldl::ChoosePanelPivot"))
    LDLPivot pivot;
    switch( pivotType )
    {
    case BUNCH_KAUFMAN_A: 
    case BUNCH_KAUFMAN_C:
        pivot = pivot::PanelBunchKaufmanA( A, X, Y, gamma ); 
        break;
    case BUNCH_KAUFMAN_D: 
        pivot = pivot::PanelBunchKaufmanD( A, X, Y, gamma ); 
        break;
    default: 
        LogicError("This pivot type not yet supported");
    }
    return pivot;
}

template<typename F>
inline LDLPivot
ChoosePanelPivot
( const DistMatrix<F>& A, 
  const DistMatrix<F,MC,STAR>& X, 
  const DistMatrix<F,MR,STAR>& Y, 
  LDLPivotType pivotType, Base<F> gamma )
{
    DEBUG_ONLY(CallStackEntry cse("ldl::ChoosePanelPivot"))
    LDLPivot pivot;
    switch( pivotType )
    {
    case BUNCH_KAUFMAN_A: 
    case BUNCH_KAUFMAN_C:
        pivot = pivot::PanelBunchKaufmanA( A, X, Y, gamma ); 
        break;
    case BUNCH_KAUFMAN_D: 
        pivot = pivot::PanelBunchKaufmanD( A, X, Y, gamma ); 
        break;
    default: 
        LogicError("This pivot type not yet supported");
    }
    return pivot;
}

// Unblocked sequential pivoted LDL
template<typename F>
inline void
UnblockedPivoted
( Matrix<F>& A, Matrix<F>& dSub, Matrix<Int>& p, bool conjugate=false,
  LDLPivotType pivotType=BUNCH_KAUFMAN_A, Base<F> gamma=0 )
{
    DEBUG_ONLY(
        CallStackEntry cse("ldl::UnblockedPivoted");
        if( A.Height() != A.Width() )
            LogicError("A must be square");
    )
    const Int n = A.Height();
    if( n == 0 )
    {
        dSub.Resize( 0, 1 );
        p.Resize( 0, 1 );
        return;
    }
    Zeros( dSub, n-1, 1 );

    // Initialize the permutation to the identity
    p.Resize( n, 1 );
    for( Int j=0; j<n; ++j )
        p.Set( j, 0, j );
     
    Matrix<F> Y21;

    Int k=0;
    while( k < n )
    {
        const IndexRange indB( k,   n   ),
                         indR( k,   n   );

        // Determine the pivot (block)
        auto ABR = View( A, indB, indR );
        if( pivotType == BUNCH_KAUFMAN_C )
        {
            LogicError("Have not yet generalized pivot storage");
            const auto diagMax = VectorMaxAbs( ABR.GetDiagonal() );
            SymmetricSwap( LOWER, A, k, k+diagMax.index, conjugate );
        }
        const LDLPivot pivot = ChoosePivot( ABR, pivotType, gamma );

        // Update trailing submatrix and store pivots
        if( pivot.nb == 1 )
        {
            const IndexRange ind1( k,   k+1 ),
                             ind2( k+1, n   );

            const Int from = k + pivot.from[0];
            SymmetricSwap( LOWER, A, k, from, conjugate );
            RowSwap( p, k, from );

            // Rank-one update: A22 -= a21 inv(delta11) a21'
            const F delta11Inv = F(1)/ABR.Get(0,0);
            auto a21 = View( A, ind2, ind1 );
            auto A22 = View( A, ind2, ind2 );
            Syr( LOWER, -delta11Inv, a21, A22, conjugate );
            Scale( delta11Inv, a21 );

            k += 1;
        }
        else
        {
            const IndexRange ind1( k,   k+2 ),
                             ind2( k+2, n   );

            const Int from0 = k + pivot.from[0];
            const Int from1 = k + pivot.from[1];
            SymmetricSwap( LOWER, A, k,   from0, conjugate );
            SymmetricSwap( LOWER, A, k+1, from1, conjugate );
            RowSwap( p, k+0, from0 );
            RowSwap( p, k+1, from1 );

            // Rank-two update: A22 -= A21 inv(D11) A21'
            auto D11 = View( A, ind1, ind1 );
            auto A21 = View( A, ind2, ind1 );
            auto A22 = View( A, ind2, ind2 );
            Y21 = A21;
            Symmetric2x2Solve( RIGHT, LOWER, D11, A21, conjugate );
            Trr2( LOWER, F(-1), A21, Y21, A22, conjugate );

            // Only leave the main diagonal of D in A, so that routines like
            // Trsm can still be used. Thus, return the subdiagonal.
            dSub.Set( k, 0, D11.Get(1,0) );
            D11.Set( 1, 0, 0 );
            k += 2;
        }
    }
}

template<typename F>
inline void
UnblockedPivoted
( AbstractDistMatrix<F>& APre, AbstractDistMatrix<F>& dSub, 
  AbstractDistMatrix<Int>& p, bool conjugate=false, 
  LDLPivotType pivotType=BUNCH_KAUFMAN_A, Base<F> gamma=0 )
{
    DEBUG_ONLY(
        CallStackEntry cse("ldl::UnblockedPivoted");
        if( APre.Height() != APre.Width() )
            LogicError("A must be square");
        AssertSameGrids( APre, dSub, p );
    )
    const Int n = APre.Height();
    const Grid& g = APre.Grid();

    Zeros( dSub, n-1, 1 );
    p.Resize( n, 1 );

    DistMatrix<F> A(g);
    Copy( APre, A, READ_WRITE_PROXY );

    // Initialize the permutation to the identity
    for( Int iLoc=0; iLoc<p.LocalHeight(); ++iLoc )
        p.SetLocal( iLoc, 0, p.GlobalRow(iLoc) );

    DistMatrix<F> Y21(g);
    DistMatrix<F,STAR,STAR> D11_STAR_STAR(g);

    Int k=0;
    while( k < n )
    {
        const IndexRange indB( k, n ),
                         indR( k, n );

        // Determine the pivot (block)
        auto ABR = View( A, indB, indR );
        if( pivotType == BUNCH_KAUFMAN_C )
        {
            LogicError("Have not yet generalized pivot storage");
            const auto diagMax = VectorMaxAbs( ABR.GetDiagonal() );
            SymmetricSwap( LOWER, A, k, k+diagMax.index, conjugate );
        }
        const LDLPivot pivot = ChoosePivot( ABR, pivotType, gamma );

        // Update trailing submatrix and store pivots
        if( pivot.nb == 1 )
        {
            const IndexRange ind1( k,   k+1 ),
                             ind2( k+1, n   );

            const Int from = k + pivot.from[0];
            SymmetricSwap( LOWER, A, k, from, conjugate );
            RowSwap( p, k, from ); 

            // Rank-one update: A22 -= a21 inv(delta11) a21'
            const F delta11Inv = F(1)/ABR.Get(0,0);
            auto a21 = View( A, ind2, ind1 );
            auto A22 = View( A, ind2, ind2 );
            Syr( LOWER, -delta11Inv, a21, A22, conjugate );
            Scale( delta11Inv, a21 );

            k += 1;
        }
        else
        {
            const IndexRange ind1( k,   k+2 ),
                             ind2( k+2, n   );

            const Int from0 = k + pivot.from[0];
            const Int from1 = k + pivot.from[1];
            SymmetricSwap( LOWER, A, k,   from0, conjugate );
            SymmetricSwap( LOWER, A, k+1, from1, conjugate );
            RowSwap( p, k+0, from0 );
            RowSwap( p, k+1, from1 );

            // Rank-two update: A22 -= A21 inv(D11) A21'
            auto D11 = View( A, ind1, ind1 );
            auto A21 = View( A, ind2, ind1 );
            auto A22 = View( A, ind2, ind2 );
            Y21 = A21;
            D11_STAR_STAR = D11;
            Symmetric2x2Solve( RIGHT, LOWER, D11_STAR_STAR, A21, conjugate );
            Trr2( LOWER, F(-1), A21, Y21, A22, conjugate );

            // Only leave the main diagonal of D in A, so that routines like
            // Trsm can still be used. Thus, return the subdiagonal.
            dSub.Set( k, 0, D11_STAR_STAR.GetLocal(1,0) );
            D11.Set( 1, 0, 0 );
            k += 2;
        }
    }
    Copy( A, APre, RESTORE_READ_WRITE_PROXY );
}

// We must use a lazy algorithm so that the symmetric pivoting does not move
// data from a fully-updated to partially-updated region (and vice-versa)
template<typename F>
inline void
PanelPivoted
( Matrix<F>& AFull, Matrix<F>& dSub, Matrix<Int>& p, 
  Matrix<F>& X, Matrix<F>& Y, Int bsize, Int off=0,
  bool conjugate=false, LDLPivotType pivotType=BUNCH_KAUFMAN_A, 
  Base<F> gamma=0 )
{
    DEBUG_ONLY(CallStackEntry cse("ldl::PanelPivoted"))
    const Int nFull = AFull.Height();
    auto A = View( AFull, IndexRange(off,nFull), IndexRange(off,nFull) );
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
        if( p.Height() != n || p.Width() != 1 )
            LogicError("permutation vector is the wrong size");
    )

    const IndexRange outerInd( 0, n );

    Int k=0;
    while( k < bsize )
    {
        const IndexRange ind0( 0, k ),
                         indB( k, n ),
                         indR( k, n );

        // Determine the pivot (block)
        auto X0 = View( X, outerInd, ind0 );
        auto Y0 = View( Y, outerInd, ind0 );
        if( pivotType == BUNCH_KAUFMAN_C )
        {
            LogicError("Have not yet generalized pivot storage");
            // TODO: Form updated diagonal and select maximum
            auto ABR = View( A, indB, indR );
            const auto diagMax = VectorMaxAbs( ABR.GetDiagonal() );
            SymmetricSwap
            ( LOWER, AFull, off+k, off+k+diagMax.index, conjugate );
            RowSwap( p,  k, k+diagMax.index );
            RowSwap( X0, k, k+diagMax.index );
            RowSwap( Y0, k, k+diagMax.index );
        }
        const auto pivot = ChoosePanelPivot( A, X0, Y0, pivotType, gamma );
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
        RowSwap( p,  to, from );
        RowSwap( X0, to, from );
        RowSwap( Y0, to, from );

        // Update the active columns and then store the new update factors
        if( pivot.nb == 1 ) 
        {
            const IndexRange ind1( k,   k+1 ),
                             ind2( k+1, n   );

            // Update A(k:end,k) -= X(k:n-1,0:k-1) Y(k,0:k-1)^T
            auto XB0 = LockedView( X, indB, ind0 ); 
            auto y10 = LockedView( Y, ind1, ind0 ); 
            auto aB1 =       View( A, indB, ind1 );
            Gemv( NORMAL, F(-1), XB0, y10, F(1), aB1 );
            if( conjugate )
                aB1.MakeReal(0,0);

            // Store x21 := a21/delta11 and y21 := a21
            const F delta11Inv = F(1)/A.Get(k,k);
            auto a21 = View( A, ind2, ind1 );
            auto x21 = View( X, ind2, ind1 ); 
            auto y21 = View( Y, ind2, ind1 ); 
            if( conjugate )
                Conjugate( a21, y21 );
            else
                y21 = a21;
            Scale( delta11Inv, a21 );
            x21 = a21;

            k += 1;
        }
        else
        {
            const IndexRange ind1( k,   k+2 ),
                             ind2( k+2, n   );

            // Update A(k:end,k:k+1) -= X(k:n-1,0:k-1) Y(k:k+1,0:k-1)^T
            // NOTE: top-right entry of AB1 is above-diagonal
            auto XB0 = LockedView( X, indB, ind0 );
            auto Y10 = LockedView( Y, ind1, ind0 );
            auto AB1 =       View( A, indB, ind1 );
            const F psi = AB1.Get(0,1);
            Gemm( NORMAL, TRANSPOSE, F(-1), XB0, Y10, F(1), AB1 );
            AB1.Set(0,1,psi);
            if( conjugate )
            {
                AB1.MakeReal(0,0);
                AB1.MakeReal(1,1);
            }

            // Store X21 := A21/D11 and Y21 := A21 or Y21 := Conj(A21)
            auto D11 = View( A, ind1, ind1 );
            auto A21 = View( A, ind2, ind1 );
            auto X21 = View( X, ind2, ind1 );
            auto Y21 = View( Y, ind2, ind1 );
            if( conjugate )
                Conjugate( A21, Y21 );
            else
                Y21 = A21;
            Symmetric2x2Solve( RIGHT, LOWER, D11, A21, conjugate );
            X21 = A21;

            // Only leave the main diagonal of D in A, so that routines like
            // Trsm can still be used. Thus, return the subdiagonal.
            dSub.Set( k, 0, D11.Get(1,0) );
            D11.Set( 1, 0, F(0) );
            k += 2;
        }
    }
}

template<typename F>
inline void
PanelPivoted
( DistMatrix<F>& AFull, 
  AbstractDistMatrix<F>& dSub, 
  AbstractDistMatrix<Int>& p, 
  DistMatrix<F,MC,STAR>& X, DistMatrix<F,MR,STAR>& Y, Int bsize, Int off=0,
  bool conjugate=false, LDLPivotType pivotType=BUNCH_KAUFMAN_A,
  Base<F> gamma=0 )
{
    DEBUG_ONLY(CallStackEntry cse("ldl::PanelPivoted"))
    const Int nFull = AFull.Height();
    auto A = View( AFull, IndexRange(off,nFull), IndexRange(off,nFull) );
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
        if( p.Height() != n || p.Width() != 1 )
            LogicError("permutation vector is the wrong size");
    )

    DistMatrix<F,STAR,STAR> D11_STAR_STAR( A.Grid() );

    const IndexRange outerInd( 0, n );

    Int k=0;
    while( k < bsize )
    {
        const IndexRange ind0( 0, k ),
                         indB( k, n ),
                         indR( k, n );

        // Determine the pivot (block)
        auto X0 = View( X, outerInd, ind0 );
        auto Y0 = View( Y, outerInd, ind0 );
        if( pivotType == BUNCH_KAUFMAN_C )
        {
            LogicError("Have not yet generalized pivot storage");
            // TODO: Form updated diagonal and select maximum
            auto ABR = View( A, indB, indR );
            const auto diagMax = VectorMaxAbs( ABR.GetDiagonal() );
            SymmetricSwap
            ( LOWER, AFull, off+k, off+k+diagMax.index, conjugate );
            RowSwap( X0, k, k+diagMax.index );
            RowSwap( Y0, k, k+diagMax.index );
            RowSwap( p,  k, k+diagMax.index );
        }
        const auto pivot = ChoosePanelPivot( A, X0, Y0, pivotType, gamma );
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
        RowSwap( p,  to, from );
        RowSwap( X0, to, from );
        RowSwap( Y0, to, from );

        // Update the active columns and then store the new update factors
        if( pivot.nb == 1 ) 
        {
            const IndexRange ind1( k,   k+1 ),
                             ind2( k+1, n   );

            // Update A(k:end,k) -= X(k:n-1,0:k-1) Y(k,0:k-1)^T
            auto aB1 = View( A, indB, ind1 );
            if( aB1.RowAlign() == aB1.RowRank() )
            {
                auto XB0 = LockedView( X, indB, ind0 );
                auto y10 = LockedView( Y, ind1, ind0 );
                LocalGemv( NORMAL, F(-1), XB0, y10, F(1), aB1 );
            }
            if( conjugate )
                aB1.MakeReal(0,0);

            // Store x21 := a21/delta11 and y21 := a21
            const F delta11Inv = F(1)/A.Get(k,k);
            auto a21 = View( A, ind2, ind1 );
            auto x21 = View( X, ind2, ind1 );
            auto y21 = View( Y, ind2, ind1 ); 
            if( conjugate )
                Conjugate( a21, y21 );
            else
                y21 = a21;
            Scale( delta11Inv, a21 );
            x21 = a21;

            k += 1;
        }
        else
        {
            const IndexRange ind1( k,   k+2 ),
                             ind2( k+2, n   );

            // Update A(k:end,k:k+1) -= X(k:end,0:k-1) Y(k:k+1,0:k-1)^T
            // NOTE: top-right entry of AB1 is above-diagonal
            auto XB0 = LockedView( X, indB, ind0 ); 
            auto Y10 = LockedView( Y, ind1, ind0 ); 
            auto AB1 =       View( A, indB, ind1 );
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
            auto D11 = View( A, ind1, ind1 );
            auto A21 = View( A, ind2, ind1 );
            auto X21 = View( X, ind2, ind1 );
            auto Y21 = View( Y, ind2, ind1 );
            if( conjugate )
                Conjugate( A21, Y21 );
            else
                Y21 = A21;
            D11_STAR_STAR = D11;
            Symmetric2x2Solve( RIGHT, LOWER, D11_STAR_STAR, A21, conjugate );
            X21 = A21;

            // Only leave the main diagonal of D in A, so that routines like
            // Trsm can still be used. Thus, return the subdiagonal.
            dSub.Set( k, 0, D11_STAR_STAR.GetLocal(1,0) );
            D11.Set( 1, 0, 0 );
            k += 2;
        }
    }
}

template<typename F>
inline void
BlockedPivoted
( Matrix<F>& A, Matrix<F>& dSub, Matrix<Int>& p, bool conjugate=false,
  LDLPivotType pivotType=BUNCH_KAUFMAN_A, Base<F> gamma=0 )
{
    DEBUG_ONLY(
        CallStackEntry cse("ldl::BlockedPivoted");
        if( A.Height() != A.Width() )
            LogicError("A must be square");
    )
    const Int n = A.Height();
    if( n == 0 )
    {
        dSub.Resize( 0, 1 );
        p.Resize( 0, 1 );
        return;
    }
    Zeros( dSub, n-1, 1 );

    // Initialize the permutation to the identity
    p.Resize( n, 1 );
    for( Int i=0; i<n; ++i )
        p.Set( i, 0, i );

    Matrix<F> XB1, YB1;
    const Int bsize = Blocksize();
    Int k=0;
    while( k < n )
    {
        const Int nbProp = Min(bsize,n-k);
        const IndexRange indB( k, n ), indBSub( k, n-1 );
        auto dSubB = View( dSub, indBSub, IndexRange(0,1) );
        auto pB    = View( p,    indB,    IndexRange(0,1) );
        PanelPivoted
        ( A, dSubB, pB, XB1, YB1, nbProp, k, conjugate, pivotType, gamma );
        const Int nb = XB1.Width();

        // Update the bottom-right panel
        const IndexRange ind2( k+nb, n ),
                         ind1Pan( 0,  nb  ),
                         ind2Pan( nb, n-k );
        auto A22 =       View( A,   ind2,    ind2    );
        auto X21 = LockedView( XB1, ind2Pan, ind1Pan );
        auto Y21 = LockedView( YB1, ind2Pan, ind1Pan );
        Trrk( LOWER, NORMAL, TRANSPOSE, F(-1), X21, Y21, F(1), A22 );

        k += nb;
    }
}

template<typename F>
inline void
BlockedPivoted
( AbstractDistMatrix<F>& APre, AbstractDistMatrix<F>& dSubPre,
  AbstractDistMatrix<Int>& pPre, bool conjugate=false, 
  LDLPivotType pivotType=BUNCH_KAUFMAN_A, Base<F> gamma=0 )
{
    DEBUG_ONLY(
        CallStackEntry cse("ldl::BlockedPivoted");
        AssertSameGrids( APre, dSubPre, pPre );
        if( APre.Height() != APre.Width() )
            LogicError("A must be square");
    )
    const Int n = APre.Height();
    pPre.Resize( n, 1 );
    if( n == 0 )
    {
        dSubPre.Resize( 0, 1 );
        return;
    }
    dSubPre.Resize( n-1, 1 );

    const Grid& g = APre.Grid();
    DistMatrix<F> A(g);
    DistMatrix<F,MC,STAR> dSub(g);
    DistMatrix<Int,VC,STAR> p(g);
    Copy( APre,    A,    READ_WRITE_PROXY );
    Copy( dSubPre, dSub, WRITE_PROXY      );
    Copy( pPre,    p,    WRITE_PROXY      );

    Zero( dSub );

    // Initialize the permutation to the identity
    if( p.IsLocalCol(0) )
        for( Int iLoc=0; iLoc<p.LocalHeight(); ++iLoc )
            p.SetLocal( iLoc, 0, p.GlobalRow(iLoc) );

    DistMatrix<F,MC,STAR> XB1(g);
    DistMatrix<F,MR,STAR> YB1(g);
    const Int bsize = Blocksize();
    Int k=0;
    while( k < n )
    {
        const Int nbProp = Min(bsize,n-k);
        const IndexRange indB( k, n ), indBSub( k, n-1 );
        auto dSubB = View( dSub, indBSub, IndexRange(0,1) );
        auto pB    = View( p,    indB,    IndexRange(0,1) );
        PanelPivoted
        ( A, dSubB, pB, XB1, YB1, nbProp, k, conjugate, pivotType, gamma );
        const Int nb = XB1.Width();

        // Update the bottom-right panel
        const IndexRange ind2( k+nb, n ),
                         ind1Pan( 0,  nb  ),
                         ind2Pan( nb, n-k );
        auto A22 =       View( A,   ind2,    ind2    );
        auto X21 = LockedView( XB1, ind2Pan, ind1Pan );
        auto Y21 = LockedView( YB1, ind2Pan, ind1Pan );
        LocalTrrk( LOWER, TRANSPOSE, F(-1), X21, Y21, F(1), A22 );

        k += nb;
    }
    Copy( A,    APre,    RESTORE_READ_WRITE_PROXY );
    Copy( dSub, dSubPre, RESTORE_WRITE_PROXY      );
    Copy( p,    pPre,    RESTORE_WRITE_PROXY      );
}

template<typename F>
inline void
Pivoted
( Matrix<F>& A, Matrix<F>& dSub, Matrix<Int>& p, bool conjugate=false,
  LDLPivotType pivotType=BUNCH_KAUFMAN_A, Base<F> gamma=0 )
{
    DEBUG_ONLY(CallStackEntry cse("ldl::Pivoted"))
    switch( pivotType )
    {
    case BUNCH_KAUFMAN_A:
    case BUNCH_KAUFMAN_C:
    case BUNCH_KAUFMAN_D:
        BlockedPivoted( A, dSub, p, conjugate, pivotType, gamma );
        break;
    default:
        UnblockedPivoted( A, dSub, p, conjugate, pivotType, gamma );
    }
}

template<typename F>
inline void
Pivoted
( AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& dSub, 
  AbstractDistMatrix<Int>& p, bool conjugate=false, 
  LDLPivotType pivotType=BUNCH_KAUFMAN_A, Base<F> gamma=0 )
{
    DEBUG_ONLY(CallStackEntry cse("ldl::Pivoted"))
    switch( pivotType )
    {
    case BUNCH_KAUFMAN_A:
    case BUNCH_KAUFMAN_C:
    case BUNCH_KAUFMAN_D:
        BlockedPivoted( A, dSub, p, conjugate, pivotType, gamma );
        break;
    default:
        UnblockedPivoted( A, dSub, p, conjugate, pivotType, gamma );
    }
}

} // namespace ldl
} // namespace El

#endif // ifndef EL_LDL_PIVOTED_HPP
