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

#include EL_ZEROS_INC

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
    const auto a21Max = VectorMaxAbs( LockedViewRange(A,1,0,n,1) );
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
    const auto leftMax   = VectorMaxAbs( LockedViewRange(A,r,  0,r+1,r  ) );
    const auto bottomMax = VectorMaxAbs( LockedViewRange(A,r+1,r,n,  r+1) );
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
    const auto a21Max = VectorMaxAbs( LockedViewRange(A,1,0,n,1) );
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
    const auto leftMax   = VectorMaxAbs( LockedViewRange(A,r,  0,r+1,r  ) );
    const auto bottomMax = VectorMaxAbs( LockedViewRange(A,r+1,r,n,  r+1) );
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
    const auto a21Max = VectorMaxAbs( LockedViewRange(A,1,0,n,1) );
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
    const auto leftMax   = VectorMaxAbs( LockedViewRange(A,r,0,r+1,r  ) );
    const auto bottomMax = VectorMaxAbs( LockedViewRange(A,r,r,n,  r+1) );
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
    const auto a21Max = VectorMaxAbs( LockedViewRange(A,1,0,n,1) );
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
    const auto leftMax   = VectorMaxAbs( LockedViewRange(A,r,0,r+1,r  ) );
    const auto bottomMax = VectorMaxAbs( LockedViewRange(A,r,r,n,  r+1) );
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

    const ValueInt<Real> diagMax = DiagonalMaxAbs( A );
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

    const ValueInt<Real> diagMax = DiagonalMaxAbs( A );
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

    auto aB1 = LockedViewRange( A, k, k, n, k+1 );
    auto zB1( aB1 );
    // A(k:n-1,k) -= X(k:n-1,0:k-1) Y(k,0:k-1)^T
    {
        auto XBL  = LockedViewRange( X, k, 0, n,   k );
        auto yRow = LockedViewRange( Y, k, 0, k+1, k );
        Gemv( NORMAL, F(-1), XBL, yRow, F(1), zB1 );
    } 

    const Real alpha11Abs = Abs(zB1.Get(0,0));
    const auto a21Max = VectorMaxAbs( LockedViewRange(zB1,1,0,n-k,1) );
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
    auto aLeft   = LockedViewRange( A, r, k, r+1, r   );
    auto aBottom = LockedViewRange( A, r, r, n,   r+1 );
        
    auto zLeft( aLeft );
    auto zBottom( aBottom );
    auto zStrictBottom = ViewRange( zBottom, 1, 0, n-r, 1 );

    //
    // Update necessary components out-of-place
    //

    // A(r,k:r-1) -= X(r,0:k-1) Y(k:r-1,0:k-1)^T
    {
        auto xMid = LockedViewRange( X, r, 0, r+1, k );
        auto YBL = LockedViewRange( Y, k, 0, r, k );
        Gemv( NORMAL, F(-1), YBL, xMid, F(1), zLeft );
    }

    // A(r:n-1,r) -= X(r:n-1,0:k-1) Y(r,0:k-1)^T
    {
        auto XBL = LockedViewRange( X, r, 0, n, k );
        auto yRow = LockedViewRange( Y, r, 0, r+1, k );
        Gemv( NORMAL, F(-1), XBL, yRow, F(1), zBottom );
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

    auto aB1 = LockedViewRange( A, k, k, n, k+1 );
    auto zB1( aB1 );
    // A(k:n-1,k) -= X(k:n-1,0:k-1) Y(k,0:k-1)^T
    if( aB1.RowAlign() == aB1.RowRank() )
    {
        auto XBL  = LockedViewRange( X, k, 0, n,   k );
        auto yRow = LockedViewRange( Y, k, 0, k+1, k );
        LocalGemv( NORMAL, F(-1), XBL, yRow, F(1), zB1 );
    } 

    const Real alpha11Abs = Abs(zB1.Get(0,0));
    const auto a21Max = VectorMaxAbs( LockedViewRange(zB1,1,0,n-k,1) );
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
    auto aLeft   = LockedViewRange( A, r, k, r+1, r   );
    auto aBottom = LockedViewRange( A, r, r, n,   r+1 );
        
    auto zLeft( aLeft );
    auto zBottom( aBottom );
    auto zStrictBottom = ViewRange( zBottom, 1, 0, n-r, 1 );

    //
    // Update necessary components out-of-place
    //

    // A(r,k:r-1) -= X(r,0:k-1) Y(k:r-1,0:k-1)^T
    if( aLeft.ColAlign() == aLeft.ColRank() )
    {
        auto xMid = LockedViewRange( X, r, 0, r+1, k );
        auto YBL = LockedViewRange( Y, k, 0, r, k );
        LocalGemv( NORMAL, F(-1), YBL, xMid, F(1), zLeft );
    }

    // A(r:n-1,r) -= X(r:n-1,0:k-1) Y(r,0:k-1)^T
    if( aBottom.RowAlign() == aBottom.RowRank() )
    {
        auto XBL = LockedViewRange( X, r, 0, n, k );
        auto yRow = LockedViewRange( Y, r, 0, r+1, k );
        LocalGemv( NORMAL, F(-1), XBL, yRow, F(1), zBottom );
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

    auto aB1 = LockedViewRange( A, k, k, n, k+1 );
    auto zB1( aB1 );
    // A(k:n-1,k) -= X(k:n-1,0:k-1) Y(k,0:k-1)^T
    {
        auto XBL  = LockedViewRange( X, k, 0, n,   k );
        auto yRow = LockedViewRange( Y, k, 0, k+1, k );
        Gemv( NORMAL, F(-1), XBL, yRow, F(1), zB1 );
    } 

    const Real alpha11Abs = Abs(zB1.Get(0,0));
    const auto a21Max = VectorMaxAbs( LockedViewRange(zB1,1,0,n-k,1) );
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
    auto aLeft   = LockedViewRange( A, r, k, r+1, r   );
    auto aBottom = LockedViewRange( A, r, r, n,   r+1 );
        
    auto zLeft( aLeft );
    auto zBottom( aBottom );

    //
    // Update necessary components out-of-place
    //

    // A(r,k:r-1) -= X(r,0:k-1) Y(k:r-1,0:k-1)^T
    {
        auto xMid = LockedViewRange( X, r, 0, r+1, k );
        auto YBL = LockedViewRange( Y, k, 0, r, k );
        Gemv( NORMAL, F(-1), YBL, xMid, F(1), zLeft );
    }

    // A(r:n-1,r) -= X(r:n-1,0:k-1) Y(r,0:k-1)^T
    {
        auto XBL = LockedViewRange( X, r, 0, n, k );
        auto yRow = LockedViewRange( Y, r, 0, r+1, k );
        Gemv( NORMAL, F(-1), XBL, yRow, F(1), zBottom );
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

    auto aB1 = LockedViewRange( A, k, k, n, k+1 );
    auto zB1( aB1 );
    // A(k:n-1,k) -= X(k:n-1,0:k-1) Y(k,0:k-1)^T
    if( aB1.RowAlign() == aB1.RowRank() )
    {
        auto XBL  = LockedViewRange( X, k, 0, n,   k );
        auto yRow = LockedViewRange( Y, k, 0, k+1, k );
        LocalGemv( NORMAL, F(-1), XBL, yRow, F(1), zB1 );
    } 

    const Real alpha11Abs = Abs(zB1.Get(0,0));
    const auto a21Max = VectorMaxAbs( LockedViewRange(zB1,1,0,n-k,1) );
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
    auto aLeft   = LockedViewRange( A, r, k, r+1, r   );
    auto aBottom = LockedViewRange( A, r, r, n,   r+1 );
        
    auto zLeft( aLeft );
    auto zBottom( aBottom );

    //
    // Update necessary components out-of-place
    //

    // A(r,k:r-1) -= X(r,0:k-1) Y(k:r-1,0:k-1)^T
    if( aLeft.ColAlign() == aLeft.ColRank() )
    {
        auto xMid = LockedViewRange( X, r, 0, r+1, k );
        auto YBL = LockedViewRange( Y, k, 0, r, k );
        LocalGemv( NORMAL, F(-1), YBL, xMid, F(1), zLeft );
    }

    // A(r:n-1,r) -= X(r:n-1,0:k-1) Y(r,0:k-1)^T
    if( aBottom.RowAlign() == aBottom.RowRank() )
    {
        auto XBL = LockedViewRange( X, r, 0, n, k );
        auto yRow = LockedViewRange( Y, r, 0, r+1, k );
        LocalGemv( NORMAL, F(-1), XBL, yRow, F(1), zBottom );
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
( Matrix<F>& A, Matrix<F>& dSub, Matrix<Int>& pPerm, bool conjugate=false,
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
        pPerm.Resize( 0, 1 );
        return;
    }
    Zeros( dSub, n-1, 1 );

    // Initialize the permutation to the identity
    pPerm.Resize( n, 1 );
    for( Int j=0; j<n; ++j )
        pPerm.Set( j, 0, j );
     
    Matrix<F> Y21;

    Int k=0;
    while( k < n )
    {
        // Determine the pivot (block)
        auto ABR = ViewRange( A, k, k, n, n );
        if( pivotType == BUNCH_KAUFMAN_C )
        {
            LogicError("Have not yet generalized pivot storage");
            const auto diagMax = DiagonalMaxAbs( ABR );
            SymmetricSwap( LOWER, A, k, k+diagMax.index, conjugate );
        }
        const LDLPivot pivot = ChoosePivot( ABR, pivotType, gamma );

        // Update trailing submatrix and store pivots
        if( pivot.nb == 1 )
        {
            const Int from = k + pivot.from[0];
            SymmetricSwap( LOWER, A, k, from, conjugate );
            RowSwap( pPerm, k, from );

            // Rank-one update: A22 -= a21 inv(delta11) a21'
            const F delta11Inv = F(1)/ABR.Get(0,0);
            auto a21 = ViewRange( ABR, 1, 0, n-k, 1   );
            auto A22 = ViewRange( ABR, 1, 1, n-k, n-k );
            Syr( LOWER, -delta11Inv, a21, A22, conjugate );
            Scale( delta11Inv, a21 );

            k += 1;
        }
        else
        {
            const Int from0 = k + pivot.from[0];
            const Int from1 = k + pivot.from[1];
            SymmetricSwap( LOWER, A, k,   from0, conjugate );
            SymmetricSwap( LOWER, A, k+1, from1, conjugate );
            RowSwap( pPerm, k+0, from0 );
            RowSwap( pPerm, k+1, from1 );

            // Rank-two update: A22 -= A21 inv(D11) A21'
            auto D11 = ViewRange( ABR, 0, 0, 2,   2   );
            auto A21 = ViewRange( ABR, 2, 0, n-k, 2   );
            auto A22 = ViewRange( ABR, 2, 2, n-k, n-k );
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

template<typename F,Dist UPerm>
inline void
UnblockedPivoted
( DistMatrix<F>& A, 
  DistMatrix<F,MD,STAR>& dSub, 
  DistMatrix<Int,UPerm,STAR>& pPerm, 
  bool conjugate=false, LDLPivotType pivotType=BUNCH_KAUFMAN_A, 
  Base<F> gamma=0 )
{
    DEBUG_ONLY(
        CallStackEntry cse("ldl::UnblockedPivoted");
        if( A.Height() != A.Width() )
            LogicError("A must be square");
        if( A.Grid() != dSub.Grid() || dSub.Grid() != pPerm.Grid() )
            LogicError("A, dSub, and pPerm must share the same grid");
    )
    const Int n = A.Height();
    dSub.SetRoot( A.DiagonalRoot(-1) );
    dSub.AlignCols( A.DiagonalAlign(-1) );
    Zeros( dSub, n-1, 1 );

    // Initialize the permutation to the identity
    pPerm.Resize( n, 1 );
    for( Int iLoc=0; iLoc<pPerm.LocalHeight(); ++iLoc )
        pPerm.SetLocal( iLoc, 0, pPerm.GlobalRow(iLoc) );

    DistMatrix<F> Y21( A.Grid() );
    DistMatrix<F,STAR,STAR> D11_STAR_STAR( A.Grid() );

    Int k=0;
    while( k < n )
    {
        // Determine the pivot (block)
        auto ABR = ViewRange( A, k, k, n, n );
        if( pivotType == BUNCH_KAUFMAN_C )
        {
            LogicError("Have not yet generalized pivot storage");
            const auto diagMax = DiagonalMaxAbs( ABR );
            SymmetricSwap( LOWER, A, k, k+diagMax.index, conjugate );
        }
        const LDLPivot pivot = ChoosePivot( ABR, pivotType, gamma );

        // Update trailing submatrix and store pivots
        if( pivot.nb == 1 )
        {
            const Int from = k + pivot.from[0];
            SymmetricSwap( LOWER, A, k, from, conjugate );
            RowSwap( pPerm, k, from ); 

            // Rank-one update: A22 -= a21 inv(delta11) a21'
            const F delta11Inv = F(1)/ABR.Get(0,0);
            auto a21 = ViewRange( ABR, 1, 0, n-k, 1   );
            auto A22 = ViewRange( ABR, 1, 1, n-k, n-k );
            Syr( LOWER, -delta11Inv, a21, A22, conjugate );
            Scale( delta11Inv, a21 );

            k += 1;
        }
        else
        {
            const Int from0 = k + pivot.from[0];
            const Int from1 = k + pivot.from[1];
            SymmetricSwap( LOWER, A, k,   from0, conjugate );
            SymmetricSwap( LOWER, A, k+1, from1, conjugate );
            RowSwap( pPerm, k+0, from0 );
            RowSwap( pPerm, k+1, from1 );

            // Rank-two update: A22 -= A21 inv(D11) A21'
            auto D11 = ViewRange( ABR, 0, 0, 2,   2   );
            auto A21 = ViewRange( ABR, 2, 0, n-k, 2   );
            auto A22 = ViewRange( ABR, 2, 2, n-k, n-k );
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
}

// We must use a lazy algorithm so that the symmetric pivoting does not move
// data from a fully-updated to partially-updated region (and vice-versa)
template<typename F>
inline void
PanelPivoted
( Matrix<F>& A, Matrix<F>& dSub, Matrix<Int>& pPerm, 
  Matrix<F>& X, Matrix<F>& Y, Int bsize, Int off=0,
  bool conjugate=false, LDLPivotType pivotType=BUNCH_KAUFMAN_A, 
  Base<F> gamma=0 )
{
    DEBUG_ONLY(CallStackEntry cse("ldl::PanelPivoted"))
    const Int n = A.Height();
    if( n == 0 )
        return;
    DEBUG_ONLY(
        if( A.Width() != n )
            LogicError("A must be square");
        if( dSub.Height() != n-1 || dSub.Width() != 1 )
            LogicError("dSub is the wrong size" );
        if( pPerm.Height() != n || pPerm.Width() != 1 )
            LogicError("permutation vector is the wrong size");
    )
    auto ABR = ViewRange( A, off, off, n, n );
    Zeros( X, n-off, bsize );
    Zeros( Y, n-off, bsize );

    Int k=0;
    while( k < bsize )
    {
        // Determine the pivot (block)
        auto X0 = ViewRange( X, 0, 0, n-off, k );
        auto Y0 = ViewRange( Y, 0, 0, n-off, k );
        if( pivotType == BUNCH_KAUFMAN_C )
        {
            LogicError("Have not yet generalized pivot storage");
            // TODO: Form updated diagonal and select maximum
            auto ABRBR = ViewRange( ABR, k, k, n-off, n-off );
            const auto diagMax = DiagonalMaxAbs( ABRBR );
            SymmetricSwap( LOWER, A, off+k, off+k+diagMax.index, conjugate );
            RowSwap( pPerm, k+off, k+off+diagMax.index );
            RowSwap( X0, k, k+diagMax.index );
            RowSwap( Y0, k, k+diagMax.index );
        }
        const auto pivot = ChoosePanelPivot( ABR, X0, Y0, pivotType, gamma );
        const Int from = off + pivot.from[pivot.nb-1];
        const Int to = (off+k) + (pivot.nb-1);
        if( k+pivot.nb > bsize )
        {
            X.Resize( n-off, bsize-1 );
            Y.Resize( n-off, bsize-1 );
            break;
        }

        // Apply the symmetric pivot
        SymmetricSwap( LOWER, A, to, from, conjugate );
        RowSwap( pPerm, to, from );
        RowSwap( X0, to-off, from-off );
        RowSwap( Y0, to-off, from-off );

        // Update the active columns and then store the new update factors
        if( pivot.nb == 1 ) 
        {
            // Update ABR(k:end,k) -= X(k:n-off-1,0:k-1) Y(k,0:k-1)^T
            auto XB0 = LockedViewRange( X,   k, 0, n-off, k   );
            auto y10 = LockedViewRange( Y,   k, 0, k+1,   k   );
            auto aB1 =       ViewRange( ABR, k, k, n-off, k+1 );
            Gemv( NORMAL, F(-1), XB0, y10, F(1), aB1 );
            if( conjugate )
                aB1.MakeReal(0,0);

            // Store x21 := a21/delta11 and y21 := a21
            const F delta11Inv = F(1)/ABR.Get(k,k);
            auto a21 = ViewRange( ABR, k+1, k, n-off, k+1 );
            auto x21 = ViewRange( X,   k+1, k, n-off, k+1 );
            auto y21 = ViewRange( Y,   k+1, k, n-off, k+1 );
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
            // Update ABR(k:end,k:k+1) -= X(k:n-off-1,0:k-1) Y(k:k+1,0:k-1)^T
            // NOTE: top-right entry of AB1 is above-diagonal
            auto XB0 = LockedViewRange( X,   k, 0, n-off, k   );
            auto Y10 = LockedViewRange( Y,   k, 0, k+2,   k   );
            auto AB1 =       ViewRange( ABR, k, k, n-off, k+2 );
            const F psi = AB1.Get(0,1);
            Gemm( NORMAL, TRANSPOSE, F(-1), XB0, Y10, F(1), AB1 );
            AB1.Set(0,1,psi);
            if( conjugate )
            {
                AB1.MakeReal(0,0);
                AB1.MakeReal(1,1);
            }

            // Store X21 := A21/D11 and Y21 := A21 or Y21 := Conj(A21)
            auto D11 = ViewRange( ABR, k,   k, k+2,   k+2 );
            auto A21 = ViewRange( ABR, k+2, k, n-off, k+2 );
            auto X21 = ViewRange( X,   k+2, k, n-off, k+2 );
            auto Y21 = ViewRange( Y,   k+2, k, n-off, k+2 );
            if( conjugate )
                Conjugate( A21, Y21 );
            else
                Y21 = A21;
            Symmetric2x2Solve( RIGHT, LOWER, D11, A21, conjugate );
            X21 = A21;

            // Only leave the main diagonal of D in A, so that routines like
            // Trsm can still be used. Thus, return the subdiagonal.
            dSub.Set( off+k, 0, D11.Get(1,0) );
            D11.Set( 1, 0, 0 );
            k += 2;
        }
    }
}

template<typename F,Dist UPerm>
inline void
PanelPivoted
( DistMatrix<F>& A, 
  DistMatrix<F,MD,STAR>& dSub, 
  DistMatrix<Int,UPerm,STAR>& pPerm, 
  DistMatrix<F,MC,STAR>& X, DistMatrix<F,MR,STAR>& Y, Int bsize, Int off=0,
  bool conjugate=false, LDLPivotType pivotType=BUNCH_KAUFMAN_A,
  Base<F> gamma=0 )
{
    DEBUG_ONLY(CallStackEntry cse("ldl::PanelPivoted"))
    const Int n = A.Height();
    if( n == 0 )
        return;
    DEBUG_ONLY(
        if( A.Width() != n )
            LogicError("A must be square");
        if( dSub.Height() != n-1 || dSub.Width() != 1 )
            LogicError("dSub is the wrong size" );
        if( pPerm.Height() != n || pPerm.Width() != 1 )
            LogicError("permutation vector is the wrong size");
    )
    auto ABR = ViewRange( A, off, off, n, n );
    X.AlignWith( ABR );
    Y.AlignWith( ABR );
    Zeros( X, n-off, bsize );
    Zeros( Y, n-off, bsize );

    DistMatrix<F,STAR,STAR> D11_STAR_STAR( A.Grid() );

    Int k=0;
    while( k < bsize )
    {
        // Determine the pivot (block)
        auto X0 = ViewRange( X, 0, 0, n-off, k );
        auto Y0 = ViewRange( Y, 0, 0, n-off, k );
        if( pivotType == BUNCH_KAUFMAN_C )
        {
            LogicError("Have not yet generalized pivot storage");
            // TODO: Form updated diagonal and select maximum
            auto ABRBR = ViewRange( ABR, k, k, n-off, n-off );
            const auto diagMax = DiagonalMaxAbs( ABRBR );
            SymmetricSwap( LOWER, A, off+k, off+k+diagMax.index, conjugate );
            RowSwap( X0, k, k+diagMax.index );
            RowSwap( Y0, k, k+diagMax.index );
            RowSwap( pPerm, off+k, off+k+diagMax.index );
        }
        const auto pivot = ChoosePanelPivot( ABR, X0, Y0, pivotType, gamma );
        const Int from = off + pivot.from[pivot.nb-1];
        const Int to = (off+k) + (pivot.nb-1);
        if( k+pivot.nb > bsize )
        {
            X.Resize( n-off, bsize-1 );
            Y.Resize( n-off, bsize-1 );
            break;
        }

        // Apply the symmetric pivot
        SymmetricSwap( LOWER, A, to, from, conjugate );
        RowSwap( pPerm, to, from );
        RowSwap( X0, to-off, from-off );
        RowSwap( Y0, to-off, from-off );

        // Update the active columns and then store the new update factors
        if( pivot.nb == 1 ) 
        {
            // Update ABR(k:end,k) -= X(k:n-off-1,0:k-1) Y(k,0:k-1)^T
            auto aB1 = ViewRange( ABR, k, k, n-off, k+1 );
            if( aB1.RowAlign() == aB1.RowRank() )
            {
                auto XB0 = LockedViewRange( X, k, 0, n-off, k );
                auto y10 = LockedViewRange( Y, k, 0, k+1,   k );
                LocalGemv( NORMAL, F(-1), XB0, y10, F(1), aB1 );
            }
            if( conjugate )
                aB1.MakeReal(0,0);

            // Store x21 := a21/delta11 and y21 := a21
            const F delta11Inv = F(1)/ABR.Get(k,k);
            auto a21 = ViewRange( ABR, k+1, k, n-off, k+1 );
            auto x21 = ViewRange( X,   k+1, k, n-off, k+1 );
            auto y21 = ViewRange( Y,   k+1, k, n-off, k+1 );
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
            // Update ABR(k:end,k:k+1) -= X(k:n-off-1,0:k-1) Y(k:k+1,0:k-1)^T
            // NOTE: top-right entry of AB1 is above-diagonal
            auto XB0 = LockedViewRange( X,   k, 0, n-off, k   );
            auto Y10 = LockedViewRange( Y,   k, 0, k+2,   k   );
            auto AB1 =       ViewRange( ABR, k, k, n-off, k+2 );
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
            auto D11 = ViewRange( ABR, k,   k, k+2,   k+2 );
            auto A21 = ViewRange( ABR, k+2, k, n-off, k+2 );
            auto X21 = ViewRange( X,   k+2, k, n-off, k+2 );
            auto Y21 = ViewRange( Y,   k+2, k, n-off, k+2 );
            if( conjugate )
                Conjugate( A21, Y21 );
            else
                Y21 = A21;
            D11_STAR_STAR = D11;
            Symmetric2x2Solve( RIGHT, LOWER, D11_STAR_STAR, A21, conjugate );
            X21 = A21;

            // Only leave the main diagonal of D in A, so that routines like
            // Trsm can still be used. Thus, return the subdiagonal.
            dSub.Set( off+k, 0, D11_STAR_STAR.GetLocal(1,0) );
            D11.Set( 1, 0, 0 );
            k += 2;
        }
    }
}

template<typename F>
inline void
BlockedPivoted
( Matrix<F>& A, Matrix<F>& dSub, Matrix<Int>& pPerm, bool conjugate=false,
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
        pPerm.Resize( 0, 1 );
        return;
    }
    Zeros( dSub, n-1, 1 );

    // Initialize the permutation to the identity
    pPerm.Resize( n, 1 );
    for( Int i=0; i<n; ++i )
        pPerm.Set( i, 0, i );

    Matrix<F> X, Y;
    const Int bsize = Blocksize();
    Int k=0;
    while( k < n )
    {
        const Int nbProp = Min(bsize,n-k);
        PanelPivoted
        ( A, dSub, pPerm, X, Y, nbProp, k, conjugate, pivotType, gamma );
        const Int nb = X.Width();

        // Update the bottom-right panel
        auto X21B  = ViewRange( X, nb,   0,    n-k, nb );
        auto Y21B  = ViewRange( Y, nb,   0,    n-k, nb );
        auto A22BR = ViewRange( A, k+nb, k+nb, n,   n  );
        Trrk( LOWER, NORMAL, TRANSPOSE, F(-1), X21B, Y21B, F(1), A22BR );

        k += nb;
    }
}

template<typename F,Dist UPerm>
inline void
BlockedPivoted
( DistMatrix<F>& A, 
  DistMatrix<F,MD,STAR>& dSub, 
  DistMatrix<Int,UPerm,STAR>& pPerm, 
  bool conjugate=false, LDLPivotType pivotType=BUNCH_KAUFMAN_A, 
  Base<F> gamma=0 )
{
    DEBUG_ONLY(
        CallStackEntry cse("ldl::BlockedPivoted");
        if( A.Height() != A.Width() )
            LogicError("A must be square");
    )
    const Grid& g = A.Grid();
    const Int n = A.Height();
    if( n == 0 )
    {
        dSub.Resize( 0, 1 );
        pPerm.Resize( 0, 1 );
        return;
    }
    dSub.SetRoot( A.DiagonalRoot(-1) );
    dSub.AlignCols( A.DiagonalAlign(-1) );
    Zeros( dSub, n-1, 1 );

    // Initialize the permutation to the identity
    pPerm.Resize( n, 1 );
    for( Int iLoc=0; iLoc<pPerm.LocalHeight(); ++iLoc )
        pPerm.SetLocal( iLoc, 0, pPerm.GlobalRow(iLoc) );

    DistMatrix<F,MC,STAR> X(g);
    DistMatrix<F,MR,STAR> Y(g);
    const Int bsize = Blocksize();
    Int k=0;
    while( k < n )
    {
        const Int nbProp = Min(bsize,n-k);
        PanelPivoted
        ( A, dSub, pPerm, X, Y, nbProp, k, conjugate, pivotType, gamma );
        const Int nb = X.Width();

        // Update the bottom-right panel
        auto X21B  = ViewRange( X, nb,   0,    n-k, nb );
        auto Y21B  = ViewRange( Y, nb,   0,    n-k, nb );
        auto A22BR = ViewRange( A, k+nb, k+nb, n,   n  );
        LocalTrrk( LOWER, TRANSPOSE, F(-1), X21B, Y21B, F(1), A22BR );

        k += nb;
    }
}

template<typename F>
inline void
Pivoted
( Matrix<F>& A, Matrix<F>& dSub, Matrix<Int>& pPerm, bool conjugate=false,
  LDLPivotType pivotType=BUNCH_KAUFMAN_A, Base<F> gamma=0 )
{
    DEBUG_ONLY(CallStackEntry cse("ldl::Pivoted"))
    switch( pivotType )
    {
    case BUNCH_KAUFMAN_A:
    case BUNCH_KAUFMAN_C:
    case BUNCH_KAUFMAN_D:
        BlockedPivoted( A, dSub, pPerm, conjugate, pivotType, gamma );
        break;
    default:
        UnblockedPivoted( A, dSub, pPerm, conjugate, pivotType, gamma );
    }
}

template<typename F,Dist UPerm>
inline void
Pivoted
( DistMatrix<F>& A, 
  DistMatrix<F,MD,STAR>& dSub, 
  DistMatrix<Int,UPerm,STAR>& pPerm, 
  bool conjugate=false, LDLPivotType pivotType=BUNCH_KAUFMAN_A, 
  Base<F> gamma=0 )
{
    DEBUG_ONLY(CallStackEntry cse("ldl::Pivoted"))
    switch( pivotType )
    {
    case BUNCH_KAUFMAN_A:
    case BUNCH_KAUFMAN_C:
    case BUNCH_KAUFMAN_D:
        BlockedPivoted( A, dSub, pPerm, conjugate, pivotType, gamma );
        break;
    default:
        UnblockedPivoted( A, dSub, pPerm, conjugate, pivotType, gamma );
    }
}

} // namespace ldl
} // namespace El

#endif // ifndef EL_LDL_PIVOTED_HPP
