/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_DIAGONALSCALETRAPEZOID_HPP
#define ELEM_DIAGONALSCALETRAPEZOID_HPP

#include ELEM_SCALE_INC

namespace elem {

template<typename T,typename TDiag>
inline void
DiagonalScaleTrapezoid
( LeftOrRight side, UpperOrLower uplo, Orientation orientation, 
  const Matrix<TDiag>& d, Matrix<T>& A, Int offset=0 )
{
    DEBUG_ONLY(
        CallStackEntry cse("DiagonalScaleTrapezoid");
        if( side==LEFT && (d.Height()!=A.Height() || d.Width()!=1) )
            LogicError("d should have been a vector of the height of A");
        if( side==RIGHT && (d.Height()!=A.Width() || d.Width()!=1) )
            LogicError("d should have been a vector of the width of A");
    )
    const Int m = A.Height();
    const Int n = A.Width();
    const Int diagLength = A.DiagonalLength(offset);
    const Int ldim = A.LDim();
    T* ABuf = A.Buffer();
    const bool conjugate = ( orientation==ADJOINT );

    const Int iOff = ( offset>=0 ? 0      : -offset );
    const Int jOff = ( offset>=0 ? offset : 0       );

    if( uplo == LOWER && side == LEFT )
    {
        // Scale from the left up to the diagonal
        for( Int i=iOff; i<m; ++i )
        {
            const Int k = i-iOff;
            const Int j = k+jOff;
            const TDiag alpha = ( conjugate ? Conj(d.Get(i,0)) : d.Get(i,0) );
            blas::Scal( Min(j+1,n), alpha, &ABuf[i], ldim );
        }
    }
    else if( uplo == UPPER && side == LEFT )
    {
        // Scale from the diagonal to the right
        for( Int i=0; i<iOff+diagLength; ++i )
        {
            const Int k = i-iOff;
            const Int j = k+jOff;
            const Int jLeft = Max(j,0);
            const TDiag alpha = ( conjugate ? Conj(d.Get(i,0)) : d.Get(i,0) );
            blas::Scal( n-jLeft, alpha, &ABuf[i+jLeft*ldim], ldim );
        }
    }
    else if( uplo == LOWER && side == RIGHT )
    {
        // Scale from the diagonal downwards
        for( Int j=0; j<jOff+diagLength; ++j )
        {
            const Int k = j-jOff;
            const Int i = k+iOff;
            const Int iTop = Max(i,0);
            const TDiag alpha = ( conjugate ? Conj(d.Get(j,0)) : d.Get(j,0) );
            blas::Scal( m-iTop, alpha, &ABuf[iTop+j*ldim], 1 );
        }
    }
    else /* uplo == UPPER && side == RIGHT */
    {
        // Scale downward to the diagonal
        for( Int j=jOff; j<n; ++j )
        {
            const Int k = j-jOff;
            const Int i = k+iOff;
            const TDiag alpha = ( conjugate ? Conj(d.Get(j,0)) : d.Get(j,0) );
            blas::Scal( Min(i+1,m), alpha, &ABuf[j*ldim], 1 );
        }
    }
}

template<typename T,typename TDiag,Dist U,Dist V,Dist W,Dist Z>
inline void
DiagonalScaleTrapezoid
( LeftOrRight side, UpperOrLower uplo, Orientation orientation,
  const DistMatrix<TDiag,U,V>& d, DistMatrix<T,W,Z>& A, Int offset=0 )
{
    DEBUG_ONLY(CallStackEntry cse("DiagonalScaleTrapezoid"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int mLoc = A.LocalHeight();
    const Int nLoc = A.LocalWidth();
    const bool conjugate = ( orientation==ADJOINT );

    const Int diagLength = A.DiagonalLength(offset);
    const Int ldim = A.LDim();
    T* ABuf = A.Buffer();

    const Int iOff = ( offset>=0 ? 0      : -offset );
    const Int jOff = ( offset>=0 ? offset : 0       );

    if( side == LEFT )
    {
        DistMatrix<TDiag,W,STAR> d_W_STAR( A.Grid() );
        if( U == W && V == STAR && d.ColAlign() == A.ColAlign() )
        {
            d_W_STAR = LockedView( d );
        }
        else
        {
            d_W_STAR.AlignWith( A );
            d_W_STAR = d;
        }

        if( uplo == LOWER )
        {
            // Scale from the left up to the diagonal
            for( Int iLoc=0; iLoc<mLoc; ++iLoc )            
            {
                const Int i = A.GlobalRow(iLoc);
                if( i >= iOff )
                {
                    const Int k = i-iOff;
                    const Int j = k+jOff;
                    const Int width = Min(j+1,n);
                    const Int localWidth = A.LocalColOffset(width);
                    const TDiag alpha = 
                        ( conjugate ? Conj(d_W_STAR.GetLocal(iLoc,0))
                                    : d_W_STAR.GetLocal(iLoc,0) );
                    blas::Scal( localWidth, alpha, &ABuf[iLoc], ldim );
                }
            }
        }
        else
        {
            // Scale from the diagonal to the right
            for( Int iLoc=0; iLoc<mLoc; ++iLoc )
            {
                const Int i = A.GlobalRow(iLoc);
                if( i < iOff+diagLength )
                {
                    const Int k = i-iOff;
                    const Int j = k+jOff;
                    const Int jLeft = Max(j,0);
                    const Int jLeftLoc = A.LocalColOffset(jLeft);
                    const TDiag alpha = 
                        ( conjugate ? Conj(d_W_STAR.GetLocal(iLoc,0))
                                    : d_W_STAR.GetLocal(iLoc,0) );
                    blas::Scal
                    ( nLoc-jLeftLoc, alpha, &ABuf[iLoc+jLeftLoc*ldim], ldim );
                }
            }
        }    
    }
    else
    {
        DistMatrix<TDiag,Z,STAR> d_Z_STAR( A.Grid() );
        if( U == Z && V == STAR && d.ColAlign() == A.RowAlign() )
        {
            d_Z_STAR = LockedView( d );
        }
        else
        {
            d_Z_STAR.AlignWith( A );
            d_Z_STAR = d;
        }

        if( uplo == LOWER )
        {
            // Scale from the diagonal downwards
            for( Int jLoc=0; jLoc<nLoc; ++jLoc )
            {
                const Int j = A.GlobalCol(jLoc);
                if( j < jOff+diagLength )
                {
                    const Int k = j-jOff;
                    const Int i = k+iOff;
                    const Int iTop = Max(i,0);
                    const Int iTopLoc = A.LocalRowOffset(iTop);
                    const TDiag alpha = 
                        ( conjugate ? Conj(d_Z_STAR.GetLocal(jLoc,0))
                                    : d_Z_STAR.GetLocal(jLoc,0) );
                    blas::Scal
                    ( mLoc-iTopLoc, alpha, &ABuf[iTopLoc+jLoc*ldim], 1 );
                }
            }
        }
        else 
        {
            // Scale downward to the diagonal
            for( Int jLoc=0; jLoc<nLoc; ++jLoc )
            {
                const Int j = A.GlobalCol(jLoc);
                if( j >= jOff )
                {
                    const Int k = j-jOff;
                    const Int i = k+iOff;
                    const Int height = Min(i+1,m);
                    const Int localHeight = A.LocalRowOffset(height);
                    const TDiag alpha = 
                        ( conjugate ? Conj(d_Z_STAR.GetLocal(jLoc,0))
                                    : d_Z_STAR.GetLocal(jLoc,0) );
                    blas::Scal( localHeight, alpha, &ABuf[jLoc*ldim], 1 );
                }
            }
        }
    }
}

} // namespace elem

#endif // ifndef ELEM_DIAGONALSCALETRAPEZOID_HPP
