/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_LU_FULL_HPP
#define ELEM_LAPACK_LU_FULL_HPP

#include "elemental/blas-like/level1/Scale.hpp"
#include "elemental/blas-like/level2/Geru.hpp"

namespace elem {
namespace lu {

template<typename F>
inline void
Full( Matrix<F>& A, Matrix<Int>& p, Matrix<Int>& q, Int pivotOffset=0 )
{
#ifndef RELEASE
    CallStackEntry entry("lu::Panel");
#endif
    const Int m = A.Height();
    const Int n = A.Width();
#ifndef RELEASE
    if( p.Height() != Min(m,n) || p.Width() != 1 )
        LogicError("p must be a vector that conforms with A");
    if( q.Height() != Min(m,n) || q.Width() != 1 )
        LogicError("q must be a vector that conforms with A");
#endif
    typedef BASE(F) Real;

    // Matrix views
    Matrix<F> 
        ATL, ATR,  A00, a01,     A02,  
        ABL, ABR,  a10, alpha11, a12,  
                   A20, a21,     A22;

    // Start the algorithm
    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    while( ATL.Height() < A.Height() && ATL.Width() < A.Width() )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ a01,     A02,
         /*************/ /**********************/
               /**/       a10, /**/ alpha11, a12,
          ABL, /**/ ABR,  A20, /**/ a21,     A22, 1 );

        //--------------------------------------------------------------------//
        // Find the index and value of the pivot candidate
        const Int offset = A00.Height();
        ValueIntPair<Real> pivot;
        pivot.value = -1;
        pivot.indices[0] = -1;
        pivot.indices[1] = -1;
        for( Int j=0; j<ABR.Width(); ++j )
        {
            for( Int i=0; i<ABR.Height(); ++i )
            {
                const Real value = FastAbs(ABR.Get(i,j));
                if( value > pivot.value )
                {
                    pivot.value = value;
                    pivot.indices[0] = offset + i + 1;
                    pivot.indices[1] = offset + j + 1;
                }
            }
        }
        p.Set( offset, 0, pivot.indices[0]+pivotOffset );
        q.Set( offset, 0, pivot.indices[1]+pivotOffset );

        // Swap the pivot row and current row
        for( Int j=0; j<n; ++j )
        {
            const F temp = A.Get( offset, j ); 
            A.Set( offset,           j, A.Get(pivot.indices[0],j) ); 
            A.Set( pivot.indices[0], j, temp                      );
        }

        // Swap the pivot column and current column
        for( Int i=0; i<m; ++i )
        {
            const F temp = A.Get( i, offset );
            A.Set( i, offset,           A.Get(i,pivot.indices[1]) );
            A.Set( i, pivot.indices[1], temp                      );
        }

        // Now we can perform the update of the current panel
        const F alpha = alpha11.Get(0,0);
        if( alpha == F(0) )
            throw SingularMatrixException();
        const F alpha11Inv = F(1) / alpha;
        Scale( alpha11Inv, a21 );
        Geru( F(-1), a21, a12, A22 );
        //--------------------------------------------------------------------//

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, a01,     /**/ A02,
               /**/       a10, alpha11, /**/ a12,
         /*************/ /**********************/
          ABL, /**/ ABR,  A20, a21,     /**/ A22 );
    }
}

template<typename F>
inline void
Full
( DistMatrix<F>& A, 
  DistMatrix<Int,VC,STAR>& p, 
  DistMatrix<Int,VC,STAR>& q,
  Int pivotOffset=0 )
{
#ifndef RELEASE
    CallStackEntry entry("lu::Panel");
    if( A.Grid() != p.Grid() || p.Grid() != q.Grid() )
        LogicError("Matrices must be distributed over the same grid");
    const Int m = A.Height();
    const Int n = A.Width();
    if( p.Height() != Min(m,n) || p.Width() != 1 )
        LogicError("p must be a vector that conforms with A");
    if( q.Height() != Min(m,n) || q.Width() != 1 )
        LogicError("q must be a vector that conforms with A");
#endif
    typedef BASE(F) Real;
    const Grid& g = A.Grid();

    // Matrix views
    DistMatrix<F> 
        ATL(g), ATR(g),  A00(g), a01(g),     A02(g),  
        ABL(g), ABR(g),  a10(g), alpha11(g), a12(g),  
                         A20(g), a21(g),     A22(g);

    // For packing rows/columns of data for pivoting
    const Int mLocal = A.LocalHeight();
    const Int nLocal = A.LocalWidth();
    std::vector<F> pivotBuffer( Max(mLocal,nLocal) );

    const Int colAlign = A.ColAlignment();
    const Int rowAlign = A.RowAlignment();
    const Int colShift = A.ColShift();
    const Int rowShift = A.RowShift();
    const Int colStride = A.ColStride();
    const Int rowStride = A.RowStride();

    // Start the algorithm
    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    while( ATL.Height() < A.Height() && ATL.Width() < A.Width() )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ a01,     A02,
         /*************/ /**********************/
               /**/       a10, /**/ alpha11, a12,
          ABL, /**/ ABR,  A20, /**/ a21,     A22, 1 );

        //--------------------------------------------------------------------//
        // Store the index/value of the local pivot candidate
        const Int offset = A00.Height();
        ValueIntPair<Real> localPivot;
        localPivot.value = -1;
        localPivot.indices[0] = -1;
        localPivot.indices[1] = -1;
        const Int mLocalBR = ABR.LocalHeight();
        const Int nLocalBR = ABR.LocalWidth();
        const Int colShiftBR = ABR.ColShift();
        const Int rowShiftBR = ABR.RowShift();
        for( Int jLoc=0; jLoc<nLocalBR; ++jLoc )
        {
            const Int j = rowShiftBR + jLoc*rowStride;
            for( Int iLoc=0; iLoc<mLocalBR; ++iLoc )
            {
                const Int i = colShiftBR + iLoc*colStride;
                const Real value = FastAbs(ABR.GetLocal(iLoc,jLoc));
                if( value > localPivot.value )
                {
                    localPivot.value = value;
                    localPivot.indices[0] = offset + i;
                    localPivot.indices[1] = offset + j;
                }
            }
        }

        // Compute and store the location of the new pivot
        const ValueIntPair<Real> pivot = 
            mpi::AllReduce( localPivot, mpi::MaxLocPairOp<Real>(), g.VCComm() );
        p.Set(offset,0,pivot.indices[0]+pivotOffset);
        q.Set(offset,0,pivot.indices[1]+pivotOffset);

        // Perform the row pivot
        // TODO: Extract this into a routine
        const Int iPiv = pivot.indices[0];
        if( iPiv != offset )
        {
            const Int curOwnerRow = (colAlign+offset) % colStride;
            const Int pivOwnerRow = (colAlign+iPiv  ) % colStride;
            if( pivOwnerRow == curOwnerRow )
            {
                if( g.Row() == curOwnerRow )
                {
                    const Int iLocCur = (offset-colShift) / colStride;
                    const Int iLocPiv = (iPiv-colShift) / colStride;
                    for( Int jLoc=0; jLoc<nLocal; ++jLoc )
                    {
                        const F temp = A.GetLocal(iLocPiv,jLoc);
                        A.SetLocal( iLocPiv, jLoc, A.GetLocal(iLocCur,jLoc) );
                        A.SetLocal( iLocCur,       jLoc, temp               );
                    }
                }
            }
            else
            {
                if( g.Row() == curOwnerRow )
                {
                    const Int iLoc = (offset-colShift) / colStride;
                    for( Int jLoc=0; jLoc<nLocal; ++jLoc )    
                        pivotBuffer[jLoc] = A.GetLocal(iLoc,jLoc);
                    mpi::SendRecv
                    ( pivotBuffer.data(), nLocal, 
                      pivOwnerRow, pivOwnerRow, g.ColComm() );
                    for( Int jLoc=0; jLoc<nLocal; ++jLoc )
                        A.SetLocal( iLoc, jLoc, pivotBuffer[jLoc] );
                } 
                else if( g.Row() == pivOwnerRow )
                {
                    const Int iLoc = (iPiv-colShift) / colStride;
                    for( Int jLoc=0; jLoc<nLocal; ++jLoc )    
                        pivotBuffer[jLoc] = A.GetLocal(iLoc,jLoc);
                    mpi::SendRecv
                    ( pivotBuffer.data(), nLocal, 
                      curOwnerRow, curOwnerRow, g.ColComm() );
                    for( Int jLoc=0; jLoc<nLocal; ++jLoc )
                        A.SetLocal( iLoc, jLoc, pivotBuffer[jLoc] );
                }
            }
        }

        // Perform the column pivot
        // TODO: Extract this into a routine
        const Int jPiv = pivot.indices[1];
        if( jPiv != offset )
        {
            const Int curOwnerCol = (rowAlign+offset) % rowStride;
            const Int pivOwnerCol = (rowAlign+jPiv  ) % rowStride;
            if( pivOwnerCol == curOwnerCol )
            {
                if( g.Col() == curOwnerCol )
                {
                    const Int jLocCur = (offset-rowShift) / rowStride;
                    const Int jLocPiv = (jPiv-rowShift) / rowStride;
                    for( Int iLoc=0; iLoc<mLocal; ++iLoc )
                    {
                        const F temp = A.GetLocal(iLoc,jLocPiv);
                        A.SetLocal( iLoc, jLocPiv, A.GetLocal(iLoc,jLocCur) );
                        A.SetLocal( iLoc, jLocCur, temp                     );
                    }
                }
            }
            else
            {
                if( g.Col() == curOwnerCol )
                {
                    const Int jLoc = (offset-rowShift) / rowStride;
                    for( Int iLoc=0; iLoc<mLocal; ++iLoc )    
                        pivotBuffer[iLoc] = A.GetLocal(iLoc,jLoc);
                    mpi::SendRecv
                    ( pivotBuffer.data(), mLocal, 
                      pivOwnerCol, pivOwnerCol, g.RowComm() );
                    for( Int iLoc=0; iLoc<mLocal; ++iLoc )
                        A.SetLocal( iLoc, jLoc, pivotBuffer[iLoc] );
                } 
                else if( g.Col() == pivOwnerCol )
                {
                    const Int jLoc = (jPiv-rowShift) / rowStride;
                    for( Int iLoc=0; iLoc<mLocal; ++iLoc )    
                        pivotBuffer[iLoc] = A.GetLocal(iLoc,jLoc);
                    mpi::SendRecv
                    ( pivotBuffer.data(), mLocal, 
                      curOwnerCol, curOwnerCol, g.RowComm() );
                    for( Int iLoc=0; iLoc<mLocal; ++iLoc )
                        A.SetLocal( iLoc, jLoc, pivotBuffer[iLoc] );
                }
            }
        }

        // Now we can perform the update of the current panel
        const F alpha = alpha11.Get(0,0);
        if( alpha == F(0) )
            throw SingularMatrixException();
        const F alpha11Inv = F(1) / alpha;
        Scale( alpha11Inv, a21 );
        Geru( F(-1), a21, a12, A22 );
        //--------------------------------------------------------------------//

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, a01,     /**/ A02,
               /**/       a10, alpha11, /**/ a12,
         /*************/ /**********************/
          ABL, /**/ ABR,  A20, a21,     /**/ A22 );
    }
}

} // namespace lu
} // namespace elem

#endif // ifndef ELEM_LAPACK_LU_FULL_HPP
