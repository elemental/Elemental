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

#include "elemental/blas-like/level1/Max.hpp"
#include "elemental/blas-like/level1/Scale.hpp"
#include "elemental/blas-like/level1/Swap.hpp"
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
    const Int minDim = Min(m,n);
    p.ResizeTo( minDim, 1 );
    q.ResizeTo( minDim, 1 );

    for( Int k=0; k<minDim; ++k )
    {
        auto ABR = ViewRange( A, k, k, m, n );

        // Find the index and value of the pivot candidate
        auto pivot = Max( ABR );
        const Int iPiv = pivot.indices[0] + k;
        const Int jPiv = pivot.indices[1] + k;
        p.Set( k, 0, iPiv+pivotOffset );
        q.Set( k, 0, jPiv+pivotOffset );

        // Swap the pivot row and current row
        if( iPiv != k )
        {
            auto aCurRow = ViewRange( A, k,      0, k+1,      n );
            auto aPivRow = ViewRange( A, iPiv, 0, iPiv+1, n );
            Swap( NORMAL, aCurRow, aPivRow );
        }

        // Swap the pivot column and current column
        if( jPiv != k )
        {
            auto aCurCol = ViewRange( A, 0, k,      m, k+1      );
            auto aPivCol = ViewRange( A, 0, jPiv, m, jPiv+1 );
            Swap( NORMAL, aCurCol, aPivCol );
        }

        // Now we can perform the update of the current panel
        const F alpha11 = A.Get(k,k);
        auto a21 = ViewRange( A, k+1, k,   m,   k+1 );
        auto a12 = ViewRange( A, k,   k+1, k+1, n   );
        auto A22 = ViewRange( A, k+1, k+1, m,   n   );
        if( alpha11 == F(0) )
            throw SingularMatrixException();
        const F alpha11Inv = F(1) / alpha11;
        Scale( alpha11Inv, a21 );
        Geru( F(-1), a21, a12, A22 );
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
#endif
    typedef BASE(F) Real;
    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);
    p.ResizeTo( minDim, 1 );
    q.ResizeTo( minDim, 1 );

    const Grid& g = A.Grid();

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

    for( Int k=0; k<minDim; ++k )
    {
        auto ABR = ViewRange( A, k, k, m, n );
        auto pivot = Max( ABR );
        const Int iPiv = pivot.indices[0] + k;
        const Int jPiv = pivot.indices[1] + k;
        p.Set( k, 0, iPiv+pivotOffset );
        q.Set( k, 0, jPiv+pivotOffset );

        // Perform the row pivot
        // TODO: Extract this into a routine
        if( iPiv != k )
        {
            const Int curOwnerRow = (colAlign+k   ) % colStride;
            const Int pivOwnerRow = (colAlign+iPiv) % colStride;
            if( pivOwnerRow == curOwnerRow )
            {
                if( g.Row() == curOwnerRow )
                {
                    const Int iLocCur = (k   -colShift) / colStride;
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
                    const Int iLoc = (k-colShift) / colStride;
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
        if( jPiv != k )
        {
            const Int curOwnerCol = (rowAlign+k   ) % rowStride;
            const Int pivOwnerCol = (rowAlign+jPiv) % rowStride;
            if( pivOwnerCol == curOwnerCol )
            {
                if( g.Col() == curOwnerCol )
                {
                    const Int jLocCur = (k   -rowShift) / rowStride;
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
                    const Int jLoc = (k-rowShift) / rowStride;
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
        const F alpha11 = A.Get(k,k);
        auto a21 = ViewRange( A, k+1, k,   m,   k+1 );
        auto a12 = ViewRange( A, k,   k+1, k+1, n   );
        auto A22 = ViewRange( A, k+1, k+1, m,   n   );
        if( alpha11 == F(0) )
            throw SingularMatrixException();
        const F alpha11Inv = F(1) / alpha11;
        Scale( alpha11Inv, a21 );
        Geru( F(-1), a21, a12, A22 );
    }
}

} // namespace lu
} // namespace elem

#endif // ifndef ELEM_LAPACK_LU_FULL_HPP
