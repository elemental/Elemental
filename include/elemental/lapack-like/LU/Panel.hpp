/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_LU_PANEL_HPP
#define ELEM_LAPACK_LU_PANEL_HPP

#include "elemental/blas-like/level1/Scale.hpp"
#include "elemental/blas-like/level2/Geru.hpp"

namespace elem {
namespace lu {

template<typename F>
inline void
Panel( Matrix<F>& A, Matrix<Int>& p, Int pivotOffset=0 )
{
#ifndef RELEASE
    CallStackEntry entry("lu::Panel");
    if( A.Width() != p.Height() || p.Width() != 1 )
        LogicError("p must be a vector that conforms with A");
#endif
    typedef BASE(F) Real;
    const Int n = A.Width();
    // Matrix views
    Matrix<F> 
        ATL, ATR,  A00, a01,     A02,  
        ABL, ABR,  a10, alpha11, a12,  
                   A20, a21,     A22;

    // Start the algorithm
    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    while( ATL.Width() < A.Width() )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ a01,     A02,
         /*************/ /**********************/
               /**/       a10, /**/ alpha11, a12,
          ABL, /**/ ABR,  A20, /**/ a21,     A22, 1 );

        //--------------------------------------------------------------------//
        // Find the index and value of the pivot candidate
        const Int offset = A00.Height();
        ValueInt<Real> pivot;
        pivot.value = FastAbs(alpha11.Get(0,0));
        pivot.index = offset;
        for( Int i=0; i<a21.Height(); ++i )
        {
            const Real value = FastAbs(a21.Get(i,0));
            if( value > pivot.value )
            {
                pivot.value = value;
                pivot.index = offset + i + 1;
            }
        }
        p.Set( offset, 0, pivot.index+pivotOffset );

        // Swap the pivot row and current row
        for( Int j=0; j<n; ++j )
        {
            const F temp = A.Get(offset,j);
            A.Set( offset,      j, A.Get(pivot.index,j) ); 
            A.Set( pivot.index, j, temp                 );
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
Panel
( DistMatrix<F,  STAR,STAR>& A, 
  DistMatrix<F,  MC,  STAR>& B, 
  DistMatrix<Int,STAR,STAR>& p, 
  Int pivotOffset=0 )
{
#ifndef RELEASE
    CallStackEntry entry("lu::Panel");
    if( A.Grid() != p.Grid() || p.Grid() != B.Grid() )
        LogicError("Matrices must be distributed over the same grid");
    if( A.Width() != B.Width() )
        LogicError("A and B must be the same width");
    if( A.Height() != p.Height() || p.Width() != 1 )
        LogicError("p must be a vector that conforms with A");
#endif
    typedef BASE(F) Real;
    const Grid& g = A.Grid();
    const Int r = g.Height();
    const Int colShift = B.ColShift();
    const Int colAlignment = B.ColAlignment();

    // Matrix views
    DistMatrix<F,STAR,STAR> 
        ATL(g), ATR(g),  A00(g), a01(g),     A02(g),  
        ABL(g), ABR(g),  a10(g), alpha11(g), a12(g),  
                         A20(g), a21(g),     A22(g);

    DistMatrix<F,MC,STAR>
        BL(g), BR(g),
        B0(g), b1(g), B2(g);

    // For packing rows of data for pivoting
    const Int n = A.Width();
    std::vector<F> pivotBuffer( n );

    // Start the algorithm
    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    PartitionRight( B, BL, BR, 0 );
    while( ATL.Height() < A.Height() )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ a01,     A02,
         /*************/ /**********************/
               /**/       a10, /**/ alpha11, a12,
          ABL, /**/ ABR,  A20, /**/ a21,     A22, 1 );

        RepartitionRight
        ( BL, /**/ BR,  
          B0, /**/ b1, B2, 1 );

        //--------------------------------------------------------------------//
        // Store the index/value of the local pivot candidate
        const Int offset = A00.Height();
        ValueInt<Real> localPivot;
        localPivot.value = FastAbs(alpha11.GetLocal(0,0));
        localPivot.index = offset;
        for( Int i=0; i<a21.Height(); ++i )
        {
            const Real value = FastAbs(a21.GetLocal(i,0));
            if( value > localPivot.value )
            {
                localPivot.value = value;
                localPivot.index = offset + i + 1;
            }
        }
        for( Int i=0; i<B.LocalHeight(); ++i )
        {
            const Real value = FastAbs(b1.GetLocal(i,0));
            if( value > localPivot.value )
            {
                localPivot.value = value;
                localPivot.index = A.Height() + colShift + i*r;
            }
        }

        // Compute and store the location of the new pivot
        const ValueInt<Real> pivot = 
            mpi::AllReduce( localPivot, mpi::MaxLocOp<Real>(), g.ColComm() );
        p.SetLocal(offset,0,pivot.index+pivotOffset);

        // Perform the pivot within this panel
        if( pivot.index < A.Height() )
        {
            // Pack pivot into temporary
            for( Int j=0; j<n; ++j )
                pivotBuffer[j] = A.GetLocal( pivot.index, j );
            // Replace pivot with current
            for( Int j=0; j<n; ++j )
                A.SetLocal( pivot.index, j, A.GetLocal(offset,j) );
        }
        else
        {
            // The owning row of the pivot row packs it into the row buffer
            // and then overwrites with the current row
            const Int relIndex = pivot.index - A.Height();
            const Int ownerRow = (colAlignment+relIndex) % r;
            if( g.Row() == ownerRow )
            {
                const int iLoc = (relIndex-colShift) / r;
                for( Int j=0; j<n; ++j )
                    pivotBuffer[j] = B.GetLocal( iLoc, j );
                for( Int j=0; j<n; ++j )
                    B.SetLocal( iLoc, j, A.GetLocal(offset,j) );
            }
            // The owning row broadcasts within process columns
            mpi::Broadcast( &pivotBuffer[0], n, ownerRow, g.ColComm() );
        }
        // Overwrite the current row with the pivot row
        for( Int j=0; j<n; ++j )
            A.SetLocal( offset, j, pivotBuffer[j] );

        // Now we can perform the update of the current panel
        const F alpha = alpha11.GetLocal(0,0);
        if( alpha == F(0) )
            throw SingularMatrixException();
        const F alpha11Inv = F(1) / alpha;
        Scale( alpha11Inv, a21 );
        Scale( alpha11Inv, b1  );
        Geru( F(-1), a21.Matrix(), a12.Matrix(), A22.Matrix() );
        Geru( F(-1), b1.Matrix(), a12.Matrix(), B2.Matrix() );
        //--------------------------------------------------------------------//

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, a01,     /**/ A02,
               /**/       a10, alpha11, /**/ a12,
         /*************/ /**********************/
          ABL, /**/ ABR,  A20, a21,     /**/ A22 );

        SlidePartitionRight
        ( BL,     /**/ BR,  
          B0, b1, /**/ B2 );
    }
}

} // namespace lu
} // namespace elem

#endif // ifndef ELEM_LAPACK_LU_PANEL_HPP
