/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_HESS_SCHUR_MULTIBULGE_TRANSFORM_HPP
#define EL_HESS_SCHUR_MULTIBULGE_TRANSFORM_HPP

namespace El {
namespace hess_schur {
namespace multibulge {

// Apply (with replacement) Z' from the left
template<typename F>
void TransformRows
( const Matrix<F>& Z,
        DistMatrix<F,MC,MR,BLOCK>& H )
{
    DEBUG_CSE
    const Int height = H.Height();
    const Grid& grid = H.Grid();

    const Int blockHeight = H.BlockHeight();
    const Int firstBlockHeight = blockHeight - H.ColCut();
    if( height <= firstBlockHeight || grid.Height() == 1 )
    {
        if( grid.Row() == H.RowOwner(0) )
        {
            // This process row can locally update its portion of H
            Matrix<F> HLocCopy( H.Matrix() );
            Gemm( ADJOINT, NORMAL, F(1), Z, HLocCopy, H.Matrix() );
        }
    }
    else if( height <= firstBlockHeight + blockHeight )
    {
        const bool firstRow = H.RowOwner( 0 );
        const bool secondRow = H.RowOwner( firstBlockHeight );
        if( grid.Row() == firstRow )
        {
            // 
            // Replace H with 
            //
            //   | ZLeft, ZRight |' | HTop    |,
            //                      | HBottom |
            //
            // where HTop is owned by this process row and HBottom by the next.
            //
            auto ZLeft = Z( ALL, IR(0,firstBlockHeight) );

            // Partition space for the combined matrix
            Matrix<F> HCombine( height, H.LocalWidth() );
            auto HTop = HCombine( IR(0,firstBlockHeight), ALL );
            auto HBottom = HCombine( IR(firstBlockHeight,END), ALL );

            // Copy our portion into the combined matrix
            HTop = H.LockedMatrix();

            // Exchange the data
            El::SendRecv( HTop, HBottom, H.ColComm(), secondRow, secondRow );
            
            // Form our portion of the result
            Gemm( ADJOINT, NORMAL, F(1), ZLeft, HCombine, H.Matrix() );
        }
        else if( grid.Row() == secondRow )
        {
            // 
            // Replace H with 
            //
            //   | ZLeft, ZRight |' | HTop    |,
            //                      | HBottom |
            //
            // where HTop is owned by the previous process row and HBottom by
            // this one.
            //
            auto ZRight = Z( ALL, IR(firstBlockHeight,END) );

            // Partition space for the combined matrix
            Matrix<F> HCombine( height, H.LocalWidth() );
            auto HTop = HCombine( IR(0,firstBlockHeight), ALL );
            auto HBottom = HCombine( IR(firstBlockHeight,END), ALL );

            // Copy our portion into the combined matrix
            HBottom = H.LockedMatrix();

            // Exchange the data
            El::SendRecv( HBottom, HTop, H.ColComm(), firstRow, firstRow );
            
            // Form our portion of the result
            Gemm( ADJOINT, NORMAL, F(1), ZRight, HCombine, H.Matrix() );
        }
    }
    else
    {
        // Fall back to the entire process column interacting.
        // TODO(poulson): Only form the subset of the result that we need.
        DistMatrix<F,STAR,MR,BLOCK> H_STAR_MR( H );
        Matrix<F> HLocCopy( H_STAR_MR.Matrix() );
        Gemm( ADJOINT, NORMAL, F(1), Z, HLocCopy, H_STAR_MR.Matrix() );
        H = H_STAR_MR;
    }
}

// Apply (with replacement) Z from the right
template<typename F>
void TransformColumns
( const Matrix<F>& Z,
        DistMatrix<F,MC,MR,BLOCK>& H )
{
    DEBUG_CSE
    const Int width = H.Width();
    const Grid& grid = H.Grid();

    const Int blockWidth = H.BlockWidth();
    const Int firstBlockWidth = blockWidth - H.RowCut();
    if( width <= firstBlockWidth || grid.Width() == 1 )
    {
        if( grid.Col() == H.ColOwner(0) )
        {
            // This process row can locally update its portion of H
            Matrix<F> HLocCopy( H.Matrix() );
            Gemm( NORMAL, NORMAL, F(1), HLocCopy, Z, H.Matrix() );
        }
    }
    else if( width <= firstBlockWidth + blockWidth )
    {
        const bool firstCol = H.ColOwner( 0 );
        const bool secondCol = H.ColOwner( firstBlockWidth );
        if( grid.Col() == firstCol )
        {
            // 
            // Replace H with 
            //
            //   | HLeft, HRight | | ZLeft, ZRight |,
            //
            // where HLeft is owned by this process column and HRight by the
            // next.
            //
            auto ZLeft = Z( ALL, IR(0,firstBlockWidth) );

            // Partition space for the combined matrix
            Matrix<F> HCombine( H.LocalHeight(), width );
            auto HLeft = HCombine( ALL, IR(0,firstBlockWidth) );
            auto HRight = HCombine( ALL, IR(firstBlockWidth,END) );

            // Copy our portion into the combined matrix
            HLeft = H.LockedMatrix();

            // Exchange the data
            El::SendRecv( HLeft, HRight, H.RowComm(), secondCol, secondCol );
            
            // Form our portion of the result
            Gemm( NORMAL, NORMAL, F(1), HCombine, ZLeft, H.Matrix() );
        }
        else if( grid.Col() == secondCol )
        {
            // 
            // Replace H with 
            //
            //   | HLeft, HRight | | ZLeft, ZRight |,
            //
            // where HLeft is owned by the previous process column and HRight
            // by this one.
            //
            auto ZRight = Z( ALL, IR(firstBlockWidth,END) );

            // Partition space for the combined matrix
            Matrix<F> HCombine( H.LocalHeight(), width );
            auto HLeft = HCombine( ALL, IR(0,firstBlockWidth) );
            auto HRight = HCombine( ALL, IR(firstBlockWidth,END) );

            // Copy our portion into the combined matrix
            HRight = H.LockedMatrix();

            // Exchange the data
            El::SendRecv( HRight, HLeft, H.RowComm(), firstCol, firstCol );
            
            // Form our portion of the result
            Gemm( NORMAL, NORMAL, F(1), HCombine, ZRight, H.Matrix() );
        }
    }
    else
    {
        // Fall back to the entire process column interacting.
        // TODO(poulson): Only form the subset of the result that we need.
        DistMatrix<F,MC,STAR,BLOCK> H_MC_STAR( H );
        Matrix<F> HLocCopy( H_MC_STAR.Matrix() );
        Gemm( NORMAL, NORMAL, F(1), HLocCopy, Z, H_MC_STAR.Matrix() );
        H = H_MC_STAR;
    }
}

} // namespace multibulge
} // namespace hess_schur
} // namespace El

#endif // ifndef EL_HESS_SCHUR_MULTIBULGE_TRANSFORM_HPP
