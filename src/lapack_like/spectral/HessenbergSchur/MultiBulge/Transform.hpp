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

// Apply (with replacement) V' from the left
// -----------------------------------------

template<typename Field>
void TransformRows
( const Matrix<Field>& V,
        Matrix<Field>& A )
{
    EL_DEBUG_CSE
    // TODO(poulson): Consider forming chunk-by-chunk to save memory
    Matrix<Field> ACopy( A );
    Gemm( ADJOINT, NORMAL, Field(1), V, ACopy, A );
}

template<typename Field>
void TransformRows
( const Matrix<Field>& V,
        DistMatrix<Field,MC,MR,BLOCK>& A )
{
    EL_DEBUG_CSE
    const Int height = A.Height();
    const Grid& grid = A.Grid();

    const Int blockHeight = A.BlockHeight();
    const Int firstBlockHeight = blockHeight - A.ColCut();
    if( height <= firstBlockHeight || grid.Height() == 1 )
    {
        if( grid.Row() == A.RowOwner(0) )
        {
            // This process row can locally update its portion of A
            TransformRows( V, A.Matrix() );
        }
    }
    else if( height <= firstBlockHeight + blockHeight )
    {
        const int firstRow = A.RowOwner( 0 );
        const int secondRow = A.RowOwner( firstBlockHeight );
        if( grid.Row() == firstRow )
        {
            //
            // Replace A with
            //
            //   | VLeft, VRight |' | ATop    |,
            //                      | ABottom |
            //
            // where ATop is owned by this process row and ABottom by the next.
            //
            auto VLeft = V( ALL, IR(0,firstBlockHeight) );

            // Partition space for the combined matrix
            Matrix<Field> ACombine( height, A.LocalWidth() );
            auto ATop = ACombine( IR(0,firstBlockHeight), ALL );
            auto ABottom = ACombine( IR(firstBlockHeight,END), ALL );

            // Copy our portion into the combined matrix
            ATop = A.LockedMatrix();

            // Exchange the data
            El::SendRecv( ATop, ABottom, A.ColComm(), secondRow, secondRow );

            // Form our portion of the result
            Gemm( ADJOINT, NORMAL, Field(1), VLeft, ACombine, A.Matrix() );
        }
        else if( grid.Row() == secondRow )
        {
            //
            // Replace A with
            //
            //   | VLeft, VRight |' | ATop    |,
            //                      | ABottom |
            //
            // where ATop is owned by the previous process row and ABottom by
            // this one.
            //
            auto VRight = V( ALL, IR(firstBlockHeight,END) );

            // Partition space for the combined matrix
            Matrix<Field> ACombine( height, A.LocalWidth() );
            auto ATop = ACombine( IR(0,firstBlockHeight), ALL );
            auto ABottom = ACombine( IR(firstBlockHeight,END), ALL );

            // Copy our portion into the combined matrix
            ABottom = A.LockedMatrix();

            // Exchange the data
            El::SendRecv( ABottom, ATop, A.ColComm(), firstRow, firstRow );

            // Form our portion of the result
            Gemm( ADJOINT, NORMAL, Field(1), VRight, ACombine, A.Matrix() );
        }
    }
    else
    {
        // Fall back to the entire process column interacting.
        // TODO(poulson): Only form the subset of the result that we need.
        DistMatrix<Field,STAR,MR,BLOCK> A_STAR_MR( A );
        Matrix<Field> ALocCopy( A_STAR_MR.Matrix() );
        Gemm( ADJOINT, NORMAL, Field(1), V, ALocCopy, A_STAR_MR.Matrix() );
        A = A_STAR_MR;
    }
}

template<typename Field>
void TransformColumns
( const Matrix<Field>& V,
        Matrix<Field>& A )
{
    EL_DEBUG_CSE
    // TODO(poulson): Consider forming chunk-by-chunk to save memory
    Matrix<Field> ACopy( A );
    Gemm( NORMAL, NORMAL, Field(1), ACopy, V, A );
}

// Apply (with replacement) V from the right
template<typename Field>
void TransformColumns
( const Matrix<Field>& V,
        DistMatrix<Field,MC,MR,BLOCK>& A )
{
    EL_DEBUG_CSE
    const Int width = A.Width();
    const Grid& grid = A.Grid();

    const Int blockWidth = A.BlockWidth();
    const Int firstBlockWidth = blockWidth - A.RowCut();
    if( width <= firstBlockWidth || grid.Width() == 1 )
    {
        if( grid.Col() == A.ColOwner(0) )
        {
            // This process row can locally update its portion of A
            TransformColumns( V, A.Matrix() );
        }
    }
    else if( width <= firstBlockWidth + blockWidth )
    {
        const int firstCol = A.ColOwner( 0 );
        const int secondCol = A.ColOwner( firstBlockWidth );
        if( grid.Col() == firstCol )
        {
            //
            // Replace A with
            //
            //   | ALeft, ARight | | VLeft, VRight |,
            //
            // where ALeft is owned by this process column and ARight by the
            // next.
            //

            // Partition space for the combined matrix
            Matrix<Field> ACombine( A.LocalHeight(), width );
            auto ALeft = ACombine( ALL, IR(0,firstBlockWidth) );
            auto ARight = ACombine( ALL, IR(firstBlockWidth,END) );

            // Copy our portion into the combined matrix
            ALeft = A.LockedMatrix();

            // Exchange the data
            El::SendRecv( ALeft, ARight, A.RowComm(), secondCol, secondCol );

            // Form our portion of the result
            auto VLeft = V( ALL, IR(0,firstBlockWidth) );
            Gemm( NORMAL, NORMAL, Field(1), ACombine, VLeft, A.Matrix() );
        }
        else if( grid.Col() == secondCol )
        {
            //
            // Replace A with
            //
            //   | ALeft, ARight | | VLeft, VRight |,
            //
            // where ALeft is owned by the previous process column and ARight
            // by this one.
            //

            // Partition space for the combined matrix
            Matrix<Field> ACombine( A.LocalHeight(), width );
            auto ALeft = ACombine( ALL, IR(0,firstBlockWidth) );
            auto ARight = ACombine( ALL, IR(firstBlockWidth,END) );

            // Copy our portion into the combined matrix
            ARight = A.LockedMatrix();

            // Exchange the data
            El::SendRecv( ARight, ALeft, A.RowComm(), firstCol, firstCol );

            // Form our portion of the result
            auto VRight = V( ALL, IR(firstBlockWidth,END) );
            Gemm( NORMAL, NORMAL, Field(1), ACombine, VRight, A.Matrix() );
        }
    }
    else
    {
        // Fall back to the entire process column interacting.
        // TODO(poulson): Only form the subset of the result that we need.
        DistMatrix<Field,MC,STAR,BLOCK> A_MC_STAR( A );
        Matrix<Field> ALocCopy( A_MC_STAR.Matrix() );
        Gemm( NORMAL, NORMAL, Field(1), ALocCopy, V, A_MC_STAR.Matrix() );
        A = A_MC_STAR;
    }
}

} // namespace multibulge
} // namespace hess_schur
} // namespace El

#endif // ifndef EL_HESS_SCHUR_MULTIBULGE_TRANSFORM_HPP
