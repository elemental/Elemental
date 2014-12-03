/*
   Copyright (c) 2009-2014, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, Stanford University, and the
   Georgia Insitute of Technology.
   All rights reserved.
 
   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_LDL_PROCESSFRONT_HPP
#define EL_LDL_PROCESSFRONT_HPP

namespace El {
namespace ldl {

template<typename F>
inline void ProcessFront( Matrix<F>& AL, Matrix<F>& ABR, bool conjugate )
{
    DEBUG_ONLY(
        CallStackEntry cse("ldl::ProcessFront");
        if( ABR.Height() != ABR.Width() )
            LogicError("ABR must be square");
        if( AL.Height() != AL.Width() + ABR.Width() )
            LogicError("AL and ABR don't have conformal dimensions");
    )
    const Int m = AL.Height();
    const Int n = AL.Width();
    const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );

    Matrix<F> d1;
    Matrix<F> S21;

    Matrix<F> S21T, S21B;
    Matrix<F> AL21T, AL21B;

    const Int bsize = Blocksize();
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);
        const Range<Int> ind1( k, k+nb ), 
                         ind2Vert( k+nb, m ), ind2Horz( k+nb, n );
        auto AL11 = AL( ind1,     ind1     );
        auto AL21 = AL( ind2Vert, ind1     );
        auto AL22 = AL( ind2Vert, ind2Horz );

        LDL( AL11, conjugate );
        AL11.GetDiagonal( d1 );

        Trsm( RIGHT, LOWER, orientation, UNIT, F(1), AL11, AL21 );

        S21 = AL21;
        DiagonalSolve( RIGHT, NORMAL, d1, AL21 );

        PartitionDown( S21, S21T, S21B, AL22.Width() );
        PartitionDown( AL21, AL21T, AL21B, AL22.Width() );
        Gemm( NORMAL, orientation, F(-1), S21, AL21T, F(1), AL22 );
        MakeTrapezoidal( LOWER, AL22 );
        Trrk( LOWER, NORMAL, orientation, F(-1), S21B, AL21B, F(1), ABR );
    }
}

template<typename F>
void ProcessFrontIntraPiv
( Matrix<F>& AL, Matrix<F>& subdiag, Matrix<Int>& piv, Matrix<F>& ABR,
  bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("ldl::ProcessFrontIntraPiv"))
    const Int n = AL.Width();
    const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );

    Matrix<F> ATL, ABL;
    PartitionDown( AL, ATL, ABL, n );

    LDL( ATL, subdiag, piv, conjugate );
    auto diag = ATL.GetDiagonal();

    PermuteCols( ABL, piv );
    Trsm( RIGHT, LOWER, orientation, UNIT, F(1), ATL, ABL );
    Matrix<F> SBL( ABL );

    QuasiDiagonalSolve( RIGHT, LOWER, diag, subdiag, ABL, conjugate );
    Trrk( LOWER, NORMAL, orientation, F(-1), SBL, ABL, F(1), ABR );
}

template<typename F> 
inline void ProcessFrontGeneral
( DistMatrix<F>& AL, DistMatrix<F>& ABR, bool conjugate=false )
{
    DEBUG_ONLY(
        CallStackEntry cse("ldl::ProcessFrontGeneral");
        if( ABR.Height() != ABR.Width() )
            LogicError("ABR must be square");
        if( AL.Height() != AL.Width()+ABR.Height() )
            LogicError("AL and ABR must have compatible dimensions");
        if( AL.Grid() != ABR.Grid() )
            LogicError("AL and ABR must use the same grid");
        if( ABR.ColAlign() !=
            (AL.ColAlign()+AL.Width()) % AL.Grid().Height() )
            LogicError("AL and ABR must have compatible col alignments");
        if( ABR.RowAlign() != 
            (AL.RowAlign()+AL.Width()) % AL.Grid().Width() )
            LogicError("AL and ABR must have compatible row alignments");
    )
    const Grid& g = AL.Grid();
    const Int m = AL.Height();
    const Int n = AL.Width();
    const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );

    DistMatrix<F,STAR,STAR> AL11_STAR_STAR(g);
    DistMatrix<F,STAR,STAR> d1_STAR_STAR(g);
    DistMatrix<F,VC,  STAR> AL21_VC_STAR(g);
    DistMatrix<F,VR,  STAR> AL21_VR_STAR(g);
    DistMatrix<F,STAR,MC  > S21Trans_STAR_MC(g);
    DistMatrix<F,STAR,MR  > AL21Trans_STAR_MR(g);

    DistMatrix<F,STAR,MC> leftL(g), leftR(g);
    DistMatrix<F,STAR,MR> rightL(g), rightR(g);
    DistMatrix<F> AL22T(g), AL22B(g);

    const Int bsize = Blocksize();
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);
        const Range<Int> ind1( k, k+nb ),
                         ind2Vert( k+nb, m ), ind2Horz( k+nb, n );
        auto AL11 = AL( ind1,     ind1     );
        auto AL21 = AL( ind2Vert, ind1     );
        auto AL22 = AL( ind2Vert, ind2Horz );

        AL11_STAR_STAR = AL11; 
        LocalLDL( AL11_STAR_STAR, conjugate );
        AL11_STAR_STAR.GetDiagonal( d1_STAR_STAR );
        AL11 = AL11_STAR_STAR;

        AL21_VC_STAR.AlignWith( AL22 );
        AL21_VC_STAR = AL21;
        LocalTrsm
        ( RIGHT, LOWER, orientation, UNIT, F(1), AL11_STAR_STAR, AL21_VC_STAR );

        S21Trans_STAR_MC.AlignWith( AL22 );
        AL21_VC_STAR.TransposePartialColAllGather( S21Trans_STAR_MC );
        DiagonalSolve( RIGHT, NORMAL, d1_STAR_STAR, AL21_VC_STAR );
        AL21_VR_STAR.AlignWith( AL22 );
        AL21_VR_STAR = AL21_VC_STAR;
        AL21Trans_STAR_MR.AlignWith( AL22 );
        AL21_VR_STAR.TransposePartialColAllGather
        ( AL21Trans_STAR_MR, conjugate );

        // Partition the update of the bottom-right corner into three pieces
        PartitionRight( S21Trans_STAR_MC, leftL, leftR, AL22.Width() );
        PartitionRight( AL21Trans_STAR_MR, rightL, rightR, AL22.Width() );
        PartitionDown( AL22, AL22T, AL22B, AL22.Width() );
        LocalTrrk
        ( LOWER, orientation, F(-1), leftL, rightL, F(1), AL22T );
        LocalGemm( orientation, NORMAL, F(-1), leftR, rightL, F(1), AL22B );
        LocalTrrk( LOWER, orientation, F(-1), leftR, rightR, F(1), ABR );

        DiagonalSolve( LEFT, NORMAL, d1_STAR_STAR, S21Trans_STAR_MC );
        AL21.TransposeRowFilterFrom( S21Trans_STAR_MC );
    }
}

template<typename F>
inline void ProcessFrontSquare
( DistMatrix<F>& AL, DistMatrix<F>& ABR, bool conjugate=false )
{
    DEBUG_ONLY(
        CallStackEntry cse("ldl::ProcessFrontSquare");
        if( ABR.Height() != ABR.Width() )
            LogicError("ABR must be square");
        if( AL.Height() != AL.Width()+ABR.Height() )
            LogicError("AL & ABR must have compatible dimensions");
        if( AL.Grid() != ABR.Grid() )
            LogicError("AL & ABR must use the same grid");
        if( ABR.ColAlign() !=
            (AL.ColAlign()+AL.Width()) % AL.Grid().Height() )
            LogicError("AL & ABR must have compatible col alignments");
        if( ABR.RowAlign() != 
            (AL.RowAlign()+AL.Width()) % AL.Grid().Width() )
            LogicError("AL & ABR must have compatible row alignments");
    )
    const Grid& g = AL.Grid();
    const Int m = AL.Height();
    const Int n = AL.Width();
    const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );
    DEBUG_ONLY(
        if( g.Height() != g.Width() )
            LogicError("This routine assumes a square process grid");
    )

    // Find the process holding our transposed data
    const int r = g.Height();
    int transposeRank;
    {
        const int colAlign = AL.ColAlign();
        const int rowAlign = AL.RowAlign();
        const int colShift = AL.ColShift();
        const int rowShift = AL.RowShift();

        const int transposeRow = (colAlign+rowShift) % r;
        const int transposeCol = (rowAlign+colShift) % r;
        transposeRank = transposeRow + r*transposeCol;
    }
    const bool onDiagonal = ( transposeRank == g.VCRank() );

    DistMatrix<F,STAR,STAR> AL11_STAR_STAR(g);
    DistMatrix<F,STAR,STAR> d1_STAR_STAR(g);
    DistMatrix<F,VC,  STAR> AL21_VC_STAR(g);
    DistMatrix<F,STAR,MC  > S21Trans_STAR_MC(g);
    DistMatrix<F,STAR,MR  > AL21Trans_STAR_MR(g);

    DistMatrix<F,STAR,MC> leftL(g), leftR(g);
    DistMatrix<F,STAR,MR> rightL(g), rightR(g);
    DistMatrix<F> AL22T(g), AL22B(g);

    const Int bsize = Blocksize();
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);
        const Range<Int> ind1( k, k+nb ),
                         ind2Vert( k+nb, m ), ind2Horz( k+nb, n );
        auto AL11 = AL( ind1,     ind1     );
        auto AL21 = AL( ind2Vert, ind1     );
        auto AL22 = AL( ind2Vert, ind2Horz );

        AL11_STAR_STAR = AL11; 
        LocalLDL( AL11_STAR_STAR, conjugate );
        AL11_STAR_STAR.GetDiagonal( d1_STAR_STAR );
        AL11 = AL11_STAR_STAR;

        AL21_VC_STAR.AlignWith( AL22 );
        AL21_VC_STAR = AL21;
        LocalTrsm
        ( RIGHT, LOWER, orientation, UNIT, F(1), AL11_STAR_STAR, AL21_VC_STAR );

        S21Trans_STAR_MC.AlignWith( AL22 );
        AL21_VC_STAR.TransposePartialColAllGather( S21Trans_STAR_MC );
        // SendRecv to form AL21^T[* ,MR] from S21^T[* ,MC], then conjugate
        // if necessary.
        AL21Trans_STAR_MR.AlignWith( AL22 );
        AL21Trans_STAR_MR.Resize( AL21.Width(), AL21.Height() );
        {
            if( onDiagonal )
            {
                const int size = AL21.LocalHeight()*AL11.Width();    
                MemCopy
                ( AL21Trans_STAR_MR.Buffer(), 
                  S21Trans_STAR_MC.Buffer(), size );
            }
            else
            {
                const int sendSize = AL21.LocalHeight()*AL11.Width();
                const int recvSize = 
                    (AL22.LocalWidth()+ABR.LocalWidth())*AL11.Height();
                // We know that the ldim is the height since we have manually
                // created both temporary matrices.
                mpi::SendRecv
                ( S21Trans_STAR_MC.Buffer(), sendSize, transposeRank,
                  AL21Trans_STAR_MR.Buffer(), 
                  recvSize, transposeRank, g.VCComm() );
            }
            DiagonalSolve( LEFT, NORMAL, d1_STAR_STAR, AL21Trans_STAR_MR );
            if( conjugate )
                Conjugate( AL21Trans_STAR_MR );
        }

        // Partition the update of the bottom-right corner into three pieces
        PartitionRight( S21Trans_STAR_MC, leftL, leftR, AL22.Width() );
        PartitionRight( AL21Trans_STAR_MR, rightL, rightR, AL22.Width() );
        PartitionDown( AL22, AL22T, AL22B, AL22.Width() );
        LocalTrrk( LOWER, orientation, F(-1), leftL, rightL, F(1), AL22T );
        LocalGemm( orientation, NORMAL, F(-1), leftR, rightL, F(1), AL22B );
        LocalTrrk( LOWER, orientation, F(-1), leftR, rightR, F(1), ABR );

        DiagonalSolve( LEFT, NORMAL, d1_STAR_STAR, S21Trans_STAR_MC );
        AL21.TransposeRowFilterFrom( S21Trans_STAR_MC );
    }
}

template<typename F> 
inline void ProcessFront
( DistMatrix<F>& AL, DistMatrix<F>& ABR, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("ldl::ProcessFront"))
    const Grid& grid = AL.Grid();
    if( grid.Height() == grid.Width() )
        ProcessFrontSquare( AL, ABR, conjugate );
    else
        ProcessFrontGeneral( AL, ABR, conjugate );
}

template<typename F>
void ProcessFrontIntraPiv
( DistMatrix<F>& AL, DistMatrix<F,MD,STAR>& subdiag, 
  DistMatrix<Int,VC,STAR>& p, DistMatrix<F>& ABR, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("ldl::ProcessFrontIntraPiv"))
    const Grid& g = AL.Grid();
    const Int n = AL.Width();
    const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );
    
    DistMatrix<F> ATL(g), ABL(g);
    PartitionDown( AL, ATL, ABL, n );

    LDL( ATL, subdiag, p, conjugate );
    auto diag = ATL.GetDiagonal();

    PermuteCols( ABL, p );
    Trsm( RIGHT, LOWER, orientation, UNIT, F(1), ATL, ABL );
    DistMatrix<F,MC,STAR> SBL_MC_STAR(g);
    SBL_MC_STAR.AlignWith( ABR );
    SBL_MC_STAR = ABL;

    QuasiDiagonalSolve( RIGHT, LOWER, diag, subdiag, ABL, conjugate );
    DistMatrix<F,VR,STAR> ABL_VR_STAR(g);
    DistMatrix<F,STAR,MR> ABLTrans_STAR_MR(g);
    ABL_VR_STAR.AlignWith( ABR );
    ABLTrans_STAR_MR.AlignWith( ABR );
    ABL_VR_STAR = ABL;
    ABL_VR_STAR.TransposePartialColAllGather( ABLTrans_STAR_MR, conjugate );
    LocalTrrk( LOWER, F(-1), SBL_MC_STAR, ABLTrans_STAR_MR, F(1), ABR );
}

} // namespace ldl
} // namespace El

#endif // ifndef EL_LDL_PROCESSFRONT_HPP
