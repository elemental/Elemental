/*
   Copyright (c) 2009-2015, Jack Poulson.
   All rights reserved.
 
   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_OPTIMIZATION_REGQSDLDL_PROCESSFRONT_HPP
#define EL_OPTIMIZATION_REGQSDLDL_PROCESSFRONT_HPP

namespace El {
namespace reg_qsd_ldl {

template<typename F>
inline void ProcessFront
( Matrix<F>& AL, Matrix<F>& ABR, Base<F> pivTol, 
  const Matrix<Base<F>>& regCand, Matrix<Base<F>>& reg,
  bool aPriori )
{
    DEBUG_ONLY(
        CallStackEntry cse("reg_qsd_ldl::ProcessFront");
        if( ABR.Height() != ABR.Width() )
            LogicError("ABR must be square");
        if( AL.Height() != AL.Width() + ABR.Width() )
            LogicError("AL and ABR don't have conformal dimensions");
    )
    const Int m = AL.Height();
    const Int n = AL.Width();

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

        auto regCand1 = regCand( ind1, IR(0,1) );
        auto reg1 = reg( ind1, IR(0,1) );

        RegularizedQSDLDL( AL11, pivTol, regCand1, reg1, aPriori );
        GetDiagonal( AL11, d1 );

        Trsm( RIGHT, LOWER, ADJOINT, UNIT, F(1), AL11, AL21 );

        S21 = AL21;
        DiagonalSolve( RIGHT, NORMAL, d1, AL21 );

        PartitionDown( S21, S21T, S21B, AL22.Width() );
        PartitionDown( AL21, AL21T, AL21B, AL22.Width() );
        Gemm( NORMAL, ADJOINT, F(-1), S21, AL21T, F(1), AL22 );
        MakeTrapezoidal( LOWER, AL22 );
        Trrk( LOWER, NORMAL, ADJOINT, F(-1), S21B, AL21B, F(1), ABR );
    }
}

template<typename F> 
inline void ProcessFrontGeneral
( DistMatrix<F>& AL, DistMatrix<F>& ABR, Base<F> pivTol,
  const DistMatrix<Base<F>,VC,STAR>& regCand, 
        DistMatrix<Base<F>,VC,STAR>& reg,
  bool aPriori )
{
    DEBUG_ONLY(
        CallStackEntry cse("reg_qsd_ldl::ProcessFrontGeneral");
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

    DistMatrix<F,STAR,STAR> AL11_STAR_STAR(g);
    DistMatrix<F,STAR,STAR> d1_STAR_STAR(g);
    DistMatrix<F,VC,  STAR> AL21_VC_STAR(g);
    DistMatrix<F,VR,  STAR> AL21_VR_STAR(g);
    DistMatrix<F,STAR,MC  > S21Trans_STAR_MC(g);
    DistMatrix<F,STAR,MR  > AL21Trans_STAR_MR(g);

    DistMatrix<F,STAR,MC> leftL(g), leftR(g);
    DistMatrix<F,STAR,MR> rightL(g), rightR(g);
    DistMatrix<F> AL22T(g), AL22B(g);

    DistMatrix<Base<F>,STAR,STAR> regCand1_STAR_STAR(g), reg1_STAR_STAR(g);

    const Int bsize = Blocksize();
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);
        const Range<Int> ind1( k, k+nb ),
                         ind2Vert( k+nb, m ), ind2Horz( k+nb, n );
        auto AL11 = AL( ind1,     ind1     );
        auto AL21 = AL( ind2Vert, ind1     );
        auto AL22 = AL( ind2Vert, ind2Horz );

        auto regCand1 = regCand( ind1, IR(0,1) );
        auto reg1 = reg( ind1, IR(0,1) );

        AL11_STAR_STAR = AL11; 
        regCand1_STAR_STAR = regCand1;
        reg1_STAR_STAR = reg1;
        RegularizedQSDLDL
        ( AL11_STAR_STAR.Matrix(), pivTol, 
          regCand1_STAR_STAR.LockedMatrix(), reg1_STAR_STAR.Matrix(),
          aPriori );
        GetDiagonal( AL11_STAR_STAR, d1_STAR_STAR );
        AL11 = AL11_STAR_STAR;
        reg1 = reg1_STAR_STAR;

        AL21_VC_STAR.AlignWith( AL22 );
        AL21_VC_STAR = AL21;
        LocalTrsm
        ( RIGHT, LOWER, ADJOINT, UNIT, F(1), AL11_STAR_STAR, AL21_VC_STAR );

        S21Trans_STAR_MC.AlignWith( AL22 );
        Transpose( AL21_VC_STAR, S21Trans_STAR_MC );
        DiagonalSolve( RIGHT, NORMAL, d1_STAR_STAR, AL21_VC_STAR );
        AL21_VR_STAR.AlignWith( AL22 );
        AL21_VR_STAR = AL21_VC_STAR;
        AL21Trans_STAR_MR.AlignWith( AL22 );
        Adjoint( AL21_VR_STAR, AL21Trans_STAR_MR );

        // Partition the update of the bottom-right corner into three pieces
        PartitionRight( S21Trans_STAR_MC, leftL, leftR, AL22.Width() );
        PartitionRight( AL21Trans_STAR_MR, rightL, rightR, AL22.Width() );
        PartitionDown( AL22, AL22T, AL22B, AL22.Width() );
        LocalTrrk( LOWER, ADJOINT, F(-1), leftL, rightL, F(1), AL22T );
        LocalGemm( ADJOINT, NORMAL, F(-1), leftR, rightL, F(1), AL22B );
        LocalTrrk( LOWER, ADJOINT, F(-1), leftR, rightR, F(1), ABR );

        DiagonalSolve( LEFT, NORMAL, d1_STAR_STAR, S21Trans_STAR_MC );
        Transpose( S21Trans_STAR_MC, AL21 );
    }
}

template<typename F>
inline void ProcessFrontSquare
( DistMatrix<F>& AL, DistMatrix<F>& ABR, Base<F> pivTol, 
  const DistMatrix<Base<F>,VC,STAR>& regCand, 
        DistMatrix<Base<F>,VC,STAR>& reg,
  bool aPriori )
{
    DEBUG_ONLY(
        CallStackEntry cse("reg_qsd_ldl::ProcessFrontSquare");
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

    DistMatrix<Base<F>,STAR,STAR> regCand1_STAR_STAR(g), reg1_STAR_STAR(g);

    const Int bsize = Blocksize();
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);
        const Range<Int> ind1( k, k+nb ),
                         ind2Vert( k+nb, m ), ind2Horz( k+nb, n );
        auto AL11 = AL( ind1,     ind1     );
        auto AL21 = AL( ind2Vert, ind1     );
        auto AL22 = AL( ind2Vert, ind2Horz );

        auto regCand1 = regCand( ind1, IR(0,1) );
        auto reg1 = reg( ind1, IR(0,1) );

        AL11_STAR_STAR = AL11; 
        regCand1_STAR_STAR = regCand1;
        reg1_STAR_STAR = reg1;
        RegularizedQSDLDL
        ( AL11_STAR_STAR.Matrix(), pivTol, 
          regCand1_STAR_STAR.LockedMatrix(), reg1_STAR_STAR.Matrix(),
          aPriori );
        GetDiagonal( AL11_STAR_STAR, d1_STAR_STAR );
        AL11 = AL11_STAR_STAR;
        reg1 = reg1_STAR_STAR;

        AL21_VC_STAR.AlignWith( AL22 );
        AL21_VC_STAR = AL21;
        LocalTrsm
        ( RIGHT, LOWER, ADJOINT, UNIT, F(1), AL11_STAR_STAR, AL21_VC_STAR );

        S21Trans_STAR_MC.AlignWith( AL22 );
        Transpose( AL21_VC_STAR, S21Trans_STAR_MC );
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
            Conjugate( AL21Trans_STAR_MR );
        }

        // Partition the update of the bottom-right corner into three pieces
        PartitionRight( S21Trans_STAR_MC, leftL, leftR, AL22.Width() );
        PartitionRight( AL21Trans_STAR_MR, rightL, rightR, AL22.Width() );
        PartitionDown( AL22, AL22T, AL22B, AL22.Width() );
        LocalTrrk( LOWER, ADJOINT, F(-1), leftL, rightL, F(1), AL22T );
        LocalGemm( ADJOINT, NORMAL, F(-1), leftR, rightL, F(1), AL22B );
        LocalTrrk( LOWER, ADJOINT, F(-1), leftR, rightR, F(1), ABR );

        DiagonalSolve( LEFT, NORMAL, d1_STAR_STAR, S21Trans_STAR_MC );
        Transpose( S21Trans_STAR_MC, AL21 );
    }
}

template<typename F> 
inline void ProcessFront
( DistMatrix<F>& AL, DistMatrix<F>& ABR, Base<F> pivTol, 
  const DistMatrix<Base<F>,VC,STAR>& regCand, 
        DistMatrix<Base<F>,VC,STAR>& reg,
  bool aPriori )
{
    DEBUG_ONLY(CallStackEntry cse("reg_qsd_ldl::ProcessFront"))
    const Grid& grid = AL.Grid();
    if( grid.Height() == grid.Width() )
        ProcessFrontSquare( AL, ABR, pivTol, regCand, reg, aPriori );
    else
        ProcessFrontGeneral( AL, ABR, pivTol, regCand, reg, aPriori );
}

} // namespace reg_qsd_ldl
} // namespace El

#endif // ifndef EL_OPTIMIZATION_REGQSDLDL_PROCESSFRONT_HPP
