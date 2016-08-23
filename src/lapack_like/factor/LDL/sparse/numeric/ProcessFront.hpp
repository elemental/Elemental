/*
   Copyright (c) 2009-2012, Jack Poulson, Lexing Ying, and 
   The University of Texas at Austin.
   All rights reserved.

   Copyright (c) 2013, Jack Poulson, Lexing Ying, and Stanford University.
   All rights reserved.

   Copyright (c) 2013-2014, Jack Poulson and 
   The Georgia Institute of Technology.
   All rights reserved.

   Copyright (c) 2014-2015, Jack Poulson and Stanford University.
   All rights reserved.
   
   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_LDL_PROCESSFRONT_HPP
#define EL_LDL_PROCESSFRONT_HPP

namespace El {
namespace ldl {

template<typename F>
void ProcessFrontVanilla( Matrix<F>& AL, Matrix<F>& ABR, bool conjugate )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( ABR.Height() != ABR.Width() )
          LogicError("ABR must be square");
      if( AL.Height() != AL.Width() + ABR.Width() )
          LogicError("AL and ABR don't have conformal dimensions");
    )
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
        const Range<Int> ind1( k, k+nb ), ind2( k+nb, END );
        auto AL11 = AL( ind1, ind1 );
        auto AL21 = AL( ind2, ind1 );
        auto AL22 = AL( ind2, ind2 );

        LDL( AL11, conjugate );
        GetDiagonal( AL11, d1 );

        Trsm( RIGHT, LOWER, orientation, UNIT, F(1), AL11, AL21 );

        S21 = AL21;
        DiagonalSolve( RIGHT, NORMAL, d1, AL21 );

        const Int ind2Size = AL22.Width();
        auto S21B = S21( IR(ind2Size,END), ALL );
        auto AL21T = AL21( IR(0,ind2Size), ALL );
        auto AL21B = AL21( IR(ind2Size,END), ALL );
        Gemm( NORMAL, orientation, F(-1), S21, AL21T, F(1), AL22 );
        MakeTrapezoidal( LOWER, AL22 );
        Trrk( LOWER, NORMAL, orientation, F(-1), S21B, AL21B, F(1), ABR );
    }
}

template<typename F>
void ProcessFrontIntraPiv
( Matrix<F>& AL,
  Matrix<F>& subdiag,
  Permutation& P,
  Matrix<F>& ABR,
  bool conjugate )
{
    DEBUG_CSE
    const Int n = AL.Width();
    const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );

    auto ATL = AL( IR(0,n  ), ALL );
    auto ABL = AL( IR(n,END), ALL );

    LDL( ATL, subdiag, P, conjugate );
    auto diag = GetDiagonal(ATL);

    P.PermuteCols( ABL );
    Trsm( RIGHT, LOWER, orientation, UNIT, F(1), ATL, ABL );
    Matrix<F> SBL( ABL );

    QuasiDiagonalSolve( RIGHT, LOWER, diag, subdiag, ABL, conjugate );
    Trrk( LOWER, NORMAL, orientation, F(-1), SBL, ABL, F(1), ABR );
}

template<typename F>
void ProcessFrontBlock
( Matrix<F>& AL,
  Matrix<F>& ABR,
  bool conjugate,
  bool intraPiv )
{
    DEBUG_CSE
    const Int n = AL.Width();

    auto ATL = AL( IR(0,n  ), ALL );
    auto ABL = AL( IR(n,END), ALL );

    // Make a copy of the original contents of ABL
    Matrix<F> BBL( ABL );

    if( intraPiv )
    {
        Permutation P;
        Matrix<F> dSub;
        // TODO: Expose the pivot type as an option?
        LDL( ATL, dSub, P, conjugate );

        // Solve against ABL and update ABR
        // NOTE: This does not exploit symmetry
        SolveAfter( ATL, dSub, P, ABL, conjugate );
        const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );
        Gemm( NORMAL, orientation, F(-1), ABL, BBL, F(1), ABR );

        // Copy the original contents of ABL back
        ABL = BBL;

        // Finish inverting ATL
        TriangularInverse( LOWER, UNIT, ATL );
        Trdtrmm( LOWER, ATL, dSub, conjugate );
        // TODO: SymmetricPermutation
        MakeSymmetric( LOWER, ATL, conjugate );
        P.PermuteRows( ATL );
        P.PermuteCols( ATL );
    }
    else
    {
        // Call the standard routine
        ProcessFrontVanilla( AL, ABR, conjugate );

        // Copy the original contents of ABL back
        ABL = BBL;

        // Finish inverting ATL
        TriangularInverse( LOWER, UNIT, ATL );
        Trdtrmm( LOWER, ATL, conjugate );
        MakeSymmetric( LOWER, ATL, conjugate );
    }
}

template<typename F>
void ProcessFront( Front<F>& front, LDLFrontType factorType )
{
    DEBUG_CSE
    front.type = factorType;
    DEBUG_ONLY(
      if( front.sparseLeaf )
          LogicError("This should not be possible");
    )
    const bool pivoted = PivotedFactorization( factorType );
    if( BlockFactorization(factorType) )
    {
        ProcessFrontBlock
        ( front.LDense,
          front.workDense,
          front.isHermitian,
          pivoted );
    }
    else if( pivoted )
    {
        ProcessFrontIntraPiv
        ( front.LDense,
          front.subdiag,
          front.p,
          front.workDense, 
          front.isHermitian );
        GetDiagonal( front.LDense, front.diag );
    }
    else
    {
        ProcessFrontVanilla
        ( front.LDense,
          front.workDense,
          front.isHermitian );
        GetDiagonal( front.LDense, front.diag );
    }
}

template<typename F> 
void ProcessFrontVanilla
( DistMatrix<F>& AL,
  DistMatrix<F>& ABR,
  bool conjugate=false )
{
    DEBUG_CSE
    DEBUG_ONLY(
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
        const Range<Int> ind1( k, k+nb ), ind2( k+nb, END );
        auto AL11 = AL( ind1, ind1 );
        auto AL21 = AL( ind2, ind1 );
        auto AL22 = AL( ind2, ind2 );

        AL11_STAR_STAR = AL11; 
        LDL( AL11_STAR_STAR, conjugate );
        GetDiagonal( AL11_STAR_STAR, d1_STAR_STAR );
        AL11 = AL11_STAR_STAR;

        AL21_VC_STAR.AlignWith( AL22 );
        AL21_VC_STAR = AL21;
        LocalTrsm
        ( RIGHT, LOWER, orientation, UNIT, F(1), AL11_STAR_STAR, AL21_VC_STAR );

        S21Trans_STAR_MC.AlignWith( AL22 );
        Transpose( AL21_VC_STAR, S21Trans_STAR_MC );
        DiagonalSolve( RIGHT, NORMAL, d1_STAR_STAR, AL21_VC_STAR );
        AL21Trans_STAR_MR.AlignWith( AL22 );
        Transpose( AL21_VC_STAR, AL21Trans_STAR_MR, conjugate );

        // Partition the update of the bottom-right corner into three pieces
        const Int ind2Size = AL22.Width();
        auto leftL = S21Trans_STAR_MC( ALL, IR(0,ind2Size) );
        auto leftR = S21Trans_STAR_MC( ALL, IR(ind2Size,END) );
        auto rightL = AL21Trans_STAR_MR( ALL, IR(0,ind2Size) );
        auto rightR = AL21Trans_STAR_MR( ALL, IR(ind2Size,END) );
        auto AL22T = AL22( IR(0,ind2Size), ALL );
        auto AL22B = AL22( IR(ind2Size,END), ALL );

        LocalTrrk( LOWER, orientation,  F(-1), leftL, rightL, F(1), AL22T );
        LocalGemm( orientation, NORMAL, F(-1), leftR, rightL, F(1), AL22B );
        LocalTrrk( LOWER, orientation,  F(-1), leftR, rightR, F(1), ABR );

        DiagonalSolve( LEFT, NORMAL, d1_STAR_STAR, S21Trans_STAR_MC );
        Transpose( S21Trans_STAR_MC, AL21 );
    }
}

template<typename F>
void ProcessFrontIntraPiv
( DistMatrix<F>& AL,
  DistMatrix<F,MD,STAR>& subdiag, 
  DistPermutation& P,
  DistMatrix<F>& ABR,
  bool conjugate )
{
    DEBUG_CSE
    const Grid& g = AL.Grid();
    const Int n = AL.Width();
    const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );
    
    auto ATL = AL( IR(0,n  ), ALL );
    auto ABL = AL( IR(n,END), ALL );

    LDL( ATL, subdiag, P, conjugate );
    auto diag = GetDiagonal(ATL);

    P.PermuteCols( ABL );
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
    Transpose( ABL_VR_STAR, ABLTrans_STAR_MR, conjugate );
    LocalTrrk( LOWER, F(-1), SBL_MC_STAR, ABLTrans_STAR_MR, F(1), ABR );
}

template<typename F>
void ProcessFrontBlock
( DistMatrix<F>& AL,
  DistMatrix<F>& ABR,
  bool conjugate,
  bool intraPiv )
{
    DEBUG_CSE
    const Int n = AL.Width();

    auto ATL = AL( IR(0,n  ), ALL );
    auto ABL = AL( IR(n,END), ALL );

    // Make a copy of the original contents of ABL
    DistMatrix<F> BBL( ABL );

    if( intraPiv )
    {
        DistPermutation P( ATL.Grid() );
        DistMatrix<F,MD,STAR> dSub( ATL.Grid() );
        // TODO: Expose the pivot type as an option?
        LDL( ATL, dSub, P, conjugate );

        // Solve against ABL and update ABR
        // NOTE: This update does not exploit symmetry
        SolveAfter( ATL, dSub, P, ABL, conjugate );
        const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );
        Gemm( NORMAL, orientation, F(-1), ABL, BBL, F(1), ABR );

        // Copy the original contents of ABL back
        ABL = BBL;

        // Finish inverting ATL
        TriangularInverse( LOWER, UNIT, ATL );
        Trdtrmm( LOWER, ATL, dSub, conjugate );
        P.InversePermuteSymmetrically( LOWER, ATL, conjugate );
    }
    else
    {
        // Call the standard routine
        ProcessFrontVanilla( AL, ABR, conjugate );

        // Copy the original contents of ABL back
        ABL = BBL;

        // Finish inverting ATL
        TriangularInverse( LOWER, UNIT, ATL );
        Trdtrmm( LOWER, ATL, conjugate );
    }
    MakeSymmetric( LOWER, ATL, conjugate );
}

template<typename F>
void ProcessFront( DistFront<F>& front, LDLFrontType factorType )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( FrontIs1D(front.type) )
          LogicError("Expected front to be in a 2D distribution");
    )
    front.type = factorType;
    const bool pivoted = PivotedFactorization( factorType );
    const Grid& grid = front.L2D.Grid();

    if( BlockFactorization(factorType) )
    {
        ProcessFrontBlock( front.L2D, front.work, front.isHermitian, pivoted );
    }
    else if( pivoted )
    {
        DistMatrix<F,MD,STAR> subdiag(grid);
        front.p.SetGrid( grid );
        ProcessFrontIntraPiv
        ( front.L2D, subdiag, front.p, front.work, front.isHermitian );

        auto diag = GetDiagonal( front.L2D );
        front.diag.SetGrid( grid );
        front.subdiag.SetGrid( grid ); 
        front.diag = diag;
        front.subdiag = subdiag;
    }
    else
    {
        ProcessFrontVanilla( front.L2D, front.work, front.isHermitian );

        auto diag = GetDiagonal( front.L2D );
        front.diag.SetGrid( grid );
        front.diag = diag;
    }
}

} // namespace ldl
} // namespace El

#endif // ifndef EL_LDL_PROCESSFRONT_HPP
