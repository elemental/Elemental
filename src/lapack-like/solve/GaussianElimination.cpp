/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

namespace El {

// Short-circuited form of LU factorization with partial pivoting
template<typename F> 
inline void
RowEchelon( Matrix<F>& A, Matrix<F>& B )
{
    DEBUG_ONLY(
        CallStackEntry cse("RowEchelon");
        if( A.Height() != B.Height() )
            LogicError("A and B must be the same height");
    )

    Matrix<Int> p1;
    Matrix<Int> p1Perm, p1InvPerm;

    const Int mA = A.Height();
    const Int nA = A.Width();
    const Int minDimA = Min(mA,nA);
    const Int nB = B.Width();
    const Int bsize = Blocksize();
    for( Int k=0; k<minDimA; k+=bsize )
    {
        const Int nb = Min(bsize,minDimA-k);
        auto A11  = ViewRange( A, k,    k,    k+nb, k+nb );
        auto A12  = ViewRange( A, k,    k+nb, k+nb, nA   );
        auto A21  = ViewRange( A, k+nb, k,    mA,   k+nb );
        auto A22  = ViewRange( A, k+nb, k+nb, mA,   nA   ); 
        auto APan = ViewRange( A, k,    k+nb, mA,   nA   );
        auto B1   = ViewRange( B, k,    0,    k+nb, nB   );
        auto B2   = ViewRange( B, k+nb, 0,    mA,   nB   );
        auto BB   = ViewRange( B, k,    0,    mA,   nB   );

        lu::Panel( APan, p1 );
        PivotsToPartialPermutation( p1, p1Perm, p1InvPerm ); 
        PermuteRows( BB, p1Perm, p1InvPerm );

        Trsm( LEFT, LOWER, NORMAL, UNIT, F(1), A11, A12 );
        Trsm( LEFT, LOWER, NORMAL, UNIT, F(1), A11, B1 );

        Gemm( NORMAL, NORMAL, F(-1), A21, A12, F(1), A22 );
        Gemm( NORMAL, NORMAL, F(-1), A21, B1,  F(1), B2 );
    }
}

// Short-circuited form of LU factorization with partial pivoting
template<typename F> 
inline void
RowEchelon( DistMatrix<F>& A, DistMatrix<F>& B )
{
    DEBUG_ONLY(
        CallStackEntry cse("RowEchelon");
        if( A.Grid() != B.Grid() )
            LogicError("{A,B} must be distributed over the same grid");
        if( A.Height() != B.Height() )
            LogicError("A and B must be the same height");
    )
    const Int mA = A.Height();
    const Int nA = A.Width();
    const Int minDimA = Min(mA,nA);
    const Int nB = B.Width();
    const Int bsize = Blocksize();
    const Grid& g = A.Grid();

    DistMatrix<F,STAR,STAR> A11_STAR_STAR(g);
    DistMatrix<F,STAR,VR  > A12_STAR_VR(g);
    DistMatrix<F,STAR,MR  > A12_STAR_MR(g);
    DistMatrix<F,MC,  STAR> A21_MC_STAR(g);
    DistMatrix<F,STAR,VR  > B1_STAR_VR(g);
    DistMatrix<F,STAR,MR  > B1_STAR_MR(g);
    DistMatrix<Int,STAR,STAR> p1_STAR_STAR(g);

    DistMatrix<Int,VC,STAR> p1Perm(g), p1InvPerm(g);

    // In case B's columns are not aligned with A's
    const bool BAligned = ( B.ColShift() == A.ColShift() );
    DistMatrix<F,MC,STAR> A21_MC_STAR_B(g);

    for( Int k=0; k<minDimA; k+=bsize )
    {
        const Int nb = Min(bsize,minDimA-k);
        auto A11  = ViewRange( A, k,    k,    k+nb, k+nb );
        auto A12  = ViewRange( A, k,    k+nb, k+nb, nA   );
        auto A21  = ViewRange( A, k+nb, k,    mA,   k+nb );
        auto A22  = ViewRange( A, k+nb, k+nb, mA,   nA   );          
        auto APan = ViewRange( A, k,    k+nb, mA,   nA   );
        auto B1   = ViewRange( B, k,    0,    k+nb, nB   );
        auto B2   = ViewRange( B, k+nb, 0,    mA,   nB   );
        auto BB   = ViewRange( B, k,    0,    mA,   nB   );

        A11_STAR_STAR = A11;
        A21_MC_STAR.AlignWith( A22 );
        A21_MC_STAR = A21;

        lu::Panel( A11_STAR_STAR, A21_MC_STAR, p1_STAR_STAR );
        PivotsToPartialPermutation( p1_STAR_STAR, p1Perm, p1InvPerm );
        PermuteRows( APan, p1Perm, p1InvPerm );
        PermuteRows( BB,   p1Perm, p1InvPerm );

        A12_STAR_VR.AlignWith( A22 );
        A12_STAR_VR = A12;
        B1_STAR_VR.AlignWith( B1 );
        B1_STAR_VR = B1;
        LocalTrsm
        ( LEFT, LOWER, NORMAL, UNIT, F(1), A11_STAR_STAR, A12_STAR_VR );
        LocalTrsm( LEFT, LOWER, NORMAL, UNIT, F(1), A11_STAR_STAR, B1_STAR_VR );

        A12_STAR_MR.AlignWith( A22 );
        A12_STAR_MR = A12_STAR_VR;
        B1_STAR_MR.AlignWith( B1 );
        B1_STAR_MR = B1_STAR_VR;
        LocalGemm( NORMAL, NORMAL, F(-1), A21_MC_STAR, A12_STAR_MR, F(1), A22 );
        if( BAligned )
        {
            LocalGemm
            ( NORMAL, NORMAL, F(-1), A21_MC_STAR, B1_STAR_MR, F(1), B2 );
        }
        else
        {
            A21_MC_STAR_B.AlignWith( B2 );
            A21_MC_STAR_B = A21_MC_STAR;
            LocalGemm
            ( NORMAL, NORMAL, F(-1), A21_MC_STAR_B, B1_STAR_MR, F(1), B2 );
        }

        A11 = A11_STAR_STAR;
        A12 = A12_STAR_MR;
        B1 = B1_STAR_MR;
    }
}

template<typename F> 
void GaussianElimination( Matrix<F>& A, Matrix<F>& B )
{
    DEBUG_ONLY(
        CallStackEntry cse("GaussianElimination");
        if( A.Height() != A.Width() )
            LogicError("A must be square");
        if( A.Height() != B.Height() )
            LogicError("A and B must be the same height");
    )
    RowEchelon( A, B );
    Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), A, B );
}

template<typename F> 
void GaussianElimination( DistMatrix<F>& A, DistMatrix<F>& B )
{
    DEBUG_ONLY(
        CallStackEntry cse("GaussianElimination");
        if( A.Grid() != B.Grid() )
            LogicError("{A,B} must be distributed over the same grid");
        if( A.Height() != A.Width() )
            LogicError("A must be square");
        if( A.Height() != B.Height() )
            LogicError("A and B must be the same height");
    )
    RowEchelon( A, B );
    Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), A, B );
}

#define PROTO(F) \
  template void GaussianElimination( Matrix<F>& A, Matrix<F>& B ); \
  template void GaussianElimination( DistMatrix<F>& A, DistMatrix<F>& B );

PROTO(float)
PROTO(double)
PROTO(Complex<float>)
PROTO(Complex<double>)

} // namespace El
