/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

#include "./LU/Local.hpp"
#include "./LU/Panel.hpp"
#include "./LU/Full.hpp"
#include "./LU/Mod.hpp"
#include "./LU/SolveAfter.hpp"

namespace El {

// Performs LU factorization without pivoting

template<typename F> 
void LU( Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("LU"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);
    const Int bsize = Blocksize();
    for( Int k=0; k<minDim; k+=bsize )
    {
        const Int nb = Min(bsize,minDim-k);
        auto A11 = ViewRange( A, k,    k,    k+nb, k+nb );   
        auto A12 = ViewRange( A, k,    k+nb, k+nb, n    );
        auto A21 = ViewRange( A, k+nb, k,    m,    k+nb );
        auto A22 = ViewRange( A, k+nb, k+nb, m,    n    );

        lu::Unb( A11 );
        Trsm( RIGHT, UPPER, NORMAL, NON_UNIT, F(1), A11, A21 );
        Trsm( LEFT, LOWER, NORMAL, UNIT, F(1), A11, A12 );
        Gemm( NORMAL, NORMAL, F(-1), A21, A12, F(1), A22 );
    }
}

template<typename F> 
void LU( DistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("LU"))
    const Grid& g = A.Grid();
    DistMatrix<F,STAR,STAR> A11_STAR_STAR(g);
    DistMatrix<F,MC,  STAR> A21_MC_STAR(g);
    DistMatrix<F,STAR,VR  > A12_STAR_VR(g);
    DistMatrix<F,STAR,MR  > A12_STAR_MR(g);

    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);
    const Int bsize = Blocksize();
    for( Int k=0; k<minDim; k+=bsize )
    {
        const Int nb = Min(bsize,minDim-k);
        auto A11 = ViewRange( A, k,    k,    k+nb, k+nb );          
        auto A12 = ViewRange( A, k,    k+nb, k+nb, n    );
        auto A21 = ViewRange( A, k+nb, k,    m,    k+nb );
        auto A22 = ViewRange( A, k+nb, k+nb, m,    n    );

        A11_STAR_STAR = A11;
        LocalLU( A11_STAR_STAR );
        A11 = A11_STAR_STAR;

        A21_MC_STAR.AlignWith( A22 );
        A21_MC_STAR = A21;
        LocalTrsm
        ( RIGHT, UPPER, NORMAL, NON_UNIT, F(1), A11_STAR_STAR, A21_MC_STAR );
        A21 = A21_MC_STAR;

        // Perhaps we should give up perfectly distributing this operation since
        // it's total contribution is only O(n^2)
        A12_STAR_VR.AlignWith( A22 );
        A12_STAR_VR = A12;
        LocalTrsm
        ( LEFT, LOWER, NORMAL, UNIT, F(1), A11_STAR_STAR, A12_STAR_VR );

        A12_STAR_MR.AlignWith( A22 );
        A12_STAR_MR = A12_STAR_VR;
        LocalGemm( NORMAL, NORMAL, F(-1), A21_MC_STAR, A12_STAR_MR, F(1), A22 );
        A12 = A12_STAR_MR;
    }
}

// Performs LU factorization with partial pivoting

template<typename F> 
void LU( Matrix<F>& A, Matrix<Int>& pPerm )
{
    DEBUG_ONLY(CallStackEntry cse("LU"))

    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);
    const Int bsize = Blocksize();

    // Initialize P to the identity matrix
    pPerm.Resize( m, 1 );
    for( Int i=0; i<m; ++i )
        pPerm.Set( i, 0, i );

    // Temporaries for accumulating partial permutations for each block
    Matrix<Int> p1;
    Matrix<Int> p1Perm, p1InvPerm;

    for( Int k=0; k<minDim; k+=bsize )
    {
        const Int nb = Min(bsize,minDim-k);
        auto A11  = ViewRange( A, k,    k,    k+nb, k+nb );
        auto A12  = ViewRange( A, k,    k+nb, k+nb, n    );
        auto A21  = ViewRange( A, k+nb, k,    m,    k+nb );
        auto A22  = ViewRange( A, k+nb, k+nb, m,    n    );
        auto ABL  = ViewRange( A, k,    0,    m,    k    );
        auto ABRL = ViewRange( A, k,    k,    m,    k+nb );
        auto ABRR = ViewRange( A, k,    k+nb, m,    n    );

        lu::Panel( ABRL, p1 );
        PivotsToPartialPermutation( p1, p1Perm, p1InvPerm );
        PermuteRows( ABL, p1Perm, p1InvPerm );
        PermuteRows( ABRR, p1Perm, p1InvPerm );

        // Update the preimage of the permutation
        auto pPermB = ViewRange( pPerm, k, 0, m, 1 ); 
        PermuteRows( pPermB, p1Perm, p1InvPerm );

        Trsm( LEFT, LOWER, NORMAL, UNIT, F(1), A11, A12 );
        Gemm( NORMAL, NORMAL, F(-1), A21, A12, F(1), A22 );
    }
}

template<typename F> 
void LU( Matrix<F>& A, Matrix<Int>& pPerm, Matrix<Int>& qPerm )
{
    DEBUG_ONLY(CallStackEntry cse("LU"))
    lu::Full( A, pPerm, qPerm );
}

template<typename F,Dist UPerm> 
void LU( DistMatrix<F>& A, DistMatrix<Int,UPerm,STAR>& pPerm )
{
    DEBUG_ONLY(
        CallStackEntry cse("LU");
        if( A.Grid() != pPerm.Grid() )
            LogicError("{A,pPerm} must be distributed over the same grid");
    )
    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);
    const Int bsize = Blocksize();
    const Grid& g = A.Grid();

    DistMatrix<F,  STAR,STAR> A11_STAR_STAR(g);
    DistMatrix<F,  MC,  STAR> A21_MC_STAR(g);
    DistMatrix<F,  STAR,VR  > A12_STAR_VR(g);
    DistMatrix<F,  STAR,MR  > A12_STAR_MR(g);
    DistMatrix<Int,STAR,STAR> p1_STAR_STAR(g);

    // Initialize the permutation to the identity
    pPerm.Resize( m, 1 );
    for( Int iLoc=0; iLoc<pPerm.LocalHeight(); ++iLoc )
        pPerm.SetLocal( iLoc, 0, pPerm.GlobalRow(iLoc) );
    DistMatrix<Int,UPerm,STAR> p1Perm(g), p1InvPerm(g);

    for( Int k=0; k<minDim; k+=bsize )
    {
        const Int nb = Min(bsize,minDim-k);
        auto A11  = ViewRange( A, k,    k,    k+nb, k+nb );
        auto A12  = ViewRange( A, k,    k+nb, k+nb, n    );
        auto A21  = ViewRange( A, k+nb, k,    m,    k+nb );
        auto A22  = ViewRange( A, k+nb, k+nb, m,    n    );
        auto AB   = ViewRange( A, k,    0,    m,    n    );
        auto ABL  = ViewRange( A, k,    0,    m,    k    );
        auto ABRL = ViewRange( A, k,    k,    m,    k+nb );
        auto ABRR = ViewRange( A, k,    k+nb, m,    n    );

        A21_MC_STAR.AlignWith( A22 );
        A21_MC_STAR = A21;
        A11_STAR_STAR = A11;

        lu::Panel( A11_STAR_STAR, A21_MC_STAR, p1_STAR_STAR );
        PivotsToPartialPermutation( p1_STAR_STAR, p1Perm, p1InvPerm );
        PermuteRows( AB, p1Perm, p1InvPerm );

        // Update the preimage of the permutation
        auto pPermB = ViewRange( pPerm, k, 0, m, 1 ); 
        PermuteRows( pPermB, p1Perm, p1InvPerm );

        // Perhaps we should give up perfectly distributing this operation since
        // it's total contribution is only O(n^2)
        A12_STAR_VR.AlignWith( A22 );
        A12_STAR_VR = A12;
        LocalTrsm
        ( LEFT, LOWER, NORMAL, UNIT, F(1), A11_STAR_STAR, A12_STAR_VR );

        A12_STAR_MR.AlignWith( A22 );
        A12_STAR_MR = A12_STAR_VR;
        LocalGemm( NORMAL, NORMAL, F(-1), A21_MC_STAR, A12_STAR_MR, F(1), A22 );

        A11 = A11_STAR_STAR;
        A12 = A12_STAR_MR;
        A21 = A21_MC_STAR;
    }
}

template<typename F,Dist UPerm> 
void LU
( DistMatrix<F>& A, 
  DistMatrix<Int,UPerm,STAR>& pPerm, DistMatrix<Int,UPerm,STAR>& qPerm )
{
    DEBUG_ONLY(CallStackEntry cse("LU"))
    lu::Full( A, pPerm, qPerm );
}

#define PROTO(F) \
  template void LU( Matrix<F>& A ); \
  template void LU( DistMatrix<F>& A ); \
  template void LU( Matrix<F>& A, Matrix<Int>& pPerm ); \
  template void LU( DistMatrix<F>& A, DistMatrix<Int,VC,STAR>& pPerm ); \
  template void LU \
  ( Matrix<F>& A, \
    Matrix<Int>& pPerm, Matrix<Int>& qPerm ); \
  template void LU \
  ( DistMatrix<F>& A, \
    DistMatrix<Int,VC,STAR>& pPerm, DistMatrix<Int,VC,STAR>& qPerm ); \
  template void LUMod \
  ( Matrix<F>& A, Matrix<Int>& perm, \
    const Matrix<F>& u, const Matrix<F>& v, bool conjugate, \
    Base<F> tau ); \
  template void LUMod \
  ( DistMatrix<F>& A, DistMatrix<Int,VC,STAR>& perm, \
    const DistMatrix<F>& u, const DistMatrix<F>& v, bool conjugate, \
    Base<F> tau ); \
  template void lu::Panel( Matrix<F>& APan, Matrix<Int>& p1 ); \
  template void lu::Panel \
  ( DistMatrix<F,  STAR,STAR>& A11, \
    DistMatrix<F,  MC,  STAR>& A21, \
    DistMatrix<Int,STAR,STAR>& p1 ); \
  template void lu::SolveAfter \
  ( Orientation orientation, const Matrix<F>& A, Matrix<F>& B ); \
  template void lu::SolveAfter \
  ( Orientation orientation, const DistMatrix<F>& A, DistMatrix<F>& B ); \
  template void lu::SolveAfter \
  ( Orientation orientation, const Matrix<F>& A, \
    const Matrix<Int>& pPerm, Matrix<F>& B ); \
  template void lu::SolveAfter \
  ( Orientation orientation, const DistMatrix<F>& A, \
    const DistMatrix<Int,VC,STAR>& pPerm, DistMatrix<F>& B ); \
  template void lu::SolveAfter \
  ( Orientation orientation, const Matrix<F>& A, \
    const Matrix<Int>& pPerm, \
    const Matrix<Int>& qPerm, Matrix<F>& B ); \
  template void lu::SolveAfter \
  ( Orientation orientation, const DistMatrix<F>& A, \
    const DistMatrix<Int,VC,STAR>& pPerm, \
    const DistMatrix<Int,VC,STAR>& qPerm, DistMatrix<F>& B );

PROTO(float)
PROTO(double)
PROTO(Complex<float>)
PROTO(Complex<double>)

} // namespace El
