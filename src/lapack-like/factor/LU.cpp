/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

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

        const IndexRange ind1( k, k+nb );
        const IndexRange ind2Vert( k+nb, m );
        const IndexRange ind2Horz( k+nb, n );

        auto A11 = View( A, ind1,     ind1     );
        auto A12 = View( A, ind1,     ind2Horz );
        auto A21 = View( A, ind2Vert, ind1     );
        auto A22 = View( A, ind2Vert, ind2Horz );

        lu::Unb( A11 );
        Trsm( RIGHT, UPPER, NORMAL, NON_UNIT, F(1), A11, A21 );
        Trsm( LEFT, LOWER, NORMAL, UNIT, F(1), A11, A12 );
        Gemm( NORMAL, NORMAL, F(-1), A21, A12, F(1), A22 );
    }
}

template<typename F> 
void LU( AbstractDistMatrix<F>& APre )
{
    DEBUG_ONLY(CallStackEntry cse("LU"))
    const Grid& g = APre.Grid();

    DistMatrix<F> A(g);
    Copy( APre, A, READ_WRITE_PROXY );

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

        const IndexRange ind1( k, k+nb );
        const IndexRange ind2Vert( k+nb, m );
        const IndexRange ind2Horz( k+nb, n );

        auto A11 = View( A, ind1,     ind1     );
        auto A12 = View( A, ind1,     ind2Horz );
        auto A21 = View( A, ind2Vert, ind1     );
        auto A22 = View( A, ind2Vert, ind2Horz );

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
    Copy( A, APre, RESTORE_READ_WRITE_PROXY );
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

        const IndexRange ind0( 0, k    );
        const IndexRange ind1( k, k+nb );
        const IndexRange indB( k, m    );
        const IndexRange ind2Vert( k+nb, m );
        const IndexRange ind2Horz( k+nb, n );

        auto A11 = View( A, ind1,     ind1     );
        auto A12 = View( A, ind1,     ind2Horz );
        auto A21 = View( A, ind2Vert, ind1     );
        auto A22 = View( A, ind2Vert, ind2Horz );

        auto AB0 = View( A, indB, ind0     );
        auto AB1 = View( A, indB, ind1     );
        auto AB2 = View( A, indB, ind2Horz );

        lu::Panel( AB1, p1 );
        PivotsToPartialPermutation( p1, p1Perm, p1InvPerm );
        PermuteRows( AB0, p1Perm, p1InvPerm );
        PermuteRows( AB2, p1Perm, p1InvPerm );

        // Update the preimage of the permutation
        auto pPermB = View( pPerm, indB, IndexRange(0,1) ); 
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

template<typename F> 
void LU( AbstractDistMatrix<F>& APre, AbstractDistMatrix<Int>& pPermPre )
{
    DEBUG_ONLY(
        CallStackEntry cse("LU");
        AssertSameGrids( APre, pPermPre );
    )
    const Int m = APre.Height();
    const Int n = APre.Width();
    const Int minDim = Min(m,n);
    const Int bsize = Blocksize();
    const Grid& g = APre.Grid();

    pPermPre.Resize( m, 1 );

    DistMatrix<F> A(g);
    DistMatrix<Int,VC,STAR> pPerm(g);
    Copy( APre, A, READ_WRITE_PROXY );
    Copy( pPermPre, pPerm, WRITE_PROXY );

    DistMatrix<F,  STAR,STAR> A11_STAR_STAR(g);
    DistMatrix<F,  MC,  STAR> A21_MC_STAR(g);
    DistMatrix<F,  STAR,VR  > A12_STAR_VR(g);
    DistMatrix<F,  STAR,MR  > A12_STAR_MR(g);
    DistMatrix<Int,STAR,STAR> p1_STAR_STAR(g);

    // Initialize the permutation to the identity
    for( Int iLoc=0; iLoc<pPerm.LocalHeight(); ++iLoc )
        pPerm.SetLocal( iLoc, 0, pPerm.GlobalRow(iLoc) );
    DistMatrix<Int,VC,STAR> p1Perm(g), p1InvPerm(g);

    const IndexRange outerInd( 0, n );

    for( Int k=0; k<minDim; k+=bsize )
    {
        const Int nb = Min(bsize,minDim-k);

        const IndexRange ind1( k, k+nb );
        const IndexRange indB( k, m    );
        const IndexRange ind2Vert( k+nb, m );
        const IndexRange ind2Horz( k+nb, n );

        auto A11 = View( A, ind1,     ind1     );
        auto A12 = View( A, ind1,     ind2Horz );
        auto A21 = View( A, ind2Vert, ind1     );
        auto A22 = View( A, ind2Vert, ind2Horz );

        auto AB  = View( A, indB, outerInd );

        A21_MC_STAR.AlignWith( A22 );
        A21_MC_STAR = A21;
        A11_STAR_STAR = A11;

        lu::Panel( A11_STAR_STAR, A21_MC_STAR, p1_STAR_STAR );
        PivotsToPartialPermutation( p1_STAR_STAR, p1Perm, p1InvPerm );
        PermuteRows( AB, p1Perm, p1InvPerm );

        // Update the preimage of the permutation
        auto pPermB = View( pPerm, indB, IndexRange(0,1) );
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
    Copy( A, APre, RESTORE_READ_WRITE_PROXY );
    Copy( pPerm, pPermPre, RESTORE_WRITE_PROXY );
}

template<typename F> 
void LU
( AbstractDistMatrix<F>& A, 
  AbstractDistMatrix<Int>& pPerm, AbstractDistMatrix<Int>& qPerm )
{
    DEBUG_ONLY(CallStackEntry cse("LU"))
    lu::Full( A, pPerm, qPerm );
}

#define PROTO(F) \
  template void LU( Matrix<F>& A ); \
  template void LU( AbstractDistMatrix<F>& A ); \
  template void LU( Matrix<F>& A, Matrix<Int>& pPerm ); \
  template void LU \
  ( AbstractDistMatrix<F>& A, AbstractDistMatrix<Int>& pPerm ); \
  template void LU \
  ( Matrix<F>& A, \
    Matrix<Int>& pPerm, Matrix<Int>& qPerm ); \
  template void LU \
  ( AbstractDistMatrix<F>& A, \
    AbstractDistMatrix<Int>& pPerm, AbstractDistMatrix<Int>& qPerm ); \
  template void LUMod \
  ( Matrix<F>& A, Matrix<Int>& perm, \
    const Matrix<F>& u, const Matrix<F>& v, bool conjugate, \
    Base<F> tau ); \
  template void LUMod \
  ( AbstractDistMatrix<F>& A, AbstractDistMatrix<Int>& perm, \
    const AbstractDistMatrix<F>& u, const AbstractDistMatrix<F>& v, \
    bool conjugate, Base<F> tau ); \
  template void lu::Panel( Matrix<F>& APan, Matrix<Int>& p1 ); \
  template void lu::Panel \
  ( DistMatrix<F,  STAR,STAR>& A11, \
    DistMatrix<F,  MC,  STAR>& A21, \
    DistMatrix<Int,STAR,STAR>& p1 ); \
  template void lu::SolveAfter \
  ( Orientation orientation, const Matrix<F>& A, Matrix<F>& B ); \
  template void lu::SolveAfter \
  ( Orientation orientation, \
    const AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& B ); \
  template void lu::SolveAfter \
  ( Orientation orientation, const Matrix<F>& A, \
    const Matrix<Int>& pPerm, Matrix<F>& B ); \
  template void lu::SolveAfter \
  ( Orientation orientation, const AbstractDistMatrix<F>& A, \
    const AbstractDistMatrix<Int>& pPerm, AbstractDistMatrix<F>& B ); \
  template void lu::SolveAfter \
  ( Orientation orientation, const Matrix<F>& A, \
    const Matrix<Int>& pPerm, \
    const Matrix<Int>& qPerm, Matrix<F>& B ); \
  template void lu::SolveAfter \
  ( Orientation orientation, const AbstractDistMatrix<F>& A, \
    const AbstractDistMatrix<Int>& pPerm, \
    const AbstractDistMatrix<Int>& qPerm, AbstractDistMatrix<F>& B );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
