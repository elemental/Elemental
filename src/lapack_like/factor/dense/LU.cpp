/*
   Copyright (c) 2009-2015, Jack Poulson
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

        const Range<Int> ind1( k, k+nb ),
                         ind2Vert( k+nb, m ),
                         ind2Horz( k+nb, n );

        auto A11 = A( ind1,     ind1     );
        auto A12 = A( ind1,     ind2Horz );
        auto A21 = A( ind2Vert, ind1     );
        auto A22 = A( ind2Vert, ind2Horz );

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

    auto APtr = ReadWriteProxy<F,MC,MR>( &APre );
    auto& A = *APtr;

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

        const Range<Int> ind1( k, k+nb ),
                         ind2Vert( k+nb, m ),
                         ind2Horz( k+nb, n );

        auto A11 = A( ind1,     ind1     );
        auto A12 = A( ind1,     ind2Horz );
        auto A21 = A( ind2Vert, ind1     );
        auto A22 = A( ind2Vert, ind2Horz );

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
void LU( Matrix<F>& A, Matrix<Int>& p )
{
    DEBUG_ONLY(CallStackEntry cse("LU"))

    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);
    const Int bsize = Blocksize();

    // Initialize P to the identity matrix
    p.Resize( m, 1 );
    for( Int i=0; i<m; ++i )
        p.Set( i, 0, i );

    // Temporaries for accumulating partial permutations for each block
    Matrix<Int> p1Piv;
    Matrix<Int> p1, p1Inv;

    for( Int k=0; k<minDim; k+=bsize )
    {
        const Int nb = Min(bsize,minDim-k);

        const Range<Int> ind0( 0, k    ),
                         ind1( k, k+nb ),
                         indB( k, m    ),
                         ind2Vert( k+nb, m ),
                         ind2Horz( k+nb, n );

        auto A11 = A( ind1,     ind1     );
        auto A12 = A( ind1,     ind2Horz );
        auto A21 = A( ind2Vert, ind1     );
        auto A22 = A( ind2Vert, ind2Horz );

        auto AB0 = A( indB, ind0     );
        auto AB1 = A( indB, ind1     );
        auto AB2 = A( indB, ind2Horz );

        lu::Panel( AB1, p1Piv );
        PivotsToPartialPermutation( p1Piv, p1, p1Inv );
        PermuteRows( AB0, p1, p1Inv );
        PermuteRows( AB2, p1, p1Inv );

        // Update the preimage of the permutation
        auto pB = p( indB, IR(0,1) ); 
        PermuteRows( pB, p1, p1Inv );

        Trsm( LEFT, LOWER, NORMAL, UNIT, F(1), A11, A12 );
        Gemm( NORMAL, NORMAL, F(-1), A21, A12, F(1), A22 );
    }
}

template<typename F> 
void LU( Matrix<F>& A, Matrix<Int>& p, Matrix<Int>& q )
{
    DEBUG_ONLY(CallStackEntry cse("LU"))
    lu::Full( A, p, q );
}

template<typename F> 
void LU( AbstractDistMatrix<F>& APre, AbstractDistMatrix<Int>& pPre )
{
    DEBUG_ONLY(
        CallStackEntry cse("LU");
        AssertSameGrids( APre, pPre );
    )

    auto APtr = ReadWriteProxy<F,MC,MR>( &APre ); auto& A = *APtr;
    auto pPtr = WriteProxy<Int,VC,STAR>( &pPre ); auto& p = *pPtr;

    const Grid& g = A.Grid();
    DistMatrix<F,  STAR,STAR> A11_STAR_STAR(g);
    DistMatrix<F,  MC,  STAR> A21_MC_STAR(g);
    DistMatrix<F,  STAR,VR  > A12_STAR_VR(g);
    DistMatrix<F,  STAR,MR  > A12_STAR_MR(g);
    DistMatrix<Int,STAR,STAR> p1Piv_STAR_STAR(g);

    // Initialize the permutation to the identity
    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);
    p.Resize( m, 1 );
    for( Int iLoc=0; iLoc<p.LocalHeight(); ++iLoc )
        p.SetLocal( iLoc, 0, p.GlobalRow(iLoc) );

    DistMatrix<Int,VC,STAR> p1(g), p1Inv(g);

    const Range<Int> outerInd( 0, n );

    const Int bsize = Blocksize();
    for( Int k=0; k<minDim; k+=bsize )
    {
        const Int nb = Min(bsize,minDim-k);

        const Range<Int> ind1( k, k+nb ),
                         indB( k, m    ),
                         ind2Vert( k+nb, m ),
                         ind2Horz( k+nb, n );

        auto A11 = A( ind1,     ind1     );
        auto A12 = A( ind1,     ind2Horz );
        auto A21 = A( ind2Vert, ind1     );
        auto A22 = A( ind2Vert, ind2Horz );

        auto AB  = A( indB, outerInd );

        A21_MC_STAR.AlignWith( A22 );
        A21_MC_STAR = A21;
        A11_STAR_STAR = A11;

        lu::Panel( A11_STAR_STAR, A21_MC_STAR, p1Piv_STAR_STAR );
        PivotsToPartialPermutation( p1Piv_STAR_STAR, p1, p1Inv );
        PermuteRows( AB, p1, p1Inv );

        // Update the preimage of the permutation
        auto pB = p( indB, IR(0,1) );
        PermuteRows( pB, p1, p1Inv );

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

template<typename F> 
void LU
( AbstractDistMatrix<F>& A, 
  AbstractDistMatrix<Int>& p, AbstractDistMatrix<Int>& q )
{
    DEBUG_ONLY(CallStackEntry cse("LU"))
    lu::Full( A, p, q );
}

#define PROTO(F) \
  template void LU( Matrix<F>& A ); \
  template void LU( AbstractDistMatrix<F>& A ); \
  template void LU( Matrix<F>& A, Matrix<Int>& p ); \
  template void LU( AbstractDistMatrix<F>& A, AbstractDistMatrix<Int>& p ); \
  template void LU( Matrix<F>& A, Matrix<Int>& p, Matrix<Int>& q ); \
  template void LU \
  ( AbstractDistMatrix<F>& A, \
    AbstractDistMatrix<Int>& p, AbstractDistMatrix<Int>& q ); \
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
    const Matrix<Int>& p, Matrix<F>& B ); \
  template void lu::SolveAfter \
  ( Orientation orientation, const AbstractDistMatrix<F>& A, \
    const AbstractDistMatrix<Int>& p, AbstractDistMatrix<F>& B ); \
  template void lu::SolveAfter \
  ( Orientation orientation, const Matrix<F>& A, \
    const Matrix<Int>& p, const Matrix<Int>& q, Matrix<F>& B ); \
  template void lu::SolveAfter \
  ( Orientation orientation, const AbstractDistMatrix<F>& A, \
    const AbstractDistMatrix<Int>& p, const AbstractDistMatrix<Int>& q, \
          AbstractDistMatrix<F>& B );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
