/*
   Copyright (c) 2009-2016, Jack Poulson
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
    DEBUG_ONLY(CSE cse("LU"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);
    const Int bsize = Blocksize();
    for( Int k=0; k<minDim; k+=bsize )
    {
        const Int nb = Min(bsize,minDim-k);
        const IR ind1( k, k+nb ), ind2( k+nb, END );

        auto A11 = A( ind1, ind1 );
        auto A12 = A( ind1, ind2 );
        auto A21 = A( ind2, ind1 );
        auto A22 = A( ind2, ind2 );

        lu::Unb( A11 );
        Trsm( RIGHT, UPPER, NORMAL, NON_UNIT, F(1), A11, A21 );
        Trsm( LEFT, LOWER, NORMAL, UNIT, F(1), A11, A12 );
        Gemm( NORMAL, NORMAL, F(-1), A21, A12, F(1), A22 );
    }
}

template<typename F> 
void LU( ElementalMatrix<F>& APre )
{
    DEBUG_ONLY(CSE cse("LU"))

    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre );
    auto& A = AProx.Get();

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
        const IR ind1( k, k+nb ), ind2( k+nb, END );

        auto A11 = A( ind1, ind1 );
        auto A12 = A( ind1, ind2 );
        auto A21 = A( ind2, ind1 );
        auto A22 = A( ind2, ind2 );

        A11_STAR_STAR = A11;
        LU( A11_STAR_STAR );
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

template<typename F> 
void LU( DistMatrix<F,STAR,STAR>& A )
{ LU( A.Matrix() ); }

// Performs LU factorization with partial pivoting

template<typename F> 
void LU( Matrix<F>& A, Permutation& P )
{
    DEBUG_ONLY(CSE cse("LU"))

    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);
    const Int bsize = Blocksize();

    P.MakeIdentity( m );
    P.ReserveSwaps( minDim );
    
    Permutation PB;
   
    // Temporaries for accumulating partial permutations for each block
    for( Int k=0; k<minDim; k+=bsize )
    {
        const Int nb = Min(bsize,minDim-k);
        const IR ind0( 0, k ), ind1( k, k+nb ), ind2( k+nb, END ), 
                 indB( k, END );

        auto A11 = A( ind1, ind1 );
        auto A12 = A( ind1, ind2 );
        auto A21 = A( ind2, ind1 );
        auto A22 = A( ind2, ind2 );

        auto AB0 = A( indB, ind0 );
        auto AB1 = A( indB, ind1 );
        auto AB2 = A( indB, ind2 );

        lu::Panel( AB1, P, PB, k );
        PB.PermuteRows( AB0 );
        PB.PermuteRows( AB2 );

        Trsm( LEFT, LOWER, NORMAL, UNIT, F(1), A11, A12 );
        Gemm( NORMAL, NORMAL, F(-1), A21, A12, F(1), A22 );
    }
}

template<typename F> 
void LU
( Matrix<F>& A,
  Permutation& P,
  Permutation& Q )
{
    DEBUG_ONLY(CSE cse("LU"))
    lu::Full( A, P, Q );
}

template<typename F> 
void LU( ElementalMatrix<F>& APre, DistPermutation& P )
{
    DEBUG_ONLY(CSE cse("LU"))

    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre );
    auto& A = AProx.Get();

    const Grid& g = A.Grid();
    DistMatrix<F,  STAR,STAR> A11_STAR_STAR(g);
    DistMatrix<F,  MC,  STAR> A21_MC_STAR(g);
    DistMatrix<F,  STAR,VR  > A12_STAR_VR(g);
    DistMatrix<F,  STAR,MR  > A12_STAR_MR(g);

    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);
    P.SetGrid( g );

    P.MakeIdentity( m );
    P.ReserveSwaps( minDim );

    DistPermutation PB(g);

    vector<F> panelBuf, pivotBuf;
    const Int bsize = Blocksize();
    for( Int k=0; k<minDim; k+=bsize )
    {
        const Int nb = Min(bsize,minDim-k);
        const IR ind1( k, k+nb ), ind2( k+nb, END ), indB( k, END );

        auto A11 = A( ind1, ind1 );
        auto A12 = A( ind1, ind2 );
        auto A21 = A( ind2, ind1 );
        auto A22 = A( ind2, ind2 );

        auto AB  = A( indB, ALL );

        const Int A21Height = A21.Height();
        const Int A21LocHeight = A21.LocalHeight();
        const Int panelLDim = nb+A21LocHeight;
        FastResize( panelBuf, panelLDim*nb );
        A11_STAR_STAR.Attach
        ( nb, nb, g, 0, 0, &panelBuf[0], panelLDim, 0 );
        A21_MC_STAR.Attach
        ( A21Height, nb, g, A21.ColAlign(), 0, &panelBuf[nb], panelLDim, 0 );
        A11_STAR_STAR = A11;
        A21_MC_STAR = A21;
        lu::Panel( A11_STAR_STAR, A21_MC_STAR, P, PB, k, pivotBuf );

        PB.PermuteRows( AB );

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
( ElementalMatrix<F>& A, 
  DistPermutation& P,
  DistPermutation& Q )
{
    DEBUG_ONLY(CSE cse("LU"))
    lu::Full( A, P, Q );
}

#define PROTO(F) \
  template void LU( Matrix<F>& A ); \
  template void LU( ElementalMatrix<F>& A ); \
  template void LU( DistMatrix<F,STAR,STAR>& A ); \
  template void LU \
  ( Matrix<F>& A, \
    Permutation& P ); \
  template void LU \
  ( ElementalMatrix<F>& A, \
    DistPermutation& P ); \
  template void LU \
  ( Matrix<F>& A, \
    Permutation& P, \
    Permutation& Q ); \
  template void LU \
  ( ElementalMatrix<F>& A, \
    DistPermutation& P, \
    DistPermutation& Q ); \
  template void LUMod \
  (       Matrix<F>& A, \
          Permutation& P, \
    const Matrix<F>& u, \
    const Matrix<F>& v, \
    bool conjugate, \
    Base<F> tau ); \
  template void LUMod \
  (       ElementalMatrix<F>& A, \
          DistPermutation& P, \
    const ElementalMatrix<F>& u, \
    const ElementalMatrix<F>& v, \
    bool conjugate, \
    Base<F> tau ); \
  template void lu::Panel \
  ( Matrix<F>& APan, \
    Permutation& P, \
    Permutation& PB, \
    Int offset ); \
  template void lu::Panel \
  ( DistMatrix<F,  STAR,STAR>& A11, \
    DistMatrix<F,  MC,  STAR>& A21, \
    DistPermutation& P, \
    DistPermutation& PB, \
    Int offset, \
    vector<F>& pivotBuf ); \
  template void lu::SolveAfter \
  ( Orientation orientation, \
    const Matrix<F>& A, \
          Matrix<F>& B ); \
  template void lu::SolveAfter \
  ( Orientation orientation, \
    const ElementalMatrix<F>& A, \
          ElementalMatrix<F>& B ); \
  template void lu::SolveAfter \
  ( Orientation orientation, \
    const Matrix<F>& A, \
    const Permutation& P, \
          Matrix<F>& B ); \
  template void lu::SolveAfter \
  ( Orientation orientation, \
    const ElementalMatrix<F>& A, \
    const DistPermutation& P, \
          ElementalMatrix<F>& B ); \
  template void lu::SolveAfter \
  ( Orientation orientation, \
    const Matrix<F>& A, \
    const Permutation& P, \
    const Permutation& Q, \
          Matrix<F>& B ); \
  template void lu::SolveAfter \
  ( Orientation orientation, \
    const ElementalMatrix<F>& A, \
    const DistPermutation& P, \
    const DistPermutation& Q, \
          ElementalMatrix<F>& B );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

} // namespace El
