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
void LU( Matrix<F>& A, Matrix<Int>& rowPiv )
{
    DEBUG_ONLY(CSE cse("LU"))

    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);
    const Int bsize = Blocksize();

    rowPiv.Resize( m, 1 );

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

        auto p1Piv = rowPiv( ind1, ALL );

        lu::Panel( AB1, p1Piv );
        ApplyRowPivots( AB0, p1Piv );
        ApplyRowPivots( AB2, p1Piv );
        Shift( p1Piv, k );

        Trsm( LEFT, LOWER, NORMAL, UNIT, F(1), A11, A12 );
        Gemm( NORMAL, NORMAL, F(-1), A21, A12, F(1), A22 );
    }
}

template<typename F> 
void LU( Matrix<F>& A, Matrix<Int>& rowPiv, Matrix<Int>& colPiv )
{
    DEBUG_ONLY(CSE cse("LU"))
    lu::Full( A, rowPiv, colPiv );
}

template<typename F> 
void LU( ElementalMatrix<F>& APre, ElementalMatrix<Int>& rowPivPre )
{
    DEBUG_ONLY(
      CSE cse("LU");
      AssertSameGrids( APre, rowPivPre );
    )

    auto APtr = ReadWriteProxy<F,MC,MR>( &APre );
    auto& A = *APtr;

    auto rowPivPtr = WriteProxy<Int,STAR,STAR>( &rowPivPre );
    auto& rowPiv = *rowPivPtr;

    const Grid& g = A.Grid();
    DistMatrix<F,  STAR,STAR> A11_STAR_STAR(g);
    DistMatrix<F,  MC,  STAR> A21_MC_STAR(g);
    DistMatrix<F,  STAR,VR  > A12_STAR_VR(g);
    DistMatrix<F,  STAR,MR  > A12_STAR_MR(g);

    // Initialize the permutation to the identity
    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);
    rowPiv.Resize( minDim, 1 );

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

        auto p1Piv = rowPiv( ind1, ALL );

        const Int A21Height = A21.Height();
        const Int A21LocHeight = A21.LocalHeight();
        const Int panelLDim = nb+A21LocHeight;
        panelBuf.reserve( panelLDim*nb );
        A11_STAR_STAR.Attach
        ( nb, nb, g, 0, 0, &panelBuf[0], panelLDim, 0 );
        A21_MC_STAR.Attach
        ( A21Height, nb, g, A21.ColAlign(), 0, &panelBuf[nb], panelLDim, 0 );
        A11_STAR_STAR = A11;
        A21_MC_STAR = A21;
        lu::Panel( A11_STAR_STAR, A21_MC_STAR, p1Piv, pivotBuf );

        ApplyRowPivots( AB, p1Piv );
        Shift( p1Piv, k );

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
  ElementalMatrix<Int>& rowPiv,
  ElementalMatrix<Int>& colPiv )
{
    DEBUG_ONLY(CSE cse("LU"))
    lu::Full( A, rowPiv, colPiv );
}

#define PROTO(F) \
  template void LU( Matrix<F>& A ); \
  template void LU( ElementalMatrix<F>& A ); \
  template void LU( DistMatrix<F,STAR,STAR>& A ); \
  template void LU \
  ( Matrix<F>& A, \
    Matrix<Int>& rowPiv ); \
  template void LU \
  ( ElementalMatrix<F>& A, \
    ElementalMatrix<Int>& rowPiv ); \
  template void LU \
  ( Matrix<F>& A, \
    Matrix<Int>& rowPiv, \
    Matrix<Int>& colPiv ); \
  template void LU \
  ( ElementalMatrix<F>& A, \
    ElementalMatrix<Int>& rowPiv, \
    ElementalMatrix<Int>& colPiv ); \
  template void LUMod \
  (       Matrix<F>& A, \
          Matrix<Int>& perm, \
    const Matrix<F>& u, \
    const Matrix<F>& v, \
    bool conjugate, Base<F> tau ); \
  template void LUMod \
  (       ElementalMatrix<F>& A, \
          ElementalMatrix<Int>& perm, \
    const ElementalMatrix<F>& u, \
    const ElementalMatrix<F>& v, \
    bool conjugate, Base<F> tau ); \
  template void lu::Panel \
  ( Matrix<F>& APan, \
    Matrix<Int>& p1 ); \
  template void lu::Panel \
  ( DistMatrix<F,  STAR,STAR>& A11, \
    DistMatrix<F,  MC,  STAR>& A21, \
    DistMatrix<Int,STAR,STAR>& p1, \
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
    const Matrix<Int>& rowPiv, \
          Matrix<F>& B ); \
  template void lu::SolveAfter \
  ( Orientation orientation, \
    const ElementalMatrix<F>& A, \
    const ElementalMatrix<Int>& rowPiv, \
          ElementalMatrix<F>& B ); \
  template void lu::SolveAfter \
  ( Orientation orientation, \
    const Matrix<F>& A, \
    const Matrix<Int>& rowPiv, \
    const Matrix<Int>& colPiv, \
          Matrix<F>& B ); \
  template void lu::SolveAfter \
  ( Orientation orientation, \
    const ElementalMatrix<F>& A, \
    const ElementalMatrix<Int>& rowPiv, \
    const ElementalMatrix<Int>& colPiv, \
          ElementalMatrix<F>& B );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
