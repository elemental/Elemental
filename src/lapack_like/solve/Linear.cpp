/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {

namespace lu {

template<typename Field>
void Panel( Matrix<Field>& APan, Permutation& P, Permutation& p1, Int offset );

template<typename Field>
void Panel
( DistMatrix<Field,STAR,STAR>& A11,
  DistMatrix<Field,MC,  STAR>& A21,
  DistPermutation& P,
  DistPermutation& PB,
  Int offset,
  vector<Field>& pivotBuf );

} // namespace lu

// Short-circuited form of LU factorization with partial pivoting
template<typename Field>
void RowEchelon( Matrix<Field>& A, Matrix<Field>& B )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( A.Height() != B.Height() )
          LogicError("A and B must be the same height");
    )

    const Int m = A.Height();
    const Int n = A.Width();
    Permutation P;
    P.MakeIdentity( m );
    P.ReserveSwaps( Min(m,n) );

    Permutation PB;

    const Int mA = A.Height();
    const Int nA = A.Width();
    const Int minDimA = Min(mA,nA);
    const Int bsize = Blocksize();
    for( Int k=0; k<minDimA; k+=bsize )
    {
        const Int nb = Min(bsize,minDimA-k);
        const Range<Int> ind1( k, k+nb ), ind2( k+nb, END ),
                         indB( k, mA   );
        auto A11 = A( ind1, ind1 );
        auto A12 = A( ind1, ind2 );
        auto A21 = A( ind2, ind1 );
        auto A22 = A( ind2, ind2 );
        auto AB1 = A( indB, ind1 );
        auto AB2 = A( indB, ind2 );
        auto B1  = B( ind1, ALL );
        auto B2  = B( ind2, ALL );
        auto BB  = B( indB, ALL );

        lu::Panel( AB1, P, PB, k );
        PB.PermuteRows( AB2 );
        PB.PermuteRows( BB );

        Trsm( LEFT, LOWER, NORMAL, UNIT, Field(1), A11, A12 );
        Trsm( LEFT, LOWER, NORMAL, UNIT, Field(1), A11, B1 );

        Gemm( NORMAL, NORMAL, Field(-1), A21, A12, Field(1), A22 );
        Gemm( NORMAL, NORMAL, Field(-1), A21, B1,  Field(1), B2 );
    }
}

// Short-circuited form of LU factorization with partial pivoting
template<typename Field>
void RowEchelon( DistMatrix<Field>& A, DistMatrix<Field>& B )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      AssertSameGrids( A, B );
      if( A.Height() != B.Height() )
          LogicError("A and B must be the same height");
    )
    const Int mA = A.Height();
    const Int nA = A.Width();
    const Int minDimA = Min(mA,nA);
    const Int bsize = Blocksize();
    const Grid& grid = A.Grid();

    DistPermutation P(grid);
    P.MakeIdentity( mA );
    P.ReserveSwaps( Min(mA,nA) );

    DistPermutation PB(grid);

    DistMatrix<Field,STAR,STAR> A11_STAR_STAR(grid);
    DistMatrix<Field,STAR,VR  > A12_STAR_VR(grid), B1_STAR_VR(grid);
    DistMatrix<Field,STAR,MR  > A12_STAR_MR(grid), B1_STAR_MR(grid);
    DistMatrix<Field,MC,  STAR> A21_MC_STAR(grid);

    // In case B's columns are not aligned with A's
    const bool BAligned = ( B.ColShift() == A.ColShift() );
    DistMatrix<Field,MC,STAR> A21_MC_STAR_B(grid);

    vector<Field> panelBuf, pivotBuf;
    for( Int k=0; k<minDimA; k+=bsize )
    {
        const Int nb = Min(bsize,minDimA-k);
        const Range<Int> ind1( k, k+nb ), ind2( k+nb, END ), indB( k, END );
        auto A11 = A( ind1, ind1 );
        auto A12 = A( ind1, ind2 );
        auto A21 = A( ind2, ind1 );
        auto A22 = A( ind2, ind2 );
        auto AB2 = A( indB, ind2 );
        auto B1  = B( ind1, ALL  );
        auto B2  = B( ind2, ALL  );
        auto BB  = B( indB, ALL  );

        const Int A21Height = A21.Height();
        const Int A21LocHeight = A21.LocalHeight();
        const Int panelLDim = nb+A21LocHeight;
        FastResize( panelBuf, panelLDim*nb );
        A11_STAR_STAR.Attach
        ( nb, nb, grid, 0, 0, &panelBuf[0], panelLDim, 0 );
        A21_MC_STAR.Attach
        ( A21Height, nb, grid, A21.ColAlign(), 0,
          &panelBuf[nb], panelLDim, 0 );
        A11_STAR_STAR = A11;
        A21_MC_STAR = A21;
        lu::Panel( A11_STAR_STAR, A21_MC_STAR, P, PB, k, pivotBuf );
        PB.PermuteRows( AB2 );
        PB.PermuteRows( BB );

        A12_STAR_VR.AlignWith( A22 );
        A12_STAR_VR = A12;
        B1_STAR_VR.AlignWith( B1 );
        B1_STAR_VR = B1;
        LocalTrsm
        ( LEFT, LOWER, NORMAL, UNIT, Field(1), A11_STAR_STAR, A12_STAR_VR );
        LocalTrsm
        ( LEFT, LOWER, NORMAL, UNIT, Field(1), A11_STAR_STAR, B1_STAR_VR );

        A12_STAR_MR.AlignWith( A22 );
        A12_STAR_MR = A12_STAR_VR;
        B1_STAR_MR.AlignWith( B1 );
        B1_STAR_MR = B1_STAR_VR;
        LocalGemm
        ( NORMAL, NORMAL, Field(-1), A21_MC_STAR, A12_STAR_MR, Field(1), A22 );
        if( BAligned )
        {
            LocalGemm
            ( NORMAL, NORMAL,
              Field(-1), A21_MC_STAR, B1_STAR_MR, Field(1), B2 );
        }
        else
        {
            A21_MC_STAR_B.AlignWith( B2 );
            A21_MC_STAR_B = A21_MC_STAR;
            LocalGemm
            ( NORMAL, NORMAL,
              Field(-1), A21_MC_STAR_B, B1_STAR_MR, Field(1), B2 );
        }

        A11 = A11_STAR_STAR;
        A12 = A12_STAR_MR;
        B1 = B1_STAR_MR;
    }
}

namespace lin_solve {

template<typename Field>
void Overwrite( Matrix<Field>& A, Matrix<Field>& B )
{
    EL_DEBUG_CSE
    // Perform Gaussian elimination
    RowEchelon( A, B );
    Trsm( LEFT, UPPER, NORMAL, NON_UNIT, Field(1), A, B );
}

template<typename Field>
void Overwrite
( AbstractDistMatrix<Field>& APre, AbstractDistMatrix<Field>& BPre )
{
    EL_DEBUG_CSE

    DistMatrixReadWriteProxy<Field,Field,MC,MR> AProx( APre ), BProx( BPre );
    auto& A = AProx.Get();
    auto& B = BProx.Get();

    const bool useFullLU = true;

    if( useFullLU )
    {
        DistPermutation P(A.Grid());
        LU( A, P );
        lu::SolveAfter( NORMAL, A, P, B );
    }
    else
    {
        RowEchelon( A, B );
        Trsm( LEFT, UPPER, NORMAL, NON_UNIT, Field(1), A, B );
    }
}

} // namespace lin_solve

template<typename Field>
void LinearSolve( const Matrix<Field>& A, Matrix<Field>& B )
{
    EL_DEBUG_CSE
    Matrix<Field> ACopy( A );
    lin_solve::Overwrite( ACopy, B );
}

namespace lin_solve {

template<typename Field,
         typename=EnableIf<IsBlasScalar<Field>>>
void ScaLAPACKHelper
( const DistMatrix<Field,MC,MR,BLOCK>& A,
        DistMatrix<Field,MC,MR,BLOCK>& B )
{
    AssertScaLAPACKSupport();
#ifdef EL_HAVE_SCALAPACK
    const int m = A.Height();
    const int n = B.Width();
    const int mb = A.BlockHeight();
    const int mLocal = A.LocalHeight();

    auto ACopy = A;
    vector<int> ipiv(mLocal+mb);

    auto descA = FillDesc( A );
    auto descB = FillDesc( B );
    scalapack::LinearSolve
    ( m, n,
      ACopy.Buffer(), descA.data(),
      ipiv.data(),
      B.Buffer(), descB.data() );
#endif
}

template<typename Field,
         typename=DisableIf<IsBlasScalar<Field>>,
         typename=void>
void ScaLAPACKHelper
( const DistMatrix<Field,MC,MR,BLOCK>& A,
        DistMatrix<Field,MC,MR,BLOCK>& B )
{
    LogicError("ScaLAPACK does not support ",TypeName<Field>());
}

} // namespace lin_solve

template<typename Field>
void LinearSolve
( const AbstractDistMatrix<Field>& A,
        AbstractDistMatrix<Field>& B,
  bool scalapack )
{
    EL_DEBUG_CSE
    if( scalapack )
    {
#ifdef EL_HAVE_SCALAPACK
        ProxyCtrl proxyCtrl;
        proxyCtrl.colConstrain = true;
        proxyCtrl.rowConstrain = true;
        proxyCtrl.blockHeight = DefaultBlockHeight();
        proxyCtrl.blockWidth = DefaultBlockWidth();
        proxyCtrl.colAlign = 0;
        proxyCtrl.rowAlign = 0;
        proxyCtrl.colCut = 0;
        proxyCtrl.rowCut = 0;
        DistMatrixReadProxy<Field,Field,MC,MR,BLOCK> AProx( A, proxyCtrl );
        DistMatrixReadWriteProxy<Field,Field,MC,MR,BLOCK> BProx( B, proxyCtrl );
        auto& ABlock = AProx.GetLocked();
        auto& BBlock = BProx.Get();
        lin_solve::ScaLAPACKHelper( ABlock, BBlock );
        return;
#else
        if( A.Grid().Rank() == 0 )
            Output
            ("WARNING: Requested a ScaLAPACK solve, "
             "but ScaLAPACK was not available");
#endif
    }
    DistMatrix<Field> ACopy( A );
    lin_solve::Overwrite( ACopy, B );
}

template<typename Field>
void LinearSolve
( const SparseMatrix<Field>& A,
        Matrix<Field>& B,
  const LeastSquaresCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    Matrix<Field> X;
    LeastSquares( NORMAL, A, B, X, ctrl );
    B = X;
}

template<typename Field>
void LinearSolve
( const DistSparseMatrix<Field>& A,
        DistMultiVec<Field>& B,
  const LeastSquaresCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    DistMultiVec<Field> X(B.Grid());
    LeastSquares( NORMAL, A, B, X, ctrl );
    B = X;
}

#define PROTO(Field) \
  template void lin_solve::Overwrite( Matrix<Field>& A, Matrix<Field>& B ); \
  template void lin_solve::Overwrite \
  ( AbstractDistMatrix<Field>& A, AbstractDistMatrix<Field>& B ); \
  template void LinearSolve( const Matrix<Field>& A, Matrix<Field>& B ); \
  template void LinearSolve \
  ( const AbstractDistMatrix<Field>& A, \
          AbstractDistMatrix<Field>& B, \
    bool scalapack ); \
  template void LinearSolve \
  ( const SparseMatrix<Field>& A, \
          Matrix<Field>& B, \
    const LeastSquaresCtrl<Base<Field>>& ctrl ); \
  template void LinearSolve \
  ( const DistSparseMatrix<Field>& A, \
          DistMultiVec<Field>& B, \
    const LeastSquaresCtrl<Base<Field>>& ctrl );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
