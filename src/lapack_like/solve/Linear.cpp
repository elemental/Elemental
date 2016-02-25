/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

namespace lu {

template<typename F>
void Panel( Matrix<F>& APan, Permutation& P, Permutation& p1, Int offset );

template<typename F>
void Panel
( DistMatrix<F,  STAR,STAR>& A11,
  DistMatrix<F,  MC,  STAR>& A21,
  DistPermutation& P,
  DistPermutation& PB,
  Int offset,
  vector<F>& pivotBuf );

} // namespace lu

// Short-circuited form of LU factorization with partial pivoting
template<typename F> 
inline void
RowEchelon( Matrix<F>& A, Matrix<F>& B )
{
    DEBUG_ONLY(
      CSE cse("RowEchelon");
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
      CSE cse("RowEchelon");
      AssertSameGrids( A, B );
      if( A.Height() != B.Height() )
          LogicError("A and B must be the same height");
    )
    const Int mA = A.Height();
    const Int nA = A.Width();
    const Int minDimA = Min(mA,nA);
    const Int bsize = Blocksize();
    const Grid& g = A.Grid();

    DistPermutation P(g);
    P.MakeIdentity( mA );
    P.ReserveSwaps( Min(mA,nA) );

    DistPermutation PB(g);

    DistMatrix<F,STAR,STAR> A11_STAR_STAR(g);
    DistMatrix<F,STAR,VR  > A12_STAR_VR(g), B1_STAR_VR(g);
    DistMatrix<F,STAR,MR  > A12_STAR_MR(g), B1_STAR_MR(g);
    DistMatrix<F,MC,  STAR> A21_MC_STAR(g);

    // In case B's columns are not aligned with A's
    const bool BAligned = ( B.ColShift() == A.ColShift() );
    DistMatrix<F,MC,STAR> A21_MC_STAR_B(g);

    vector<F> panelBuf, pivotBuf;
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
        ( nb, nb, g, 0, 0, &panelBuf[0], panelLDim, 0 );
        A21_MC_STAR.Attach
        ( A21Height, nb, g, A21.ColAlign(), 0, &panelBuf[nb], panelLDim, 0 );
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

namespace lin_solve {

template<typename F> 
void Overwrite( Matrix<F>& A, Matrix<F>& B )
{
    DEBUG_ONLY(CSE cse("lin_solve::Overwrite"))
    // Perform Gaussian elimination
    RowEchelon( A, B );
    Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), A, B );
}

template<typename F> 
void Overwrite
( ElementalMatrix<F>& APre, ElementalMatrix<F>& BPre )
{
    DEBUG_ONLY(CSE cse("lin_solve::Overwrite"))

    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre ), BProx( BPre );
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
        Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), A, B );
    }
}

} // namespace lin_solve

template<typename F> 
void LinearSolve( const Matrix<F>& A, Matrix<F>& B )
{
    DEBUG_ONLY(CSE cse("LinearSolve"))
    Matrix<F> ACopy( A );
    lin_solve::Overwrite( ACopy, B );
}

template<typename F> 
void LinearSolve
( const ElementalMatrix<F>& A,
        ElementalMatrix<F>& B )
{
    DEBUG_ONLY(CSE cse("LinearSolve"))
    DistMatrix<F> ACopy( A );
    lin_solve::Overwrite( ACopy, B );
}

namespace lin_solve {

template<typename F,typename=EnableIf<IsBlasScalar<F>>>
void ScaLAPACKHelper
( const DistMatrix<F,MC,MR,BLOCK>& A,
        DistMatrix<F,MC,MR,BLOCK>& B )
{
    AssertScaLAPACKSupport();
#ifdef EL_HAVE_SCALAPACK
    const int m = A.Height();
    const int n = B.Width();
    const int mb = A.BlockHeight();
    const int mLocal = A.LocalHeight();

    auto ACopy = A;
    vector<int> ipiv(mLocal+mb);

    const int bHandle = blacs::Handle( A );
    const int context = blacs::GridInit( bHandle, A );
    auto descA = FillDesc( A, context );
    auto descB = FillDesc( B, context );

    scalapack::LinearSolve
    ( m, n,
      ACopy.Buffer(), descA.data(),
      ipiv.data(),
      B.Buffer(), descB.data() );

    // TODO: Cache context, handle, and exit BLACS during El::Finalize()
    blacs::FreeGrid( context );
    blacs::FreeHandle( bHandle );
#endif
}

template<typename F,typename=DisableIf<IsBlasScalar<F>>,typename=void>
void ScaLAPACKHelper
( const DistMatrix<F,MC,MR,BLOCK>& A,
        DistMatrix<F,MC,MR,BLOCK>& B )
{
    LogicError("ScaLAPACK does not support this datatype");
}

} // namespace lin_solve

template<typename F>
void LinearSolve
( const DistMatrix<F,MC,MR,BLOCK>& A,
        DistMatrix<F,MC,MR,BLOCK>& B )
{
    DEBUG_ONLY(CSE cse("LinearSolve"))
    lin_solve::ScaLAPACKHelper( A, B );
}

template<typename F>
void LinearSolve
( const SparseMatrix<F>& A, Matrix<F>& B, 
  const LeastSquaresCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("LinearSolve"))
    Matrix<F> X;
    LeastSquares( NORMAL, A, B, X, ctrl );
    B = X;
}

template<typename F>
void LinearSolve
( const DistSparseMatrix<F>& A, DistMultiVec<F>& B, 
  const LeastSquaresCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("LinearSolve"))
    DistMultiVec<F> X;
    X.SetComm( B.Comm() );
    LeastSquares( NORMAL, A, B, X, ctrl );
    B = X;
}

#define PROTO(F) \
  template void lin_solve::Overwrite( Matrix<F>& A, Matrix<F>& B ); \
  template void lin_solve::Overwrite \
  ( ElementalMatrix<F>& A, ElementalMatrix<F>& B ); \
  template void LinearSolve( const Matrix<F>& A, Matrix<F>& B ); \
  template void LinearSolve \
  ( const ElementalMatrix<F>& A, \
          ElementalMatrix<F>& B ); \
  template void LinearSolve \
  ( const DistMatrix<F,MC,MR,BLOCK>& A, \
          DistMatrix<F,MC,MR,BLOCK>& B ); \
  template void LinearSolve \
  ( const SparseMatrix<F>& A, Matrix<F>& B, \
    const LeastSquaresCtrl<Base<F>>& ctrl ); \
  template void LinearSolve \
  ( const DistSparseMatrix<F>& A, DistMultiVec<F>& B, \
    const LeastSquaresCtrl<Base<F>>& ctrl );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

} // namespace El
