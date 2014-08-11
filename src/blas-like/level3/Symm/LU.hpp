/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace El {
namespace symm {

template<typename T>
void LocalAccumulateLU
( Orientation orientation, T alpha,
  const DistMatrix<T>& A,
  const DistMatrix<T,MC,  STAR>& B_MC_STAR,
  const DistMatrix<T,STAR,MR  >& BTrans_STAR_MR,
        DistMatrix<T,MC,  STAR>& Z_MC_STAR,
        DistMatrix<T,MR,  STAR>& Z_MR_STAR )
{
    DEBUG_ONLY(
        CallStackEntry cse("symm::LocalAccumulateLU");
        AssertSameGrids( A, B_MC_STAR, BTrans_STAR_MR, Z_MC_STAR, Z_MR_STAR );
        if( A.Height() != A.Width() ||
            A.Height() != B_MC_STAR.Height() ||
            A.Height() != BTrans_STAR_MR.Width() ||
            A.Height() != Z_MC_STAR.Height() ||
            A.Height() != Z_MR_STAR.Height() ||
            B_MC_STAR.Width() != BTrans_STAR_MR.Height() ||
            BTrans_STAR_MR.Height() != Z_MC_STAR.Width() ||
            Z_MC_STAR.Width() != Z_MR_STAR.Width() )
            LogicError
            ("Nonconformal:\n",
             DimsString(A,"A"),"\n",
             DimsString(B_MC_STAR,"B[MC,* ]"),"\n",
             DimsString(BTrans_STAR_MR,"B'[* ,MR]"),"\n",
             DimsString(Z_MC_STAR,"Z[MC,* ]"),"\n",
             DimsString(Z_MR_STAR,"Z[MR,* ]"));
        if( B_MC_STAR.ColAlign() != A.ColAlign() ||
            BTrans_STAR_MR.RowAlign() != A.RowAlign() ||
            Z_MC_STAR.ColAlign() != A.ColAlign() ||
            Z_MR_STAR.ColAlign() != A.RowAlign() )
            LogicError("Partial matrix distributions are misaligned");
    )
    const Int m = B_MC_STAR.Height();
    const Int n = B_MC_STAR.Width();
    const Grid& g = A.Grid();
    const Int ratio = Max( g.Height(), g.Width() );
    const Int bsize = ratio*Blocksize();

    DistMatrix<T> D11(g);

    for( Int k=0; k<m; k+=bsize )
    {
        const Int nb = Min(bsize,m-k);

        auto A11 = LockedViewRange( A, k, k,    k+nb, k+nb );
        auto A12 = LockedViewRange( A, k, k+nb, k+nb, m    );

        auto B1_MC_STAR = LockedViewRange( B_MC_STAR, k, 0, k+nb, n );

        auto B1Trans_STAR_MR = 
            LockedViewRange( BTrans_STAR_MR, 0, k,    n, k+nb );
        auto B2Trans_STAR_MR = 
            LockedViewRange( BTrans_STAR_MR, 0, k+nb, n, m    );

        auto Z1_MC_STAR = ViewRange( Z_MC_STAR, k,    0, k+nb, n );

        auto Z1_MR_STAR = ViewRange( Z_MR_STAR, k,    0, k+nb, n );
        auto Z2_MR_STAR = ViewRange( Z_MR_STAR, k+nb, 0, m,    n );

        D11.AlignWith( A11 );
        D11 = A11;
        MakeTriangular( UPPER, D11 );
        LocalGemm
        ( NORMAL, orientation, alpha, D11, B1Trans_STAR_MR, T(1), Z1_MC_STAR );
        SetDiagonal( D11, T(0) );

        LocalGemm
        ( orientation, NORMAL, alpha, D11, B1_MC_STAR, T(1), Z1_MR_STAR );

        LocalGemm
        ( NORMAL, orientation, alpha, A12, B2Trans_STAR_MR, T(1), Z1_MC_STAR );

        LocalGemm
        ( orientation, NORMAL, alpha, A12, B1_MC_STAR, T(1), Z2_MR_STAR );
    }
}

template<typename T>
inline void
LUA
( T alpha, const AbstractDistMatrix<T>& APre, const AbstractDistMatrix<T>& BPre,
  T beta,        AbstractDistMatrix<T>& CPre,
  bool conjugate=false )
{
    DEBUG_ONLY(
        CallStackEntry cse("symm::LUA");
        AssertSameGrids( APre, BPre, CPre );
    )
    const Int m = CPre.Height();
    const Int n = CPre.Width();
    const Int bsize = Blocksize();
    const Grid& g = APre.Grid();
    const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );

    // Force 'A', 'B', and 'C' to be in [MC,MR] distributions
    DistMatrix<T> A(g), B(g), C(g);
    Copy( APre, A, READ_PROXY );
    Copy( BPre, B, READ_PROXY );
    Copy( CPre, C, READ_WRITE_PROXY );

    DistMatrix<T,MC,STAR> B1_MC_STAR(g);
    DistMatrix<T,VR,STAR> B1_VR_STAR(g);
    DistMatrix<T,STAR,MR> B1Trans_STAR_MR(g);
    DistMatrix<T> Z1(g);
    DistMatrix<T,MC,STAR> Z1_MC_STAR(g);
    DistMatrix<T,MR,STAR> Z1_MR_STAR(g);
    DistMatrix<T,MR,MC  > Z1_MR_MC(g);

    B1_MC_STAR.AlignWith( A );
    B1_VR_STAR.AlignWith( A );
    B1Trans_STAR_MR.AlignWith( A );
    Z1_MC_STAR.AlignWith( A );
    Z1_MR_STAR.AlignWith( A );

    Scale( beta, C );
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);

        auto B1 = LockedView( B, 0, k, m, nb );
        auto C1 = LockedView( C, 0, k, m, nb );

        B1_MC_STAR = B1;
        B1_VR_STAR = B1_MC_STAR;
        B1_VR_STAR.TransposePartialColAllGather( B1Trans_STAR_MR, conjugate );
        Zeros( Z1_MC_STAR, m, nb );
        Zeros( Z1_MR_STAR, m, nb );
        LocalAccumulateLU
        ( orientation,
          alpha, A, B1_MC_STAR, B1Trans_STAR_MR, Z1_MC_STAR, Z1_MR_STAR );

        Z1_MR_MC.RowSumScatterFrom( Z1_MR_STAR );
        Z1.AlignWith( C1 );
        Z1 = Z1_MR_MC;
        Z1.RowSumScatterUpdate( T(1), Z1_MC_STAR );
        Axpy( T(1), Z1, C1 );
    }

    Copy( C, CPre, RESTORE_READ_WRITE_PROXY );
}

template<typename T>
inline void
LUC
( T alpha, const AbstractDistMatrix<T>& APre, const AbstractDistMatrix<T>& BPre,
  T beta,        AbstractDistMatrix<T>& CPre,
  bool conjugate=false )
{
    DEBUG_ONLY(
        CallStackEntry cse("symm::LUC");
        AssertSameGrids( APre, BPre, CPre );
    )
    const Int m = CPre.Height();
    const Int n = CPre.Width();
    const Int bsize = Blocksize();
    const Grid& g = APre.Grid();
    const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );

    // Force 'A', 'B', and 'C' to be in [MC,MR] distributions
    DistMatrix<T> A(g), B(g), C(g);
    Copy( APre, A, READ_PROXY );
    Copy( BPre, B, READ_PROXY );
    Copy( CPre, C, READ_WRITE_PROXY );

    // Temporary distributions
    DistMatrix<T,MC,  STAR> AT1_MC_STAR(g);
    DistMatrix<T,STAR,MC  > A1R_STAR_MC(g);
    DistMatrix<T,MR,  STAR> B1Trans_MR_STAR(g);

    B1Trans_MR_STAR.AlignWith( C );

    Scale( beta, C );
    for( Int k=0; k<m; k+=bsize )
    {
        const Int nb = Min(bsize,m-k);

        auto A1R = LockedViewRange( A, k, k, k+nb, m    );
        auto AT1 = LockedViewRange( A, 0, k, k+nb, k+nb );

        auto B1 = LockedViewRange( B, k, 0, k+nb, n );

        auto CAbove = ViewRange( C, 0, 0, k+nb, n );
        auto CBelow = ViewRange( C, k, 0, m,    n );

        AT1_MC_STAR.AlignWith( CAbove );
        A1R_STAR_MC.AlignWith( CBelow );
        AT1_MC_STAR = AT1;
        A1R_STAR_MC = A1R;
        MakeTrapezoidal( UPPER, AT1_MC_STAR, -k );
        MakeTrapezoidal( UPPER, A1R_STAR_MC, 1 );

        B1.TransposeColAllGather( B1Trans_MR_STAR );

        LocalGemm
        ( NORMAL, TRANSPOSE, 
          alpha, AT1_MC_STAR, B1Trans_MR_STAR, T(1), CAbove );

        LocalGemm
        ( orientation, TRANSPOSE, 
          alpha, A1R_STAR_MC, B1Trans_MR_STAR, T(1), CBelow );
    }

    Copy( C, CPre, RESTORE_READ_WRITE_PROXY );
}

template<typename T>
inline void
LU
( T alpha, const AbstractDistMatrix<T>& A, const AbstractDistMatrix<T>& B,
  T beta,        AbstractDistMatrix<T>& C, bool conjugate=false )
{
    DEBUG_ONLY(CallStackEntry cse("symm::LU"))
    // TODO: Come up with a better routing mechanism
    if( A.Height() > 5*B.Width() )
        symm::LUA( alpha, A, B, beta, C, conjugate );
    else
        symm::LUC( alpha, A, B, beta, C, conjugate );
}

} // namespace symm
} // namespace El
