/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_HERMITIANTRIDIAG_U_HPP
#define EL_HERMITIANTRIDIAG_U_HPP

#include "./UPan.hpp"

namespace El {
namespace herm_tridiag {

template<typename F>
void U( Matrix<F>& A, Matrix<F>& t )
{
    DEBUG_ONLY(
      CSE cse("herm_tridiag::U");
      if( A.Height() != A.Width() )
          LogicError("A must be square");
    )
    const Int n = A.Height();
    if( n == 0 )
    {
        t.Resize( 0, 1 );
        return;
    }
    t.Resize( n-1, 1 );

    Matrix<F> w01;
    for( Int k=n-1; k>0; --k )
    {
        const Range<Int> ind0( 0, k   ),
                         ind1( k, k+1 );

        auto A00      = A( ind0, ind0 );
        auto a01      = A( ind0, ind1 );

        auto a01T     = A( IR(0,k-1), ind1 ); 
        auto alpha01B = A( IR(k-1),   ind1 );

        const F tau = LeftReflector( alpha01B, a01T );
        const Base<F> epsilon1 = alpha01B.GetRealPart(0,0);
        t.Set( k-1, 0, tau );
        alpha01B.Set(0,0,F(1));

        Zeros( w01, k, 1 );
        Hemv( UPPER, Conj(tau), A00, a01, F(0), w01 );
        const F alpha = -Conj(tau)*Dot( w01, a01 )/F(2);
        Axpy( alpha, a01, w01 );
        Her2( UPPER, F(-1), a01, w01, A00 );
        alpha01B.Set(0,0,epsilon1);
    }
}

// TODO: If there is only a single MPI process, fall down to the sequential
//       implementation.
template<typename F>
void U
( ElementalMatrix<F>& APre,
  ElementalMatrix<F>& tPre,
  const SymvCtrl<F>& ctrl )
{
    DEBUG_ONLY(
      CSE cse("herm_tridiag::U");
      AssertSameGrids( APre, tPre );
      if( APre.Height() != APre.Width() )
          LogicError("A must be square");
    )

    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre );
    DistMatrixWriteProxy<F,F,STAR,STAR> tProx( tPre );
    auto& A = AProx.Get();
    auto& t = tProx.Get();

    const Grid& g = A.Grid();
    const Int n = A.Height();
    if( n == 0 )
    {
        t.Resize( 0, 1 );
        return;
    }
    DistMatrix<F,MD,STAR> tDiag(g);
    tDiag.SetRoot( A.DiagonalRoot(1) );
    tDiag.AlignCols( A.DiagonalAlign(1) );
    tDiag.Resize( n-1, 1 );

    DistMatrix<F> WPan(g);
    DistMatrix<F,STAR,STAR> A11_STAR_STAR(g), t1_STAR_STAR(g);
    DistMatrix<F,MC,  STAR> APan_MC_STAR(g), WPan_MC_STAR(g);
    DistMatrix<F,MR,  STAR> APan_MR_STAR(g), WPan_MR_STAR(g);
    
    const Int bsize = Blocksize();
    const Int kLast = LastOffset( n, bsize );
    for( Int k=kLast; k>=0; k-=bsize )
    {
        const Int nb = Min(bsize,n-k);

        const Range<Int> ind0( 0, k    ),
                         indT( 0, k+nb ), indL( 0, k+nb ),
                         ind1( k, k+nb );

        auto A00 = A( ind0, ind0 );
        auto A01 = A( ind0, ind1 );
        auto A11 = A( ind1, ind1 );
        auto ATL = A( indT, indL );
        
        if( k > 0 )
        {
            auto t1 = tDiag( IR(k-1,k+nb-1), ALL );
            WPan.AlignWith( A01 );
            WPan.Resize( k+nb, nb );
            APan_MC_STAR.AlignWith( A00 );
            APan_MC_STAR.Resize( k+nb, nb );
            WPan_MC_STAR.AlignWith( A00 );
            WPan_MC_STAR.Resize( k+nb, nb );
            APan_MR_STAR.AlignWith( A00 );
            APan_MR_STAR.Resize( k+nb, nb );
            WPan_MR_STAR.AlignWith( A00 );
            WPan_MR_STAR.Resize( k+nb, nb );

            UPan
            ( ATL, WPan, t1,
              APan_MC_STAR, APan_MR_STAR, 
              WPan_MC_STAR, WPan_MR_STAR, ctrl );

            auto A01_MC_STAR = APan_MC_STAR( ind0, ind1-k );
            auto A01_MR_STAR = APan_MR_STAR( ind0, ind1-k );
            auto W01_MC_STAR = WPan_MC_STAR( ind0, ind1-k );
            auto W01_MR_STAR = WPan_MR_STAR( ind0, ind1-k );

            LocalTrr2k
            ( UPPER, NORMAL, ADJOINT, NORMAL, ADJOINT,
              F(-1), A01_MC_STAR, W01_MR_STAR,
              F(-1), W01_MC_STAR, A01_MR_STAR,
              F(1),  A00 );
        }
        else
        {
            auto t1 = tDiag( IR(0,nb-1), ALL );
            A11_STAR_STAR = A11;
            t1_STAR_STAR.Resize( nb-1, 1 );
            HermitianTridiag
            ( UPPER, A11_STAR_STAR.Matrix(), t1_STAR_STAR.Matrix() );
            A11 = A11_STAR_STAR;
            t1 = t1_STAR_STAR;
        }
    }
    // Redistribute from matrix-diagonal form to fully replicated
    t = tDiag;
}

} // namespace herm_tridiag
} // namespace El

#endif // ifndef EL_HERMITIANTRIDIAG_U_HPP
