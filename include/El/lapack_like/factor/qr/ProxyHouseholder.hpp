/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_QR_PROXY_HOUSEHOLDER_HPP
#define EL_QR_PROXY_HOUSEHOLDER_HPP

namespace El {
namespace qr {

// The following implementation is based upon the algorithm in Fig. 3 from
// 
//   P.G. Martinsson,
//   "Blocked rank-revealing QR factorizations: How randomized sampling can
//    be used to avoid single-vector pivoting",
//   Available from http://arxiv.org/abs/1505.08115, May 2015.
//
// While I am currently unhappy with the forfeiting of optimal pivots, the
// performance improvements should be substantial once Elemental is moved
// back from arbitrary permutations back to pivot sequences.

template<typename F>
class StandardProxy
{
private:
    Int numPower_, numOversample_;

public:
    StandardProxy( Int numPower=2, Int numOversample=5 )
    : numPower_(numPower), numOversample_(numOversample)
    { }

    void operator()
    ( const Matrix<F>& A,
            Matrix<Int>& perm,
            Int numPivots,
            bool smallestFirst=false ) const
    {
        const Int m = A.Height();

        // Generate a Gaussian random matrix 
        Matrix<F> Omega; 
        Gaussian( Omega, numPivots+numOversample_, m );
        // TODO: Force the row norms to be one?
        /*
        Matrix<Base<F>> rowNorms;
        RowTwoNorms( Omega, rowNorms );
        DiagonalSolve( LEFT, NORMAL, rowNorms, Omega );
        */

        // Form  Omega A (A^H A)^q 
        Matrix<F> Y, Z; 
        Gemm( NORMAL, NORMAL, F(1), Omega, A, Y );
        for( Int powerIter=0; powerIter<numPower_; ++powerIter )
        {
            Gemm( NORMAL, ADJOINT, F(1), Y, A, Z );
            Gemm( NORMAL, NORMAL, F(1), Z, A, Y );
        }

        QRCtrl<Base<F>> ctrl;
        ctrl.boundRank = true;
        ctrl.maxRank = numPivots;
        ctrl.smallestFirst = smallestFirst;
        Matrix<F> t, d;
        QR( Y, t, d, perm, ctrl );
    }

    void operator()
    ( const ElementalMatrix<F>& APre,
            ElementalMatrix<Int>& perm,
            Int numPivots,
            bool smallestFirst=false ) const
    {
        const Int m = APre.Height();
        const Grid& g = APre.Grid();

        DistMatrixReadProxy<F,F,MC,MR> AProxy( APre ); 
        auto& A = AProxy.GetLocked();

        // Generate a Gaussian random matrix
        DistMatrix<F> Omega(g);
        Gaussian( Omega, numPivots+numOversample_, m );
        // TODO: Force the row norms to be one?
        /*
        DistMatrix<Base<F>,MC,STAR> rowNorms(g);
        RowTwoNorms( Omega, rowNorms );
        DiagonalSolve( LEFT, NORMAL, rowNorms, Omega );
        */

        // Form  Omega A (A^H A)^q 
        DistMatrix<F> Y(g), Z(g);
        Gemm( NORMAL, NORMAL, F(1), Omega, A, Y );
        for( Int powerIter=0; powerIter<numPower_; ++powerIter )
        {
            Gemm( NORMAL, ADJOINT, F(1), Y, A, Z );
            Gemm( NORMAL, NORMAL, F(1), Z, A, Y );
        }

        QRCtrl<Base<F>> ctrl;
        ctrl.boundRank = true;
        ctrl.maxRank = numPivots;
        ctrl.smallestFirst = smallestFirst;
        DistMatrix<F,MD,STAR> t(g), d(g);
        QR( Y, t, d, perm, ctrl );
    }
};

template<typename F,class ProxyType> 
inline void
ProxyHouseholder
( Matrix<F>& A,
  Matrix<F>& t,
  Matrix<Base<F>>& d,
  Matrix<Int>& p,
  const ProxyType& proxy,
  bool usePanelPerm=true,
  bool smallestFirst=false )
{
    DEBUG_ONLY(CSE cse("qr::ProxyHouseholder"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);
    t.Resize( minDim, 1 );
    d.Resize( minDim, 1 );

    Matrix<Int> pInv;
    pInv.Resize( n, 1 );
    for( Int j=0; j<n; ++j )
        pInv.Set( j, 0, j );

    Matrix<Int> proxyPerm, panelPerm;

    QRCtrl<Base<F>> ctrl;
    ctrl.smallestFirst = smallestFirst;

    const Int bsize = Blocksize();
    for( Int k=0; k<minDim; k+=bsize )
    {
        const Int nb = Min(bsize,minDim-k);

        const Range<Int> ind1( k,    k+nb ),
                         indB( k,    END  ),
                         indR( k,    END  ),
                         ind2( k+nb, END  );

        // Decide this set of pivots using a proxy
        auto ABR = A(indB,indR);
        proxy( ABR, proxyPerm, nb, smallestFirst ); 

        auto AR = A(ALL,indR);
        PermuteCols( AR, proxyPerm );

        auto pInvB = pInv(indB,ALL);
        PermuteRows( pInvB, proxyPerm );

        auto AB1 = A( indB, ind1 );
        auto AB2 = A( indB, ind2 );
        auto t1 = t( ind1, ALL );
        auto d1 = d( ind1, ALL );

        if( usePanelPerm )
        {
            QR( AB1, t1, d1, panelPerm, ctrl );
            auto pInv1 = pInv(ind1,ALL);
            PermuteRows( pInv1, panelPerm );
        }
        else
        {
            QR( AB1, t1, d1 );
        }
        ApplyQ( LEFT, ADJOINT, AB1, t1, d1, AB2 );
    }
    InvertPermutation( pInv, p );
}

template<typename F,class ProxyType> 
inline void
ProxyHouseholder
( ElementalMatrix<F>& APre,
  ElementalMatrix<F>& tPre, 
  ElementalMatrix<Base<F>>& dPre,
  ElementalMatrix<Int>& p,
  const ProxyType& proxy,
  bool usePanelPerm=true,
  bool smallestFirst=false )
{
    DEBUG_ONLY(
      CSE cse("qr::ProxyHouseholder");
      AssertSameGrids( APre, tPre, dPre, p );
    )
    const Int m = APre.Height();
    const Int n = APre.Width();
    const Int minDim = Min(m,n);
    const Grid& g = APre.Grid();

    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre );
    DistMatrixWriteProxy<F,F,MD,STAR> tProx( tPre );
    DistMatrixWriteProxy<Base<F>,Base<F>,MD,STAR> dProx( dPre );
    auto& A = AProx.Get();
    auto& t = tProx.Get();
    auto& d = dProx.Get();

    t.Resize( minDim, 1 );
    d.Resize( minDim, 1 );

    DistMatrix<Int,VC,STAR> pInv(g);
    pInv.Resize( n, 1 );
    for( Int jLoc=0; jLoc<pInv.LocalHeight(); ++jLoc )
        pInv.SetLocal( jLoc, 0, pInv.GlobalRow(jLoc) );

    DistMatrix<Int,VC,STAR> proxyPerm(g), panelPerm(g);

    QRCtrl<Base<F>> ctrl;
    ctrl.smallestFirst = smallestFirst;

    const Int bsize = Blocksize();
    for( Int k=0; k<minDim; k+=bsize )
    {
        const Int nb = Min(bsize,minDim-k);

        const Range<Int> ind1( k,    k+nb ),
                         indB( k,    END  ),
                         indR( k,    END  ),
                         ind2( k+nb, END  );

        // Decide this set of pivots using a proxy
        auto ABR = A(indB,indR);
        proxy( ABR, proxyPerm, nb, smallestFirst ); 

        auto AR = A(ALL,indR);
        PermuteCols( AR, proxyPerm );

        auto pInvB = pInv(indB,ALL);
        PermuteRows( pInvB, proxyPerm );

        auto AB1 = A( indB, ind1 );
        auto AB2 = A( indB, ind2 );
        auto t1 = t( ind1, ALL );
        auto d1 = d( ind1, ALL );

        if( usePanelPerm )
        {
            QR( AB1, t1, d1, panelPerm, ctrl );
            auto pInv1 = pInv(ind1,ALL);
            PermuteRows( pInv1, panelPerm );
        }
        else
        {
            QR( AB1, t1, d1 );
        }
        ApplyQ( LEFT, ADJOINT, AB1, t1, d1, AB2 );
    }
    InvertPermutation( pInv, p );
}

} // namespace qr
} // namespace El

#endif // ifndef EL_QR_PROXY_HOUSEHOLDER_HPP
