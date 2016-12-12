/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
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
// It is worth noting that multiple research groups have been actively working
// on nearly identical algorithms, with the relevant timeline being:
//
// 1) May, 2014;
//    Jed Duersch's final presentation for Math 273 at UC Berkeley
//    introduced the fundamental algorithm and its ideas.
//
// 2) May, 2015;
//    P.G. Martinsson, "Blocked rank-revealing QR factorizations:
//    How randomized sampling can be used to avoid single-vector pivoting",
//    arXiv:1505.08115. Available from http://arxiv.org/abs/1505.08115
//
// 3) September, 2015;
//    Jed Duersch and Ming Gu, "True BLAS-3 Performance QRCP using Random
//    Sampling", arXiv:1509.06820.
//    Available from http://arxiv.org/abs/1509.06820
//
// 4) December, 2015:
//    P.G. Martinsson, Gregorio Quintana-Orti, Nathan Heavner, and
//    Robert van de Geijn,"Householder QR Factorization: Adding Randomization
//    for Column-Pivoting", arXiv:1512.02671
//    Available from http://arxiv.org/abs/1512.02671
//

template<typename Field>
class StandardProxy
{
private:
    Int numPower_, numOversample_;

public:
    StandardProxy( Int numPower=1, Int numOversample=10 )
    : numPower_(numPower), numOversample_(numOversample)
    { }

    void operator()
    ( const Matrix<Field>& A,
            Permutation& Omega,
            Int numPivots,
            bool smallestFirst=false ) const
    {
        const Int m = A.Height();

        // Generate a Gaussian random matrix
        Matrix<Field> G;
        Gaussian( G, numPivots+numOversample_, m );
        // TODO(poulson): Force the row norms to be one?
        /*
        Matrix<Base<Field>> rowNorms;
        RowTwoNorms( G, rowNorms );
        DiagonalSolve( LEFT, NORMAL, rowNorms, G );
        */

        // Form  G A (A^H A)^q
        Matrix<Field> Y, Z;
        Gemm( NORMAL, NORMAL, Field(1), G, A, Y );
        for( Int powerIter=0; powerIter<numPower_; ++powerIter )
        {
            Gemm( NORMAL, ADJOINT, Field(1), Y, A, Z );
            Gemm( NORMAL, NORMAL, Field(1), Z, A, Y );
        }

        QRCtrl<Base<Field>> ctrl;
        ctrl.boundRank = true;
        ctrl.maxRank = numPivots;
        ctrl.smallestFirst = smallestFirst;
        Matrix<Field> phase;
        Matrix<Base<Field>> signature;
        QR( Y, phase, signature, Omega, ctrl );
    }

    void operator()
    ( const AbstractDistMatrix<Field>& APre,
            DistPermutation& Omega,
            Int numPivots,
            bool smallestFirst=false ) const
    {
        const Int m = APre.Height();
        const Grid& g = APre.Grid();

        DistMatrixReadProxy<Field,Field,MC,MR> AProxy( APre );
        auto& A = AProxy.GetLocked();

        // Generate a Gaussian random matrix
        DistMatrix<Field> G(g);
        Gaussian( G, numPivots+numOversample_, m );
        // TODO(poulson): Force the row norms to be one?
        /*
        DistMatrix<Base<Field>,MC,STAR> rowNorms(g);
        RowTwoNorms( G, rowNorms );
        DiagonalSolve( LEFT, NORMAL, rowNorms, G );
        */

        // Form  G A (A^H A)^q
        DistMatrix<Field> Y(g), Z(g);
        Gemm( NORMAL, NORMAL, Field(1), G, A, Y );
        for( Int powerIter=0; powerIter<numPower_; ++powerIter )
        {
            Gemm( NORMAL, ADJOINT, Field(1), Y, A, Z );
            Gemm( NORMAL, NORMAL, Field(1), Z, A, Y );
        }

        QRCtrl<Base<Field>> ctrl;
        ctrl.boundRank = true;
        ctrl.maxRank = numPivots;
        ctrl.smallestFirst = smallestFirst;
        DistMatrix<Field,MD,STAR> phase(g);
        DistMatrix<Base<Field>,MD,STAR> signature(g);
        QR( Y, phase, signature, Omega, ctrl );
    }
};

template<typename Field,class ProxyType>
void
ProxyHouseholder
( Matrix<Field>& A,
  Matrix<Field>& phase,
  Matrix<Base<Field>>& signature,
  Permutation& Omega,
  const ProxyType& proxy,
  bool usePanelPerm=false,
  bool smallestFirst=false )
{
    EL_DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);
    phase.Resize( minDim, 1 );
    signature.Resize( minDim, 1 );

    Omega.MakeIdentity( n );
    if( usePanelPerm )
        Omega.ReserveSwaps( 2*n );
    else
        Omega.ReserveSwaps( n );

    Permutation proxyPerm, panelPerm;

    QRCtrl<Base<Field>> ctrl;
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
        proxyPerm.PermuteCols( AR );
        Omega.SwapSequence( proxyPerm, k );

        auto AB1 = A( indB, ind1 );
        auto AB2 = A( indB, ind2 );
        auto phase1 = phase( ind1, ALL );
        auto sig1 = signature( ind1, ALL );

        if( usePanelPerm )
        {
            QR( AB1, phase1, sig1, panelPerm, ctrl );
            auto AT1 = A( IR(0,k), ind1 );
            panelPerm.PermuteCols( AT1 );
            Omega.SwapSequence( panelPerm, k );
        }
        else
        {
            QR( AB1, phase1, sig1 );
        }
        ApplyQ( LEFT, ADJOINT, AB1, phase1, sig1, AB2 );
    }
}

template<typename Field,class ProxyType>
void
ProxyHouseholder
( AbstractDistMatrix<Field>& APre,
  AbstractDistMatrix<Field>& phasePre,
  AbstractDistMatrix<Base<Field>>& signaturePre,
  DistPermutation& Omega,
  const ProxyType& proxy,
  bool usePanelPerm=false,
  bool smallestFirst=false )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(AssertSameGrids( APre, phasePre, signaturePre ))
    const Int m = APre.Height();
    const Int n = APre.Width();
    const Int minDim = Min(m,n);
    const Grid& g = APre.Grid();

    DistMatrixReadWriteProxy<Field,Field,MC,MR> AProx( APre );
    DistMatrixWriteProxy<Field,Field,MD,STAR> phaseProx( phasePre );
    DistMatrixWriteProxy<Base<Field>,Base<Field>,MD,STAR>
      signatureProx( signaturePre );
    auto& A = AProx.Get();
    auto& phase = phaseProx.Get();
    auto& signature = signatureProx.Get();

    phase.Resize( minDim, 1 );
    signature.Resize( minDim, 1 );

    Omega.MakeIdentity( n );
    if( usePanelPerm )
        Omega.ReserveSwaps( 2*n );
    else
        Omega.ReserveSwaps( n );

    DistPermutation proxyPerm(g), panelPerm(g);

    QRCtrl<Base<Field>> ctrl;
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
        proxyPerm.PermuteCols( AR );
        Omega.SwapSequence( proxyPerm, k );

        auto AB1 = A( indB, ind1 );
        auto AB2 = A( indB, ind2 );
        auto phase1 = phase( ind1, ALL );
        auto sig1 = signature( ind1, ALL );

        if( usePanelPerm )
        {
            QR( AB1, phase1, sig1, panelPerm, ctrl );
            auto AT1 = A( IR(0,k), ind1 );
            panelPerm.PermuteCols( AT1 );
            Omega.SwapSequence( panelPerm, k );
        }
        else
        {
            QR( AB1, phase1, sig1 );
        }
        ApplyQ( LEFT, ADJOINT, AB1, phase1, sig1, AB2 );
    }
}

} // namespace qr
} // namespace El

#endif // ifndef EL_QR_PROXY_HOUSEHOLDER_HPP
