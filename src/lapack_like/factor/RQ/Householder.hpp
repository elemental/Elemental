/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_RQ_HOUSEHOLDER_HPP
#define EL_RQ_HOUSEHOLDER_HPP

#include "./ApplyQ.hpp"
#include "./PanelHouseholder.hpp"

namespace El {
namespace rq {

template<typename F> 
void
Householder
( Matrix<F>& A,
  Matrix<F>& phase,
  Matrix<Base<F>>& signature )
{
    DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);
    const Int iOff = m-minDim;
    const Int jOff = n-minDim;

    phase.Resize( minDim, 1 );
    signature.Resize( minDim, 1 );

    const Int bsize = Blocksize();
    const Int kLast = LastOffset( minDim, bsize );
    for( Int k=kLast; k>=0; k-=bsize )
    {
        const Int nb = Min(bsize,minDim-k);

        const Int ki = k + iOff;
        const Int kj = k + jOff;

        const Range<Int> ind0Vert( 0,  ki    ),
                         ind1(     k,  k+nb  ),
                         ind1Vert( ki, ki+nb ),
                         indL( 0, kj+nb );

        auto A0L = A( ind0Vert, indL );
        auto A1L = A( ind1Vert, indL );
        auto phase1 = phase( ind1, ALL );
        auto sig1 = signature( ind1, ALL );

        PanelHouseholder( A1L, phase1, sig1 );
        ApplyQ( RIGHT, ADJOINT, A1L, phase1, sig1, A0L );
    }
}

template<typename F> 
void
Householder
( ElementalMatrix<F>& APre,
  ElementalMatrix<F>& phasePre, 
  ElementalMatrix<Base<F>>& signaturePre )
{
    DEBUG_CSE
    DEBUG_ONLY(AssertSameGrids( APre, phasePre, signaturePre ))
    const Int m = APre.Height();
    const Int n = APre.Width();
    const Int minDim = Min(m,n);
    const Int iOff = m-minDim;
    const Int jOff = n-minDim;

    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre );
    DistMatrixWriteProxy<F,F,MD,STAR> phaseProx( phasePre );
    DistMatrixWriteProxy<Base<F>,Base<F>,MD,STAR> signatureProx( signaturePre );
    auto& A = AProx.Get();
    auto& phase = phaseProx.Get();
    auto& signature = signatureProx.Get();
    phase.Resize( minDim, 1 );
    signature.Resize( minDim, 1 );

    const Int bsize = Blocksize();
    const Int kLast = LastOffset( minDim, bsize );
    for( Int k=kLast; k>=0; k-=bsize )
    {
        const Int nb = Min(bsize,minDim-k);

        const Int ki = k + iOff;
        const Int kj = k + jOff;

        const Range<Int> ind0Vert( 0,  ki    ),
                         ind1(     k,  k+nb  ),
                         ind1Vert( ki, ki+nb ),
                         indL( 0, kj+nb );

        auto A0L = A( ind0Vert, indL );
        auto A1L = A( ind1Vert, indL );
        auto phase1 = phase( ind1, ALL );
        auto sig1 = signature( ind1, ALL );

        PanelHouseholder( A1L, phase1, sig1 );
        ApplyQ( RIGHT, ADJOINT, A1L, phase1, sig1, A0L );
    }
}

} // namespace rq
} // namespace El

#endif // ifndef EL_RQ_HOUSEHOLDER_HPP
