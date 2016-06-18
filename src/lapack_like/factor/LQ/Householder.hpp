/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_LQ_HOUSEHOLDER_HPP
#define EL_LQ_HOUSEHOLDER_HPP

#include "./ApplyQ.hpp"
#include "./PanelHouseholder.hpp"

namespace El {
namespace lq {

template<typename F> 
void
Householder( Matrix<F>& A, Matrix<F>& phase, Matrix<Base<F>>& signature )
{
    DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);
    phase.Resize( minDim, 1 );
    signature.Resize( minDim, 1 );

    const Int bsize = Blocksize();
    for( Int k=0; k<minDim; k+=bsize )
    {
        const Int nb = Min(bsize,minDim-k);

        const Range<Int> ind1( k,    k+nb ), 
                         indR( k,    END  ),
                         ind2( k+nb, END  );

        auto A1R = A( ind1, indR );
        auto A2R = A( ind2, indR );
        auto phase1 = phase( ind1, ALL  );
        auto sig1 = signature( ind1, ALL  );

        PanelHouseholder( A1R, phase1, sig1 );
        ApplyQ( RIGHT, ADJOINT, A1R, phase1, sig1, A2R );
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

    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre );
    DistMatrixWriteProxy<F,F,MD,STAR> phaseProx( phasePre );
    DistMatrixWriteProxy<Base<F>,Base<F>,MD,STAR> signatureProx( signaturePre );
    auto& A = AProx.Get();
    auto& phase = phaseProx.Get();
    auto& signature = signatureProx.Get();
    phase.Resize( minDim, 1 );
    signature.Resize( minDim, 1 );

    const Int bsize = Blocksize();
    for( Int k=0; k<minDim; k+=bsize )
    {
        const Int nb = Min(bsize,minDim-k);

        const Range<Int> ind1( k, k+nb ),
                         indR( k, END  ),
                         ind2( k+nb, END );

        auto A1R = A( ind1, indR );
        auto A2R = A( ind2, indR );
        auto phase1 = phase( ind1, ALL  );
        auto sig1 = signature( ind1, ALL  );

        PanelHouseholder( A1R, phase1, sig1 );
        ApplyQ( RIGHT, ADJOINT, A1R, phase1, sig1, A2R );
    }
}

} // namespace lq
} // namespace El

#endif // ifndef EL_LQ_HOUSEHOLDER_HPP
