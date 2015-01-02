/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_LQ_HOUSEHOLDER_HPP
#define EL_LQ_HOUSEHOLDER_HPP

#include "./ApplyQ.hpp"
#include "./PanelHouseholder.hpp"

namespace El {
namespace lq {

template<typename F> 
inline void
Householder( Matrix<F>& A, Matrix<F>& t, Matrix<Base<F>>& d )
{
    DEBUG_ONLY(CallStackEntry cse("lq::Householder"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);
    t.Resize( minDim, 1 );
    d.Resize( minDim, 1 );

    const Int bsize = Blocksize();
    for( Int k=0; k<minDim; k+=bsize )
    {
        const Int nb = Min(bsize,minDim-k);

        const Range<Int> ind1( k, k+nb ), 
                         indR( k, n    ),
                         ind2Vert( k+nb, m );

        auto A1R = A( ind1,     indR    );
        auto A2R = A( ind2Vert, indR    );
        auto t1  = t( ind1,     IR(0,1) );
        auto d1  = d( ind1,     IR(0,1) );

        PanelHouseholder( A1R, t1, d1 );
        ApplyQ( RIGHT, ADJOINT, A1R, t1, d1, A2R );
    }
}

template<typename F> 
inline void
Householder
( AbstractDistMatrix<F>& APre, AbstractDistMatrix<F>& tPre, 
  AbstractDistMatrix<Base<F>>& dPre )
{
    DEBUG_ONLY(
        CallStackEntry cse("Householder");
        AssertSameGrids( APre, tPre, dPre );
    )
    const Int m = APre.Height();
    const Int n = APre.Width();
    const Int minDim = Min(m,n);

    auto APtr = ReadWriteProxy<F,MC,MR>( &APre );
    auto& A = *APtr;

    auto tPtr = WriteProxy<F,      MD,STAR>( &tPre ); auto& t = *tPtr;
    auto dPtr = WriteProxy<Base<F>,MD,STAR>( &dPre ); auto& d = *dPtr;
    t.Resize( minDim, 1 );
    d.Resize( minDim, 1 );

    const Int bsize = Blocksize();
    for( Int k=0; k<minDim; k+=bsize )
    {
        const Int nb = Min(bsize,minDim-k);

        const Range<Int> ind1( k, k+nb ),
                         indR( k, n    ),
                         ind2Vert( k+nb, m );

        auto A1R = A( ind1,     indR    );
        auto A2R = A( ind2Vert, indR    );
        auto t1  = t( ind1,     IR(0,1) );
        auto d1  = d( ind1,     IR(0,1) );

        PanelHouseholder( A1R, t1, d1 );
        ApplyQ( RIGHT, ADJOINT, A1R, t1, d1, A2R );
    }
}

} // namespace lq
} // namespace El

#endif // ifndef EL_LQ_HOUSEHOLDER_HPP
