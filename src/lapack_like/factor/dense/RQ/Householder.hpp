/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_RQ_HOUSEHOLDER_HPP
#define EL_RQ_HOUSEHOLDER_HPP

#include "./ApplyQ.hpp"
#include "./PanelHouseholder.hpp"

namespace El {
namespace rq {

template<typename F> 
inline void
Householder( Matrix<F>& A, Matrix<F>& t, Matrix<Base<F>>& d )
{
    DEBUG_ONLY(CallStackEntry cse("rq::Householder"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);
    const Int iOff = m-minDim;
    const Int jOff = n-minDim;

    t.Resize( minDim, 1 );
    d.Resize( minDim, 1 );

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
        auto t1 = t( ind1, IR(0,1) );
        auto d1 = d( ind1, IR(0,1) );

        PanelHouseholder( A1L, t1, d1 );
        ApplyQ( RIGHT, ADJOINT, A1L, t1, d1, A0L );
    }
}

template<typename F> 
inline void
Householder
( AbstractDistMatrix<F>& APre, AbstractDistMatrix<F>& tPre, 
  AbstractDistMatrix<Base<F>>& dPre )
{
    DEBUG_ONLY(
        CallStackEntry cse("rq::Householder");
        AssertSameGrids( APre, tPre, dPre );
    )
    const Int m = APre.Height();
    const Int n = APre.Width();
    const Int minDim = Min(m,n);
    const Int iOff = m-minDim;
    const Int jOff = n-minDim;

    auto APtr = ReadWriteProxy<F,MC,MR>( &APre );
    auto& A = *APtr;

    auto tPtr = WriteProxy<F,      MD,STAR>( &tPre ); auto& t = *tPtr;
    auto dPtr = WriteProxy<Base<F>,MD,STAR>( &dPre ); auto& d = *dPtr;
    t.Resize( minDim, 1 );
    d.Resize( minDim, 1 );

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
        auto t1 = t( ind1, IR(0,1) );
        auto d1 = d( ind1, IR(0,1) );

        PanelHouseholder( A1L, t1, d1 );
        ApplyQ( RIGHT, ADJOINT, A1L, t1, d1, A0L );
    }
}

} // namespace rq
} // namespace El

#endif // ifndef EL_RQ_HOUSEHOLDER_HPP
