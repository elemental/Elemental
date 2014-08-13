/*
   Copyright (c) 2009-2014, Jack Poulson
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

        const IndexRange ind0Vert( 0,  ki    ),
                         ind1(     k,  k+nb  ),
                         ind1Vert( ki, ki+nb ),
                         indL( 0, kj+nb );

        auto A0L = View( A, ind0Vert, indL );
        auto A1L = View( A, ind1Vert, indL );
        auto t1 = View( t, ind1, IndexRange(0,1) );
        auto d1 = View( d, ind1, IndexRange(0,1) );

        PanelHouseholder( A1L, t1, d1 );
        ApplyQ( RIGHT, ADJOINT, A1L, t1, d1, A0L );
    }
}

template<typename F> 
inline void
Householder( Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("rq::Householder"))
    Matrix<F> t;
    Matrix<Base<F>> d;
    Householder( A, t, d );
    MakeTrapezoidal( UPPER, A, Min(A.Height(),A.Width()) );
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
    const Int offset = n-m;
    const Int iOff = m-minDim;
    const Int jOff = n-minDim;

    const Grid& g = APre.Grid();
    DistMatrix<F> A(g);
    Copy( APre, A, READ_WRITE_PROXY );

    tPre.Resize( minDim, 1 );
    dPre.Resize( minDim, 1 );
    DistMatrix<F,      MD,STAR> t(g);
    DistMatrix<Base<F>,MD,STAR> d(g);
    t.SetRoot( A.DiagonalRoot(offset) );
    d.SetRoot( A.DiagonalRoot(offset) );
    t.AlignCols( A.DiagonalAlign(offset) );
    d.AlignCols( A.DiagonalAlign(offset) );
    Copy( tPre, t, WRITE_PROXY );
    Copy( dPre, d, WRITE_PROXY );

    const Int bsize = Blocksize();
    const Int kLast = LastOffset( minDim, bsize );
    for( Int k=kLast; k>=0; k-=bsize )
    {
        const Int nb = Min(bsize,minDim-k);

        const Int ki = k + iOff;
        const Int kj = k + jOff;

        const IndexRange ind0Vert( 0,  ki    ),
                         ind1(     k,  k+nb  ),
                         ind1Vert( ki, ki+nb ),
                         indL( 0, kj+nb );

        auto A0L = View( A, ind0Vert, indL );
        auto A1L = View( A, ind1Vert, indL );
        auto t1 = View( t, ind1, IndexRange(0,1) );
        auto d1 = View( d, ind1, IndexRange(0,1) );

        PanelHouseholder( A1L, t1, d1 );
        ApplyQ( RIGHT, ADJOINT, A1L, t1, d1, A0L );
    }
    Copy( A, APre, RESTORE_READ_WRITE_PROXY );
    Copy( t, tPre, RESTORE_WRITE_PROXY      );
    Copy( d, dPre, RESTORE_WRITE_PROXY      );
}

template<typename F> 
inline void
Householder( AbstractDistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("rq::Householder"))
    DistMatrix<F,MD,STAR> t(A.Grid());
    DistMatrix<Base<F>,MD,STAR> d(A.Grid());
    Householder( A, t, d );
    MakeTrapezoidal( UPPER, A, Min(A.Height(),A.Width()) );
}

} // namespace rq
} // namespace El

#endif // ifndef EL_RQ_HOUSEHOLDER_HPP
