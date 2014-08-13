/*
   Copyright (c) 2009-2014, Jack Poulson
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

        const IndexRange ind1( k, k+nb ), 
                         indR( k, n    ),
                         ind2Vert( k+nb, m );

        auto A1R = View( A, ind1,     indR            );
        auto A2R = View( A, ind2Vert, indR            );
        auto t1  = View( t, ind1,     IndexRange(0,1) );
        auto d1  = View( d, ind1,     IndexRange(0,1) );

        PanelHouseholder( A1R, t1, d1 );
        ApplyQ( RIGHT, ADJOINT, A1R, t1, d1, A2R );
    }
}

template<typename F> 
inline void
Householder( Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("lq::Householder"))
    Matrix<F> t;
    Matrix<Base<F>> d;
    Householder( A, t, d );
    MakeTriangular( LOWER, A );
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
    const Grid& g = APre.Grid();

    DistMatrix<F> A(g);
    Copy( APre, A, READ_WRITE_PROXY );

    tPre.Resize( minDim, 1 );
    dPre.Resize( minDim, 1 );
    DistMatrix<F,      MD,STAR> t(g);
    DistMatrix<Base<F>,MD,STAR> d(g);
    t.SetRoot( A.DiagonalRoot() );
    d.SetRoot( A.DiagonalRoot() );
    t.AlignCols( A.DiagonalAlign() );
    d.AlignCols( A.DiagonalAlign() );
    Copy( tPre, t, WRITE_PROXY );
    Copy( dPre, d, WRITE_PROXY );

    const Int bsize = Blocksize();
    for( Int k=0; k<minDim; k+=bsize )
    {
        const Int nb = Min(bsize,minDim-k);

        const IndexRange ind1( k, k+nb ),
                         indR( k, n    ),
                         ind2Vert( k+nb, m );

        auto A1R = View( A, ind1,     indR            );
        auto A2R = View( A, ind2Vert, indR            );
        auto t1  = View( t, ind1,     IndexRange(0,1) );
        auto d1  = View( d, ind1,     IndexRange(0,1) );

        PanelHouseholder( A1R, t1, d1 );
        ApplyQ( RIGHT, ADJOINT, A1R, t1, d1, A2R );
    }
    Copy( A, APre, RESTORE_READ_WRITE_PROXY );
    Copy( t, tPre, RESTORE_WRITE_PROXY      );
    Copy( d, dPre, RESTORE_WRITE_PROXY      );
}

template<typename F> 
inline void
Householder( AbstractDistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("Householder"))
    DistMatrix<F,MD,STAR> t(A.Grid());
    DistMatrix<Base<F>,MD,STAR> d(A.Grid());
    Householder( A, t, d );
    MakeTriangular( LOWER, A );
}

} // namespace lq
} // namespace El

#endif // ifndef EL_LQ_HOUSEHOLDER_HPP
