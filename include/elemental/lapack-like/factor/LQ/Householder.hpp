/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LQ_HOUSEHOLDER_HPP
#define ELEM_LQ_HOUSEHOLDER_HPP

#include "./ApplyQ.hpp"
#include "./PanelHouseholder.hpp"

namespace elem {
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
        auto ATopPan    = ViewRange( A, k,    k, k+nb, n );
        auto ABottomPan = ViewRange( A, k+nb, k, m,    n );
        auto t1 = View( t, k, 0, nb, 1 );
        auto d1 = View( d, k, 0, nb, 1 );

        PanelHouseholder( ATopPan, t1, d1 );
        ApplyQ( RIGHT, ADJOINT, ATopPan, t1, d1, ABottomPan );
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
( DistMatrix<F>& A, DistMatrix<F,MD,STAR>& t, DistMatrix<Base<F>,MD,STAR>& d )
{
    DEBUG_ONLY(
        CallStackEntry cse("Householder");
        if( A.Grid() != t.Grid() || t.Grid() != d.Grid() )
            LogicError("{A,t,d} must be distributed over the same grid");
    )
    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);

    t.SetRoot( A.DiagonalRoot() );
    d.SetRoot( A.DiagonalRoot() );
    t.AlignCols( A.DiagonalAlign() );
    d.AlignCols( A.DiagonalAlign() );
    t.Resize( minDim, 1 );
    d.Resize( minDim, 1 );

    const Int bsize = Blocksize();
    for( Int k=0; k<minDim; k+=bsize )
    {
        const Int nb = Min(bsize,minDim-k);
        auto ATopPan    = ViewRange( A, k,    k, k+nb, n );
        auto ABottomPan = ViewRange( A, k+nb, k, m,    n );
        auto t1 = View( t, k, 0, nb, 1 );
        auto d1 = View( d, k, 0, nb, 1 );

        PanelHouseholder( ATopPan, t1, d1 );
        ApplyQ( RIGHT, ADJOINT, ATopPan, t1, d1, ABottomPan );
    }
}

template<typename F> 
inline void
Householder( DistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("Householder"))
    DistMatrix<F,MD,STAR> t(A.Grid());
    DistMatrix<Base<F>,MD,STAR> d(A.Grid());
    Householder( A, t, d );
    MakeTriangular( LOWER, A );
}

} // namespace lq
} // namespace elem

#endif // ifndef ELEM_LQ_HOUSEHOLDER_HPP
