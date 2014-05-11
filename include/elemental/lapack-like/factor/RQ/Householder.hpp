/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_RQ_HOUSEHOLDER_HPP
#define ELEM_RQ_HOUSEHOLDER_HPP

#include ELEM_MAKETRAPEZOIDAL_INC

#include "./ApplyQ.hpp"
#include "./PanelHouseholder.hpp"

namespace elem {
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
        auto ATopPan    = View( A, 0,  0, ki, kj+nb ); 
        auto ABottomPan = View( A, ki, 0, nb, kj+nb );
        auto t1 = View( t, k, 0, nb, 1 );
        auto d1 = View( d, k, 0, nb, 1 );

        PanelHouseholder( ABottomPan, t1, d1 );
        ApplyQ( RIGHT, ADJOINT, ABottomPan, t1, d1, ATopPan );
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
( DistMatrix<F>& A, DistMatrix<F,MD,STAR>& t, DistMatrix<Base<F>,MD,STAR>& d )
{
    DEBUG_ONLY(
        CallStackEntry cse("rq::Householder");
        if( A.Grid() != t.Grid() || t.Grid() != d.Grid() )
            LogicError("{A,t,d} must be distributed over the same grid");
    )
    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);
    const Int offset = n-m;
    const Int iOff = m-minDim;
    const Int jOff = n-minDim;

    t.SetRoot( A.DiagonalRoot(offset) );
    d.SetRoot( A.DiagonalRoot(offset) );
    t.AlignCols( A.DiagonalAlign(offset) );
    d.AlignCols( A.DiagonalAlign(offset) );
    t.Resize( minDim, 1 );
    d.Resize( minDim, 1 );

    const Int bsize = Blocksize();
    const Int kLast = LastOffset( minDim, bsize );
    for( Int k=kLast; k>=0; k-=bsize )
    {
        const Int nb = Min(bsize,minDim-k);
        const Int ki = k + iOff;
        const Int kj = k + jOff;
        auto ATopPan    = View( A, 0,  0, ki, kj+nb ); 
        auto ABottomPan = View( A, ki, 0, nb, kj+nb );
        auto t1 = View( t, k, 0, nb, 1 );
        auto d1 = View( d, k, 0, nb, 1 );

        PanelHouseholder( ABottomPan, t1, d1 );
        ApplyQ( RIGHT, ADJOINT, ABottomPan, t1, d1, ATopPan );
    }
}

template<typename F> 
inline void
Householder( DistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("rq::Householder"))
    DistMatrix<F,MD,STAR> t(A.Grid());
    DistMatrix<Base<F>,MD,STAR> d(A.Grid());
    Householder( A, t, d );
    MakeTrapezoidal( UPPER, A, Min(A.Height(),A.Width()) );
}

} // namespace rq
} // namespace elem

#endif // ifndef ELEM_RQ_HOUSEHOLDER_HPP
