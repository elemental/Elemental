/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_QR_HOUSEHOLDER_HPP
#define EL_QR_HOUSEHOLDER_HPP

#include "./ApplyQ.hpp"
#include "./PanelHouseholder.hpp"

namespace El {
namespace qr {

template<typename F> 
inline void
Householder( Matrix<F>& A, Matrix<F>& t, Matrix<Base<F>>& d )
{
    DEBUG_ONLY(CallStackEntry cse("qr::Householder"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);
    t.Resize( minDim, 1 );
    d.Resize( minDim, 1 );

    const Int bsize = Blocksize();
    for( Int k=0; k<minDim; k+=bsize )
    {
        const Int nb = Min(bsize,minDim-k);
        auto AB1 = ViewRange( A, k, k,    m, k+nb );
        auto AB2 = ViewRange( A, k, k+nb, m, n    ); 
        auto t1 = View( t, k, 0, nb, 1 );
        auto d1 = View( d, k, 0, nb, 1 );

        PanelHouseholder( AB1, t1, d1 );
        ApplyQ( LEFT, ADJOINT, AB1, t1, d1, AB2 );
    }
}

template<typename F> 
inline void
Householder( Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("qr::Householder"))
    Matrix<F> t;
    Matrix<Base<F>> d;
    Householder( A, t, d );
    MakeTriangular( UPPER, A );
}

template<typename F> 
inline void
Householder
( DistMatrix<F>& A, DistMatrix<F,MD,STAR>& t, DistMatrix<Base<F>,MD,STAR>& d )
{
    DEBUG_ONLY(
        CallStackEntry cse("qr::Householder");
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
        auto AB1 = ViewRange( A, k, k,    m, k+nb );
        auto AB2 = ViewRange( A, k, k+nb, m, n    ); 
        auto t1 = View( t, k, 0, nb, 1 );
        auto d1 = View( d, k, 0, nb, 1 );

        PanelHouseholder( AB1, t1, d1 );
        ApplyQ( LEFT, ADJOINT, AB1, t1, d1, AB2 );
    }
}

template<typename F> 
inline void
Householder( DistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("qr::Householder"))
    DistMatrix<F,MD,STAR> t(A.Grid());
    DistMatrix<Base<F>,MD,STAR> d(A.Grid());
    Householder( A, t, d );
    MakeTriangular( UPPER, A );
}

} // namespace qr
} // namespace El

#endif // ifndef EL_QR_HOUSEHOLDER_HPP
