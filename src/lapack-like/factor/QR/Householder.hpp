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

        const Range<Int> ind1(     k,    k+nb ),
                         indB(     k,    m    ),
                         ind2Horz( k+nb, n    );

        auto AB1 = A( indB, ind1     );
        auto AB2 = A( indB, ind2Horz );
        auto t1 = t( ind1, IR(0,1) );
        auto d1 = d( ind1, IR(0,1) );

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
( AbstractDistMatrix<F>& APre, AbstractDistMatrix<F>& tPre, 
  AbstractDistMatrix<Base<F>>& dPre )
{
    DEBUG_ONLY(
        CallStackEntry cse("qr::Householder");
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

        const Range<Int> ind1(     k,    k+nb ),
                         indB(     k,    m    ),
                         ind2Horz( k+nb, n    );

        auto AB1 = A( indB, ind1     );
        auto AB2 = A( indB, ind2Horz );
        auto t1 = t( ind1, IR(0,1) );
        auto d1 = d( ind1, IR(0,1) );

        PanelHouseholder( AB1, t1, d1 );
        ApplyQ( LEFT, ADJOINT, AB1, t1, d1, AB2 );
    }
    Copy( A, APre, RESTORE_READ_WRITE_PROXY );
    Copy( t, tPre, RESTORE_WRITE_PROXY      );
    Copy( d, dPre, RESTORE_WRITE_PROXY      );
}

template<typename F> 
inline void
Householder( AbstractDistMatrix<F>& A )
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
