/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_LQ_HOUSEHOLDER_HPP
#define ELEM_LAPACK_LQ_HOUSEHOLDER_HPP

#include "elemental/lapack-like/LQ/ApplyQ.hpp"
#include "elemental/lapack-like/LQ/PanelHouseholder.hpp"

namespace elem {
namespace lq {

// On exit, the lower triangle of A is overwritten by L, and the Householder
// transforms that determine Q are stored above the diagonal of A with an 
// implicit one on the diagonal. 
//
// In the complex case, the column-vector t stores the unit-magnitude complex 
// rotations that map the norms of the implicit Householder vectors to their
// coefficient:  
//                psi_j = 2 tau_j / ( u_j^H u_j ),
// where tau_j is the j'th entry of t and u_j is the j'th unscaled Householder
// reflector.

template<typename F> 
inline void
Householder( Matrix<F>& A, Matrix<F>& t )
{
    DEBUG_ONLY(CallStackEntry cse("lq::Householder"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);
    t.Resize( minDim, 1 );

    const Int bsize = Blocksize();
    for( Int k=0; k<minDim; k+=bsize )
    {
        const Int nb = Min(bsize,minDim-k);
        auto ATopPan    = ViewRange( A, k,    k, k+nb, n );
        auto ABottomPan = ViewRange( A, k+nb, k, m,    n );
        auto t1 = View( t, k, 0, nb, 1 );

        PanelHouseholder( ATopPan, t1 );
        ApplyQ( RIGHT, ADJOINT, ATopPan, t1, ABottomPan );
    }
}

template<typename F> 
inline void
Householder( Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("lq::Householder"))
    Matrix<F> t;
    Householder( A, t );
}

template<typename F> 
inline void
Householder( DistMatrix<F>& A, DistMatrix<F,MD,STAR>& t )
{
    DEBUG_ONLY(
        CallStackEntry cse("Householder");
        if( A.Grid() != t.Grid() )
            LogicError("{A,t} must be distributed over the same grid");
    )
    t.AlignWithDiagonal( A );
    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);
    t.Resize( minDim, 1 );

    const Int bsize = Blocksize();
    for( Int k=0; k<minDim; k+=bsize )
    {
        const Int nb = Min(bsize,minDim-k);
        auto ATopPan    = ViewRange( A, k,    k, k+nb, n );
        auto ABottomPan = ViewRange( A, k+nb, k, m,    n );
        auto t1 = View( t, k, 0, nb, 1 );

        PanelHouseholder( ATopPan, t1 );
        ApplyQ( RIGHT, ADJOINT, ATopPan, t1, ABottomPan );
    }
}

template<typename F> 
inline void
Householder( DistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("Householder"))
    DistMatrix<F,MD,STAR> t(A.Grid());
    Householder( A, t );
}

} // namespace lq
} // namespace elem

#endif // ifndef ELEM_LAPACK_LQ_HOUSEHOLDER_HPP
