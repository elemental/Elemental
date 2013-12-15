/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_BIDIAG_HPP
#define ELEM_LAPACK_BIDIAG_HPP

#include "./Bidiag/Apply.hpp"
#include "./Bidiag/PanelL.hpp"
#include "./Bidiag/PanelU.hpp"
#include "./Bidiag/L.hpp"
#include "./Bidiag/LUnb.hpp"
#include "./Bidiag/U.hpp"
#include "./Bidiag/UUnb.hpp"

namespace elem {

template<typename F>
inline void Bidiag( Matrix<F>& A, Matrix<F>& tP, Matrix<F>& tQ )
{
    DEBUG_ONLY(CallStackEntry cse("Bidiag"))
    if( A.Height() >= A.Width() )
        bidiag::U( A, tP, tQ );
    else
        bidiag::L( A, tP, tQ );
}

template<typename F> 
inline void Bidiag
( DistMatrix<F>& A, DistMatrix<F,STAR,STAR>& tP, DistMatrix<F,STAR,STAR>& tQ )
{
    DEBUG_ONLY(CallStackEntry cse("Bidiag"))
    if( A.Height() >= A.Width() )
        bidiag::U( A, tP, tQ );
    else
        bidiag::L( A, tP, tQ );
}

template<typename F>
inline void Bidiag( Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("Bidiag"))
    Matrix<F> tP, tQ;
    Bidiag( A, tP, tQ );
}

template<typename F> 
inline void Bidiag( DistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("Bidiag"))
    DistMatrix<F,STAR,STAR> tP(A.Grid()), tQ(A.Grid());
    Bidiag( A, tP, tQ );
}

} // namespace elem

#endif // ifndef ELEM_LAPACK_BIDIAG_HPP
