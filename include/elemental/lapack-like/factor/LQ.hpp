/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LQ_HPP
#define ELEM_LQ_HPP

#include "./LQ/ApplyQ.hpp"
#include "./LQ/Householder.hpp"
#include "./LQ/SolveAfter.hpp"
#include "./LQ/Explicit.hpp"

namespace elem {

template<typename F> 
inline void
LQ( Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("LQ"))
    lq::Householder( A );
}

template<typename F> 
inline void
LQ( DistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("LQ"))
    lq::Householder( A );
}

template<typename F> 
inline void
LQ( Matrix<F>& A, Matrix<F>& t, Matrix<Base<F>>& d )
{
    DEBUG_ONLY(CallStackEntry cse("LQ"))
    lq::Householder( A, t, d );
}

template<typename F> 
inline void
LQ( DistMatrix<F>& A, DistMatrix<F,MD,STAR>& t, DistMatrix<Base<F>,MD,STAR>& d )
{
    DEBUG_ONLY(CallStackEntry cse("LQ"))
    lq::Householder( A, t, d );
}

// Variants which perform (Businger-Golub) row-pivoting
// ====================================================
// TODO

} // namespace elem

#endif // ifndef ELEM_LQ_HPP
