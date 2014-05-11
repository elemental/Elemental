/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_RQ_HPP
#define ELEM_RQ_HPP

#include "./RQ/ApplyQ.hpp"
#include "./RQ/Cholesky.hpp"
#include "./RQ/Householder.hpp"
#include "./RQ/SolveAfter.hpp"

namespace elem {

template<typename F> 
inline void
RQ( Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("RQ"))
    rq::Householder( A );
}

template<typename F> 
inline void
RQ( DistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("RQ"))
    rq::Householder( A );
}

template<typename F> 
inline void
RQ( Matrix<F>& A, Matrix<F>& t, Matrix<Base<F>>& d )
{
    DEBUG_ONLY(CallStackEntry cse("RQ"))
    rq::Householder( A, t, d );
}

template<typename F> 
inline void
RQ( DistMatrix<F>& A, DistMatrix<F,MD,STAR>& t, DistMatrix<Base<F>,MD,STAR>& d )
{
    DEBUG_ONLY(CallStackEntry cse("RQ"))
    rq::Householder( A, t, d );
}

} // namespace elem

#endif // ifndef ELEM_RQ_HPP
