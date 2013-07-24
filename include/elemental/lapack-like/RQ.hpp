/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_RQ_HPP
#define ELEM_LAPACK_RQ_HPP

#include "elemental/lapack-like/RQ/ApplyQ.hpp"
#include "elemental/lapack-like/RQ/Cholesky.hpp"
#include "elemental/lapack-like/RQ/Householder.hpp"

namespace elem {

template<typename F> 
inline void
RQ( Matrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("RQ");
#endif
    rq::Householder( A );
}

template<typename F> 
inline void
RQ( DistMatrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("RQ");
#endif
    rq::Householder( A );
}

template<typename F> 
inline void
RQ( Matrix<F>& A, Matrix<F>& t )
{
#ifndef RELEASE
    CallStackEntry entry("RQ");
#endif
    rq::Householder( A, t );
}

template<typename F> 
inline void
RQ( DistMatrix<F>& A, DistMatrix<F,MD,STAR>& t )
{
#ifndef RELEASE
    CallStackEntry entry("RQ");
#endif
    rq::Householder( A, t );
}

// TODO: BusingerGolub pivoting?

} // namespace elem

#endif // ifndef ELEM_LAPACK_RQ_HPP
