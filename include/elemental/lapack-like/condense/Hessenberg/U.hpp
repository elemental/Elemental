/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_HESSENBERG_U_HPP
#define ELEM_HESSENBERG_U_HPP

#include "./UUnb.hpp"

namespace elem {
namespace hessenberg {

template<typename F>
inline void U( Matrix<F>& A, Matrix<F>& t )
{
    DEBUG_ONLY(CallStackEntry cse("hessenberg::U"))
    // TODO: Blocked algorithm
    UUnb( A, t );
}

template<typename F> 
inline void U( DistMatrix<F>& A, DistMatrix<F,STAR,STAR>& t )
{
    DEBUG_ONLY(CallStackEntry cse("hessenberg::U"))
    // TODO: Blocked algorithm
    UUnb( A, t );
}

} // namespace hessenberg
} // namespace elem

#endif // ifndef ELEM_HESSENBERG_U_HPP
