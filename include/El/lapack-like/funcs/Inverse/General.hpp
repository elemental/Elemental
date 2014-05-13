/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_INVERSE_GENERAL_HPP
#define EL_INVERSE_GENERAL_HPP

#include "./General/LUPartialPiv.hpp"

namespace El {

template<typename F> 
inline void
Inverse( Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("Inverse"))
    inverse::LUPartialPiv( A );
}

template<typename F> 
inline void
Inverse( DistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("Inverse"))
    inverse::LUPartialPiv( A );
}

template<typename F>
inline void
LocalInverse( DistMatrix<F,STAR,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("LocalInverse"))
    Inverse( A.Matrix() );
}

} // namespace El

#endif // ifndef EL_INVERSE_GENERAL_HPP
