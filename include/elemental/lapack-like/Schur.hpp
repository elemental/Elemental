/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_SCHUR_HPP
#define ELEM_LAPACK_SCHUR_HPP

#include "elemental/lapack-like/Schur/InverseFreeSDC.hpp"
#include "elemental/lapack-like/Schur/SDC.hpp"

namespace elem {

template<typename F>
inline void
Schur( Matrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry cse("Schur");
#endif
    throw std::logic_error("This routine not yet written");
    // TODO: Call LAPACK...
}

template<typename F>
inline int
Schur( DistMatrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry cse("Schur");
#endif
    throw std::logic_error("This routine not yet written");
    // TODO: Spectral D&C
}

} // namespace elem

#endif // ifndef ELEM_LAPACK_SCHUR_HPP
