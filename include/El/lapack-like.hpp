/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_LAPACK_HPP
#define EL_LAPACK_HPP

namespace El {

// Ones
// ====
template<typename T>
void Ones( Matrix<T>& A, Int m, Int n );
template<typename T>
void Ones( AbstractDistMatrix<T>& A, Int m, Int n );
template<typename T>
void Ones( AbstractBlockDistMatrix<T>& A, Int m, Int n );

// Zeros
// =====
template<typename T>
void Zeros( Matrix<T>& A, Int m, Int n );
template<typename T>
void Zeros( AbstractDistMatrix<T>& A, Int m, Int n );
template<typename T>
void Zeros( AbstractBlockDistMatrix<T>& A, Int m, Int n );

} // namespace El

#include "./lapack-like/util.hpp"
#include "./lapack-like/perm.hpp"

#include "./lapack-like/factor.hpp"
#include "./lapack-like/condense.hpp"
#include "./lapack-like/funcs.hpp"

#include "./lapack-like/decomp.hpp"
#include "./lapack-like/solve.hpp"

#include "./lapack-like/props.hpp"

#endif // ifndef EL_LAPACK_HPP
