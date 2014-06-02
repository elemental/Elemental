/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_UNIFORM_HPP
#define EL_UNIFORM_HPP

namespace El {

// Draw each entry from a uniform PDF over a closed ball.

template<typename T>
void MakeUniform( Matrix<T>& A, T center=0, Base<T> radius=1 );

template<typename T>
void MakeUniform( AbstractDistMatrix<T>& A, T center=0, Base<T> radius=1 );

template<typename T>
void MakeUniform( AbstractBlockDistMatrix<T>& A, T center=0, Base<T> radius=1 );

template<typename T>
void Uniform( Matrix<T>& A, Int m, Int n, T center=0, Base<T> radius=1 );

template<typename T>
void Uniform
( AbstractDistMatrix<T>& A, Int m, Int n, T center=0, Base<T> radius=1 );

template<typename T>
void Uniform
( AbstractBlockDistMatrix<T>& A, Int m, Int n, T center=0, Base<T> radius=1 );

} // namespace El

#endif // ifndef EL_UNIFORM_HPP
