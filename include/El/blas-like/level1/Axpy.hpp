/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_AXPY_HPP
#define EL_AXPY_HPP

namespace El {

template<typename T>
void Axpy( T alpha, const Matrix<T>& X, Matrix<T>& Y );

template<typename T,Dist U,Dist V,Dist W,Dist Z>
void Axpy( T alpha, const DistMatrix<T,U,V>& X, DistMatrix<T,W,Z>& Y );

template<typename T>
void Axpy( T alpha, const AbstractDistMatrix<T>& X, AbstractDistMatrix<T>& Y );

template<typename T>
inline void Axpy( Base<T> alpha, const Matrix<T>& X, Matrix<T>& Y )
{ Axpy( T(alpha), X, Y ); }

template<typename T,Dist U,Dist V,Dist W,Dist Z>
inline void Axpy
( Base<T> alpha, const DistMatrix<T,U,V>& X, DistMatrix<T,W,Z>& Y )
{ Axpy( T(alpha), X, Y ); }

template<typename T>
inline void Axpy
( Base<T> alpha, const AbstractDistMatrix<T>& X, AbstractDistMatrix<T>& Y )
{ Axpy( T(alpha), X, Y ); }

} // namespace El

#endif // ifndef EL_AXPY_HPP
