/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_AXPYTRIANGLE_HPP
#define EL_AXPYTRIANGLE_HPP

namespace El {

template<typename T>
void AxpyTriangle
( UpperOrLower uplo, T alpha, const Matrix<T>& X, Matrix<T>& Y );

template<typename T,Dist U,Dist V>
void AxpyTriangle
( UpperOrLower uplo, T alpha, 
  const DistMatrix<T,U,V>& X, DistMatrix<T,U,V>& Y );

template<typename T>
void AxpyTriangle
( UpperOrLower uplo, T alpha, 
  const AbstractDistMatrix<T>& X, AbstractDistMatrix<T>& Y );

template<typename T>
inline void AxpyTriangle
( UpperOrLower uplo, Base<T> alpha, const Matrix<T>& X, Matrix<T>& Y )
{ AxpyTriangle( uplo, T(alpha), X, Y ); }

template<typename T,Dist U,Dist V>
inline void AxpyTriangle
( UpperOrLower uplo, Base<T> alpha, 
  const DistMatrix<T,U,V>& X, DistMatrix<T,U,V>& Y )
{ AxpyTriangle( uplo, T(alpha), X, Y ); }

template<typename T>
inline void AxpyTriangle
( UpperOrLower uplo, Base<T> alpha, 
  const AbstractDistMatrix<T>& X, AbstractDistMatrix<T>& Y )
{ AxpyTriangle( uplo, T(alpha), X, Y ); }

} // namespace El

#endif // ifndef EL_AXPYTRIANGLE_HPP
