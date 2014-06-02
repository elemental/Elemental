/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_CONJUGATE_HPP
#define EL_CONJUGATE_HPP

namespace El {

template<typename T>
inline void Conjugate( Matrix<T>& A ) { }

template<typename T>
void Conjugate( Matrix<Complex<T>>& A );

template<typename T>
void Conjugate( const Matrix<T>& A, Matrix<T>& B );

template<typename T>
void Conjugate( AbstractDistMatrix<T>& A );

template<typename T,Dist U,Dist V,Dist W,Dist Z>
void Conjugate( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B );
template<typename T>
void Conjugate( const AbstractDistMatrix<T>& A, AbstractDistMatrix<T>& B );

} // namespace El

#endif // ifndef EL_CONJUGATE_HPP
