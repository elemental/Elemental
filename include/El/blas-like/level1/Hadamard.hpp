/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_HADAMARD_HPP
#define EL_HADAMARD_HPP

// C(i,j) := A(i,j) B(i,j)

namespace El {

template<typename T> 
void Hadamard( const Matrix<T>& A, const Matrix<T>& B, Matrix<T>& C );

template<typename T,Dist U,Dist V> 
void Hadamard
( const DistMatrix<T,U,V>& A, 
  const DistMatrix<T,U,V>& B, 
        DistMatrix<T,U,V>& C );

template<typename T>
void Hadamard
( const AbstractDistMatrix<T>& A,
  const AbstractDistMatrix<T>& B,
        AbstractDistMatrix<T>& C );

} // namespace El

#endif // ifndef EL_HADAMARD_HPP
