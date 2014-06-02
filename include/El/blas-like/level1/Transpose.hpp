/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_TRANSPOSE_HPP
#define EL_TRANSPOSE_HPP

namespace El {

template<typename T>
void Transpose( const Matrix<T>& A, Matrix<T>& B, bool conjugate=false );

template<typename T,Dist U,Dist V,Dist W,Dist Z>
void Transpose
( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B, bool conjugate=false );

template<typename T,Dist U,Dist V,Dist W,Dist Z>
void Transpose
( const BlockDistMatrix<T,U,V>& A, BlockDistMatrix<T,W,Z>& B, 
  bool conjugate=false );

template<typename T>
void Transpose
( const AbstractDistMatrix<T>& A, AbstractDistMatrix<T>& B, 
  bool conjugate=false );
template<typename T>
void Transpose
( const AbstractBlockDistMatrix<T>& A, AbstractBlockDistMatrix<T>& B, 
  bool conjugate=false );

} // namespace El

#endif // ifndef EL_TRANSPOSE_HPP
