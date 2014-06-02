/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_COPY_HPP
#define EL_COPY_HPP

namespace El {

template<typename T>
void Copy( const Matrix<T>& A, Matrix<T>& B );

template<typename Real>
void Copy( const Matrix<Real>& A, Matrix<Complex<Real>>& B );

template<typename T,Dist U,Dist V,Dist W,Dist Z>
void Copy( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B );

template<typename Real,Dist U,Dist V,Dist W,Dist Z>
void Copy( const DistMatrix<Real,U,V>& A, DistMatrix<Complex<Real>,W,Z>& B );

template<typename T,Dist U,Dist V,Dist W,Dist Z>
void Copy( const BlockDistMatrix<T,U,V>& A, BlockDistMatrix<T,W,Z>& B );

template<typename Real,Dist U,Dist V,Dist W,Dist Z>
void Copy
( const BlockDistMatrix<Real,U,V>& A, BlockDistMatrix<Complex<Real>,W,Z>& B );

template<typename T>
void Copy( const AbstractDistMatrix<T>& A, AbstractDistMatrix<T>& B );

} // namespace El

#endif // ifndef EL_COPY_HPP
