/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_SPY_HPP
#define EL_SPY_HPP

namespace El {

template<typename T>
void Spy( const Matrix<T>& A, std::string title="Default", Base<T> tol=0 );

template<typename T,Distribution U,Distribution V>
void Spy
( const DistMatrix<T,U,V>& A, std::string title="Default", Base<T> tol=0 );

template<typename T,Dist U,Dist V>
void Spy
( const BlockDistMatrix<T,U,V>& A, std::string title="Default", Base<T> tol=0 );

template<typename T>
void Spy
( const AbstractDistMatrix<T>& A, std::string title="Default", Base<T> tol=0 );

template<typename T>
void Spy
( const AbstractBlockDistMatrix<T>& A, 
  std::string title="Default", Base<T> tol=0 );

} // namespace El

#endif // ifndef EL_SPY_HPP
