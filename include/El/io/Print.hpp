/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_PRINT_HPP
#define EL_PRINT_HPP

namespace El {

template<typename T>
void Print
( const Matrix<T>& A, std::string title="", std::ostream& os=std::cout );

template<typename T,Dist U,Dist V>
void Print
( const DistMatrix<T,U,V>& A, 
  std::string title="", std::ostream& os=std::cout );

template<typename T,Dist U,Dist V>
void Print
( const BlockDistMatrix<T,U,V>& A, 
  std::string title="", std::ostream& os=std::cout );

template<typename T>
void Print
( const AbstractDistMatrix<T>& AAbs, std::string title="", 
  std::ostream& os=std::cout );

template<typename T>
void Print
( const AbstractBlockDistMatrix<T>& AAbs, std::string title="", 
  std::ostream& os=std::cout );

// Utilities
// =========
template<typename T>
void Print
( const std::vector<T>& x, std::string title="", std::ostream& os=std::cout );

} // namespace El

#endif // ifndef EL_PRINT_HPP
