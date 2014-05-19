/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_DISPLAY_HPP
#define EL_DISPLAY_HPP

namespace El {

void ProcessEvents( int numMsecs );

template<typename T>
void Display( const Matrix<T>& A, std::string title="Default" );

template<typename T>
void Display( const Matrix<Complex<T>>& A, std::string title="Default" );

template<typename T,Dist U,Dist V>
void Display( const DistMatrix<T,U,V>& A, std::string title="Default" );

template<typename T,Dist U,Dist V>
void Display( const BlockDistMatrix<T,U,V>& A, std::string title="Default" );

template<typename T>
void Display( const AbstractDistMatrix<T>& AAbs, std::string title="" );

template<typename T>
void Display( const AbstractBlockDistMatrix<T>& AAbs, std::string title="" );

} // namespace El

#endif // ifndef EL_DISPLAY_HPP
