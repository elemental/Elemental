/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_ADJOINT_HPP
#define EL_ADJOINT_HPP

namespace El {

template<typename T>
void Adjoint( const Matrix<T>& A, Matrix<T>& B );

template<typename T>
void Adjoint( const AbstractDistMatrix<T>& A, AbstractDistMatrix<T>& B );

template<typename T>
void Adjoint
( const AbstractBlockDistMatrix<T>& A, AbstractBlockDistMatrix<T>& B );

} // namespace El

#endif // ifndef EL_ADJOINT_HPP
