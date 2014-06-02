/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_MAKESYMMETRIC_HPP
#define EL_MAKESYMMETRIC_HPP

namespace El {

template<typename T>
void MakeSymmetric( UpperOrLower uplo, Matrix<T>& A, bool conjugate=false );

template<typename T,Dist U,Dist V>
void MakeSymmetric
( UpperOrLower uplo, DistMatrix<T,U,V>& A, bool conjugate=false );

template<typename F>
void MakeSymmetric
( UpperOrLower uplo, AbstractDistMatrix<F>& A, bool conjugate=false );

} // namespace El

#endif // ifndef EL_MAKESYMMETRIC_HPP
