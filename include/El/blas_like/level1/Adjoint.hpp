/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_BLAS_ADJOINT_HPP
#define EL_BLAS_ADJOINT_HPP

namespace El {

template<typename T>
void Adjoint( const Matrix<T>& A, Matrix<T>& B )
{
    DEBUG_ONLY(CSE cse("Adjoint"))
    Transpose( A, B, true );
}

template<typename T>
void Adjoint( const ElementalMatrix<T>& A, ElementalMatrix<T>& B )
{
    DEBUG_ONLY(CSE cse("Adjoint"))
    Transpose( A, B, true );
}

template<typename T>
void Adjoint
( const BlockMatrix<T>& A, BlockMatrix<T>& B )
{
    DEBUG_ONLY(CSE cse("Adjoint"))
    Transpose( A, B, true );
}

template<typename T>
void Adjoint
( const AbstractDistMatrix<T>& A, AbstractDistMatrix<T>& B )
{
    DEBUG_ONLY(CSE cse("Adjoint"))
    Transpose( A, B, true );
}

template<typename T>
void Adjoint( const SparseMatrix<T>& A, SparseMatrix<T>& B )
{
    DEBUG_ONLY(CSE cse("Adjoint"))
    Transpose( A, B, true );
}

template<typename T>
void Adjoint( const DistSparseMatrix<T>& A, DistSparseMatrix<T>& B )
{
    DEBUG_ONLY(CSE cse("Adjoint"))
    Transpose( A, B, true );
}

} // namespace El

#endif // ifndef EL_BLAS_ADJOINT_HPP
