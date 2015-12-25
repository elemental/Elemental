/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_BLAS_ADJOINTCONTRACT_HPP
#define EL_BLAS_ADJOINTCONTRACT_HPP

namespace El {

template<typename T>
void AdjointContract( const ElementalMatrix<T>& A, ElementalMatrix<T>& B )
{
    DEBUG_ONLY(CSE cse("AdjointContract"))
    TransposeContract( A, B, true );
}

template<typename T>
void AdjointContract
( const BlockMatrix<T>& A, 
        BlockMatrix<T>& B )
{
    DEBUG_ONLY(CSE cse("AdjointContract"))
    TransposeContract( A, B, true );
}

} // namespace El

#endif // ifndef EL_BLAS_ADJOINTCONTRACT_HPP
