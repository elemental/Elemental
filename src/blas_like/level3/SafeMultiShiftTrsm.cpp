/*
   Copyright (c) 2009-2016, Jack Poulson and Tim Moon
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

#include "./SafeMultiShiftTrsm/Overflow.hpp"
#include "./SafeMultiShiftTrsm/LUN.hpp"

namespace El {

template<typename F>
void SafeMultiShiftTrsm
( LeftOrRight side, UpperOrLower uplo, Orientation orientation,
  F alpha, Matrix<F>& A, const Matrix<F>& shifts,
  Matrix<F>& B, Matrix<F>& scales )
{
    DEBUG_ONLY(CSE cse("SafeMultiShiftTrsm"))
    B *= alpha;
    if( side == LEFT && uplo == UPPER && orientation == NORMAL)
    {
        safemstrsm::LUN( A, shifts, B, scales );
    }
    else
    {
        LogicError("This option is not yet supported");
    }
}

template<typename F>
void SafeMultiShiftTrsm
( LeftOrRight side, UpperOrLower uplo, Orientation orientation,
  F alpha, const ElementalMatrix<F>& A, const ElementalMatrix<F>& shifts, 
  ElementalMatrix<F>& B, ElementalMatrix<F>& scales )
{
    DEBUG_ONLY(CSE cse("SafeMultiShiftTrsm"))
    B *= alpha;
    if( side == LEFT && uplo == UPPER && orientation == NORMAL)
    {
        safemstrsm::LUN( A, shifts, B, scales );
    }
    else
        LogicError("This option is not yet supported");
}
  
#define PROTO(F) \
  template void SafeMultiShiftTrsm \
  ( LeftOrRight side, UpperOrLower uplo, Orientation orientation, \
    F alpha, Matrix<F>& A, const Matrix<F>& shifts, \
    Matrix<F>& B, Matrix<F>& scales ); \
  template void SafeMultiShiftTrsm \
  ( LeftOrRight side, UpperOrLower uplo, Orientation orientation, \
    F alpha, const ElementalMatrix<F>& A, const ElementalMatrix<F>& shifts, \
    ElementalMatrix<F>& B, ElementalMatrix<F>& scales );
  
#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
