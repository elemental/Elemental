/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_MULTISHIFTTRSM_HPP
#define ELEM_MULTISHIFTTRSM_HPP

#include "./MultiShiftTrsm/LUN.hpp"
#include "./MultiShiftTrsm/LUT.hpp"

namespace elem {

template<typename F>
inline void
MultiShiftTrsm
( LeftOrRight side, UpperOrLower uplo, Orientation orientation,
  F alpha, const Matrix<F>& U, const Matrix<F>& shifts, Matrix<F>& X )
{
    DEBUG_ONLY(CallStackEntry cse("MultiShiftTrsm"))
    if( side == LEFT && uplo == UPPER )
    {
        if( orientation == NORMAL )
            mstrsm::LUN( alpha, U, shifts, X );
        else
            mstrsm::LUT( orientation, alpha, U, shifts, X );
    }
    else
        LogicError("This option is not yet supported");
}

template<typename F>
inline void
MultiShiftTrsm
( LeftOrRight side, UpperOrLower uplo, Orientation orientation,
  F alpha, const DistMatrix<F>& U, const DistMatrix<F,VR,STAR>& shifts, 
  DistMatrix<F>& X )
{
    DEBUG_ONLY(CallStackEntry cse("MultiShiftTrsm"))
    if( side == LEFT && uplo == UPPER )
    {
        if( orientation == NORMAL )
            mstrsm::LUN( alpha, U, shifts, X );
        else
            mstrsm::LUT( orientation, alpha, U, shifts, X );
    }
    else
        LogicError("This option is not yet supported");
}

} // namespace elem

#endif // ifndef ELEM_MULTISHIFTTRSM_HPP
