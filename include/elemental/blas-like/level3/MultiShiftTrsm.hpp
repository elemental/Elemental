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

namespace elem {

template<typename F,Dist shiftColDist,
                    Dist     XColDist,Dist XRowDist>
inline void
LocalMultiShiftTrsm
( LeftOrRight side, UpperOrLower uplo, Orientation orientation,
  F alpha, const DistMatrix<F,STAR,STAR>& A,
           const DistMatrix<F,shiftColDist,STAR    >& shifts,
                 DistMatrix<F,    XColDist,XRowDist>& X )
{
    DEBUG_ONLY(
        CallStackEntry cse("LocalMultiShiftTrsm");
        if( (side == LEFT &&  (     XColDist != STAR ||
                                shiftColDist != XRowDist) ) ||
            (side == RIGHT && (     XRowDist != STAR ||
                                shiftColDist != XColDist) ) )
            LogicError
            ("Dist of RHS and shifts must conform with that of triangle");
    )
    // NOTE: Is this prototype available yet?!?
    MultiShiftTrsm
    ( side, uplo, orientation, 
      alpha, A.LockedMatrix(), shifts.LockedMatrix(), X.Matrix() );
}

} // namespace elem

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
