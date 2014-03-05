/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_MULTISHIFTHESSSOLVE_HPP
#define ELEM_MULTISHIFTHESSSOLVE_HPP

namespace elem {
namespace mshs {

template<typename F>
inline void
UN( F alpha, Matrix<F>& H, const Matrix<F>& shifts, Matrix<F>& X ) 
{
    DEBUG_ONLY(CallStackEntry cse("mshs::UN"))
    Scale( alpha, X );
    // TODO
    LogicError("This routine is not yet written");
}

// NOTE: We pick a [VC,* ] distribution for the upper-Hessenberg matrix H
//       since whole columns will need to be formed on every process and
//       this distribution will keep the communication balanced.
template<typename F>
inline void
UN
( F alpha, const DistMatrix<F,VC,STAR>& H, const DistMatrix<F,VR,STAR>& shifts,
  DistMatrix<F>& X ) 
{
    DEBUG_ONLY(CallStackEntry cse("mshs::UN"))
    Scale( alpha, X );

    const Grid& g = H.Grid();
    // TODO
    LogicError("This routine is not yet written");
}

// TODO: UT, LN, and LT

} // namespace mshs

template<typename F>
inline void
MultiShiftHessSolve
( UpperOrLower uplo, Orientation orientation,
  F alpha, const Matrix<F>& H, const Matrix<F>& shifts, Matrix<F>& X )
{
    DEBUG_ONLY(CallStackEntry cse("MultiShiftHessSolve"))
    if( uplo == UPPER )
    {
        if( orientation == NORMAL )
            mshs::UN( alpha, H, shifts, X );
        else
            LogicError("This option is not yet supported");
    }
    else
        LogicError("This option is not yet supported");
}

template<typename F>
inline void
MultiShiftHessSolve
( UpperOrLower uplo, Orientation orientation,
  F alpha, const DistMatrix<F,VC,STAR>& H, const DistMatrix<F,VR,STAR>& shifts, 
  DistMatrix<F>& X )
{
    DEBUG_ONLY(CallStackEntry cse("MultiShiftHessSolve"))
    if( uplo == UPPER )
    {
        if( orientation == NORMAL )
            mshs::UN( alpha, H, shifts, X );
        else
            LogicError("This option is not yet supported");
    }
    else
        LogicError("This option is not yet supported");
}

} // namespace elem

#endif // ifndef ELEM_MULTISHIFTHESSSOLVE_HPP
