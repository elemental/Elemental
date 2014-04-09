/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_MULTISHIFTQUASITRSM_HPP
#define ELEM_MULTISHIFTQUASITRSM_HPP

namespace elem {

template<typename F,Dist shiftColDist,
                    Dist     XColDist,Dist XRowDist>
inline void
LocalMultiShiftQuasiTrsm
( LeftOrRight side, UpperOrLower uplo, Orientation orientation, 
  F alpha, const DistMatrix<F,STAR,STAR>& A,
           const DistMatrix<F,shiftColDist,STAR    >& shifts,
                 DistMatrix<F,    XColDist,XRowDist>& X )
{
    DEBUG_ONLY(
        CallStackEntry cse("LocalMultiShiftQuasiTrsm");
        if( (side == LEFT &&  (     XColDist != STAR || 
                                shiftColDist != XRowDist) ) ||
            (side == RIGHT && (     XRowDist != STAR || 
                                shiftColDist != XColDist) ) )
            LogicError
            ("Dist of RHS and shifts must conform with that of triangle");
    )
    // NOTE: Is this prototype available yet?!?
    MultiShiftQuasiTrsm
    ( side, uplo, orientation, 
      alpha, A.LockedMatrix(), shifts.LockedMatrix(), X.Matrix() );
}

} // namespace elem

#include "./MultiShiftQuasiTrsm/LLN.hpp"
#include "./MultiShiftQuasiTrsm/LLT.hpp"
#include "./MultiShiftQuasiTrsm/LUN.hpp"
#include "./MultiShiftQuasiTrsm/LUT.hpp"
//#include "./MultiShiftQuasiTrsm/RLN.hpp"
//#include "./MultiShiftQuasiTrsm/RLT.hpp"
//#include "./MultiShiftQuasiTrsm/RUN.hpp"
//#include "./MultiShiftQuasiTrsm/RUT.hpp"

namespace elem {

template<typename F>
inline void
MultiShiftQuasiTrsm
( LeftOrRight side, UpperOrLower uplo, Orientation orientation, 
  F alpha, const Matrix<F>& A, const Matrix<F>& shifts, Matrix<F>& B )
{
    DEBUG_ONLY(
        CallStackEntry cse("MultiShiftQuasiTrsm");
        if( A.Height() != A.Width() )
            LogicError("A must be square");
        if( side == LEFT )
        {
            if( A.Height() != B.Height() )
                LogicError("Nonconformal Trsm");
        }
        else
        {
            if( A.Height() != B.Width() )
                LogicError("Nonconformal Trsm");
        }
    )
    Scale( alpha, B );
    // TODO: Call the single right-hand side algorithm if appropriate

    if( side == LEFT && uplo == LOWER )
    {
        if( orientation == NORMAL )
            msquasitrsm::LLN( A, shifts, B );
        else
            msquasitrsm::LLT( orientation, A, shifts, B );
    }
    else if( side == LEFT && uplo == UPPER )
    {
        if( orientation == NORMAL )
            msquasitrsm::LUN( A, shifts, B );
        else
            msquasitrsm::LUT( orientation, A, shifts, B );
    }
    else if( side == RIGHT && uplo == LOWER )
    {
        if( orientation == NORMAL )
            //quasitrsm::RLN( A, B, checkIfSingular );
            LogicError("This case not yet handled");
        else
            //quasitrsm::RLT( orientation, A, B, checkIfSingular );
            LogicError("This case not yet handled");
    }
    else if( side == RIGHT && uplo == UPPER )
    {
        if( orientation == NORMAL )
            //quasitrsm::RUN( A, B, checkIfSingular );
            LogicError("This case not yet handled");
        else
            //quasitrsm::RUT( orientation, A, B, checkIfSingular );
            LogicError("This case not yet handled");
    }
}

template<typename F>
inline void
MultiShiftQuasiTrsm
( LeftOrRight side, UpperOrLower uplo, Orientation orientation, 
  F alpha, const DistMatrix<F>& A, const DistMatrix<F,VR,STAR>& shifts, 
  DistMatrix<F>& B )
{
    DEBUG_ONLY(
        CallStackEntry cse("MultiShiftQuasiTrsm");
        if( A.Grid() != B.Grid() )
            LogicError("A and B must use the same grid");
        if( A.Height() != A.Width() )
            LogicError("A must be square");
        if( side == LEFT )
        {
            if( A.Height() != B.Height() )
                LogicError("Nonconformal Trsm");
        }
        else
        {
            if( A.Height() != B.Width() )
                LogicError("Nonconformal Trsm");
        }
    )
    Scale( alpha, B );
    // TODO: Call the single right-hand side algorithm if appropriate

    const Int p = B.Grid().Size();
    if( side == LEFT && uplo == LOWER )
    {
        if( orientation == NORMAL )
        {
            if( B.Width() > 5*p )
                msquasitrsm::LLNLarge( A, shifts, B );
            else
                msquasitrsm::LLNMedium( A, shifts, B );
        }
        else
        {
            if( B.Width() > 5*p )
                msquasitrsm::LLTLarge( orientation, A, shifts, B );
            else
                msquasitrsm::LLTMedium( orientation, A, shifts, B );
        }
    }
    else if( side == LEFT && uplo == UPPER )
    {
        if( orientation == NORMAL )
        {
            if( B.Width() > 5*p )
                msquasitrsm::LUNLarge( A, shifts, B );
            else
                msquasitrsm::LUNMedium( A, shifts, B );
        }
        else
        {
            if( B.Width() > 5*p )
                msquasitrsm::LUTLarge( orientation, A, shifts, B );
            else
                msquasitrsm::LUTMedium( orientation, A, shifts, B );
        }
    }
    else if( side == RIGHT && uplo == LOWER )
    {
        if( orientation == NORMAL )
            //msquasitrsm::RLN( A, shifts, B );
            LogicError("This case not yet handled");
        else
            //msquasitrsm::RLT( orientation, A, shifts, B );
            LogicError("This case not yet handled");
    }
    else if( side == RIGHT && uplo == UPPER )
    {
        if( orientation == NORMAL )
            //msquasitrsm::RUN( A, shifts, B );
            LogicError("This case not yet handled");
        else
            //msquasitrsm::RUT( orientation, A, shifts, B );
            LogicError("This case not yet handled");
    }
}

} // namespace elem

#endif // ifndef ELEM_MULTISHIFTQUASITRSM_HPP
