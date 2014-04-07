/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_QUASITRSM_HPP
#define ELEM_QUASITRSM_HPP

#include ELEM_QUASITRSV_INC

namespace elem {

template<typename F,Dist XColDist,Dist XRowDist>
inline void
LocalQuasiTrsm
( LeftOrRight side, UpperOrLower uplo,
  Orientation orientation, 
  F alpha, const DistMatrix<F,STAR,STAR>& A,
                 DistMatrix<F,XColDist,XRowDist>& X,
  bool checkIfSingular=false )
{
    DEBUG_ONLY(
        CallStackEntry cse("LocalQuasiTrsm");
        if( (side == LEFT && XColDist != STAR) ||
            (side == RIGHT && XRowDist != STAR) )
            LogicError
            ("Dist of RHS must conform with that of triangle");
    )
    QuasiTrsm
    ( side, uplo, orientation, 
      alpha, A.LockedMatrix(), X.Matrix(), checkIfSingular );
}

} // namespace elem

#include "./QuasiTrsm/LLN.hpp"
#include "./QuasiTrsm/LLT.hpp"
#include "./QuasiTrsm/LUN.hpp"
#include "./QuasiTrsm/LUT.hpp"
//#include "./QuasiTrsm/RLN.hpp"
//#include "./QuasiTrsm/RLT.hpp"
//#include "./QuasiTrsm/RUN.hpp"
//#include "./QuasiTrsm/RUT.hpp"

namespace elem {

template<typename F>
inline void
QuasiTrsm
( LeftOrRight side, UpperOrLower uplo, Orientation orientation, 
  F alpha, const Matrix<F>& A, Matrix<F>& B,
  bool checkIfSingular=false )
{
    DEBUG_ONLY(
        CallStackEntry cse("QuasiTrsm");
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
    // Call the single right-hand side algorithm if appropriate
    if( side == LEFT && B.Width() == 1 )
    {
        QuasiTrsv( uplo, orientation, A, B );
        return;
    }
    // TODO: Compute appropriate transpose/conjugation options to convert
    //       to Trsv.
    /*
    else if( side == RIGHT && B.Height() == 1 )
    {
        QuasiTrsv( uplo, orientation, A, B );
        return;
    }
    */

    if( side == LEFT && uplo == LOWER )
    {
        if( orientation == NORMAL )
            quasitrsm::LLN( A, B, checkIfSingular );
        else
            quasitrsm::LLT( orientation, A, B, checkIfSingular );
    }
    else if( side == LEFT && uplo == UPPER )
    {
        if( orientation == NORMAL )
            quasitrsm::LUN( A, B, checkIfSingular );
        else
            quasitrsm::LUT( orientation, A, B, checkIfSingular );
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
QuasiTrsm
( LeftOrRight side, UpperOrLower uplo, Orientation orientation, 
  F alpha, const DistMatrix<F>& A, DistMatrix<F>& B,
  bool checkIfSingular=false )
{
    DEBUG_ONLY(
        CallStackEntry cse("QuasiTrsm");
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
    // Call the single right-hand side algorithm if appropriate
    if( side == LEFT && B.Width() == 1 )
    {
        QuasiTrsv( uplo, orientation, A, B );
        return;
    }
    // TODO: Compute appropriate transpose/conjugation options to convert
    //       to Trsv.
    /*
    else if( side == RIGHT && B.Height() == 1 )
    {
        QuasiTrsv( uplo, orientation, A, B );
        return;
    }
    */

    const Int p = B.Grid().Size();
    if( side == LEFT && uplo == LOWER )
    {
        if( orientation == NORMAL )
        {
            if( B.Width() > 5*p )
                quasitrsm::LLNLarge( A, B, checkIfSingular );
            else
                quasitrsm::LLNMedium( A, B, checkIfSingular );
        }
        else
        {
            if( B.Width() > 5*p )
                quasitrsm::LLTLarge( orientation, A, B, checkIfSingular );
            else
                quasitrsm::LLTMedium( orientation, A, B, checkIfSingular );
        }
    }
    else if( side == LEFT && uplo == UPPER )
    {
        if( orientation == NORMAL )
        {
            if( B.Width() > 5*p )
                quasitrsm::LUNLarge( A, B, checkIfSingular );
            else
                quasitrsm::LUNMedium( A, B, checkIfSingular );
        }
        else
        {
            if( B.Width() > 5*p )
                quasitrsm::LUTLarge( orientation, A, B, checkIfSingular );
            else
                quasitrsm::LUTMedium( orientation, A, B, checkIfSingular );
        }
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

} // namespace elem

#endif // ifndef ELEM_QUASITRSM_HPP
