/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_TRSTRM_HPP
#define ELEM_TRSTRM_HPP

namespace elem {

template<typename F>
inline void
LocalTrstrm
( LeftOrRight side, UpperOrLower uplo,
  Orientation orientation, UnitOrNonUnit diag,
  F alpha, const DistMatrix<F,STAR,STAR>& A,
                 DistMatrix<F,STAR,STAR>& X,
  bool checkIfSingular=true )
{
    DEBUG_ONLY(CallStackEntry cse("LocalTrstrm"))
    Trstrm
    ( side, uplo, orientation, diag,
      alpha, A.LockedMatrix(), X.Matrix(), checkIfSingular );
}

} // namespace elem

#include "./Trstrm/LLN.hpp"

namespace elem {

template<typename F>
inline void
Trstrm
( LeftOrRight side, UpperOrLower uplo,
  Orientation orientation, UnitOrNonUnit diag,
  F alpha, const Matrix<F>& A, Matrix<F>& X,
  bool checkIfSingular=true )
{
    DEBUG_ONLY(
        CallStackEntry cse("Trstrm");
        if( A.Height() != A.Width() || X.Height() != X.Width() )
            LogicError("Triangular matrices must be square");
        if( A.Height() != X.Height() )
            LogicError("Nonconformal Trstrm");
    )
    if( side == LEFT && uplo == LOWER )
    {
        if( orientation == NORMAL )
            internal::TrstrmLLN( diag, alpha, A, X, checkIfSingular );
        else
            LogicError("This option not yet implemented");
        /*
            internal::TrstrmLLT
            ( orientation, diag, alpha, A, X, checkIfSingular );
        */
    }
    else
        LogicError("This option not yet implemented");
    /*
    else if( side == LEFT && uplo == UPPER )
    {
        if( orientation == NORMAL )
            internal::TrstrmLUN( diag, alpha, A, X, checkIfSingular );
        else
            internal::TrstrmLUT
            ( orientation, diag, alpha, A, X, checkIfSingular );
    }
    else if( side == RIGHT && uplo == LOWER )
    {
        if( orientation == NORMAL )
            internal::TrstrmRLN( diag, alpha, A, X, checkIfSingular );
        else
            internal::TrstrmRLT
            ( orientation, diag, alpha, A, X, checkIfSingular );
    }
    else if( side == RIGHT && uplo == UPPER )
    {
        if( orientation == NORMAL )
            internal::TrstrmRUN( diag, alpha, A, X, checkIfSingular );
        else
            internal::TrstrmRUT
            ( orientation, diag, alpha, A, X, checkIfSingular );
    }
    */
}

template<typename F>
inline void
Trstrm
( LeftOrRight side, UpperOrLower uplo, 
  Orientation orientation, UnitOrNonUnit diag,
  F alpha, const DistMatrix<F>& A, DistMatrix<F>& X,
  bool checkIfSingular=true )
{
    DEBUG_ONLY(CallStackEntry cse("Trstrm"))
    if( side == LEFT && uplo == LOWER )
    {
        if( orientation == NORMAL )
            internal::TrstrmLLN( diag, alpha, A, X, checkIfSingular );
        else
            LogicError("This option not yet implemented");
        /*
            internal::TrstrmLLT
            ( orientation, diag, alpha, A, X, checkIfSingular );
        */
    }
    else
        LogicError("This option not yet implemented");
    /*
    else if( side == LEFT && uplo == UPPER )
    {
        if( orientation == NORMAL )
            internal::TrstrmLUN( diag, alpha, A, X, checkIfSingular );
        else
            internal::TrstrmLUT
            ( orientation, diag, alpha, A, X, checkIfSingular );
    }
    else if( side == RIGHT && uplo == LOWER )
    {
        if( orientation == NORMAL )
            internal::TrstrmRLN( diag, alpha, A, X, checkIfSingular );
        else
            internal::TrstrmRLT
            ( orientation, diag, alpha, A, X, checkIfSingular );
    }
    else if( side == RIGHT && uplo == UPPER )
    {
        if( orientation == NORMAL )
            internal::TrstrmRUN( diag, alpha, A, X, checkIfSingular );
        else
            internal::TrstrmRUT
            ( orientation, diag, alpha, A, X, checkIfSingular );
    }
    */
}

} // namespace elem

#endif // ifndef ELEM_TRSTRM_HPP
