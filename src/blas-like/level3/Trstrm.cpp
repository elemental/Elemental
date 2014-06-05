/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

#include "./Trstrm/LLN.hpp"

namespace El {

template<typename F>
void Trstrm
( LeftOrRight side, UpperOrLower uplo,
  Orientation orientation, UnitOrNonUnit diag,
  F alpha, const Matrix<F>& A, Matrix<F>& X,
  bool checkIfSingular )
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
            trstrm::LLN( diag, alpha, A, X, checkIfSingular );
        else
            LogicError("This option not yet implemented");
        /*
            trstrm::LLT
            ( orientation, diag, alpha, A, X, checkIfSingular );
        */
    }
    else
        LogicError("This option not yet implemented");
    /*
    else if( side == LEFT && uplo == UPPER )
    {
        if( orientation == NORMAL )
            trstrm::LUN( diag, alpha, A, X, checkIfSingular );
        else
            trstrm::LUT
            ( orientation, diag, alpha, A, X, checkIfSingular );
    }
    else if( side == RIGHT && uplo == LOWER )
    {
        if( orientation == NORMAL )
            trstrm::RLN( diag, alpha, A, X, checkIfSingular );
        else
            trstrm::RLT
            ( orientation, diag, alpha, A, X, checkIfSingular );
    }
    else if( side == RIGHT && uplo == UPPER )
    {
        if( orientation == NORMAL )
            trstrm::RUN( diag, alpha, A, X, checkIfSingular );
        else
            trstrm::RUT
            ( orientation, diag, alpha, A, X, checkIfSingular );
    }
    */
}

template<typename F>
void Trstrm
( LeftOrRight side, UpperOrLower uplo, 
  Orientation orientation, UnitOrNonUnit diag,
  F alpha, const DistMatrix<F>& A, DistMatrix<F>& X,
  bool checkIfSingular )
{
    DEBUG_ONLY(CallStackEntry cse("Trstrm"))
    if( side == LEFT && uplo == LOWER )
    {
        if( orientation == NORMAL )
            trstrm::LLN( diag, alpha, A, X, checkIfSingular );
        else
            LogicError("This option not yet implemented");
        /*
            trstrm::LLT
            ( orientation, diag, alpha, A, X, checkIfSingular );
        */
    }
    else
        LogicError("This option not yet implemented");
    /*
    else if( side == LEFT && uplo == UPPER )
    {
        if( orientation == NORMAL )
            trstrm::LUN( diag, alpha, A, X, checkIfSingular );
        else
            trstrm::LUT
            ( orientation, diag, alpha, A, X, checkIfSingular );
    }
    else if( side == RIGHT && uplo == LOWER )
    {
        if( orientation == NORMAL )
            trstrm::RLN( diag, alpha, A, X, checkIfSingular );
        else
            trstrm::RLT
            ( orientation, diag, alpha, A, X, checkIfSingular );
    }
    else if( side == RIGHT && uplo == UPPER )
    {
        if( orientation == NORMAL )
            trstrm::RUN( diag, alpha, A, X, checkIfSingular );
        else
            trstrm::RUT
            ( orientation, diag, alpha, A, X, checkIfSingular );
    }
    */
}

#define PROTO(F) \
  template void Trstrm \
  ( LeftOrRight side, UpperOrLower uplo, \
    Orientation orientation, UnitOrNonUnit diag, \
    F alpha, const Matrix<F>& A, Matrix<F>& X, bool checkIfSingular ); \
  template void Trstrm \
  ( LeftOrRight side, UpperOrLower uplo, \
    Orientation orientation, UnitOrNonUnit diag, \
    F alpha, const DistMatrix<F>& A, DistMatrix<F>& X, bool checkIfSingular );

PROTO(float)
PROTO(double)
PROTO(Complex<float>)
PROTO(Complex<double>)

} // namespace El
