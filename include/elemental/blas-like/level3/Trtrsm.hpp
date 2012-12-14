/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

#include "./Trtrsm/LLN.hpp"

namespace elem {

template<typename F>
inline void
Trtrsm
( LeftOrRight side, UpperOrLower uplo,
  Orientation orientation, UnitOrNonUnit diag,
  F alpha, const Matrix<F>& A, Matrix<F>& X,
  bool checkIfSingular )
{
#ifndef RELEASE
    PushCallStack("Trtrsm");
    if( A.Height() != A.Width() || X.Height() != X.Width() )
        throw std::logic_error("Triangular matrices must be square");
    if( A.Height() != X.Height() )
        throw std::logic_error("Nonconformal Trtrsm");
#endif
    if( side == LEFT && uplo == LOWER )
    {
        if( orientation == NORMAL )
            internal::TrtrsmLLN( diag, alpha, A, X, checkIfSingular );
        else
            throw std::logic_error("This option not yet implemented");
        /*
            internal::TrtrsmLLT
            ( orientation, diag, alpha, A, X, checkIfSingular );
        */
    }
    else
        throw std::logic_error("This option not yet implemented");
    /*
    else if( side == LEFT && uplo == UPPER )
    {
        if( orientation == NORMAL )
            internal::TrtrsmLUN( diag, alpha, A, X, checkIfSingular );
        else
            internal::TrtrsmLUT
            ( orientation, diag, alpha, A, X, checkIfSingular );
    }
    else if( side == RIGHT && uplo == LOWER )
    {
        if( orientation == NORMAL )
            internal::TrtrsmRLN( diag, alpha, A, X, checkIfSingular );
        else
            internal::TrtrsmRLT
            ( orientation, diag, alpha, A, X, checkIfSingular );
    }
    else if( side == RIGHT && uplo == UPPER )
    {
        if( orientation == NORMAL )
            internal::TrtrsmRUN( diag, alpha, A, X, checkIfSingular );
        else
            internal::TrtrsmRUT
            ( orientation, diag, alpha, A, X, checkIfSingular );
    }
    */
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
Trtrsm
( LeftOrRight side, UpperOrLower uplo, 
  Orientation orientation, UnitOrNonUnit diag,
  F alpha, const DistMatrix<F>& A, DistMatrix<F>& X,
  bool checkIfSingular )
{
#ifndef RELEASE
    PushCallStack("Trtrsm");
#endif
    if( side == LEFT && uplo == LOWER )
    {
        if( orientation == NORMAL )
            internal::TrtrsmLLN( diag, alpha, A, X, checkIfSingular );
        else
            throw std::logic_error("This option not yet implemented");
        /*
            internal::TrtrsmLLT
            ( orientation, diag, alpha, A, X, checkIfSingular );
        */
    }
    else
        throw std::logic_error("This option not yet implemented");
    /*
    else if( side == LEFT && uplo == UPPER )
    {
        if( orientation == NORMAL )
            internal::TrtrsmLUN( diag, alpha, A, X, checkIfSingular );
        else
            internal::TrtrsmLUT
            ( orientation, diag, alpha, A, X, checkIfSingular );
    }
    else if( side == RIGHT && uplo == LOWER )
    {
        if( orientation == NORMAL )
            internal::TrtrsmRLN( diag, alpha, A, X, checkIfSingular );
        else
            internal::TrtrsmRLT
            ( orientation, diag, alpha, A, X, checkIfSingular );
    }
    else if( side == RIGHT && uplo == UPPER )
    {
        if( orientation == NORMAL )
            internal::TrtrsmRUN( diag, alpha, A, X, checkIfSingular );
        else
            internal::TrtrsmRUT
            ( orientation, diag, alpha, A, X, checkIfSingular );
    }
    */
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
