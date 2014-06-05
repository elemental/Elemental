/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

#include "./Trdtrmm/Unblocked.hpp"
#include "./Trdtrmm/LVar1.hpp"
#include "./Trdtrmm/UVar1.hpp"

namespace El {

template<typename F>
void Trdtrmm( UpperOrLower uplo, Matrix<F>& A, bool conjugate )
{
    DEBUG_ONLY(
        CallStackEntry cse("Trdtrdmm");
        if( A.Height() != A.Width() )
            LogicError("A must be square");
    )
    if( uplo == LOWER )
        trdtrmm::LVar1( A, conjugate );
    else
        trdtrmm::UVar1( A, conjugate );
}

template<typename F>
void Trdtrmm
( UpperOrLower uplo, Matrix<F>& A, const Matrix<F>& dOff, bool conjugate )
{
    DEBUG_ONLY(
        CallStackEntry cse("Trdtrdmm");
        if( A.Height() != A.Width() )
            LogicError("A must be square");
    )
    if( uplo == LOWER )
        trdtrmm::LVar1( A, dOff, conjugate );
    else
        LogicError("Not yet written");
}

template<typename F>
void Trdtrmm( UpperOrLower uplo, DistMatrix<F>& A, bool conjugate )
{
    DEBUG_ONLY(
        CallStackEntry cse("Trdtrmm");
        if( A.Height() != A.Width() )
            LogicError("A must be square");
    )
    if( uplo == LOWER )
        trdtrmm::LVar1( A, conjugate );
    else
        trdtrmm::UVar1( A, conjugate );
}

template<typename F>
void Trdtrmm
( UpperOrLower uplo, 
  DistMatrix<F>& A, const DistMatrix<F,MD,STAR>& dOff, bool conjugate )
{
    DEBUG_ONLY(
        CallStackEntry cse("Trdtrmm");
        if( A.Height() != A.Width() )
            LogicError("A must be square");
    )
    if( uplo == LOWER )
        trdtrmm::LVar1( A, dOff, conjugate );
    else
        LogicError("Not yet written");
}

#define PROTO(F) \
  template void Trdtrmm( UpperOrLower uplo, Matrix<F>& A, bool conjugate ); \
  template void Trdtrmm \
  ( UpperOrLower uplo, Matrix<F>& A, const Matrix<F>& dOff, bool conjugate ); \
  template void Trdtrmm \
  ( UpperOrLower uplo, DistMatrix<F>& A, bool conjugate ); \
  template void Trdtrmm \
  ( UpperOrLower uplo, DistMatrix<F>& A, const DistMatrix<F,MD,STAR>& dOff, \
    bool conjugate );

PROTO(float)
PROTO(double)
PROTO(Complex<float>)
PROTO(Complex<double>)

} // namespace El
