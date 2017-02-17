/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El-lite.hpp>
#include <El/blas_like/level1.hpp>
#include <El/blas_like/level3.hpp>

#include "./Trdtrmm/Unblocked.hpp"
#include "./Trdtrmm/LVar1.hpp"
#include "./Trdtrmm/UVar1.hpp"

namespace El {

template<typename F>
void Trdtrmm( UpperOrLower uplo, Matrix<F>& A, bool conjugate )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
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
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( A.Height() != A.Width() )
          LogicError("A must be square");
    )
    if( uplo == LOWER )
        trdtrmm::LVar1( A, dOff, conjugate );
    else
        LogicError("Not yet written");
}

template<typename F>
void Trdtrmm( UpperOrLower uplo, AbstractDistMatrix<F>& A, bool conjugate )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
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
  AbstractDistMatrix<F>& A, const AbstractDistMatrix<F>& dOff, 
  bool conjugate )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( A.Height() != A.Width() )
          LogicError("A must be square");
    )
    if( uplo == LOWER )
        trdtrmm::LVar1( A, dOff, conjugate );
    else
        LogicError("Not yet written");
}

template<typename F>
void Trdtrmm
( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A, bool conjugate )
{ Trdtrmm( uplo, A.Matrix(), conjugate ); }

template<typename F>
void Trdtrmm
( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A, 
  const DistMatrix<F,STAR,STAR>& dOff, bool conjugate )
{ Trdtrmm( uplo, A.Matrix(), dOff.LockedMatrix(), conjugate ); }

#define PROTO(F) \
  template void Trdtrmm( UpperOrLower uplo, Matrix<F>& A, bool conjugate ); \
  template void Trdtrmm \
  ( UpperOrLower uplo, Matrix<F>& A, const Matrix<F>& dOff, bool conjugate ); \
  template void Trdtrmm \
  ( UpperOrLower uplo, AbstractDistMatrix<F>& A, bool conjugate ); \
  template void Trdtrmm \
  ( UpperOrLower uplo, \
    AbstractDistMatrix<F>& A, const AbstractDistMatrix<F>& dOff, \
    bool conjugate ); \
  template void Trdtrmm \
  ( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A, bool conjugate ); \
  template void Trdtrmm \
  ( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A, \
    const DistMatrix<F,STAR,STAR>& dOff, bool conjugate );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
