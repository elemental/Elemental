/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El-lite.hpp>
#include <El/blas_like/level2.hpp>
#include <El/blas_like/level3.hpp>

#include "./Symm/LL.hpp"
#include "./Symm/LU.hpp"
#include "./Symm/RL.hpp"
#include "./Symm/RU.hpp"

namespace El {

template<typename T>
void Symm
( LeftOrRight side, UpperOrLower uplo,
  T alpha, const Matrix<T>& A,
           const Matrix<T>& B,
  T beta,        Matrix<T>& C,
  bool conjugate )
{
    DEBUG_CSE
    if( side == LEFT && B.Width() == 1 )
    {
        Symv( uplo, alpha, A, B, beta, C , conjugate );
        return;
    }

    const char sideChar = LeftOrRightToChar( side );
    const char uploChar = UpperOrLowerToChar( uplo );
    if( conjugate )
    {
        blas::Hemm
        ( sideChar, uploChar, C.Height(), C.Width(),
          alpha, A.LockedBuffer(), A.LDim(),
                 B.LockedBuffer(), B.LDim(),
          beta,  C.Buffer(),       C.LDim() );
    }
    else
    {
        blas::Symm
        ( sideChar, uploChar, C.Height(), C.Width(),
          alpha, A.LockedBuffer(), A.LDim(),
                 B.LockedBuffer(), B.LDim(),
          beta,  C.Buffer(),       C.LDim() );
    }
}

template<typename T>
void Symm
( LeftOrRight side, UpperOrLower uplo,
  T alpha, const AbstractDistMatrix<T>& A,
           const AbstractDistMatrix<T>& B,
  T beta,        AbstractDistMatrix<T>& C,
  bool conjugate )
{
    DEBUG_CSE
    if( side == LEFT && B.Width() == 1 )
    {
        Symv( uplo, alpha, A, B, beta, C , conjugate );
        return;
    }

    C *= beta;
    if( side == LEFT && uplo == LOWER )
        symm::LL( alpha, A, B, C, conjugate );
    else if( side == LEFT )
        symm::LU( alpha, A, B, C, conjugate );
    else if( uplo == LOWER )
        symm::RL( alpha, A, B, C, conjugate );
    else
        symm::RU( alpha, A, B, C, conjugate );
}

#define PROTO(T) \
  template void Symm \
  ( LeftOrRight side, UpperOrLower uplo, \
    T alpha, const Matrix<T>& A, \
             const Matrix<T>& B, \
    T beta,        Matrix<T>& C, \
    bool conjugate ); \
  template void Symm \
  ( LeftOrRight side, UpperOrLower uplo, \
    T alpha, const AbstractDistMatrix<T>& A, \
             const AbstractDistMatrix<T>& B, \
    T beta,        AbstractDistMatrix<T>& C, \
    bool conjugate ); \
  template void symm::LocalAccumulateLL \
  ( Orientation orientation, T alpha, \
    const DistMatrix<T>& A, \
    const DistMatrix<T,MC,  STAR>& B_MC_STAR, \
    const DistMatrix<T,STAR,MR  >& BTrans_STAR_MR, \
          DistMatrix<T,MC,  STAR>& Z_MC_STAR, \
          DistMatrix<T,MR,  STAR>& Z_MR_STAR ); \
  template void symm::LocalAccumulateLU \
  ( Orientation orientation, T alpha, \
    const DistMatrix<T>& A, \
    const DistMatrix<T,MC,  STAR>& B_MC_STAR, \
    const DistMatrix<T,STAR,MR  >& BTrans_STAR_MR, \
          DistMatrix<T,MC,  STAR>& Z_MC_STAR, \
          DistMatrix<T,MR,  STAR>& Z_MR_STAR ); \
  template void symm::LocalAccumulateRL \
  ( Orientation orientation, T alpha, \
    const DistMatrix<T>& A, \
    const DistMatrix<T,STAR,MC  >& B_STAR_MC, \
    const DistMatrix<T,MR,  STAR>& BTrans_MR_STAR, \
          DistMatrix<T,MC,  STAR>& ZTrans_MC_STAR, \
          DistMatrix<T,MR,  STAR>& ZTrans_MR_STAR ); \
  template void symm::LocalAccumulateRU \
  ( Orientation orientation, T alpha, \
    const DistMatrix<T,MC,  MR  >& A, \
    const DistMatrix<T,STAR,MC  >& B_STAR_MC, \
    const DistMatrix<T,MR,  STAR>& BTrans_MR_STAR, \
          DistMatrix<T,MC,  STAR>& ZTrans_MC_STAR, \
          DistMatrix<T,MR,  STAR>& ZTrans_MR_STAR );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
