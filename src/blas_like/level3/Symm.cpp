/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

#include "./Symm/LL.hpp"
#include "./Symm/LU.hpp"
#include "./Symm/RL.hpp"
#include "./Symm/RU.hpp"

namespace El {

template<typename Ring>
void Symm
( LeftOrRight side, UpperOrLower uplo,
  Ring alpha, const Matrix<Ring>& A, 
              const Matrix<Ring>& B, 
  Ring beta,        Matrix<Ring>& C,
  bool conjugate )
{
    DEBUG_ONLY(CSE cse("Symm"))
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

template<typename Ring>
void Symm
( LeftOrRight side, UpperOrLower uplo,
  Ring alpha, const AbstractDistMatrix<Ring>& A, 
              const AbstractDistMatrix<Ring>& B,
  Ring beta,        AbstractDistMatrix<Ring>& C, 
  bool conjugate )
{
    DEBUG_ONLY(CSE cse("Symm"))
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

#define PROTO(Ring) \
  template void Symm \
  ( LeftOrRight side, UpperOrLower uplo, \
    Ring alpha, const Matrix<Ring>& A, \
                const Matrix<Ring>& B, \
    Ring beta,        Matrix<Ring>& C, bool conjugate ); \
  template void Symm \
  ( LeftOrRight side, UpperOrLower uplo, \
    Ring alpha, const AbstractDistMatrix<Ring>& A, \
                const AbstractDistMatrix<Ring>& B, \
    Ring beta,        AbstractDistMatrix<Ring>& C, bool conjugate ); \
  template void symm::LocalAccumulateLL \
  ( Orientation orientation, Ring alpha, \
    const DistMatrix<Ring>& A, \
    const DistMatrix<Ring,MC,  STAR>& B_MC_STAR, \
    const DistMatrix<Ring,STAR,MR  >& BTrans_STAR_MR, \
          DistMatrix<Ring,MC,  STAR>& Z_MC_STAR, \
          DistMatrix<Ring,MR,  STAR>& Z_MR_STAR ); \
  template void symm::LocalAccumulateLU \
  ( Orientation orientation, Ring alpha, \
    const DistMatrix<Ring>& A, \
    const DistMatrix<Ring,MC,  STAR>& B_MC_STAR, \
    const DistMatrix<Ring,STAR,MR  >& BTrans_STAR_MR, \
          DistMatrix<Ring,MC,  STAR>& Z_MC_STAR, \
          DistMatrix<Ring,MR,  STAR>& Z_MR_STAR ); \
  template void symm::LocalAccumulateRL \
  ( Orientation orientation, Ring alpha, \
    const DistMatrix<Ring>& A, \
    const DistMatrix<Ring,STAR,MC  >& B_STAR_MC, \
    const DistMatrix<Ring,MR,  STAR>& BTrans_MR_STAR, \
          DistMatrix<Ring,MC,  STAR>& ZTrans_MC_STAR, \
          DistMatrix<Ring,MR,  STAR>& ZTrans_MR_STAR ); \
  template void symm::LocalAccumulateRU \
  ( Orientation orientation, Ring alpha, \
    const DistMatrix<Ring,MC,  MR  >& A, \
    const DistMatrix<Ring,STAR,MC  >& B_STAR_MC, \
    const DistMatrix<Ring,MR,  STAR>& BTrans_MR_STAR, \
          DistMatrix<Ring,MC,  STAR>& ZTrans_MC_STAR, \
          DistMatrix<Ring,MR,  STAR>& ZTrans_MR_STAR );

#define EL_ENABLE_QUAD
#include "El/macros/Instantiate.h"

} // namespace El
