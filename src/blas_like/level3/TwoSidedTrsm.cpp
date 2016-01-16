/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

#include "./TwoSidedTrsm/Unblocked.hpp"
#include "./TwoSidedTrsm/LVar4.hpp"
#include "./TwoSidedTrsm/UVar4.hpp"

namespace El {

template<typename F> 
void TwoSidedTrsm
( UpperOrLower uplo, UnitOrNonUnit diag,
        Matrix<F>& A,
  const Matrix<F>& B )
{
    DEBUG_ONLY(CSE cse("TwoSidedTrsm"))
    if( uplo == LOWER )
        twotrsm::LVar4( diag, A, B );
    else
        twotrsm::UVar4( diag, A, B );
}

template<typename F> 
void TwoSidedTrsm
( UpperOrLower uplo, UnitOrNonUnit diag, 
        ElementalMatrix<F>& A,
  const ElementalMatrix<F>& B )
{
    DEBUG_ONLY(CSE cse("TwoSidedTrsm"))
    if( uplo == LOWER )
        twotrsm::LVar4( diag, A, B );
    else
        twotrsm::UVar4( diag, A, B );
}

template<typename F>
void TwoSidedTrsm
( UpperOrLower uplo, UnitOrNonUnit diag,
        DistMatrix<F,STAR,STAR>& A,
  const DistMatrix<F,STAR,STAR>& B )
{ TwoSidedTrsm( uplo, diag, A.Matrix(), B.LockedMatrix() ); }

namespace twotrsm {

template<typename F,typename=EnableIf<IsBlasScalar<F>>>
void ScaLAPACKHelper
( UpperOrLower uplo, UnitOrNonUnit diag,
        DistMatrix<F,MC,MR,BLOCK>& A,
  const DistMatrix<F,MC,MR,BLOCK>& B )
{
    if( diag == UNIT )
        LogicError("ScaLAPACK does not support unit-diagonal two-sided TRSM");
    // NOTE: ScaLAPACK additionally assumes that the diagonal of the triangular
    //       matrix is real and positive.
    AssertScaLAPACKSupport();
#ifdef EL_HAVE_SCALAPACK
    const Int n = A.Height();
    const int bHandle = blacs::Handle( A );
    const int context = blacs::GridInit( bHandle, A );
    auto descA = FillDesc( A, context );
    auto descB = FillDesc( B, context );
    const char uploChar = UpperOrLowerToChar( uplo );
    scalapack::TwoSidedTrsm
    ( uploChar, n, A.Buffer(), descA.data(), B.LockedBuffer(), descB.data() ); 
    blacs::FreeGrid( context );
    blacs::FreeHandle( bHandle );
#endif
}

template<typename F,typename=DisableIf<IsBlasScalar<F>>,typename=void>
void ScaLAPACKHelper
( UpperOrLower uplo, UnitOrNonUnit diag,
        DistMatrix<F,MC,MR,BLOCK>& A,
  const DistMatrix<F,MC,MR,BLOCK>& B )
{
    LogicError("ScaLAPACK does not support this datatype");
}

} // namespace twotrsm

template<typename F>
void TwoSidedTrsm
( UpperOrLower uplo, UnitOrNonUnit diag,
        DistMatrix<F,MC,MR,BLOCK>& A,
  const DistMatrix<F,MC,MR,BLOCK>& B )
{
    DEBUG_ONLY(CSE cse("TwoSidedTrsm"))
    twotrsm::ScaLAPACKHelper( uplo, diag, A, B );
}

#define PROTO(F) \
  template void TwoSidedTrsm \
  ( UpperOrLower uplo, UnitOrNonUnit diag, \
          Matrix<F>& A, \
    const Matrix<F>& B ); \
  template void TwoSidedTrsm \
  ( UpperOrLower uplo, UnitOrNonUnit diag, \
          ElementalMatrix<F>& A, \
    const ElementalMatrix<F>& B ); \
  template void TwoSidedTrsm \
  ( UpperOrLower uplo, UnitOrNonUnit diag, \
          DistMatrix<F,STAR,STAR>& A, \
    const DistMatrix<F,STAR,STAR>& B ); \
  template void TwoSidedTrsm \
  ( UpperOrLower uplo, UnitOrNonUnit diag, \
          DistMatrix<F,MC,MR,BLOCK>& A, \
    const DistMatrix<F,MC,MR,BLOCK>& B );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

} // namespace El
