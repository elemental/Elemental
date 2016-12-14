/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El-lite.hpp>
#include <El/blas_like/level3.hpp>

#include "./TwoSidedTrmm/Unblocked.hpp"
#include "./TwoSidedTrmm/LVar4.hpp"
#include "./TwoSidedTrmm/UVar4.hpp"

namespace El {

template<typename T> 
void TwoSidedTrmm
( UpperOrLower uplo, UnitOrNonUnit diag,
        Matrix<T>& A,
  const Matrix<T>& B )
{
    EL_DEBUG_CSE
    if( uplo == LOWER )
        twotrmm::LVar4( diag, A, B );
    else
        twotrmm::UVar4( diag, A, B );
}

template<typename T> 
void TwoSidedTrmm
( UpperOrLower uplo, UnitOrNonUnit diag, 
        AbstractDistMatrix<T>& A,
  const AbstractDistMatrix<T>& B )
{
    EL_DEBUG_CSE
    if( uplo == LOWER )
        twotrmm::LVar4( diag, A, B );
    else
        twotrmm::UVar4( diag, A, B );
}

template<typename T>
void TwoSidedTrmm
( UpperOrLower uplo, UnitOrNonUnit diag,
        DistMatrix<T,STAR,STAR>& A,
  const DistMatrix<T,STAR,STAR>& B )
{ TwoSidedTrmm( uplo, diag, A.Matrix(), B.LockedMatrix() ); }

namespace twotrmm {

template<typename T,typename=EnableIf<IsBlasScalar<T>>>
void ScaLAPACKHelper
( UpperOrLower uplo, UnitOrNonUnit diag,
        DistMatrix<T,MC,MR,BLOCK>& A,
  const DistMatrix<T,MC,MR,BLOCK>& B )
{
    if( diag == UNIT )
        LogicError("ScaLAPACK does not support unit-diagonal two-sided TRMM");
    // NOTE: ScaLAPACK additionally assumes that the diagonal of the triangular
    //       matrix is real and positive.
    AssertScaLAPACKSupport();
#ifdef EL_HAVE_SCALAPACK
    const Int n = A.Height();
    const char uploChar = UpperOrLowerToChar( uplo );

    auto descA = FillDesc( A );
    auto descB = FillDesc( B );
    scalapack::TwoSidedTrmm
    ( uploChar, n, A.Buffer(), descA.data(), B.LockedBuffer(), descB.data() );
#endif
}

template<typename T,typename=DisableIf<IsBlasScalar<T>>,typename=void>
void ScaLAPACKHelper
( UpperOrLower uplo, UnitOrNonUnit diag,
        DistMatrix<T,MC,MR,BLOCK>& A,
  const DistMatrix<T,MC,MR,BLOCK>& B )
{
    LogicError("ScaLAPACK does not support this datatype");
}

} // namespace twotrmm

template<typename T>
void TwoSidedTrmm
( UpperOrLower uplo, UnitOrNonUnit diag,
        DistMatrix<T,MC,MR,BLOCK>& A,
  const DistMatrix<T,MC,MR,BLOCK>& B )
{
    EL_DEBUG_CSE
    twotrmm::ScaLAPACKHelper( uplo, diag, A, B );
}

#define PROTO(T) \
  template void TwoSidedTrmm \
  ( UpperOrLower uplo, UnitOrNonUnit diag, \
          Matrix<T>& A, \
    const Matrix<T>& B ); \
  template void TwoSidedTrmm \
  ( UpperOrLower uplo, UnitOrNonUnit diag, \
          AbstractDistMatrix<T>& A, \
    const AbstractDistMatrix<T>& B ); \
  template void TwoSidedTrmm \
  ( UpperOrLower uplo, UnitOrNonUnit diag, \
          DistMatrix<T,STAR,STAR>& A, \
    const DistMatrix<T,STAR,STAR>& B ); \
  template void TwoSidedTrmm \
  ( UpperOrLower uplo, UnitOrNonUnit diag, \
          DistMatrix<T,MC,MR,BLOCK>& A, \
    const DistMatrix<T,MC,MR,BLOCK>& B );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
