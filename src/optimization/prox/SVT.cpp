/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

#include "./SVT/Normal.hpp"
#include "./SVT/Cross.hpp"
#include "./SVT/PivotedQR.hpp"
#include "./SVT/TSQR.hpp"

namespace El {

template<typename Field>
Int SVT( Matrix<Field>& A, const Base<Field>& tau, bool relative )
{
    EL_DEBUG_CSE
    return svt::Normal( A, tau, relative );
}

template<typename Field>
Int SVT( AbstractDistMatrix<Field>& A, const Base<Field>& tau, bool relative )
{
    EL_DEBUG_CSE
    // NOTE: This should be less accurate (but faster) than svt::Normal
    return svt::Cross( A, tau, relative );
}

template<typename Field>
Int SVT
( Matrix<Field>& A, const Base<Field>& tau, Int relaxedRank, bool relative )
{
    EL_DEBUG_CSE
    // Preprocess with numSteps iterations of pivoted QR factorization
    return svt::PivotedQR( A, tau, relaxedRank, relative );
}

template<typename Field>
Int SVT
( AbstractDistMatrix<Field>& A, const Base<Field>& tau,
  Int relaxedRank, bool relative )
{
    EL_DEBUG_CSE
    // Preprocess with numSteps iterations of pivoted QR factorization
    return svt::PivotedQR( A, tau, relaxedRank, relative );
}

// Singular-value soft-thresholding based on TSQR
template<typename Field,Dist U>
Int SVT( DistMatrix<Field,U,STAR>& A, const Base<Field>& tau, bool relative )
{
    EL_DEBUG_CSE
    return svt::TSQR( A, tau, relative );
}

#define PROTO_DIST(Field,U) \
  template Int SVT \
  ( DistMatrix<Field,U,STAR>& A, const Base<Field>& tau, bool relative );

#define PROTO(Field) \
  template Int SVT \
  ( Matrix<Field>& A, const Base<Field>& tau, bool relative ); \
  template Int SVT \
  ( AbstractDistMatrix<Field>& A, const Base<Field>& tau, bool relative ); \
  template Int SVT \
  ( Matrix<Field>& A, const Base<Field>& tau, \
    Int relaxedRank, bool relative ); \
  template Int SVT \
  ( AbstractDistMatrix<Field>& A, const Base<Field>& tau, \
    Int relaxedRank, bool relative ); \
  template Int svt::Cross \
  ( Matrix<Field>& A, const Base<Field>& tau, bool relative ); \
  template Int svt::Cross \
  ( AbstractDistMatrix<Field>& A, const Base<Field>& tau, bool relative ); \
  template Int svt::Cross \
  ( DistMatrix<Field,VC,STAR>& A, const Base<Field>& tau, bool relative ); \
  template Int svt::PivotedQR \
  ( Matrix<Field>& A, const Base<Field>& tau, Int numSteps, bool relative ); \
  template Int svt::PivotedQR \
  ( AbstractDistMatrix<Field>& A, const Base<Field>& tau, Int numSteps, \
    bool relative ); \
  template Int svt::TSQR \
  ( AbstractDistMatrix<Field>& A, const Base<Field>& tau, bool relative ); \
  PROTO_DIST(Field,MC  ) \
  PROTO_DIST(Field,MD  ) \
  PROTO_DIST(Field,MR  ) \
  PROTO_DIST(Field,STAR) \
  PROTO_DIST(Field,VC  ) \
  PROTO_DIST(Field,VR  )

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
