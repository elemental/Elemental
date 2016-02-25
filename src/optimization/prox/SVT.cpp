/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

#include "./SVT/Normal.hpp"
#include "./SVT/Cross.hpp"
#include "./SVT/PivotedQR.hpp"
#include "./SVT/TSQR.hpp"

namespace El {

template<typename F>
Int SVT( Matrix<F>& A, Base<F> tau, bool relative )
{
    DEBUG_ONLY(CSE cse("SVT"))
    return svt::Normal( A, tau, relative );
}

template<typename F>
Int SVT( ElementalMatrix<F>& A, Base<F> tau, bool relative )
{
    DEBUG_ONLY(CSE cse("SVT"))
    // NOTE: This should be less accurate (but faster) than svt::Normal
    return svt::Cross( A, tau, relative );
}

template<typename F>
Int SVT( Matrix<F>& A, Base<F> tau, Int relaxedRank, bool relative )
{
    DEBUG_ONLY(CSE cse("SVT"))
    // Preprocess with numSteps iterations of pivoted QR factorization
    return svt::PivotedQR( A, tau, relaxedRank, relative );
}

template<typename F>
Int SVT( ElementalMatrix<F>& A, Base<F> tau, Int relaxedRank, bool relative )
{
    DEBUG_ONLY(CSE cse("SVT"))
    // Preprocess with numSteps iterations of pivoted QR factorization
    return svt::PivotedQR( A, tau, relaxedRank, relative );
}

// Singular-value soft-thresholding based on TSQR
template<typename F,Dist U>
Int SVT( DistMatrix<F,U,STAR>& A, Base<F> tau, bool relative )
{
    DEBUG_ONLY(CSE cse("SVT"))
    return svt::TSQR( A, tau, relative );
}

#define PROTO_DIST(F,U) \
  template Int SVT( DistMatrix<F,U,STAR>& A, Base<F> tau, bool relative );

#define PROTO(F) \
  template Int SVT( Matrix<F>& A, Base<F> tau, bool relative ); \
  template Int SVT( ElementalMatrix<F>& A, Base<F> tau, bool relative ); \
  template Int SVT \
  ( Matrix<F>& A, Base<F> tau, Int relaxedRank, bool relative ); \
  template Int SVT \
  ( ElementalMatrix<F>& A, Base<F> tau, Int relaxedRank, bool relative ); \
  template Int svt::Cross \
  ( Matrix<F>& A, Base<F> tau, bool relative ); \
  template Int svt::Cross \
  ( ElementalMatrix<F>& A, Base<F> tau, bool relative ); \
  template Int svt::Cross \
  ( DistMatrix<F,VC,STAR>& A, Base<F> tau, bool relative ); \
  template Int svt::PivotedQR \
  ( Matrix<F>& A, Base<F> tau, Int numSteps, bool relative ); \
  template Int svt::PivotedQR \
  ( ElementalMatrix<F>& A, Base<F> tau, Int numSteps, bool relative ); \
  template Int svt::TSQR \
  ( ElementalMatrix<F>& A, Base<F> tau, bool relative ); \
  PROTO_DIST(F,MC  ) \
  PROTO_DIST(F,MD  ) \
  PROTO_DIST(F,MR  ) \
  PROTO_DIST(F,STAR) \
  PROTO_DIST(F,VC  ) \
  PROTO_DIST(F,VR  )

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
