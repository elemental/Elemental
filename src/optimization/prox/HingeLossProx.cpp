/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename Real>
void HingeLossProx( Matrix<Real>& A, Real tau )
{
    DEBUG_ONLY(CSE cse("HingeLossProx"))
    auto hingeProx = 
      [=]( Real alpha ) -> Real
      { if( alpha < 1 ) { return Min(alpha+1/tau,Real(1)); }
        else            { return alpha;                    } };
    EntrywiseMap( A, function<Real(Real)>(hingeProx) );
}

template<typename Real>
void HingeLossProx( AbstractDistMatrix<Real>& A, Real tau )
{
    DEBUG_ONLY(CSE cse("HingeLossProx"))
    auto hingeProx = 
      [=]( Real alpha ) -> Real
      { if( alpha < 1 ) { return Min(alpha+1/tau,Real(1)); }
        else            { return alpha;                    } };
    EntrywiseMap( A, function<Real(Real)>(hingeProx) );
}

#define PROTO(Real) \
  template void HingeLossProx( Matrix<Real>& A, Real tau ); \
  template void HingeLossProx( AbstractDistMatrix<Real>& A, Real tau );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

} // namespace El
