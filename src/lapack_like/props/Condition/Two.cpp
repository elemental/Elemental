/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {

template<typename F> 
Base<F> TwoCondition( const Matrix<F>& A )
{
    DEBUG_CSE
    typedef Base<F> Real;
    Matrix<Real> s;
    SVD( A, s );

    Real cond = 1;
    const Int numVals = s.Height();
    if( numVals > 0 )
        cond = s(0) / s(numVals-1);
    return cond;
}

template<typename F> 
Base<F> TwoCondition( const ElementalMatrix<F>& A )
{
    DEBUG_CSE
    typedef Base<F> Real;
    DistMatrix<Real,VR,STAR> s( A.Grid() );
    SVD( A, s );

    Real cond = 1;
    const Int numVals = s.Height();
    if( numVals > 0 )
        cond = s.Get(0,0) / s.Get(numVals-1,0);
    return cond;
}

#define PROTO(F) \
  template Base<F> TwoCondition( const Matrix<F>& A ); \
  template Base<F> TwoCondition( const ElementalMatrix<F>& A );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
