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
Base<F> LatticeLogPotential( const Matrix<F>& R )
{
    EL_DEBUG_CSE
    typedef Base<F> Real;
    const Int m = R.Height();
    const Int n = R.Width();
    const Int minDim = Min(m,n);

    // TODO: Carefully check this
    Real logPotential=0;
    for( Int j=0; j<minDim; ++j )
        logPotential += 2*(n-j)*Log(Abs(R(j,j)));
    return logPotential;
}

#define PROTO(F) template Base<F> LatticeLogPotential( const Matrix<F>& R );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
