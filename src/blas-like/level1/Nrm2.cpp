/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

namespace El {

template<typename F>
Base<F> Nrm2( const Matrix<F>& x )
{
    DEBUG_ONLY(
        CallStackEntry cse("Nrm2");
        if( x.Height() != 1 && x.Width() != 1 )
            LogicError("Expected vector input");
    )
    Base<F> norm;
    if( x.Width() == 1 )
        norm = blas::Nrm2( x.Height(), x.LockedBuffer(), 1 );
    else
        norm = blas::Nrm2( x.Width(), x.LockedBuffer(), x.LDim() );
    return norm;
}

template<typename F>
Base<F> Nrm2( const AbstractDistMatrix<F>& x )
{
    DEBUG_ONLY(
        CallStackEntry cse("Nrm2");
        if( x.Height() != 1 && x.Width() != 1 )
            LogicError("x must be a vector");
    )
    return FrobeniusNorm( x );
}

#define PROTO(F) \
  template Base<F> Nrm2( const Matrix<F>& x ); \
  template Base<F> Nrm2( const AbstractDistMatrix<F>& x ); 

PROTO(float)
PROTO(double)
PROTO(Complex<float>)
PROTO(Complex<double>)

} // namespace El
