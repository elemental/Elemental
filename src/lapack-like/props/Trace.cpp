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
F Trace( const Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("Trace"))
    if( A.Height() != A.Width() )
        LogicError("Cannot compute trace of nonsquare matrix");

    Matrix<F> d;
    A.GetDiagonal( d );
    F trace = 0;
    const Int n = A.Height();
    for( Int i=0; i<n; ++i )
        trace += d.Get(i,0);
    return trace;
}

template<typename F> 
F Trace( const DistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("Trace"))
    if( A.Height() != A.Width() )
        LogicError("Cannot compute trace of nonsquare matrix");
    const Grid& g = A.Grid();

    DistMatrix<F,MD,STAR> d(g);
    A.GetDiagonal( d );
    F localTrace = 0;
    if( d.Participating() )
    {
        const Int nLocalDiag = d.LocalHeight();
        for( Int iLoc=0; iLoc<nLocalDiag; ++iLoc )
            localTrace += d.GetLocal(iLoc,0);
    }
    return mpi::AllReduce( localTrace, g.VCComm() );
}

#define PROTO(F) \
  template F Trace( const Matrix<F>& A ); \
  template F Trace( const DistMatrix<F>& A );

PROTO(float)
PROTO(double)
PROTO(Complex<float>)
PROTO(Complex<double>)

} // namespace El
