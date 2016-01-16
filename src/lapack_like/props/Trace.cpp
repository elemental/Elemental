/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename T>
T Trace( const Matrix<T>& A )
{
    DEBUG_ONLY(CSE cse("Trace"))
    if( A.Height() != A.Width() )
        LogicError("Cannot compute trace of nonsquare matrix");

    Matrix<T> d;
    GetDiagonal( A, d );
    T trace = 0;
    const Int n = A.Height();
    for( Int i=0; i<n; ++i )
        trace += d.Get(i,0);
    return trace;
}

template<typename T> 
T Trace( const AbstractDistMatrix<T>& A )
{
    DEBUG_ONLY(CSE cse("Trace"))
    if( A.Height() != A.Width() )
        LogicError("Cannot compute trace of nonsquare matrix");

    T trace = 0;
    if( A.Participating() )
    {
        T localTrace = 0;
        for( Int jLoc=0; jLoc<A.LocalWidth(); ++jLoc )
        {
            const Int j = A.GlobalCol(jLoc);
            if( A.IsLocalRow(j) )
            {
                const Int iLoc = A.LocalRow(j);
                localTrace += A.GetLocal(iLoc,jLoc); 
            }
        }
        trace = mpi::AllReduce( localTrace, A.DistComm() ); 
    }
    mpi::Broadcast( trace, A.Root(), A.CrossComm() );
    return trace;
}

#define PROTO(T) \
  template T Trace( const Matrix<T>& A ); \
  template T Trace( const AbstractDistMatrix<T>& A );

#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

} // namespace El
