/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_NORM_ZERO_HPP
#define LAPACK_NORM_ZERO_HPP

// The number of nonzeros in a matrix isn't really a norm...but it's useful

namespace elem {

template<typename F>
inline int
ZeroNorm( const Matrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("ZeroNorm");
#endif
    int numNonzeros = 0;
    const int height = A.Height();
    const int width = A.Width();
    for( int j=0; j<width; ++j )
        for( int i=0; i<height; ++i )
            if( Abs(A.Get(i,j)) > 0 )
                ++numNonzeros;
    return numNonzeros;
}

template<typename F,Distribution U,Distribution V>
inline int
ZeroNorm( const DistMatrix<F,U,V>& A )
{
#ifndef RELEASE
    CallStackEntry entry("ZeroNorm");
#endif
    const int numLocalNonzeros = ZeroNorm( A.LockedMatrix() );
    mpi::Comm comm = ReduceComm<U,V>( A.Grid() );
    int numNonzeros;
    mpi::AllReduce( &numLocalNonzeros, &numNonzeros, 1, mpi::SUM, comm );
    return numNonzeros;
}

} // namespace elem

#endif // ifndef LAPACK_NORM_ZERO_HPP
