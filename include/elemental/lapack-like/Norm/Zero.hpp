/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_NORM_ZERO_HPP
#define ELEM_LAPACK_NORM_ZERO_HPP

// The number of nonzeros in a matrix isn't really a norm...but it's useful

namespace elem {

template<typename F>
inline Int
ZeroNorm( const Matrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("ZeroNorm");
#endif
    Int numNonzeros = 0;
    const Int height = A.Height();
    const Int width = A.Width();
    for( Int j=0; j<width; ++j )
        for( Int i=0; i<height; ++i )
            if( Abs(A.Get(i,j)) > 0 )
                ++numNonzeros;
    return numNonzeros;
}

template<typename F,Distribution U,Distribution V>
inline Int
ZeroNorm( const DistMatrix<F,U,V>& A )
{
#ifndef RELEASE
    CallStackEntry entry("ZeroNorm");
#endif
    Int numNonzeros;
    if( A.Participating() )
    {
        const Int numLocalNonzeros = ZeroNorm( A.LockedMatrix() );
        numNonzeros = mpi::AllReduce( numLocalNonzeros, A.DistComm() );
    }
    mpi::Broadcast( numNonzeros, A.Root(), A.CrossComm() );
    return numNonzeros;
}

} // namespace elem

#endif // ifndef ELEM_LAPACK_NORM_ZERO_HPP
