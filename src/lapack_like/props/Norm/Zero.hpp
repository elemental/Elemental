/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_NORM_ZERO_HPP
#define EL_NORM_ZERO_HPP

// The number of nonzeros in a matrix isn't really a norm...
// but the terminology is common

namespace El {

template<typename T>
Int ZeroNorm( const Matrix<T>& A, Base<T> tol )
{
    DEBUG_ONLY(CallStackEntry cse("ZeroNorm"))
    Int numNonzeros = 0;
    const Int height = A.Height();
    const Int width = A.Width();
    for( Int j=0; j<width; ++j )
        for( Int i=0; i<height; ++i )
            if( Abs(A.Get(i,j)) > tol )
                ++numNonzeros;
    return numNonzeros;
}

template<typename T>
Int ZeroNorm( const SparseMatrix<T>& A, Base<T> tol )
{
    DEBUG_ONLY(CallStackEntry cse("ZeroNorm"))
    Int numNonzeros = 0;
    const Int numEntries = A.NumEntries();
    for( Int k=0; k<numEntries; ++k )
        if( Abs(A.Value(k)) > tol )
            ++numNonzeros;
    return numNonzeros;
}

template<typename T>
Int ZeroNorm( const AbstractDistMatrix<T>& A, Base<T> tol )
{
    DEBUG_ONLY(CallStackEntry cse("ZeroNorm"))
    Int numNonzeros;
    if( A.Participating() )
    {
        const Int numLocalNonzeros = ZeroNorm( A.LockedMatrix(), tol );
        numNonzeros = mpi::AllReduce( numLocalNonzeros, A.DistComm() );
    }
    mpi::Broadcast( numNonzeros, A.Root(), A.CrossComm() );
    return numNonzeros;
}

template<typename T>
Int ZeroNorm( const DistSparseMatrix<T>& A, Base<T> tol )
{
    DEBUG_ONLY(CallStackEntry cse("ZeroNorm"))
    Int numNonzeros = 0;
    const Int numLocalEntries = A.NumLocalEntries();
    for( Int k=0; k<numLocalEntries; ++k )
        if( Abs(A.Value(k)) > tol )
            ++numNonzeros;
    return mpi::AllReduce( numNonzeros, A.Comm() );
}

} // namespace El

#endif // ifndef EL_NORM_ZERO_HPP
