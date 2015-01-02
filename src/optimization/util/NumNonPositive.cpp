/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename Real>
Int NumNonPositive( const Matrix<Real>& A )
{
    DEBUG_ONLY(CallStackEntry cse("NumNonPositive"))
    Int numNonPos = 0;
    const Int height = A.Height();
    const Int width = A.Width();
    for( Int j=0; j<width; ++j )
        for( Int i=0; i<height; ++i )
            if( A.Get(i,j) <= Real(0) )
                ++numNonPos;
    return numNonPos;
}

template<typename Real>
Int NumNonPositive( const SparseMatrix<Real>& A )
{
    DEBUG_ONLY(CallStackEntry cse("NumNonPositive"))
    Int numNonPos = 0;
    const Int numEntries = A.NumEntries();
    for( Int k=0; k<numEntries; ++k )
        if( A.Value(k) <= Real(0) )
            ++numNonPos;
    return numNonPos;
}


template<typename Real>
Int NumNonPositive( const AbstractDistMatrix<Real>& A )
{
    DEBUG_ONLY(CallStackEntry cse("NumNonPositive"))
    Int numNonPos = 0;
    if( A.Participating() )
    {
        const Int numLocalNonPos = NumNonPositive( A.LockedMatrix() );
        numNonPos = mpi::AllReduce( numLocalNonPos, A.DistComm() );
    }
    mpi::Broadcast( numNonPos, A.Root(), A.CrossComm() );
    return numNonPos;
}

template<typename Real>
Int NumNonPositive( const DistSparseMatrix<Real>& A )
{
    DEBUG_ONLY(CallStackEntry cse("NumNonPositive"))

    Int numLocalNonPos = 0;
    const Int numLocalEntries = A.NumLocalEntries(); 
    for( Int k=0; k<numLocalEntries; ++k )
        if( A.Value(k) <= Real(0) )
            ++numLocalNonPos;

    return mpi::AllReduce( numLocalNonPos, A.Comm() );
}

template<typename Real>
Int NumNonPositive( const DistMultiVec<Real>& A )
{
    DEBUG_ONLY(CallStackEntry cse("NumNonPositive"))

    Int numLocalNonPos = 0;
    const Int localHeight = A.LocalHeight();
    const Int width = A.Width();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        for( Int j=0; j<width; ++j )
            if( A.GetLocal(iLoc,j) <= Real(0) )
                ++numLocalNonPos;

    return mpi::AllReduce( numLocalNonPos, A.Comm() );
}


#define PROTO(Real) \
  template Int NumNonPositive( const Matrix<Real>& A ); \
  template Int NumNonPositive( const SparseMatrix<Real>& A ); \
  template Int NumNonPositive( const AbstractDistMatrix<Real>& A ); \
  template Int NumNonPositive( const DistSparseMatrix<Real>& A ); \
  template Int NumNonPositive( const DistMultiVec<Real>& A );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
