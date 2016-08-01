/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {
namespace pos_orth {

template<typename Real,typename>
Int NumOutside( const Matrix<Real>& A )
{
    DEBUG_CSE
    const Int height = A.Height();
    const Int width = A.Width();

    Int numNonPos = 0;
    for( Int j=0; j<width; ++j )
        for( Int i=0; i<height; ++i )
            if( A(i,j) <= Real(0) )
                ++numNonPos;
    return numNonPos;
}

template<typename Real,typename>
Int NumOutside( const SparseMatrix<Real>& A )
{
    DEBUG_CSE
    const Int numEntries = A.NumEntries();
    const Real* valBuf = A.LockedValueBuffer();

    Int numNonPos = 0;
    for( Int k=0; k<numEntries; ++k )
        if( valBuf[k] <= Real(0) )
            ++numNonPos;
    return numNonPos;
}

template<typename Real,typename>
Int NumOutside( const AbstractDistMatrix<Real>& A )
{
    DEBUG_CSE
    Int numNonPos = 0;
    if( A.Participating() )
    {
        const Int numLocalNonPos = NumOutside( A.LockedMatrix() );
        numNonPos = mpi::AllReduce( numLocalNonPos, A.DistComm() );
    }
    mpi::Broadcast( numNonPos, A.Root(), A.CrossComm() );
    return numNonPos;
}

template<typename Real,typename>
Int NumOutside( const DistSparseMatrix<Real>& A )
{
    DEBUG_CSE
    const Int numLocalEntries = A.NumLocalEntries(); 
    const Real* valBuf = A.LockedValueBuffer();

    Int numLocalNonPos = 0;
    for( Int k=0; k<numLocalEntries; ++k )
        if( valBuf[k] <= Real(0) )
            ++numLocalNonPos;

    return mpi::AllReduce( numLocalNonPos, A.Comm() );
}

template<typename Real,typename>
Int NumOutside( const DistMultiVec<Real>& A )
{
    DEBUG_CSE
    const Int localHeight = A.LocalHeight();
    const Int width = A.Width();
    const Real* ABuf = A.LockedMatrix().LockedBuffer();
    const Int ALDim = A.LockedMatrix().LDim();

    Int numLocalNonPos = 0;
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        for( Int j=0; j<width; ++j )
            if( ABuf[iLoc+j*ALDim] <= Real(0) )
                ++numLocalNonPos;

    return mpi::AllReduce( numLocalNonPos, A.Comm() );
}

#define PROTO(Real) \
  template Int NumOutside( const Matrix<Real>& A ); \
  template Int NumOutside( const SparseMatrix<Real>& A ); \
  template Int NumOutside( const AbstractDistMatrix<Real>& A ); \
  template Int NumOutside( const DistSparseMatrix<Real>& A ); \
  template Int NumOutside( const DistMultiVec<Real>& A );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace pos_orth
} // namespace El
