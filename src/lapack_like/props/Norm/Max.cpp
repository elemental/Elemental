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
Base<T> MaxNorm( const Matrix<T>& A )
{
    DEBUG_ONLY(CSE cse("MaxNorm"))
    typedef Base<T> Real;
    const Int height = A.Height();
    const Int width = A.Width();
    const T* ABuf = A.LockedBuffer();
    const Int ALDim = A.LDim();

    Real maxAbs = 0;
    for( Int j=0; j<width; ++j )
        for( Int i=0; i<height; ++i )
            maxAbs = Max( maxAbs, Abs(ABuf[i+j*ALDim]) );

    return maxAbs;
}

template<typename T> 
Base<T> MaxNorm( const SparseMatrix<T>& A )
{
    DEBUG_ONLY(CSE cse("MaxNorm"))
    typedef Base<T> Real;
    const Int numEntries = A.NumEntries();
    const T* AValBuf = A.LockedValueBuffer();

    Real maxAbs = 0;
    for( Int k=0; k<numEntries; ++k )
        maxAbs = Max( maxAbs, Abs(AValBuf[k]) );

    return maxAbs;
}

template<typename T>
Base<T> HermitianMaxNorm( UpperOrLower uplo, const Matrix<T>& A )
{
    DEBUG_ONLY(CSE cse("HermitianMaxNorm"))
    if( A.Height() != A.Width() )
        LogicError("Hermitian matrices must be square.");

    typedef Base<T> Real;
    const Int height = A.Height();
    const Int width = A.Width();
    const T* ABuf = A.LockedBuffer();
    const Int ALDim = A.LDim();

    Real maxAbs = 0;
    if( uplo == UPPER )
    {
        for( Int j=0; j<width; ++j )
            for( Int i=0; i<=j; ++i )
                maxAbs = Max( maxAbs, Abs(ABuf[i+j*ALDim]) );
    }
    else
    {
        for( Int j=0; j<width; ++j )
            for( Int i=j; i<height; ++i )
                maxAbs = Max( maxAbs, Abs(ABuf[i+j*ALDim]) );
    }
    return maxAbs;
}

template<typename T> 
Base<T> HermitianMaxNorm( UpperOrLower uplo, const SparseMatrix<T>& A )
{
    DEBUG_ONLY(CSE cse("HermitianMaxNorm"))
    if( A.Height() != A.Width() )
        LogicError("Hermitian matrices must be square.");

    typedef Base<T> Real;
    const Int numEntries = A.NumEntries();
    const Int* ARowBuf = A.LockedSourceBuffer();
    const Int* AColBuf = A.LockedTargetBuffer();
    const T* AValBuf = A.LockedValueBuffer();

    Real maxAbs = 0;
    for( Int k=0; k<numEntries; ++k )
    {
        const Int i = ARowBuf[k];
        const Int j = AColBuf[k];
        if( (uplo==UPPER && i<=j) || (uplo==LOWER && i>=j) )
            maxAbs = Max( maxAbs, Abs(AValBuf[k]) );
    }
    return maxAbs;
}

template<typename T>
Base<T> SymmetricMaxNorm( UpperOrLower uplo, const Matrix<T>& A )
{
    DEBUG_ONLY(CSE cse("SymmetricMaxNorm"))
    return HermitianMaxNorm( uplo, A );
}

template<typename T>
Base<T> SymmetricMaxNorm( UpperOrLower uplo, const SparseMatrix<T>& A )
{
    DEBUG_ONLY(CSE cse("SymmetricMaxNorm"))
    return HermitianMaxNorm( uplo, A );
}

template<typename T>
Base<T> MaxNorm( const AbstractDistMatrix<T>& A )
{
    DEBUG_ONLY(CSE cse("MaxNorm"))
    Base<T> norm=0;
    if( A.Participating() )
    {
        Base<T> localMaxAbs = MaxNorm( A.LockedMatrix() );
        norm = mpi::AllReduce( localMaxAbs, mpi::MAX, A.DistComm() );
    }
    mpi::Broadcast( norm, A.Root(), A.CrossComm() );
    return norm;
}

template<typename T>
Base<T> MaxNorm( const DistSparseMatrix<T>& A )
{
    DEBUG_ONLY(CSE cse("MaxNorm"))
    const Int numLocalEntries = A.NumLocalEntries();
    const T* AValBuf = A.LockedValueBuffer();

    Base<T> localNorm = 0;
    for( Int k=0; k<numLocalEntries; ++k )
        localNorm = Max( localNorm, Abs(AValBuf[k]) );

    return mpi::AllReduce( localNorm, mpi::MAX, A.Comm() );
}

template<typename T>
Base<T> MaxNorm( const DistMultiVec<T>& A )
{
    DEBUG_ONLY(CSE cse("MaxNorm"))
    Base<T> localNorm = MaxNorm( A.LockedMatrix() );
    return mpi::AllReduce( localNorm, mpi::MAX, A.Comm() );
}

template<typename T>
Base<T> HermitianMaxNorm( UpperOrLower uplo, const AbstractDistMatrix<T>& A )
{
    DEBUG_ONLY(CSE cse("HermitianMaxNorm"))
    if( A.Height() != A.Width() )
        LogicError("Hermitian matrices must be square.");
    typedef Base<T> Real;

    Real norm;
    if( A.Participating() )
    {
        const Int localWidth = A.LocalWidth();
        const Int localHeight = A.LocalHeight();
        const T* ABuf = A.LockedBuffer();
        const Int ALDim = A.LDim();

        Real localMaxAbs = 0;
        if( uplo == UPPER )
        {
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const Int j = A.GlobalCol(jLoc);
                const Int numUpperRows = A.LocalRowOffset(j+1);
                for( Int iLoc=0; iLoc<numUpperRows; ++iLoc )
                    localMaxAbs = Max(localMaxAbs,Abs(ABuf[iLoc+jLoc*ALDim]));
            }
        }
        else
        {
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const Int j = A.GlobalCol(jLoc);
                const Int numStrictlyUpperRows = A.LocalRowOffset(j);
                for( Int iLoc=numStrictlyUpperRows; iLoc<localHeight; ++iLoc )
                    localMaxAbs = Max(localMaxAbs,Abs(ABuf[iLoc+jLoc*ALDim]));
            }
        }
        norm = mpi::AllReduce( localMaxAbs, mpi::MAX, A.DistComm() );
    }
    mpi::Broadcast( norm, A.Root(), A.CrossComm() );
    return norm;
}

template<typename T>
Base<T> HermitianMaxNorm( UpperOrLower uplo, const DistSparseMatrix<T>& A )
{
    DEBUG_ONLY(CSE cse("HermitianMaxNorm"))
    if( A.Height() != A.Width() )
        LogicError("Hermitian matrices must be square.");
    const Int numLocalEntries = A.NumLocalEntries();
    const T* AValBuf = A.LockedValueBuffer();
    const Int* ARowBuf = A.LockedSourceBuffer();
    const Int* AColBuf = A.LockedTargetBuffer();

    Base<T> localNorm = 0;
    for( Int k=0; k<numLocalEntries; ++k )
    {
        const Int i = ARowBuf[k];
        const Int j = AColBuf[k];
        if( (uplo==UPPER && i<=j) || (uplo==LOWER && i>=j) )
            localNorm = Max( localNorm, Abs(AValBuf[k]) );
    }

    return mpi::AllReduce( localNorm, mpi::MAX, A.Comm() );
}

template<typename T>
Base<T> SymmetricMaxNorm( UpperOrLower uplo, const AbstractDistMatrix<T>& A )
{
    DEBUG_ONLY(CSE cse("SymmetricMaxNorm"))
    return HermitianMaxNorm( uplo, A );
}

template<typename T>
Base<T> SymmetricMaxNorm( UpperOrLower uplo, const DistSparseMatrix<T>& A )
{
    DEBUG_ONLY(CSE cse("SymmetricMaxNorm"))
    return HermitianMaxNorm( uplo, A );
}

#define PROTO(T) \
  template Base<T> MaxNorm( const Matrix<T>& A ); \
  template Base<T> MaxNorm ( const AbstractDistMatrix<T>& A ); \
  template Base<T> MaxNorm( const SparseMatrix<T>& A ); \
  template Base<T> MaxNorm( const DistSparseMatrix<T>& A ); \
  template Base<T> MaxNorm( const DistMultiVec<T>& A ); \
  template Base<T> HermitianMaxNorm \
  ( UpperOrLower uplo, const Matrix<T>& A ); \
  template Base<T> HermitianMaxNorm \
  ( UpperOrLower uplo, const AbstractDistMatrix<T>& A ); \
  template Base<T> HermitianMaxNorm \
  ( UpperOrLower uplo, const SparseMatrix<T>& A ); \
  template Base<T> HermitianMaxNorm \
  ( UpperOrLower uplo, const DistSparseMatrix<T>& A ); \
  template Base<T> SymmetricMaxNorm \
  ( UpperOrLower uplo, const Matrix<T>& A ); \
  template Base<T> SymmetricMaxNorm \
  ( UpperOrLower uplo, const AbstractDistMatrix<T>& A ); \
  template Base<T> SymmetricMaxNorm \
  ( UpperOrLower uplo, const SparseMatrix<T>& A ); \
  template Base<T> SymmetricMaxNorm \
  ( UpperOrLower uplo, const DistSparseMatrix<T>& A );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

} // namespace El
