/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_NORM_MAX_HPP
#define EL_NORM_MAX_HPP

namespace El {

template<typename T> 
Base<T> MaxNorm( const Matrix<T>& A )
{
    DEBUG_ONLY(CallStackEntry cse("MaxNorm"))
    typedef Base<T> Real;
    Real maxAbs = 0;
    const Int height = A.Height();
    const Int width = A.Width();
    for( Int j=0; j<width; ++j )
        for( Int i=0; i<height; ++i )
            maxAbs = Max( maxAbs, Abs(A.Get(i,j)) );
    return maxAbs;
}

template<typename T> 
Base<T> MaxNorm( const SparseMatrix<T>& A )
{
    DEBUG_ONLY(CallStackEntry cse("MaxNorm"))
    typedef Base<T> Real;
    Real maxAbs = 0;
    const Int numEntries = A.NumEntries();
    for( Int k=0; k<numEntries; ++k )
        maxAbs = Max( maxAbs, Abs(A.Value(k)) );
    return maxAbs;
}

template<typename T>
Base<T> HermitianMaxNorm( UpperOrLower uplo, const Matrix<T>& A )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianMaxNorm"))
    if( A.Height() != A.Width() )
        LogicError("Hermitian matrices must be square.");

    typedef Base<T> Real;
    Real maxAbs = 0;
    const Int height = A.Height();
    const Int width = A.Width();
    if( uplo == UPPER )
    {
        for( Int j=0; j<width; ++j )
            for( Int i=0; i<=j; ++i )
                maxAbs = Max( maxAbs, Abs(A.Get(i,j)) );
    }
    else
    {
        for( Int j=0; j<width; ++j )
            for( Int i=j; i<height; ++i )
                maxAbs = Max( maxAbs, Abs(A.Get(i,j)) );
    }
    return maxAbs;
}

template<typename T> 
Base<T> HermitianMaxNorm( UpperOrLower uplo, const SparseMatrix<T>& A )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianMaxNorm"))
    if( A.Height() != A.Width() )
        LogicError("Hermitian matrices must be square.");

    typedef Base<T> Real;
    Real maxAbs = 0;
    const Int numEntries = A.NumEntries();
    for( Int k=0; k<numEntries; ++k )
    {
        const Int i = A.Row(k);
        const Int j = A.Col(k);
        if( (uplo==UPPER && i<=j) || (uplo==LOWER && i>=j) )
            maxAbs = Max( maxAbs, Abs(A.Value(k)) );
    }
    return maxAbs;
}

template<typename T>
Base<T> SymmetricMaxNorm( UpperOrLower uplo, const Matrix<T>& A )
{
    DEBUG_ONLY(CallStackEntry cse("SymmetricMaxNorm"))
    return HermitianMaxNorm( uplo, A );
}

template<typename T>
Base<T> SymmetricMaxNorm( UpperOrLower uplo, const SparseMatrix<T>& A )
{
    DEBUG_ONLY(CallStackEntry cse("SymmetricMaxNorm"))
    return HermitianMaxNorm( uplo, A );
}

template<typename T>
Base<T> MaxNorm( const AbstractDistMatrix<T>& A )
{
    DEBUG_ONLY(CallStackEntry cse("MaxNorm"))
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
    DEBUG_ONLY(CallStackEntry cse("MaxNorm"))

    Base<T> localNorm = 0;
    const Int numLocalEntries = A.NumLocalEntries();
    for( Int k=0; k<numLocalEntries; ++k )
        localNorm = Max(localNorm,Abs(A.Value(k)));

    return mpi::AllReduce( localNorm, mpi::MAX, A.Comm() );
}

template<typename T>
Base<T> MaxNorm( const DistMultiVec<T>& A )
{
    DEBUG_ONLY(CallStackEntry cse("MaxNorm"))
    
    Base<T> localNorm = 0;
    const Int localHeight = A.LocalHeight();
    const Int width = A.Width();
    for( Int j=0; j<width; ++j )
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
            localNorm = Max(localNorm,Abs(A.GetLocal(iLoc,j)));
    
    return mpi::AllReduce( localNorm, mpi::MAX, A.Comm() );
}

template<typename T>
Base<T> HermitianMaxNorm( UpperOrLower uplo, const AbstractDistMatrix<T>& A )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianMaxNorm"))
    if( A.Height() != A.Width() )
        LogicError("Hermitian matrices must be square.");

    typedef Base<T> Real;
    Real norm;
    if( A.Participating() )
    {
        Real localMaxAbs = 0;
        const Int localWidth = A.LocalWidth();
        const Int localHeight = A.LocalHeight();
        if( uplo == UPPER )
        {
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const Int j = A.GlobalCol(jLoc);
                const Int numUpperRows = A.LocalRowOffset(j+1);
                for( Int iLoc=0; iLoc<numUpperRows; ++iLoc )
                    localMaxAbs = Max(localMaxAbs,Abs(A.GetLocal(iLoc,jLoc)));
            }
        }
        else
        {
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const Int j = A.GlobalCol(jLoc);
                const Int numStrictlyUpperRows = A.LocalRowOffset(j);
                for( Int iLoc=numStrictlyUpperRows; iLoc<localHeight; ++iLoc )
                    localMaxAbs = Max(localMaxAbs,Abs(A.GetLocal(iLoc,jLoc)));
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
    DEBUG_ONLY(CallStackEntry cse("HermitianMaxNorm"))
    if( A.Height() != A.Width() )
        LogicError("Hermitian matrices must be square.");

    Base<T> localNorm = 0;
    const Int numLocalEntries = A.NumLocalEntries();
    for( Int k=0; k<numLocalEntries; ++k )
    {
        const Int i = A.Row(k);
        const Int j = A.Col(k);
        if( (uplo==UPPER && i<=j) || (uplo==LOWER && i>=j) )
            localNorm = Max(localNorm,Abs(A.Value(k)));
    }

    return mpi::AllReduce( localNorm, mpi::MAX, A.Comm() );
}

template<typename T>
Base<T> SymmetricMaxNorm( UpperOrLower uplo, const AbstractDistMatrix<T>& A )
{
    DEBUG_ONLY(CallStackEntry cse("SymmetricMaxNorm"))
    return HermitianMaxNorm( uplo, A );
}

template<typename T>
Base<T> SymmetricMaxNorm( UpperOrLower uplo, const DistSparseMatrix<T>& A )
{
    DEBUG_ONLY(CallStackEntry cse("SymmetricMaxNorm"))
    return HermitianMaxNorm( uplo, A );
}

} // namespace El

#endif // ifndef EL_NORM_MAX_HPP
