/*
   Copyright (c) 2009-2014, Jack Poulson
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
    {
        for( Int i=0; i<height; ++i )
        {
            const Real thisAbs = Abs(A.Get(i,j));
            maxAbs = std::max( maxAbs, thisAbs );
        }
    }
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
        {
            for( Int i=0; i<=j; ++i )
            {
                const Real thisAbs = Abs(A.Get(i,j));
                maxAbs = std::max( maxAbs, thisAbs );
            }
        }
    }
    else
    {
        for( Int j=0; j<width; ++j )
        {
            for( Int i=j; i<height; ++i )
            {
                const Real thisAbs = Abs(A.Get(i,j));
                maxAbs = std::max( maxAbs, thisAbs );
            }
        }
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
        if( uplo == UPPER )
        {
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const Int j = A.GlobalCol(jLoc);
                const Int numUpperRows = A.LocalRowOffset(j+1);
                for( Int iLoc=0; iLoc<numUpperRows; ++iLoc )
                {
                    const Real thisAbs = Abs(A.GetLocal(iLoc,jLoc));
                    localMaxAbs = std::max( localMaxAbs, thisAbs );
                }
            }
        }
        else
        {
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const Int j = A.GlobalCol(jLoc);
                const Int numStrictlyUpperRows = A.LocalRowOffset(j);
                for( Int iLoc=numStrictlyUpperRows;
                     iLoc<A.LocalHeight(); ++iLoc )
                {
                    const Real thisAbs = Abs(A.GetLocal(iLoc,jLoc));
                    localMaxAbs = std::max( localMaxAbs, thisAbs );
                }
            }
        }
        norm = mpi::AllReduce( localMaxAbs, mpi::MAX, A.DistComm() );
    }
    mpi::Broadcast( norm, A.Root(), A.CrossComm() );
    return norm;
}

template<typename T>
Base<T> SymmetricMaxNorm( UpperOrLower uplo, const AbstractDistMatrix<T>& A )
{
    DEBUG_ONLY(CallStackEntry cse("SymmetricMaxNorm"))
    return HermitianMaxNorm( uplo, A );
}

} // namespace El

#endif // ifndef EL_NORM_MAX_HPP
