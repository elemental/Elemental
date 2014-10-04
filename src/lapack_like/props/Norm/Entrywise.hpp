/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_NORM_ENTRYWISE_HPP
#define EL_NORM_ENTRYWISE_HPP

namespace El {

template<typename F> 
Base<F> EntrywiseNorm( const Matrix<F>& A, Base<F> p )
{
    DEBUG_ONLY(CallStackEntry cse("EntrywiseNorm"))
    // TODO: Make this more numerically stable
    typedef Base<F> R;
    R sum = 0;
    const Int width = A.Width();
    const Int height = A.Height();
    for( Int j=0; j<width; ++j )
        for( Int i=0; i<height; ++i )
            sum += Pow( Abs(A.Get(i,j)), p );
    return Pow( sum, 1/p );
}

template<typename F>
Base<F> HermitianEntrywiseNorm
( UpperOrLower uplo, const Matrix<F>& A, Base<F> p )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianEntrywiseNorm"))
    if( A.Height() != A.Width() )
        LogicError("Hermitian matrices must be square.");

    // TODO: make this more numerically stable
    typedef Base<F> R;
    R sum = 0;
    const Int height = A.Height();
    const Int width = A.Width();
    if( uplo == UPPER )
    {
        for( Int j=0; j<width; ++j )
        {
            for( Int i=0; i<j; ++i )
            {
                const R term = Pow( Abs(A.Get(i,j)), p );
                if( i ==j )
                    sum += term;
                else
                    sum += 2*term;
            }
        }
    }
    else
    {
        for( Int j=0; j<width; ++j )
        {
            for( Int i=j+1; i<height; ++i )
            {
                const R term = Pow( Abs(A.Get(i,j)), p );
                if( i ==j )
                    sum += term;
                else
                    sum += 2*term;
            }
        }
    }
    return Pow( sum, 1/p );
}

template<typename F>
Base<F> SymmetricEntrywiseNorm
( UpperOrLower uplo, const Matrix<F>& A, Base<F> p )
{
    DEBUG_ONLY(CallStackEntry cse("SymmetricEntrywiseNorm"))
    return HermitianEntrywiseNorm( uplo, A, p );
}

template<typename F> 
Base<F> EntrywiseNorm( const AbstractDistMatrix<F>& A, Base<F> p )
{
    DEBUG_ONLY(CallStackEntry cse("EntrywiseNorm"))
    typedef Base<F> Real;
    Real norm;
    if( A.Participating() )
    {
        Real localSum = 0;
        const Int localHeight = A.LocalHeight();
        const Int localWidth = A.LocalWidth();
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                localSum += Pow( Abs(A.GetLocal(iLoc,jLoc)), p ); 
        const Real sum = mpi::AllReduce( localSum, A.DistComm() );
        norm = Pow( sum, 1/p );
    }
    mpi::Broadcast( norm, A.Root(), A.CrossComm() );
    return norm;
}

template<typename F>
Base<F> HermitianEntrywiseNorm
( UpperOrLower uplo, const AbstractDistMatrix<F>& A, Base<F> p )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianEntrywiseNorm"))
    if( A.Height() != A.Width() )
        LogicError("Hermitian matrices must be square.");

    typedef Base<F> Real;
    Real sum;
    if( A.Participating() )
    {
        Real localSum = 0;
        const Int localWidth = A.LocalWidth();
        if( uplo == UPPER )
        {
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const Int j = A.GlobalCol(jLoc);
                const Int numUpperRows = A.LocalRowOffset(j+1);
                for( Int iLoc=0; iLoc<numUpperRows; ++iLoc )
                {
                    const Int i = A.GlobalRow(iLoc);
                    const Real term = Pow( Abs(A.GetLocal(iLoc,jLoc)), p );
                    if( i ==j )
                        localSum += term;
                    else
                        localSum += 2*term;
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
                    const Int i = A.GlobalRow(iLoc);
                    const Real term = Pow( Abs(A.GetLocal(iLoc,jLoc)), p );
                    if( i ==j )
                        localSum += term;
                    else
                        localSum += 2*term;
                }
            }
        }
        sum = mpi::AllReduce( localSum, A.DistComm() );
    }
    mpi::Broadcast( sum, A.Root(), A.CrossComm() );
    return Pow( sum, 1/p );
}

template<typename F>
Base<F> SymmetricEntrywiseNorm
( UpperOrLower uplo, const AbstractDistMatrix<F>& A, Base<F> p )
{
    DEBUG_ONLY(CallStackEntry cse("SymmetricEntrywiseNorm"))
    return HermitianEntrywiseNorm( uplo, A, p );
}

} // namespace El

#endif // ifndef EL_NORM_ENTRYWISE_HPP
