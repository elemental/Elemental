/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_NORM_ENTRYWISE_HPP
#define ELEM_LAPACK_NORM_ENTRYWISE_HPP

namespace elem {

template<typename F> 
inline BASE(F)
EntrywiseNorm( const Matrix<F>& A, BASE(F) p )
{
#ifndef RELEASE
    CallStackEntry entry("EntrywiseNorm");
#endif
    // TODO: Make this more numerically stable
    typedef BASE(F) R;
    R sum = 0;
    const int width = A.Width();
    const int height = A.Height();
    for( int j=0; j<width; ++j )
        for( int i=0; i<height; ++i )
            sum += Pow( Abs(A.Get(i,j)), p );
    return Pow( sum, 1/p );
}

template<typename F>
inline BASE(F)
HermitianEntrywiseNorm
( UpperOrLower uplo, const Matrix<F>& A, BASE(F) p )
{
#ifndef RELEASE
    CallStackEntry entry("HermitianEntrywiseNorm");
#endif
    if( A.Height() != A.Width() )
        throw std::logic_error("Hermitian matrices must be square.");

    // TODO: make this more numerically stable
    typedef BASE(F) R;
    R sum = 0;
    const int height = A.Height();
    const int width = A.Width();
    if( uplo == UPPER )
    {
        for( int j=0; j<width; ++j )
        {
            for( int i=0; i<j; ++i )
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
        for( int j=0; j<width; ++j )
        {
            for( int i=j+1; i<height; ++i )
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
inline BASE(F)
SymmetricEntrywiseNorm
( UpperOrLower uplo, const Matrix<F>& A, BASE(F) p )
{
#ifndef RELEASE
    CallStackEntry entry("SymmetricEntrywiseNorm");
#endif
    return HermitianEntrywiseNorm( uplo, A, p );
}

template<typename F,Distribution U,Distribution V> 
inline BASE(F)
EntrywiseNorm( const DistMatrix<F,U,V>& A, BASE(F) p )
{
#ifndef RELEASE
    CallStackEntry entry("EntrywiseNorm");
#endif
    typedef BASE(F) R;
    R localSum = 0;
    const int localHeight = A.LocalHeight();
    const int localWidth = A.LocalWidth();
    for( int jLoc=0; jLoc<localWidth; ++jLoc )
        for( int iLoc=0; iLoc<localHeight; ++iLoc )
            localSum += Pow( Abs(A.GetLocal(iLoc,jLoc)), p ); 

    R sum;
    mpi::Comm comm = ReduceComm<U,V>( A.Grid() );
    mpi::AllReduce( &localSum, &sum, 1, mpi::SUM, comm );
    return Pow( sum, 1/p );
}

template<typename F>
inline BASE(F)
HermitianEntrywiseNorm
( UpperOrLower uplo, const DistMatrix<F>& A, BASE(F) p )
{
#ifndef RELEASE
    CallStackEntry entry("HermitianEntrywiseNorm");
#endif
    if( A.Height() != A.Width() )
        throw std::logic_error("Hermitian matrices must be square.");

    const int r = A.Grid().Height();
    const int c = A.Grid().Width();
    const int colShift = A.ColShift();
    const int rowShift = A.RowShift();

    typedef BASE(F) R;
    R localSum = 0;
    const int localWidth = A.LocalWidth();
    if( uplo == UPPER )
    {
        for( int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            int j = rowShift + jLoc*c;
            int numUpperRows = Length(j+1,colShift,r);
            for( int iLoc=0; iLoc<numUpperRows; ++iLoc )
            {
                int i = colShift + iLoc*r;
                const R term = Pow( Abs(A.GetLocal(iLoc,jLoc)), p );
                if( i ==j )
                    localSum += term;
                else
                    localSum += 2*term;
            }
        }
    }
    else
    {
        for( int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            int j = rowShift + jLoc*c;
            int numStrictlyUpperRows = Length(j,colShift,r);
            for( int iLoc=numStrictlyUpperRows;
                 iLoc<A.LocalHeight(); ++iLoc )
            {
                int i = colShift + iLoc*r;
                const R term = Pow( Abs(A.GetLocal(iLoc,jLoc)), p );
                if( i ==j )
                    localSum += term;
                else
                    localSum += 2*term;
            }
        }
    }

    R sum;
    mpi::AllReduce( &localSum, &sum, 1, mpi::SUM, A.Grid().VCComm() );
    return Pow( sum, 1/p );
}

template<typename F,Distribution U,Distribution V>
inline BASE(F)
SymmetricEntrywiseNorm
( UpperOrLower uplo, const DistMatrix<F,U,V>& A, BASE(F) p )
{
#ifndef RELEASE
    CallStackEntry entry("SymmetricEntrywiseNorm");
#endif
    return HermitianEntrywiseNorm( uplo, A, p );
}

} // namespace elem

#endif // ifndef ELEM_LAPACK_NORM_ENTRYWISE_HPP
