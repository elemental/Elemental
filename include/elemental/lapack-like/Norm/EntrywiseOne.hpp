/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_NORM_ENTRYWISEONE_HPP
#define LAPACK_NORM_ENTRYWISEONE_HPP

namespace elem {

template<typename F> 
inline BASE(F)
EntrywiseOneNorm( const Matrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("EntrywiseOneNorm");
#endif
    typedef BASE(F) R;
    R norm = 0;
    const int width = A.Width();
    const int height = A.Height();
    for( int j=0; j<width; ++j )
        for( int i=0; i<height; ++i )
            norm += Abs(A.Get(i,j));
    return norm;
}

template<typename F>
inline BASE(F)
HermitianEntrywiseOneNorm( UpperOrLower uplo, const Matrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("HermitianEntrywiseOneNorm");
#endif
    if( A.Height() != A.Width() )
        throw std::logic_error("Hermitian matrices must be square.");

    typedef BASE(F) R;
    R norm = 0;
    const int height = A.Height();
    const int width = A.Width();
    if( uplo == UPPER )
    {
        for( int j=0; j<width; ++j )
        {
            for( int i=0; i<j; ++i )
            {
                const R alpha = Abs(A.Get(i,j));
                if( i ==j )
                    norm += alpha;
                else
                    norm += 2*alpha;
            }
        }
    }
    else
    {
        for( int j=0; j<width; ++j )
        {
            for( int i=j+1; i<height; ++i )
            {
                const R alpha = Abs(A.Get(i,j));
                if( i ==j )
                    norm += alpha;
                else
                    norm += 2*alpha;
            }
        }
    }
    return norm;
}

template<typename F>
inline BASE(F)
SymmetricEntrywiseOneNorm( UpperOrLower uplo, const Matrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("SymmetricEntrywiseOneNorm");
#endif
    return HermitianEntrywiseOneNorm( uplo, A );
}

template<typename F,Distribution U,Distribution V> 
inline BASE(F)
EntrywiseOneNorm( const DistMatrix<F,U,V>& A )
{
#ifndef RELEASE
    CallStackEntry entry("EntrywiseOneNorm");
#endif
    typedef BASE(F) R;
    R localSum = 0;
    const int localHeight = A.LocalHeight();
    const int localWidth = A.LocalWidth();
    for( int jLoc=0; jLoc<localWidth; ++jLoc )
        for( int iLoc=0; iLoc<localHeight; ++iLoc )
            localSum += Abs(A.GetLocal(iLoc,jLoc)); 

    R norm;
    mpi::Comm comm = ReduceComm<U,V>( A.Grid() );
    mpi::AllReduce( &localSum, &norm, 1, mpi::SUM, comm );
    return norm;
}

template<typename F>
inline BASE(F)
HermitianEntrywiseOneNorm( UpperOrLower uplo, const DistMatrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("HermitianEntrywiseOneNorm");
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
                const R alpha = Abs(A.GetLocal(iLoc,jLoc));
                if( i ==j )
                    localSum += alpha;
                else
                    localSum += 2*alpha;
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
                const R alpha = Abs(A.GetLocal(iLoc,jLoc));
                if( i ==j )
                    localSum += alpha;
                else
                    localSum += 2*alpha;
            }
        }
    }

    R norm;
    mpi::AllReduce( &localSum, &norm, 1, mpi::SUM, A.Grid().VCComm() );
    return norm;
}

template<typename F,Distribution U,Distribution V>
inline BASE(F)
SymmetricEntrywiseOneNorm( UpperOrLower uplo, const DistMatrix<F,U,V>& A )
{
#ifndef RELEASE
    CallStackEntry entry("SymmetricEntrywiseOneNorm");
#endif
    return HermitianEntrywiseOneNorm( uplo, A );
}

} // namespace elem

#endif // ifndef LAPACK_NORM_ENTRYWISEONE_HPP
