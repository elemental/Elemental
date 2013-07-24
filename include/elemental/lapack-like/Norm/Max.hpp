/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_NORM_MAX_HPP
#define ELEM_LAPACK_NORM_MAX_HPP

namespace elem {

template<typename F> 
inline BASE(F)
MaxNorm( const Matrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("MaxNorm");
#endif
    typedef BASE(F) R;
    R maxAbs = 0;
    const int height = A.Height();
    const int width = A.Width();
    for( int j=0; j<width; ++j )
    {
        for( int i=0; i<height; ++i )
        {
            const R thisAbs = Abs(A.Get(i,j));
            maxAbs = std::max( maxAbs, thisAbs );
        }
    }
    return maxAbs;
}

template<typename F>
inline BASE(F)
HermitianMaxNorm( UpperOrLower uplo, const Matrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("HermitianMaxNorm");
#endif
    if( A.Height() != A.Width() )
        throw std::logic_error("Hermitian matrices must be square.");

    typedef BASE(F) R;
    R maxAbs = 0;
    const int height = A.Height();
    const int width = A.Width();
    if( uplo == UPPER )
    {
        for( int j=0; j<width; ++j )
        {
            for( int i=0; i<=j; ++i )
            {
                const R thisAbs = Abs(A.Get(i,j));
                maxAbs = std::max( maxAbs, thisAbs );
            }
        }
    }
    else
    {
        for( int j=0; j<width; ++j )
        {
            for( int i=j; i<height; ++i )
            {
                const R thisAbs = Abs(A.Get(i,j));
                maxAbs = std::max( maxAbs, thisAbs );
            }
        }
    }
    return maxAbs;
}

template<typename F>
inline BASE(F)
SymmetricMaxNorm( UpperOrLower uplo, const Matrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("SymmetricMaxNorm");
#endif
    return HermitianMaxNorm( uplo, A );
}

template<typename F,Distribution U,Distribution V>
inline BASE(F)
MaxNorm( const DistMatrix<F,U,V>& A )
{
#ifndef RELEASE
    CallStackEntry entry("MaxNorm");
#endif
    typedef BASE(F) R;
    R localMaxAbs = 0;
    const int localHeight = A.LocalHeight();
    const int localWidth = A.LocalWidth();
    for( int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        for( int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const R thisAbs = Abs(A.GetLocal(iLoc,jLoc));
            localMaxAbs = std::max( localMaxAbs, thisAbs );
        }
    }

    R maxAbs;
    mpi::Comm reduceComm = ReduceComm<U,V>( A.Grid() );
    mpi::AllReduce( &localMaxAbs, &maxAbs, 1, mpi::MAX, reduceComm );
    return maxAbs;
}

template<typename F>
inline BASE(F)
HermitianMaxNorm( UpperOrLower uplo, const DistMatrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("HermitianMaxNorm");
#endif
    if( A.Height() != A.Width() )
        throw std::logic_error("Hermitian matrices must be square.");

    const int r = A.Grid().Height();
    const int c = A.Grid().Width();
    const int colShift = A.ColShift();
    const int rowShift = A.RowShift();

    typedef BASE(F) R;
    R localMaxAbs = 0;
    const int localWidth = A.LocalWidth();
    if( uplo == UPPER )
    {
        for( int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            int j = rowShift + jLoc*c;
            int numUpperRows = Length(j+1,colShift,r);
            for( int iLoc=0; iLoc<numUpperRows; ++iLoc )
            {
                const R thisAbs = Abs(A.GetLocal(iLoc,jLoc));
                localMaxAbs = std::max( localMaxAbs, thisAbs );
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
                const R thisAbs = Abs(A.GetLocal(iLoc,jLoc));
                localMaxAbs = std::max( localMaxAbs, thisAbs );
            }
        }
    }

    R maxAbs;
    mpi::AllReduce( &localMaxAbs, &maxAbs, 1, mpi::MAX, A.Grid().VCComm() );
    return maxAbs;
}

template<typename F>
inline BASE(F)
SymmetricMaxNorm( UpperOrLower uplo, const DistMatrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("SymmetricMaxNorm");
#endif
    return HermitianMaxNorm( uplo, A );
}

} // namespace elem

#endif // ifndef ELEM_LAPACK_NORM_MAX_HPP
