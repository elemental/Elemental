/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_NORM_ONE_HPP
#define ELEM_LAPACK_NORM_ONE_HPP

namespace elem {

template<typename F>
inline BASE(F)
OneNorm( const Matrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("OneNorm");
#endif
    typedef BASE(F) R;
    R maxColSum = 0;
    const int height = A.Height();
    const int width = A.Width();
    for( int j=0; j<width; ++j )
    {
        R colSum = 0;
        for( int i=0; i<height; ++i )
            colSum += Abs(A.Get(i,j));
        maxColSum = std::max( maxColSum, colSum );
    }
    return maxColSum;
}

template<typename F>
inline BASE(F)
HermitianOneNorm( UpperOrLower uplo, const Matrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("HermitianOneNorm");
#endif
    typedef BASE(F) R;
    if( A.Height() != A.Width() )
        throw std::runtime_error("Hermitian matrices must be square.");
    R maxColSum = 0;
    const int height = A.Height();
    if( uplo == UPPER )
    {
        for( int j=0; j<height; ++j )
        {
            R colSum = 0;
            for( int i=0; i<=j; ++i )
                colSum += Abs(A.Get(i,j));
            for( int i=j+1; i<height; ++i )
                colSum += Abs(A.Get(j,i));
            maxColSum = std::max( maxColSum, colSum );
        }
    }
    else
    {
        for( int j=0; j<height; ++j )
        {
            R colSum = 0;
            for( int i=0; i<j; ++i )
                colSum += Abs(A.Get(j,i));
            for( int i=j; i<height; ++i )
                colSum += Abs(A.Get(i,j));
            maxColSum = std::max( maxColSum, colSum );
        }
    }
    return maxColSum;
}

template<typename F>
inline BASE(F)
SymmetricOneNorm( UpperOrLower uplo, const Matrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("SymmetricOneNorm");
#endif
    return HermitianOneNorm( uplo, A );
}

template<typename F,Distribution U,Distribution V>
inline BASE(F)
OneNorm( const DistMatrix<F,U,V>& A )
{
#ifndef RELEASE
    CallStackEntry entry("OneNorm");
#endif
    // Compute the partial column sums defined by our local matrix, A[U,V]
    typedef BASE(F) R;
    const int localHeight = A.LocalHeight();
    const int localWidth = A.LocalWidth();
    std::vector<R> myPartialColSums( localWidth );
    for( int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        myPartialColSums[jLoc] = 0;
        for( int iLoc=0; iLoc<localHeight; ++iLoc )
            myPartialColSums[jLoc] += Abs(A.GetLocal(iLoc,jLoc));
    }

    // Sum our partial column sums to get the column sums over A[* ,V]
    std::vector<R> myColSums( localWidth );
    mpi::Comm colComm = ReduceColComm<U,V>( A.Grid() );
    mpi::AllReduce
    ( &myPartialColSums[0], &myColSums[0], localWidth, mpi::SUM, colComm );

    // Find the maximum out of the column sums
    R myMaxColSum = 0;
    for( int jLoc=0; jLoc<localWidth; ++jLoc )
        myMaxColSum = std::max( myMaxColSum, myColSums[jLoc] );

    // Find the global maximum column sum by searching over the MR team
    R maxColSum = 0;
    mpi::Comm rowComm = ReduceRowComm<U,V>( A.Grid() );
    mpi::AllReduce( &myMaxColSum, &maxColSum, 1, mpi::MAX, rowComm );
    return maxColSum;
}

template<typename F>
inline BASE(F)
HermitianOneNorm( UpperOrLower uplo, const DistMatrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("HermitianOneNorm");
#endif
    typedef BASE(F) R;

    if( A.Height() != A.Width() )
        throw std::runtime_error("Hermitian matrices must be square.");
    const int height = A.Height();

    // For now, we take the 'easy' approach to exploiting the implicit symmetry
    // by storing all of the column sums of the triangular matrix and the 
    // row sums of the strictly triangular matrix. We can then add them.

    int r = A.Grid().Height();
    int c = A.Grid().Width();
    int rowShift = A.RowShift();
    int colShift = A.ColShift();

    R maxColSum = 0;
    const int localHeight = A.LocalHeight();
    const int localWidth = A.LocalWidth();
    if( uplo == UPPER )
    {
        std::vector<R> myPartialUpperColSums( localWidth );
        std::vector<R> myPartialStrictlyUpperRowSums( localHeight );
        for( int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            int j = rowShift + jLoc*c;
            int numUpperRows = Length(j+1,colShift,r);
            myPartialUpperColSums[jLoc] = 0;
            for( int iLoc=0; iLoc<numUpperRows; ++iLoc )
                myPartialUpperColSums[jLoc] +=
                    Abs(A.GetLocal(iLoc,jLoc));
        }
        for( int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            int i = colShift + iLoc*r;
            int numLowerCols = Length(i+1,rowShift,c);
            myPartialStrictlyUpperRowSums[iLoc] = 0;
            for( int jLoc=numLowerCols; jLoc<localWidth; ++jLoc )
                myPartialStrictlyUpperRowSums[iLoc] +=
                    Abs(A.GetLocal(iLoc,jLoc));
        }

        // Just place the sums into their appropriate places in a vector an 
        // AllReduce sum to get the results. This isn't optimal, but it should
        // be good enough.
        std::vector<R> partialColSums( height, 0 );
        for( int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            int j = rowShift + jLoc*c;
            partialColSums[j] = myPartialUpperColSums[jLoc];
        }
        for( int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            int i = colShift + iLoc*r;
            partialColSums[i] += myPartialStrictlyUpperRowSums[iLoc];
        }
        std::vector<R> colSums( height );
        mpi::AllReduce
        ( &partialColSums[0], &colSums[0], height, mpi::SUM,
          A.Grid().VCComm() );

        // Find the maximum sum
        for( int j=0; j<height; ++j )
            maxColSum = std::max( maxColSum, colSums[j] );
    }
    else
    {
        std::vector<R> myPartialLowerColSums( localWidth );
        std::vector<R> myPartialStrictlyLowerRowSums( localHeight );
        for( int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            int j = rowShift + jLoc*c;
            int numStrictlyUpperRows = Length(j,colShift,r);
            myPartialLowerColSums[jLoc] = 0;
            for( int iLoc=numStrictlyUpperRows; iLoc<localHeight; ++iLoc )
                myPartialLowerColSums[jLoc] +=
                    Abs(A.GetLocal(iLoc,jLoc));
        }
        for( int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            int i = colShift + iLoc*r;
            int numStrictlyLowerCols = Length(i,rowShift,c);
            myPartialStrictlyLowerRowSums[iLoc] = 0;
            for( int jLoc=0; jLoc<numStrictlyLowerCols; ++jLoc )
                myPartialStrictlyLowerRowSums[iLoc] +=
                    Abs(A.GetLocal(iLoc,jLoc));
        }

        // Just place the sums into their appropriate places in a vector an 
        // AllReduce sum to get the results. This isn't optimal, but it should
        // be good enough.
        std::vector<R> partialColSums( height, 0 );
        for( int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            int j = rowShift + jLoc*c;
            partialColSums[j] = myPartialLowerColSums[jLoc];
        }
        for( int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            int i = colShift + iLoc*r;
            partialColSums[i] += myPartialStrictlyLowerRowSums[iLoc];
        }
        std::vector<R> colSums( height );
        mpi::AllReduce
        ( &partialColSums[0], &colSums[0], height, mpi::SUM,
          A.Grid().VCComm() );

        // Find the maximum sum
        for( int j=0; j<height; ++j )
            maxColSum = std::max( maxColSum, colSums[j] );
    }
    return maxColSum;
}

template<typename F>
inline BASE(F)
SymmetricOneNorm( UpperOrLower uplo, const DistMatrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("SymmetricOneNorm");
#endif
    return HermitianOneNorm( uplo, A );
}

} // namespace elem

#endif // ifndef ELEM_LAPACK_NORM_ONE_HPP
