/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_NORM_ONE_HPP
#define LAPACK_NORM_ONE_HPP

namespace elem {

template<typename F>
inline typename Base<F>::type
OneNorm( const Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("OneNorm");
#endif
    typedef typename Base<F>::type R;
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
#ifndef RELEASE
    PopCallStack();
#endif
    return maxColSum;
}

template<typename F>
inline typename Base<F>::type
HermitianOneNorm( UpperOrLower uplo, const Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("HermitianOneNorm");
#endif
    typedef typename Base<F>::type R;
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
#ifndef RELEASE
    PopCallStack();
#endif
    return maxColSum;
}

template<typename F>
inline typename Base<F>::type
SymmetricOneNorm( UpperOrLower uplo, const Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("SymmetricOneNorm");
#endif
    typedef typename Base<F>::type R;
    const R norm = HermitianOneNorm( uplo, A );
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

template<typename F,Distribution U,Distribution V>
inline typename Base<F>::type
OneNorm( const DistMatrix<F,U,V>& A )
{
#ifndef RELEASE
    PushCallStack("OneNorm");
#endif
    // Compute the partial column sums defined by our local matrix, A[U,V]
    typedef typename Base<F>::type R;
    const int localHeight = A.LocalHeight();
    const int localWidth = A.LocalWidth();
    std::vector<R> myPartialColSums( localWidth );
    for( int jLocal=0; jLocal<localWidth; ++jLocal )
    {
        myPartialColSums[jLocal] = 0;
        for( int iLocal=0; iLocal<localHeight; ++iLocal )
            myPartialColSums[jLocal] += Abs(A.GetLocal(iLocal,jLocal));
    }

    // Sum our partial column sums to get the column sums over A[* ,V]
    std::vector<R> myColSums( localWidth );
    mpi::Comm colComm = ReduceColComm<U,V>( A.Grid() );
    mpi::AllReduce
    ( &myPartialColSums[0], &myColSums[0], localWidth, mpi::SUM, colComm );

    // Find the maximum out of the column sums
    R myMaxColSum = 0;
    for( int jLocal=0; jLocal<localWidth; ++jLocal )
        myMaxColSum = std::max( myMaxColSum, myColSums[jLocal] );

    // Find the global maximum column sum by searching over the MR team
    R maxColSum = 0;
    mpi::Comm rowComm = ReduceRowComm<U,V>( A.Grid() );
    mpi::AllReduce( &myMaxColSum, &maxColSum, 1, mpi::MAX, rowComm );
#ifndef RELEASE
    PopCallStack();
#endif
    return maxColSum;
}

template<typename F>
inline typename Base<F>::type
HermitianOneNorm( UpperOrLower uplo, const DistMatrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("HermitianOneNorm");
#endif
    typedef typename Base<F>::type R;

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
        for( int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            int j = rowShift + jLocal*c;
            int numUpperRows = Length(j+1,colShift,r);
            myPartialUpperColSums[jLocal] = 0;
            for( int iLocal=0; iLocal<numUpperRows; ++iLocal )
                myPartialUpperColSums[jLocal] +=
                    Abs(A.GetLocal(iLocal,jLocal));
        }
        for( int iLocal=0; iLocal<localHeight; ++iLocal )
        {
            int i = colShift + iLocal*r;
            int numLowerCols = Length(i+1,rowShift,c);
            myPartialStrictlyUpperRowSums[iLocal] = 0;
            for( int jLocal=numLowerCols; jLocal<localWidth; ++jLocal )
                myPartialStrictlyUpperRowSums[iLocal] +=
                    Abs(A.GetLocal(iLocal,jLocal));
        }

        // Just place the sums into their appropriate places in a vector an 
        // AllReduce sum to get the results. This isn't optimal, but it should
        // be good enough.
        std::vector<R> partialColSums( height, 0 );
        for( int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            int j = rowShift + jLocal*c;
            partialColSums[j] = myPartialUpperColSums[jLocal];
        }
        for( int iLocal=0; iLocal<localHeight; ++iLocal )
        {
            int i = colShift + iLocal*r;
            partialColSums[i] += myPartialStrictlyUpperRowSums[iLocal];
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
        for( int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            int j = rowShift + jLocal*c;
            int numStrictlyUpperRows = Length(j,colShift,r);
            myPartialLowerColSums[jLocal] = 0;
            for( int iLocal=numStrictlyUpperRows; iLocal<localHeight; ++iLocal )
                myPartialLowerColSums[jLocal] +=
                    Abs(A.GetLocal(iLocal,jLocal));
        }
        for( int iLocal=0; iLocal<localHeight; ++iLocal )
        {
            int i = colShift + iLocal*r;
            int numStrictlyLowerCols = Length(i,rowShift,c);
            myPartialStrictlyLowerRowSums[iLocal] = 0;
            for( int jLocal=0; jLocal<numStrictlyLowerCols; ++jLocal )
                myPartialStrictlyLowerRowSums[iLocal] +=
                    Abs(A.GetLocal(iLocal,jLocal));
        }

        // Just place the sums into their appropriate places in a vector an 
        // AllReduce sum to get the results. This isn't optimal, but it should
        // be good enough.
        std::vector<R> partialColSums( height, 0 );
        for( int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            int j = rowShift + jLocal*c;
            partialColSums[j] = myPartialLowerColSums[jLocal];
        }
        for( int iLocal=0; iLocal<localHeight; ++iLocal )
        {
            int i = colShift + iLocal*r;
            partialColSums[i] += myPartialStrictlyLowerRowSums[iLocal];
        }
        std::vector<R> colSums( height );
        mpi::AllReduce
        ( &partialColSums[0], &colSums[0], height, mpi::SUM,
          A.Grid().VCComm() );

        // Find the maximum sum
        for( int j=0; j<height; ++j )
            maxColSum = std::max( maxColSum, colSums[j] );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return maxColSum;
}

template<typename F>
inline typename Base<F>::type
SymmetricOneNorm( UpperOrLower uplo, const DistMatrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("SymmetricOneNorm");
#endif
    typedef typename Base<F>::type R;
    const R norm = HermitianOneNorm( uplo, A );
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

} // namespace elem

#endif // ifndef LAPACK_NORM_ONE_HPP
