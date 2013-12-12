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
    DEBUG_ONLY(CallStackEntry cse("OneNorm"))
    typedef Base<F> R;
    R maxColSum = 0;
    const Int height = A.Height();
    const Int width = A.Width();
    for( Int j=0; j<width; ++j )
    {
        R colSum = 0;
        for( Int i=0; i<height; ++i )
            colSum += Abs(A.Get(i,j));
        maxColSum = std::max( maxColSum, colSum );
    }
    return maxColSum;
}

template<typename F>
inline BASE(F)
HermitianOneNorm( UpperOrLower uplo, const Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianOneNorm"))
    typedef Base<F> R;
    if( A.Height() != A.Width() )
        RuntimeError("Hermitian matrices must be square.");
    R maxColSum = 0;
    const Int height = A.Height();
    if( uplo == UPPER )
    {
        for( Int j=0; j<height; ++j )
        {
            R colSum = 0;
            for( Int i=0; i<=j; ++i )
                colSum += Abs(A.Get(i,j));
            for( Int i=j+1; i<height; ++i )
                colSum += Abs(A.Get(j,i));
            maxColSum = std::max( maxColSum, colSum );
        }
    }
    else
    {
        for( Int j=0; j<height; ++j )
        {
            R colSum = 0;
            for( Int i=0; i<j; ++i )
                colSum += Abs(A.Get(j,i));
            for( Int i=j; i<height; ++i )
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
    DEBUG_ONLY(CallStackEntry cse("SymmetricOneNorm"))
    return HermitianOneNorm( uplo, A );
}

template<typename F,Distribution U,Distribution V>
inline BASE(F)
OneNorm( const DistMatrix<F,U,V>& A )
{
    DEBUG_ONLY(CallStackEntry cse("OneNorm"))
    typedef Base<F> Real;
    Real norm;
    if( A.Participating() )
    {
        // Compute the partial column sums defined by our local matrix, A[U,V]
        const Int localHeight = A.LocalHeight();
        const Int localWidth = A.LocalWidth();
        std::vector<Real> myPartialColSums( localWidth );
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            myPartialColSums[jLoc] = 0;
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                myPartialColSums[jLoc] += Abs(A.GetLocal(iLoc,jLoc));
        }

        // Sum our partial column sums to get the column sums over A[* ,V]
        std::vector<Real> myColSums( localWidth );
        mpi::AllReduce
        ( myPartialColSums.data(), myColSums.data(), localWidth, A.ColComm() );

        // Find the maximum out of the column sums
        Real myMaxColSum = 0;
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            myMaxColSum = std::max( myMaxColSum, myColSums[jLoc] );

        // Find the global maximum column sum by searching the row team
        norm = mpi::AllReduce( myMaxColSum, A.RowComm() );
    }
    mpi::Broadcast( norm, A.Root(), A.CrossComm() );
    return norm;
}

template<typename F>
inline BASE(F)
HermitianOneNorm( UpperOrLower uplo, const DistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianOneNorm"))
    typedef Base<F> R;
    if( A.Height() != A.Width() )
        RuntimeError("Hermitian matrices must be square.");
    const Int height = A.Height();

    // For now, we take the 'easy' approach to exploiting the implicit symmetry
    // by storing all of the column sums of the triangular matrix and the 
    // row sums of the strictly triangular matrix. We can then add them.

    Int r = A.Grid().Height();
    Int c = A.Grid().Width();
    Int rowShift = A.RowShift();
    Int colShift = A.ColShift();

    R maxColSum = 0;
    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();
    if( uplo == UPPER )
    {
        std::vector<R> myPartialUpperColSums( localWidth ),
                       myPartialStrictlyUpperRowSums( localHeight );
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            Int j = rowShift + jLoc*c;
            Int numUpperRows = Length(j+1,colShift,r);
            myPartialUpperColSums[jLoc] = 0;
            for( Int iLoc=0; iLoc<numUpperRows; ++iLoc )
                myPartialUpperColSums[jLoc] +=
                    Abs(A.GetLocal(iLoc,jLoc));
        }
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            Int i = colShift + iLoc*r;
            Int numLowerCols = Length(i+1,rowShift,c);
            myPartialStrictlyUpperRowSums[iLoc] = 0;
            for( Int jLoc=numLowerCols; jLoc<localWidth; ++jLoc )
                myPartialStrictlyUpperRowSums[iLoc] +=
                    Abs(A.GetLocal(iLoc,jLoc));
        }

        // Just place the sums into their appropriate places in a vector an 
        // AllReduce sum to get the results. This isn't optimal, but it should
        // be good enough.
        std::vector<R> partialColSums( height, 0 );
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            Int j = rowShift + jLoc*c;
            partialColSums[j] = myPartialUpperColSums[jLoc];
        }
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            Int i = colShift + iLoc*r;
            partialColSums[i] += myPartialStrictlyUpperRowSums[iLoc];
        }
        std::vector<R> colSums( height );
        mpi::AllReduce
        ( partialColSums.data(), colSums.data(), height, A.Grid().VCComm() );

        // Find the maximum sum
        for( Int j=0; j<height; ++j )
            maxColSum = std::max( maxColSum, colSums[j] );
    }
    else
    {
        std::vector<R> myPartialLowerColSums( localWidth ),
                       myPartialStrictlyLowerRowSums( localHeight );
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            Int j = rowShift + jLoc*c;
            Int numStrictlyUpperRows = Length(j,colShift,r);
            myPartialLowerColSums[jLoc] = 0;
            for( Int iLoc=numStrictlyUpperRows; iLoc<localHeight; ++iLoc )
                myPartialLowerColSums[jLoc] += Abs(A.GetLocal(iLoc,jLoc));
        }
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            Int i = colShift + iLoc*r;
            Int numStrictlyLowerCols = Length(i,rowShift,c);
            myPartialStrictlyLowerRowSums[iLoc] = 0;
            for( Int jLoc=0; jLoc<numStrictlyLowerCols; ++jLoc )
                myPartialStrictlyLowerRowSums[iLoc] +=
                    Abs(A.GetLocal(iLoc,jLoc));
        }

        // Just place the sums into their appropriate places in a vector an 
        // AllReduce sum to get the results. This isn't optimal, but it should
        // be good enough.
        std::vector<R> partialColSums( height, 0 );
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            Int j = rowShift + jLoc*c;
            partialColSums[j] = myPartialLowerColSums[jLoc];
        }
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            Int i = colShift + iLoc*r;
            partialColSums[i] += myPartialStrictlyLowerRowSums[iLoc];
        }
        std::vector<R> colSums( height );
        mpi::AllReduce
        ( partialColSums.data(), colSums.data(), height, A.Grid().VCComm() );

        // Find the maximum sum
        for( Int j=0; j<height; ++j )
            maxColSum = std::max( maxColSum, colSums[j] );
    }
    return maxColSum;
}

template<typename F>
inline BASE(F)
SymmetricOneNorm( UpperOrLower uplo, const DistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("SymmetricOneNorm"))
    return HermitianOneNorm( uplo, A );
}

} // namespace elem

#endif // ifndef ELEM_LAPACK_NORM_ONE_HPP
