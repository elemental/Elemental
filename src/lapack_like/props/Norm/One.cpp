/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {

template<typename T>
Base<T> OneNorm( const Matrix<T>& A )
{
    DEBUG_CSE
    typedef Base<T> Real;
    const Int height = A.Height();
    const Int width = A.Width();

    Real maxColSum = 0;
    for( Int j=0; j<width; ++j )
    {
        Real colSum = 0;
        for( Int i=0; i<height; ++i )
            colSum += Abs(A(i,j));
        maxColSum = Max( maxColSum, colSum );
    }
    return maxColSum;
}

template<typename T>
Base<T> HermitianOneNorm( UpperOrLower uplo, const Matrix<T>& A )
{
    DEBUG_CSE
    typedef Base<T> Real;
    const Int height = A.Height();

    if( height != A.Width() )
        RuntimeError("Hermitian matrices must be square.");

    Real maxColSum = 0;
    if( uplo == UPPER )
    {
        for( Int j=0; j<height; ++j )
        {
            Real colSum = 0;
            for( Int i=0; i<=j; ++i )
                colSum += Abs(A(i,j));
            for( Int i=j+1; i<height; ++i )
                colSum += Abs(A(j,i));
            maxColSum = Max( maxColSum, colSum );
        }
    }
    else
    {
        for( Int j=0; j<height; ++j )
        {
            Real colSum = 0;
            for( Int i=0; i<j; ++i )
                colSum += Abs(A(j,i));
            for( Int i=j; i<height; ++i )
                colSum += Abs(A(i,j));
            maxColSum = Max( maxColSum, colSum );
        }
    }
    return maxColSum;
}

template<typename T>
Base<T> SymmetricOneNorm( UpperOrLower uplo, const Matrix<T>& A )
{
    DEBUG_CSE
    return HermitianOneNorm( uplo, A );
}

template<typename T>
Base<T> OneNorm( const AbstractDistMatrix<T>& A )
{
    DEBUG_CSE
    typedef Base<T> Real;
    Real norm;
    if( A.Participating() )
    {
        // Compute the partial column sums defined by our local matrix, A[U,V]
        const Int localHeight = A.LocalHeight();
        const Int localWidth = A.LocalWidth();
        const Matrix<T>& ALoc = A.LockedMatrix();

        vector<Real> myPartialColSums( localWidth );
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            myPartialColSums[jLoc] = 0;
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                myPartialColSums[jLoc] += Abs(ALoc(iLoc,jLoc));
        }

        // Sum our partial column sums to get the column sums over A[* ,V]
        vector<Real> myColSums( localWidth );
        mpi::AllReduce
        ( myPartialColSums.data(), myColSums.data(), localWidth, A.ColComm() );

        // Find the maximum out of the column sums
        Real myMaxColSum = 0;
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            myMaxColSum = Max( myMaxColSum, myColSums[jLoc] );

        // Find the global maximum column sum by searching the row team
        norm = mpi::AllReduce( myMaxColSum, mpi::MAX, A.RowComm() );
    }
    mpi::Broadcast( norm, A.Root(), A.CrossComm() );
    return norm;
}

template<typename T>
Base<T> HermitianOneNorm( UpperOrLower uplo, const AbstractDistMatrix<T>& A )
{
    DEBUG_CSE
    typedef Base<T> Real;
    if( A.Height() != A.Width() )
        RuntimeError("Hermitian matrices must be square.");
    const Int height = A.Height();

    // For now, we take the 'easy' approach to exploiting the implicit symmetry
    // by storing all of the column sums of the triangular matrix and the 
    // row sums of the strictly triangular matrix. We can then add them.

    Real maxColSum = 0;
    if( A.Participating() )
    {
        const Int localHeight = A.LocalHeight();
        const Int localWidth = A.LocalWidth();
        const Matrix<T>& ALoc = A.LockedMatrix();

        if( uplo == UPPER )
        {
            vector<Real> myPartialUpperColSums( localWidth ),
                         myPartialStrictlyUpperRowSums( localHeight );
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const Int j = A.GlobalCol(jLoc);
                const Int numUpperRows = A.LocalRowOffset(j+1);
                myPartialUpperColSums[jLoc] = 0;
                for( Int iLoc=0; iLoc<numUpperRows; ++iLoc )
                    myPartialUpperColSums[jLoc] += Abs(ALoc(iLoc,jLoc));
            }
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
            {
                const Int i = A.GlobalRow(iLoc);
                const Int numLowerCols = A.LocalColOffset(i+1);
                myPartialStrictlyUpperRowSums[iLoc] = 0;
                for( Int jLoc=numLowerCols; jLoc<localWidth; ++jLoc )
                    myPartialStrictlyUpperRowSums[iLoc] += Abs(ALoc(iLoc,jLoc));
            }

            // Just place the sums into their appropriate places in a vector an 
            // AllReduce sum to get the results. This isn't optimal, but it 
            // should be good enough.
            vector<Real> partialColSums( height, 0 );
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const Int j = A.GlobalCol(jLoc);
                partialColSums[j] = myPartialUpperColSums[jLoc];
            }
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
            {
                const Int i = A.GlobalRow(iLoc);
                partialColSums[i] += myPartialStrictlyUpperRowSums[iLoc];
            }
            vector<Real> colSums( height );
            mpi::AllReduce
            ( partialColSums.data(), colSums.data(), height, A.DistComm() );

            // Find the maximum sum
            for( Int j=0; j<height; ++j )
                maxColSum = Max( maxColSum, colSums[j] );
        }
        else
        {
            vector<Real> myPartialLowerColSums( localWidth ),
                         myPartialStrictlyLowerRowSums( localHeight );
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const Int j = A.GlobalCol(jLoc);
                const Int numStrictlyUpperRows = A.LocalRowOffset(j);
                myPartialLowerColSums[jLoc] = 0;
                for( Int iLoc=numStrictlyUpperRows; iLoc<localHeight; ++iLoc )
                    myPartialLowerColSums[jLoc] += Abs(ALoc(iLoc,jLoc));
            }
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
            {
                const Int i = A.GlobalRow(iLoc);
                const Int numStrictlyLowerCols = A.LocalColOffset(i);
                myPartialStrictlyLowerRowSums[iLoc] = 0;
                for( Int jLoc=0; jLoc<numStrictlyLowerCols; ++jLoc )
                    myPartialStrictlyLowerRowSums[iLoc] += Abs(ALoc(iLoc,jLoc));
            }

            // Just place the sums into their appropriate places in a vector an 
            // AllReduce sum to get the results. This isn't optimal, but it 
            // should be good enough.
            vector<Real> partialColSums( height, 0 );
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const Int j = A.GlobalCol(jLoc);
                partialColSums[j] = myPartialLowerColSums[jLoc];
            }
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
            {
                const Int i = A.GlobalRow(iLoc);
                partialColSums[i] += myPartialStrictlyLowerRowSums[iLoc];
            }
            vector<Real> colSums( height );
            mpi::AllReduce
            ( partialColSums.data(), colSums.data(), height, A.DistComm() );

            // Find the maximum sum
            for( Int j=0; j<height; ++j )
                maxColSum = Max( maxColSum, colSums[j] );
        }
    }
    mpi::Broadcast( maxColSum, A.Root(), A.CrossComm() );
    return maxColSum;
}

template<typename T>
Base<T> SymmetricOneNorm( UpperOrLower uplo, const AbstractDistMatrix<T>& A )
{
    DEBUG_CSE
    return HermitianOneNorm( uplo, A );
}

template<typename T>
Base<T> OneNorm( const SparseMatrix<T>& A )
{
    DEBUG_CSE
    SparseMatrix<T> ATrans;
    Transpose( A, ATrans );
    return InfinityNorm( ATrans );
}

template<typename T>
Base<T> OneNorm( const DistSparseMatrix<T>& A )
{
    DEBUG_CSE
    DistSparseMatrix<T> ATrans(A.Comm());
    Transpose( A, ATrans );
    return InfinityNorm( ATrans );
}

template<typename T>
Base<T> HermitianTridiagOneNorm
( const Matrix<Base<T>>& d, const Matrix<T>& e )
{
    DEBUG_CSE
    typedef Base<T> Real;
    const Int n = d.Height();

    // | d(0) conj(e(0))     0       ... |
    // | e(0)    d(1)    conj(e(1))  ... |
    // |   0     e(1)       d(2)     ... |
    // |   0      0         e(2)     ... |
    // |   .      .          .        .  |

    Real maxColSum = 0;
    for( Int j=0; j<n; ++j )
    {
        Real colSum = Abs(d(j));
        if( j > 0 )
            colSum += Abs(e(j-1));
        if( j < n-1 )
            colSum += Abs(e(j));
        maxColSum = Max( maxColSum, colSum );
    }
    return maxColSum;
}

#define PROTO(T) \
  template Base<T> OneNorm( const Matrix<T>& A ); \
  template Base<T> OneNorm ( const AbstractDistMatrix<T>& A ); \
  template Base<T> HermitianOneNorm \
  ( UpperOrLower uplo, const Matrix<T>& A ); \
  template Base<T> HermitianOneNorm \
  ( UpperOrLower uplo, const AbstractDistMatrix<T>& A ); \
  template Base<T> SymmetricOneNorm \
  ( UpperOrLower uplo, const Matrix<T>& A ); \
  template Base<T> SymmetricOneNorm \
  ( UpperOrLower uplo, const AbstractDistMatrix<T>& A ); \
  template Base<T> OneNorm( const SparseMatrix<T>& A ); \
  template Base<T> OneNorm( const DistSparseMatrix<T>& A ); \
  template Base<T> HermitianTridiagOneNorm \
  ( const Matrix<Base<T>>& d, const Matrix<T>& e );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
