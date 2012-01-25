/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

    - Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

    - Neither the name of the owner nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/

namespace elem {

template<typename F> 
inline typename Base<F>::type
internal::HermitianOneNorm( UpperOrLower uplo, const Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("internal::HermitianOneNorm");
#endif
    typedef typename Base<F>::type R;

    if( A.Height() != A.Width() )
        throw std::runtime_error("Hermitian matrices must be square.");
    R maxColSum = 0;
    if( uplo == UPPER )
    {
        for( int j=0; j<A.Width(); ++j )
        {
            R colSum = 0;
            for( int i=0; i<=j; ++i )
                colSum += Abs(A.Get(i,j));
            for( int i=j+1; i<A.Height(); ++i )
                colSum += Abs(A.Get(j,i));
            maxColSum = std::max( maxColSum, colSum );
        }
    }
    else
    {
        for( int j=0; j<A.Width(); ++j )
        {
            R colSum = 0;
            for( int i=0; i<j; ++i )
                colSum += Abs(A.Get(j,i));
            for( int i=j; i<A.Height(); ++i )
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
internal::HermitianOneNorm( UpperOrLower uplo, const DistMatrix<F,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("internal::HermitianOneNorm");
#endif
    typedef typename Base<F>::type R;

    if( A.Height() != A.Width() )
        throw std::runtime_error("Hermitian matrices must be square.");

    // For now, we take the 'easy' approach to exploiting the implicit symmetry
    // by storing all of the column sums of the triangular matrix and the 
    // row sums of the strictly triangular matrix. We can then add them.

    int r = A.Grid().Height();
    int c = A.Grid().Width();
    int rowShift = A.RowShift();
    int colShift = A.ColShift();

    R maxColSum = 0;
    if( uplo == UPPER )
    {
        std::vector<R> myPartialUpperColSums(A.LocalWidth());
        std::vector<R> myPartialStrictlyUpperRowSums(A.LocalHeight());
        for( int jLocal=0; jLocal<A.LocalWidth(); ++jLocal )
        {
            int j = rowShift + jLocal*c;
            int numUpperRows = LocalLength(j+1,colShift,r);
            myPartialUpperColSums[jLocal] = 0;
            for( int iLocal=0; iLocal<numUpperRows; ++iLocal )
                myPartialUpperColSums[jLocal] += 
                    Abs(A.GetLocalEntry(iLocal,jLocal));
        } 
        for( int iLocal=0; iLocal<A.LocalHeight(); ++iLocal )
        {
            int i = colShift + iLocal*r;
            int numLowerCols = LocalLength(i+1,rowShift,c);
            myPartialStrictlyUpperRowSums[iLocal] = 0;
            for( int jLocal=numLowerCols; jLocal<A.LocalWidth(); ++jLocal )
                myPartialStrictlyUpperRowSums[iLocal] += 
                    Abs(A.GetLocalEntry(iLocal,jLocal));
        }

        // Just place the sums into their appropriate places in a vector an 
        // AllReduce sum to get the results. This isn't optimal, but it should
        // be good enough.
        std::vector<R> partialColSums(A.Width(),0);
        for( int jLocal=0; jLocal<A.LocalWidth(); ++jLocal )
        {
            int j = rowShift + jLocal*c;
            partialColSums[j] = myPartialUpperColSums[jLocal];
        }
        for( int iLocal=0; iLocal<A.LocalHeight(); ++iLocal )
        {
            int i = colShift + iLocal*r;
            partialColSums[i] += myPartialStrictlyUpperRowSums[iLocal];
        }
        std::vector<R> colSums(A.Width());
        mpi::AllReduce
        ( &partialColSums[0], &colSums[0], A.Width(), mpi::SUM, 
          A.Grid().VCComm() );

        // Find the maximum sum
        for( int j=0; j<A.Width(); ++j )
            maxColSum = std::max( maxColSum, colSums[j] );
    }
    else
    {
        std::vector<R> myPartialLowerColSums(A.LocalWidth());
        std::vector<R> myPartialStrictlyLowerRowSums(A.LocalHeight());
        for( int jLocal=0; jLocal<A.LocalWidth(); ++jLocal )
        {
            int j = rowShift + jLocal*c;
            int numStrictlyUpperRows = LocalLength(j,colShift,r);
            myPartialLowerColSums[jLocal] = 0;
            for( int iLocal=numStrictlyUpperRows; 
                 iLocal<A.LocalHeight(); ++iLocal )
                myPartialLowerColSums[jLocal] += 
                    Abs(A.GetLocalEntry(iLocal,jLocal));
        } 
        for( int iLocal=0; iLocal<A.LocalHeight(); ++iLocal )
        {
            int i = colShift + iLocal*r;
            int numStrictlyLowerCols = LocalLength(i,rowShift,c);
            myPartialStrictlyLowerRowSums[iLocal] = 0;
            for( int jLocal=0; jLocal<numStrictlyLowerCols; ++jLocal )
                myPartialStrictlyLowerRowSums[iLocal] += 
                    Abs(A.GetLocalEntry(iLocal,jLocal));
        }

        // Just place the sums into their appropriate places in a vector an 
        // AllReduce sum to get the results. This isn't optimal, but it should
        // be good enough.
        std::vector<R> partialColSums(A.Width(),0);
        for( int jLocal=0; jLocal<A.LocalWidth(); ++jLocal )
        {
            int j = rowShift + jLocal*c;
            partialColSums[j] = myPartialLowerColSums[jLocal];
        }
        for( int iLocal=0; iLocal<A.LocalHeight(); ++iLocal )
        {
            int i = colShift + iLocal*r;
            partialColSums[i] += myPartialStrictlyLowerRowSums[iLocal];
        }
        std::vector<R> colSums(A.Width());
        mpi::AllReduce
        ( &partialColSums[0], &colSums[0], A.Width(), mpi::SUM, 
          A.Grid().VCComm() );

        // Find the maximum sum
        for( int j=0; j<A.Width(); ++j )
            maxColSum = std::max( maxColSum, colSums[j] );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return maxColSum;
}

} // namespace elem
