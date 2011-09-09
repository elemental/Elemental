/*
   Copyright (c) 2009-2011, Jack Poulson
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

template<typename R> // represents a real number
inline R
advanced::internal::OneNorm( const Matrix<R>& A )
{
#ifndef RELEASE
    PushCallStack("advanced::internal::OneNorm");
#endif
    R maxColSum = 0;
    for( int j=0; j<A.Width(); ++j )
    {
        R colSum = 0;
        for( int i=0; i<A.Height(); ++i )
            colSum += Abs(A.Get(i,j));
        maxColSum = std::max( maxColSum, colSum );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return maxColSum;
}

#ifndef WITHOUT_COMPLEX
template<typename R> // representation of a real number
inline R
advanced::internal::OneNorm( const Matrix<std::complex<R> >& A )
{
#ifndef RELEASE
    PushCallStack("advanced::internal::OneNorm");
#endif
    R maxColSum = 0;
    for( int j=0; j<A.Width(); ++j )
    {
        R colSum = 0;
        for( int i=0; i<A.Height(); ++i )
            colSum += Abs(A.Get(i,j));
        maxColSum = std::max( maxColSum, colSum );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return maxColSum;
}
#endif // WITHOUT_COMPLEX

template<typename R> // representation of a real number
inline R
advanced::internal::OneNorm( const DistMatrix<R,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("advanced::internal::OneNorm");
#endif
    // Compute the partial column sums defined by our local matrix, A[MC,MR]
    std::vector<R> myPartialColSums(A.LocalWidth());
    for( int j=0; j<A.LocalWidth(); ++j )
    {
        myPartialColSums[j] = 0;
        for( int i=0; i<A.LocalHeight(); ++i )
            myPartialColSums[j] += Abs(A.GetLocalEntry(i,j));
    }

    // Sum our partial column sums to get the column sums over A[* ,MR]
    std::vector<R> myColSums(A.LocalWidth());
    mpi::AllReduce
    ( &myPartialColSums[0], &myColSums[0], A.LocalWidth(), mpi::SUM,
      A.Grid().MCComm() );

    // Find the maximum out of the column sums
    R myMaxColSum = 0;
    for( int j=0; j<A.LocalWidth(); ++j )
        myMaxColSum = std::max( myMaxColSum, myColSums[j] );

    // Find the global maximum column sum by searching over the MR team
    R maxColSum = 0;
    mpi::AllReduce
    ( &myMaxColSum, &maxColSum, 1, mpi::MAX, A.Grid().MRComm() );
#ifndef RELEASE
    PopCallStack();
#endif
    return maxColSum;
}

#ifndef WITHOUT_COMPLEX
template<typename R> // representation of a real number
inline R
advanced::internal::OneNorm( const DistMatrix<std::complex<R>,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("advanced::internal::OneNorm");
#endif
    // Compute the partial column sums defined by our local matrix, A[MC,MR]
    std::vector<R> myPartialColSums(A.LocalWidth());
    for( int j=0; j<A.LocalWidth(); ++j )
    {
        myPartialColSums[j] = 0;
        for( int i=0; i<A.LocalHeight(); ++i )
            myPartialColSums[j] += Abs(A.GetLocalEntry(i,j));
    }

    // Sum our partial column sums to get the column sums over A[* ,MR]
    std::vector<R> myColSums(A.LocalWidth());
    mpi::AllReduce
    ( &myPartialColSums[0], &myColSums[0], A.LocalWidth(), mpi::SUM,
      A.Grid().MCComm() );

    // Find the maximum out of the column sums
    R myMaxColSum = 0;
    for( int j=0; j<A.LocalWidth(); ++j )
        myMaxColSum = std::max( myMaxColSum, myColSums[j] );

    // Find the global maximum column sum by searching over the MR team
    R maxColSum = 0;
    mpi::AllReduce
    ( &myMaxColSum, &maxColSum, 1, mpi::MAX, A.Grid().MRComm() );
#ifndef RELEASE
    PopCallStack();
#endif
    return maxColSum;
}
#endif // WITHOUT_COMPLEX
