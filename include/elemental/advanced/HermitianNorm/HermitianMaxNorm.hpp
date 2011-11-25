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

template<typename R> // representation of a real number
inline R
elemental::advanced::internal::HermitianMaxNorm
( UpperOrLower uplo, const Matrix<R>& A )
{
#ifndef RELEASE
    PushCallStack("advanced::internal::HermitianMaxNorm");
#endif
    if( A.Height() != A.Width() )
        throw std::logic_error("Hermitian matrices must be square.");

    R maxAbs = 0;
    if( uplo == UPPER )
    {
        for( int j=0; j<A.Width(); ++j )
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
        for( int j=0; j<A.Width(); ++j )
        {
            for( int i=j; i<A.Height(); ++i )
            {
                const R thisAbs = Abs(A.Get(i,j));
                maxAbs = std::max( maxAbs, thisAbs );
            }
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return maxAbs;
}

#ifndef WITHOUT_COMPLEX
template<typename R> // representation of a real number
inline R
elemental::advanced::internal::HermitianMaxNorm
( UpperOrLower uplo, const Matrix<std::complex<R> >& A )
{
#ifndef RELEASE
    PushCallStack("advanced::internal::HermitianMaxNorm");
#endif
    if( A.Height() != A.Width() )
        throw std::logic_error("Hermitian matrices must be square.");

    R maxAbs = 0;
    if( uplo == UPPER )
    {
        for( int j=0; j<A.Width(); ++j )
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
        for( int j=0; j<A.Width(); ++j )
        {
            for( int i=j; i<A.Height(); ++i )
            {
                const R thisAbs = Abs(A.Get(i,j));
                maxAbs = std::max( maxAbs, thisAbs );
            }
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return maxAbs;
}
#endif // WITHOUT_COMPLEX

template<typename R> // representation of a real number
inline R
elemental::advanced::internal::HermitianMaxNorm
( UpperOrLower uplo, const DistMatrix<R,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("advanced::internal::HermitianMaxNorm");
#endif
    if( A.Height() != A.Width() )
        throw std::logic_error("Hermitian matrices must be square.");

    const int r = A.Grid().Height();
    const int c = A.Grid().Width();
    const int colShift = A.ColShift();
    const int rowShift = A.RowShift();

    R localMaxAbs = 0;
    if( uplo == UPPER )
    {
        for( int jLocal=0; jLocal<A.LocalWidth(); ++jLocal )
        {
            int j = rowShift + jLocal*c;
            int numUpperRows = LocalLength(j+1,colShift,r);
            for( int iLocal=0; iLocal<numUpperRows; ++iLocal )
            {
                const R thisAbs = Abs(A.GetLocalEntry(iLocal,jLocal));
                localMaxAbs = std::max( localMaxAbs, thisAbs );
            }
        }
    }
    else
    {
        for( int jLocal=0; jLocal<A.LocalWidth(); ++jLocal )
        {
            int j = rowShift + jLocal*c;
            int numStrictlyUpperRows = LocalLength(j,colShift,r);
            for( int iLocal=numStrictlyUpperRows; 
                 iLocal<A.LocalHeight(); ++iLocal )
            {
                const R thisAbs = Abs(A.GetLocalEntry(iLocal,jLocal));
                localMaxAbs = std::max( localMaxAbs, thisAbs );
            }
        }
    }

    R maxAbs;
    mpi::AllReduce
    ( &localMaxAbs, &maxAbs, 1, mpi::MAX, A.Grid().VCComm() );
#ifndef RELEASE
    PopCallStack();
#endif
    return maxAbs;
}

#ifndef WITHOUT_COMPLEX
template<typename R> // representation of a real number
inline R
elemental::advanced::internal::HermitianMaxNorm
( UpperOrLower uplo, const DistMatrix<std::complex<R>,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("advanced::internal::HermitianMaxNorm");
#endif
    if( A.Height() != A.Width() )
        throw std::logic_error("Hermitian matrices must be square.");

    const int r = A.Grid().Height();
    const int c = A.Grid().Width();
    const int colShift = A.ColShift();
    const int rowShift = A.RowShift();

    R localMaxAbs = 0;
    if( uplo == UPPER )
    {
        for( int jLocal=0; jLocal<A.LocalWidth(); ++jLocal )
        {
            int j = rowShift + jLocal*c;
            int numUpperRows = LocalLength(j+1,colShift,r);
            for( int iLocal=0; iLocal<numUpperRows; ++iLocal )
            {
                const R thisAbs = Abs(A.GetLocalEntry(iLocal,jLocal));
                localMaxAbs = std::max( localMaxAbs, thisAbs );
            }
        }
    }
    else
    {
        for( int jLocal=0; jLocal<A.LocalWidth(); ++jLocal )
        {
            int j = rowShift + jLocal*c;
            int numStrictlyUpperRows = LocalLength(j,colShift,r);
            for( int iLocal=numStrictlyUpperRows; 
                 iLocal<A.LocalHeight(); ++iLocal )
            {
                const R thisAbs = Abs(A.GetLocalEntry(iLocal,jLocal));
                localMaxAbs = std::max( localMaxAbs, thisAbs );
            }
        }
    }

    R maxAbs;
    mpi::AllReduce
    ( &localMaxAbs, &maxAbs, 1, mpi::MAX, A.Grid().VCComm() );
#ifndef RELEASE
    PopCallStack();
#endif
    return maxAbs;
}
#endif // WITHOUT_COMPLEX
