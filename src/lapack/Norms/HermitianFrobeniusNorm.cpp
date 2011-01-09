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
#include "elemental/lapack.hpp"
using namespace elemental;

template<typename R> // representation of a real number
R
lapack::HermitianFrobeniusNorm
( Shape shape, const Matrix<R>& A )
{
#ifndef RELEASE
    PushCallStack("lapack::HermitianFrobeniusNorm");
#endif
    if( A.Height() != A.Width() )
        throw std::logic_error("Hermitian matrices must be square.");

    R normSquared = 0;
    if( shape == Upper )
    {
        for( int j=0; j<A.Width(); ++j )
        {
            for( int i=0; i<j; ++j )
            {
                R alpha = A.Get(i,j);
                normSquared += 2*alpha*alpha;
            }
            R alpha = A.Get(j,j);
            normSquared += alpha*alpha;
        }
    }
    else
    {
        for( int j=0; j<A.Width(); ++j )
        {
            for( int i=j+1; i<A.Height(); ++i )
            {
                R alpha = A.Get(i,j);
                normSquared += 2*alpha*alpha;
            }
            R alpha = A.Get(j,j);
            normSquared += alpha*alpha;
        }
    }

    R norm = sqrt(normSquared);
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

#ifndef WITHOUT_COMPLEX
template<typename R> // representation of a real number
R
lapack::HermitianFrobeniusNorm
( Shape shape, const Matrix< std::complex<R> >& A )
{
#ifndef RELEASE
    PushCallStack("lapack::HermitianFrobeniusNorm");
#endif
    if( A.Height() != A.Width() )
        throw std::logic_error("Hermitian matrices must be square.");

    // The std::norm function is a field norm rather than a vector norm.

    R normSquared = 0;
    if( shape == Upper )
    {
        for( int j=0; j<A.Width(); ++j )
        {
            for( int i=0; i<j; ++j )
            {
                std::complex<R> alpha = A.Get(i,j);
                normSquared += 2*norm(alpha);
            }
            std::complex<R> alpha = A.Get(j,j);
            normSquared += norm(alpha);
        }
    }
    else
    {
        for( int j=0; j<A.Width(); ++j )
        {
            for( int i=j+1; i<A.Height(); ++i )
            {
                std::complex<R> alpha = A.Get(i,j);
                normSquared += 2*norm(alpha);
            }
            std::complex<R> alpha = A.Get(j,j);
            normSquared += norm(alpha);
        }
    }

    R norm = sqrt(normSquared);
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}
#endif

template<typename R> // representation of a real number
R
lapack::HermitianFrobeniusNorm
( Shape shape, const DistMatrix<R,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("lapack::HermitianFrobeniusNorm");
#endif
    if( A.Height() != A.Width() )
        throw std::logic_error("Hermitian matrices must be square.");

    const int r = A.Grid().Height();
    const int c = A.Grid().Width();
    const int colShift = A.ColShift();
    const int rowShift = A.RowShift();

    R localNormSquared = 0;
    if( shape == Upper )
    {
        for( int jLocal=0; jLocal<A.LocalWidth(); ++jLocal )
        {
            int j = rowShift + jLocal*c;
            int numUpperRows = utilities::LocalLength(j+1,colShift,r);
            for( int iLocal=0; iLocal<numUpperRows; ++iLocal )
            {
                int i = colShift + iLocal*r;
                R alpha = A.GetLocalEntry(iLocal,jLocal);
                if( i != j )
                    localNormSquared += 2*alpha*alpha;
                else
                    localNormSquared += alpha*alpha;
            }
        }
    }
    else
    {
        for( int jLocal=0; jLocal<A.LocalWidth(); ++jLocal )
        {
            int j = rowShift + jLocal*c;
            int numStrictlyUpperRows = utilities::LocalLength(j,colShift,r);
            for( int iLocal=numStrictlyUpperRows; 
                 iLocal<A.LocalHeight(); ++iLocal )
            {
                int i = colShift + iLocal*r;
                R alpha = A.GetLocalEntry(iLocal,jLocal);
                if( i != j )
                    localNormSquared += 2*alpha*alpha;
                else
                    localNormSquared += alpha*alpha;
            }
        }
    }

    // Sum the local contributions
    R normSquared;
    import::mpi::AllReduce
    ( &localNormSquared, &normSquared, 1, MPI_SUM, A.Grid().VCComm() );

    R norm = sqrt(normSquared);
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

#ifndef WITHOUT_COMPLEX
template<typename R> // representation of a real number
R
lapack::HermitianFrobeniusNorm
( Shape shape, const DistMatrix<std::complex<R>,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("lapack::HermitianFrobeniusNorm");
#endif
    if( A.Height() != A.Width() )
        throw std::logic_error("Hermitian matrices must be square.");

    const int r = A.Grid().Height();
    const int c = A.Grid().Width();
    const int colShift = A.ColShift();
    const int rowShift = A.RowShift();

    // The std::norm function is a field norm rather than a vector norm.

    R localNormSquared = 0;
    if( shape == Upper )
    {
        for( int jLocal=0; jLocal<A.LocalWidth(); ++jLocal )
        {
            int j = rowShift + jLocal*c;
            int numUpperRows = utilities::LocalLength(j+1,colShift,r);
            for( int iLocal=0; iLocal<numUpperRows; ++iLocal )
            {
                int i = colShift + iLocal*r;
                std::complex<R> alpha = A.GetLocalEntry(iLocal,jLocal);
                if( i != j )
                    localNormSquared += 2*norm(alpha);
                else
                    localNormSquared += norm(alpha);
            }
        }
    }
    else
    {
        for( int jLocal=0; jLocal<A.LocalWidth(); ++jLocal )
        {
            int j = rowShift + jLocal*c;
            int numStrictlyUpperRows = utilities::LocalLength(j,colShift,r);
            for( int iLocal=numStrictlyUpperRows; 
                 iLocal<A.LocalHeight(); ++iLocal )
            {
                int i = colShift + iLocal*r;
                std::complex<R> alpha = A.GetLocalEntry(iLocal,jLocal);
                if( i != j )
                    localNormSquared += 2*norm(alpha);
                else
                    localNormSquared += norm(alpha);
            }
        }
    }

    // Sum the local contributions
    R normSquared;
    import::mpi::AllReduce
    ( &localNormSquared, &normSquared, 1, MPI_SUM, A.Grid().VCComm() );

    R norm = sqrt(normSquared);
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}
#endif

template float elemental::lapack::HermitianFrobeniusNorm
( Shape shape, const Matrix<float>& A );
template double elemental::lapack::HermitianFrobeniusNorm
( Shape shape, const Matrix<double>& A );
#ifndef WITHOUT_COMPLEX
template float elemental::lapack::HermitianFrobeniusNorm
( Shape shape, const Matrix< std::complex<float> >& A );
template double elemental::lapack::HermitianFrobeniusNorm
( Shape shape, const Matrix< std::complex<double> >& A );
#endif

template float elemental::lapack::HermitianFrobeniusNorm
( Shape shape, const DistMatrix<float,MC,MR>& A );
template double elemental::lapack::HermitianFrobeniusNorm
( Shape shape, const DistMatrix<double,MC,MR>& A );
#ifndef WITHOUT_COMPLEX
template float elemental::lapack::HermitianFrobeniusNorm
( Shape shape, const DistMatrix<std::complex<float>,MC,MR>& A );
template double elemental::lapack::HermitianFrobeniusNorm
( Shape shape, const DistMatrix<std::complex<double>,MC,MR>& A );
#endif
