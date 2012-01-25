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

// TODO: Switch to non-naive methods which are less likely to overflow

template<typename R> 
inline R
internal::FrobeniusNorm( const Matrix<R>& A )
{
#ifndef RELEASE
    PushCallStack("internal::FrobeniusNorm");
#endif
    R normSquared = 0;
    for( int j=0; j<A.Width(); ++j )
    {
        for( int i=0; i<A.Height(); ++i )
        {
            const R alpha = A.Get(i,j);
            normSquared += alpha*alpha;
        }
    }
    R norm = sqrt(normSquared);
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

template<typename R>
inline R
internal::FrobeniusNorm( const Matrix<Complex<R> >& A )
{
#ifndef RELEASE
    PushCallStack("internal::FrobeniusNorm");
#endif
    R normSquared = 0;
    for( int j=0; j<A.Width(); ++j )
    {
        for( int i=0; i<A.Height(); ++i )
        {
            const Complex<R> alpha = A.Get(i,j);
            normSquared += Abs(alpha)*Abs(alpha);
        }
    }
    R norm = sqrt(normSquared);
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

template<typename R> 
inline R
internal::FrobeniusNorm( const DistMatrix<R,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("internal::FrobeniusNorm");
#endif
    R localNormSquared = 0;
    for( int jLocal=0; jLocal<A.LocalWidth(); ++jLocal )
    {
        for( int iLocal=0; iLocal<A.LocalHeight(); ++iLocal )
        {
            const R alpha = A.GetLocalEntry(iLocal,jLocal);
            localNormSquared += alpha*alpha;
        }
    }

    // The norm squared is simply the sum of the local contributions
    R normSquared;
    mpi::AllReduce
    ( &localNormSquared, &normSquared, 1, mpi::SUM, A.Grid().VCComm() );

    R norm = sqrt(normSquared);
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

template<typename R> 
inline R
internal::FrobeniusNorm( const DistMatrix<Complex<R>,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("internal::FrobeniusNorm");
#endif
    R localNormSquared = 0;
    for( int jLocal=0; jLocal<A.LocalWidth(); ++jLocal )
    {
        for( int iLocal=0; iLocal<A.LocalHeight(); ++iLocal )
        {
            const Complex<R> alpha = A.GetLocalEntry(iLocal,jLocal);
            localNormSquared += Abs(alpha)*Abs(alpha);
        }
    }

    // The norm squared is simply the sum of the local contributions
    R normSquared;
    mpi::AllReduce
    ( &localNormSquared, &normSquared, 1, mpi::SUM, A.Grid().VCComm() );

    R norm = sqrt(normSquared);
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

} // namespace elem
