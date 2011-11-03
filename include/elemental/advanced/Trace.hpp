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

template<typename F> 
inline F advanced::Trace( const DistMatrix<F,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("advanced::Trace");
#endif
    if( A.Height() != A.Width() )
        throw std::logic_error("Cannot compute trace of nonsquare matrix");
    const Grid& g = A.Grid();

    DistMatrix<F,MD,STAR> d(g);
    A.GetDiagonal( d );
    F localTrace = 0;
    if( d.InDiagonal() )
    {
        const int nLocalDiag = d.LocalHeight();
        for( int iLocal=0; iLocal<nLocalDiag; ++iLocal )
            localTrace += d.GetLocalEntry(iLocal,0);
    }
    F trace;
    mpi::AllReduce( &localTrace, &trace, 1, mpi::SUM, g.VCComm() );
#ifndef RELEASE
    PopCallStack();
#endif
    return trace;
}

template<typename F>
inline F advanced::Trace( const Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("advanced::Trace");
#endif
    if( A.Height() != A.Width() )
        throw std::logic_error("Cannot compute trace of nonsquare matrix");

    Matrix<F> d;
    A.GetDiagonal( d );
    F trace = 0;
    const int n = A.Height();
    for( int i=0; i<n; ++i )
        trace += d.Get(i,0);
#ifndef RELEASE
    PopCallStack();
#endif
    return trace;
}
