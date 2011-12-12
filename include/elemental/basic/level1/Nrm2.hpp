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

namespace elemental {

template<typename R>
inline R basic::Nrm2( const DistMatrix<R,MC,MR>& x )
{
#ifndef RELEASE
    PushCallStack("basic::Nrm2");
    if( x.Height() != 1 && x.Width() != 1 )
        throw std::logic_error("x must be a vector");
#endif
    R norm;
    const Grid& g = x.Grid();

    if( x.Width() == 1 )
    {
        const int ownerCol = x.RowAlignment();
        if( g.MRRank() == ownerCol )
        {
            R localNorm = Nrm2( x.LockedLocalMatrix() ); 
            
            const int r = g.Height();
            std::vector<R> localNorms(r);
            R* localNormsPtr = &localNorms[0];
            mpi::AllGather( &localNorm, 1, localNormsPtr, 1, g.MCComm() );
            norm = blas::Nrm2( r, localNormsPtr, 1 );
        }
        mpi::Broadcast( &norm, 1, ownerCol, g.MRComm() );
    }
    else
    {
        const int ownerRow = x.ColAlignment();
        if( g.MCRank() == ownerRow )
        {
            R localNorm = Nrm2( x.LockedLocalMatrix() );

            const int c = g.Width();
            std::vector<R> localNorms(c);
            R* localNormsPtr = &localNorms[0];
            mpi::AllGather( &localNorm, 1, localNormsPtr, 1, g.MRComm() );
            norm = blas::Nrm2( c, localNormsPtr, 1 );
        }
        mpi::Broadcast( &norm, 1, ownerRow, g.MCComm() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

template<typename R>
inline R basic::Nrm2( const DistMatrix<std::complex<R>, MC, MR >& x )
{
#ifndef RELEASE
    PushCallStack("basic::Nrm2");
    if( x.Height() != 1 && x.Width() != 1 )
        throw std::logic_error("x must be a vector");
#endif
    R norm;
    const Grid& g = x.Grid();

    if( x.Width() == 1 )
    {
        const int ownerCol = x.RowAlignment();
        if( g.MRRank() == ownerCol )
        {
            R localNorm = Nrm2( x.LockedLocalMatrix() ); 
            
            const int r = g.Height();
            std::vector<R> localNorms(r);
            R* localNormsPtr = &localNorms[0];
            mpi::AllGather( &localNorm, 1, localNormsPtr, 1, g.MCComm() );
            norm = blas::Nrm2( r, localNormsPtr, 1 );
        }
        mpi::Broadcast( &norm, 1, ownerCol, g.MRComm() );
    }
    else
    {
        const int ownerRow = x.ColAlignment();
        if( g.MCRank() == ownerRow )
        {
            R localNorm = Nrm2( x.LockedLocalMatrix() );

            const int c = g.Width();
            std::vector<R> localNorms(c);
            R* localNormsPtr = &localNorms[0];
            mpi::AllGather( &localNorm, 1, localNormsPtr, 1, g.MRComm() );
            norm = blas::Nrm2( c, localNormsPtr, 1 );
        }
        mpi::Broadcast( &norm, 1, ownerRow, g.MCComm() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

} // namespace elemental
