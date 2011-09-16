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
using namespace std;

template<typename T>
template<typename Z>
inline void
DistMatrix<T,MR,STAR>::SetToRandomHermitianHelper<Z>::Func
( DistMatrix<Z,MR,STAR>& parent )
{
#ifndef RELEASE
    PushCallStack("[MR,* ]::SetToRandomHermitian");
    parent.AssertNotLockedView();
    if( parent.Height() != parent.Width() )
        throw logic_error( "Hermitian matrices must be square." );
#endif
    parent.SetToRandom();
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
template<typename T>
template<typename Z>
inline void
DistMatrix<T,MR,STAR>::SetToRandomHermitianHelper<complex<Z> >::Func
( DistMatrix<complex<Z>,MR,STAR>& parent )
{
#ifndef RELEASE
    PushCallStack("[MR,* ]::SetToRandomHermitian");
    parent.AssertNotLockedView();
    if( parent.Height() != parent.Width() )
        throw logic_error( "Hermitian matrices must be square." );
#endif
    const int width = parent.Width();
    const int localHeight = parent.LocalHeight();
    const int c = parent.Grid().Width();
    const int colShift = parent.ColShift();

    parent.SetToRandom();

    complex<Z>* thisLocalBuffer = parent.LocalBuffer();
    const int thisLDim = parent.LocalLDim();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( int iLocal=0; iLocal<localHeight; ++iLocal )
    {
        const int i = colShift + iLocal*c;
        if( i < width )
        {
            const Z value = real(thisLocalBuffer[iLocal+i*thisLDim]);
            thisLocalBuffer[iLocal+i*thisLDim] = value;
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX

template<typename T>
template<typename Z>
inline void
DistMatrix<T,MR,STAR>::SetToRandomHPDHelper<Z>::Func
( DistMatrix<Z,MR,STAR>& parent )
{
#ifndef RELEASE
    PushCallStack("[MR,* ]::SetToRandomHPD");
    parent.AssertNotLockedView();
    if( parent.Height() != parent.Width() )
        throw logic_error( "Positive-definite matrices must be square." );
#endif
    const int width = parent.Width();
    const int localHeight = parent.LocalHeight();
    const int c = parent.Grid().Width();
    const int colShift = parent.ColShift();

    parent.SetToRandom();

    Z* thisLocalBuffer = parent.LocalBuffer();
    const int thisLDim = parent.LocalLDim();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( int iLocal=0; iLocal<localHeight; ++iLocal )
    {
        const int i = colShift + iLocal*c;
        if( i < width )
            thisLocalBuffer[iLocal+i*thisLDim] += width;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
template<typename T>
template<typename Z>
inline void
DistMatrix<T,MR,STAR>::SetToRandomHPDHelper<complex<Z> >::Func
( DistMatrix<complex<Z>,MR,STAR>& parent )
{
#ifndef RELEASE
    PushCallStack("[MR,* ]::SetToRandomHPD");
    parent.AssertNotLockedView();
    if( parent.Height() != parent.Width() )
        throw logic_error( "Positive-definite matrices must be square." );
#endif
    const int width = parent.Width();
    const int localHeight = parent.LocalHeight();
    const int c = parent.Grid().Width();
    const int colShift = parent.ColShift();

    parent.SetToRandom();

    complex<Z>* thisLocalBuffer = parent.LocalBuffer();
    const int thisLDim = parent.LocalLDim();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( int iLocal=0; iLocal<localHeight; ++iLocal )
    {
        const int i = colShift + iLocal*c;
        if( i < width )
        {
            const Z value = real(thisLocalBuffer[iLocal+i*thisLDim]);
            thisLocalBuffer[iLocal+i*thisLDim] = value + width;
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
template<typename Z>
inline Z
DistMatrix<T,MR,STAR>::GetRealHelper<complex<Z> >::Func
( const DistMatrix<complex<Z>,MR,STAR>& parent, int i, int j ) 
{
#ifndef RELEASE
    PushCallStack("[MR,* ]::GetReal");
    parent.AssertValidEntry( i, j );
#endif
    // We will determine the owner column of entry (i,j) and broadcast from that
    // columns within each process row
    const elemental::Grid& g = parent.Grid();
    const int ownerCol = (i + parent.ColAlignment()) % g.Width();

    Z u;
    if( g.MRRank() == ownerCol )
    {
        const int iLoc = (i-parent.ColShift()) / g.Width();
        u = parent.GetRealLocalEntry(iLoc,j);
    }
    mpi::Broadcast( &u, 1, ownerCol, g.MRComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T>
template<typename Z>
inline Z
DistMatrix<T,MR,STAR>::GetImagHelper<complex<Z> >::Func
( const DistMatrix<complex<Z>,MR,STAR>& parent, int i, int j ) 
{
#ifndef RELEASE
    PushCallStack("[MR,* ]::GetImag");
    parent.AssertValidEntry( i, j );
#endif
    // We will determine the owner column of entry (i,j) and broadcast from that
    // columns within each process row
    const elemental::Grid& g = parent.Grid();
    const int ownerCol = (i + parent.ColAlignment()) % g.Width();

    Z u;
    if( g.MRRank() == ownerCol )
    {
        const int iLoc = (i-parent.ColShift()) / g.Width();
        u = parent.GetImagLocalEntry(iLoc,j);
    }
    mpi::Broadcast( &u, 1, ownerCol, g.MRComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T>
template<typename Z>
inline void
DistMatrix<T,MR,STAR>::SetRealHelper<complex<Z> >::Func
( DistMatrix<complex<Z>,MR,STAR>& parent, int i, int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[MR,* ]::SetReal");
    parent.AssertValidEntry( i, j );
#endif
    const elemental::Grid& g = parent.Grid();
    const int ownerCol = (i + parent.ColAlignment()) % g.Width();

    if( g.MRRank() == ownerCol )
    {
        const int iLoc = (i-parent.ColShift()) / g.Width();
        parent.SetRealLocalEntry(iLoc,j,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
template<typename Z>
inline void
DistMatrix<T,MR,STAR>::SetImagHelper<complex<Z> >::Func
( DistMatrix<complex<Z>,MR,STAR>& parent, int i, int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[MR,* ]::SetImag");
    parent.AssertValidEntry( i, j );
#endif
    const elemental::Grid& g = parent.Grid();
    const int ownerCol = (i + parent.ColAlignment()) % g.Width();

    if( g.MRRank() == ownerCol )
    {
        const int iLoc = (i-parent.ColShift()) / g.Width();
        parent.SetImagLocalEntry(iLoc,j,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
template<typename Z>
inline void
DistMatrix<T,MR,STAR>::UpdateRealHelper<complex<Z> >::Func
( DistMatrix<complex<Z>,MR,STAR>& parent, int i, int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[MR,* ]::UpdateReal");
    parent.AssertValidEntry( i, j );
#endif
    const elemental::Grid& g = parent.Grid();
    const int ownerCol = (i + parent.ColAlignment()) % g.Width();

    if( g.MRRank() == ownerCol )
    {
        const int iLoc = (i-parent.ColShift()) / g.Width();
        parent.UpdateRealLocalEntry(iLoc,j,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
template<typename Z>
inline void
DistMatrix<T,MR,STAR>::UpdateImagHelper<complex<Z> >::Func
( DistMatrix<complex<Z>,MR,STAR>& parent, int i, int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[MR,* ]::UpdateImag");
    parent.AssertValidEntry( i, j );
#endif
    const elemental::Grid& g = parent.Grid();
    const int ownerCol = (i + parent.ColAlignment()) % g.Width();

    if( g.MRRank() == ownerCol )
    {
        const int iLoc = (i-parent.ColShift()) / g.Width();
        parent.UpdateImagLocalEntry(iLoc,j,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX

} // namespace elemental
