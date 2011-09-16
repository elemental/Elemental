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
DistMatrix<T,STAR,VR>::SetToRandomHermitianHelper<Z>::Func
( DistMatrix<Z,STAR,VR>& parent )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::SetToRandomHermitian");
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
DistMatrix<T,STAR,VR>::SetToRandomHermitianHelper<complex<Z> >::Func
( DistMatrix<complex<Z>,STAR,VR>& parent )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::SetToRandomHermitian");
    parent.AssertNotLockedView();
    if( parent.Height() != parent.Width() )
        throw logic_error( "Hermitian matrices must be square." );
#endif
    const elemental::Grid& g = parent.Grid();
    const int height = parent.Height();
    const int localWidth = parent.LocalWidth();
    const int p = g.Size();
    const int rowShift = parent.RowShift();

    parent.SetToRandom();

    complex<Z>* thisLocalBuffer = parent.LocalBuffer();
    const int thisLDim = parent.LocalLDim();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( int jLocal=0; jLocal<localWidth; ++jLocal )
    {
        const int j = rowShift + jLocal*p;
        if( j < height )
        {
            const Z value = real(thisLocalBuffer[j+jLocal*thisLDim]);
            thisLocalBuffer[j+jLocal*thisLDim] = value;
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
DistMatrix<T,STAR,VR>::SetToRandomHPDHelper<Z>::Func
( DistMatrix<Z,STAR,VR>& parent )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::SetToRandomHPD");
    parent.AssertNotLockedView();
    if( parent.Height() != parent.Width() )
        throw logic_error( "Positive-definite matrices must be square." );
#endif
    const elemental::Grid& g = parent.Grid();
    const int height = parent.Height();
    const int width = parent.Width();
    const int localWidth = parent.LocalWidth();
    const int p = g.Size();
    const int rowShift = parent.RowShift();

    parent.SetToRandom();

    Z* thisLocalBuffer = parent.LocalBuffer();
    const int thisLDim = parent.LocalLDim();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( int jLocal=0; jLocal<localWidth; ++jLocal )
    {
        const int j = rowShift + jLocal*p;
        if( j < height )
            thisLocalBuffer[j+jLocal*thisLDim] += width;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
template<typename T>
template<typename Z>
inline void
DistMatrix<T,STAR,VR>::SetToRandomHPDHelper<complex<Z> >::Func
( DistMatrix<complex<Z>,STAR,VR>& parent )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::SetToRandomHPD");
    parent.AssertNotLockedView();
    if( parent.Height() != parent.Width() )
        throw logic_error( "Positive-definite matrices must be square." );
#endif
    const elemental::Grid& g = parent.Grid();
    const int height = parent.Height();
    const int width = parent.Width();
    const int localWidth = parent.LocalWidth();
    const int p = g.Size();
    const int rowShift = parent.RowShift();

    parent.SetToRandom();

    complex<Z>* thisLocalBuffer = parent.LocalBuffer();
    const int thisLDim = parent.LocalLDim();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( int jLocal=0; jLocal<localWidth; ++jLocal )
    {
        const int j = rowShift + jLocal*p;
        if( j < height )
        {
            const Z value = real(thisLocalBuffer[j+jLocal*thisLDim]);
            thisLocalBuffer[j+jLocal*thisLDim] = value + width;
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
template<typename Z>
inline Z
DistMatrix<T,STAR,VR>::GetRealHelper<complex<Z> >::Func
( const DistMatrix<complex<Z>,STAR,VR>& parent, int i, int j ) 
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::GetReal");
    parent.AssertValidEntry( i, j );
#endif
    // We will determine the owner rank of entry (i,j) and broadcast from that
    // process over the entire g
    const elemental::Grid& g = parent.Grid();
    const int ownerRank = (j + parent.RowAlignment()) % g.Size();

    Z u;
    if( g.VRRank() == ownerRank )
    {
        const int jLoc = (j-parent.RowShift()) / g.Size();
        u = parent.GetRealLocalEntry(i,jLoc);
    }
    mpi::Broadcast( &u, 1, ownerRank, g.VRComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T>
template<typename Z>
inline Z
DistMatrix<T,STAR,VR>::GetImagHelper<complex<Z> >::Func
( const DistMatrix<complex<Z>,STAR,VR>& parent, int i, int j ) 
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::GetImag");
    parent.AssertValidEntry( i, j );
#endif
    // We will determine the owner rank of entry (i,j) and broadcast from that
    // process over the entire g
    const elemental::Grid& g = parent.Grid();
    const int ownerRank = (j + parent.RowAlignment()) % g.Size();

    Z u;
    if( g.VRRank() == ownerRank )
    {
        const int jLoc = (j-parent.RowShift()) / g.Size();
        u = parent.GetImagLocalEntry(i,jLoc);
    }
    mpi::Broadcast( &u, 1, ownerRank, g.VRComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T>
template<typename Z>
inline void
DistMatrix<T,STAR,VR>::SetRealHelper<complex<Z> >::Func
( DistMatrix<complex<Z>,STAR,VR>& parent, int i, int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::SetReal");
    parent.AssertValidEntry( i, j );
#endif
    const elemental::Grid& g = parent.Grid();
    const int ownerRank = (j + parent.RowAlignment()) % g.Size();

    if( g.VRRank() == ownerRank )
    {
        const int jLoc = (j-parent.RowShift()) / g.Size();
        parent.SetRealLocalEntry(i,jLoc,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
template<typename Z>
inline void
DistMatrix<T,STAR,VR>::SetImagHelper<complex<Z> >::Func
( DistMatrix<complex<Z>,STAR,VR>& parent, int i, int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::SetImag");
    parent.AssertValidEntry( i, j );
#endif
    const elemental::Grid& g = parent.Grid();
    const int ownerRank = (j + parent.RowAlignment()) % g.Size();

    if( g.VRRank() == ownerRank )
    {
        const int jLoc = (j-parent.RowShift()) / g.Size();
        parent.SetImagLocalEntry(i,jLoc,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
template<typename Z>
inline void
DistMatrix<T,STAR,VR>::UpdateRealHelper<complex<Z> >::Func
( DistMatrix<complex<Z>,STAR,VR>& parent, int i, int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::UpdateReal");
    parent.AssertValidEntry( i, j );
#endif
    const elemental::Grid& g = parent.Grid();
    const int ownerRank = (j + parent.RowAlignment()) % g.Size();

    if( g.VRRank() == ownerRank )
    {
        const int jLoc = (j-parent.RowShift()) / g.Size();
        parent.UpdateRealLocalEntry(i,jLoc,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
template<typename Z>
inline void
DistMatrix<T,STAR,VR>::UpdateImagHelper<complex<Z> >::Func
( DistMatrix<complex<Z>,STAR,VR>& parent, int i, int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[* ,VR]::UpdateImag");
    parent.AssertValidEntry( i, j );
#endif
    const elemental::Grid& g = parent.Grid();
    const int ownerRank = (j + parent.RowAlignment()) % g.Size();

    if( g.VRRank() == ownerRank )
    {
        const int jLoc = (j-parent.RowShift()) / g.Size();
        parent.UpdateImagLocalEntry(i,jLoc,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX

} // namespace elemental
