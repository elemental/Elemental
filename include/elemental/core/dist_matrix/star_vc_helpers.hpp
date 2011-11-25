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

template<typename T>
template<typename Z>
inline void
DistMatrix<T,STAR,VC>::SetToRandomHermitianHelper<Z>::Func
( DistMatrix<Z,STAR,VC>& parent )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::SetToRandomHermitian");
    parent.AssertNotLockedView();
    if( parent.Height() != parent.Width() )
        throw std::logic_error("Hermitian matrices must be square");
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
DistMatrix<T,STAR,VC>::SetToRandomHermitianHelper<std::complex<Z> >::Func
( DistMatrix<std::complex<Z>,STAR,VC>& parent )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::SetToRandomHermitian");
    parent.AssertNotLockedView();
    if( parent.Height() != parent.Width() )
        throw std::logic_error("Hermitian matrices must be square");
#endif
    const int height     = parent.Height();
    const int localWidth = parent.LocalWidth();
    const int p          = parent.Grid().Size();
    const int rowShift   = parent.RowShift();

    parent.SetToRandom();

    std::complex<Z>* thisLocalBuffer = parent.LocalBuffer();
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
DistMatrix<T,STAR,VC>::SetToRandomHPDHelper<Z>::Func
( DistMatrix<Z,STAR,VC>& parent )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::SetToRandomHPD");
    parent.AssertNotLockedView();
    if( parent.Height() != parent.Width() )
        throw std::logic_error("Positive-definite matrices must be square");
#endif
    const int height     = parent.Height();
    const int width      = parent.Width();
    const int localWidth = parent.LocalWidth();
    const int p          = parent.Grid().Size();
    const int rowShift   = parent.RowShift();

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
DistMatrix<T,STAR,VC>::SetToRandomHPDHelper<std::complex<Z> >::Func
( DistMatrix<std::complex<Z>,STAR,VC>& parent )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::SetToRandomHPD");
    parent.AssertNotLockedView();
    if( parent.Height() != parent.Width() )
        throw std::logic_error("Positive-definite matrices must be square");
#endif
    const int height     = parent.Height();
    const int width      = parent.Width();
    const int localWidth = parent.LocalWidth();
    const int p          = parent.Grid().Size();
    const int rowShift   = parent.RowShift();

    parent.SetToRandom();

    std::complex<Z>* thisLocalBuffer = parent.LocalBuffer();
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
DistMatrix<T,STAR,VC>::GetRealHelper<std::complex<Z> >::Func
( const DistMatrix<std::complex<Z>,STAR,VC>& parent, int i, int j )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::GetReal");
    parent.AssertValidEntry( i, j );
#endif
    // We will determine the owner rank of entry (i,j) and broadcast from that
    // process over the entire g
    const elemental::Grid& g = parent.Grid();
    const int ownerRank = (j + parent.RowAlignment()) % g.Size();

    Z u;
    if( g.VCRank() == ownerRank )
    {
        const int jLoc = (j-parent.RowShift()) / g.Size();
        u = parent.GetRealLocalEntry(i,jLoc);
    }
    mpi::Broadcast( &u, 1, ownerRank, g.VCComm() );
    
#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T>
template<typename Z>
inline Z
DistMatrix<T,STAR,VC>::GetImagHelper<std::complex<Z> >::Func
( const DistMatrix<std::complex<Z>,STAR,VC>& parent, int i, int j ) 
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::GetImag");
    parent.AssertValidEntry( i, j );
#endif
    // We will determine the owner rank of entry (i,j) and broadcast from that
    // process over the entire g
    const elemental::Grid& g = parent.Grid();
    const int ownerRank = (j + parent.RowAlignment()) % g.Size();

    Z u;
    if( g.VCRank() == ownerRank )
    {
        const int jLoc = (j-parent.RowShift()) / g.Size();
        u = parent.GetImagLocalEntry(i,jLoc);
    }
    mpi::Broadcast( &u, 1, ownerRank, g.VCComm() );
    
#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T>
template<typename Z>
inline void
DistMatrix<T,STAR,VC>::SetRealHelper<std::complex<Z> >::Func
( DistMatrix<std::complex<Z>,STAR,VC>& parent, int i, int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::SetReal");
    parent.AssertValidEntry( i, j );
#endif
    const elemental::Grid& g = parent.Grid();
    const int ownerRank = (j + parent.RowAlignment()) % g.Size();

    if( g.VCRank() == ownerRank )
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
DistMatrix<T,STAR,VC>::SetImagHelper<std::complex<Z> >::Func
( DistMatrix<std::complex<Z>,STAR,VC>& parent, int i, int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::SetImag");
    parent.AssertValidEntry( i, j );
#endif
    const elemental::Grid& g = parent.Grid();
    const int ownerRank = (j + parent.RowAlignment()) % g.Size();

    if( g.VCRank() == ownerRank )
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
DistMatrix<T,STAR,VC>::UpdateRealHelper<std::complex<Z> >::Func
( DistMatrix<std::complex<Z>,STAR,VC>& parent, int i, int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::UpdateReal");
    parent.AssertValidEntry( i, j );
#endif
    const elemental::Grid& g = parent.Grid();
    const int ownerRank = (j + parent.RowAlignment()) % g.Size();

    if( g.VCRank() == ownerRank )
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
DistMatrix<T,STAR,VC>::UpdateImagHelper<std::complex<Z> >::Func
( DistMatrix<std::complex<Z>,STAR,VC>& parent, int i, int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[* ,VC]::UpdateImag");
    parent.AssertValidEntry( i, j );
#endif
    const elemental::Grid& g = parent.Grid();
    const int ownerRank = (j + parent.RowAlignment()) % g.Size();

    if( g.VCRank() == ownerRank )
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
