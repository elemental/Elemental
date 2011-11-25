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
DistMatrix<T,MD,STAR>::SetToRandomHermitianHelper<Z>::Func
( DistMatrix<Z,MD,STAR>& parent )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::SetToRandomHermitian");
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
DistMatrix<T,MD,STAR>::SetToRandomHermitianHelper<std::complex<Z> >::Func
( DistMatrix<std::complex<Z>,MD,STAR>& parent )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::SetToRandomHermitian");
    parent.AssertNotLockedView();
    if( parent.Height() != parent.Width() )
        throw std::logic_error("Hermitian matrices must be square");
#endif
    parent.SetToRandom();

    if( parent.InDiagonal() )
    {
        const int width = parent.Width();
        const int localHeight = parent.LocalHeight();
        const int lcm = parent.Grid().LCM();
        const int colShift = parent.ColShift();

        std::complex<Z>* thisLocalBuffer = parent.LocalBuffer();
        const int thisLDim = parent.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int iLocal=0; iLocal<localHeight; ++iLocal )
        {
            const int i = colShift + iLocal*lcm;
            if( i < width )
            {
                const Z value = real(thisLocalBuffer[iLocal+i*thisLDim]);
                thisLocalBuffer[iLocal+i*thisLDim] = value;
            }
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
DistMatrix<T,MD,STAR>::SetToRandomHPDHelper<Z>::Func
( DistMatrix<Z,MD,STAR>& parent )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::SetToRandomHPD");
    parent.AssertNotLockedView();
    if( parent.Height() != parent.Width() )
        throw std::logic_error("Positive-definite matrices must be square");
#endif
    parent.SetToRandom();

    if( parent.InDiagonal() )
    {
        const int width = parent.Width();
        const int localHeight = parent.LocalHeight();
        const int lcm = parent.Grid().LCM();
        const int colShift = parent.ColShift();

        Z* thisLocalBuffer = parent.LocalBuffer();
        const int thisLDim = parent.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int iLocal=0; iLocal<localHeight; ++iLocal )
        {
            const int i = colShift + iLocal*lcm;
            if( i < width )
                thisLocalBuffer[iLocal+i*thisLDim] += width;
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
template<typename T>
template<typename Z>
inline void
DistMatrix<T,MD,STAR>::SetToRandomHPDHelper<std::complex<Z> >::Func
( DistMatrix<std::complex<Z>,MD,STAR>& parent )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::SetToRandomHPD");
    parent.AssertNotLockedView();
    if( parent.Height() != parent.Width() )
        throw std::logic_error("Positive-definite matrices must be square");
#endif
    parent.SetToRandom();

    if( parent.InDiagonal() )
    {
        const int width = parent.Width();
        const int localHeight = parent.LocalHeight();
        const int lcm = parent.Grid().LCM();
        const int colShift = parent.ColShift();

        std::complex<Z>* thisLocalBuffer = parent.LocalBuffer();
        const int thisLDim = parent.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int iLocal=0; iLocal<localHeight; ++iLocal )
        {
            const int i = colShift + iLocal*lcm;
            if( i < width )
            {
                const Z value = real(thisLocalBuffer[iLocal+i*thisLDim]);
                thisLocalBuffer[iLocal+i*thisLDim] = value + width;
            }
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
template<typename Z>
inline Z
DistMatrix<T,MD,STAR>::GetRealHelper<std::complex<Z> >::Func
( const DistMatrix<std::complex<Z>,MD,STAR>& parent, int i, int j ) 
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::GetReal");
    parent.AssertValidEntry( i, j );
#endif
    // We will determine the owner of entry (i,j) and broadcast from it
    int ownerRank;
    const elemental::Grid& g = parent.Grid();
    {
        const int r = g.Height();
        const int c = g.Width();
        const int alignmentRank = parent.ColAlignment();
        const int alignmentRow = alignmentRank % r;
        const int alignmentCol = alignmentRank / r;
        const int ownerRow = (alignmentRow + i) % r;
        const int ownerCol = (alignmentCol + i) % c;
        ownerRank = ownerRow + r*ownerCol;
    }

    Z u;
    if( g.VCRank() == ownerRank )
    {
        const int iLoc = (i-parent.ColShift()) / g.LCM();
        u = parent.GetRealLocalEntry(iLoc,j);
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
DistMatrix<T,MD,STAR>::GetImagHelper<std::complex<Z> >::Func
( const DistMatrix<std::complex<Z>,MD,STAR>& parent, int i, int j ) 
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::GetImag");
    parent.AssertValidEntry( i, j );
#endif
    // We will determine the owner of entry (i,j) and broadcast from it
    int ownerRank;
    const elemental::Grid& g = parent.Grid();
    {
        const int r = g.Height();
        const int c = g.Width();
        const int alignmentRank = parent.ColAlignment();
        const int alignmentRow = alignmentRank % r;
        const int alignmentCol = alignmentRank / r;
        const int ownerRow = (alignmentRow + i) % r;
        const int ownerCol = (alignmentCol + i) % c;
        ownerRank = ownerRow + r*ownerCol;
    }

    Z u;
    if( g.VCRank() == ownerRank )
    {
        const int iLoc = (i-parent.ColShift()) / g.LCM();
        u = parent.GetImagLocalEntry(iLoc,j);
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
DistMatrix<T,MD,STAR>::SetRealHelper<std::complex<Z> >::Func
( DistMatrix<std::complex<Z>,MD,STAR>& parent, int i, int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::SetReal");
    parent.AssertValidEntry( i, j );
#endif
    int ownerRank;
    const elemental::Grid& g = parent.Grid();
    {
        const int r = g.Height();
        const int c = g.Width();
        const int alignmentRank = parent.ColAlignment();
        const int alignmentRow = alignmentRank % r;
        const int alignmentCol = alignmentRank / r;
        const int ownerRow = (alignmentRow + i) % r;
        const int ownerCol = (alignmentCol + i) % c;
        ownerRank = ownerRow + r*ownerCol;
    }

    if( g.VCRank() == ownerRank )
    {
        const int iLoc = (i-parent.ColShift()) / g.LCM();
        parent.SetRealLocalEntry(iLoc,j,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
template<typename Z>
inline void
DistMatrix<T,MD,STAR>::SetImagHelper<std::complex<Z> >::Func
( DistMatrix<std::complex<Z>,MD,STAR>& parent, int i, int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::SetImag");
    parent.AssertValidEntry( i, j );
#endif
    int ownerRank;
    const elemental::Grid& g = parent.Grid();
    {
        const int r = g.Height();
        const int c = g.Width();
        const int alignmentRank = parent.ColAlignment();
        const int alignmentRow = alignmentRank % r;
        const int alignmentCol = alignmentRank / r;
        const int ownerRow = (alignmentRow + i) % r;
        const int ownerCol = (alignmentCol + i) % c;
        ownerRank = ownerRow + r*ownerCol;
    }

    if( g.VCRank() == ownerRank )
    {
        const int iLoc = (i-parent.ColShift()) / g.LCM();
        parent.SetImagLocalEntry(iLoc,j,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
template<typename Z>
inline void
DistMatrix<T,MD,STAR>::UpdateRealHelper<std::complex<Z> >::Func
( DistMatrix<std::complex<Z>,MD,STAR>& parent, int i, int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::UpdateReal");
    parent.AssertValidEntry( i, j );
#endif
    int ownerRank;
    const elemental::Grid& g = parent.Grid();
    {
        const int r = g.Height();
        const int c = g.Width();
        const int alignmentRank = parent.ColAlignment();
        const int alignmentRow = alignmentRank % r;
        const int alignmentCol = alignmentRank / r;
        const int ownerRow = (alignmentRow + i) % r;
        const int ownerCol = (alignmentCol + i) % c;
        ownerRank = ownerRow + r*ownerCol;
    }

    if( g.VCRank() == ownerRank )
    {
        const int iLoc = (i-parent.ColShift()) / g.LCM();
        parent.UpdateRealLocalEntry(iLoc,j,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
template<typename Z>
inline void
DistMatrix<T,MD,STAR>::UpdateImagHelper<std::complex<Z> >::Func
( DistMatrix<std::complex<Z>,MD,STAR>& parent, int i, int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::UpdateImag");
    parent.AssertValidEntry( i, j );
#endif
    int ownerRank;
    const elemental::Grid& g = parent.Grid();
    {
        const int r = g.Height();
        const int c = g.Width();
        const int alignmentRank = parent.ColAlignment();
        const int alignmentRow = alignmentRank % r;
        const int alignmentCol = alignmentRank / r;
        const int ownerRow = (alignmentRow + i) % r;
        const int ownerCol = (alignmentCol + i) % c;
        ownerRank = ownerRow + r*ownerCol;
    }

    if( g.VCRank() == ownerRank )
    {
        const int iLoc = (i-parent.ColShift()) / g.LCM();
        parent.UpdateImagLocalEntry(iLoc,j,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX

} // namespace elemental
