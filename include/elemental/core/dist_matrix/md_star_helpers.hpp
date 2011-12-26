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

template<typename T,typename Int>
inline void
DistMatrix<T,MD,STAR,Int>::SetToRandomHermitian()
{ SetToRandomHermitianHelper<T>::Func( *this ); }

template<typename T,typename Int>
inline void
DistMatrix<T,MD,STAR,Int>::SetToRandomHPD()
{ SetToRandomHPDHelper<T>::Func( *this ); }

template<typename T,typename Int>
inline typename RealBase<T>::type
DistMatrix<T,MD,STAR,Int>::GetReal( Int i, Int j ) const
{ return GetRealHelper<T>::Func( *this, i, j ); }

template<typename T,typename Int>
template<typename Z>
inline Z
DistMatrix<T,MD,STAR,Int>::GetRealHelper<Z>::Func
( const DistMatrix<Z,MD,STAR,Int>& parent, Int i, Int j )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::GetRealHelper");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T,typename Int>
inline typename RealBase<T>::type
DistMatrix<T,MD,STAR,Int>::GetImag( Int i, Int j ) const
{ return GetImagHelper<T>::Func( *this, i, j ); }

template<typename T,typename Int>
template<typename Z>
inline Z
DistMatrix<T,MD,STAR,Int>::GetImagHelper<Z>::Func
( const DistMatrix<Z,MD,STAR,Int>& parent, Int i, Int j )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::GetImag");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T,typename Int>
inline void
DistMatrix<T,MD,STAR,Int>::SetReal
( Int i, Int j, typename RealBase<T>::type alpha )
{ SetRealHelper<T>::Func( *this, i, j, alpha ); }

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,MD,STAR,Int>::SetRealHelper<Z>::Func
( DistMatrix<Z,MD,STAR,Int>& parent, Int i, Int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::SetReal");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T,typename Int>
inline void
DistMatrix<T,MD,STAR,Int>::SetImag
( Int i, Int j, typename RealBase<T>::type alpha )
{ SetImagHelper<T>::Func( *this, i, j, alpha ); }

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,MD,STAR,Int>::SetImagHelper<Z>::Func
( DistMatrix<Z,MD,STAR,Int>& parent, Int i, Int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::SetImag");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T,typename Int>
inline void
DistMatrix<T,MD,STAR,Int>::UpdateReal
( Int i, Int j, typename RealBase<T>::type alpha )
{ UpdateRealHelper<T>::Func( *this, i, j, alpha ); }

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,MD,STAR,Int>::UpdateRealHelper<Z>::Func
( DistMatrix<Z,MD,STAR,Int>& parent, Int i, Int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::UpdateReal");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T,typename Int>
inline void
DistMatrix<T,MD,STAR,Int>::UpdateImag
( Int i, Int j, typename RealBase<T>::type alpha )
{ UpdateImagHelper<T>::Func( *this, i, j, alpha ); }

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,MD,STAR,Int>::UpdateImagHelper<Z>::Func
( DistMatrix<Z,MD,STAR,Int>& parent, Int i, Int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::UpdateImag");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}


template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,MD,STAR,Int>::SetToRandomHermitianHelper<Z>::Func
( DistMatrix<Z,MD,STAR,Int>& parent )
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

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,MD,STAR,Int>::SetToRandomHermitianHelper<std::complex<Z> >::Func
( DistMatrix<std::complex<Z>,MD,STAR,Int>& parent )
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
        const Int width = parent.Width();
        const Int localHeight = parent.LocalHeight();
        const Int lcm = parent.Grid().LCM();
        const Int colShift = parent.ColShift();

        std::complex<Z>* thisLocalBuffer = parent.LocalBuffer();
        const Int thisLDim = parent.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( Int iLocal=0; iLocal<localHeight; ++iLocal )
        {
            const Int i = colShift + iLocal*lcm;
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

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,MD,STAR,Int>::SetToRandomHPDHelper<Z>::Func
( DistMatrix<Z,MD,STAR,Int>& parent )
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
        const Int width = parent.Width();
        const Int localHeight = parent.LocalHeight();
        const Int lcm = parent.Grid().LCM();
        const Int colShift = parent.ColShift();

        Z* thisLocalBuffer = parent.LocalBuffer();
        const Int thisLDim = parent.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( Int iLocal=0; iLocal<localHeight; ++iLocal )
        {
            const Int i = colShift + iLocal*lcm;
            if( i < width )
                thisLocalBuffer[iLocal+i*thisLDim] += width;
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,MD,STAR,Int>::SetToRandomHPDHelper<std::complex<Z> >::Func
( DistMatrix<std::complex<Z>,MD,STAR,Int>& parent )
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
        const Int width = parent.Width();
        const Int localHeight = parent.LocalHeight();
        const Int lcm = parent.Grid().LCM();
        const Int colShift = parent.ColShift();

        std::complex<Z>* thisLocalBuffer = parent.LocalBuffer();
        const Int thisLDim = parent.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( Int iLocal=0; iLocal<localHeight; ++iLocal )
        {
            const Int i = colShift + iLocal*lcm;
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

template<typename T,typename Int>
template<typename Z>
inline Z
DistMatrix<T,MD,STAR,Int>::GetRealHelper<std::complex<Z> >::Func
( const DistMatrix<std::complex<Z>,MD,STAR,Int>& parent, Int i, Int j ) 
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::GetReal");
    parent.AssertValidEntry( i, j );
#endif
    // We will determine the owner of entry (i,j) and broadcast from it
    Int ownerRank;
    const elemental::Grid& g = parent.Grid();
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int alignmentRank = parent.ColAlignment();
        const Int alignmentRow = alignmentRank % r;
        const Int alignmentCol = alignmentRank / r;
        const Int ownerRow = (alignmentRow + i) % r;
        const Int ownerCol = (alignmentCol + i) % c;
        ownerRank = ownerRow + r*ownerCol;
    }

    Z u;
    if( g.VCRank() == ownerRank )
    {
        const Int iLoc = (i-parent.ColShift()) / g.LCM();
        u = parent.GetRealLocalEntry(iLoc,j);
    }
    mpi::Broadcast( &u, 1, ownerRank, g.VCComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T,typename Int>
template<typename Z>
inline Z
DistMatrix<T,MD,STAR,Int>::GetImagHelper<std::complex<Z> >::Func
( const DistMatrix<std::complex<Z>,MD,STAR,Int>& parent, Int i, Int j ) 
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::GetImag");
    parent.AssertValidEntry( i, j );
#endif
    // We will determine the owner of entry (i,j) and broadcast from it
    Int ownerRank;
    const elemental::Grid& g = parent.Grid();
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int alignmentRank = parent.ColAlignment();
        const Int alignmentRow = alignmentRank % r;
        const Int alignmentCol = alignmentRank / r;
        const Int ownerRow = (alignmentRow + i) % r;
        const Int ownerCol = (alignmentCol + i) % c;
        ownerRank = ownerRow + r*ownerCol;
    }

    Z u;
    if( g.VCRank() == ownerRank )
    {
        const Int iLoc = (i-parent.ColShift()) / g.LCM();
        u = parent.GetImagLocalEntry(iLoc,j);
    }
    mpi::Broadcast( &u, 1, ownerRank, g.VCComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,MD,STAR,Int>::SetRealHelper<std::complex<Z> >::Func
( DistMatrix<std::complex<Z>,MD,STAR,Int>& parent, Int i, Int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::SetReal");
    parent.AssertValidEntry( i, j );
#endif
    Int ownerRank;
    const elemental::Grid& g = parent.Grid();
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int alignmentRank = parent.ColAlignment();
        const Int alignmentRow = alignmentRank % r;
        const Int alignmentCol = alignmentRank / r;
        const Int ownerRow = (alignmentRow + i) % r;
        const Int ownerCol = (alignmentCol + i) % c;
        ownerRank = ownerRow + r*ownerCol;
    }

    if( g.VCRank() == ownerRank )
    {
        const Int iLoc = (i-parent.ColShift()) / g.LCM();
        parent.SetRealLocalEntry(iLoc,j,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,MD,STAR,Int>::SetImagHelper<std::complex<Z> >::Func
( DistMatrix<std::complex<Z>,MD,STAR,Int>& parent, Int i, Int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::SetImag");
    parent.AssertValidEntry( i, j );
#endif
    Int ownerRank;
    const elemental::Grid& g = parent.Grid();
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int alignmentRank = parent.ColAlignment();
        const Int alignmentRow = alignmentRank % r;
        const Int alignmentCol = alignmentRank / r;
        const Int ownerRow = (alignmentRow + i) % r;
        const Int ownerCol = (alignmentCol + i) % c;
        ownerRank = ownerRow + r*ownerCol;
    }

    if( g.VCRank() == ownerRank )
    {
        const Int iLoc = (i-parent.ColShift()) / g.LCM();
        parent.SetImagLocalEntry(iLoc,j,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,MD,STAR,Int>::UpdateRealHelper<std::complex<Z> >::Func
( DistMatrix<std::complex<Z>,MD,STAR,Int>& parent, Int i, Int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::UpdateReal");
    parent.AssertValidEntry( i, j );
#endif
    Int ownerRank;
    const elemental::Grid& g = parent.Grid();
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int alignmentRank = parent.ColAlignment();
        const Int alignmentRow = alignmentRank % r;
        const Int alignmentCol = alignmentRank / r;
        const Int ownerRow = (alignmentRow + i) % r;
        const Int ownerCol = (alignmentCol + i) % c;
        ownerRank = ownerRow + r*ownerCol;
    }

    if( g.VCRank() == ownerRank )
    {
        const Int iLoc = (i-parent.ColShift()) / g.LCM();
        parent.UpdateRealLocalEntry(iLoc,j,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,MD,STAR,Int>::UpdateImagHelper<std::complex<Z> >::Func
( DistMatrix<std::complex<Z>,MD,STAR,Int>& parent, Int i, Int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::UpdateImag");
    parent.AssertValidEntry( i, j );
#endif
    Int ownerRank;
    const elemental::Grid& g = parent.Grid();
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int alignmentRank = parent.ColAlignment();
        const Int alignmentRow = alignmentRank % r;
        const Int alignmentCol = alignmentRank / r;
        const Int ownerRow = (alignmentRow + i) % r;
        const Int ownerCol = (alignmentCol + i) % c;
        ownerRank = ownerRow + r*ownerCol;
    }

    if( g.VCRank() == ownerRank )
    {
        const Int iLoc = (i-parent.ColShift()) / g.LCM();
        parent.UpdateImagLocalEntry(iLoc,j,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elemental
