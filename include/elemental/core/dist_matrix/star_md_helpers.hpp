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

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,MD,Int>::SetToRandomHermitian()
{ SetToRandomHermitianHelper<T>::Func( *this ); }

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,MD,Int>::SetToRandomHPD()
{ SetToRandomHPDHelper<T>::Func( *this ); }

template<typename T,typename Int>
inline typename Base<T>::type
DistMatrix<T,STAR,MD,Int>::GetReal( Int i, Int j ) const
{ return GetRealHelper<T>::Func( *this, i, j ); }

template<typename T,typename Int>
template<typename Z>
inline Z
DistMatrix<T,STAR,MD,Int>::GetRealHelper<Z>::Func
( const DistMatrix<Z,STAR,MD,Int>& parent, Int i, Int j )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::GetRealHelper");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T,typename Int>
inline typename Base<T>::type
DistMatrix<T,STAR,MD,Int>::GetImag( Int i, Int j ) const
{ return GetImagHelper<T>::Func( *this, i, j ); }

template<typename T,typename Int>
template<typename Z>
inline Z
DistMatrix<T,STAR,MD,Int>::GetImagHelper<Z>::Func
( const DistMatrix<Z,STAR,MD,Int>& parent, Int i, Int j )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::GetImag");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,MD,Int>::SetReal
( Int i, Int j, typename Base<T>::type alpha )
{ SetRealHelper<T>::Func( *this, i, j, alpha ); }

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,STAR,MD,Int>::SetRealHelper<Z>::Func
( DistMatrix<Z,STAR,MD,Int>& parent, Int i, Int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::SetReal");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,MD,Int>::SetImag
( Int i, Int j, typename Base<T>::type alpha )
{ SetImagHelper<T>::Func( *this, i, j, alpha ); }

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,STAR,MD,Int>::SetImagHelper<Z>::Func
( DistMatrix<Z,STAR,MD,Int>& parent, Int i, Int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::SetImag");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,MD,Int>::UpdateReal
( Int i, Int j, typename Base<T>::type alpha )
{ UpdateRealHelper<T>::Func( *this, i, j, alpha ); }

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,STAR,MD,Int>::UpdateRealHelper<Z>::Func
( DistMatrix<Z,STAR,MD,Int>& parent, Int i, Int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::UpdateReal");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,MD,Int>::UpdateImag
( Int i, Int j, typename Base<T>::type alpha )
{ UpdateImagHelper<T>::Func( *this, i, j, alpha ); }

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,STAR,MD,Int>::UpdateImagHelper<Z>::Func
( DistMatrix<Z,STAR,MD,Int>& parent, Int i, Int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::UpdateImag");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,STAR,MD,Int>::SetToRandomHermitianHelper<Z>::Func
( DistMatrix<Z,STAR,MD,Int>& parent )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::SetToRandomHermitian");
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
DistMatrix<T,STAR,MD,Int>::SetToRandomHermitianHelper<Complex<Z> >::Func
( DistMatrix<Complex<Z>,STAR,MD,Int>& parent )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::SetToRandomHermitian");
    parent.AssertNotLockedView();
    if( parent.Height() != parent.Width() )
        throw std::logic_error("Hermitian matrices must be square");
#endif
    parent.SetToRandom();

    if( parent.InDiagonal() )
    {
        const Int height = parent.Height();
        const Int localWidth = parent.LocalWidth();
        const Int lcm = parent.Grid().LCM();
        const Int rowShift = parent.RowShift();

        Complex<Z>* thisLocalBuffer = parent.LocalBuffer();
        const Int thisLDim = parent.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( Int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            const Int j = rowShift + jLocal*lcm;
            if( j < height )
            {
                const Z value = thisLocalBuffer[j+jLocal*thisLDim].real;
                thisLocalBuffer[j+jLocal*thisLDim] = value;
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
DistMatrix<T,STAR,MD,Int>::SetToRandomHPDHelper<Z>::Func
( DistMatrix<Z,STAR,MD,Int>& parent )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::SetToRandomHPD");
    parent.AssertNotLockedView();
    if( parent.Height() != parent.Width() )
        throw std::logic_error("Positive-definite matrices must be square");
#endif
    parent.SetToRandom();

    if( parent.InDiagonal() )
    {
        const Int height = parent.Height();
        const Int width = parent.Width();
        const Int localWidth = parent.LocalWidth();
        const Int lcm = parent.Grid().LCM();
        const Int rowShift = parent.RowShift();

        Z* thisLocalBuffer = parent.LocalBuffer();
        const Int thisLDim = parent.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( Int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            const Int j = rowShift + jLocal*lcm;
            if( j < height )
                thisLocalBuffer[j+jLocal*thisLDim] += width;
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,STAR,MD,Int>::SetToRandomHPDHelper<Complex<Z> >::Func
( DistMatrix<Complex<Z>,STAR,MD,Int>& parent )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::SetToRandomHPD");
    parent.AssertNotLockedView();
    if( parent.Height() != parent.Width() )
        throw std::logic_error("Positive-definite matrices must be square");
#endif
    parent.SetToRandom();

    if( parent.InDiagonal() )
    {
        const Int height = parent.Height();
        const Int width = parent.Width();
        const Int localWidth = parent.LocalWidth();
        const Int lcm = parent.Grid().LCM();
        const Int rowShift = parent.RowShift();

        Complex<Z>* thisLocalBuffer = parent.LocalBuffer();
        const Int thisLDim = parent.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( Int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            const Int j = rowShift + jLocal*lcm;
            if( j < height )
            {
                const Z value = thisLocalBuffer[j+jLocal*thisLDim].real;
                thisLocalBuffer[j+jLocal*thisLDim] = value + width;
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
DistMatrix<T,STAR,MD,Int>::GetRealHelper<Complex<Z> >::Func
( const DistMatrix<Complex<Z>,STAR,MD,Int>& parent, Int i, Int j ) 
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::GetReal");
    parent.AssertValidEntry( i, j );
#endif
    // We will determine the owner of entry (i,j) and broadcast from it
    Int ownerRank;
    const elem::Grid& g = parent.Grid();
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int alignmentRank = parent.RowAlignment();
        const Int alignmentRow = alignmentRank % r;
        const Int alignmentCol = alignmentRank / r;
        const Int ownerRow = (alignmentRow + j) % r;
        const Int ownerCol = (alignmentCol + j) % c;
        ownerRank = ownerRow + r*ownerCol;
    }

    Z u;
    if( g.VCRank() == ownerRank )
    {
        const Int jLoc = (j-parent.RowShift()) / g.LCM();
        u = parent.GetRealLocalEntry(i,jLoc);
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
DistMatrix<T,STAR,MD,Int>::GetImagHelper<Complex<Z> >::Func
( const DistMatrix<Complex<Z>,STAR,MD,Int>& parent, Int i, Int j )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::GetImag");
    parent.AssertValidEntry( i, j );
#endif
    // We will determine the owner of entry (i,j) and broadcast from it
    Int ownerRank;
    const elem::Grid& g = parent.Grid();
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int alignmentRank = parent.RowAlignment();
        const Int alignmentRow = alignmentRank % r;
        const Int alignmentCol = alignmentRank / r;
        const Int ownerRow = (alignmentRow + j) % r;
        const Int ownerCol = (alignmentCol + j) % c;
        ownerRank = ownerRow + r*ownerCol;
    }

    Z u;
    if( g.VCRank() == ownerRank )
    {
        const Int jLoc = (j-parent.RowShift()) / g.LCM();
        u = parent.GetImagLocalEntry(i,jLoc);
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
DistMatrix<T,STAR,MD,Int>::SetRealHelper<Complex<Z> >::Func
( DistMatrix<Complex<Z>,STAR,MD,Int>& parent, Int i, Int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::SetReal");
    parent.AssertValidEntry( i, j );
#endif
    Int ownerRank;
    const elem::Grid& g = parent.Grid();
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int alignmentRank = parent.RowAlignment();
        const Int alignmentRow = alignmentRank % r;
        const Int alignmentCol = alignmentRank / r;
        const Int ownerRow = (alignmentRow + j) % r;
        const Int ownerCol = (alignmentCol + j) % c;
        ownerRank = ownerRow + r*ownerCol;
    }

    if( g.VCRank() == ownerRank )
    {
        const Int jLoc = (j-parent.RowShift()) / g.LCM();
        parent.SetRealLocalEntry(i,jLoc,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,STAR,MD,Int>::SetImagHelper<Complex<Z> >::Func
( DistMatrix<Complex<Z>,STAR,MD,Int>& parent, Int i, Int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::SetImag");
    parent.AssertValidEntry( i, j );
#endif
    Int ownerRank;
    const elem::Grid& g = parent.Grid();
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int alignmentRank = parent.RowAlignment();
        const Int alignmentRow = alignmentRank % r;
        const Int alignmentCol = alignmentRank / r;
        const Int ownerRow = (alignmentRow + j) % r;
        const Int ownerCol = (alignmentCol + j) % c;
        ownerRank = ownerRow + r*ownerCol;
    }

    if( g.VCRank() == ownerRank )
    {
        const Int jLoc = (j-parent.RowShift()) / g.LCM();
        parent.SetImagLocalEntry(i,jLoc,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,STAR,MD,Int>::UpdateRealHelper<Complex<Z> >::Func
( DistMatrix<Complex<Z>,STAR,MD,Int>& parent, Int i, Int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::UpdateReal");
    parent.AssertValidEntry( i, j );
#endif
    Int ownerRank;
    const elem::Grid& g = parent.Grid();
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int alignmentRank = parent.RowAlignment();
        const Int alignmentRow = alignmentRank % r;
        const Int alignmentCol = alignmentRank / r;
        const Int ownerRow = (alignmentRow + j) % r;
        const Int ownerCol = (alignmentCol + j) % c;
        ownerRank = ownerRow + r*ownerCol;
    }

    if( g.VCRank() == ownerRank )
    {
        const Int jLoc = (j-parent.RowShift()) / g.LCM();
        parent.UpdateRealLocalEntry(i,jLoc,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,STAR,MD,Int>::UpdateImagHelper<Complex<Z> >::Func
( DistMatrix<Complex<Z>,STAR,MD,Int>& parent, Int i, Int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::UpdateImag");
    parent.AssertValidEntry( i, j );
#endif
    Int ownerRank;
    const elem::Grid& g = parent.Grid();
    {
        const Int r = g.Height();
        const Int c = g.Width();
        const Int alignmentRank = parent.RowAlignment();
        const Int alignmentRow = alignmentRank % r;
        const Int alignmentCol = alignmentRank / r;
        const Int ownerRow = (alignmentRow + j) % r;
        const Int ownerCol = (alignmentCol + j) % c;
        ownerRank = ownerRow + r*ownerCol;
    }

    if( g.VCRank() == ownerRank )
    {
        const Int jLoc = (j-parent.RowShift()) / g.LCM();
        parent.UpdateImagLocalEntry(i,jLoc,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
