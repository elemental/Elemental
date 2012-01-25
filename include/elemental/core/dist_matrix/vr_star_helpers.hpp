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
DistMatrix<T,VR,STAR,Int>::SetToRandomHermitian()
{ SetToRandomHermitianHelper<T>::Func( *this ); }

template<typename T,typename Int>
inline void
DistMatrix<T,VR,STAR,Int>::SetToRandomHPD()
{ SetToRandomHPDHelper<T>::Func( *this ); }

template<typename T,typename Int>
inline typename Base<T>::type
DistMatrix<T,VR,STAR,Int>::GetReal( Int i, Int j ) const
{ return GetRealHelper<T>::Func( *this, i, j ); }

template<typename T,typename Int>
template<typename Z>
inline Z
DistMatrix<T,VR,STAR,Int>::GetRealHelper<Z>::Func
( const DistMatrix<Z,VR,STAR,Int>& parent, Int i, Int j )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::GetRealHelper");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T,typename Int>
inline typename Base<T>::type
DistMatrix<T,VR,STAR,Int>::GetImag( Int i, Int j ) const
{ return GetImagHelper<T>::Func( *this, i, j ); }

template<typename T,typename Int>
template<typename Z>
inline Z
DistMatrix<T,VR,STAR,Int>::GetImagHelper<Z>::Func
( const DistMatrix<Z,VR,STAR,Int>& parent, Int i, Int j )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::GetImag");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T,typename Int>
inline void
DistMatrix<T,VR,STAR,Int>::SetReal( Int i, Int j, typename Base<T>::type alpha )
{ SetRealHelper<T>::Func( *this, i, j, alpha ); }

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,VR,STAR,Int>::SetRealHelper<Z>::Func
( DistMatrix<Z,VR,STAR,Int>& parent, Int i, Int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::SetReal");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T,typename Int>
inline void
DistMatrix<T,VR,STAR,Int>::SetImag( Int i, Int j, typename Base<T>::type alpha )
{ SetImagHelper<T>::Func( *this, i, j, alpha ); }

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,VR,STAR,Int>::SetImagHelper<Z>::Func
( DistMatrix<Z,VR,STAR,Int>& parent, Int i, Int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::SetImag");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T,typename Int>
inline void
DistMatrix<T,VR,STAR,Int>::UpdateReal
( Int i, Int j, typename Base<T>::type alpha )
{ UpdateRealHelper<T>::Func( *this, i, j, alpha ); }

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,VR,STAR,Int>::UpdateRealHelper<Z>::Func
( DistMatrix<Z,VR,STAR,Int>& parent, Int i, Int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::UpdateReal");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T,typename Int>
inline void
DistMatrix<T,VR,STAR,Int>::UpdateImag
( Int i, Int j, typename Base<T>::type alpha )
{ UpdateImagHelper<T>::Func( *this, i, j, alpha ); }

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,VR,STAR,Int>::UpdateImagHelper<Z>::Func
( DistMatrix<Z,VR,STAR,Int>& parent, Int i, Int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::UpdateImag");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}


template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,VR,STAR,Int>::SetToRandomHermitianHelper<Z>::Func
( DistMatrix<Z,VR,STAR,Int>& parent )
{   
#ifndef RELEASE
    PushCallStack("[VR,* ]::SetToRandomHermitian");
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
DistMatrix<T,VR,STAR,Int>::SetToRandomHermitianHelper<Complex<Z> >::Func
( DistMatrix<Complex<Z>,VR,STAR,Int>& parent )
{   
#ifndef RELEASE
    PushCallStack("[VR,* ]::SetToRandomHermitian");
    parent.AssertNotLockedView();
    if( parent.Height() != parent.Width() )
        throw std::logic_error("Hermitian matrices must be square");
#endif
    const Int width = parent.Width();
    const Int localHeight = parent.LocalHeight();
    const Int p = parent.Grid().Size();
    const Int colShift = parent.ColShift();

    parent.SetToRandom();

    Complex<Z>* thisLocalBuffer = parent.LocalBuffer();
    const Int thisLDim = parent.LocalLDim();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( Int iLocal=0; iLocal<localHeight; ++iLocal )
    {
        const Int i = colShift + iLocal*p;
        if( i < width )
        {
            const Z value = thisLocalBuffer[iLocal+i*thisLDim].real;
            thisLocalBuffer[iLocal+i*thisLDim] = value;
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,VR,STAR,Int>::SetToRandomHPDHelper<Z>::Func
( DistMatrix<Z,VR,STAR,Int>& parent )
{   
#ifndef RELEASE
    PushCallStack("[VR,* ]::SetToRandomHPD");
    parent.AssertNotLockedView();
    if( parent.Height() != parent.Width() )
        throw std::logic_error("Positive-definite matrices must be square");
#endif
    const Int width = parent.Width();
    const Int localHeight = parent.LocalHeight();
    const Int p = parent.Grid().Size();
    const Int colShift = parent.ColShift();

    parent.SetToRandom();

    Z* thisLocalBuffer = parent.LocalBuffer();
    const Int thisLDim = parent.LocalLDim();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( Int iLocal=0; iLocal<localHeight; ++iLocal )
    {
        const Int i = colShift + iLocal*p;
        if( i < width )
            thisLocalBuffer[iLocal+i*thisLDim] += width;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,VR,STAR,Int>::SetToRandomHPDHelper<Complex<Z> >::Func
( DistMatrix<Complex<Z>,VR,STAR,Int>& parent )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::SetToRandomHPD");
    parent.AssertNotLockedView();
    if( parent.Height() != parent.Width() )
        throw std::logic_error("Positive-definite matrices must be square");
#endif
    const Int width = parent.Width();
    const Int localHeight = parent.LocalHeight();
    const Int p = parent.Grid().Size();
    const Int colShift = parent.ColShift();

    parent.SetToRandom();

    Complex<Z>* thisLocalBuffer = parent.LocalBuffer();
    const Int thisLDim = parent.LocalLDim();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( Int iLocal=0; iLocal<localHeight; ++iLocal )
    {
        const Int i = colShift + iLocal*p;
        if( i < width )
        {
            const Z value = thisLocalBuffer[iLocal+i*thisLDim].real;
            thisLocalBuffer[iLocal+i*thisLDim] = value + width;
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename Z>
inline Z
DistMatrix<T,VR,STAR,Int>::GetRealHelper<Complex<Z> >::Func
( const DistMatrix<Complex<Z>,VR,STAR,Int>& parent, Int i, Int j ) 
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::GetReal");
    parent.AssertValidEntry( i, j );
#endif
    // We will determine the owner rank of entry (i,j) and broadcast from that
    // process over the entire g
    const elem::Grid& g = parent.Grid();
    const Int ownerRank = (i + parent.ColAlignment()) % g.Size();

    Z u;
    if( g.VRRank() == ownerRank )
    {
        const Int iLoc = (i-parent.ColShift()) / g.Size();
        u = parent.GetRealLocalEntry(iLoc,j);
    }
    mpi::Broadcast( &u, 1, ownerRank, g.VRComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T,typename Int>
template<typename Z>
inline Z
DistMatrix<T,VR,STAR,Int>::GetImagHelper<Complex<Z> >::Func
( const DistMatrix<Complex<Z>,VR,STAR,Int>& parent, Int i, Int j ) 
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::GetImag");
    parent.AssertValidEntry( i, j );
#endif
    // We will determine the owner rank of entry (i,j) and broadcast from that
    // process over the entire g
    const elem::Grid& g = parent.Grid();
    const Int ownerRank = (i + parent.ColAlignment()) % g.Size();

    Z u;
    if( g.VRRank() == ownerRank )
    {
        const Int iLoc = (i-parent.ColShift()) / g.Size();
        u = parent.GetImagLocalEntry(iLoc,j);
    }
    mpi::Broadcast( &u, 1, ownerRank, g.VRComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,VR,STAR,Int>::SetRealHelper<Complex<Z> >::Func
( DistMatrix<Complex<Z>,VR,STAR,Int>& parent, Int i, Int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::SetReal");
    parent.AssertValidEntry( i, j );
#endif
    const elem::Grid& g = parent.Grid();
    const Int ownerRank = (i + parent.ColAlignment()) % g.Size();

    if( g.VRRank() == ownerRank )
    {
        const Int iLoc = (i-parent.ColShift()) / g.Size();
        parent.SetRealLocalEntry(iLoc,j,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,VR,STAR,Int>::SetImagHelper<Complex<Z> >::Func
( DistMatrix<Complex<Z>,VR,STAR,Int>& parent, Int i, Int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::SetImag");
    parent.AssertValidEntry( i, j );
#endif
    const elem::Grid& g = parent.Grid();
    const Int ownerRank = (i + parent.ColAlignment()) % g.Size();

    if( g.VRRank() == ownerRank )
    {
        const Int iLoc = (i-parent.ColShift()) / g.Size();
        parent.SetImagLocalEntry(iLoc,j,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,VR,STAR,Int>::UpdateRealHelper<Complex<Z> >::Func
( DistMatrix<Complex<Z>,VR,STAR,Int>& parent, Int i, Int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::UpdateReal");
    parent.AssertValidEntry( i, j );
#endif
    const elem::Grid& g = parent.Grid();
    const Int ownerRank = (i + parent.ColAlignment()) % g.Size();

    if( g.VRRank() == ownerRank )
    {
        const Int iLoc = (i-parent.ColShift()) / g.Size();
        parent.UpdateRealLocalEntry(iLoc,j,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,VR,STAR,Int>::UpdateImagHelper<Complex<Z> >::Func
( DistMatrix<Complex<Z>,VR,STAR,Int>& parent, Int i, Int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[VR,* ]::UpdateImag");
    parent.AssertValidEntry( i, j );
#endif
    const elem::Grid& g = parent.Grid();
    const Int ownerRank = (i + parent.ColAlignment()) % g.Size();

    if( g.VRRank() == ownerRank )
    {
        const Int iLoc = (i-parent.ColShift()) / g.Size();
        parent.UpdateImagLocalEntry(iLoc,j,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
