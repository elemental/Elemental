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
DistMatrix<T,STAR,MR,Int>::SetToRandomHermitian()
{ SetToRandomHermitianHelper<T>::Func( *this ); }

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,MR,Int>::SetToRandomHPD()
{ SetToRandomHPDHelper<T>::Func( *this ); }

template<typename T,typename Int>
inline typename Base<T>::type
DistMatrix<T,STAR,MR,Int>::GetReal( Int i, Int j ) const
{ return GetRealHelper<T>::Func( *this, i, j ); }

template<typename T,typename Int>
template<typename Z>
inline Z
DistMatrix<T,STAR,MR,Int>::GetRealHelper<Z>::Func
( const DistMatrix<Z,STAR,MR,Int>& parent, Int i, Int j )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::GetRealHelper");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T,typename Int>
inline typename Base<T>::type
DistMatrix<T,STAR,MR,Int>::GetImag( Int i, Int j ) const
{ return GetImagHelper<T>::Func( *this, i, j ); }

template<typename T,typename Int>
template<typename Z>
inline Z
DistMatrix<T,STAR,MR,Int>::GetImagHelper<Z>::Func
( const DistMatrix<Z,STAR,MR,Int>& parent, Int i, Int j )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::GetImag");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,MR,Int>::SetReal
( Int i, Int j, typename Base<T>::type alpha )
{ SetRealHelper<T>::Func( *this, i, j, alpha ); }

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,STAR,MR,Int>::SetRealHelper<Z>::Func
( DistMatrix<Z,STAR,MR,Int>& parent, Int i, Int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::SetReal");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,MR,Int>::SetImag
( Int i, Int j, typename Base<T>::type alpha )
{ SetImagHelper<T>::Func( *this, i, j, alpha ); }

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,STAR,MR,Int>::SetImagHelper<Z>::Func
( DistMatrix<Z,STAR,MR,Int>& parent, Int i, Int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::SetImag");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,MR,Int>::UpdateReal
( Int i, Int j, typename Base<T>::type alpha )
{ UpdateRealHelper<T>::Func( *this, i, j, alpha ); }

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,STAR,MR,Int>::UpdateRealHelper<Z>::Func
( DistMatrix<Z,STAR,MR,Int>& parent, Int i, Int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::UpdateReal");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,MR,Int>::UpdateImag
( Int i, Int j, typename Base<T>::type alpha )
{ UpdateImagHelper<T>::Func( *this, i, j, alpha ); }

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,STAR,MR,Int>::UpdateImagHelper<Z>::Func
( DistMatrix<Z,STAR,MR,Int>& parent, Int i, Int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::UpdateImag");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}


template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,STAR,MR,Int>::SetToRandomHermitianHelper<Z>::Func
( DistMatrix<Z,STAR,MR,Int>& parent )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::SetToRandomHermitian");
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
DistMatrix<T,STAR,MR,Int>::SetToRandomHermitianHelper<Complex<Z> >::Func
( DistMatrix<Complex<Z>,STAR,MR,Int>& parent )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::SetToRandomHermitian");
    parent.AssertNotLockedView();
    if( parent.Height() != parent.Width() )
        throw std::logic_error("Hermitian matrices must be square");
#endif
    const Int height     = parent.Height();
    const Int localWidth = parent.LocalWidth();
    const Int c          = parent.Grid().Width();
    const Int rowShift   = parent.RowShift();

    parent.SetToRandom();

    Complex<Z>* thisLocalBuffer = parent.LocalBuffer();
    const Int thisLDim = parent.LocalLDim();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( Int jLocal=0; jLocal<localWidth; ++jLocal )
    {
        const Int j = rowShift + jLocal*c;
        if( j < height )
        {
            const Z value = thisLocalBuffer[j+jLocal*thisLDim].real;
            thisLocalBuffer[j+jLocal*thisLDim] = value;
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,STAR,MR,Int>::SetToRandomHPDHelper<Z>::Func
( DistMatrix<Z,STAR,MR,Int>& parent )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::SetToRandomHPD");
    parent.AssertNotLockedView();
    if( parent.Height() != parent.Width() )
        throw std::logic_error("Positive-definite matrices must be square");
#endif
    const Int height     = parent.Height();
    const Int width      = parent.Width();
    const Int localWidth = parent.LocalWidth();
    const Int c          = parent.Grid().Width();
    const Int rowShift   = parent.RowShift();

    parent.SetToRandom();

    Z* thisLocalBuffer = parent.LocalBuffer();
    const Int thisLDim = parent.LocalLDim();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( Int jLocal=0; jLocal<localWidth; ++jLocal )
    {
        const Int j = rowShift + jLocal*c;
        if( j < height )
            thisLocalBuffer[j+jLocal*thisLDim] += width;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,STAR,MR,Int>::SetToRandomHPDHelper<Complex<Z> >::Func
( DistMatrix<Complex<Z>,STAR,MR,Int>& parent )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::SetToRandomHPD");
    parent.AssertNotLockedView();
    if( parent.Height() != parent.Width() )
        throw std::logic_error("Positive-definite matrices must be square");
#endif
    const Int height     = parent.Height();
    const Int width      = parent.Width();
    const Int localWidth = parent.LocalWidth();
    const Int c          = parent.Grid().Width();
    const Int rowShift   = parent.RowShift();

    parent.SetToRandom();

    Complex<Z>* thisLocalBuffer = parent.LocalBuffer();
    const Int thisLDim = parent.LocalLDim();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( Int jLocal=0; jLocal<localWidth; ++jLocal )
    {
        const Int j = rowShift + jLocal*c;
        if( j < height )
        {
            const Z value = thisLocalBuffer[j+jLocal*thisLDim].real;
            thisLocalBuffer[j+jLocal*thisLDim] = value + width;
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename Z>
inline Z
DistMatrix<T,STAR,MR,Int>::GetRealHelper<Complex<Z> >::Func
( const DistMatrix<Complex<Z>,STAR,MR,Int>& parent, Int i, Int j ) 
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::GetReal");
    parent.AssertValidEntry( i, j );
#endif
    // We will determine the owner column of entry (i,j) and broadcast from that
    // column within each process row
    const elem::Grid& g = parent.Grid();
    const Int ownerCol = (j + parent.RowAlignment()) % g.Width();

    Z u;
    if( g.MRRank() == ownerCol )
    {
        const Int jLoc = (j-parent.RowShift()) / g.Width();
        u = parent.GetRealLocalEntry(i,jLoc);
    }
    mpi::Broadcast( &u, 1, ownerCol, g.MRComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T,typename Int>
template<typename Z>
inline Z
DistMatrix<T,STAR,MR,Int>::GetImagHelper<Complex<Z> >::Func
( const DistMatrix<Complex<Z>,STAR,MR,Int>& parent, Int i, Int j ) 
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::GetImag");
    parent.AssertValidEntry( i, j );
#endif
    // We will determine the owner column of entry (i,j) and broadcast from that
    // column within each process row
    const elem::Grid& g = parent.Grid();
    const Int ownerCol = (j + parent.RowAlignment()) % g.Width();

    Z u;
    if( g.MRRank() == ownerCol )
    {
        const Int jLoc = (j-parent.RowShift()) / g.Width();
        u = parent.GetImagLocalEntry(i,jLoc);
    }
    mpi::Broadcast( &u, 1, ownerCol, g.MRComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,STAR,MR,Int>::SetRealHelper<Complex<Z> >::Func
( DistMatrix<Complex<Z>,STAR,MR,Int>& parent, Int i, Int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::SetReal");
    parent.AssertValidEntry( i, j );
#endif
    const elem::Grid& g = parent.Grid();
    const Int ownerCol = (j + parent.RowAlignment()) % g.Width();

    if( g.MRRank() == ownerCol )
    {
        const Int jLoc = (j-parent.RowShift()) / g.Width();
        parent.SetRealLocalEntry(i,jLoc,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,STAR,MR,Int>::SetImagHelper<Complex<Z> >::Func
( DistMatrix<Complex<Z>,STAR,MR,Int>& parent, Int i, Int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::SetImag");
    parent.AssertValidEntry( i, j );
#endif
    const elem::Grid& g = parent.Grid();
    const Int ownerCol = (j + parent.RowAlignment()) % g.Width();

    if( g.MRRank() == ownerCol )
    {
        const Int jLoc = (j-parent.RowShift()) / g.Width();
        parent.SetImagLocalEntry(i,jLoc,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,STAR,MR,Int>::UpdateRealHelper<Complex<Z> >::Func
( DistMatrix<Complex<Z>,STAR,MR,Int>& parent, Int i, Int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::UpdateReal");
    parent.AssertValidEntry( i, j );
#endif
    const elem::Grid& g = parent.Grid();
    const Int ownerCol = (j + parent.RowAlignment()) % g.Width();

    if( g.MRRank() == ownerCol )
    {
        const Int jLoc = (j-parent.RowShift()) / g.Width();
        parent.UpdateRealLocalEntry(i,jLoc,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,STAR,MR,Int>::UpdateImagHelper<Complex<Z> >::Func
( DistMatrix<Complex<Z>,STAR,MR,Int>& parent, Int i, Int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::UpdateImag");
    parent.AssertValidEntry( i, j );
#endif
    const elem::Grid& g = parent.Grid();
    const Int ownerCol = (j + parent.RowAlignment()) % g.Width();

    if( g.MRRank() == ownerCol )
    {
        const Int jLoc = (j-parent.RowShift()) / g.Width();
        parent.UpdateImagLocalEntry(i,jLoc,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
