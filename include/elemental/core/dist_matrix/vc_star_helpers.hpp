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
DistMatrix<T,VC,STAR,Int>::SetToRandomHermitian()
{ SetToRandomHermitianHelper<T>::Func( *this ); }

template<typename T,typename Int>
inline void
DistMatrix<T,VC,STAR,Int>::SetToRandomHPD()
{ SetToRandomHPDHelper<T>::Func( *this ); }

template<typename T,typename Int>
inline typename Base<T>::type
DistMatrix<T,VC,STAR,Int>::GetReal( Int i, Int j ) const
{ return GetRealHelper<T>::Func( *this, i, j ); }

template<typename T,typename Int>
template<typename Z>
inline Z
DistMatrix<T,VC,STAR,Int>::GetRealHelper<Z>::Func
( const DistMatrix<Z,VC,STAR,Int>& parent, Int i, Int j )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::GetRealHelper");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T,typename Int>
inline typename Base<T>::type
DistMatrix<T,VC,STAR,Int>::GetImag( Int i, Int j ) const
{ return GetImagHelper<T>::Func( *this, i, j ); }

template<typename T,typename Int>
template<typename Z>
inline Z
DistMatrix<T,VC,STAR,Int>::GetImagHelper<Z>::Func
( const DistMatrix<Z,VC,STAR,Int>& parent, Int i, Int j )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::GetImag");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T,typename Int>
inline void
DistMatrix<T,VC,STAR,Int>::SetReal( Int i, Int j, typename Base<T>::type alpha )
{ SetRealHelper<T>::Func( *this, i, j, alpha ); }

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,VC,STAR,Int>::SetRealHelper<Z>::Func
( DistMatrix<Z,VC,STAR,Int>& parent, Int i, Int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::SetReal");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T,typename Int>
inline void
DistMatrix<T,VC,STAR,Int>::SetImag( Int i, Int j, typename Base<T>::type alpha )
{ SetImagHelper<T>::Func( *this, i, j, alpha ); }

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,VC,STAR,Int>::SetImagHelper<Z>::Func
( DistMatrix<Z,VC,STAR,Int>& parent, Int i, Int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::SetImag");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T,typename Int>
inline void
DistMatrix<T,VC,STAR,Int>::UpdateReal
( Int i, Int j, typename Base<T>::type alpha )
{ UpdateRealHelper<T>::Func( *this, i, j, alpha ); }

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,VC,STAR,Int>::UpdateRealHelper<Z>::Func
( DistMatrix<Z,VC,STAR,Int>& parent, Int i, Int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::UpdateReal");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T,typename Int>
inline void
DistMatrix<T,VC,STAR,Int>::UpdateImag
( Int i, Int j, typename Base<T>::type alpha )
{ UpdateImagHelper<T>::Func( *this, i, j, alpha ); }

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,VC,STAR,Int>::UpdateImagHelper<Z>::Func
( DistMatrix<Z,VC,STAR,Int>& parent, Int i, Int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::UpdateImag");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}


template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,VC,STAR,Int>::SetToRandomHermitianHelper<Z>::Func
( DistMatrix<Z,VC,STAR,Int>& parent )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::SetToRandomHermitian");
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
DistMatrix<T,VC,STAR,Int>::SetToRandomHermitianHelper<Complex<Z> >::Func
( DistMatrix<Complex<Z>,VC,STAR,Int>& parent )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::SetToRandomHermitian");
    parent.AssertNotLockedView();
    if( parent.Height() != parent.Width() )
        throw std::logic_error("Hermitian matrices must be square");
#endif
    const Int width       = parent.Width();
    const Int localHeight = parent.LocalHeight();
    const Int p           = parent.Grid().Size();
    const Int colShift    = parent.ColShift();

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
DistMatrix<T,VC,STAR,Int>::SetToRandomHPDHelper<Z>::Func
( DistMatrix<Z,VC,STAR,Int>& parent )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::SetToRandomHPD");
    parent.AssertNotLockedView();
    if( parent.Height() != parent.Width() )
        throw std::logic_error("Positive-definite matrices must be square");
#endif
    const Int width       = parent.Width();
    const Int localHeight = parent.LocalHeight();
    const Int p           = parent.Grid().Size();
    const Int colShift    = parent.ColShift();

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
DistMatrix<T,VC,STAR,Int>::SetToRandomHPDHelper<Complex<Z> >::Func
( DistMatrix<Complex<Z>,VC,STAR,Int>& parent )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::SetToRandomHPD");
    parent.AssertNotLockedView();
    if( parent.Height() != parent.Width() )
        throw std::logic_error("Positive-definite matrices must be square");
#endif
    const Int width       = parent.Width();
    const Int localHeight = parent.LocalHeight();
    const Int p           = parent.Grid().Size();
    const Int colShift    = parent.ColShift();

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
DistMatrix<T,VC,STAR,Int>::GetRealHelper<Complex<Z> >::Func
( const DistMatrix<Complex<Z>,VC,STAR,Int>& parent, Int i, Int j ) 
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::GetReal");
    parent.AssertValidEntry( i, j );
#endif
    // We will determine the owner rank of entry (i,j) and broadcast from that
    // process over the entire g
    const elem::Grid& g = parent.Grid();
    const Int ownerRank = (i + parent.ColAlignment()) % g.Size();

    Z u;
    if( g.VCRank() == ownerRank )
    {
        const Int iLoc = (i-parent.ColShift()) / g.Size();
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
DistMatrix<T,VC,STAR,Int>::GetImagHelper<Complex<Z> >::Func
( const DistMatrix<Complex<Z>,VC,STAR,Int>& parent, Int i, Int j ) 
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::GetImag");
    parent.AssertValidEntry( i, j );
#endif
    // We will determine the owner rank of entry (i,j) and broadcast from that
    // process over the entire g
    const elem::Grid& g = parent.Grid();
    const Int ownerRank = (i + parent.ColAlignment()) % g.Size();

    Z u;
    if( g.VCRank() == ownerRank )
    {
        const Int iLoc = (i-parent.ColShift()) / g.Size();
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
DistMatrix<T,VC,STAR,Int>::SetRealHelper<Complex<Z> >::Func
( DistMatrix<Complex<Z>,VC,STAR,Int>& parent, Int i, Int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::SetReal");
    parent.AssertValidEntry( i, j );
#endif
    const elem::Grid& g = parent.Grid();
    const Int ownerRank = (i + parent.ColAlignment()) % g.Size();

    if( g.VCRank() == ownerRank )
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
DistMatrix<T,VC,STAR,Int>::SetImagHelper<Complex<Z> >::Func
( DistMatrix<Complex<Z>,VC,STAR,Int>& parent, Int i, Int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::SetImag");
    parent.AssertValidEntry( i, j );
#endif
    const elem::Grid& g = parent.Grid();
    const Int ownerRank = (i + parent.ColAlignment()) % g.Size();

    if( g.VCRank() == ownerRank )
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
DistMatrix<T,VC,STAR,Int>::UpdateRealHelper<Complex<Z> >::Func
( DistMatrix<Complex<Z>,VC,STAR,Int>& parent, Int i, Int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::UpdateReal");
    parent.AssertValidEntry( i, j );
#endif
    const elem::Grid& g = parent.Grid();
    const Int ownerRank = (i + parent.ColAlignment()) % g.Size();

    if( g.VCRank() == ownerRank )
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
void
DistMatrix<T,VC,STAR,Int>::UpdateImagHelper<Complex<Z> >::Func
( DistMatrix<Complex<Z>,VC,STAR,Int>& parent, Int i, Int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::UpdateImag");
    parent.AssertValidEntry( i, j );
#endif
    const elem::Grid& g = parent.Grid();
    const Int ownerRank = (i + parent.ColAlignment()) % g.Size();

    if( g.VCRank() == ownerRank )
    {
        const Int iLoc = (i-parent.ColShift()) / g.Size();
        parent.UpdateImagLocalEntry(iLoc,j,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,VC,STAR,Int>::GetRealDiagonalHelper<Complex<Z> >::Func
( const DistMatrix<Complex<Z>,VC,STAR,Int>& parent,
        DistMatrix<Z,VC,STAR,Int>& d, Int offset )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::GetRealDiagonal");
    if( d.Viewing() )
        parent.AssertSameGrid( d );
#endif
    const Int length = parent.DiagonalLength(offset);
#ifndef RELEASE
    if( d.Viewing() && (length != d.Height() || d.Width() != 1) )
    {
        std::ostringstream msg;
        msg << "d is not a column vec of the same length as the diagonal:\n"
            << "  A ~ " << parent.Height() << " x " << parent.Width() << "\n"
            << "  d ~ " << d.Height() << " x " << d.Width() << "\n"
            << "  A diag length: " << length << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( ( d.Viewing() || d.ConstrainedColAlignment() ) &&
        !d.AlignedWithDiagonal( parent, offset ) )
        throw std::logic_error("d must be aligned with the offset diag");
#endif
    const elem::Grid& g = parent.Grid();
    if( !d.Viewing() )
    {
        d.SetGrid( g );
        if( !d.ConstrainedColAlignment() )
            d.AlignWithDiagonal( parent, offset );
        d.ResizeTo( length, 1 );
    }

    if( g.InGrid() )
    {
        const Int p = g.Size();
        const Int colShift = parent.ColShift();
        const Int diagShift = d.ColShift();

        Int iStart, jStart;
        if( offset >= 0 )
        {
            iStart = diagShift;
            jStart = diagShift+offset;
        }
        else
        {
            iStart = diagShift-offset;
            jStart = diagShift;
        }

        const Int iLocalStart = (iStart-colShift) / p;
        const Int localDiagLength = d.LocalHeight();
        const Complex<Z>* thisLocalBuffer = parent.LockedLocalBuffer();
        const Int thisLDim = parent.LocalLDim();
        Z* dLocalBuffer = d.LocalBuffer();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( Int k=0; k<localDiagLength; ++k )
        {
            const Int iLocal = iLocalStart+k;
            const Int jLocal = jStart+k*p;
            dLocalBuffer[k] = thisLocalBuffer[iLocal+jLocal*thisLDim].real;
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,VC,STAR,Int>::GetImagDiagonalHelper<Complex<Z> >::Func
( const DistMatrix<Complex<Z>,VC,STAR,Int>& parent,
        DistMatrix<Z,VC,STAR,Int>& d, Int offset )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::GetImagDiagonal");
    if( d.Viewing() )
        parent.AssertSameGrid( d );
#endif
    const Int length = parent.DiagonalLength(offset);
#ifndef RELEASE
    if( d.Viewing() && length != d.Height() )
    {
        std::ostringstream msg;
        msg << "d is not of the same length as the diagonal:\n"
            << "  A ~ " << parent.Height() << " x " << parent.Width() << "\n"
            << "  d ~ " << d.Height() << " x " << d.Width() << "\n"
            << "  A diag length: " << length << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( ( d.Viewing() || d.ConstrainedColAlignment() ) &&
        !d.AlignedWithDiagonal( parent, offset ) )
        throw std::logic_error("d must be aligned with the offset diag");
#endif
    const elem::Grid& g = parent.Grid();
    if( !d.Viewing() )
    {
        d.SetGrid( g );
        if( !d.ConstrainedColAlignment() )
            d.AlignWithDiagonal( parent, offset );
        d.ResizeTo( length, 1 );
    }

    if( g.InGrid() )
    {
        const Int p = g.Size();
        const Int colShift = parent.ColShift();
        const Int diagShift = d.ColShift();

        Int iStart, jStart;
        if( offset >= 0 )
        {
            iStart = diagShift;
            jStart = diagShift+offset;
        }
        else
        {
            iStart = diagShift-offset;
            jStart = diagShift;
        }

        const Int iLocalStart = (iStart-colShift) / p;
        const Int localDiagLength = d.LocalHeight();
        const Complex<Z>* thisLocalBuffer = parent.LockedLocalBuffer();
        const Int thisLDim = parent.LocalLDim();
        Z* dLocalBuffer = d.LocalBuffer();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( Int k=0; k<localDiagLength; ++k )
        {
            const Int iLocal = iLocalStart+k;
            const Int jLocal = jStart+k*p;
            dLocalBuffer[k] = thisLocalBuffer[iLocal+jLocal*thisLDim].imag;
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,VC,STAR,Int>::GetRealDiagonalHelper<Complex<Z> >::Func
( const DistMatrix<Complex<Z>,VC,STAR,Int>& parent,
        DistMatrix<Z,STAR,VC,Int>& d, Int offset )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::GetRealDiagonal");
    if( d.Viewing() )
        parent.AssertSameGrid( d );
#endif
    const Int length = parent.DiagonalLength(offset);
#ifndef RELEASE
    if( d.Viewing() && length != d.Width() )
    {
        std::ostringstream msg;
        msg << "d is not of the same length as the diagonal:\n"
            << "  A ~ " << parent.Height() << " x " << parent.Width() << "\n"
            << "  d ~ " << d.Height() << " x " << d.Width() << "\n"
            << "  A diag length: " << length << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( ( d.Viewing() || d.ConstrainedRowAlignment() ) &&
        !d.AlignedWithDiagonal( parent, offset ) )
        throw std::logic_error("d must be aligned with the offset diag");
#endif
    const elem::Grid& g = parent.Grid();
    if( !d.Viewing() )
    {
        d.SetGrid( g );
        if( !d.ConstrainedRowAlignment() )
            d.AlignWithDiagonal( parent, offset );
        d.ResizeTo( 1, length );
    }

    if( g.InGrid() )
    {
        const Int p = g.Size();
        const Int colShift = parent.ColShift();
        const Int diagShift = d.RowShift();

        Int iStart, jStart;
        if( offset >= 0 )
        {
            iStart = diagShift;
            jStart = diagShift+offset;
        }
        else
        {
            iStart = diagShift-offset;
            jStart = diagShift;
        }

        const Int iLocalStart = (iStart-colShift) / p;
        const Int localDiagLength = d.LocalWidth();
        const Complex<Z>* thisLocalBuffer = parent.LockedLocalBuffer();
        const Int thisLDim = parent.LocalLDim();
        Z* dLocalBuffer = d.LocalBuffer();
        const Int dLDim = d.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( Int k=0; k<localDiagLength; ++k )
        {
            const Int iLocal = iLocalStart+k;
            const Int jLocal = jStart+k*p;
            dLocalBuffer[k*dLDim] =thisLocalBuffer[iLocal+jLocal*thisLDim].real;
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,VC,STAR,Int>::GetImagDiagonalHelper<Complex<Z> >::Func
( const DistMatrix<Complex<Z>,VC,STAR,Int>& parent,
        DistMatrix<Z,STAR,VC,Int>& d, Int offset )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::GetImagDiagonal");
    if( d.Viewing() )
        parent.AssertSameGrid( d );
#endif
    const Int length = parent.DiagonalLength(offset);
#ifndef RELEASE
    if( d.Viewing() && length != d.Width() )
    {
        std::ostringstream msg;
        msg << "d is not of the same length as the diagonal:\n"
            << "  A ~ " << parent.Height() << " x " << parent.Width() << "\n"
            << "  d ~ " << d.Height() << " x " << d.Width() << "\n"
            << "  A diag length: " << length << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( ( d.Viewing() || d.ConstrainedRowAlignment() ) &&
        !d.AlignedWithDiagonal( parent, offset ) )
        throw std::logic_error("d must be aligned with the offset diag");
#endif
    const elem::Grid& g = parent.Grid();
    if( !d.Viewing() )
    {
        d.SetGrid( g );
        if( !d.ConstrainedRowAlignment() )
            d.AlignWithDiagonal( parent, offset );
        d.ResizeTo( 1, length );
    }

    if( g.InGrid() )
    {
        const Int p = g.Size();
        const Int colShift = parent.ColShift();
        const Int diagShift = d.RowShift();

        Int iStart, jStart;
        if( offset >= 0 )
        {
            iStart = diagShift;
            jStart = diagShift+offset;
        }
        else
        {
            iStart = diagShift-offset;
            jStart = diagShift;
        }

        const Int iLocalStart = (iStart-colShift) / p;
        const Int localDiagLength = d.LocalWidth();
        const Complex<Z>* thisLocalBuffer = parent.LockedLocalBuffer();
        const Int thisLDim = parent.LocalLDim();
        Z* dLocalBuffer = d.LocalBuffer();
        const Int dLDim = d.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( Int k=0; k<localDiagLength; ++k )
        {
            const Int iLocal = iLocalStart+k;
            const Int jLocal = jStart+k*p;
            dLocalBuffer[k*dLDim] =thisLocalBuffer[iLocal+jLocal*thisLDim].imag;
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,VC,STAR,Int>::SetRealDiagonalHelper<Complex<Z> >::Func
(       DistMatrix<Complex<Z>,VC,STAR,Int>& parent,
  const DistMatrix<Z,VC,STAR,Int>& d, Int offset )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::SetRealDiagonal");
    parent.AssertSameGrid( d );
    if( d.Width() != 1 )
        throw std::logic_error("d must be a column vector");
    const Int length = parent.DiagonalLength(offset);
    if( length != d.Height() )
    {
        std::ostringstream msg;
        msg << "d is not of the same length as the diagonal:\n"
            << "  A ~ " << parent.Height() << " x " << parent.Width() << "\n"
            << "  d ~ " << d.Height() << " x " << d.Width() << "\n"
            << "  A diag length: " << length << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( !d.AlignedWithDiagonal( parent, offset ) )
        throw std::logic_error("d must be aligned with the 'offset' diagonal");
#endif
    const elem::Grid& g = parent.Grid();
    if( g.InGrid() )
    {
        const Int p = g.Size();
        const Int colShift = parent.ColShift();
        const Int diagShift = d.ColShift();

        Int iStart, jStart;
        if( offset >= 0 )
        {
            iStart = diagShift;
            jStart = diagShift+offset;
        }
        else
        {
            iStart = diagShift-offset;
            jStart = diagShift;
        }

        const Int iLocalStart = (iStart-colShift)/p;
        const Int localDiagLength = d.LocalHeight();

        const Z* dLocalBuffer = d.LockedLocalBuffer();
        Complex<Z>* thisLocalBuffer = parent.LocalBuffer();
        const Int thisLDim = parent.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( Int k=0; k<localDiagLength; ++k )
        {
            const Int iLocal = iLocalStart+k;
            const Int jLocal = jStart+k*p;
            const Z u = dLocalBuffer[k];
            const Z v = thisLocalBuffer[iLocal+jLocal*thisLDim].imag;
            thisLocalBuffer[iLocal+jLocal*thisLDim] = Complex<Z>(u,v);
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,VC,STAR,Int>::SetImagDiagonalHelper<Complex<Z> >::Func
(       DistMatrix<Complex<Z>,VC,STAR,Int>& parent,
  const DistMatrix<Z,VC,STAR,Int>& d, Int offset )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::SetImagDiagonal");
    parent.AssertSameGrid( d );
    if( d.Width() != 1 )
        throw std::logic_error("d must be a column vector");
    const Int length = parent.DiagonalLength(offset);
    if( length != d.Height() )
    {
        std::ostringstream msg;
        msg << "d is not of the same length as the diagonal:\n"
            << "  A ~ " << parent.Height() << " x " << parent.Width() << "\n"
            << "  d ~ " << d.Height() << " x " << d.Width() << "\n"
            << "  A diag length: " << length << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( !d.AlignedWithDiagonal( parent, offset ) )
        throw std::logic_error("d must be aligned with the 'offset' diagonal");
#endif
    const elem::Grid& g = parent.Grid();
    if( g.InGrid() )
    {
        const Int p = g.Size();
        const Int colShift = parent.ColShift();
        const Int diagShift = d.ColShift();

        Int iStart, jStart;
        if( offset >= 0 )
        {
            iStart = diagShift;
            jStart = diagShift+offset;
        }
        else
        {
            iStart = diagShift-offset;
            jStart = diagShift;
        }

        const Int iLocalStart = (iStart-colShift)/p;
        const Int localDiagLength = d.LocalHeight();

        const Z* dLocalBuffer = d.LockedLocalBuffer();
        Complex<Z>* thisLocalBuffer = parent.LocalBuffer();
        const Int thisLDim = parent.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( Int k=0; k<localDiagLength; ++k )
        {
            const Int iLocal = iLocalStart+k;
            const Int jLocal = jStart+k*p;
            const Z u = thisLocalBuffer[iLocal+jLocal*thisLDim].real;
            const Z v = dLocalBuffer[k];
            thisLocalBuffer[iLocal+jLocal*thisLDim] = Complex<Z>(u,v);
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,VC,STAR,Int>::SetRealDiagonalHelper<Complex<Z> >::Func
(       DistMatrix<Complex<Z>,VC,STAR,Int>& parent,
  const DistMatrix<Z,STAR,VC,Int>& d, Int offset )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::SetRealDiagonal");
    parent.AssertSameGrid( d );
    if( d.Height() != 1 )
        throw std::logic_error("d must be a row vector");
    const Int length = parent.DiagonalLength(offset);
    if( length != d.Width() )
    {
        std::ostringstream msg;
        msg << "d is not of the same length as the diagonal:\n"
            << "  A ~ " << parent.Height() << " x " << parent.Width() << "\n"
            << "  d ~ " << d.Height() << " x " << d.Width() << "\n"
            << "  A diag length: " << length << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( !d.AlignedWithDiagonal( parent, offset ) )
        throw std::logic_error("d must be aligned with the 'offset' diagonal");
#endif
    const elem::Grid& g = parent.Grid();
    if( g.InGrid() )
    {
        const Int p = g.Size();
        const Int colShift = parent.ColShift();
        const Int diagShift = d.RowShift();

        Int iStart, jStart;
        if( offset >= 0 )
        {
            iStart = diagShift;
            jStart = diagShift+offset;
        }
        else
        {
            iStart = diagShift-offset;
            jStart = diagShift;
        }

        const Int iLocalStart = (iStart-colShift)/p;
        const Int localDiagLength = d.LocalWidth();

        const Z* dLocalBuffer = d.LockedLocalBuffer();
        Complex<Z>* thisLocalBuffer = parent.LocalBuffer();
        const Int dLDim = d.LocalLDim();
        const Int thisLDim = parent.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( Int k=0; k<localDiagLength; ++k )
        {
            const Int iLocal = iLocalStart+k;
            const Int jLocal = jStart+k*p;
            const Z u = dLocalBuffer[k*dLDim];
            const Z v = thisLocalBuffer[iLocal+jLocal*thisLDim].imag;
            thisLocalBuffer[iLocal+jLocal*thisLDim] = Complex<Z>(u,v);
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,VC,STAR,Int>::SetImagDiagonalHelper<Complex<Z> >::Func
(       DistMatrix<Complex<Z>,VC,STAR,Int>& parent,
  const DistMatrix<Z,STAR,VC,Int>& d, Int offset )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::SetImagDiagonal");
    parent.AssertSameGrid( d );
    if( d.Height() != 1 )
        throw std::logic_error("d must be a row vector");
    const Int length = parent.DiagonalLength(offset);
    if( length != d.Width() )
    {
        std::ostringstream msg;
        msg << "d is not of the same length as the diagonal:\n"
            << "  A ~ " << parent.Height() << " x " << parent.Width() << "\n"
            << "  d ~ " << d.Height() << " x " << d.Width() << "\n"
            << "  A diag length: " << length << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( !d.AlignedWithDiagonal( parent, offset ) )
        throw std::logic_error("d must be aligned with the 'offset' diagonal");
#endif
    const elem::Grid& g = parent.Grid();
    if( g.InGrid() )
    {
        const Int p = g.Size();
        const Int colShift = parent.ColShift();
        const Int diagShift = d.RowShift();

        Int iStart, jStart;
        if( offset >= 0 )
        {
            iStart = diagShift;
            jStart = diagShift+offset;
        }
        else
        {
            iStart = diagShift-offset;
            jStart = diagShift;
        }

        const Int iLocalStart = (iStart-colShift)/p;
        const Int localDiagLength = d.LocalWidth();

        const Z* dLocalBuffer = d.LockedLocalBuffer();
        Complex<Z>* thisLocalBuffer = parent.LocalBuffer();
        const Int dLDim = d.LocalLDim();
        const Int thisLDim = parent.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( Int k=0; k<localDiagLength; ++k )
        {
            const Int iLocal = iLocalStart+k;
            const Int jLocal = jStart+k*p;
            const Z u = thisLocalBuffer[iLocal+jLocal*thisLDim].real;
            const Z v = dLocalBuffer[k*dLDim];
            thisLocalBuffer[iLocal+jLocal*thisLDim] = Complex<Z>(u,v);
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
