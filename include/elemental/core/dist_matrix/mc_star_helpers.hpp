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
DistMatrix<T,MC,STAR,Int>::SetToRandomHermitian()
{ SetToRandomHermitianHelper<T>::Func( *this ); }

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,MC,STAR,Int>::SetToRandomHermitianHelper<Z>::Func
( DistMatrix<Z,MC,STAR,Int>& parent )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::SetToRandomHermitian");
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
DistMatrix<T,MC,STAR,Int>::SetToRandomHermitianHelper<std::complex<Z> >::Func
( DistMatrix<std::complex<Z>,MC,STAR,Int>& parent )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::SetToRandomHermitian");
    parent.AssertNotLockedView();
    if( parent.Height() != parent.Width() )
        throw std::logic_error("Hermitian matrices must be square");
#endif
    const Int width = parent.Width();
    const Int localHeight = parent.LocalHeight();
    const Int r = parent.Grid().Height();
    const Int colShift = parent.ColShift();

    parent.SetToRandom();

    std::complex<Z>* thisLocalBuffer = parent.LocalBuffer();
    const Int thisLDim = parent.LocalLDim();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( Int iLocal=0; iLocal<localHeight; ++iLocal )
    {
        const Int i = colShift + iLocal*r;
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

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::SetToRandomHPD()
{ SetToRandomHPDHelper<T>::Func( *this ); }

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,MC,STAR,Int>::SetToRandomHPDHelper<Z>::Func
( DistMatrix<Z,MC,STAR,Int>& parent )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::SetToRandomHPD");
    parent.AssertNotLockedView();
    if( parent.Height() != parent.Width() )
        throw std::logic_error("Positive-definite matrices must be square");
#endif
    const Int width = parent.Width();
    const Int localHeight = parent.LocalHeight();
    const Int r = parent.Grid().Height();
    const Int colShift = parent.ColShift();

    parent.SetToRandom();

    Z* thisLocalBuffer = parent.LocalBuffer();
    const Int thisLDim = parent.LocalLDim();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( Int iLocal=0; iLocal<localHeight; ++iLocal )
    {
        const Int i = colShift + iLocal*r;
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
DistMatrix<T,MC,STAR,Int>::SetToRandomHPDHelper<std::complex<Z> >::Func
( DistMatrix<std::complex<Z>,MC,STAR,Int>& parent )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::SetToRandomHPD");
    parent.AssertNotLockedView();
    if( parent.Height() != parent.Width() )
        throw std::logic_error("Positive-definite matrices must be square");
#endif
    const Int width = parent.Width();
    const Int localHeight = parent.LocalHeight();
    const Int r = parent.Grid().Height();
    const Int colShift = parent.ColShift();

    parent.SetToRandom();

    std::complex<Z>* thisLocalBuffer = parent.LocalBuffer();
    const Int thisLDim = parent.LocalLDim();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( Int iLocal=0; iLocal<localHeight; ++iLocal )
    {
        const Int i = colShift + iLocal*r;
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

template<typename T,typename Int>
inline typename RealBase<T>::type
DistMatrix<T,MC,STAR,Int>::GetReal( Int i, Int j ) const
{ return GetRealHelper<T>::Func( *this, i, j ); }

template<typename T,typename Int>
template<typename Z>
inline Z
DistMatrix<T,MC,STAR,Int>::GetRealHelper<Z>::Func
( const DistMatrix<Z,MC,STAR,Int>& parent, Int i, Int j )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::GetRealHelper");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T,typename Int>
template<typename Z>
inline Z
DistMatrix<T,MC,STAR,Int>::GetRealHelper<std::complex<Z> >::Func
( const DistMatrix<std::complex<Z>,MC,STAR,Int>& parent, Int i, Int j )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::GetReal");
    parent.AssertValidEntry( i, j );
#endif
    // We will determine the owner row of entry (i,j) and broadcast from that 
    // row within each process column
    const elemental::Grid& g = parent.Grid();
    const Int ownerRow = (i + parent.ColAlignment()) % g.Height();

    Z u;
    if( g.MCRank() == ownerRow )
    {
        const Int iLoc = (i-parent.ColShift()) / g.Height();
        u = parent.GetRealLocalEntry(iLoc,j);
    }
    mpi::Broadcast( &u, 1, ownerRow, g.MCComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T,typename Int>
inline typename RealBase<T>::type
DistMatrix<T,MC,STAR,Int>::GetImag( Int i, Int j ) const
{ return GetImagHelper<T>::Func( *this, i, j ); }

template<typename T,typename Int>
template<typename Z>
inline Z
DistMatrix<T,MC,STAR,Int>::GetImagHelper<Z>::Func
( const DistMatrix<Z,MC,STAR,Int>& parent, Int i, Int j )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::GetImag");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T,typename Int>
template<typename Z>
inline Z
DistMatrix<T,MC,STAR,Int>::GetImagHelper<std::complex<Z> >::Func
( const DistMatrix<std::complex<Z>,MC,STAR,Int>& parent, Int i, Int j ) 
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::GetImag");
    parent.AssertValidEntry( i, j );
#endif
    // We will determine the owner row of entry (i,j) and broadcast from that 
    // row within each process column
    const elemental::Grid& g = parent.Grid();
    const Int ownerRow = (i + parent.ColAlignment()) % g.Height();

    Z u;
    if( g.MCRank() == ownerRow )
    {
        const Int iLoc = (i-parent.ColShift()) / g.Height();
        u = parent.GetImagLocalEntry(iLoc,j);
    }
    mpi::Broadcast( &u, 1, ownerRow, g.MCComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::SetReal
( Int i, Int j, typename RealBase<T>::type alpha )
{ SetRealHelper<T>::Func( *this, i, j, alpha ); }

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,MC,STAR,Int>::SetRealHelper<Z>::Func
( DistMatrix<Z,MC,STAR,Int>& parent, Int i, Int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::SetReal");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,MC,STAR,Int>::SetRealHelper<std::complex<Z> >::Func
( DistMatrix<std::complex<Z>,MC,STAR,Int>& parent, Int i, Int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::SetReal");
    parent.AssertValidEntry( i, j );
#endif
    const elemental::Grid& g = parent.Grid();
    const Int ownerRow = (i + parent.ColAlignment()) % g.Height();

    if( g.MCRank() == ownerRow )
    {
        const Int iLoc = (i-parent.ColShift()) / g.Height();
        parent.SetRealLocalEntry(iLoc,j,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::SetImag
( Int i, Int j, typename RealBase<T>::type alpha )
{ SetImagHelper<T>::Func( *this, i, j, alpha ); }

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,MC,STAR,Int>::SetImagHelper<Z>::Func
( DistMatrix<Z,MC,STAR,Int>& parent, Int i, Int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::SetImag");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,MC,STAR,Int>::SetImagHelper<std::complex<Z> >::Func
( DistMatrix<std::complex<Z>,MC,STAR,Int>& parent, Int i, Int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::SetImag");
    parent.AssertValidEntry( i, j );
#endif
    const elemental::Grid& g = parent.Grid();
    const Int ownerRow = (i + parent.ColAlignment()) % g.Height();

    if( g.MCRank() == ownerRow )
    {
        const Int iLoc = (i-parent.ColShift()) / g.Height();
        parent.SetImagLocalEntry(iLoc,j,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::UpdateReal
( Int i, Int j, typename RealBase<T>::type alpha )
{ UpdateRealHelper<T>::Func( *this, i, j, alpha ); }

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,MC,STAR,Int>::UpdateRealHelper<Z>::Func
( DistMatrix<Z,MC,STAR,Int>& parent, Int i, Int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::UpdateReal");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,MC,STAR,Int>::UpdateRealHelper<std::complex<Z> >::Func
( DistMatrix<std::complex<Z>,MC,STAR,Int>& parent, Int i, Int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::UpdateReal");
    parent.AssertValidEntry( i, j );
#endif
    const elemental::Grid& g = parent.Grid();
    const Int ownerRow = (i + parent.ColAlignment()) % g.Height();

    if( g.MCRank() == ownerRow )
    {
        const Int iLoc = (i-parent.ColShift()) / g.Height();
        parent.UpdateRealLocalEntry(iLoc,j,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::UpdateImag
( Int i, Int j, typename RealBase<T>::type alpha )
{ UpdateImagHelper<T>::Func( *this, i, j, alpha ); }

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,MC,STAR,Int>::UpdateImagHelper<Z>::Func
( DistMatrix<Z,MC,STAR,Int>& parent, Int i, Int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::UpdateImag");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,MC,STAR,Int>::UpdateImagHelper<std::complex<Z> >::Func
( DistMatrix<std::complex<Z>,MC,STAR,Int>& parent, Int i, Int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::UpdateImag");
    parent.AssertValidEntry( i, j );
#endif
    const elemental::Grid& g = parent.Grid();
    const Int ownerRow = (i + parent.ColAlignment()) % g.Height();

    if( g.MCRank() == ownerRow )
    {
        const Int iLoc = (i-parent.ColShift()) / g.Height();
        parent.UpdateImagLocalEntry(iLoc,j,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::GetRealDiagonal
( DistMatrix<typename RealBase<T>::type,MC,STAR,Int>& d, Int offset ) const
{ GetRealDiagonalHelper<T>::Func( *this, d, offset ); }

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,MC,STAR,Int>::GetRealDiagonalHelper<Z>::Func
( const DistMatrix<Z,MC,STAR,Int>& parent, 
        DistMatrix<Z,MC,STAR,Int>& d, 
  Int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::GetRealDiagonal");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,MC,STAR,Int>::GetRealDiagonalHelper<std::complex<Z> >::Func
( const DistMatrix<std::complex<Z>,MC,STAR,Int>& parent,
        DistMatrix<Z,MC,STAR,Int>& d, Int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::GetRealDiagonal");
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
    const elemental::Grid& g = parent.Grid();
    if( !d.Viewing() )
    {
        d.SetGrid( g );
        if( !d.ConstrainedColAlignment() )
            d.AlignWithDiagonal( parent, offset );
        d.ResizeTo( length, 1 );
    }

    if( g.InGrid() )
    {
        const Int r = g.Height();
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

        const Int iLocalStart = (iStart-colShift) / r;
        const Int localDiagLength = d.LocalHeight();
        const std::complex<Z>* thisLocalBuffer = parent.LockedLocalBuffer();
        const Int thisLDim = parent.LocalLDim();
        Z* dLocalBuffer = d.LocalBuffer();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( Int k=0; k<localDiagLength; ++k )
        {
            const Int iLocal = iLocalStart+k;
            const Int jLocal = jStart+k*r;
            dLocalBuffer[k] = real(thisLocalBuffer[iLocal+jLocal*thisLDim]);
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::GetRealDiagonal
( DistMatrix<typename RealBase<T>::type,STAR,MC,Int>& d, Int offset ) const
{ GetRealDiagonalHelper<T>::Func( *this, d, offset ); }

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,MC,STAR,Int>::GetRealDiagonalHelper<Z>::Func
( const DistMatrix<Z,MC,STAR,Int>& parent, 
        DistMatrix<Z,STAR,MC,Int>& d, Int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::GetRealDiagonal");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,MC,STAR,Int>::GetRealDiagonalHelper<std::complex<Z> >::Func
( const DistMatrix<std::complex<Z>,MC,STAR,Int>& parent,
        DistMatrix<Z,STAR,MC,Int>& d, Int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::GetRealDiagonal");
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
    const elemental::Grid& g = parent.Grid();
    if( !d.Viewing() )
    {
        d.SetGrid( g );
        if( !d.ConstrainedRowAlignment() )
            d.AlignWithDiagonal( parent, offset );
        d.ResizeTo( 1, length );
    }

    if( g.InGrid() )
    {
        const Int r = g.Height();
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

        const Int iLocalStart = (iStart-colShift) / r;
        const Int localDiagLength = d.LocalWidth();
        const std::complex<Z>* thisLocalBuffer = parent.LockedLocalBuffer();
        const Int thisLDim = parent.LocalLDim();
        Z* dLocalBuffer = d.LocalBuffer();
        const Int dLDim = d.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( Int k=0; k<localDiagLength; ++k )
        {
            const Int iLocal = iLocalStart+k;
            const Int jLocal = jStart+k*r;
            dLocalBuffer[k*dLDim] = 
                real(thisLocalBuffer[iLocal+jLocal*thisLDim]);
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::GetImagDiagonal
( DistMatrix<typename RealBase<T>::type,MC,STAR,Int>& d, Int offset ) const
{ GetImagDiagonalHelper<T>::Func( *this, d, offset ); }

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,MC,STAR,Int>::GetImagDiagonalHelper<Z>::Func
( const DistMatrix<Z,MC,STAR,Int>& parent, 
        DistMatrix<Z,MC,STAR,Int>& d, Int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::GetImagDiagonal");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,MC,STAR,Int>::GetImagDiagonalHelper<std::complex<Z> >::Func
( const DistMatrix<std::complex<Z>,MC,STAR,Int>& parent,
        DistMatrix<Z,MC,STAR,Int>& d, Int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::GetImagDiagonal");
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
    const elemental::Grid& g = parent.Grid();
    if( !d.Viewing() )
    {
        d.SetGrid( g );
        if( !d.ConstrainedColAlignment() )
            d.AlignWithDiagonal( parent, offset );
        d.ResizeTo( length, 1 );
    }

    if( g.InGrid() )
    {
        const Int r = g.Height();
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

        const Int iLocalStart = (iStart-colShift) / r;
        const Int localDiagLength = d.LocalHeight();
        const std::complex<Z>* thisLocalBuffer = parent.LockedLocalBuffer();
        const Int thisLDim = parent.LocalLDim();
        Z* dLocalBuffer = d.LocalBuffer();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( Int k=0; k<localDiagLength; ++k )
        {
            const Int iLocal = iLocalStart+k;
            const Int jLocal = jStart+k*r;
            dLocalBuffer[k] = imag(thisLocalBuffer[iLocal+jLocal*thisLDim]);
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::GetImagDiagonal
( DistMatrix<typename RealBase<T>::type,STAR,MC,Int>& d, Int offset ) const
{ GetImagDiagonalHelper<T>::Func( *this, d, offset ); }

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,MC,STAR,Int>::GetImagDiagonalHelper<Z>::Func
( const DistMatrix<Z,MC,STAR,Int>& parent, 
        DistMatrix<Z,STAR,MC,Int>& d, Int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::GetImagDiagonal");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,MC,STAR,Int>::GetImagDiagonalHelper<std::complex<Z> >::Func
( const DistMatrix<std::complex<Z>,MC,STAR,Int>& parent,
        DistMatrix<Z,STAR,MC,Int>& d, Int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::GetImagDiagonal");
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
    const elemental::Grid& g = parent.Grid();
    if( !d.Viewing() )
    {
        d.SetGrid( g );
        if( !d.ConstrainedRowAlignment() )
            d.AlignWithDiagonal( parent, offset );
        d.ResizeTo( 1, length );
    }

    if( g.InGrid() )
    {
        const Int r = g.Height();
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

        const Int iLocalStart = (iStart-colShift) / r;
        const Int localDiagLength = d.LocalWidth();
        const std::complex<Z>* thisLocalBuffer = parent.LockedLocalBuffer();
        const Int thisLDim = parent.LocalLDim();
        Z* dLocalBuffer = d.LocalBuffer();
        const Int dLDim = d.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( Int k=0; k<localDiagLength; ++k )
        {
            const Int iLocal = iLocalStart+k;
            const Int jLocal = jStart+k*r;
            dLocalBuffer[k*dLDim] = 
                imag(thisLocalBuffer[iLocal+jLocal*thisLDim]);
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::SetRealDiagonal
( const DistMatrix<typename RealBase<T>::type,MC,STAR,Int>& d, Int offset )
{ SetRealDiagonalHelper<T>::Func( *this, d, offset ); }

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,MC,STAR,Int>::SetRealDiagonalHelper<Z>::Func
(       DistMatrix<Z,MC,STAR,Int>& parent, 
  const DistMatrix<Z,MC,STAR,Int>& d, Int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::SetRealDiagonal");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,MC,STAR,Int>::SetRealDiagonalHelper<std::complex<Z> >::Func
(       DistMatrix<std::complex<Z>,MC,STAR,Int>& parent,
  const DistMatrix<Z,MC,STAR,Int>& d, Int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::SetRealDiagonal");
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
    const elemental::Grid& g = parent.Grid();
    if( g.InGrid() )
    {
        const Int r = g.Height();
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

        const Int iLocalStart = (iStart-colShift)/r;
        const Int localDiagLength = d.LocalHeight();

        const Z* dLocalBuffer = d.LockedLocalBuffer();
        std::complex<Z>* thisLocalBuffer = parent.LocalBuffer();
        const Int thisLDim = parent.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( Int k=0; k<localDiagLength; ++k )
        {
            const Int iLocal = iLocalStart+k;
            const Int jLocal = jStart+k*r;
            const Z u = dLocalBuffer[k];
            const Z v = imag(thisLocalBuffer[iLocal+jLocal*thisLDim]);
            thisLocalBuffer[iLocal+jLocal*thisLDim] = std::complex<Z>(u,v);
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::SetRealDiagonal
( const DistMatrix<typename RealBase<T>::type,STAR,MC,Int>& d, Int offset )
{ SetRealDiagonalHelper<T>::Func( *this, d, offset ); }

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,MC,STAR,Int>::SetRealDiagonalHelper<Z>::Func
(       DistMatrix<Z,MC,STAR,Int>& parent, 
  const DistMatrix<Z,STAR,MC,Int>& d, Int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::SetRealDiagonal");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,MC,STAR,Int>::SetRealDiagonalHelper<std::complex<Z> >::Func
(       DistMatrix<std::complex<Z>,MC,STAR,Int>& parent,
  const DistMatrix<Z,STAR,MC,Int>& d, Int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::SetRealDiagonal");
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
    const elemental::Grid& g = parent.Grid();
    if( g.InGrid() )
    {
        const Int r = g.Height();
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

        const Int iLocalStart = (iStart-colShift)/r;
        const Int localDiagLength = d.LocalWidth();

        const Z* dLocalBuffer = d.LockedLocalBuffer();
        std::complex<Z>* thisLocalBuffer = parent.LocalBuffer();
        const Int dLDim = d.LocalLDim();
        const Int thisLDim = parent.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( Int k=0; k<localDiagLength; ++k )
        {
            const Int iLocal = iLocalStart+k;
            const Int jLocal = jStart+k*r;
            const Z u = dLocalBuffer[k*dLDim];
            const Z v = imag(thisLocalBuffer[iLocal+jLocal*thisLDim]);
            thisLocalBuffer[iLocal+jLocal*thisLDim] = std::complex<Z>(u,v);
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::SetImagDiagonal
( const DistMatrix<typename RealBase<T>::type,MC,STAR,Int>& d, Int offset )
{ SetImagDiagonalHelper<T>::Func( *this, d, offset ); }

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,MC,STAR,Int>::SetImagDiagonalHelper<Z>::Func
(       DistMatrix<Z,MC,STAR,Int>& parent, 
  const DistMatrix<Z,MC,STAR,Int>& d, Int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::SetImagDiagonal");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,MC,STAR,Int>::SetImagDiagonalHelper<std::complex<Z> >::Func
(       DistMatrix<std::complex<Z>,MC,STAR,Int>& parent,
  const DistMatrix<Z,MC,STAR,Int>& d, Int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::SetImagDiagonal");
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
    const elemental::Grid& g = parent.Grid();
    if( g.InGrid() )
    {
        const Int r = g.Height();
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

        const Int iLocalStart = (iStart-colShift)/r;
        const Int localDiagLength = d.LocalHeight();

        const Z* dLocalBuffer = d.LockedLocalBuffer();
        std::complex<Z>* thisLocalBuffer = parent.LocalBuffer();
        const Int thisLDim = parent.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( Int k=0; k<localDiagLength; ++k )
        {
            const Int iLocal = iLocalStart+k;
            const Int jLocal = jStart+k*r;
            const Z u = real(thisLocalBuffer[iLocal+jLocal*thisLDim]);
            const Z v = dLocalBuffer[k];
            thisLocalBuffer[iLocal+jLocal*thisLDim] = std::complex<Z>(u,v);
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::SetImagDiagonal
( const DistMatrix<typename RealBase<T>::type,STAR,MC,Int>& d, Int offset )
{ SetImagDiagonalHelper<T>::Func( *this, d, offset ); }

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,MC,STAR,Int>::SetImagDiagonalHelper<Z>::Func
(       DistMatrix<Z,MC,STAR,Int>& parent, 
  const DistMatrix<Z,STAR,MC,Int>& d, Int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::SetImagDiagonal");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,MC,STAR,Int>::SetImagDiagonalHelper<std::complex<Z> >::Func
(       DistMatrix<std::complex<Z>,MC,STAR,Int>& parent,
  const DistMatrix<Z,STAR,MC,Int>& d, Int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::SetImagDiagonal");
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
    const elemental::Grid& g = parent.Grid();
    if( g.InGrid() )
    {
        const Int r = g.Height();
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

        const Int iLocalStart = (iStart-colShift)/r;
        const Int localDiagLength = d.LocalWidth();

        const Z* dLocalBuffer = d.LockedLocalBuffer();
        std::complex<Z>* thisLocalBuffer = parent.LocalBuffer();
        const Int dLDim = d.LocalLDim();
        const Int thisLDim = parent.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( Int k=0; k<localDiagLength; ++k )
        {
            const Int iLocal = iLocalStart+k;
            const Int jLocal = jStart+k*r;
            const Z u = real(thisLocalBuffer[iLocal+jLocal*thisLDim]);
            const Z v = dLocalBuffer[k*dLDim];
            thisLocalBuffer[iLocal+jLocal*thisLDim] = std::complex<Z>(u,v);
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elemental
