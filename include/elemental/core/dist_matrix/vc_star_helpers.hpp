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
DistMatrix<T,VC,STAR>::SetToRandomHermitianHelper<Z>::Func
( DistMatrix<Z,VC,STAR>& parent )
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

#ifndef WITHOUT_COMPLEX
template<typename T>
template<typename Z>
inline void
DistMatrix<T,VC,STAR>::SetToRandomHermitianHelper<std::complex<Z> >::Func
( DistMatrix<std::complex<Z>,VC,STAR>& parent )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::SetToRandomHermitian");
    parent.AssertNotLockedView();
    if( parent.Height() != parent.Width() )
        throw std::logic_error("Hermitian matrices must be square");
#endif
    const int width       = parent.Width();
    const int localHeight = parent.LocalHeight();
    const int p           = parent.Grid().Size();
    const int colShift    = parent.ColShift();

    parent.SetToRandom();

    std::complex<Z>* thisLocalBuffer = parent.LocalBuffer();
    const int thisLDim = parent.LocalLDim();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( int iLocal=0; iLocal<localHeight; ++iLocal )
    {
        const int i = colShift + iLocal*p;
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
DistMatrix<T,VC,STAR>::SetToRandomHPDHelper<Z>::Func
( DistMatrix<Z,VC,STAR>& parent )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::SetToRandomHPD");
    parent.AssertNotLockedView();
    if( parent.Height() != parent.Width() )
        throw std::logic_error("Positive-definite matrices must be square");
#endif
    const int width       = parent.Width();
    const int localHeight = parent.LocalHeight();
    const int p           = parent.Grid().Size();
    const int colShift    = parent.ColShift();

    parent.SetToRandom();

    Z* thisLocalBuffer = parent.LocalBuffer();
    const int thisLDim = parent.LocalLDim();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( int iLocal=0; iLocal<localHeight; ++iLocal )
    {
        const int i = colShift + iLocal*p;
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
DistMatrix<T,VC,STAR>::SetToRandomHPDHelper<std::complex<Z> >::Func
( DistMatrix<std::complex<Z>,VC,STAR>& parent )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::SetToRandomHPD");
    parent.AssertNotLockedView();
    if( parent.Height() != parent.Width() )
        throw std::logic_error("Positive-definite matrices must be square");
#endif
    const int width       = parent.Width();
    const int localHeight = parent.LocalHeight();
    const int p           = parent.Grid().Size();
    const int colShift    = parent.ColShift();

    parent.SetToRandom();

    std::complex<Z>* thisLocalBuffer = parent.LocalBuffer();
    const int thisLDim = parent.LocalLDim();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( int iLocal=0; iLocal<localHeight; ++iLocal )
    {
        const int i = colShift + iLocal*p;
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
DistMatrix<T,VC,STAR>::GetRealHelper<std::complex<Z> >::Func
( const DistMatrix<std::complex<Z>,VC,STAR>& parent, int i, int j ) 
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::GetReal");
    parent.AssertValidEntry( i, j );
#endif
    // We will determine the owner rank of entry (i,j) and broadcast from that
    // process over the entire g
    const elemental::Grid& g = parent.Grid();
    const int ownerRank = (i + parent.ColAlignment()) % g.Size();

    Z u;
    if( g.VCRank() == ownerRank )
    {
        const int iLoc = (i-parent.ColShift()) / g.Size();
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
DistMatrix<T,VC,STAR>::GetImagHelper<std::complex<Z> >::Func
( const DistMatrix<std::complex<Z>,VC,STAR>& parent, int i, int j ) 
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::GetImag");
    parent.AssertValidEntry( i, j );
#endif
    // We will determine the owner rank of entry (i,j) and broadcast from that
    // process over the entire g
    const elemental::Grid& g = parent.Grid();
    const int ownerRank = (i + parent.ColAlignment()) % g.Size();

    Z u;
    if( g.VCRank() == ownerRank )
    {
        const int iLoc = (i-parent.ColShift()) / g.Size();
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
DistMatrix<T,VC,STAR>::SetRealHelper<std::complex<Z> >::Func
( DistMatrix<std::complex<Z>,VC,STAR>& parent, int i, int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::SetReal");
    parent.AssertValidEntry( i, j );
#endif
    const elemental::Grid& g = parent.Grid();
    const int ownerRank = (i + parent.ColAlignment()) % g.Size();

    if( g.VCRank() == ownerRank )
    {
        const int iLoc = (i-parent.ColShift()) / g.Size();
        parent.SetRealLocalEntry(iLoc,j,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
template<typename Z>
inline void
DistMatrix<T,VC,STAR>::SetImagHelper<std::complex<Z> >::Func
( DistMatrix<std::complex<Z>,VC,STAR>& parent, int i, int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::SetImag");
    parent.AssertValidEntry( i, j );
#endif
    const elemental::Grid& g = parent.Grid();
    const int ownerRank = (i + parent.ColAlignment()) % g.Size();

    if( g.VCRank() == ownerRank )
    {
        const int iLoc = (i-parent.ColShift()) / g.Size();
        parent.SetImagLocalEntry(iLoc,j,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
template<typename Z>
inline void
DistMatrix<T,VC,STAR>::UpdateRealHelper<std::complex<Z> >::Func
( DistMatrix<std::complex<Z>,VC,STAR>& parent, int i, int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::UpdateReal");
    parent.AssertValidEntry( i, j );
#endif
    const elemental::Grid& g = parent.Grid();
    const int ownerRank = (i + parent.ColAlignment()) % g.Size();

    if( g.VCRank() == ownerRank )
    {
        const int iLoc = (i-parent.ColShift()) / g.Size();
        parent.UpdateRealLocalEntry(iLoc,j,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
template<typename Z>
void
DistMatrix<T,VC,STAR>::UpdateImagHelper<std::complex<Z> >::Func
( DistMatrix<std::complex<Z>,VC,STAR>& parent, int i, int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::UpdateImag");
    parent.AssertValidEntry( i, j );
#endif
    const elemental::Grid& g = parent.Grid();
    const int ownerRank = (i + parent.ColAlignment()) % g.Size();

    if( g.VCRank() == ownerRank )
    {
        const int iLoc = (i-parent.ColShift()) / g.Size();
        parent.UpdateImagLocalEntry(iLoc,j,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
template<typename Z>
inline void
DistMatrix<T,VC,STAR>::GetRealDiagonalHelper<std::complex<Z> >::Func
( const DistMatrix<std::complex<Z>,VC,STAR>& parent,
        DistMatrix<Z,VC,STAR>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::GetRealDiagonal");
    if( d.Viewing() )
        parent.AssertSameGrid( d );
#endif
    const int length = parent.DiagonalLength(offset);
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
        const int p = g.Size();
        const int colShift = parent.ColShift();
        const int diagShift = d.ColShift();

        int iStart, jStart;
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

        const int iLocalStart = (iStart-colShift) / p;
        const int localDiagLength = d.LocalHeight();
        const std::complex<Z>* thisLocalBuffer = parent.LockedLocalBuffer();
        const int thisLDim = parent.LocalLDim();
        Z* dLocalBuffer = d.LocalBuffer();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<localDiagLength; ++k )
        {
            const int iLocal = iLocalStart+k;
            const int jLocal = jStart+k*p;
            dLocalBuffer[k] = real(thisLocalBuffer[iLocal+jLocal*thisLDim]);
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
template<typename Z>
inline void
DistMatrix<T,VC,STAR>::GetImagDiagonalHelper<std::complex<Z> >::Func
( const DistMatrix<std::complex<Z>,VC,STAR>& parent,
        DistMatrix<Z,VC,STAR>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::GetImagDiagonal");
    if( d.Viewing() )
        parent.AssertSameGrid( d );
#endif
    const int length = parent.DiagonalLength(offset);
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
        const int p = g.Size();
        const int colShift = parent.ColShift();
        const int diagShift = d.ColShift();

        int iStart, jStart;
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

        const int iLocalStart = (iStart-colShift) / p;
        const int localDiagLength = d.LocalHeight();
        const std::complex<Z>* thisLocalBuffer = parent.LockedLocalBuffer();
        const int thisLDim = parent.LocalLDim();
        Z* dLocalBuffer = d.LocalBuffer();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<localDiagLength; ++k )
        {
            const int iLocal = iLocalStart+k;
            const int jLocal = jStart+k*p;
            dLocalBuffer[k] = imag(thisLocalBuffer[iLocal+jLocal*thisLDim]);
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
template<typename Z>
inline void
DistMatrix<T,VC,STAR>::GetRealDiagonalHelper<std::complex<Z> >::Func
( const DistMatrix<std::complex<Z>,VC,STAR>& parent,
        DistMatrix<Z,STAR,VC>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::GetRealDiagonal");
    if( d.Viewing() )
        parent.AssertSameGrid( d );
#endif
    const int length = parent.DiagonalLength(offset);
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
        const int p = g.Size();
        const int colShift = parent.ColShift();
        const int diagShift = d.RowShift();

        int iStart, jStart;
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

        const int iLocalStart = (iStart-colShift) / p;
        const int localDiagLength = d.LocalWidth();
        const std::complex<Z>* thisLocalBuffer = parent.LockedLocalBuffer();
        const int thisLDim = parent.LocalLDim();
        Z* dLocalBuffer = d.LocalBuffer();
        const int dLDim = d.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<localDiagLength; ++k )
        {
            const int iLocal = iLocalStart+k;
            const int jLocal = jStart+k*p;
            dLocalBuffer[k*dLDim] =
                real(thisLocalBuffer[iLocal+jLocal*thisLDim]);
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
template<typename Z>
inline void
DistMatrix<T,VC,STAR>::GetImagDiagonalHelper<std::complex<Z> >::Func
( const DistMatrix<std::complex<Z>,VC,STAR>& parent,
        DistMatrix<Z,STAR,VC>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::GetImagDiagonal");
    if( d.Viewing() )
        parent.AssertSameGrid( d );
#endif
    const int length = parent.DiagonalLength(offset);
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
        const int p = g.Size();
        const int colShift = parent.ColShift();
        const int diagShift = d.RowShift();

        int iStart, jStart;
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

        const int iLocalStart = (iStart-colShift) / p;
        const int localDiagLength = d.LocalWidth();
        const std::complex<Z>* thisLocalBuffer = parent.LockedLocalBuffer();
        const int thisLDim = parent.LocalLDim();
        Z* dLocalBuffer = d.LocalBuffer();
        const int dLDim = d.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<localDiagLength; ++k )
        {
            const int iLocal = iLocalStart+k;
            const int jLocal = jStart+k*p;
            dLocalBuffer[k*dLDim] =
                imag(thisLocalBuffer[iLocal+jLocal*thisLDim]);
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
template<typename Z>
inline void
DistMatrix<T,VC,STAR>::SetRealDiagonalHelper<std::complex<Z> >::Func
(      DistMatrix<std::complex<Z>,VC,STAR>& parent,
 const DistMatrix<Z,VC,STAR>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::SetRealDiagonal");
    parent.AssertSameGrid( d );
    if( d.Width() != 1 )
        throw std::logic_error("d must be a column vector");
    const int length = parent.DiagonalLength(offset);
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
        const int p = g.Size();
        const int colShift = parent.ColShift();
        const int diagShift = d.ColShift();

        int iStart, jStart;
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

        const int iLocalStart = (iStart-colShift)/p;
        const int localDiagLength = d.LocalHeight();

        const Z* dLocalBuffer = d.LockedLocalBuffer();
        std::complex<Z>* thisLocalBuffer = parent.LocalBuffer();
        const int thisLDim = parent.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<localDiagLength; ++k )
        {
            const int iLocal = iLocalStart+k;
            const int jLocal = jStart+k*p;
            const Z u = dLocalBuffer[k];
            const Z v = imag(thisLocalBuffer[iLocal+jLocal*thisLDim]);
            thisLocalBuffer[iLocal+jLocal*thisLDim] = std::complex<Z>(u,v);
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
template<typename Z>
inline void
DistMatrix<T,VC,STAR>::SetImagDiagonalHelper<std::complex<Z> >::Func
(      DistMatrix<std::complex<Z>,VC,STAR>& parent,
 const DistMatrix<Z,VC,STAR>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::SetImagDiagonal");
    parent.AssertSameGrid( d );
    if( d.Width() != 1 )
        throw std::logic_error("d must be a column vector");
    const int length = parent.DiagonalLength(offset);
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
        const int p = g.Size();
        const int colShift = parent.ColShift();
        const int diagShift = d.ColShift();

        int iStart, jStart;
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

        const int iLocalStart = (iStart-colShift)/p;
        const int localDiagLength = d.LocalHeight();

        const Z* dLocalBuffer = d.LockedLocalBuffer();
        std::complex<Z>* thisLocalBuffer = parent.LocalBuffer();
        const int thisLDim = parent.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<localDiagLength; ++k )
        {
            const int iLocal = iLocalStart+k;
            const int jLocal = jStart+k*p;
            const Z u = real(thisLocalBuffer[iLocal+jLocal*thisLDim]);
            const Z v = dLocalBuffer[k];
            thisLocalBuffer[iLocal+jLocal*thisLDim] = std::complex<Z>(u,v);
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
template<typename Z>
inline void
DistMatrix<T,VC,STAR>::SetRealDiagonalHelper<std::complex<Z> >::Func
(      DistMatrix<std::complex<Z>,VC,STAR>& parent,
 const DistMatrix<Z,STAR,VC>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::SetRealDiagonal");
    parent.AssertSameGrid( d );
    if( d.Height() != 1 )
        throw std::logic_error("d must be a row vector");
    const int length = parent.DiagonalLength(offset);
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
        const int p = g.Size();
        const int colShift = parent.ColShift();
        const int diagShift = d.RowShift();

        int iStart, jStart;
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

        const int iLocalStart = (iStart-colShift)/p;
        const int localDiagLength = d.LocalWidth();

        const Z* dLocalBuffer = d.LockedLocalBuffer();
        std::complex<Z>* thisLocalBuffer = parent.LocalBuffer();
        const int dLDim = d.LocalLDim();
        const int thisLDim = parent.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<localDiagLength; ++k )
        {
            const int iLocal = iLocalStart+k;
            const int jLocal = jStart+k*p;
            const Z u = dLocalBuffer[k*dLDim];
            const Z v = imag(thisLocalBuffer[iLocal+jLocal*thisLDim]);
            thisLocalBuffer[iLocal+jLocal*thisLDim] = std::complex<Z>(u,v);
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
template<typename Z>
inline void
DistMatrix<T,VC,STAR>::SetImagDiagonalHelper<std::complex<Z> >::Func
(      DistMatrix<std::complex<Z>,VC,STAR>& parent,
 const DistMatrix<Z,STAR,VC>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::SetImagDiagonal");
    parent.AssertSameGrid( d );
    if( d.Height() != 1 )
        throw std::logic_error("d must be a row vector");
    const int length = parent.DiagonalLength(offset);
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
        const int p = g.Size();
        const int colShift = parent.ColShift();
        const int diagShift = d.RowShift();

        int iStart, jStart;
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

        const int iLocalStart = (iStart-colShift)/p;
        const int localDiagLength = d.LocalWidth();

        const Z* dLocalBuffer = d.LockedLocalBuffer();
        std::complex<Z>* thisLocalBuffer = parent.LocalBuffer();
        const int dLDim = d.LocalLDim();
        const int thisLDim = parent.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<localDiagLength; ++k )
        {
            const int iLocal = iLocalStart+k;
            const int jLocal = jStart+k*p;
            const Z u = real(thisLocalBuffer[iLocal+jLocal*thisLDim]);
            const Z v = dLocalBuffer[k*dLDim];
            thisLocalBuffer[iLocal+jLocal*thisLDim] = std::complex<Z>(u,v);
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX

} // namespace elemental
