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
DistMatrix<T,MC,STAR>::SetToRandomHermitianHelper<Z>::Func
( DistMatrix<Z,MC,STAR>& parent )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::SetToRandomHermitian");
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
DistMatrix<T,MC,STAR>::SetToRandomHermitianHelper<complex<Z> >::Func
( DistMatrix<complex<Z>,MC,STAR>& parent )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::SetToRandomHermitian");
    parent.AssertNotLockedView();
    if( parent.Height() != parent.Width() )
        throw logic_error( "Hermitian matrices must be square." );
#endif
    const int width = parent.Width();
    const int localHeight = parent.LocalHeight();
    const int r = parent.Grid().Height();
    const int colShift = parent.ColShift();

    parent.SetToRandom();

    complex<Z>* thisLocalBuffer = parent.LocalBuffer();
    const int thisLDim = parent.LocalLDim();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( int iLocal=0; iLocal<localHeight; ++iLocal )
    {
        const int i = colShift + iLocal*r;
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
DistMatrix<T,MC,STAR>::SetToRandomHPDHelper<Z>::Func
( DistMatrix<Z,MC,STAR>& parent )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::SetToRandomHPD");
    parent.AssertNotLockedView();
    if( parent.Height() != parent.Width() )
        throw logic_error( "Positive-definite matrices must be square." );
#endif
    const int width = parent.Width();
    const int localHeight = parent.LocalHeight();
    const int r = parent.Grid().Height();
    const int colShift = parent.ColShift();

    parent.SetToRandom();

    Z* thisLocalBuffer = parent.LocalBuffer();
    const int thisLDim = parent.LocalLDim();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( int iLocal=0; iLocal<localHeight; ++iLocal )
    {
        const int i = colShift + iLocal*r;
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
DistMatrix<T,MC,STAR>::SetToRandomHPDHelper<complex<Z> >::Func
( DistMatrix<complex<Z>,MC,STAR>& parent )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::SetToRandomHPD");
    parent.AssertNotLockedView();
    if( parent.Height() != parent.Width() )
        throw logic_error( "Positive-definite matrices must be square." );
#endif
    const int width = parent.Width();
    const int localHeight = parent.LocalHeight();
    const int r = parent.Grid().Height();
    const int colShift = parent.ColShift();

    parent.SetToRandom();

    complex<Z>* thisLocalBuffer = parent.LocalBuffer();
    const int thisLDim = parent.LocalLDim();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( int iLocal=0; iLocal<localHeight; ++iLocal )
    {
        const int i = colShift + iLocal*r;
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
DistMatrix<T,MC,STAR>::GetRealHelper<complex<Z> >::Func
( const DistMatrix<complex<Z>,MC,STAR>& parent, int i, int j )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::GetReal");
    parent.AssertValidEntry( i, j );
#endif
    // We will determine the owner row of entry (i,j) and broadcast from that 
    // row within each process column
    const elemental::Grid& g = parent.Grid();
    const int ownerRow = (i + parent.ColAlignment()) % g.Height();

    Z u;
    if( g.MCRank() == ownerRow )
    {
        const int iLoc = (i-parent.ColShift()) / g.Height();
        u = parent.GetRealLocalEntry(iLoc,j);
    }
    mpi::Broadcast( &u, 1, ownerRow, g.MCComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T>
template<typename Z>
inline Z
DistMatrix<T,MC,STAR>::GetImagHelper<complex<Z> >::Func
( const DistMatrix<complex<Z>,MC,STAR>& parent, int i, int j ) 
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::GetImag");
    parent.AssertValidEntry( i, j );
#endif
    // We will determine the owner row of entry (i,j) and broadcast from that 
    // row within each process column
    const elemental::Grid& g = parent.Grid();
    const int ownerRow = (i + parent.ColAlignment()) % g.Height();

    Z u;
    if( g.MCRank() == ownerRow )
    {
        const int iLoc = (i-parent.ColShift()) / g.Height();
        u = parent.GetImagLocalEntry(iLoc,j);
    }
    mpi::Broadcast( &u, 1, ownerRow, g.MCComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T>
template<typename Z>
inline void
DistMatrix<T,MC,STAR>::SetRealHelper<complex<Z> >::Func
( DistMatrix<complex<Z>,MC,STAR>& parent, int i, int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::SetReal");
    parent.AssertValidEntry( i, j );
#endif
    const elemental::Grid& g = parent.Grid();
    const int ownerRow = (i + parent.ColAlignment()) % g.Height();

    if( g.MCRank() == ownerRow )
    {
        const int iLoc = (i-parent.ColShift()) / g.Height();
        parent.SetRealLocalEntry(iLoc,j,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
template<typename Z>
inline void
DistMatrix<T,MC,STAR>::SetImagHelper<complex<Z> >::Func
( DistMatrix<complex<Z>,MC,STAR>& parent, int i, int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::SetImag");
    parent.AssertValidEntry( i, j );
#endif
    const elemental::Grid& g = parent.Grid();
    const int ownerRow = (i + parent.ColAlignment()) % g.Height();

    if( g.MCRank() == ownerRow )
    {
        const int iLoc = (i-parent.ColShift()) / g.Height();
        parent.SetImagLocalEntry(iLoc,j,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
template<typename Z>
inline void
DistMatrix<T,MC,STAR>::UpdateRealHelper<complex<Z> >::Func
( DistMatrix<complex<Z>,MC,STAR>& parent, int i, int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::UpdateReal");
    parent.AssertValidEntry( i, j );
#endif
    const elemental::Grid& g = parent.Grid();
    const int ownerRow = (i + parent.ColAlignment()) % g.Height();

    if( g.MCRank() == ownerRow )
    {
        const int iLoc = (i-parent.ColShift()) / g.Height();
        parent.UpdateRealLocalEntry(iLoc,j,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
template<typename Z>
inline void
DistMatrix<T,MC,STAR>::UpdateImagHelper<complex<Z> >::Func
( DistMatrix<complex<Z>,MC,STAR>& parent, int i, int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::UpdateImag");
    parent.AssertValidEntry( i, j );
#endif
    const elemental::Grid& g = parent.Grid();
    const int ownerRow = (i + parent.ColAlignment()) % g.Height();

    if( g.MCRank() == ownerRow )
    {
        const int iLoc = (i-parent.ColShift()) / g.Height();
        parent.UpdateImagLocalEntry(iLoc,j,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
template<typename Z>
inline void
DistMatrix<T,MC,STAR>::GetRealDiagonalHelper<complex<Z> >::Func
( const DistMatrix<complex<Z>,MC,STAR>& parent,
        DistMatrix<Z,MC,STAR>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::GetRealDiagonal");
    if( d.Viewing() )
        parent.AssertSameGrid( d );
#endif
    const int length = parent.DiagonalLength(offset);
#ifndef RELEASE
    if( d.Viewing() && (length != d.Height() || d.Width() != 1) )
    {
        ostringstream msg; 
        msg << "d is not a column vec of the same length as the diagonal:\n"
            << "  A ~ " << parent.Height() << " x " << parent.Width() << "\n"
            << "  d ~ " << d.Height() << " x " << d.Width() << "\n"
            << "  A diag length: " << length << "\n";
        throw logic_error( msg.str().c_str() );
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
        const int r = g.Height();
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

        const int iLocalStart = (iStart-colShift) / r;
        const int localDiagLength = d.LocalHeight();
        const complex<Z>* thisLocalBuffer = parent.LockedLocalBuffer();
        const int thisLDim = parent.LocalLDim();
        Z* dLocalBuffer = d.LocalBuffer();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<localDiagLength; ++k )
        {
            const int iLocal = iLocalStart+k;
            const int jLocal = jStart+k*r;
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
DistMatrix<T,MC,STAR>::GetImagDiagonalHelper<complex<Z> >::Func
( const DistMatrix<complex<Z>,MC,STAR>& parent,
        DistMatrix<Z,MC,STAR>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::GetImagDiagonal");
    if( d.Viewing() )
        parent.AssertSameGrid( d );
#endif
    const int length = parent.DiagonalLength(offset);
#ifndef RELEASE
    if( d.Viewing() && length != d.Height() )
    {
        ostringstream msg; 
        msg << "d is not of the same length as the diagonal:\n"
            << "  A ~ " << parent.Height() << " x " << parent.Width() << "\n"
            << "  d ~ " << d.Height() << " x " << d.Width() << "\n"
            << "  A diag length: " << length << "\n";
        throw logic_error( msg.str().c_str() );
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
        const int r = g.Height();
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

        const int iLocalStart = (iStart-colShift) / r;
        const int localDiagLength = d.LocalHeight();
        const complex<Z>* thisLocalBuffer = parent.LockedLocalBuffer();
        const int thisLDim = parent.LocalLDim();
        Z* dLocalBuffer = d.LocalBuffer();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<localDiagLength; ++k )
        {
            const int iLocal = iLocalStart+k;
            const int jLocal = jStart+k*r;
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
DistMatrix<T,MC,STAR>::GetRealDiagonalHelper<complex<Z> >::Func
( const DistMatrix<complex<Z>,MC,STAR>& parent,
        DistMatrix<Z,STAR,MC>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::GetRealDiagonal");
    if( d.Viewing() )
        parent.AssertSameGrid( d );
#endif
    const int length = parent.DiagonalLength(offset);
#ifndef RELEASE
    if( d.Viewing() && length != d.Width() )
    {
        ostringstream msg; 
        msg << "d is not of the same length as the diagonal:\n"
            << "  A ~ " << parent.Height() << " x " << parent.Width() << "\n"
            << "  d ~ " << d.Height() << " x " << d.Width() << "\n"
            << "  A diag length: " << length << "\n";
        throw logic_error( msg.str().c_str() );
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
        const int r = g.Height();
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

        const int iLocalStart = (iStart-colShift) / r;
        const int localDiagLength = d.LocalWidth();
        const complex<Z>* thisLocalBuffer = parent.LockedLocalBuffer();
        const int thisLDim = parent.LocalLDim();
        Z* dLocalBuffer = d.LocalBuffer();
        const int dLDim = d.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<localDiagLength; ++k )
        {
            const int iLocal = iLocalStart+k;
            const int jLocal = jStart+k*r;
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
DistMatrix<T,MC,STAR>::GetImagDiagonalHelper<complex<Z> >::Func
( const DistMatrix<complex<Z>,MC,STAR>& parent,
        DistMatrix<Z,STAR,MC>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::GetImagDiagonal");
    if( d.Viewing() )
        parent.AssertSameGrid( d );
#endif
    const int length = parent.DiagonalLength(offset);
#ifndef RELEASE
    if( d.Viewing() && length != d.Width() )
    {
        ostringstream msg; 
        msg << "d is not of the same length as the diagonal:\n"
            << "  A ~ " << parent.Height() << " x " << parent.Width() << "\n"
            << "  d ~ " << d.Height() << " x " << d.Width() << "\n"
            << "  A diag length: " << length << "\n";
        throw logic_error( msg.str().c_str() );
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
        const int r = g.Height();
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

        const int iLocalStart = (iStart-colShift) / r;
        const int localDiagLength = d.LocalWidth();
        const complex<Z>* thisLocalBuffer = parent.LockedLocalBuffer();
        const int thisLDim = parent.LocalLDim();
        Z* dLocalBuffer = d.LocalBuffer();
        const int dLDim = d.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<localDiagLength; ++k )
        {
            const int iLocal = iLocalStart+k;
            const int jLocal = jStart+k*r;
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
DistMatrix<T,MC,STAR>::SetRealDiagonalHelper<complex<Z> >::Func
(      DistMatrix<complex<Z>,MC,STAR>& parent,
 const DistMatrix<Z,MC,STAR>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::SetRealDiagonal");
    parent.AssertSameGrid( d );
    if( d.Width() != 1 )
        throw std::logic_error("d must be a column vector");
    const int length = parent.DiagonalLength(offset);
    if( length != d.Height() )
    {
        ostringstream msg;
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
        const int r = g.Height();
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

        const int iLocalStart = (iStart-colShift)/r;
        const int localDiagLength = d.LocalHeight();

        const Z* dLocalBuffer = d.LockedLocalBuffer();
        complex<Z>* thisLocalBuffer = parent.LocalBuffer();
        const int thisLDim = parent.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<localDiagLength; ++k )
        {
            const int iLocal = iLocalStart+k;
            const int jLocal = jStart+k*r;
            const Z u = dLocalBuffer[k];
            const Z v = imag(thisLocalBuffer[iLocal+jLocal*thisLDim]);
            thisLocalBuffer[iLocal+jLocal*thisLDim] = complex<Z>(u,v);
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
template<typename Z>
inline void
DistMatrix<T,MC,STAR>::SetImagDiagonalHelper<complex<Z> >::Func
(      DistMatrix<complex<Z>,MC,STAR>& parent,
 const DistMatrix<Z,MC,STAR>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::SetImagDiagonal");
    parent.AssertSameGrid( d );
    if( d.Width() != 1 )
        throw std::logic_error("d must be a column vector");
    const int length = parent.DiagonalLength(offset);
    if( length != d.Height() )
    {
        ostringstream msg;
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
        const int r = g.Height();
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

        const int iLocalStart = (iStart-colShift)/r;
        const int localDiagLength = d.LocalHeight();

        const Z* dLocalBuffer = d.LockedLocalBuffer();
        complex<Z>* thisLocalBuffer = parent.LocalBuffer();
        const int thisLDim = parent.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<localDiagLength; ++k )
        {
            const int iLocal = iLocalStart+k;
            const int jLocal = jStart+k*r;
            const Z u = real(thisLocalBuffer[iLocal+jLocal*thisLDim]);
            const Z v = dLocalBuffer[k];
            thisLocalBuffer[iLocal+jLocal*thisLDim] = complex<Z>(u,v);
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
template<typename Z>
inline void
DistMatrix<T,MC,STAR>::SetRealDiagonalHelper<complex<Z> >::Func
(      DistMatrix<complex<Z>,MC,STAR>& parent,
 const DistMatrix<Z,STAR,MC>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::SetRealDiagonal");
    parent.AssertSameGrid( d );
    if( d.Height() != 1 )
        throw std::logic_error("d must be a row vector");
    const int length = parent.DiagonalLength(offset);
    if( length != d.Width() )
    {
        ostringstream msg;
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
        const int r = g.Height();
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

        const int iLocalStart = (iStart-colShift)/r;
        const int localDiagLength = d.LocalWidth();

        const Z* dLocalBuffer = d.LockedLocalBuffer();
        complex<Z>* thisLocalBuffer = parent.LocalBuffer();
        const int dLDim = d.LocalLDim();
        const int thisLDim = parent.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<localDiagLength; ++k )
        {
            const int iLocal = iLocalStart+k;
            const int jLocal = jStart+k*r;
            const Z u = dLocalBuffer[k*dLDim];
            const Z v = imag(thisLocalBuffer[iLocal+jLocal*thisLDim]);
            thisLocalBuffer[iLocal+jLocal*thisLDim] = complex<Z>(u,v);
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
template<typename Z>
inline void
DistMatrix<T,MC,STAR>::SetImagDiagonalHelper<complex<Z> >::Func
(      DistMatrix<complex<Z>,MC,STAR>& parent,
 const DistMatrix<Z,STAR,MC>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::SetImagDiagonal");
    parent.AssertSameGrid( d );
    if( d.Height() != 1 )
        throw std::logic_error("d must be a row vector");
    const int length = parent.DiagonalLength(offset);
    if( length != d.Width() )
    {
        ostringstream msg;
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
        const int r = g.Height();
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

        const int iLocalStart = (iStart-colShift)/r;
        const int localDiagLength = d.LocalWidth();

        const Z* dLocalBuffer = d.LockedLocalBuffer();
        complex<Z>* thisLocalBuffer = parent.LocalBuffer();
        const int dLDim = d.LocalLDim();
        const int thisLDim = parent.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<localDiagLength; ++k )
        {
            const int iLocal = iLocalStart+k;
            const int jLocal = jStart+k*r;
            const Z u = real(thisLocalBuffer[iLocal+jLocal*thisLDim]);
            const Z v = dLocalBuffer[k*dLDim];
            thisLocalBuffer[iLocal+jLocal*thisLDim] = complex<Z>(u,v);
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX

} // namespace elemental
