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
DistMatrix<T,MC,MR>::SetToRandomHermitianHelper<Z>::Func
( DistMatrix<Z,MC,MR>& parent )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::SetToRandomHermitian");
    parent.AssertNotLockedView();
    if( parent._height != parent._width )
        throw logic_error("Hermitian matrices must be square");
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
DistMatrix<T,MC,MR>::SetToRandomHermitianHelper<complex<Z> >::Func
( DistMatrix<complex<Z>,MC,MR>& parent )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::SetToRandomHermitian");
    parent.AssertNotLockedView();
    if( parent._height != parent._width )
        throw logic_error("Hermitian matrices must be square");
#endif
    const int r = parent.Grid().Height();
    const int c = parent.Grid().Width();

    const int localHeight = parent.LocalHeight();
    const int localWidth = parent.LocalWidth();
    const int colShift = parent.ColShift();
    const int rowShift = parent.RowShift();

    parent.SetToRandom();

    complex<Z>* thisLocalBuffer = parent.LocalBuffer();
    const int thisLDim = parent.LocalLDim();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( int iLocal=0; iLocal<localHeight; ++iLocal )
    {
        const int i = colShift + iLocal*r;                
        if( i % c == rowShift )
        {
            const int jLocal = (i-rowShift) / c;
            if( jLocal < localWidth )
            {
                const Z value = real(thisLocalBuffer[iLocal+jLocal*thisLDim]);
                thisLocalBuffer[iLocal+jLocal*thisLDim] = value;
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
DistMatrix<T,MC,MR>::SetToRandomHPDHelper<Z>::Func
( DistMatrix<Z,MC,MR>& parent )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::SetToRandomHPD");
    parent.AssertNotLockedView();
    if( parent._height != parent._width )
        throw logic_error("Positive-definite matrices must be square");
#endif
    const int r = parent.Grid().Height();
    const int c = parent.Grid().Width();

    const int width = parent._width;
    const int localHeight = parent.LocalHeight();
    const int localWidth = parent.LocalWidth();
    const int colShift = parent.ColShift();
    const int rowShift = parent.RowShift();

    parent.SetToRandom();

    Z* thisLocalBuffer = parent.LocalBuffer();
    const int thisLDim = parent.LocalLDim();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( int iLocal=0; iLocal<localHeight; ++iLocal )
    {
        const int i = colShift + iLocal*r;                
        if( i % c == rowShift )
        {
            const int jLocal = (i-rowShift) / c;
            if( jLocal < localWidth )
                thisLocalBuffer[iLocal+jLocal*thisLDim] += width;
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
DistMatrix<T,MC,MR>::SetToRandomHPDHelper<complex<Z> >::Func
( DistMatrix<complex<Z>,MC,MR>& parent )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::SetToRandomHPD");
    parent.AssertNotLockedView();
    if( parent.Height() != parent.Width() )
        throw logic_error("Positive-definite matrices must be square");
#endif
    const int r = parent.Grid().Height();
    const int c = parent.Grid().Width();

    const int width = parent.Width();
    const int localHeight = parent.LocalHeight();
    const int localWidth = parent.LocalWidth();
    const int colShift = parent.ColShift();
    const int rowShift = parent.RowShift();

    parent.SetToRandom();

    complex<Z>* thisLocalBuffer = parent.LocalBuffer();
    const int thisLDim = parent.LocalLDim();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( int iLocal=0; iLocal<localHeight; ++iLocal )
    {
        const int i = colShift + iLocal*r;                
        if( i % c == rowShift )
        {
            const int jLocal = (i-rowShift) / c;
            if( jLocal < localWidth )
            {
                const Z value = real(thisLocalBuffer[iLocal+jLocal*thisLDim]);
                thisLocalBuffer[iLocal+jLocal*thisLDim] = value + width;
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
DistMatrix<T,MC,MR>::GetRealHelper<complex<Z> >::Func
( const DistMatrix<complex<Z>,MC,MR>& parent, int i, int j ) 
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::GetReal");
    parent.AssertValidEntry( i, j );
#endif
    // We will determine the owner of the (i,j) entry and have him Broadcast
    // throughout the entire process grid
    const elemental::Grid& g = parent.Grid();
    const int ownerRow = (i + parent.ColAlignment()) % g.Height();
    const int ownerCol = (j + parent.RowAlignment()) % g.Width();
    const int ownerRank = ownerRow + ownerCol * g.Height();

    Z u;
    if( g.VCRank() == ownerRank )
    {
        const int iLocal = (i-parent.ColShift()) / g.Height();
        const int jLocal = (j-parent.RowShift()) / g.Width();
        u = parent.GetRealLocalEntry(iLocal,jLocal);
    }
    mpi::Broadcast( &u, 1, g.VCToViewingMap(ownerRank), g.ViewingComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T>
template<typename Z>
inline Z
DistMatrix<T,MC,MR>::GetImagHelper<complex<Z> >::Func
( const DistMatrix<complex<Z>,MC,MR>& parent, int i, int j )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::GetImag");
    parent.AssertValidEntry( i, j );
#endif
    // We will determine the owner of the (i,j) entry and have him Broadcast
    // throughout the entire process grid
    const elemental::Grid& g = parent.Grid();
    const int ownerRow = (i + parent.ColAlignment()) % g.Height();
    const int ownerCol = (j + parent.RowAlignment()) % g.Width();
    const int ownerRank = ownerRow + ownerCol * g.Height();

    Z u;
    if( g.VCRank() == ownerRank )
    {
        const int iLocal = (i-parent.ColShift()) / g.Height();
        const int jLocal = (j-parent.RowShift()) / g.Width();
        u = parent.GetImagLocalEntry(iLocal,jLocal);
    }
    mpi::Broadcast( &u, 1, g.VCToViewingMap(ownerRank), g.ViewingComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T>
template<typename Z>
inline void
DistMatrix<T,MC,MR>::SetRealHelper<complex<Z> >::Func
( DistMatrix<complex<Z>,MC,MR>& parent, int i, int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::SetReal");
    parent.AssertValidEntry( i, j );
#endif
    const elemental::Grid& g = parent.Grid();
    const int ownerRow = (i + parent.ColAlignment()) % g.Height();
    const int ownerCol = (j + parent.RowAlignment()) % g.Width();
    const int ownerRank = ownerRow + ownerCol * g.Height();

    if( g.VCRank() == ownerRank )
    {
        const int iLocal = (i-parent.ColShift()) / g.Height();
        const int jLocal = (j-parent.RowShift()) / g.Width();
        parent.SetRealLocalEntry(iLocal,jLocal,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
template<typename Z>
inline void
DistMatrix<T,MC,MR>::SetImagHelper<complex<Z> >::Func
( DistMatrix<complex<Z>,MC,MR>& parent, int i, int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::SetImag");
    parent.AssertValidEntry( i, j );
#endif
    const elemental::Grid& g = parent.Grid();
    const int ownerRow = (i + parent.ColAlignment()) % g.Height();
    const int ownerCol = (j + parent.RowAlignment()) % g.Width();
    const int ownerRank = ownerRow + ownerCol * g.Height();

    if( g.VCRank() == ownerRank )
    {
        const int iLocal = (i-parent.ColShift()) / g.Height();
        const int jLocal = (j-parent.RowShift()) / g.Width();
        parent.SetImagLocalEntry(iLocal,jLocal,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
template<typename Z>
inline void
DistMatrix<T,MC,MR>::UpdateRealHelper<complex<Z> >::Func
( DistMatrix<complex<Z>,MC,MR>& parent, int i, int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::UpdateReal");
    parent.AssertValidEntry( i, j );
#endif
    const elemental::Grid& g = parent.Grid();
    const int ownerRow = (i + parent.ColAlignment()) % g.Height();
    const int ownerCol = (j + parent.RowAlignment()) % g.Width();
    const int ownerRank = ownerRow + ownerCol * g.Height();

    if( g.VCRank() == ownerRank )
    {
        const int iLocal = (i-parent.ColShift()) / g.Height();
        const int jLocal = (j-parent.RowShift()) / g.Width();
        parent.UpdateRealLocalEntry(iLocal,jLocal,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
template<typename Z>
inline void
DistMatrix<T,MC,MR>::UpdateImagHelper<complex<Z> >::Func
( DistMatrix<complex<Z>,MC,MR>& parent, int i, int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::UpdateImag");
    parent.AssertValidEntry( i, j );
#endif
    const elemental::Grid& g = parent.Grid();
    const int ownerRow = (i + parent.ColAlignment()) % g.Height();
    const int ownerCol = (j + parent.RowAlignment()) % g.Width();
    const int ownerRank = ownerRow + ownerCol * g.Height();

    if( g.VCRank() == ownerRank )
    {
        const int iLocal = (i-parent.ColShift()) / g.Height();
        const int jLocal = (j-parent.RowShift()) / g.Width();
        parent.UpdateImagLocalEntry(iLocal,jLocal,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
template<typename Z>
inline void
DistMatrix<T,MC,MR>::GetRealDiagonalHelper<complex<Z> >::Func
( const DistMatrix<complex<Z>,MC,MR>& parent, 
        DistMatrix<Z,MD,STAR>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::GetRealDiagonal");
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
#endif
    const elemental::Grid& g = parent.Grid();
    if( !d.Viewing() )
    {
        d.SetGrid( g );
        if( !d.ConstrainedColAlignment() )
            d.AlignWithDiagonal( parent, offset );
        d.ResizeTo( length, 1 );
    }

    if( d.InDiagonal() )
    {
        const int r = g.Height();
        const int c = g.Width();
        const int lcm = g.LCM();
        const int colShift = parent.ColShift();
        const int rowShift = parent.RowShift();
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
        const int jLocalStart = (jStart-rowShift) / c;

        const int localDiagLength = d.LocalHeight();

        const complex<Z>* thisLocalBuffer = parent.LockedLocalBuffer();
        const int thisLDim = parent.LocalLDim();
        Z* dLocalBuffer = d.LocalBuffer();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<localDiagLength; ++k )
        {
            const int iLocal = iLocalStart + k*(lcm/r);
            const int jLocal = jLocalStart + k*(lcm/c);
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
DistMatrix<T,MC,MR>::GetImagDiagonalHelper<complex<Z> >::Func
( const DistMatrix<complex<Z>,MC,MR>& parent,
        DistMatrix<Z,MD,STAR>& d, int offset ) 
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::GetImagDiagonal");
    if( d.Viewing() )
        parent.AssertSameGrid( d );
#endif
    int width = parent.Width();
    int height = parent.Height();
    int length;
    if( offset > 0 )
    {
        const int remainingWidth = max(width-offset,0);
        length = min(height,remainingWidth);
    }
    else
    {
        const int remainingHeight = max(height+offset,0);
        length = min(remainingHeight,width);
    }
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
#endif

    const elemental::Grid& g = parent.Grid();
    if( !d.Viewing() )
    {
        d.SetGrid( g );
        if( !d.ConstrainedColAlignment() )
            d.AlignWithDiagonal( parent, offset );
        d.ResizeTo( length, 1 );
    }

    if( d.InDiagonal() )
    {
        const int r = g.Height();
        const int c = g.Width();
        const int lcm = g.LCM();
        const int colShift = parent.ColShift();
        const int rowShift = parent.RowShift();
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
        const int jLocalStart = (jStart-rowShift) / c;

        const int localDiagLength = d.LocalHeight();

        const complex<Z>* thisLocalBuffer = parent.LockedLocalBuffer();
        const int thisLDim = parent.LocalLDim();
        Z* dLocalBuffer = d.LocalBuffer();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<localDiagLength; ++k )
        {
            const int iLocal = iLocalStart + k*(lcm/r);
            const int jLocal = jLocalStart + k*(lcm/c);
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
DistMatrix<T,MC,MR>::GetRealDiagonalHelper<complex<Z> >::Func
( const DistMatrix<complex<Z>,MC,MR>& parent,
        DistMatrix<Z,STAR,MD>& d, int offset ) 
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::GetRealDiagonal");
    if( d.Viewing() )
        parent.AssertSameGrid( d );
#endif
    int height = parent.Height();
    int width = parent.Width();
    int length;
    if( offset > 0 )
    {
        const int remainingWidth = max(width-offset,0);
        length = min(height,remainingWidth);
    }
    else
    {
        const int remainingHeight = max(height+offset,0);
        length = min(remainingHeight,width);
    }
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
#endif

    const elemental::Grid& g = parent.Grid();
    if( !d.Viewing() )
    {
        d.SetGrid( g );
        if( !d.ConstrainedRowAlignment() )
            d.AlignWithDiagonal( parent, offset );
        d.ResizeTo( 1, length );
    }

    if( d.InDiagonal() )
    {
        const int r = g.Height();
        const int c = g.Width();
        const int lcm = g.LCM();
        const int colShift = parent.ColShift();
        const int rowShift = parent.RowShift();
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
        const int jLocalStart = (jStart-rowShift) / c;

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
            const int iLocal = iLocalStart + k*(lcm/r);
            const int jLocal = jLocalStart + k*(lcm/c);
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
DistMatrix<T,MC,MR>::GetImagDiagonalHelper<complex<Z> >::Func
( const DistMatrix<complex<Z>,MC,MR>& parent,
        DistMatrix<Z,STAR,MD>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::GetImagDiagonal");
    if( d.Viewing() )
        parent.AssertSameGrid( d );
#endif
    int height = parent.Height();
    int width = parent.Width();
    int length;
    if( offset > 0 )
    {
        const int remainingWidth = max(width-offset,0);
        length = min(height,remainingWidth);
    }
    else
    {
        const int remainingHeight = max(height+offset,0);
        length = min(remainingHeight,width);
    }
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
#endif

    const elemental::Grid& g = parent.Grid();
    if( !d.Viewing() )
    {
        d.SetGrid( g );
        if( !d.ConstrainedRowAlignment() )
            d.AlignWithDiagonal( parent, offset );
        d.ResizeTo( 1, length );
    }

    if( d.InDiagonal() )
    {
        const int r = g.Height();
        const int c = g.Width();
        const int lcm = g.LCM();
        const int colShift = parent.ColShift();
        const int rowShift = parent.RowShift();
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
        const int jLocalStart = (jStart-rowShift) / c;

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
            const int iLocal = iLocalStart + k*(lcm/r);
            const int jLocal = jLocalStart + k*(lcm/c);
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
DistMatrix<T,MC,MR>::SetRealDiagonalHelper<complex<Z> >::Func
(       DistMatrix<complex<Z>,MC,MR>& parent,
  const DistMatrix<Z,MD,STAR>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::SetRealDiagonal");
    parent.AssertSameGrid( d );
    if( d.Width() != 1 )
        throw logic_error("d must be a column vector");
    const int length = parent.DiagonalLength(offset);
    if( length != d.Height() )
    {
        ostringstream msg;
        msg << "d is not of the same length as the diagonal:\n"
            << "  A ~ " << parent.Height() << " x " << parent.Width() << "\n"
            << "  d ~ " << d.Height() << " x " << d.Width() << "\n"
            << "  A diag length: " << length << "\n";
        throw logic_error( msg.str().c_str() );
    }
#endif
    if( d.InDiagonal() )
    {
        const elemental::Grid& g = parent.Grid();
        const int r = g.Height();
        const int c = g.Width();
        const int lcm = g.LCM();
        const int colShift = parent.ColShift();
        const int rowShift = parent.RowShift();
        const int diagShift = d.ColShift();

        int iStart,jStart;
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
        const int jLocalStart = (jStart-rowShift) / c;

        const int localDiagLength = d.LocalHeight();

        const Z* dLocalBuffer = d.LockedLocalBuffer(); 
        complex<Z>* thisLocalBuffer = parent.LocalBuffer();
        const int thisLDim = parent.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<localDiagLength; ++k )
        {
            const int iLocal = iLocalStart + k*(lcm/r);
            const int jLocal = jLocalStart + k*(lcm/c);
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
DistMatrix<T,MC,MR>::SetImagDiagonalHelper<complex<Z> >::Func
(       DistMatrix<complex<Z>,MC,MR>& parent,
  const DistMatrix<Z,MD,STAR>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::SetImagDiagonal");
    parent.AssertSameGrid( d );
    if( d.Width() != 1 )
        throw logic_error("d must be a column vector");
    const int length = parent.DiagonalLength(offset);
    if( length != d.Height() )
    {
        ostringstream msg;
        msg << "d is not of the same length as the diagonal:\n"
            << "  A ~ " << parent.Height() << " x " << parent.Width() << "\n"
            << "  d ~ " << d.Height() << " x " << d.Width() << "\n"
            << "  A diag length: " << length << "\n";
        throw logic_error( msg.str().c_str() );
    }
#endif
    if( d.InDiagonal() )
    {
        const elemental::Grid& g = parent.Grid();
        const int r = g.Height();
        const int c = g.Width();
        const int lcm = g.LCM();
        const int colShift = parent.ColShift();
        const int rowShift = parent.RowShift();
        const int diagShift = d.ColShift();

        int iStart,jStart;
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
        const int jLocalStart = (jStart-rowShift) / c;

        const int localDiagLength = d.LocalHeight();

        const Z* dLocalBuffer = d.LockedLocalBuffer();
        complex<Z>* thisLocalBuffer = parent.LocalBuffer();
        const int thisLDim = parent.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<localDiagLength; ++k )
        {
            const int iLocal = iLocalStart + k*(lcm/r);
            const int jLocal = jLocalStart + k*(lcm/c);
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
DistMatrix<T,MC,MR>::SetRealDiagonalHelper<complex<Z> >::Func
(       DistMatrix<complex<Z>,MC,MR>& parent,
  const DistMatrix<Z,STAR,MD>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::SetRealDiagonal");
    parent.AssertSameGrid( d );
    if( d.Height() != 1 )
        throw logic_error("d must be a row vector");
    const int length = parent.DiagonalLength(offset);
    if( length != d.Width() )
    {
        ostringstream msg;
        msg << "d is not of the same length as the diagonal:\n"
            << "  A ~ " << parent.Height() << " x " << parent.Width() << "\n"
            << "  d ~ " << d.Height() << " x " << d.Width() << "\n"
            << "  A diag length: " << length << "\n";
        throw logic_error( msg.str().c_str() );
    }
#endif
    if( d.InDiagonal() )
    {
        const elemental::Grid& g = parent.Grid();
        const int r = g.Height();
        const int c = g.Width();
        const int lcm = g.LCM();
        const int colShift = parent.ColShift();
        const int rowShift = parent.RowShift();
        const int diagShift = d.RowShift();

        int iStart,jStart;
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
        const int jLocalStart = (jStart-rowShift) / c;

        const int localDiagLength = d.LocalWidth();

        const Z* dLocalBuffer = d.LockedLocalBuffer();
        const int dLDim = d.LocalLDim();
        complex<Z>* thisLocalBuffer = parent.LocalBuffer();
        const int thisLDim = parent.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<localDiagLength; ++k )
        {
            const int iLocal = iLocalStart + k*(lcm/r);
            const int jLocal = jLocalStart + k*(lcm/c);
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
DistMatrix<T,MC,MR>::SetImagDiagonalHelper<complex<Z> >::Func
(       DistMatrix<complex<Z>,MC,MR>& parent,
  const DistMatrix<Z,STAR,MD>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::SetImagDiagonal");
    parent.AssertSameGrid( d );
    if( d.Height() != 1 )
        throw logic_error("d must be a row vector");
    const int length = parent.DiagonalLength(offset);
    if( length != d.Width() )
    {
        ostringstream msg;
        msg << "d is not of the same length as the diagonal:\n"
            << "  A ~ " << parent.Height() << " x " << parent.Width() << "\n"
            << "  d ~ " << d.Height() << " x " << d.Width() << "\n"
            << "  A diag length: " << length << "\n";
        throw logic_error( msg.str().c_str() );
    }
#endif
    if( d.InDiagonal() )
    {
        const elemental::Grid& g = parent.Grid();
        const int r = g.Height();
        const int c = g.Width();
        const int lcm = g.LCM();
        const int colShift = parent.ColShift();
        const int rowShift = parent.RowShift();
        const int diagShift = d.RowShift();

        int iStart,jStart;
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
        const int jLocalStart = (jStart-rowShift) / c;

        const int localDiagLength = d.LocalWidth();

        const Z* dLocalBuffer = d.LockedLocalBuffer();
        const int dLDim = d.LocalLDim();
        complex<Z>* thisLocalBuffer = parent.LocalBuffer();
        const int thisLDim = parent.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<localDiagLength; ++k )
        {
            const int iLocal = iLocalStart + k*(lcm/r);
            const int jLocal = jLocalStart + k*(lcm/c);
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
