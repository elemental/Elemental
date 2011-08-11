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
#include "elemental/dist_matrix.hpp"
using namespace std;
using namespace elemental;
using namespace elemental::imports;
using namespace elemental::utilities;

// Template conventions:
//   G: general datatype
//
//   T: any ring, e.g., the (Gaussian) integers and the real/complex numbers
//   Z: representation of a real ring, e.g., the integers or real numbers
//   std::complex<Z>: representation of a complex ring, e.g. Gaussian integers
//                    or complex numbers
//
//   F: representation of real or complex number
//   R: representation of real number
//   std::complex<R>: representation of complex number

template<typename Z>
void
elemental::DistMatrix<Z,MR,MC>::SetToRandomHermitian()
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::SetToRandomHermitian");
    this->AssertNotLockedView();
    if( this->Height() != this->Width() )
        throw logic_error( "Hermitian matrices must be square." );
#endif
    this->SetToRandom();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename Z>
void
elemental::DistMatrix<Z,MR,MC>::SetToRandomHPD()
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::SetToRandomHPD");
    this->AssertNotLockedView();
    if( this->Height() != this->Width() )
        throw logic_error( "Positive-definite matrices must be square." );
#endif
    const elemental::Grid& g = this->Grid();

    const int width = this->Width();
    const int localHeight = this->LocalHeight();
    const int localWidth = this->LocalWidth();
    const int r = g.Height();
    const int c = g.Width();
    const int colShift = this->ColShift();
    const int rowShift = this->RowShift();

    this->SetToRandom();

    Z* thisLocalBuffer = this->LocalBuffer();
    const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( int iLocal=0; iLocal<localHeight; ++iLocal )
    {
        const int i = colShift + iLocal*c;
        if( i % r == rowShift )
        {
            const int jLocal = (i-rowShift) / r;
            if( jLocal < localWidth )
                thisLocalBuffer[iLocal+jLocal*thisLDim] += width;
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
template<typename Z>
void
elemental::DistMatrix<complex<Z>,MR,MC>::SetToRandomHermitian()
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::SetToRandomHermitian");
    this->AssertNotLockedView();
    if( this->Height() != this->Width() )
        throw logic_error( "Hermitian matrices must be square." );
#endif
    const elemental::Grid& g = this->Grid();

    const int localHeight = this->LocalHeight();
    const int localWidth = this->LocalWidth();
    const int r = g.Height();
    const int c = g.Width();
    const int colShift = this->ColShift();
    const int rowShift = this->RowShift();

    this->SetToRandom();

    complex<Z>* thisLocalBuffer = this->LocalBuffer();
    const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( int iLocal=0; iLocal<localHeight; ++iLocal )
    {
        const int i = colShift + iLocal*c;
        if( i % r == rowShift )
        {
            const int jLocal = (i-rowShift) / r;
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

template<typename Z>
void
elemental::DistMatrix<complex<Z>,MR,MC>::SetToRandomHPD()
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::SetToRandomHPD");
    this->AssertNotLockedView();
    if( this->Height() != this->Width() )
        throw logic_error( "Positive-definite matrices must be square." );
#endif
    const elemental::Grid& g = this->Grid();

    const int width = this->Width();
    const int localHeight = this->LocalHeight();
    const int localWidth = this->LocalWidth();
    const int r = g.Height();
    const int c = g.Width();
    const int colShift = this->ColShift();
    const int rowShift = this->RowShift();

    this->SetToRandom();

    complex<Z>* thisLocalBuffer = this->LocalBuffer();
    const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( int iLocal=0; iLocal<localHeight; ++iLocal )
    {
        const int i = colShift + iLocal*c;
        if( i % r == rowShift )
        {
            const int jLocal = (i-rowShift) / r;
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

template<typename Z>
Z
elemental::DistMatrix<complex<Z>,MR,MC>::GetReal
( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::GetReal");
    this->AssertValidEntry( i, j );
#endif
    // We will determine the owner of the (i,j) entry and have him Broadcast
    // throughout the entire process grid
    const elemental::Grid& g = this->Grid();
    const int ownerRow = (j + this->RowAlignment()) % g.Height();
    const int ownerCol = (i + this->ColAlignment()) % g.Width();
    const int ownerRank = ownerRow + ownerCol * g.Height();

    Z u;
    if( g.VCRank() == ownerRank )
    {
        const int iLoc = (i-this->ColShift()) / g.Width();
        const int jLoc = (j-this->RowShift()) / g.Height();
        u = this->GetRealLocalEntry(iLoc,jLoc);
    }
    mpi::Broadcast( &u, 1, ownerRank, g.VCComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename Z>
Z
elemental::DistMatrix<complex<Z>,MR,MC>::GetImag
( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::GetImag");
    this->AssertValidEntry( i, j );
#endif
    // We will determine the owner of the (i,j) entry and have him Broadcast
    // throughout the entire process grid
    const elemental::Grid& g = this->Grid();
    const int ownerRow = (j + this->RowAlignment()) % g.Height();
    const int ownerCol = (i + this->ColAlignment()) % g.Width();
    const int ownerRank = ownerRow + ownerCol * g.Height();

    Z u;
    if( g.VCRank() == ownerRank )
    {
        const int iLoc = (i-this->ColShift()) / g.Width();
        const int jLoc = (j-this->RowShift()) / g.Height();
        u = this->GetImagLocalEntry(iLoc,jLoc);
    }
    mpi::Broadcast( &u, 1, ownerRank, g.VCComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename Z>
void
elemental::DistMatrix<complex<Z>,MR,MC>::SetReal
( int i, int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::SetReal");
    this->AssertValidEntry( i, j );
#endif
    const elemental::Grid& g = this->Grid();
    const int ownerRow = (j + this->RowAlignment()) % g.Height();
    const int ownerCol = (i + this->ColAlignment()) % g.Width();
    const int ownerRank = ownerRow + ownerCol * g.Height();

    if( g.VCRank() == ownerRank )
    {
        const int iLoc = (i-this->ColShift()) / g.Width();
        const int jLoc = (j-this->RowShift()) / g.Height();
        this->SetRealLocalEntry(iLoc,jLoc,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename Z>
void
elemental::DistMatrix<complex<Z>,MR,MC>::SetImag
( int i, int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::SetImag");
    this->AssertValidEntry( i, j );
#endif
    const elemental::Grid& g = this->Grid();
    const int ownerRow = (j + this->RowAlignment()) % g.Height();
    const int ownerCol = (i + this->ColAlignment()) % g.Width();
    const int ownerRank = ownerRow + ownerCol * g.Height();

    if( g.VCRank() == ownerRank )
    {
        const int iLoc = (i-this->ColShift()) / g.Width();
        const int jLoc = (j-this->RowShift()) / g.Height();
        this->SetImagLocalEntry(iLoc,jLoc,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename Z>
void
elemental::DistMatrix<complex<Z>,MR,MC>::UpdateReal
( int i, int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::UpdateReal");
    this->AssertValidEntry( i, j );
#endif
    const elemental::Grid& g = this->Grid();
    const int ownerRow = (j + this->RowAlignment()) % g.Height();
    const int ownerCol = (i + this->ColAlignment()) % g.Width();
    const int ownerRank = ownerRow + ownerCol * g.Height();

    if( g.VCRank() == ownerRank )
    {
        const int iLoc = (i-this->ColShift()) / g.Width();
        const int jLoc = (j-this->RowShift()) / g.Height();
        this->UpdateRealLocalEntry(iLoc,jLoc,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename Z>
void
elemental::DistMatrix<complex<Z>,MR,MC>::UpdateImag
( int i, int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::UpdateImag");
    this->AssertValidEntry( i, j );
#endif
    const elemental::Grid& g = this->Grid();
    const int ownerRow = (j + this->RowAlignment()) % g.Height();
    const int ownerCol = (i + this->ColAlignment()) % g.Width();
    const int ownerRank = ownerRow + ownerCol * g.Height();

    if( g.VCRank() == ownerRank )
    {
        const int iLoc = (i-this->ColShift()) / g.Width();
        const int jLoc = (j-this->RowShift()) / g.Height();
        this->UpdateImagLocalEntry(iLoc,jLoc,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename Z>
void
elemental::DistMatrix<complex<Z>,MR,MC>::GetRealDiagonal
( DistMatrix<Z,MD,STAR>& d, int offset ) const
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::GetRealDiagonal([MD,* ])");
    this->AssertNotLockedView();
#endif
    int height = this->Height();
    int width = this->Width();
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
        throw logic_error( "d is not of the correct length." );
#endif
    if( !d.Viewing() )
    {
        if( !d.ConstrainedColAlignment() )
        {
            d.AlignWithDiag( *this, offset );
        }
        d.ResizeTo( length, 1 );
    }

    if( d.InDiagonal() )
    {
        const elemental::Grid& g = this->Grid();
        const int r = g.Height();
        const int c = g.Width();
        const int lcm = g.LCM();
        const int colShift = this->ColShift();
        const int rowShift = this->RowShift();
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

        const int iLocalStart = (iStart-colShift) / c;
        const int jLocalStart = (jStart-rowShift) / r;

        const int localDiagLength = d.LocalHeight();

        Z* dLocalBuffer = d.LocalBuffer();
        const complex<Z>* thisLocalBuffer = this->LockedLocalBuffer();
        const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<localDiagLength; ++k )
        {
            const int iLocal = iLocalStart + k*(lcm/c);
            const int jLocal = jLocalStart + k*(lcm/r);
            dLocalBuffer[k] = real(thisLocalBuffer[iLocal+jLocal*thisLDim]);
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename Z>
void
elemental::DistMatrix<complex<Z>,MR,MC>::GetImagDiagonal
( DistMatrix<Z,MD,STAR>& d, int offset ) const
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::GetImagDiagonal([MD,* ])");
    this->AssertNotLockedView();
#endif
    int height = this->Height();
    int width = this->Width();
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
        throw logic_error( "d is not of the correct length." );
#endif
    if( !d.Viewing() )
    {
        if( !d.ConstrainedColAlignment() )
        {
            d.AlignWithDiag( *this, offset );
        }
        d.ResizeTo( length, 1 );
    }

    if( d.InDiagonal() )
    {
        const elemental::Grid& g = this->Grid();
        const int r = g.Height();
        const int c = g.Width();
        const int lcm = g.LCM();
        const int colShift = this->ColShift();
        const int rowShift = this->RowShift();
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

        const int iLocalStart = (iStart-colShift) / c;
        const int jLocalStart = (jStart-rowShift) / r;

        const int localDiagLength = d.LocalHeight();

        Z* dLocalBuffer = d.LocalBuffer();
        const complex<Z>* thisLocalBuffer = this->LockedLocalBuffer();
        const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<localDiagLength; ++k )
        {
            const int iLocal = iLocalStart + k*(lcm/c);
            const int jLocal = jLocalStart + k*(lcm/r);
            dLocalBuffer[k] = imag(thisLocalBuffer[iLocal+jLocal*thisLDim]);
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename Z>
void
elemental::DistMatrix<complex<Z>,MR,MC>::GetRealDiagonal
( DistMatrix<Z,STAR,MD>& d, int offset ) const
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::GetRealDiagonal([* ,MD])");
    this->AssertNotLockedView();
#endif
    int height = this->Height();
    int width = this->Width();
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
        throw logic_error( "d is not of the correct length." );
#endif

    if( !d.Viewing() )
    {
        if( !d.ConstrainedRowAlignment() )
        {
            d.AlignWithDiag( *this, offset );
        }
        d.ResizeTo( 1, length );
    }

    if( d.InDiagonal() )
    {
        const elemental::Grid& g = this->Grid();
        const int r = g.Height();
        const int c = g.Width();
        const int lcm = g.LCM();
        const int colShift = this->ColShift();
        const int rowShift = this->RowShift();
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

        const int iLocalStart = (iStart-colShift) / c;
        const int jLocalStart = (jStart-rowShift) / r;

        const int localDiagLength = d.LocalWidth();

        Z* dLocalBuffer = d.LocalBuffer();
        const int dLDim = d.LocalLDim();
        const complex<Z>* thisLocalBuffer = this->LockedLocalBuffer();
        const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<localDiagLength; ++k )
        {
            const int iLocal = iLocalStart + k*(lcm/c);
            const int jLocal = jLocalStart + k*(lcm/r);
            dLocalBuffer[k*dLDim] = 
                real(thisLocalBuffer[iLocal+jLocal*thisLDim]);
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename Z>
void
elemental::DistMatrix<complex<Z>,MR,MC>::GetImagDiagonal
( DistMatrix<Z,STAR,MD>& d, int offset ) const
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::GetImagDiagonal([* ,MD])");
    this->AssertNotLockedView();
#endif
    int height = this->Height();
    int width = this->Width();
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
        throw logic_error( "d is not of the correct length." );
#endif

    if( !d.Viewing() )
    {
        if( !d.ConstrainedRowAlignment() )
        {
            d.AlignWithDiag( *this, offset );
        }
        d.ResizeTo( 1, length );
    }

    if( d.InDiagonal() )
    {
        const elemental::Grid& g = this->Grid();
        const int r = g.Height();
        const int c = g.Width();
        const int lcm = g.LCM();
        const int colShift = this->ColShift();
        const int rowShift = this->RowShift();
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

        const int iLocalStart = (iStart-colShift) / c;
        const int jLocalStart = (jStart-rowShift) / r;

        const int localDiagLength = d.LocalWidth();

        Z* dLocalBuffer = d.LocalBuffer();
        const int dLDim = d.LocalLDim();
        const complex<Z>* thisLocalBuffer = this->LockedLocalBuffer();
        const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<localDiagLength; ++k )
        {
            const int iLocal = iLocalStart + k*(lcm/c);
            const int jLocal = jLocalStart + k*(lcm/r);
            dLocalBuffer[k*dLDim] = 
                imag(thisLocalBuffer[iLocal+jLocal*thisLDim]);
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename Z>
void
elemental::DistMatrix<complex<Z>,MR,MC>::SetDiagonal
( const DistMatrixBase<Z,MD,STAR>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::SetDiagonal([MD,* ])");
    if( d.Width() != 1 )
        throw logic_error( "d must be a column vector." );
    {
        int height = this->Height();
        int width = this->Width();
        int length;
        if( offset >= 0 )
        {
            const int remainingWidth = max(width-offset,0);
            length = min(remainingWidth,height);
        }
        else
        {
            const int remainingHeight = max(height+offset,0);
            length = min(remainingHeight,width);
        }
        if( length != d.Height() )
        {
            ostringstream msg;
            msg << "d is not of the same length as the diagonal:" << endl
                << "  A ~ " << this->Height() << " x " << this->Width() << endl
                << "  d ~ " << d.Height() << " x " << d.Width() << endl
                << "  A diag length: " << length << endl;
            throw logic_error( msg.str() );
        }
    }
#endif
    if( d.InDiagonal() )
    {
        const elemental::Grid& g = this->Grid();
        const int r = g.Height();
        const int c = g.Width();
        const int lcm = g.LCM();
        const int colShift = this->ColShift();
        const int rowShift = this->RowShift();
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

        const int iLocalStart = (iStart-colShift) / c;
        const int jLocalStart = (jStart-rowShift) / r;

        const int localDiagLength = d.LocalHeight();

        const Z* dLocalBuffer = d.LockedLocalBuffer();
        complex<Z>* thisLocalBuffer = this->LocalBuffer();
        const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<localDiagLength; ++k )
        {
            const int iLocal = iLocalStart + k*(lcm/c);
            const int jLocal = jLocalStart + k*(lcm/r);
            thisLocalBuffer[iLocal+jLocal*thisLDim] = dLocalBuffer[k];
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename Z>
void
elemental::DistMatrix<complex<Z>,MR,MC>::SetRealDiagonal
( const DistMatrixBase<Z,MD,STAR>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::SetRealDiagonal([MD,* ])");
    if( d.Width() != 1 )
        throw logic_error( "d must be a column vector." );
    {
        int height = this->Height();
        int width = this->Width();
        int length;
        if( offset >= 0 )
        {
            const int remainingWidth = max(width-offset,0);
            length = min(remainingWidth,height);
        }
        else
        {
            const int remainingHeight = max(height+offset,0);
            length = min(remainingHeight,width);
        }
        if( length != d.Height() )
        {
            ostringstream msg;
            msg << "d is not of the same length as the diagonal:" << endl
                << "  A ~ " << this->Height() << " x " << this->Width() << endl
                << "  d ~ " << d.Height() << " x " << d.Width() << endl
                << "  A diag length: " << length << endl;
            throw logic_error( msg.str() );
        }
    }
#endif
    if( d.InDiagonal() )
    {
        const elemental::Grid& g = this->Grid();
        const int r = g.Height();
        const int c = g.Width();
        const int lcm = g.LCM();
        const int colShift = this->ColShift();
        const int rowShift = this->RowShift();
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

        const int iLocalStart = (iStart-colShift) / c;
        const int jLocalStart = (jStart-rowShift) / r;

        const int localDiagLength = d.LocalHeight();

        complex<Z>* thisLocalBuffer = this->LocalBuffer();
        const int thisLDim = this->LocalLDim();
        const Z* dLocalBuffer = d.LockedLocalBuffer();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<localDiagLength; ++k )
        {
            const int iLocal = iLocalStart + k*(lcm/c);
            const int jLocal = jLocalStart + k*(lcm/r);
            const Z u = dLocalBuffer[k];
            const Z v = imag(thisLocalBuffer[iLocal+jLocal*thisLDim]);
            thisLocalBuffer[iLocal+jLocal*thisLDim] = complex<Z>(u,v);
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename Z>
void
elemental::DistMatrix<complex<Z>,MR,MC>::SetImagDiagonal
( const DistMatrixBase<Z,MD,STAR>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::SetImagDiagonal([MD,* ])");
    if( d.Width() != 1 )
        throw logic_error( "d must be a column vector." );
    {
        int height = this->Height();
        int width = this->Width();
        int length;
        if( offset >= 0 )
        {
            const int remainingWidth = max(width-offset,0);
            length = min(remainingWidth,height);
        }
        else
        {
            const int remainingHeight = max(height+offset,0);
            length = min(remainingHeight,width);
        }
        if( length != d.Height() )
        {
            ostringstream msg;
            msg << "d is not of the same length as the diagonal:" << endl
                << "  A ~ " << this->Height() << " x " << this->Width() << endl
                << "  d ~ " << d.Height() << " x " << d.Width() << endl
                << "  A diag length: " << length << endl;
            throw logic_error( msg.str() );
        }
    }
#endif
    if( d.InDiagonal() )
    {
        const elemental::Grid& g = this->Grid();
        const int r = g.Height();
        const int c = g.Width();
        const int lcm = g.LCM();
        const int colShift = this->ColShift();
        const int rowShift = this->RowShift();
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

        const int iLocalStart = (iStart-colShift) / c;
        const int jLocalStart = (jStart-rowShift) / r;

        const int localDiagLength = d.LocalHeight();

        complex<Z>* thisLocalBuffer = this->LocalBuffer();
        const int thisLDim = this->LocalLDim();
        const Z* dLocalBuffer = d.LockedLocalBuffer();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<localDiagLength; ++k )
        {
            const int iLocal = iLocalStart + k*(lcm/c);
            const int jLocal = jLocalStart + k*(lcm/r);
            const Z u = real(thisLocalBuffer[iLocal+jLocal*thisLDim]);
            const Z v = dLocalBuffer[k];
            thisLocalBuffer[iLocal+jLocal*thisLDim] = complex<Z>(u,v);
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename Z>
void
elemental::DistMatrix<complex<Z>,MR,MC>::SetDiagonal
( const DistMatrixBase<Z,STAR,MD>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::SetDiagonal([* ,MD])");
    if( d.Height() != 1 )
        throw logic_error( "d must be a row vector." );
    {
        int height = this->Height();
        int width = this->Width();
        int length;
        if( offset >= 0 )
        {
            const int remainingWidth = max(width-offset,0);
            length = min(remainingWidth,height);
        }
        else
        {
            const int remainingHeight = max(height+offset,0);
            length = min(remainingHeight,width);
        }
        if( length != d.Width() )
        {
            ostringstream msg;
            msg << "d is not of the same length as the diagonal:" << endl
                << "  A ~ " << this->Height() << " x " << this->Width() << endl
                << "  d ~ " << d.Height() << " x " << d.Width() << endl
                << "  A diag length: " << length << endl;
            throw logic_error( msg.str() );
        }
    }
#endif
    if( d.InDiagonal() )
    {
        const elemental::Grid& g = this->Grid();
        const int r = g.Height();
        const int c = g.Width();
        const int lcm = g.LCM();
        const int colShift = this->ColShift();
        const int rowShift = this->RowShift();
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

        const int iLocalStart = (iStart-colShift) / c;
        const int jLocalStart = (jStart-rowShift) / r;

        const int localDiagLength = d.LocalWidth();

        complex<Z>* thisLocalBuffer = this->LocalBuffer();
        const int thisLDim = this->LocalLDim();
        const Z* dLocalBuffer = d.LockedLocalBuffer();
        const int dLDim = d.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<localDiagLength; ++k )
        {
            const int iLocal = iLocalStart + k*(lcm/c);
            const int jLocal = jLocalStart + k*(lcm/r);
            thisLocalBuffer[iLocal+jLocal*thisLDim] = dLocalBuffer[k*dLDim];
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename Z>
void
elemental::DistMatrix<complex<Z>,MR,MC>::SetRealDiagonal
( const DistMatrixBase<Z,STAR,MD>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::SetRealDiagonal([* ,MD])");
    if( d.Height() != 1 )
        throw logic_error( "d must be a row vector." );
    {
        int height = this->Height();
        int width = this->Width();
        int length;
        if( offset >= 0 )
        {
            const int remainingWidth = max(width-offset,0);
            length = min(remainingWidth,height);
        }
        else
        {
            const int remainingHeight = max(height+offset,0);
            length = min(remainingHeight,width);
        }
        if( length != d.Width() )
        {
            ostringstream msg;
            msg << "d is not of the same length as the diagonal:" << endl
                << "  A ~ " << this->Height() << " x " << this->Width() << endl
                << "  d ~ " << d.Height() << " x " << d.Width() << endl
                << "  A diag length: " << length << endl;
            throw logic_error( msg.str() );
        }
    }
#endif
    if( d.InDiagonal() )
    {
        const elemental::Grid& g = this->Grid();
        const int r = g.Height();
        const int c = g.Width();
        const int lcm = g.LCM();
        const int colShift = this->ColShift();
        const int rowShift = this->RowShift();
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

        const int iLocalStart = (iStart-colShift) / c;
        const int jLocalStart = (jStart-rowShift) / r;

        const int localDiagLength = d.LocalWidth();

        complex<Z>* thisLocalBuffer = this->LocalBuffer();
        const int thisLDim = this->LocalLDim();
        const Z* dLocalBuffer = d.LockedLocalBuffer();
        const int dLDim = d.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<localDiagLength; ++k )
        {
            const int iLocal = iLocalStart + k*(lcm/c);
            const int jLocal = jLocalStart + k*(lcm/r);
            const Z u = dLocalBuffer[k*dLDim];
            const Z v = imag(thisLocalBuffer[iLocal+jLocal*thisLDim]);
            thisLocalBuffer[iLocal+jLocal*thisLDim] = complex<Z>(u,v);
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename Z>
void
elemental::DistMatrix<complex<Z>,MR,MC>::SetImagDiagonal
( const DistMatrixBase<Z,STAR,MD>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::SetImagDiagonal([* ,MD])");
    if( d.Height() != 1 )
        throw logic_error( "d must be a row vector." );
    {
        int height = this->Height();
        int width = this->Width();
        int length;
        if( offset >= 0 )
        {
            const int remainingWidth = max(width-offset,0);
            length = min(remainingWidth,height);
        }
        else
        {
            const int remainingHeight = max(height+offset,0);
            length = min(remainingHeight,width);
        }
        if( length != d.Width() )
        {
            ostringstream msg;
            msg << "d is not of the same length as the diagonal:" << endl
                << "  A ~ " << this->Height() << " x " << this->Width() << endl
                << "  d ~ " << d.Height() << " x " << d.Width() << endl
                << "  A diag length: " << length << endl;
            throw logic_error( msg.str() );
        }
    }
#endif
    if( d.InDiagonal() )
    {
        const elemental::Grid& g = this->Grid();
        const int r = g.Height();
        const int c = g.Width();
        const int lcm = g.LCM();
        const int colShift = this->ColShift();
        const int rowShift = this->RowShift();
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

        const int iLocalStart = (iStart-colShift) / c;
        const int jLocalStart = (jStart-rowShift) / r;

        const int localDiagLength = d.LocalWidth();

        complex<Z>* thisLocalBuffer = this->LocalBuffer();
        const int thisLDim = this->LocalLDim();
        const Z* dLocalBuffer = d.LockedLocalBuffer();
        const int dLDim = d.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<localDiagLength; ++k )
        {
            const int iLocal = iLocalStart + k*(lcm/c);
            const int jLocal = jLocalStart + k*(lcm/r);
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

template class elemental::DistMatrix<int,   MR,MC>;
template class elemental::DistMatrix<float, MR,MC>;
template class elemental::DistMatrix<double,MR,MC>;
#ifndef WITHOUT_COMPLEX
template class elemental::DistMatrix<scomplex,MR,MC>;
template class elemental::DistMatrix<dcomplex,MR,MC>;
#endif

