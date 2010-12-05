/*
   Copyright (c) 2009-2010, Jack Poulson
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
using namespace elemental::utilities;
using namespace elemental::wrappers::mpi;

template<typename R>
void
elemental::DistMatrix<R,MC,MR>::SetToRandomHPD()
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::SetToRandomHPD");
    this->AssertNotLockedView();
    if( this->Height() != this->Width() )
        throw logic_error( "Positive-definite matrices must be square." );
#endif
    const int r = this->Grid().Height();
    const int c = this->Grid().Width();

    const int localHeight = this->LocalHeight();
    const int localWidth = this->LocalWidth();
    const int colShift = this->ColShift();
    const int rowShift = this->RowShift();

    this->SetToRandom();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const int i = colShift + iLoc*r;                
        if( i % c == rowShift )
        {
            const int jLoc = (i-rowShift) / c;
            if( jLoc < localWidth )
            {
                const R value = this->GetLocalEntry(iLoc,jLoc);
                this->SetLocalEntry(iLoc,jLoc,value+this->Width());
            }
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
template<typename R>
void
elemental::DistMatrix<complex<R>,MC,MR>::SetToRandomHPD()
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::SetToRandomHPD");
    this->AssertNotLockedView();
    if( this->Height() != this->Width() )
        throw logic_error( "Positive-definite matrices must be square." );
#endif
    const int r = this->Grid().Height();
    const int c = this->Grid().Width();

    const int localHeight = this->LocalHeight();
    const int localWidth = this->LocalWidth();
    const int colShift = this->ColShift();
    const int rowShift = this->RowShift();

    this->SetToRandom();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const int i = colShift + iLoc*r;                
        if( i % c == rowShift )
        {
            const int jLoc = (i-rowShift) / c;
            if( jLoc < localWidth )
            {
                const R value = real(this->GetLocalEntry(iLoc,jLoc));
                this->SetLocalEntry(iLoc,jLoc,value+this->Width());
            }
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
R
elemental::DistMatrix<complex<R>,MC,MR>::GetReal
( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::GetReal");
    this->AssertValidEntry( i, j );
#endif
    // We will determine the owner of the (i,j) entry and have him Broadcast
    // throughout the entire process grid
    const Grid& g = this->Grid();
    const int ownerRow = (i + this->ColAlignment()) % g.Height();
    const int ownerCol = (j + this->RowAlignment()) % g.Width();
    const int ownerRank = ownerRow + ownerCol * g.Height();

    R u;
    if( g.VCRank() == ownerRank )
    {
        const int iLoc = (i-this->ColShift()) / g.Height();
        const int jLoc = (j-this->RowShift()) / g.Width();
        u = real(this->GetLocalEntry(iLoc,jLoc));
    }
    Broadcast( &u, 1, ownerRank, g.VCComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename R>
R
elemental::DistMatrix<complex<R>,MC,MR>::GetImag
( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::GetImag");
    this->AssertValidEntry( i, j );
#endif
    // We will determine the owner of the (i,j) entry and have him Broadcast
    // throughout the entire process grid
    const Grid& g = this->Grid();
    const int ownerRow = (i + this->ColAlignment()) % g.Height();
    const int ownerCol = (j + this->RowAlignment()) % g.Width();
    const int ownerRank = ownerRow + ownerCol * g.Height();

    R u;
    if( g.VCRank() == ownerRank )
    {
        const int iLoc = (i-this->ColShift()) / g.Height();
        const int jLoc = (j-this->RowShift()) / g.Width();
        u = imag(this->GetLocalEntry(iLoc,jLoc));
    }
    Broadcast( &u, 1, ownerRank, g.VCComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename R>
void
elemental::DistMatrix<complex<R>,MC,MR>::SetReal
( int i, int j, R u )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::SetReal");
    this->AssertValidEntry( i, j );
#endif
    const Grid& g = this->Grid();
    const int ownerRow = (i + this->ColAlignment()) % g.Height();
    const int ownerCol = (j + this->RowAlignment()) % g.Width();
    const int ownerRank = ownerRow + ownerCol * g.Height();

    if( g.VCRank() == ownerRank )
    {
        const int iLoc = (i-this->ColShift()) / g.Height();
        const int jLoc = (j-this->RowShift()) / g.Width();
        const R v = imag(this->GetLocalEntry(iLoc,jLoc));
        this->SetLocalEntry(iLoc,jLoc,complex<R>(u,v));
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void
elemental::DistMatrix<complex<R>,MC,MR>::SetImag
( int i, int j, R v )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::SetImag");
    this->AssertValidEntry( i, j );
#endif
    const Grid& g = this->Grid();
    const int ownerRow = (i + this->ColAlignment()) % g.Height();
    const int ownerCol = (j + this->RowAlignment()) % g.Width();
    const int ownerRank = ownerRow + ownerCol * g.Height();

    if( g.VCRank() == ownerRank )
    {
        const int iLoc = (i-this->ColShift()) / g.Height();
        const int jLoc = (j-this->RowShift()) / g.Width();
        R u = real(this->GetLocalEntry(iLoc,jLoc));
        this->SetLocalEntry(iLoc,jLoc,complex<R>(u,v));
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void
elemental::DistMatrix<complex<R>,MC,MR>::GetRealDiagonal
( DistMatrix<R,MD,Star>& d, int offset ) const
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::GetRealDiagonal");
    this->AssertNotLockedView();
#endif
    int width = this->Width();
    int height = this->Height();
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
        msg << "d is not of the same length as the diagonal:" << endl
            << "  A ~ " << this->Height() << " x " << this->Width() << endl
            << "  d ~ " << d.Height() << " x " << d.Width() << endl
            << "  A diag length: " << length << endl;
        throw logic_error( msg.str() );
    }
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
        const Grid& g = this->Grid();
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

        const int iLocStart = (iStart-colShift) / r;
        const int jLocStart = (jStart-rowShift) / c;

        const int localDiagLength = d.LocalHeight();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<localDiagLength; ++k )
        {
            const R u = real(this->GetLocalEntry(iLocStart+k*(lcm/r),
                                                 jLocStart+k*(lcm/c)));
            d.SetLocalEntry(k,0,u);
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void
elemental::DistMatrix<complex<R>,MC,MR>::GetImagDiagonal
( DistMatrix<R,MD,Star>& d, int offset ) const
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::GetImagDiagonal");
    this->AssertNotLockedView();
#endif
    int width = this->Width();
    int height = this->Height();
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
        msg << "d is not of the same length as the diagonal:" << endl
            << "  A ~ " << this->Height() << " x " << this->Width() << endl
            << "  d ~ " << d.Height() << " x " << d.Width() << endl
            << "  A diag length: " << length << endl;
        throw logic_error( msg.str() );
    }
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
        const Grid& g = this->Grid();
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

        const int iLocStart = (iStart-colShift) / r;
        const int jLocStart = (jStart-rowShift) / c;

        const int localDiagLength = d.LocalHeight();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<localDiagLength; ++k )
        {
            const R v = imag(this->GetLocalEntry(iLocStart+k*(lcm/r),
                                                 jLocStart+k*(lcm/c)));
            d.SetLocalEntry(k,0,v);
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void
elemental::DistMatrix<complex<R>,MC,MR>::GetRealDiagonal
( DistMatrix<R,Star,MD>& d, int offset ) const
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::GetRealDiagonal");
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
    {
        ostringstream msg;
        msg << "d is not of the same length as the diagonal:" << endl
            << "  A ~ " << this->Height() << " x " << this->Width() << endl
            << "  d ~ " << d.Height() << " x " << d.Width() << endl
            << "  A diag length: " << length << endl;
        throw logic_error( msg.str() );
    }
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
        const Grid& g = this->Grid();
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

        const int iLocStart = (iStart-colShift) / r;
        const int jLocStart = (jStart-rowShift) / c;

        const int localDiagLength = d.LocalWidth();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<localDiagLength; ++k )
        {
            const R u = real(this->GetLocalEntry(iLocStart+k*(lcm/r),
                                                 jLocStart+k*(lcm/c)));
            d.SetLocalEntry(0,k,u);
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void
elemental::DistMatrix<complex<R>,MC,MR>::GetImagDiagonal
( DistMatrix<R,Star,MD>& d, int offset ) const
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::GetImagDiagonal");
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
    {
        ostringstream msg;
        msg << "d is not of the same length as the diagonal:" << endl
            << "  A ~ " << this->Height() << " x " << this->Width() << endl
            << "  d ~ " << d.Height() << " x " << d.Width() << endl
            << "  A diag length: " << length << endl;
        throw logic_error( msg.str() );
    }
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
        const Grid& g = this->Grid();
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

        const int iLocStart = (iStart-colShift) / r;
        const int jLocStart = (jStart-rowShift) / c;

        const int localDiagLength = d.LocalWidth();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<localDiagLength; ++k )
        {
            const R v = imag(this->GetLocalEntry(iLocStart+k*(lcm/r),
                                                 jLocStart+k*(lcm/c)));
            d.SetLocalEntry(0,k,v);
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void
elemental::DistMatrix<complex<R>,MC,MR>::SetDiagonal
( const DistMatrixBase<R,MD,Star>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::SetDiagonal");
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
        const Grid& g = this->Grid();
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

        const int iLocStart = (iStart-colShift) / r;
        const int jLocStart = (jStart-rowShift) / c;

        const int localDiagLength = d.LocalHeight();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<localDiagLength; ++k )
            this->SetLocalEntry
                (iLocStart+k*(lcm/r),jLocStart+k*(lcm/c),d.GetLocalEntry(k,0));
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void
elemental::DistMatrix<complex<R>,MC,MR>::SetRealDiagonal
( const DistMatrixBase<R,MD,Star>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::SetRealDiagonal");
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
        const Grid& g = this->Grid();
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

        const int iLocStart = (iStart-colShift) / r;
        const int jLocStart = (jStart-rowShift) / c;

        const int localDiagLength = d.LocalHeight();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<localDiagLength; ++k )
        {
            const R v = imag(this->GetLocalEntry(iLocStart+k*(lcm/r),
                                                 jLocStart+k*(lcm/c)));
            this->SetLocalEntry
                (iLocStart+k*(lcm/r),jLocStart+k*(lcm/c),
                 complex<R>(d.GetLocalEntry(k,0),v));
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void
elemental::DistMatrix<complex<R>,MC,MR>::SetImagDiagonal
( const DistMatrixBase<R,MD,Star>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::SetImagDiagonal");
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
        const Grid& g = this->Grid();
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

        const int iLocStart = (iStart-colShift) / r;
        const int jLocStart = (jStart-rowShift) / c;

        const int localDiagLength = d.LocalHeight();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<localDiagLength; ++k )
        {
            const R u = real(this->GetLocalEntry(iLocStart+k*(lcm/r),
                                                 jLocStart+k*(lcm/c)));
            this->SetLocalEntry
                (iLocStart+k*(lcm/r),jLocStart+k*(lcm/c),
                 complex<R>(u,d.GetLocalEntry(k,0)));
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void
elemental::DistMatrix<complex<R>,MC,MR>::SetDiagonal
( const DistMatrixBase<R,Star,MD>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::SetDiagonal");
    if( d.Height() != 1 )
        throw logic_error( "d must be a row vector." );
    {
        const int height = this->Height();
        const int width = this->Width();

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
        const Grid& g = this->Grid();
        const int r = g.Height();
        const int c = g.Width();
        const int lcm = g.LCM();
        const int colShift = this->ColShift();
        const int rowShift = this->RowShift();
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

        const int iLocStart = (iStart-colShift) / r;
        const int jLocStart = (jStart-rowShift) / c;

        const int localDiagLength = d.LocalWidth();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<localDiagLength; ++k )
            this->SetLocalEntry
                (iLocStart+k*(lcm/r),jLocStart+k*(lcm/c),d.GetLocalEntry(0,k));
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void
elemental::DistMatrix<complex<R>,MC,MR>::SetRealDiagonal
( const DistMatrixBase<R,Star,MD>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::SetRealDiagonal");
    if( d.Height() != 1 )
        throw logic_error( "d must be a row vector." );
    {
        const int height = this->Height();
        const int width = this->Width();

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
        const Grid& g = this->Grid();
        const int r = g.Height();
        const int c = g.Width();
        const int lcm = g.LCM();
        const int colShift = this->ColShift();
        const int rowShift = this->RowShift();
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

        const int iLocStart = (iStart-colShift) / r;
        const int jLocStart = (jStart-rowShift) / c;

        const int localDiagLength = d.LocalWidth();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<localDiagLength; ++k )
        {
            const R v = imag(this->GetLocalEntry(iLocStart+k*(lcm/r),
                                                 jLocStart+k*(lcm/c)));
            this->SetLocalEntry
                (iLocStart+k*(lcm/r),jLocStart+k*(lcm/c),
                 complex<R>(d.GetLocalEntry(0,k),v));
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void
elemental::DistMatrix<complex<R>,MC,MR>::SetImagDiagonal
( const DistMatrixBase<R,Star,MD>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::SetImagDiagonal");
    if( d.Height() != 1 )
        throw logic_error( "d must be a row vector." );
    {
        const int height = this->Height();
        const int width = this->Width();

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
        const Grid& g = this->Grid();
        const int r = g.Height();
        const int c = g.Width();
        const int lcm = g.LCM();
        const int colShift = this->ColShift();
        const int rowShift = this->RowShift();
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

        const int iLocStart = (iStart-colShift) / r;
        const int jLocStart = (jStart-rowShift) / c;

        const int localDiagLength = d.LocalWidth();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int k=0; k<localDiagLength; ++k )
        {
            const R u = real(this->GetLocalEntry(iLocStart+k*(lcm/r),
                                                 jLocStart+k*(lcm/c)));
            this->SetLocalEntry
                (iLocStart+k*(lcm/r),jLocStart+k*(lcm/c),
                 complex<R>(u,d.GetLocalEntry(0,k)));
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX

template class elemental::DistMatrix<int,   MC,MR>;
template class elemental::DistMatrix<float, MC,MR>;
template class elemental::DistMatrix<double,MC,MR>;
#ifndef WITHOUT_COMPLEX
template class elemental::DistMatrix<scomplex,MC,MR>;
template class elemental::DistMatrix<dcomplex,MC,MR>;
#endif

