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
bool
elemental::DistMatrix<R,Star,MD>::AlignedWithDiag
( const DistMatrixBase<R,MC,MR>& A, int offset ) const
{ return DMB::AlignedWithDiag( A, offset ); }

template<typename R>
void
elemental::DistMatrix<R,Star,MD>::AlignWithDiag
( const DistMatrixBase<R,MC,MR>& A, int offset )
{ DMB::AlignWithDiag( A, offset ); }

template<typename R>
bool
elemental::DistMatrix<R,Star,MD>::AlignedWithDiag
( const DistMatrixBase<R,MR,MC>& A, int offset ) const
{ return DMB::AlignedWithDiag( A, offset ); }

template<typename R>
void
elemental::DistMatrix<R,Star,MD>::AlignWithDiag
( const DistMatrixBase<R,MR,MC>& A, int offset )
{ DMB::AlignWithDiag( A, offset ); }

template<typename R>
void
elemental::DistMatrix<R,Star,MD>::SetToRandomHPD()
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::SetToRandomHPD");
    this->AssertNotLockedView();
    if( this->Height() != this->Width() )
        throw logic_error( "Positive-definite matrices must be square." );
#endif
    this->SetToRandom();

    if( this->InDiagonal() )
    {
        const int height = this->Height();
        const int localWidth = this->LocalWidth();
        const int lcm = this->Grid().LCM();
        const int rowShift = this->RowShift();

#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const int j = rowShift + jLoc*lcm;
            if( j < height )
            {
                const R value = this->GetLocalEntry(j,jLoc);
                this->SetLocalEntry(j,jLoc,value+this->Width());
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
elemental::DistMatrix<complex<R>,Star,MD>::SetToRandomHPD()
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::SetToRandomHPD");
    this->AssertNotLockedView();
    if( this->Height() != this->Width() )
        throw logic_error( "Positive-definite matrices must be square." );
#endif
    this->SetToRandom();

    if( this->InDiagonal() )
    {
        const int height = this->Height();
        const int localWidth = this->LocalWidth();
        const int lcm = this->Grid().LCM();
        const int rowShift = this->RowShift();

#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const int j = rowShift + jLoc*lcm;
            if( j < height )
            {
                const R value = real(this->GetLocalEntry(j,jLoc));
                this->SetLocalEntry(j,jLoc,value+this->Width());
            }
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
R
elemental::DistMatrix<complex<R>,Star,MD>::GetReal
( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::GetReal");
    this->AssertValidEntry( i, j );
#endif
    // We will determine the owner of entry (i,j) and broadcast from it
    int ownerRank;
    const Grid& g = this->Grid();
    {
        const int r = g.Height();
        const int c = g.Width();
        const int alignmentRank = this->RowAlignment();
        const int alignmentRow = alignmentRank % r;
        const int alignmentCol = alignmentRank / r;
        const int ownerRow = (alignmentRow + j) % r;
        const int ownerCol = (alignmentCol + j) % c;
        ownerRank = ownerRow + r*ownerCol;
    }

    R u;
    if( g.VCRank() == ownerRank )
    {
        const int jLoc = (j-this->RowShift()) / g.LCM();
        u = real(this->GetLocalEntry(i,jLoc));
    }
    Broadcast( &u, 1, ownerRank, g.VCComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename R>
R
elemental::DistMatrix<complex<R>,Star,MD>::GetImag
( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::GetImag");
    this->AssertValidEntry( i, j );
#endif
    // We will determine the owner of entry (i,j) and broadcast from it
    int ownerRank;
    const Grid& g = this->Grid();
    {
        const int r = g.Height();
        const int c = g.Width();
        const int alignmentRank = this->RowAlignment();
        const int alignmentRow = alignmentRank % r;
        const int alignmentCol = alignmentRank / r;
        const int ownerRow = (alignmentRow + j) % r;
        const int ownerCol = (alignmentCol + j) % c;
        ownerRank = ownerRow + r*ownerCol;
    }

    R u;
    if( g.VCRank() == ownerRank )
    {
        const int jLoc = (j-this->RowShift()) / g.LCM();
        u = imag(this->GetLocalEntry(i,jLoc));
    }
    Broadcast( &u, 1, ownerRank, g.VCComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename R>
void
elemental::DistMatrix<complex<R>,Star,MD>::SetReal
( int i, int j, R u )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::SetReal");
    this->AssertValidEntry( i, j );
#endif
    int ownerRank;
    const Grid& g = this->Grid();
    {
        const int r = g.Height();
        const int c = g.Width();
        const int alignmentRank = this->RowAlignment();
        const int alignmentRow = alignmentRank % r;
        const int alignmentCol = alignmentRank / r;
        const int ownerRow = (alignmentRow + j) % r;
        const int ownerCol = (alignmentCol + j) % c;
        ownerRank = ownerRow + r*ownerCol;
    }

    if( g.VCRank() == ownerRank )
    {
        const int jLoc = (j-this->RowShift()) / g.LCM();
        const R v = imag(this->GetLocalEntry(i,jLoc));
        this->SetLocalEntry(i,jLoc,complex<R>(u,v));
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void
elemental::DistMatrix<complex<R>,Star,MD>::SetImag
( int i, int j, R v )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::SetImag");
    this->AssertValidEntry( i, j );
#endif
    int ownerRank;
    const Grid& g = this->Grid();
    {
        const int r = g.Height();
        const int c = g.Width();
        const int alignmentRank = this->RowAlignment();
        const int alignmentRow = alignmentRank % r;
        const int alignmentCol = alignmentRank / r;
        const int ownerRow = (alignmentRow + j) % r;
        const int ownerCol = (alignmentCol + j) % c;
        ownerRank = ownerRow + r*ownerCol;
    }

    if( g.VCRank() == ownerRank )
    {
        const int jLoc = (j-this->RowShift()) / g.LCM();
        const R u = real(this->GetLocalEntry(i,jLoc));
        this->SetLocalEntry(i,jLoc,complex<R>(u,v));
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
bool
elemental::DistMatrix<R,Star,MD>::AlignedWithDiag
( const DistMatrixBase<complex<R>,MC,MR>& A, int offset ) const
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::AlignWithedDiag([MC,MR])");
    this->AssertSameGrid( A );
#endif
    const Grid& g = this->Grid();
    const int r = g.Height();
    const int c = g.Width();
    const int colAlignment = A.ColAlignment();
    const int rowAlignment = A.RowAlignment();
    bool aligned;

    if( offset >= 0 )
    {
        const int ownerRow = colAlignment;
        const int ownerCol = (rowAlignment + offset) % c;
        aligned = ( this->RowAlignment() == ownerRow + r*ownerCol );
    }
    else
    {
        const int ownerRow = (colAlignment-offset) % r;
        const int ownerCol = rowAlignment;
        aligned = ( this->RowAlignment() == ownerRow + r*ownerCol );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return aligned;
}

template<typename R>
void
elemental::DistMatrix<R,Star,MD>::AlignWithDiag
( const DistMatrixBase<complex<R>,MC,MR>& A, int offset )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::AlignWithDiag([MC,MR])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    const Grid& g = this->Grid();
    const int r = g.Height();
    const int c = g.Width();
    const int lcm = g.LCM();
    const int colAlignment = A.ColAlignment();
    const int rowAlignment = A.RowAlignment();

    if( offset >= 0 )
    {
        const int ownerRow = colAlignment;
        const int ownerCol = (rowAlignment + offset) % c;
        this->_rowAlignment = ownerRow + r*ownerCol;
        this->_inDiagonal = 
            ( g.DiagPath() == g.DiagPath( this->RowAlignment() ) );
    }
    else
    {
        const int ownerRow = (colAlignment-offset) % r;
        const int ownerCol = rowAlignment;
        this->_rowAlignment = ownerRow + r*ownerCol;
        this->_inDiagonal = 
            ( g.DiagPath() == g.DiagPath( this->RowAlignment() ) );
    }
    if( this->InDiagonal() )
    {
        this->_rowShift = 
            ( g.DiagPathRank() + lcm - 
              g.DiagPathRank( this->RowAlignment() ) ) % lcm;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
bool
elemental::DistMatrix<R,Star,MD>::AlignedWithDiag
( const DistMatrixBase<complex<R>,MR,MC>& A, int offset ) const
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::AlignedWithDiag([MR,MC])");
    this->AssertSameGrid( A );
#endif
    const Grid& g = this->Grid();
    const int r = g.Height();
    const int c = g.Width();
    const int colAlignment = A.ColAlignment();
    const int rowAlignment = A.RowAlignment();
    bool aligned;

    if( offset >= 0 )
    {
        const int ownerRow = rowAlignment;
        const int ownerCol = (colAlignment + offset) % c;
        aligned = ( this->RowAlignment() == ownerRow + r*ownerCol );
    }
    else
    {
        const int ownerRow = (rowAlignment-offset) % r;
        const int ownerCol = colAlignment;
        aligned = ( this->RowAlignment() == ownerRow + r*ownerCol );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return aligned;
}

template<typename R>
void
elemental::DistMatrix<R,Star,MD>::AlignWithDiag
( const DistMatrixBase<complex<R>,MR,MC>& A, int offset )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::AlignWithDiag([MR,MC])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    const Grid& g = this->Grid();
    const int r = g.Height();
    const int c = g.Width();
    const int lcm = g.LCM();
    const int colAlignment = A.ColAlignment();
    const int rowAlignment = A.RowAlignment();

    if( offset >= 0 )
    {
        const int ownerRow = rowAlignment;
        const int ownerCol = (colAlignment + offset) % c;
        this->_rowAlignment = ownerRow + r*ownerCol;
        this->_inDiagonal = 
            ( g.DiagPath() == g.DiagPath( this->RowAlignment() ) );
    }
    else
    {
        const int ownerRow = (rowAlignment-offset) % r;
        const int ownerCol = colAlignment;
        this->_rowAlignment = ownerRow + r*ownerCol;
        this->_inDiagonal = 
            ( g.DiagPath() == g.DiagPath( this->RowAlignment() ) );
    }
    if( this->InDiagonal() )
    {
        this->_rowShift = 
            ( g.DiagPathRank() + lcm -
              g.DiagPathRank( this->RowAlignment() ) ) % lcm;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX

template class elemental::DistMatrix<int,   Star,MD>;
template class elemental::DistMatrix<float, Star,MD>;
template class elemental::DistMatrix<double,Star,MD>;
#ifndef WITHOUT_COMPLEX
template class elemental::DistMatrix<scomplex,Star,MD>;
template class elemental::DistMatrix<dcomplex,Star,MD>;
#endif

