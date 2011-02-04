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
using namespace elemental::utilities;
using namespace elemental::imports::mpi;

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
bool
elemental::DistMatrix<Z,MD,Star>::AlignedWithDiag
( const DistMatrixBase<Z,MC,MR>& A, int offset ) const
{ return DMB::AlignedWithDiag( A, offset ); }

template<typename Z>
void
elemental::DistMatrix<Z,MD,Star>::AlignWithDiag
( const DistMatrixBase<Z,MC,MR>& A, int offset )
{ DMB::AlignWithDiag( A, offset ); }

template<typename Z>
bool
elemental::DistMatrix<Z,MD,Star>::AlignedWithDiag
( const DistMatrixBase<Z,MR,MC>& A, int offset ) const
{ return DMB::AlignedWithDiag( A, offset ); }

template<typename Z>
void
elemental::DistMatrix<Z,MD,Star>::AlignWithDiag
( const DistMatrixBase<Z,MR,MC>& A, int offset )
{ DMB::AlignWithDiag( A, offset ); }

template<typename Z>
void
elemental::DistMatrix<Z,MD,Star>::SetToRandomHPD()
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::SetToRandomHPD");
    this->AssertNotLockedView();
    if( this->Height() != this->Width() )
        throw logic_error( "Positive-definite matrices must be square." );
#endif
    this->SetToRandom();

    if( this->InDiagonal() )
    {
        const int width = this->Width();
        const int localHeight = this->LocalHeight();
        const int lcm = this->Grid().LCM();
        const int colShift = this->ColShift();

        Z* thisLocalBuffer = this->LocalBuffer();
        const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int iLocal=0; iLocal<localHeight; ++iLocal )
        {
            const int i = colShift + iLocal*lcm;
            if( i < width )
                thisLocalBuffer[iLocal+i*thisLDim] += width;
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
template<typename Z>
void
elemental::DistMatrix<complex<Z>,MD,Star>::SetToRandomHPD()
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::SetToRandomHPD");
    this->AssertNotLockedView();
    if( this->Height() != this->Width() )
        throw logic_error( "Positive-definite matrices must be square." );
#endif
    this->SetToRandom();

    if( this->InDiagonal() )
    {
        const int width = this->Width();
        const int localHeight = this->LocalHeight();
        const int lcm = this->Grid().LCM();
        const int colShift = this->ColShift();

        complex<Z>* thisLocalBuffer = this->LocalBuffer();
        const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int iLocal=0; iLocal<localHeight; ++iLocal )
        {
            const int i = colShift + iLocal*lcm;
            if( i < width )
            {
                const Z value = real(thisLocalBuffer[iLocal+i*thisLDim]);
                thisLocalBuffer[iLocal+i*thisLDim] = value + width;
            }
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename Z>
Z
elemental::DistMatrix<complex<Z>,MD,Star>::GetReal
( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::GetReal");
    this->AssertValidEntry( i, j );
#endif
    // We will determine the owner of entry (i,j) and broadcast from it
    int ownerRank;
    const Grid& g = this->Grid();
    {
        const int r = g.Height();
        const int c = g.Width();
        const int alignmentRank = this->ColAlignment();
        const int alignmentRow = alignmentRank % r;
        const int alignmentCol = alignmentRank / r;
        const int ownerRow = (alignmentRow + i) % r;
        const int ownerCol = (alignmentCol + i) % c;
        ownerRank = ownerRow + r*ownerCol;
    }

    Z u;
    if( g.VCRank() == ownerRank )
    {
        const int iLoc = (i-this->ColShift()) / g.LCM();
        u = this->GetRealLocalEntry(iLoc,j);
    }
    Broadcast( &u, 1, ownerRank, g.VCComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename Z>
Z
elemental::DistMatrix<complex<Z>,MD,Star>::GetImag
( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::GetImag");
    this->AssertValidEntry( i, j );
#endif
    // We will determine the owner of entry (i,j) and broadcast from it
    int ownerRank;
    const Grid& g = this->Grid();
    {
        const int r = g.Height();
        const int c = g.Width();
        const int alignmentRank = this->ColAlignment();
        const int alignmentRow = alignmentRank % r;
        const int alignmentCol = alignmentRank / r;
        const int ownerRow = (alignmentRow + i) % r;
        const int ownerCol = (alignmentCol + i) % c;
        ownerRank = ownerRow + r*ownerCol;
    }

    Z u;
    if( g.VCRank() == ownerRank )
    {
        const int iLoc = (i-this->ColShift()) / g.LCM();
        u = this->GetImagLocalEntry(iLoc,j);
    }
    Broadcast( &u, 1, ownerRank, g.VCComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename Z>
void
elemental::DistMatrix<complex<Z>,MD,Star>::SetReal
( int i, int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::SetReal");
    this->AssertValidEntry( i, j );
#endif
    int ownerRank;
    const Grid& g = this->Grid();
    {
        const int r = g.Height();
        const int c = g.Width();
        const int alignmentRank = this->ColAlignment();
        const int alignmentRow = alignmentRank % r;
        const int alignmentCol = alignmentRank / r;
        const int ownerRow = (alignmentRow + i) % r;
        const int ownerCol = (alignmentCol + i) % c;
        ownerRank = ownerRow + r*ownerCol;
    }

    if( g.VCRank() == ownerRank )
    {
        const int iLoc = (i-this->ColShift()) / g.LCM();
        this->SetRealLocalEntry(iLoc,j,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename Z>
void
elemental::DistMatrix<complex<Z>,MD,Star>::SetImag
( int i, int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::SetImag");
    this->AssertValidEntry( i, j );
#endif
    int ownerRank;
    const Grid& g = this->Grid();
    {
        const int r = g.Height();
        const int c = g.Width();
        const int alignmentRank = this->ColAlignment();
        const int alignmentRow = alignmentRank % r;
        const int alignmentCol = alignmentRank / r;
        const int ownerRow = (alignmentRow + i) % r;
        const int ownerCol = (alignmentCol + i) % c;
        ownerRank = ownerRow + r*ownerCol;
    }

    if( g.VCRank() == ownerRank )
    {
        const int iLoc = (i-this->ColShift()) / g.LCM();
        this->SetImagLocalEntry(iLoc,j,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename Z>
void
elemental::DistMatrix<complex<Z>,MD,Star>::UpdateReal
( int i, int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::UpdateReal");
    this->AssertValidEntry( i, j );
#endif
    int ownerRank;
    const Grid& g = this->Grid();
    {
        const int r = g.Height();
        const int c = g.Width();
        const int alignmentRank = this->ColAlignment();
        const int alignmentRow = alignmentRank % r;
        const int alignmentCol = alignmentRank / r;
        const int ownerRow = (alignmentRow + i) % r;
        const int ownerCol = (alignmentCol + i) % c;
        ownerRank = ownerRow + r*ownerCol;
    }

    if( g.VCRank() == ownerRank )
    {
        const int iLoc = (i-this->ColShift()) / g.LCM();
        this->UpdateRealLocalEntry(iLoc,j,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename Z>
void
elemental::DistMatrix<complex<Z>,MD,Star>::UpdateImag
( int i, int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::UpdateImag");
    this->AssertValidEntry( i, j );
#endif
    int ownerRank;
    const Grid& g = this->Grid();
    {
        const int r = g.Height();
        const int c = g.Width();
        const int alignmentRank = this->ColAlignment();
        const int alignmentRow = alignmentRank % r;
        const int alignmentCol = alignmentRank / r;
        const int ownerRow = (alignmentRow + i) % r;
        const int ownerCol = (alignmentCol + i) % c;
        ownerRank = ownerRow + r*ownerCol;
    }

    if( g.VCRank() == ownerRank )
    {
        const int iLoc = (i-this->ColShift()) / g.LCM();
        this->UpdateImagLocalEntry(iLoc,j,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename Z>
bool
elemental::DistMatrix<Z,MD,Star>::AlignedWithDiag
( const DistMatrixBase<complex<Z>,MC,MR>& A, int offset ) const
{ 
#ifndef RELEASE
    PushCallStack("[MD,* ]::AlignedWithDiag([MC,MR])");
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
        aligned = ( this->ColAlignment() == ownerRow + r*ownerCol );
    }
    else
    {
        const int ownerRow = (colAlignment-offset) % r;
        const int ownerCol = rowAlignment;
        aligned = ( this->ColAlignment() == ownerRow + r*ownerCol );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return aligned;
}

template<typename Z>
void
elemental::DistMatrix<Z,MD,Star>::AlignWithDiag
( const DistMatrixBase<complex<Z>,MC,MR>& A, int offset )
{ 
#ifndef RELEASE
    PushCallStack("[MD,* ]::AlignWithDiag([MC,MR])");
    this->AssertFreeColAlignment();
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
        this->_colAlignment = ownerRow + r*ownerCol;
        this->_inDiagonal = 
            ( g.DiagPath() == g.DiagPath( this->ColAlignment() ) );
    }
    else
    {
        const int ownerRow = (colAlignment-offset) % r;
        const int ownerCol = rowAlignment;
        this->_colAlignment = ownerRow + r*ownerCol;
        this->_inDiagonal = 
            ( g.DiagPath() == g.DiagPath( this->ColAlignment() ) );
    }
    if( this->InDiagonal() )
    {
        this->_colShift = 
            ( g.DiagPathRank() + lcm - 
              g.DiagPathRank( this->ColAlignment() ) ) % lcm;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename Z>
bool
elemental::DistMatrix<Z,MD,Star>::AlignedWithDiag
( const DistMatrixBase<complex<Z>,MR,MC>& A, int offset ) const
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::AlignedWithDiag([MR,MC])");
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
        aligned = ( this->ColAlignment() == ownerRow + r*ownerCol );
    }
    else
    {
        const int ownerRow = (rowAlignment-offset) % r;
        const int ownerCol = colAlignment;
        aligned = ( this->ColAlignment() == ownerRow + r*ownerCol );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return aligned;
}
template<typename Z>
void
elemental::DistMatrix<Z,MD,Star>::AlignWithDiag
( const DistMatrixBase<complex<Z>,MR,MC>& A, int offset )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::AlignWithDiag([MR,MC])");
    this->AssertFreeColAlignment();
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
        this->_colAlignment = ownerRow + r*ownerCol;
        this->_inDiagonal = 
            ( g.DiagPath() == g.DiagPath( this->ColAlignment() ) );
    }
    else
    {
        const int ownerRow = (rowAlignment-offset) % r;
        const int ownerCol = colAlignment;
        this->_colAlignment = ownerRow + r*ownerCol;
        this->_inDiagonal = 
            ( g.DiagPath() == g.DiagPath( this->ColAlignment() ) );
    }
    if( this->InDiagonal() )
    {
        this->_colShift = 
            ( g.DiagPathRank() + lcm -
              g.DiagPathRank( this->ColAlignment() ) ) % lcm;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX

template class elemental::DistMatrix<int,   MD,Star>;
template class elemental::DistMatrix<float, MD,Star>;
template class elemental::DistMatrix<double,MD,Star>;
#ifndef WITHOUT_COMPLEX
template class elemental::DistMatrix<scomplex,MD,Star>;
template class elemental::DistMatrix<dcomplex,MD,Star>;
#endif

