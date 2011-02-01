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
void
elemental::DistMatrix<Z,Star,MR>::SetToRandomHPD()
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::SetToRandomHPD");
    this->AssertNotLockedView();
    if( this->Height() != this->Width() )
        throw logic_error( "Positive-definite matrices must be square." );
#endif
    const int height     = this->Height();
    const int localWidth = this->LocalWidth();
    const int c          = this->Grid().Width();
    const int rowShift   = this->RowShift();

    this->SetToRandom();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const int j = rowShift + jLoc*c;
        if( j < height )
        {
            const Z value = this->GetLocalEntry(j,jLoc);
            this->SetLocalEntry(j,jLoc,value+this->Width());
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
template<typename Z>
void
elemental::DistMatrix<complex<Z>,Star,MR>::SetToRandomHPD()
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::SetToRandomHPD");
    this->AssertNotLockedView();
    if( this->Height() != this->Width() )
        throw logic_error( "Positive-definite matrices must be square." );
#endif
    const int height     = this->Height();
    const int localWidth = this->LocalWidth();
    const int c          = this->Grid().Width();
    const int rowShift   = this->RowShift();

    this->SetToRandom();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const int j = rowShift + jLoc*c;
        if( j < height )
        {
            const Z value = real(this->GetLocalEntry(j,jLoc));
            this->SetLocalEntry(j,jLoc,value+this->Width());
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename Z>
Z
elemental::DistMatrix<complex<Z>,Star,MR>::GetReal
( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::GetReal");
    this->AssertValidEntry( i, j );
#endif
    // We will determine the owner column of entry (i,j) and broadcast from that
    // column within each process row
    const Grid& g = this->Grid();
    const int ownerCol = (j + this->RowAlignment()) % g.Width();

    Z u;
    if( g.MRRank() == ownerCol )
    {
        const int jLoc = (j-this->RowShift()) / g.Width();
        u = this->GetRealLocalEntry(i,jLoc);
    }
    Broadcast( &u, 1, ownerCol, g.MRComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename Z>
Z
elemental::DistMatrix<complex<Z>,Star,MR>::GetImag
( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::GetImag");
    this->AssertValidEntry( i, j );
#endif
    // We will determine the owner column of entry (i,j) and broadcast from that
    // column within each process row
    const Grid& g = this->Grid();
    const int ownerCol = (j + this->RowAlignment()) % g.Width();

    Z u;
    if( g.MRRank() == ownerCol )
    {
        const int jLoc = (j-this->RowShift()) / g.Width();
        u = this->GetImagLocalEntry(i,jLoc);
    }
    Broadcast( &u, 1, ownerCol, g.MRComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename Z>
void
elemental::DistMatrix<complex<Z>,Star,MR>::SetReal
( int i, int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::SetReal");
    this->AssertValidEntry( i, j );
#endif
    const Grid& g = this->Grid();
    const int ownerCol = (j + this->RowAlignment()) % g.Width();

    if( g.MRRank() == ownerCol )
    {
        const int jLoc = (j-this->RowShift()) / g.Width();
        this->SetRealLocalEntry(i,jLoc,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename Z>
void
elemental::DistMatrix<complex<Z>,Star,MR>::SetImag
( int i, int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::SetImag");
    this->AssertValidEntry( i, j );
#endif
    const Grid& g = this->Grid();
    const int ownerCol = (j + this->RowAlignment()) % g.Width();

    if( g.MRRank() == ownerCol )
    {
        const int jLoc = (j-this->RowShift()) / g.Width();
        this->SetImagLocalEntry(i,jLoc,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename Z>
void
elemental::DistMatrix<complex<Z>,Star,MR>::UpdateReal
( int i, int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::UpdateReal");
    this->AssertValidEntry( i, j );
#endif
    const Grid& g = this->Grid();
    const int ownerCol = (j + this->RowAlignment()) % g.Width();

    if( g.MRRank() == ownerCol )
    {
        const int jLoc = (j-this->RowShift()) / g.Width();
        this->UpdateRealLocalEntry(i,jLoc,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename Z>
void
elemental::DistMatrix<complex<Z>,Star,MR>::UpdateImag
( int i, int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::UpdateImag");
    this->AssertValidEntry( i, j );
#endif
    const Grid& g = this->Grid();
    const int ownerCol = (j + this->RowAlignment()) % g.Width();

    if( g.MRRank() == ownerCol )
    {
        const int jLoc = (j-this->RowShift()) / g.Width();
        this->UpdateImagLocalEntry(i,jLoc,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX

template class elemental::DistMatrix<int,   Star,MR>;
template class elemental::DistMatrix<float, Star,MR>;
template class elemental::DistMatrix<double,Star,MR>;
#ifndef WITHOUT_COMPLEX
template class elemental::DistMatrix<scomplex,Star,MR>;
template class elemental::DistMatrix<dcomplex,Star,MR>;
#endif

