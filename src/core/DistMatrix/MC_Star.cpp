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
elemental::DistMatrix<R,MC,Star>::SetToRandomHPD()
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::SetToRandomHPD");
    this->AssertNotLockedView();
    if( this->Height() != this->Width() )
        throw logic_error( "Positive-definite matrices must be square." );
#endif
    const int width = this->Width();
    const int localHeight = this->LocalHeight();
    const int r = this->GetGrid().Height();
    const int colShift = this->ColShift();

    this->SetToRandom();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const int i = colShift + iLoc*r;
        if( i < width )
        {
            const R value = this->GetLocalEntry(iLoc,i);
            this->SetLocalEntry(iLoc,i,value+this->Width());
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
template<typename R>
void
elemental::DistMatrix<complex<R>,MC,Star>::SetToRandomHPD()
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::SetToRandomHPD");
    this->AssertNotLockedView();
    if( this->Height() != this->Width() )
        throw logic_error( "Positive-definite matrices must be square." );
#endif
    const int width = this->Width();
    const int localHeight = this->LocalHeight();
    const int r = this->GetGrid().Height();
    const int colShift = this->ColShift();

    this->SetToRandom();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const int i = colShift + iLoc*r;
        if( i < width )
        {
            const R value = real(this->GetLocalEntry(iLoc,i));
            this->SetLocalEntry(iLoc,i,value+this->Width());
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
R
elemental::DistMatrix<complex<R>,MC,Star>::GetReal
( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::GetReal");
    this->AssertValidEntry( i, j );
#endif
    // We will determine the owner row of entry (i,j) and broadcast from that 
    // row within each process column
    const Grid& g = this->GetGrid();
    const int ownerRow = (i + this->ColAlignment()) % g.Height();

    R u;
    if( g.MCRank() == ownerRow )
    {
        const int iLoc = (i-this->ColShift()) / g.Height();
        u = real(this->GetLocalEntry(iLoc,j));
    }
    Broadcast( &u, 1, ownerRow, g.MCComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename R>
R
elemental::DistMatrix<complex<R>,MC,Star>::GetImag
( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::GetImag");
    this->AssertValidEntry( i, j );
#endif
    // We will determine the owner row of entry (i,j) and broadcast from that 
    // row within each process column
    const Grid& g = this->GetGrid();
    const int ownerRow = (i + this->ColAlignment()) % g.Height();

    R u;
    if( g.MCRank() == ownerRow )
    {
        const int iLoc = (i-this->ColShift()) / g.Height();
        u = imag(this->GetLocalEntry(iLoc,j));
    }
    Broadcast( &u, 1, ownerRow, g.MCComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename R>
void
elemental::DistMatrix<complex<R>,MC,Star>::SetReal
( int i, int j, R u )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::SetReal");
    this->AssertValidEntry( i, j );
#endif
    const Grid& g = this->GetGrid();
    const int ownerRow = (i + this->ColAlignment()) % g.Height();

    if( g.MCRank() == ownerRow )
    {
        const int iLoc = (i-this->ColShift()) / g.Height();
        const R v = imag(this->GetLocalEntry(iLoc,j));
        this->SetLocalEntry(iLoc,j,complex<R>(u,v));
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void
elemental::DistMatrix<complex<R>,MC,Star>::SetImag
( int i, int j, R v )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::SetImag");
    this->AssertValidEntry( i, j );
#endif
    const Grid& g = this->GetGrid();
    const int ownerRow = (i + this->ColAlignment()) % g.Height();

    if( g.MCRank() == ownerRow )
    {
        const int iLoc = (i-this->ColShift()) / g.Height();
        const R u = real(this->GetLocalEntry(iLoc,j)); 
        this->SetLocalEntry(iLoc,j,complex<R>(u,v));
    }
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX

template class elemental::DistMatrix<int,   MC,Star>;
template class elemental::DistMatrix<float, MC,Star>;
template class elemental::DistMatrix<double,MC,Star>;
#ifndef WITHOUT_COMPLEX
template class elemental::DistMatrix<scomplex,MC,Star>;
template class elemental::DistMatrix<dcomplex,MC,Star>;
#endif

