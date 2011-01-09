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
using namespace elemental::import::mpi;

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
elemental::DistMatrix<Z,Star,Star>::SetToRandomHPD()
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::SetToRandomHPD");
    this->AssertNotLockedView();
    if( this->Height() != this->Width() )
        throw logic_error( "Positive-definite matrices must be square." );
#endif
    if( this->Grid().InGrid() )
    {
        const int height = this->Height();
        const int width = this->Width();

        this->SetToRandom();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int j=0; j<min(height,width); ++j )
        {
            const Z value = this->GetLocalEntry(j,j);
            this->SetLocalEntry(j,j,value+this->Width());
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
template<typename Z>
void
elemental::DistMatrix<complex<Z>,Star,Star>::SetToRandomHPD()
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::SetToRandomHPD");
    this->AssertNotLockedView();
    if( this->Height() != this->Width() )
        throw logic_error( "Positive-definite matrices must be square." );
#endif
    if( this->Grid().InGrid() )
    {
        const int height = this->Height();
        const int width = this->Width();

        this->SetToRandom();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int j=0; j<min(height,width); ++j )
        {
            const Z value = real(this->GetLocalEntry(j,j));
            this->SetLocalEntry(j,j,value+this->Width());
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename Z>
Z
elemental::DistMatrix<complex<Z>,Star,Star>::GetReal
( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::GetReal");
    this->AssertValidEntry( i, j );
#endif
    int viewingSize, owningSize;
    MPI_Comm_size( this->Grid().ViewingComm(), &viewingSize );
    MPI_Group_size( this->Grid().OwningGroup(), &owningSize );
    Z u;
    if( viewingSize == owningSize )
    {
        // Everyone can just grab their own copy of the data
        u = real(this->GetLocalEntry(i,j));
    }
    else
    {
        // Have the root broadcast its data
        if( this->Grid().VCRank() == 0 )
            u = real(this->GetLocalEntry(i,j));
        Broadcast
        ( &u, 1, this->Grid().OwningToViewingMap(0),
          this->Grid().ViewingComm() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename Z>
Z
elemental::DistMatrix<complex<Z>,Star,Star>::GetImag
( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::GetImag");
    this->AssertValidEntry( i, j );
#endif
    int viewingSize, owningSize;
    MPI_Comm_size( this->Grid().ViewingComm(), &viewingSize );
    MPI_Group_size( this->Grid().OwningGroup(), &owningSize );
    Z u;
    if( viewingSize == owningSize )
    {
        // Everyone can just grab their own copy of the data
        u = imag(this->GetLocalEntry(i,j));
    }
    else
    { 
        // Have the root broadcast its data
        if( this->Grid().VCRank() == 0 )
            u = imag(this->GetLocalEntry(i,j));
        Broadcast
        ( &u, 1, this->Grid().OwningToViewingMap(0),
          this->Grid().ViewingComm() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename Z>
void
elemental::DistMatrix<complex<Z>,Star,Star>::SetReal
( int i, int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::SetReal");
    this->AssertValidEntry( i, j );
#endif
    if( this->Grid().InGrid() )
    {
        const Z v = imag(this->GetLocalEntry(i,j));
        this->SetLocalEntry(i,j,complex<Z>(u,v));
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename Z>
void
elemental::DistMatrix<complex<Z>,Star,Star>::SetImag
( int i, int j, Z v )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::SetImag");
    this->AssertValidEntry( i, j );
#endif
    if( this->Grid().InGrid() )
    {
        const Z u = real(this->GetLocalEntry(i,j));
        this->SetLocalEntry(i,j,complex<Z>(u,v));
    }
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX

template class elemental::DistMatrix<int,   Star,Star>;
template class elemental::DistMatrix<float, Star,Star>;
template class elemental::DistMatrix<double,Star,Star>;
#ifndef WITHOUT_COMPLEX
template class elemental::DistMatrix<scomplex,Star,Star>;
template class elemental::DistMatrix<dcomplex,Star,Star>;
#endif

