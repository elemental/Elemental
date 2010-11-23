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
elemental::DistMatrix<R,Star,Star>::SetToRandomHPD()
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::SetToRandomHPD");
    this->AssertNotLockedView();
    if( this->Height() != this->Width() )
        throw logic_error( "Positive-definite matrices must be square." );
#endif
    const int height = this->Height();
    const int width = this->Width();

    this->SetToRandom();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( int j=0; j<min(height,width); ++j )
    {
        const R value = this->GetLocalEntry(j,j);
        this->SetLocalEntry(j,j,value+this->Width());
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
template<typename R>
void
elemental::DistMatrix<complex<R>,Star,Star>::SetToRandomHPD()
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::SetToRandomHPD");
    this->AssertNotLockedView();
    if( this->Height() != this->Width() )
        throw logic_error( "Positive-definite matrices must be square." );
#endif
    const int height = this->Height();
    const int width = this->Width();

    this->SetToRandom();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( int j=0; j<min(height,width); ++j )
    {
        const R value = real(this->GetLocalEntry(j,j));
        this->SetLocalEntry(j,j,value+this->Width());
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
R
elemental::DistMatrix<complex<R>,Star,Star>::GetReal
( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::GetReal");
    this->AssertValidEntry( i, j );
#endif
    R u = real(this->GetLocalEntry(i,j));
#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename R>
R
elemental::DistMatrix<complex<R>,Star,Star>::GetImag
( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::GetImag");
    this->AssertValidEntry( i, j );
#endif
    R u = imag(this->GetLocalEntry(i,j));
#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename R>
void
elemental::DistMatrix<complex<R>,Star,Star>::SetReal
( int i, int j, R u )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::SetReal");
    this->AssertValidEntry( i, j );
#endif
    const R v = imag(this->GetLocalEntry(i,j));
    this->SetLocalEntry(i,j,complex<R>(u,v));
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void
elemental::DistMatrix<complex<R>,Star,Star>::SetImag
( int i, int j, R v )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::SetImag");
    this->AssertValidEntry( i, j );
#endif
    const R u = real(this->GetLocalEntry(i,j));
    this->SetLocalEntry(i,j,complex<R>(u,v));
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

