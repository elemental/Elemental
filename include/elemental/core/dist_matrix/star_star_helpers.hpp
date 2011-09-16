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
DistMatrix<T,STAR,STAR>::SetToRandomHermitianHelper<Z>::Func
( DistMatrix<Z,STAR,STAR>& parent )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::SetToRandomHermitian");
    parent.AssertNotLockedView();
    if( parent.Height() != parent.Width() )
        throw logic_error( "Hermitian matrices must be square." );
#endif
    if( parent.Grid().InGrid() )
        parent.SetToRandom();
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
template<typename T>
template<typename Z>
inline void
DistMatrix<T,STAR,STAR>::SetToRandomHermitianHelper<complex<Z> >::
Func( DistMatrix<complex<Z>,STAR,STAR>& parent )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::SetToRandomHermitian");
    parent.AssertNotLockedView();
    if( parent.Height() != parent.Width() )
        throw logic_error( "Hermitian matrices must be square." );
#endif
    if( parent.Grid().InGrid() )
    {
        const int height = parent.Height();
        const int width = parent.Width();

        parent.SetToRandom();

        complex<Z>* thisLocalBuffer = parent.LocalBuffer();
        const int thisLDim = parent.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int j=0; j<min(height,width); ++j )
        {
            const Z value = real(thisLocalBuffer[j+j*thisLDim]);
            thisLocalBuffer[j+j*thisLDim] = value;
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
DistMatrix<T,STAR,STAR>::SetToRandomHPDHelper<Z>::Func
( DistMatrix<Z,STAR,STAR>& parent )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::SetToRandomHPD");
    parent.AssertNotLockedView();
    if( parent.Height() != parent.Width() )
        throw logic_error( "Positive-definite matrices must be square." );
#endif
    if( parent.Grid().InGrid() )
    {
        const int height = parent.Height();
        const int width = parent.Width();

        parent.SetToRandom();

        Z* thisLocalBuffer = parent.LocalBuffer();
        const int thisLDim = parent.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int j=0; j<min(height,width); ++j )
            thisLocalBuffer[j+j*thisLDim] += width;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
template<typename T>
template<typename Z>
inline void
DistMatrix<T,STAR,STAR>::SetToRandomHPDHelper<complex<Z> >::Func
( DistMatrix<complex<Z>,STAR,STAR>& parent )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::SetToRandomHPD");
    parent.AssertNotLockedView();
    if( parent.Height() != parent.Width() )
        throw logic_error( "Positive-definite matrices must be square." );
#endif
    if( parent.Grid().InGrid() )
    {
        const int height = parent.Height();
        const int width = parent.Width();

        parent.SetToRandom();

        complex<Z>* thisLocalBuffer = parent.LocalBuffer();
        const int thisLDim = parent.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( int j=0; j<min(height,width); ++j )
        {
            const Z value = real(thisLocalBuffer[j+j*thisLDim]);
            thisLocalBuffer[j+j*thisLDim] = value + width;
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
template<typename Z>
inline Z
DistMatrix<T,STAR,STAR>::GetRealHelper<complex<Z> >::Func
( const DistMatrix<complex<Z>,STAR,STAR>& parent, int i, int j ) 
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::GetReal");
    parent.AssertValidEntry( i, j );
#endif
    const int viewingSize = mpi::CommSize( parent.Grid().ViewingComm() );
    const int owningSize = mpi::GroupSize( parent.Grid().OwningGroup() );
    Z u;
    if( viewingSize == owningSize )
    {
        // Everyone can just grab their own copy of the data
        u = parent.GetRealLocalEntry(i,j);
    }
    else
    {
        // Have the root broadcast its data
        if( parent.Grid().VCRank() == 0 )
            u = parent.GetRealLocalEntry(i,j);
        mpi::Broadcast
        ( &u, 1, parent.Grid().VCToViewingMap(0),
          parent.Grid().ViewingComm() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T>
template<typename Z>
inline Z
DistMatrix<T,STAR,STAR>::GetImagHelper<complex<Z> >::Func
( const DistMatrix<complex<Z>,STAR,STAR>& parent, int i, int j ) 
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::GetImag");
    parent.AssertValidEntry( i, j );
#endif
    const int viewingSize = mpi::CommSize( parent.Grid().ViewingComm() );
    const int owningSize = mpi::GroupSize( parent.Grid().OwningGroup() );
    Z u;
    if( viewingSize == owningSize )
    {
        // Everyone can just grab their own copy of the data
        u = parent.GetImagLocalEntry(i,j);
    }
    else
    { 
        // Have the root broadcast its data
        if( parent.Grid().VCRank() == 0 )
            u = parent.GetImagLocalEntry(i,j);
        mpi::Broadcast
        ( &u, 1, parent.Grid().VCToViewingMap(0),
          parent.Grid().ViewingComm() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename T>
template<typename Z>
inline void
DistMatrix<T,STAR,STAR>::SetRealHelper<complex<Z> >::Func
( DistMatrix<complex<Z>,STAR,STAR>& parent, int i, int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::SetReal");
    parent.AssertValidEntry( i, j );
#endif
    if( parent.Grid().InGrid() )
        parent.SetRealLocalEntry(i,j,u);
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
template<typename Z>
inline void
DistMatrix<T,STAR,STAR>::SetImagHelper<complex<Z> >::Func
( DistMatrix<complex<Z>,STAR,STAR>& parent, int i, int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::SetImag");
    parent.AssertValidEntry( i, j );
#endif
    if( parent.Grid().InGrid() )
        parent.SetImagLocalEntry(i,j,u);
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
template<typename Z>
inline void
DistMatrix<T,STAR,STAR>::UpdateRealHelper<complex<Z> >::Func
( DistMatrix<complex<Z>,STAR,STAR>& parent, int i, int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::UpdateReal");
    parent.AssertValidEntry( i, j );
#endif
    if( parent.Grid().InGrid() )
        parent.UpdateRealLocalEntry(i,j,u);
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
template<typename Z>
inline void
DistMatrix<T,STAR,STAR>::UpdateImagHelper<complex<Z> >::Func
( DistMatrix<complex<Z>,STAR,STAR>& parent, int i, int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::UpdateImag");
    parent.AssertValidEntry( i, j );
#endif
    if( parent.Grid().InGrid() )
        parent.UpdateImagLocalEntry(i,j,u);
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX

} // namespace elemental
