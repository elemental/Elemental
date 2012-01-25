/*
   Copyright (c) 2009-2012, Jack Poulson
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

namespace elem {

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,STAR,Int>::SetToRandomHermitian()
{ SetToRandomHermitianHelper<T>::Func( *this ); }

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,STAR,Int>::SetToRandomHPD()
{ SetToRandomHPDHelper<T>::Func( *this ); }

template<typename T,typename Int>
inline typename Base<T>::type
DistMatrix<T,STAR,STAR,Int>::GetReal( Int i, Int j ) const
{ return GetRealHelper<T>::Func( *this, i, j ); }

template<typename T,typename Int>
template<typename Z>
inline Z
DistMatrix<T,STAR,STAR,Int>::GetRealHelper<Z>::Func
( const DistMatrix<Z,STAR,STAR,Int>& parent, Int i, Int j )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::GetRealHelper");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T,typename Int>
inline typename Base<T>::type
DistMatrix<T,STAR,STAR,Int>::GetImag( Int i, Int j ) const
{ return GetImagHelper<T>::Func( *this, i, j ); }

template<typename T,typename Int>
template<typename Z>
inline Z
DistMatrix<T,STAR,STAR,Int>::GetImagHelper<Z>::Func
( const DistMatrix<Z,STAR,STAR,Int>& parent, Int i, Int j )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::GetImag");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,STAR,Int>::SetReal
( Int i, Int j, typename Base<T>::type alpha )
{ SetRealHelper<T>::Func( *this, i, j, alpha ); }

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,STAR,STAR,Int>::SetRealHelper<Z>::Func
( DistMatrix<Z,STAR,STAR,Int>& parent, Int i, Int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::SetReal");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,STAR,Int>::SetImag
( Int i, Int j, typename Base<T>::type alpha )
{ SetImagHelper<T>::Func( *this, i, j, alpha ); }

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,STAR,STAR,Int>::SetImagHelper<Z>::Func
( DistMatrix<Z,STAR,STAR,Int>& parent, Int i, Int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::SetImag");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,STAR,Int>::UpdateReal
( Int i, Int j, typename Base<T>::type alpha )
{ UpdateRealHelper<T>::Func( *this, i, j, alpha ); }

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,STAR,STAR,Int>::UpdateRealHelper<Z>::Func
( DistMatrix<Z,STAR,STAR,Int>& parent, Int i, Int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::UpdateReal");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,STAR,Int>::UpdateImag
( Int i, Int j, typename Base<T>::type alpha )
{ UpdateImagHelper<T>::Func( *this, i, j, alpha ); }

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,STAR,STAR,Int>::UpdateImagHelper<Z>::Func
( DistMatrix<Z,STAR,STAR,Int>& parent, Int i, Int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::UpdateImag");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}


template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,STAR,STAR,Int>::SetToRandomHermitianHelper<Z>::Func
( DistMatrix<Z,STAR,STAR,Int>& parent )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::SetToRandomHermitian");
    parent.AssertNotLockedView();
    if( parent.Height() != parent.Width() )
        throw std::logic_error("Hermitian matrices must be square");
#endif
    if( parent.Grid().InGrid() )
        parent.SetToRandom();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,STAR,STAR,Int>::SetToRandomHermitianHelper<Complex<Z> >::
Func( DistMatrix<Complex<Z>,STAR,STAR,Int>& parent )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::SetToRandomHermitian");
    parent.AssertNotLockedView();
    if( parent.Height() != parent.Width() )
        throw std::logic_error("Hermitian matrices must be square");
#endif
    if( parent.Grid().InGrid() )
    {
        const Int height = parent.Height();
        const Int width = parent.Width();

        parent.SetToRandom();

        Complex<Z>* thisLocalBuffer = parent.LocalBuffer();
        const Int thisLDim = parent.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( Int j=0; j<std::min(height,width); ++j )
        {
            const Z value = thisLocalBuffer[j+j*thisLDim].real;
            thisLocalBuffer[j+j*thisLDim] = value;
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,STAR,STAR,Int>::SetToRandomHPDHelper<Z>::Func
( DistMatrix<Z,STAR,STAR,Int>& parent )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::SetToRandomHPD");
    parent.AssertNotLockedView();
    if( parent.Height() != parent.Width() )
        throw std::logic_error("Positive-definite matrices must be square");
#endif
    if( parent.Grid().InGrid() )
    {
        const Int height = parent.Height();
        const Int width = parent.Width();

        parent.SetToRandom();

        Z* thisLocalBuffer = parent.LocalBuffer();
        const Int thisLDim = parent.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( Int j=0; j<std::min(height,width); ++j )
            thisLocalBuffer[j+j*thisLDim] += width;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,STAR,STAR,Int>::SetToRandomHPDHelper<Complex<Z> >::Func
( DistMatrix<Complex<Z>,STAR,STAR,Int>& parent )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::SetToRandomHPD");
    parent.AssertNotLockedView();
    if( parent.Height() != parent.Width() )
        throw std::logic_error("Positive-definite matrices must be square");
#endif
    if( parent.Grid().InGrid() )
    {
        const Int height = parent.Height();
        const Int width = parent.Width();

        parent.SetToRandom();

        Complex<Z>* thisLocalBuffer = parent.LocalBuffer();
        const Int thisLDim = parent.LocalLDim();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for( Int j=0; j<std::min(height,width); ++j )
        {
            const Z value = thisLocalBuffer[j+j*thisLDim].real;
            thisLocalBuffer[j+j*thisLDim] = value + width;
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename Z>
inline Z
DistMatrix<T,STAR,STAR,Int>::GetRealHelper<Complex<Z> >::Func
( const DistMatrix<Complex<Z>,STAR,STAR,Int>& parent, Int i, Int j ) 
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::GetReal");
    parent.AssertValidEntry( i, j );
#endif
    const Int viewingSize = mpi::CommSize( parent.Grid().ViewingComm() );
    const Int owningSize = mpi::GroupSize( parent.Grid().OwningGroup() );
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

template<typename T,typename Int>
template<typename Z>
inline Z
DistMatrix<T,STAR,STAR,Int>::GetImagHelper<Complex<Z> >::Func
( const DistMatrix<Complex<Z>,STAR,STAR,Int>& parent, Int i, Int j ) 
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::GetImag");
    parent.AssertValidEntry( i, j );
#endif
    const Int viewingSize = mpi::CommSize( parent.Grid().ViewingComm() );
    const Int owningSize = mpi::GroupSize( parent.Grid().OwningGroup() );
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

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,STAR,STAR,Int>::SetRealHelper<Complex<Z> >::Func
( DistMatrix<Complex<Z>,STAR,STAR,Int>& parent, Int i, Int j, Z u )
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

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,STAR,STAR,Int>::SetImagHelper<Complex<Z> >::Func
( DistMatrix<Complex<Z>,STAR,STAR,Int>& parent, Int i, Int j, Z u )
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

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,STAR,STAR,Int>::UpdateRealHelper<Complex<Z> >::Func
( DistMatrix<Complex<Z>,STAR,STAR,Int>& parent, Int i, Int j, Z u )
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

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,STAR,STAR,Int>::UpdateImagHelper<Complex<Z> >::Func
( DistMatrix<Complex<Z>,STAR,STAR,Int>& parent, Int i, Int j, Z u )
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

} // namespace elem
