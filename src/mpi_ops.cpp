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
#include "elemental.hpp"

namespace {
bool createdPivotOpFloat = false;
bool createdPivotOpDouble = false;
bool createdPivotOpScomplex = false;
bool createdPivotOpDcomplex = false;
elem::mpi::Op pivotOpFloat;
elem::mpi::Op pivotOpDouble;
elem::mpi::Op pivotOpScomplex;
elem::mpi::Op pivotOpDcomplex;
}   

namespace elem {

template<typename T>
void
internal::PivotFunc
( void* inData, void* outData, int* length, mpi::Datatype* datatype )
{           
    if( *length == 0 )
        return;
        
    const int rowSize = (*length - sizeof(int))/sizeof(T);

    const T a = ((T*)inData)[0];
    const T b = ((T*)outData)[0];

    const int inIdx = *((int*)(((T*)inData)+rowSize));
    const int outIdx = *((int*)(((T*)outData)+rowSize));

    if( FastAbs(a) > FastAbs(b) ||
        ( FastAbs(a) == FastAbs(b) && inIdx < outIdx ) )
    {
        // Copy the row of T's
        for( int j=0; j<rowSize; ++j )
            ((T*)outData)[j] = ((T*)inData)[j];

        // Copy the integer
        *((int*)(((T*)outData)+rowSize)) = inIdx;
    }
}

template<>
void internal::CreatePivotOp<float>()
{
#ifndef RELEASE
    PushCallStack("internal::CreatePivotOp<float>");
    if( ::createdPivotOpFloat )
        throw std::logic_error("Already created pivot op");
#endif
    mpi::OpCreate
    ( (mpi::UserFunction*)PivotFunc<float>, true, ::pivotOpFloat );
    ::createdPivotOpFloat = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<>
void internal::CreatePivotOp<double>()
{
#ifndef RELEASE
    PushCallStack("internal::CreatePivotOp<double>");
    if( ::createdPivotOpDouble )
        throw std::logic_error("Already created pivot op");
#endif  
    mpi::OpCreate
    ( (mpi::UserFunction*)PivotFunc<double>, true, ::pivotOpDouble );
    ::createdPivotOpDouble = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<>
void internal::CreatePivotOp<Complex<float> >()
{
#ifndef RELEASE
    PushCallStack("internal::CreatePivotOp");
    if( ::createdPivotOpScomplex )
        throw std::logic_error("Alread created pivot op");
#endif
    mpi::OpCreate
    ( (mpi::UserFunction*)PivotFunc<Complex<float> >, true, 
      ::pivotOpScomplex );
    ::createdPivotOpScomplex = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<>
void internal::CreatePivotOp<Complex<double> >()
{
#ifndef RELEASE
    PushCallStack("internal::CreatePivotOp");
    if( ::createdPivotOpDcomplex )
        throw std::logic_error("Already created pivot op");
#endif
    mpi::OpCreate
    ( (mpi::UserFunction*)PivotFunc<Complex<double> >, true, 
      ::pivotOpDcomplex );
    ::createdPivotOpDcomplex = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<>
void internal::DestroyPivotOp<float>()
{
#ifndef RELEASE
    PushCallStack("internal::DestroyPivotOp<float>");
    if( ! ::createdPivotOpFloat )
        throw std::logic_error("Have not created this pivot op");
#endif
    if( ::createdPivotOpFloat )
        mpi::OpFree( ::pivotOpFloat );
    ::createdPivotOpFloat = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<>
void internal::DestroyPivotOp<double>()
{
#ifndef RELEASE
    PushCallStack("internal::DestroyPivotOp<double>");
    if( ! ::createdPivotOpDouble )
        throw std::logic_error("Have not created ths pivot op");
#endif
    if( ::createdPivotOpDouble )
        mpi::OpFree( ::pivotOpDouble );
    ::createdPivotOpDouble = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<>
void internal::DestroyPivotOp<Complex<float> >()
{
#ifndef RELEASE
    PushCallStack("internal::DestroyPivotOp");
    if( ! ::createdPivotOpScomplex )
        throw std::logic_error("Have not created this pivot op");
#endif
    if( ::createdPivotOpScomplex )
        mpi::OpFree( ::pivotOpScomplex );
    ::createdPivotOpScomplex = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<>
void internal::DestroyPivotOp<Complex<double> >()
{
#ifndef RELEASE
    PushCallStack("internal::DestroyPivotOp");
    if( ! ::createdPivotOpDcomplex )
        throw std::logic_error("Have not created this pivot op");
#endif
    if( ::createdPivotOpDcomplex )
        mpi::OpFree( ::pivotOpDcomplex );
    ::createdPivotOpDcomplex = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<>
mpi::Op internal::PivotOp<float>()
{
#ifndef RELEASE
    PushCallStack("internal::PivotOp<float>");
    if( ! ::createdPivotOpFloat )
        throw std::logic_error("Tried to return uncreated pivot op");
    PopCallStack();
#endif
    return ::pivotOpFloat;
}

template<>
mpi::Op internal::PivotOp<double>()
{
#ifndef RELEASE
    PushCallStack("internal::PivotOp<double>");
    if( ! ::createdPivotOpDouble )
        throw std::logic_error("Tried to return uncreated pivot op");
    PopCallStack();
#endif
    return ::pivotOpDouble;
}

template<>
mpi::Op internal::PivotOp<Complex<float> >()
{
#ifndef RELEASE
    PushCallStack("internal::PivotOp");
    if( ! ::createdPivotOpScomplex )
        throw std::logic_error("Tried to return uncreated pivot op");
    PopCallStack();
#endif
    return ::pivotOpScomplex;
}

template<>
mpi::Op internal::PivotOp<Complex<double> >()
{
#ifndef RELEASE
    PushCallStack("internal::PivotOp");
    if( ! ::createdPivotOpDcomplex )
        throw std::logic_error("Tried to return uncreated pivot op");
    PopCallStack();
#endif
    return ::pivotOpDcomplex;
}

} // namespace elem

template void
elem::internal::PivotFunc<float>
( void* inData, void* outData, int* length, mpi::Datatype* datatype );

template void
elem::internal::PivotFunc<double>
( void* inData, void* outData, int* length, mpi::Datatype* datatype );

template void
elem::internal::PivotFunc<elem::Complex<float> >
( void* inData, void* outData, int* length, mpi::Datatype* datatype );

template void
elem::internal::PivotFunc<elem::Complex<double> >
( void* inData, void* outData, int* length, mpi::Datatype* datatype );
