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
#include "elemental.hpp"

namespace {
bool createdPivotOpFloat = false;
bool createdPivotOpDouble = false;
#ifndef WITHOUT_COMPLEX
bool createdPivotOpScomplex = false;
bool createdPivotOpDcomplex = false;
#endif
elemental::mpi::Op pivotOpFloat;
elemental::mpi::Op pivotOpDouble;
#ifndef WITHOUT_COMPLEX
elemental::mpi::Op pivotOpScomplex;
elemental::mpi::Op pivotOpDcomplex;
#endif
}   
template<typename T> // represents a real or complex ring
void
elemental::advanced::internal::PivotFunc
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
void
elemental::advanced::internal::CreatePivotOp<float>()
{
#ifndef RELEASE
    PushCallStack("advanced::internal::CreatePivotOp<float>");
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
void
elemental::advanced::internal::CreatePivotOp<double>()
{
#ifndef RELEASE
    PushCallStack("advanced::internal::CreatePivotOp<double>");
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

#ifndef WITHOUT_COMPLEX
template<>
void
elemental::advanced::internal::CreatePivotOp<elemental::scomplex>()
{
#ifndef RELEASE
    PushCallStack("advanced::internal::CreatePivotOp<scomplex>");
    if( ::createdPivotOpScomplex )
        throw std::logic_error("Alread created pivot op");
#endif
    mpi::OpCreate
    ( (mpi::UserFunction*)PivotFunc<scomplex>, true, ::pivotOpScomplex );
    ::createdPivotOpScomplex = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<>
void
elemental::advanced::internal::CreatePivotOp<elemental::dcomplex>()
{
#ifndef RELEASE
    PushCallStack("advanced::internal::CreatePivotOp<dcomplex>");
    if( ::createdPivotOpDcomplex )
        throw std::logic_error("Already created pivot op");
#endif
    mpi::OpCreate
    ( (mpi::UserFunction*)PivotFunc<dcomplex>, true, ::pivotOpDcomplex );
    ::createdPivotOpDcomplex = true;
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif

template<>
void
elemental::advanced::internal::DestroyPivotOp<float>()
{
#ifndef RELEASE
    PushCallStack("advanced::internal::DestroyPivotOp<float>");
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
void
elemental::advanced::internal::DestroyPivotOp<double>()
{
#ifndef RELEASE
    PushCallStack("advanced::internal::DestroyPivotOp<double>");
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

#ifndef WITHOUT_COMPLEX
template<>
void
elemental::advanced::internal::DestroyPivotOp<elemental::scomplex>()
{
#ifndef RELEASE
    PushCallStack("advanced::internal::DestroyPivotOp<scomplex>");
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
void
elemental::advanced::internal::DestroyPivotOp<elemental::dcomplex>()
{
#ifndef RELEASE
    PushCallStack("advanced::internal::DestroyPivotOp<dcomplex>");
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
#endif

template<>
MPI_Op
elemental::advanced::internal::PivotOp<float>()
{
#ifndef RELEASE
    PushCallStack("advanced::internal::PivotOp<float>");
    if( ! ::createdPivotOpFloat )
        throw std::logic_error("Tried to return uncreated pivot op");
    PopCallStack();
#endif
    return ::pivotOpFloat;
}

template<>
elemental::mpi::Op
elemental::advanced::internal::PivotOp<double>()
{
#ifndef RELEASE
    PushCallStack("advanced::internal::PivotOp<double>");
    if( ! ::createdPivotOpDouble )
        throw std::logic_error("Tried to return uncreated pivot op");
    PopCallStack();
#endif
    return ::pivotOpDouble;
}

#ifndef WITHOUT_COMPLEX
template<>
elemental::mpi::Op
elemental::advanced::internal::PivotOp<elemental::scomplex>()
{
#ifndef RELEASE
    PushCallStack("advanced::internal::PivotOp<scomplex>");
    if( ! ::createdPivotOpScomplex )
        throw std::logic_error("Tried to return uncreated pivot op");
    PopCallStack();
#endif
    return ::pivotOpScomplex;
}

template<>
elemental::mpi::Op
elemental::advanced::internal::PivotOp<elemental::dcomplex>()
{
#ifndef RELEASE
    PushCallStack("advanced::internal::PivotOp<dcomplex>");
    if( ! ::createdPivotOpDcomplex )
        throw std::logic_error("Tried to return uncreated pivot op");
    PopCallStack();
#endif
    return ::pivotOpDcomplex;
}
#endif // WITHOUT_COMPLEX

template void
elemental::advanced::internal::PivotFunc<float>
( void* inData, void* outData, int* length, mpi::Datatype* datatype );

template void
elemental::advanced::internal::PivotFunc<double>
( void* inData, void* outData, int* length, mpi::Datatype* datatype );

#ifndef WITHOUT_COMPLEX
template void
elemental::advanced::internal::PivotFunc<elemental::scomplex>
( void* inData, void* outData, int* length, mpi::Datatype* datatype );

template void
elemental::advanced::internal::PivotFunc<elemental::dcomplex>
( void* inData, void* outData, int* length, mpi::Datatype* datatype );
#endif
