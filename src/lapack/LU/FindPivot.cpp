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
#include "elemental/lapack_internal.hpp"
using namespace std;
using namespace elemental;
using namespace elemental::blas;
using namespace elemental::lapack::internal;
using namespace elemental::import::mpi;

namespace {
#ifndef RELEASE
bool createdPivotOpFloat = false;
bool createdPivotOpDouble = false;
#ifndef WITHOUT_COMPLEX
bool createdPivotOpScomplex = false;
bool createdPivotOpDcomplex = false;
#endif
#endif
MPI_Op pivotOpFloat;
MPI_Op pivotOpDouble;
#ifndef WITHOUT_COMPLEX
MPI_Op pivotOpScomplex;
MPI_Op pivotOpDcomplex;
#endif
}
template<typename T> // represents a real or complex ring
void
elemental::lapack::internal::PivotFunc
( void* inData, void* outData, int* length, MPI_Datatype* datatype )
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
elemental::lapack::internal::CreatePivotOp<float>()
{
#ifndef RELEASE
    PushCallStack("lapack::internal::CreatePivotOp<float>");
    if( ::createdPivotOpFloat )
        throw logic_error( "Already created pivot op." );
#endif
    MPI_Op_create
    ( (MPI_User_function*)PivotFunc<float>, 1, &::pivotOpFloat );
#ifndef RELEASE
    ::createdPivotOpFloat = true;
    PopCallStack();
#endif
}

template<>
void
elemental::lapack::internal::CreatePivotOp<double>()
{
#ifndef RELEASE
    PushCallStack("lapack::internal::CreatePivotOp<double>");
    if( ::createdPivotOpDouble )
        throw logic_error( "Already created pivot op." );
#endif
    MPI_Op_create
    ( (MPI_User_function*)PivotFunc<double>, 1, &::pivotOpDouble );
#ifndef RELEASE
    ::createdPivotOpDouble = true;
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
template<>
void
elemental::lapack::internal::CreatePivotOp<scomplex>()
{
#ifndef RELEASE
    PushCallStack("lapack::internal::CreatePivotOp<scomplex>");
    if( ::createdPivotOpScomplex )
        throw logic_error( "Alread created pivot op." );
#endif
    MPI_Op_create
    ( (MPI_User_function*)PivotFunc<scomplex>, 1, &::pivotOpScomplex );
#ifndef RELEASE
    ::createdPivotOpScomplex = true;
    PopCallStack();
#endif
}

template<>
void
elemental::lapack::internal::CreatePivotOp<dcomplex>()
{
#ifndef RELEASE
    PushCallStack("lapack::internal::CreatePivotOp<dcomplex>");
    if( ::createdPivotOpDcomplex )
        throw logic_error( "Already created pivot op." );
#endif
    MPI_Op_create
    ( (MPI_User_function*)PivotFunc<dcomplex>, 1, &::pivotOpDcomplex );
#ifndef RELEASE
    ::createdPivotOpDcomplex = true;
    PopCallStack();
#endif
}
#endif

template<>
void
elemental::lapack::internal::DestroyPivotOp<float>()
{
#ifndef RELEASE
    PushCallStack("lapack::internal::DestroyPivotOp<float>");
    if( ! ::createdPivotOpFloat )
        throw logic_error( "Have not created this pivot op." );
#endif
    MPI_Op_free( &::pivotOpFloat );
#ifndef RELEASE
    ::createdPivotOpFloat = false;
    PopCallStack();
#endif
}

template<>
void
elemental::lapack::internal::DestroyPivotOp<double>()
{
#ifndef RELEASE
    PushCallStack("lapack::internal::DestroyPivotOp<double>");
    if( ! ::createdPivotOpDouble )
        throw logic_error( "Have not created ths pivot op." );
#endif
    MPI_Op_free( &::pivotOpDouble );
#ifndef RELEASE
    ::createdPivotOpDouble = false;
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
template<>
void
elemental::lapack::internal::DestroyPivotOp<scomplex>()
{
#ifndef RELEASE
    PushCallStack("lapack::internal::DestroyPivotOp<scomplex>");
    if( ! ::createdPivotOpScomplex )
        throw logic_error( "Have not created this pivot op." );
#endif
    MPI_Op_free( &::pivotOpScomplex );
#ifndef RELEASE
    ::createdPivotOpScomplex = false;
    PopCallStack();
#endif
}

template<>
void
elemental::lapack::internal::DestroyPivotOp<dcomplex>()
{
#ifndef RELEASE
    PushCallStack("lapack::internal::DestroyPivotOp<dcomplex>");
    if( ! ::createdPivotOpDcomplex )
        throw logic_error( "Have not created this pivot op." );
#endif
    MPI_Op_free( &::pivotOpDcomplex );
#ifndef RELEASE
    ::createdPivotOpDcomplex = false;
    PopCallStack();
#endif
}
#endif

template<>
MPI_Op
elemental::lapack::internal::PivotOp<float>()
{ 
#ifndef RELEASE
    PushCallStack("lapack::internal::PivotOp<float>");
    if( ! ::createdPivotOpFloat )
        throw logic_error( "Tried to return uncreated pivot op." );
    PopCallStack();
#endif
    return ::pivotOpFloat; 
}

template<>
MPI_Op
elemental::lapack::internal::PivotOp<double>()
{ 
#ifndef RELEASE
    PushCallStack("lapack::internal::PivotOp<double>");
    if( ! ::createdPivotOpDouble )
        throw logic_error( "Tried to return uncreated pivot op." );
    PopCallStack();
#endif
    return ::pivotOpDouble; 
}

#ifndef WITHOUT_COMPLEX
template<>
MPI_Op
elemental::lapack::internal::PivotOp<scomplex>()
{ 
#ifndef RELEASE
    PushCallStack("lapack::internal::PivotOp<scomplex>");
    if( ! ::createdPivotOpScomplex )
        throw logic_error( "Tried to return uncreated pivot op." );
    PopCallStack();
#endif
    return ::pivotOpScomplex; 
}

template<>
MPI_Op
elemental::lapack::internal::PivotOp<dcomplex>()
{ 
#ifndef RELEASE
    PushCallStack("lapack::internal::PivotOp<dcomplex>");
    if( ! ::createdPivotOpDcomplex )
        throw logic_error( "Tried to return uncreated pivot op." );
    PopCallStack();
#endif
    return ::pivotOpDcomplex; 
}
#endif

template void 
elemental::lapack::internal::PivotFunc<float>
( void* inData, void* outData, int* length, MPI_Datatype* datatype );

template void
elemental::lapack::internal::PivotFunc<double>
( void* inData, void* outData, int* length, MPI_Datatype* datatype );

#ifndef WITHOUT_COMPLEX
template void
elemental::lapack::internal::PivotFunc<scomplex>
( void* inData, void* outData, int* length, MPI_Datatype* datatype );

template void
elemental::lapack::internal::PivotFunc<dcomplex>
( void* inData, void* outData, int* length, MPI_Datatype* datatype );
#endif

