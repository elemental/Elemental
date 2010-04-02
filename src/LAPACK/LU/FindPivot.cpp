/*
   Copyright 2009-2010 Jack Poulson

   This file is part of Elemental.

   Elemental is free software: you can redistribute it and/or modify it under
   the terms of the GNU Lesser General Public License as published by the
   Free Software Foundation; either version 3 of the License, or 
   (at your option) any later version.

   Elemental is distributed in the hope that it will be useful, but 
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with Elemental. If not, see <http://www.gnu.org/licenses/>.
*/
#include "Elemental/LAPACKInternal.hpp"
using namespace std;
using namespace Elemental;
using namespace Elemental::BLAS;
using namespace Elemental::LAPACK::Internal;
using namespace Elemental::wrappers::MPI;

#ifndef RELEASE
static bool createdPivotOpFloat = false;
static bool createdPivotOpDouble = false;
#ifndef WITHOUT_COMPLEX
static bool createdPivotOpScomplex = false;
static bool createdPivotOpDcomplex = false;
#endif
#endif
static MPI_Op pivotOpFloat;
static MPI_Op pivotOpDouble;
#ifndef WITHOUT_COMPLEX
static MPI_Op pivotOpScomplex;
static MPI_Op pivotOpDcomplex;
#endif

template<typename T>
void
Elemental::LAPACK::Internal::PivotFunc
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
Elemental::LAPACK::Internal::CreatePivotOp<float>()
{
#ifndef RELEASE
    PushCallStack("LAPACK::Internal::CreatePivotOp<float>");
    if( ::createdPivotOpFloat )
        throw "Already created pivot op.";
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
Elemental::LAPACK::Internal::CreatePivotOp<double>()
{
#ifndef RELEASE
    PushCallStack("LAPACK::Internal::CreatePivotOp<double>");
    if( ::createdPivotOpDouble )
        throw "Already created pivot op.";
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
Elemental::LAPACK::Internal::CreatePivotOp<scomplex>()
{
#ifndef RELEASE
    PushCallStack("LAPACK::Internal::CreatePivotOp<scomplex>");
    if( ::createdPivotOpScomplex )
        throw "Alread created pivot op.";
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
Elemental::LAPACK::Internal::CreatePivotOp<dcomplex>()
{
#ifndef RELEASE
    PushCallStack("LAPACK::Internal::CreatePivotOp<dcomplex>");
    if( ::createdPivotOpDcomplex )
        throw "Already created pivot op.";
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
Elemental::LAPACK::Internal::DestroyPivotOp<float>()
{
#ifndef RELEASE
    PushCallStack("LAPACK::Internal::DestroyPivotOp<float>");
    if( ! ::createdPivotOpFloat )
        throw "Have not created this pivot op.";
#endif
    MPI_Op_free( &::pivotOpFloat );
#ifndef RELEASE
    ::createdPivotOpFloat = false;
    PopCallStack();
#endif
}

template<>
void
Elemental::LAPACK::Internal::DestroyPivotOp<double>()
{
#ifndef RELEASE
    PushCallStack("LAPACK::Internal::DestroyPivotOp<double>");
    if( ! ::createdPivotOpDouble )
        throw "Have not created ths pivot op.";
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
Elemental::LAPACK::Internal::DestroyPivotOp<scomplex>()
{
#ifndef RELEASE
    PushCallStack("LAPACK::Internal::DestroyPivotOp<scomplex>");
    if( ! ::createdPivotOpScomplex )
        throw "Have not created this pivot op.";
#endif
    MPI_Op_free( &::pivotOpScomplex );
#ifndef RELEASE
    ::createdPivotOpScomplex = false;
    PopCallStack();
#endif
}

template<>
void
Elemental::LAPACK::Internal::DestroyPivotOp<dcomplex>()
{
#ifndef RELEASE
    PushCallStack("LAPACK::Internal::DestroyPivotOp<dcomplex>");
    if( ! ::createdPivotOpDcomplex )
        throw "Have not created this pivot op.";
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
Elemental::LAPACK::Internal::PivotOp<float>()
{ 
#ifndef RELEASE
    PushCallStack("LAPACK::Internal::PivotOp<float>");
    if( ! ::createdPivotOpFloat )
        throw "Tried to return uncreated pivot op.";
    PopCallStack();
#endif
    return ::pivotOpFloat; 
}

template<>
MPI_Op
Elemental::LAPACK::Internal::PivotOp<double>()
{ 
#ifndef RELEASE
    PushCallStack("LAPACK::Internal::PivotOp<double>");
    if( ! ::createdPivotOpDouble )
        throw "Tried to return uncreated pivot op.";
    PopCallStack();
#endif
    return ::pivotOpDouble; 
}

#ifndef WITHOUT_COMPLEX
template<>
MPI_Op
Elemental::LAPACK::Internal::PivotOp<scomplex>()
{ 
#ifndef RELEASE
    PushCallStack("LAPACK::Internal::PivotOp<scomplex>");
    if( ! ::createdPivotOpScomplex )
        throw "Tried to return uncreated pivot op.";
    PopCallStack();
#endif
    return ::pivotOpScomplex; 
}

template<>
MPI_Op
Elemental::LAPACK::Internal::PivotOp<dcomplex>()
{ 
#ifndef RELEASE
    PushCallStack("LAPACK::Internal::PivotOp<dcomplex>");
    if( ! ::createdPivotOpDcomplex )
        throw "Tried to return uncreated pivot op.";
    PopCallStack();
#endif
    return ::pivotOpDcomplex; 
}
#endif

template void 
Elemental::LAPACK::Internal::PivotFunc<float>
( void* inData, void* outData, int* length, MPI_Datatype* datatype );

template void
Elemental::LAPACK::Internal::PivotFunc<double>
( void* inData, void* outData, int* length, MPI_Datatype* datatype );

#ifndef WITHOUT_COMPLEX
template void
Elemental::LAPACK::Internal::PivotFunc<scomplex>
( void* inData, void* outData, int* length, MPI_Datatype* datatype );

template void
Elemental::LAPACK::Internal::PivotFunc<dcomplex>
( void* inData, void* outData, int* length, MPI_Datatype* datatype );
#endif

