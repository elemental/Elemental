/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
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
PivotFunc( void* inData, void* outData, int* length, mpi::Datatype* datatype )
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
void CreatePivotOp<float>()
{
#ifndef RELEASE
    PushCallStack("CreatePivotOp<float>");
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
void CreatePivotOp<double>()
{
#ifndef RELEASE
    PushCallStack("CreatePivotOp<double>");
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
void CreatePivotOp<Complex<float> >()
{
#ifndef RELEASE
    PushCallStack("CreatePivotOp");
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
void CreatePivotOp<Complex<double> >()
{
#ifndef RELEASE
    PushCallStack("CreatePivotOp");
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
void DestroyPivotOp<float>()
{
#ifndef RELEASE
    PushCallStack("DestroyPivotOp<float>");
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
void DestroyPivotOp<double>()
{
#ifndef RELEASE
    PushCallStack("DestroyPivotOp<double>");
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
void DestroyPivotOp<Complex<float> >()
{
#ifndef RELEASE
    PushCallStack("DestroyPivotOp");
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
void DestroyPivotOp<Complex<double> >()
{
#ifndef RELEASE
    PushCallStack("DestroyPivotOp");
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
mpi::Op PivotOp<float>()
{
#ifndef RELEASE
    PushCallStack("PivotOp<float>");
    if( ! ::createdPivotOpFloat )
        throw std::logic_error("Tried to return uncreated pivot op");
    PopCallStack();
#endif
    return ::pivotOpFloat;
}

template<>
mpi::Op PivotOp<double>()
{
#ifndef RELEASE
    PushCallStack("PivotOp<double>");
    if( ! ::createdPivotOpDouble )
        throw std::logic_error("Tried to return uncreated pivot op");
    PopCallStack();
#endif
    return ::pivotOpDouble;
}

template<>
mpi::Op PivotOp<Complex<float> >()
{
#ifndef RELEASE
    PushCallStack("PivotOp");
    if( ! ::createdPivotOpScomplex )
        throw std::logic_error("Tried to return uncreated pivot op");
    PopCallStack();
#endif
    return ::pivotOpScomplex;
}

template<>
mpi::Op PivotOp<Complex<double> >()
{
#ifndef RELEASE
    PushCallStack("PivotOp");
    if( ! ::createdPivotOpDcomplex )
        throw std::logic_error("Tried to return uncreated pivot op");
    PopCallStack();
#endif
    return ::pivotOpDcomplex;
}

template void
PivotFunc<float>
( void* inData, void* outData, int* length, mpi::Datatype* datatype );

template void
PivotFunc<double>
( void* inData, void* outData, int* length, mpi::Datatype* datatype );

template void
PivotFunc<elem::Complex<float> >
( void* inData, void* outData, int* length, mpi::Datatype* datatype );

template void
PivotFunc<elem::Complex<double> >
( void* inData, void* outData, int* length, mpi::Datatype* datatype );

} // namespace elem
