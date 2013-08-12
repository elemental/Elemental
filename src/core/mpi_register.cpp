/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "elemental.hpp"

namespace {
elem::mpi::Op pivotFloatOp;
elem::mpi::Op pivotDoubleOp;
elem::mpi::Op pivotScomplexOp;
elem::mpi::Op pivotDcomplexOp;

elem::mpi::Datatype typeIntInt;
elem::mpi::Datatype typeFloatInt;
elem::mpi::Datatype typeDoubleInt;

elem::mpi::Op maxLocIntOp;
elem::mpi::Op maxLocFloatOp;
elem::mpi::Op maxLocDoubleOp;
}   

namespace elem {
namespace mpi {

// TODO: Decide how to avoid problems with AllReduce being implemented as 
//       ReduceScatter + AllGather, which would potentially artificially 
//       split the datatypes up and confuse this algorithm.
template<typename T>
void
PivotFunc( void* inData, void* outData, int* length, Datatype* datatype )
{           
    if( *length == 0 )
        return;
        
    const Int rowSize = (*length - sizeof(Int))/sizeof(T);

    const T a = ((T*)inData)[0];
    const T b = ((T*)outData)[0];

    const Int inIdx = *((Int*)(((T*)inData)+rowSize));
    const Int outIdx = *((Int*)(((T*)outData)+rowSize));

    if( FastAbs(a) > FastAbs(b) ||
        ( FastAbs(a) == FastAbs(b) && inIdx < outIdx ) )
    {
        // Copy the row of T's
        for( Int j=0; j<rowSize; ++j )
            ((T*)outData)[j] = ((T*)inData)[j];

        // Copy the integer
        *((Int*)(((T*)outData)+rowSize)) = inIdx;
    }
}

template<typename T>
void
MaxLocFunc
( ValueInt<T>* inData, ValueInt<T>* outData, int* length, Datatype* datatype )
{           
    for( int j=0; j<*length; ++j )
        if( inData[j].value > outData[j].value )
            outData[j] = inData[j];
}

template<>
Datatype& ValueIntType<Int>()
{ return ::typeIntInt; }

template<>
Datatype& ValueIntType<float>()
{ return ::typeFloatInt; }

template<>
Datatype& ValueIntType<double>()
{ return ::typeDoubleInt; }

template<typename T>
void CreateValueIntType()
{
#ifndef RELEASE
    CallStackEntry cse("CreateValueIntType<Int>");
#endif
    Datatype& type = ValueIntType<T>();
#ifndef RELEASE
    if( type != 0 )
        LogicError("type was already created");
#endif
    Datatype typeList[2];
    typeList[0] = TypeMap<T>();
    typeList[1] = TypeMap<Int>();
    
    int blockLengths[2];
    blockLengths[0] = 1;
    blockLengths[1] = 1; 

    ValueInt<T> v;
    MPI_Aint startAddr, valueAddr, indexAddr;
    MPI_Address( &v,       &startAddr );
    MPI_Address( &v.value, &valueAddr );
    MPI_Address( &v.index, &indexAddr );

    MPI_Aint displs[2];
    displs[0] = valueAddr - startAddr;
    displs[1] = indexAddr - startAddr;

    MPI_Type_create_struct( 2, blockLengths, displs, typeList, &type );
    MPI_Type_commit( &type );
}
template void CreateValueIntType<Int>();
template void CreateValueIntType<float>();
template void CreateValueIntType<double>();

template<>
void CreatePivotOp<float>()
{
#ifndef RELEASE
    CallStackEntry cse("CreatePivotOp<float>");
#endif
    OpCreate( (UserFunction*)PivotFunc<float>, true, ::pivotFloatOp );
}

template<>
void CreatePivotOp<double>()
{
#ifndef RELEASE
    CallStackEntry cse("CreatePivotOp<double>");
#endif  
    OpCreate( (UserFunction*)PivotFunc<double>, true, ::pivotDoubleOp );
}

template<>
void CreatePivotOp<Complex<float> >()
{
#ifndef RELEASE
    CallStackEntry cse("CreatePivotOp");
#endif
    OpCreate
    ( (UserFunction*)PivotFunc<Complex<float> >, true, ::pivotScomplexOp );
}

template<>
void CreatePivotOp<Complex<double> >()
{
#ifndef RELEASE
    CallStackEntry cse("CreatePivotOp");
#endif
    OpCreate
    ( (UserFunction*)PivotFunc<Complex<double> >, true, ::pivotDcomplexOp );
}

template<>
void CreateMaxLocOp<Int>()
{
#ifndef RELEASE
    CallStackEntry cse("CreateMaxLocOp<Int>");
#endif
    OpCreate( (UserFunction*)MaxLocFunc<Int>, true, ::maxLocIntOp );
}

template<>
void CreateMaxLocOp<float>()
{
#ifndef RELEASE
    CallStackEntry cse("CreateMaxLocOp<float>");
#endif
    OpCreate( (UserFunction*)MaxLocFunc<float>, true, ::maxLocFloatOp );
}

template<>
void CreateMaxLocOp<double>()
{
#ifndef RELEASE
    CallStackEntry cse("CreateMaxLocOp<double>");
#endif
    OpCreate( (UserFunction*)MaxLocFunc<double>, true, ::maxLocDoubleOp );
}

template<>
void DestroyPivotOp<float>()
{
#ifndef RELEASE
    CallStackEntry cse("DestroyPivotOp<float>");
#endif
    OpFree( ::pivotFloatOp );
}

template<>
void DestroyPivotOp<double>()
{
#ifndef RELEASE
    CallStackEntry cse("DestroyPivotOp<double>");
#endif
    OpFree( ::pivotDoubleOp );
}

template<>
void DestroyPivotOp<Complex<float> >()
{
#ifndef RELEASE
    CallStackEntry cse("DestroyPivotOp");
#endif
    OpFree( ::pivotScomplexOp );
}

template<>
void DestroyPivotOp<Complex<double> >()
{
#ifndef RELEASE
    CallStackEntry cse("DestroyPivotOp");
#endif
    OpFree( ::pivotDcomplexOp );
}

template<>
void DestroyMaxLocOp<Int>()
{
#ifndef RELEASE
    CallStackEntry cse("DestroyMaxLocOp<Int>");
#endif
    OpFree( ::maxLocIntOp );
}

template<>
void DestroyMaxLocOp<float>()
{
#ifndef RELEASE
    CallStackEntry cse("DestroyMaxLocOp<float>");
#endif
    OpFree( ::maxLocFloatOp );
}

template<>
void DestroyMaxLocOp<double>()
{
#ifndef RELEASE
    CallStackEntry cse("DestroyMaxLocOp<double>");
#endif
    OpFree( ::maxLocDoubleOp );
}

template<>
Op PivotOp<float>()
{
#ifndef RELEASE
    CallStackEntry cse("PivotOp<float>");
#endif
    return ::pivotFloatOp;
}

template<>
Op PivotOp<double>()
{
#ifndef RELEASE
    CallStackEntry cse("PivotOp<double>");
#endif
    return ::pivotDoubleOp;
}

template<>
Op PivotOp<Complex<float> >()
{
#ifndef RELEASE
    CallStackEntry cse("PivotOp");
#endif
    return ::pivotScomplexOp;
}

template<>
Op PivotOp<Complex<double> >()
{
#ifndef RELEASE
    CallStackEntry cse("PivotOp");
#endif
    return ::pivotDcomplexOp;
}

template<>
Op MaxLocOp<Int>()
{
#ifndef RELEASE
    CallStackEntry cse("MaxLocOp<Int>");
#endif
    return ::maxLocIntOp;
}

template<>
Op MaxLocOp<float>()
{
#ifndef RELEASE
    CallStackEntry cse("MaxLocOp<float>");
#endif
    return ::maxLocFloatOp;
}

template<>
Op MaxLocOp<double>()
{
#ifndef RELEASE
    CallStackEntry cse("MaxLocOp<double>");
#endif
    return ::maxLocDoubleOp;
}

template void
PivotFunc<float>
( void* inData, void* outData, int* length, Datatype* datatype );
template void
PivotFunc<double>
( void* inData, void* outData, int* length, Datatype* datatype );
template void
PivotFunc<Complex<float> >
( void* inData, void* outData, int* length, Datatype* datatype );
template void
PivotFunc<Complex<double> >
( void* inData, void* outData, int* length, Datatype* datatype );

template void
MaxLocFunc<Int>
( ValueInt<Int>* inData, ValueInt<Int>* outData, int* length, 
  Datatype* datatype );
template void
MaxLocFunc<float>
( ValueInt<float>* inData, ValueInt<float>* outData, int* length, 
  Datatype* datatype );
template void
MaxLocFunc<double>
( ValueInt<double>* inData, ValueInt<double>* outData, int* length, 
  Datatype* datatype );

} // namespace mpi
} // namespace elem
