/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "elemental.hpp"

namespace {
elem::mpi::Datatype typeIntInt;
elem::mpi::Datatype typeFloatInt;
elem::mpi::Datatype typeDoubleInt;

elem::mpi::Op maxLocIntOp;
elem::mpi::Op maxLocFloatOp;
elem::mpi::Op maxLocDoubleOp;
}   

namespace elem {
namespace mpi {

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
    CallStackEntry cse("CreateValueIntType");
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

    Datatype& type = ValueIntType<T>();
    MPI_Type_create_struct( 2, blockLengths, displs, typeList, &type );
    MPI_Type_commit( &type );
}
template void CreateValueIntType<Int>();
template void CreateValueIntType<float>();
template void CreateValueIntType<double>();

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
