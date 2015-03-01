/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace {

// Datatypes
// =========
El::mpi::Datatype QuadType, QuadComplexType;

El::mpi::Datatype IntIntType, floatIntType, doubleIntType;
#ifdef EL_HAVE_QUAD
El::mpi::Datatype QuadIntType;
#endif

El::mpi::Datatype IntIntPairType, floatIntPairType, doubleIntPairType;
#ifdef EL_HAVE_QUAD
El::mpi::Datatype QuadIntPairType;
#endif

// Operations
// ==========
#ifdef EL_HAVE_QUAD
El::mpi::Op minQuadOp;
El::mpi::Op maxQuadOp;
El::mpi::Op sumQuadOp, sumQuadComplexOp;
#endif

El::mpi::Op maxLocIntOp, maxLocFloatOp, maxLocDoubleOp;
#ifdef EL_HAVE_QUAD
El::mpi::Op maxLocQuadOp;
#endif

El::mpi::Op maxLocPairIntOp, maxLocPairFloatOp, maxLocPairDoubleOp;
#ifdef EL_HAVE_QUAD
El::mpi::Op maxLocPairQuadOp;
#endif

El::mpi::Op minLocIntOp, minLocFloatOp, minLocDoubleOp;
#ifdef EL_HAVE_QUAD
El::mpi::Op minLocQuadOp;
#endif

El::mpi::Op minLocPairIntOp, minLocPairFloatOp, minLocPairDoubleOp;
#ifdef EL_HAVE_QUAD
El::mpi::Op minLocPairQuadOp;
#endif

} // anonymouse namespace   

namespace El {
namespace mpi {

#ifdef EL_HAVE_QUAD
static void
MaxQuad( void* inVoid, void* outVoid, int* length, Datatype* datatype )
{
    const Quad* inData = static_cast<Quad*>(inVoid);
    Quad* outData = static_cast<Quad*>(outVoid);
    for( int j=0; j<*length; ++j )
    {
        if( inData[j] > outData[j] )
            outData[j] = inData[j];
    }
}

static void
MinQuad( void* inVoid, void* outVoid, int* length, Datatype* datatype )
{
    const Quad* inData = static_cast<Quad*>(inVoid);
    Quad* outData = static_cast<Quad*>(outVoid);
    for( int j=0; j<*length; ++j )
    {
        if( inData[j] < outData[j] )
            outData[j] = inData[j];
    }
}

static void
SumQuad( void* inVoid, void* outVoid, int* length, Datatype* datatype )
{
    const Quad* inData = static_cast<Quad*>(inVoid);
    Quad* outData = static_cast<Quad*>(outVoid);
    for( int j=0; j<*length; ++j )
        outData[j] += inData[j];
}

static void
SumQuadComplex( void* inVoid, void* outVoid, int* length, Datatype* datatype )
{
    const Complex<Quad>* inData = static_cast<Complex<Quad>*>(inVoid);
    Complex<Quad>* outData = static_cast<Complex<Quad>*>(outVoid);
    for( int j=0; j<*length; ++j )
        outData[j] += inData[j];
}
#endif 

template<typename T>
static void
MaxLocFunc( void* inVoid, void* outVoid, int* length, Datatype* datatype )
{           
    const ValueInt<T>* inData = static_cast<ValueInt<T>*>(inVoid);
    ValueInt<T>* outData = static_cast<ValueInt<T>*>(outVoid);
    for( int j=0; j<*length; ++j )
    {
        const T inVal = inData[j].value;
        const T outVal = outData[j].value;
        const Int inInd = inData[j].index;
        const Int outInd = outData[j].index; 
        if( inVal > outVal || (inVal == outVal && inInd < outInd) )
            outData[j] = inData[j];
    }
}
template void
MaxLocFunc<Int>( void* in, void* out, int* length, Datatype* datatype );
template void
MaxLocFunc<float>( void* in, void* out, int* length, Datatype* datatype );
template void
MaxLocFunc<double>( void* in, void* out, int* length, Datatype* datatype );
#ifdef EL_HAVE_QUAD
template void
MaxLocFunc<Quad>( void* in, void* out, int* length, Datatype* datatype );
#endif

template<typename T>
static void
MaxLocPairFunc( void* inVoid, void* outVoid, int* length, Datatype* datatype )
{           
    const ValueIntPair<T>* inData = static_cast<ValueIntPair<T>*>(inVoid);
    ValueIntPair<T>* outData = static_cast<ValueIntPair<T>*>(outVoid);
    for( int j=0; j<*length; ++j )
    {
        const T inVal = inData[j].value;
        const T outVal = outData[j].value;
        const Int inInd0 = inData[j].indices[0];
        const Int inInd1 = inData[j].indices[1];
        const Int outInd0 = outData[j].indices[0];
        const Int outInd1 = outData[j].indices[1];
        const bool inIndLess = 
            ( inInd0 < outInd0 || (inInd0 == outInd0 && inInd1 < outInd1) );
        if( inVal > outVal || (inVal == outVal && inIndLess) )
            outData[j] = inData[j];
    }
}
template void
MaxLocPairFunc<Int>( void* in, void* out, int* length, Datatype* datatype );
template void
MaxLocPairFunc<float>( void* in, void* out, int* length, Datatype* datatype );
template void
MaxLocPairFunc<double>( void* in, void* out, int* length, Datatype* datatype );
#ifdef EL_HAVE_QUAD
template void
MaxLocPairFunc<Quad>( void* in, void* out, int* length, Datatype* datatype );
#endif

template<typename T>
static void
MinLocFunc( void* inVoid, void* outVoid, int* length, Datatype* datatype )
{           
    const ValueInt<T>* inData = static_cast<ValueInt<T>*>(inVoid);
    ValueInt<T>* outData = static_cast<ValueInt<T>*>(outVoid);
    for( int j=0; j<*length; ++j )
    {
        const T inVal = inData[j].value;
        const T outVal = outData[j].value;
        const Int inInd = inData[j].index;
        const Int outInd = outData[j].index; 
        if( inVal < outVal || (inVal == outVal && inInd < outInd) )
            outData[j] = inData[j];
    }
}
template void
MinLocFunc<Int>( void* in, void* out, int* length, Datatype* datatype );
template void
MinLocFunc<float>( void* in, void* out, int* length, Datatype* datatype );
template void
MinLocFunc<double>( void* in, void* out, int* length, Datatype* datatype );
#ifdef EL_HAVE_QUAD
template void
MinLocFunc<Quad>( void* in, void* out, int* length, Datatype* datatype );
#endif

template<typename T>
static void
MinLocPairFunc( void* inVoid, void* outVoid, int* length, Datatype* datatype )
{           
    const ValueIntPair<T>* inData = static_cast<ValueIntPair<T>*>(inVoid);
    ValueIntPair<T>* outData = static_cast<ValueIntPair<T>*>(outVoid);
    for( int j=0; j<*length; ++j )
    {
        const T inVal = inData[j].value;
        const T outVal = outData[j].value;
        const Int inInd0 = inData[j].indices[0];
        const Int inInd1 = inData[j].indices[1];
        const Int outInd0 = outData[j].indices[0];
        const Int outInd1 = outData[j].indices[1];
        const bool inIndLess = 
            ( inInd0 < outInd0 || (inInd0 == outInd0 && inInd1 < outInd1) );
        if( inVal < outVal || (inVal == outVal && inIndLess) )
            outData[j] = inData[j];
    }
}
template void
MinLocPairFunc<Int>( void* in, void* out, int* length, Datatype* datatype );
template void
MinLocPairFunc<float>( void* in, void* out, int* length, Datatype* datatype );
template void
MinLocPairFunc<double>( void* in, void* out, int* length, Datatype* datatype );
#ifdef EL_HAVE_QUAD
template void
MinLocPairFunc<Quad>( void* in, void* out, int* length, Datatype* datatype );
#endif

template<typename R> static Datatype& ValueIntType();
template<>
Datatype& ValueIntType<Int>()    { return ::IntIntType; }
template<>
Datatype& ValueIntType<float>()  { return ::floatIntType; }
template<>
Datatype& ValueIntType<double>() { return ::doubleIntType; }
#ifdef EL_HAVE_QUAD
template<>
Datatype& ValueIntType<Quad>()   { return ::QuadIntType; }
#endif

template<typename R> static Datatype& ValueIntPairType();
template<>
Datatype& ValueIntPairType<Int>()    { return ::IntIntPairType; }
template<>
Datatype& ValueIntPairType<float>()  { return ::floatIntPairType; }
template<>
Datatype& ValueIntPairType<double>() { return ::doubleIntPairType; }
#ifdef EL_HAVE_QUAD
template<>
Datatype& ValueIntPairType<Quad>()   { return ::QuadIntPairType; }
#endif

template<> Datatype TypeMap<byte>()          { return MPI_UNSIGNED_CHAR; }
template<> Datatype TypeMap<int>()           { return MPI_INT; }
template<> Datatype TypeMap<unsigned>()      { return MPI_UNSIGNED; }
template<> Datatype TypeMap<long int>()      { return MPI_LONG_INT; }
template<> Datatype TypeMap<long unsigned>() { return MPI_UNSIGNED_LONG; }
template<> Datatype TypeMap<long long int>()
{
#ifdef EL_HAVE_MPI_LONG_LONG
    return MPI_LONG_LONG_INT;
#else
    throw std::runtime_error("MPI_LONG_LONG_INT does not exist");
    return 0;
#endif
}
template<>
Datatype TypeMap<unsigned long long>()
{
#ifdef EL_HAVE_MPI_LONG_LONG
    return MPI_UNSIGNED_LONG_LONG;
#else
    throw std::runtime_error("MPI_UNSIGNED_LONG_LONG does not exist");
    return 0;
#endif
}

template<> Datatype TypeMap<float>()  { return MPI_FLOAT; }
template<> Datatype TypeMap<double>() { return MPI_DOUBLE; }
#ifdef EL_HAVE_QUAD
template<> Datatype TypeMap<Quad>()   { return ::QuadType; }
#endif

template<> Datatype TypeMap<Complex<float>>()  { return MPI_COMPLEX; }
template<> Datatype TypeMap<Complex<double>>() { return MPI_DOUBLE_COMPLEX; }
#ifdef EL_HAVE_QUAD
template<> Datatype TypeMap<Complex<Quad>>()   { return ::QuadComplexType; }
#endif

template<> Datatype TypeMap<ValueInt<Int>>()
{ return ValueIntType<Int>(); }
template<> Datatype TypeMap<ValueInt<float>>()
{ return ValueIntType<float>(); }
template<> Datatype TypeMap<ValueInt<double>>()
{ return ValueIntType<double>(); }
#ifdef EL_HAVE_QUAD
template<> Datatype TypeMap<ValueInt<Quad>>()
{ return ValueIntType<Quad>(); }
#endif

template<> Datatype TypeMap<ValueIntPair<Int>>()
{ return ValueIntPairType<Int>(); }
template<> Datatype TypeMap<ValueIntPair<float>>()
{ return ValueIntPairType<float>(); }
template<> Datatype TypeMap<ValueIntPair<double>>()
{ return ValueIntPairType<double>(); }
#ifdef EL_HAVE_QUAD
template<> Datatype TypeMap<ValueIntPair<Quad>>()
{ return ValueIntPairType<Quad>(); }
#endif

template<typename T>
static void CreateValueIntType()
{
    DEBUG_ONLY(CallStackEntry cse("CreateValueIntType"))
    Datatype typeList[2];
    typeList[0] = TypeMap<T>();
    typeList[1] = TypeMap<Int>();
    
    int blockLengths[2];
    blockLengths[0] = 1;
    blockLengths[1] = 1; 

    ValueInt<T> v;
    MPI_Aint startAddr, valueAddr, indexAddr;
    MPI_Get_address( &v,       &startAddr );
    MPI_Get_address( &v.value, &valueAddr );
    MPI_Get_address( &v.index, &indexAddr );

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
#ifdef EL_HAVE_QUAD
template void CreateValueIntType<Quad>();
#endif

template<typename T>
static void CreateValueIntPairType()
{
    DEBUG_ONLY(CallStackEntry cse("CreateValueIntPairType"))
    Datatype typeList[2];
    typeList[0] = TypeMap<T>();
    typeList[1] = TypeMap<Int>();
    
    int blockLengths[2];
    blockLengths[0] = 1;
    blockLengths[1] = 2; 

    ValueIntPair<T> v;
    MPI_Aint startAddr, valueAddr, indexAddr;
    MPI_Get_address( &v,        &startAddr );
    MPI_Get_address( &v.value,  &valueAddr );
    MPI_Get_address( v.indices, &indexAddr );

    MPI_Aint displs[2];
    displs[0] = valueAddr - startAddr;
    displs[1] = indexAddr - startAddr;

    Datatype& type = ValueIntPairType<T>();
    MPI_Type_create_struct( 2, blockLengths, displs, typeList, &type );
    MPI_Type_commit( &type );
}
template void CreateValueIntPairType<Int>();
template void CreateValueIntPairType<float>();
template void CreateValueIntPairType<double>();
#ifdef EL_HAVE_QUAD
template void CreateValueIntPairType<Quad>();
#endif

void CreateCustom()
{
    // Create the necessary types
    // ==========================
#ifdef EL_HAVE_QUAD
    // Create an MPI type for Quad
    // ---------------------------
    MPI_Type_contiguous( 2, MPI_DOUBLE, &::QuadType );
    MPI_Type_commit( &::QuadType );

    // Create an MPI type for Complex<Quad>
    // ------------------------------------
    MPI_Type_contiguous( 4, MPI_DOUBLE, &::QuadComplexType );
    MPI_Type_commit( &::QuadComplexType );
#endif
    // A value and an integer
    // ----------------------
    mpi::CreateValueIntType<Int>();
    mpi::CreateValueIntType<float>();
    mpi::CreateValueIntType<double>();
#ifdef EL_HAVE_QUAD
    mpi::CreateValueIntType<Quad>();
#endif
    // A triplet of a value and a pair of integers
    // -------------------------------------------
    mpi::CreateValueIntPairType<Int>();
    mpi::CreateValueIntPairType<float>();
    mpi::CreateValueIntPairType<double>();
#ifdef EL_HAVE_QUAD
    mpi::CreateValueIntPairType<Quad>();
#endif

    // Create the necessary MPI operations
    // ===================================
    // Functions for scalar types
    // --------------------------
#ifdef EL_HAVE_QUAD
    Create( (UserFunction*)MaxQuad, true, ::maxQuadOp );
    Create( (UserFunction*)MinQuad, true, ::minQuadOp );
    Create( (UserFunction*)SumQuad, true, ::sumQuadOp );
    Create( (UserFunction*)SumQuadComplex, true, ::sumQuadComplexOp );
#endif
    // Functions for the value and integer
    // -----------------------------------
    Create( (UserFunction*)MaxLocFunc<Int>,    true, ::maxLocIntOp    );
    Create( (UserFunction*)MinLocFunc<Int>,    true, ::minLocIntOp    );
    Create( (UserFunction*)MaxLocFunc<float>,  true, ::maxLocFloatOp  );
    Create( (UserFunction*)MinLocFunc<float>,  true, ::minLocFloatOp  );
    Create( (UserFunction*)MaxLocFunc<double>, true, ::maxLocDoubleOp );
    Create( (UserFunction*)MinLocFunc<double>, true, ::minLocDoubleOp );
#ifdef EL_HAVE_QUAD
    Create( (UserFunction*)MaxLocFunc<Quad>,   true, ::maxLocQuadOp   );
    Create( (UserFunction*)MinLocFunc<Quad>,   true, ::minLocQuadOp   );
#endif
    // Functions for the triplet of a value and a pair of integers
    // -----------------------------------------------------------
    Create( (UserFunction*)MaxLocPairFunc<Int>,    true, ::maxLocPairIntOp    );
    Create( (UserFunction*)MinLocPairFunc<Int>,    true, ::minLocPairIntOp    );
    Create( (UserFunction*)MaxLocPairFunc<float>,  true, ::maxLocPairFloatOp  );
    Create( (UserFunction*)MinLocPairFunc<float>,  true, ::minLocPairFloatOp  );
    Create( (UserFunction*)MaxLocPairFunc<double>, true, ::maxLocPairDoubleOp );
    Create( (UserFunction*)MinLocPairFunc<double>, true, ::minLocPairDoubleOp );
#ifdef EL_HAVE_QUAD
    Create( (UserFunction*)MaxLocPairFunc<Quad>,   true, ::maxLocPairQuadOp   );
    Create( (UserFunction*)MinLocPairFunc<Quad>,   true, ::minLocPairQuadOp   );
#endif
}

void DestroyCustom()
{
    // Destroy the created types
    // =========================
#ifdef EL_HAVE_QUAD
    Free( ::QuadType );
    Free( ::QuadComplexType );
#endif

    Free( ValueIntType<Int>() );
    Free( ValueIntType<float>() );
    Free( ValueIntType<double>() );
#ifdef EL_HAVE_QUAD
    Free( ValueIntType<Quad>() );
#endif

    Free( ValueIntPairType<Int>() );
    Free( ValueIntPairType<float>() );
    Free( ValueIntPairType<double>() );
#ifdef EL_HAVE_QUAD
    Free( ValueIntPairType<Quad>() );
#endif

    // Destroy the created operations
    // ==============================
#ifdef EL_HAVE_QUAD
    Free( ::maxQuadOp );
    Free( ::minQuadOp );
    Free( ::sumQuadOp );
    Free( ::sumQuadComplexOp );
#endif

    Free( ::maxLocIntOp );
    Free( ::maxLocFloatOp );
    Free( ::maxLocDoubleOp );
#ifdef EL_HAVE_QUAD
    Free( ::maxLocQuadOp );
#endif

    Free( ::maxLocPairIntOp );
    Free( ::maxLocPairFloatOp );
    Free( ::maxLocPairDoubleOp );
#ifdef EL_HAVE_QUAD
    Free( ::maxLocPairQuadOp );
#endif

    Free( ::minLocIntOp );
    Free( ::minLocFloatOp );
    Free( ::minLocDoubleOp );
#ifdef EL_HAVE_QUAD
    Free( ::minLocQuadOp );
#endif

    Free( ::minLocPairIntOp );
    Free( ::minLocPairFloatOp );
    Free( ::minLocPairDoubleOp );
#ifdef EL_HAVE_QUAD
    Free( ::minLocPairQuadOp );
#endif
}

#ifdef EL_HAVE_QUAD
template<> Op MaxOp<Quad>() { return ::maxQuadOp; }
template<> Op MinOp<Quad>() { return ::maxQuadOp; }

template<> Op SumOp<Quad>()          { return ::sumQuadOp; }
template<> Op SumOp<Complex<Quad>>() { return ::sumQuadComplexOp; }
#endif

template<> Op MaxLocOp<Int>()    { return ::maxLocIntOp; }
template<> Op MaxLocOp<float>()  { return ::maxLocFloatOp; }
template<> Op MaxLocOp<double>() { return ::maxLocDoubleOp; }
#ifdef EL_HAVE_QUAD
template<> Op MaxLocOp<Quad>()   { return ::maxLocQuadOp; }
#endif

template<> Op MaxLocPairOp<Int>()    { return ::maxLocPairIntOp; }
template<> Op MaxLocPairOp<float>()  { return ::maxLocPairFloatOp; }
template<> Op MaxLocPairOp<double>() { return ::maxLocPairDoubleOp; }
#ifdef EL_HAVE_QUAD
template<> Op MaxLocPairOp<Quad>()   { return ::maxLocPairQuadOp; }
#endif

template<> Op MinLocOp<Int>()    { return ::minLocIntOp; }
template<> Op MinLocOp<float>()  { return ::minLocFloatOp; }
template<> Op MinLocOp<double>() { return ::minLocDoubleOp; }
#ifdef EL_HAVE_QUAD
template<> Op MinLocOp<Quad>()   { return ::minLocQuadOp; }
#endif

template<> Op MinLocPairOp<Int>()    { return ::minLocPairIntOp; }
template<> Op MinLocPairOp<float>()  { return ::minLocPairFloatOp; }
template<> Op MinLocPairOp<double>() { return ::minLocPairDoubleOp; }
#ifdef EL_HAVE_QUAD
template<> Op MinLocPairOp<Quad>()   { return ::minLocPairQuadOp; }
#endif

} // namespace mpi
} // namespace El
