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
#ifdef EL_HAVE_QUAD
El::mpi::Datatype QuadType, QuadComplexType;
#endif

El::mpi::Datatype IntIntType, floatIntType, doubleIntType,
                  floatComplexIntType, doubleComplexIntType;
#ifdef EL_HAVE_QUAD
El::mpi::Datatype QuadIntType, QuadComplexIntType;
#endif

El::mpi::Datatype IntEntryType, floatEntryType, doubleEntryType,
                  floatComplexEntryType, doubleComplexEntryType;
#ifdef EL_HAVE_QUAD
El::mpi::Datatype QuadEntryType, QuadComplexEntryType;
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
    const Entry<T>* inData = static_cast<Entry<T>*>(inVoid);
    Entry<T>* outData = static_cast<Entry<T>*>(outVoid);
    for( int k=0; k<*length; ++k )
    {
        const Entry<T>& in = inData[k];
        Entry<T>& out = outData[k];
        bool inIndLess = ( in.i < out.i || (in.i == out.i && in.j < out.j) );
        if( in.value > out.value || (in.value == out.value && inIndLess) )
            out = in;
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
    const Entry<T>* inData = static_cast<Entry<T>*>(inVoid);
    Entry<T>* outData = static_cast<Entry<T>*>(outVoid);
    for( int k=0; k<*length; ++k )
    {
        const Entry<T>& in = inData[k];
        Entry<T>& out = outData[k];
        bool inIndLess = ( in.i < out.i || (in.i == out.i && in.j < out.j) );
        if( in.value < out.value || (in.value == out.value && inIndLess) )
            out = in;
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
template<>
Datatype& ValueIntType<Complex<float>>()  { return ::floatComplexIntType; }
template<>
Datatype& ValueIntType<Complex<double>>() { return ::doubleComplexIntType; }
#ifdef EL_HAVE_QUAD
template<>
Datatype& ValueIntType<Complex<Quad>>()   { return ::QuadComplexIntType; }
#endif

template<typename R> static Datatype& EntryType();
template<>
Datatype& EntryType<Int>()    { return ::IntEntryType; }
template<>
Datatype& EntryType<float>()  { return ::floatEntryType; }
template<>
Datatype& EntryType<double>() { return ::doubleEntryType; }
#ifdef EL_HAVE_QUAD
template<>
Datatype& EntryType<Quad>()   { return ::QuadEntryType; }
#endif
template<>
Datatype& EntryType<Complex<float>>() { return ::floatComplexEntryType; }
template<>
Datatype& EntryType<Complex<double>>() { return ::doubleComplexEntryType; }
#ifdef EL_HAVE_QUAD
template<>
Datatype& EntryType<Complex<Quad>>() { return ::QuadComplexEntryType; }
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

/* I'm not sure of whether it is better to manually implement these
   or not. MPI_COMPLEX and MPI_DOUBLE_COMPLEX are dangerous since it 
   appears that recent versions of MPICH leave them as NULL when 
   compiling with Clang. 

   It also appears that certain versions of OpenMPI do not support 
   MPI_C_FLOAT_COMPLEX and MPI_C_DOUBLE_COMPLEX, and so we will, for now,
   use these by default and fall back to MPI_COMPLEX and 
   MPI_DOUBLE_COMPLEX otherwise. */
template<> Datatype TypeMap<Complex<float>>()  
{ 
#ifdef EL_HAVE_MPI_C_COMPLEX
    return MPI_C_FLOAT_COMPLEX; 
#else
    return MPI_COMPLEX;
#endif
}
template<> Datatype TypeMap<Complex<double>>() 
{ 
#ifdef EL_HAVE_MPI_C_COMPLEX
    return MPI_C_DOUBLE_COMPLEX; 
#else
    return MPI_DOUBLE_COMPLEX;
#endif
}
#ifdef EL_HAVE_QUAD
template<> Datatype TypeMap<Complex<Quad>>() { return ::QuadComplexType; }
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
template<> Datatype TypeMap<ValueInt<Complex<float>>>()
{ return ValueIntType<Complex<float>>(); }
template<> Datatype TypeMap<ValueInt<Complex<double>>>()
{ return ValueIntType<Complex<double>>(); }
#ifdef EL_HAVE_QUAD
template<> Datatype TypeMap<ValueInt<Complex<Quad>>>()
{ return ValueIntType<Complex<Quad>>(); }
#endif

template<> Datatype TypeMap<Entry<Int>>()
{ return EntryType<Int>(); }
template<> Datatype TypeMap<Entry<float>>()
{ return EntryType<float>(); }
template<> Datatype TypeMap<Entry<double>>()
{ return EntryType<double>(); }
#ifdef EL_HAVE_QUAD
template<> Datatype TypeMap<Entry<Quad>>()
{ return EntryType<Quad>(); }
#endif
template<> Datatype TypeMap<Entry<Complex<float>>>()
{ return EntryType<Complex<float>>(); }
template<> Datatype TypeMap<Entry<Complex<double>>>()
{ return EntryType<Complex<double>>(); }
#ifdef EL_HAVE_QUAD
template<> Datatype TypeMap<Entry<Complex<Quad>>>()
{ return EntryType<Complex<Quad>>(); }
#endif

template<typename T>
static void CreateValueIntType()
{
    DEBUG_ONLY(CSE cse("CreateValueIntType"))
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
template void CreateValueIntType<Complex<float>>();
template void CreateValueIntType<Complex<double>>();
#ifdef EL_HAVE_QUAD
template void CreateValueIntType<Complex<Quad>>();
#endif

template<typename T>
static void CreateEntryType()
{
    DEBUG_ONLY(CSE cse("CreateEntryType"))
    Datatype typeList[3];
    typeList[0] = TypeMap<Int>();
    typeList[1] = TypeMap<Int>();
    typeList[2] = TypeMap<T>();
    
    int blockLengths[3];
    blockLengths[0] = 1;
    blockLengths[1] = 1; 
    blockLengths[2] = 1; 

    Entry<T> v;
    MPI_Aint startAddr, iAddr, jAddr, valueAddr;
    MPI_Get_address( &v,       &startAddr );
    MPI_Get_address( &v.i,     &iAddr );
    MPI_Get_address( &v.j,     &jAddr );
    MPI_Get_address( &v.value, &valueAddr );

    MPI_Aint displs[3];
    displs[0] = iAddr - startAddr;
    displs[1] = jAddr - startAddr;
    displs[2] = valueAddr - startAddr;

    Datatype& type = EntryType<T>();
    MPI_Type_create_struct( 3, blockLengths, displs, typeList, &type );
    MPI_Type_commit( &type );
}
template void CreateEntryType<Int>();
template void CreateEntryType<float>();
template void CreateEntryType<double>();
#ifdef EL_HAVE_QUAD
template void CreateEntryType<Quad>();
#endif
template void CreateEntryType<Complex<float>>();
template void CreateEntryType<Complex<double>>();
#ifdef EL_HAVE_QUAD
template void CreateEntryType<Complex<Quad>>();
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
    mpi::CreateValueIntType<Complex<float>>();
    mpi::CreateValueIntType<Complex<double>>();
#ifdef EL_HAVE_QUAD
    mpi::CreateValueIntType<Complex<Quad>>();
#endif
    // A triplet of a value and a pair of integers
    // -------------------------------------------
    mpi::CreateEntryType<Int>();
    mpi::CreateEntryType<float>();
    mpi::CreateEntryType<double>();
#ifdef EL_HAVE_QUAD
    mpi::CreateEntryType<Quad>();
#endif
    mpi::CreateEntryType<Complex<float>>();
    mpi::CreateEntryType<Complex<double>>();
#ifdef EL_HAVE_QUAD
    mpi::CreateEntryType<Complex<Quad>>();
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
    Free( ValueIntType<Complex<float>>() );
    Free( ValueIntType<Complex<double>>() );
#ifdef EL_HAVE_QUAD
    Free( ValueIntType<Complex<Quad>>() );
#endif

    Free( EntryType<Int>() );
    Free( EntryType<float>() );
    Free( EntryType<double>() );
#ifdef EL_HAVE_QUAD
    Free( EntryType<Quad>() );
#endif
    Free( EntryType<Complex<float>>() );
    Free( EntryType<Complex<double>>() );
#ifdef EL_HAVE_QUAD
    Free( EntryType<Complex<Quad>>() );
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
template<> Op MinOp<Quad>() { return ::minQuadOp; }

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
