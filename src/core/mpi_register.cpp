/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El-lite.hpp>
using std::function;

namespace El {

namespace mpi {

// Datatypes
// =========

// Provide default datatypes (though they should not end up being used)

template<typename T>
Datatype Types<T>::type = MPI_UNSIGNED_CHAR;

template<typename T> Op Types<T>::sumOp = MPI_SUM;
template<typename T> Op Types<T>::maxOp = MPI_MAX;
template<typename T> Op Types<T>::minOp = MPI_MIN;
template<typename T> Op Types<T>::userOp = MPI_SUM;
template<typename T> Op Types<T>::userCommOp = MPI_SUM;

template<typename T> function<T(const T&,const T&)> Types<T>::userFunc;
template<typename T> function<T(const T&,const T&)> Types<T>::userCommFunc;

template<> Datatype Types<byte>::type = MPI_UNSIGNED_CHAR;
template<> Datatype Types<short>::type = MPI_SHORT;
template<> Datatype Types<int>::type = MPI_INT;
template<> Datatype Types<unsigned>::type = MPI_UNSIGNED;
template<> Datatype Types<long int>::type = MPI_LONG_INT;
template<> Datatype Types<unsigned long>::type = MPI_UNSIGNED_LONG;
#ifdef EL_HAVE_MPI_LONG_LONG
template<> Datatype Types<long long int>::type = MPI_LONG_LONG_INT;
template<> Datatype Types<unsigned long long>::type = MPI_UNSIGNED_LONG_LONG;
#endif
template<> Datatype Types<float>::type = MPI_FLOAT;
template<> Datatype Types<double>::type = MPI_DOUBLE;
/* I'm not sure of whether it is better to manually implement these
   or not. MPI_COMPLEX and MPI_DOUBLE_COMPLEX are dangerous since it 
   appears that recent versions of MPICH leave them as NULL when 
   compiling with Clang. 

   It also appears that certain versions of OpenMPI do not support 
   MPI_C_FLOAT_COMPLEX and MPI_C_DOUBLE_COMPLEX, and so we will, for now,
   use these by default and fall back to MPI_COMPLEX and 
   MPI_DOUBLE_COMPLEX otherwise. */
#ifdef EL_HAVE_MPI_C_COMPLEX
template<> Datatype Types<Complex<float>>::type = MPI_C_FLOAT_COMPLEX;
template<> Datatype Types<Complex<double>>::type = MPI_C_DOUBLE_COMPLEX;
#else
template<> Datatype Types<Complex<float>>::type = MPI_COMPLEX;
template<> Datatype Types<Complex<double>>::type = MPI_DOUBLE_COMPLEX;
#endif

template struct Types<byte>;
template struct Types<short>;
template struct Types<unsigned>;
template struct Types<unsigned long>;
#ifdef EL_USE_64BIT_INTS
template struct Types<int>; // Avoid conflict with Int
#endif
template struct Types<long int>;
#ifdef EL_HAVE_MPI_LONG_LONG
template struct Types<unsigned long long>;
#ifndef EL_USE_64BIT_INTS
template struct Types<long long int>; // Avoid conflict with Int
#endif
#endif

#define PROTO(T) \
  template struct Types<T>; \
  template struct Types<ValueInt<T>>; \
  template struct Types<Entry<T>>;

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

// TODO: ValueInt<Real> user functions and ops
// TODO: ValueIntPair<Real> user functions and ops

namespace {

// The first stage; to be finalized via 'CreateResized'.
// TODO: Decide whether the intermediate type (i.e., 'newType') needs to be
//       committed and/or freed.
void CreateStruct
( int numItems,
  int* blockLengths,
  MPI_Aint* displs,
  Datatype* typeList,
  Datatype& newType )
{
    int err =
      MPI_Type_create_struct
      ( numItems, blockLengths, displs, typeList, &newType );
    if( err != MPI_SUCCESS )
        RuntimeError("MPI_Type_create_struct returned with err=",err);
}

void CreateResized
( Datatype& type,
  MPI_Aint lowerBound,
  MPI_Aint extent,
  Datatype& newType )
{
    int err = MPI_Type_create_resized( type, 0, extent, &newType );
    if( err != MPI_SUCCESS )
        RuntimeError("MPI_Type_create_resized returned with err=",err);

    err = MPI_Type_commit( &newType );
    if( err != MPI_SUCCESS )
        RuntimeError("MPI_Type_commit returned with err=",err);
}

void CreateContiguous
( int numItems,
  Datatype type,
  Datatype& newType )
{
    int err = MPI_Type_contiguous( numItems, type, &newType );
    if( err != MPI_SUCCESS )
        RuntimeError("MPI_Type_contiguous returned with err=",err);

    err = MPI_Type_commit( &newType );
    if( err != MPI_SUCCESS )
        RuntimeError("MPI_Type_commit returned with err=",err);
}

} // anonymous namespace

#ifdef EL_HAVE_MPC

void CreateBigIntType()
{
    BigInt alpha;
    const auto packedSize = alpha.SerializedSize();
    const int numLimbs = mpfr::NumIntLimbs();

    Datatype typeList[2];
    typeList[0] = TypeMap<int>();
    typeList[1] = TypeMap<mp_limb_t>();
    
    int blockLengths[2];
    blockLengths[0] = 1;
    blockLengths[1] = numLimbs;

    MPI_Aint displs[2];
    displs[0] = 0;
    displs[1] = sizeof(int);
     
    MPI_Datatype tempType;
    CreateStruct( 2, blockLengths, displs, typeList, tempType );
    CreateResized( tempType, 0, packedSize, TypeMap<BigInt>() );
}

void CreateBigFloatType()
{
    BigFloat alpha;
    const auto packedSize = alpha.SerializedSize();
    const auto numLimbs = alpha.NumLimbs();

    Datatype typeList[4];
    typeList[0] = TypeMap<mpfr_prec_t>();
    typeList[1] = TypeMap<mpfr_sign_t>();
    typeList[2] = TypeMap<mpfr_exp_t>();
    typeList[3] = TypeMap<mp_limb_t>();
    
    int blockLengths[4];
    blockLengths[0] = 1;
    blockLengths[1] = 1; 
    blockLengths[2] = 1;
    blockLengths[3] = numLimbs;

    MPI_Aint displs[4];
    displs[0] = 0;
    displs[1] = sizeof(mpfr_prec_t);
    displs[2] = sizeof(mpfr_prec_t) + sizeof(mpfr_sign_t);
    displs[3] = sizeof(mpfr_prec_t) + sizeof(mpfr_sign_t) + sizeof(mpfr_exp_t);
    
    MPI_Datatype tempType;
    CreateStruct( 4, blockLengths, displs, typeList, tempType );
    CreateResized( tempType, 0, packedSize, TypeMap<BigFloat>() );
}
#endif // ifdef EL_HAVE_MPC

// TODO: Expose hooks for these routines so that the user could run
//
//  AllReduce
//  ( value,
//    []( const T& alpha, const T& beta ) { return Min(alpha,beta); }
//    comm )
//

template<typename T,typename=EnableIf<IsPacked<T>>>
static void
UserReduce
( void* inVoid, void* outVoid, int* lengthPtr, Datatype* datatype )
EL_NO_EXCEPT
{
    auto inData  = static_cast<const T*>(inVoid);
    auto outData = static_cast<      T*>(outVoid);
    const int length = *lengthPtr;
    auto func = Types<T>::userFunc;
    for( int j=0; j<length; ++j )
        outData[j] = func(inData[j],outData[j]);
}
template<typename T,typename=DisableIf<IsPacked<T>>,typename=void>
static void
UserReduce
( void* inVoid, void* outVoid, int* lengthPtr, Datatype* datatype )
EL_NO_EXCEPT
{
    T a, b;
    auto inData  = static_cast<const byte*>(inVoid);
    auto outData = static_cast<      byte*>(outVoid);
    const int length = *lengthPtr;
    auto func = Types<T>::userFunc;
    for( int j=0; j<length; ++j )
    {
        inData = a.Deserialize(inData);
        b.Deserialize(outData);

        b = func(a,b);
        outData = b.Serialize(outData); 
    }
}

template<typename T,typename=EnableIf<IsPacked<T>>>
static void
UserReduceComm
( void* inVoid, void* outVoid, int* lengthPtr, Datatype* datatype )
EL_NO_EXCEPT
{
    auto inData  = static_cast<const T*>(inVoid);
    auto outData = static_cast<      T*>(outVoid);
    const int length = *lengthPtr;
    auto func = Types<T>::userCommFunc;
    for( int j=0; j<length; ++j )
        outData[j] = func(inData[j],outData[j]);
}
template<typename T,typename=DisableIf<IsPacked<T>>,typename=void>
static void
UserReduceComm
( void* inVoid, void* outVoid, int* lengthPtr, Datatype* datatype )
EL_NO_EXCEPT
{
    T a, b;
    auto inData  = static_cast<const byte*>(inVoid);
    auto outData = static_cast<      byte*>(outVoid);
    const int length = *lengthPtr;
    auto func = Types<T>::userCommFunc;
    for( int j=0; j<length; ++j )
    {
        inData = a.Deserialize(inData);
        b.Deserialize(outData);

        b = func(a,b);
        outData = b.Serialize(outData); 
    }
}

template<typename T,typename=EnableIf<IsPacked<T>>>
static void
MaxFunc( void* inVoid, void* outVoid, int* lengthPtr, Datatype* datatype )
EL_NO_EXCEPT
{
    auto inData  = static_cast<const T*>(inVoid);
    auto outData = static_cast<      T*>(outVoid);
    const int length = *lengthPtr;
    for( int j=0; j<length; ++j )
    {
        if( inData[j] > outData[j] )
            outData[j] = inData[j];
    }
}
template<typename T,typename=DisableIf<IsPacked<T>>,typename=void>
static void
MaxFunc( void* inVoid, void* outVoid, int* lengthPtr, Datatype* datatype )
EL_NO_EXCEPT
{
    T a, b;
    auto inData  = static_cast<const byte*>(inVoid);
    auto outData = static_cast<      byte*>(outVoid);
    const int length = *lengthPtr;
    for( int j=0; j<length; ++j )
    {
        inData = a.Deserialize(inData);
        auto bAfter = b.Deserialize(outData);

        if( a > b )
            a.Serialize(outData);
        outData = bAfter;
    }
}

template<typename T,typename=EnableIf<IsPacked<T>>>
static void
MinFunc( void* inVoid, void* outVoid, int* lengthPtr, Datatype* datatype )
EL_NO_EXCEPT
{
    auto inData  = static_cast<const T*>(inVoid);
    auto outData = static_cast<      T*>(outVoid);
    const int length = *lengthPtr;
    for( int j=0; j<length; ++j )
    {
        if( inData[j] < outData[j] )
            outData[j] = inData[j];
    }
}
template<typename T,typename=DisableIf<IsPacked<T>>,typename=void>
static void
MinFunc( void* inVoid, void* outVoid, int* lengthPtr, Datatype* datatype )
EL_NO_EXCEPT
{
    T a, b;
    auto inData  = static_cast<const byte*>(inVoid);
    auto outData = static_cast<      byte*>(outVoid);
    const int length = *lengthPtr;
    for( int j=0; j<length; ++j )
    {
        inData = a.Deserialize(inData);
        auto bAfter = b.Deserialize(outData);

        if( a < b )
            a.Serialize(outData);
        outData = bAfter;
    }
}

template<typename T,typename=EnableIf<IsPacked<T>>>
static void
SumFunc( void* inVoid, void* outVoid, int* lengthPtr, Datatype* datatype )
EL_NO_EXCEPT
{
    auto inData  = static_cast<const T*>(inVoid);
    auto outData = static_cast<      T*>(outVoid);
    const int length = *lengthPtr;
    for( int j=0; j<length; ++j )
        outData[j] += inData[j];
}
template<typename T,typename=DisableIf<IsPacked<T>>,typename=void>
static void
SumFunc( void* inVoid, void* outVoid, int* lengthPtr, Datatype* datatype )
EL_NO_EXCEPT
{
    T a, b;
    auto inData  = static_cast<const byte*>(inVoid);
    auto outData = static_cast<      byte*>(outVoid);
    const int length = *lengthPtr;
    for( int j=0; j<length; ++j )
    {
        inData = a.Deserialize(inData);
        b.Deserialize(outData);

        b += a;
        outData = b.Serialize(outData); 
    }
}

template<typename T,typename=EnableIf<IsPacked<T>>>
static void
MaxLocFunc( void* inVoid, void* outVoid, int* lengthPtr, Datatype* datatype )
EL_NO_EXCEPT
{           
    auto inData  = static_cast<const ValueInt<T>*>(inVoid);
    auto outData = static_cast<      ValueInt<T>*>(outVoid);
    const int length = *lengthPtr;
    for( int j=0; j<length; ++j )
    {
        const T inVal = inData[j].value;
        const T outVal = outData[j].value;
        const Int inInd = inData[j].index;
        const Int outInd = outData[j].index; 
        if( inVal > outVal || (inVal == outVal && inInd < outInd) )
            outData[j] = inData[j];
    }
}
template<typename T,typename=DisableIf<IsPacked<T>>,typename=void>
static void
MaxLocFunc( void* inVoid, void* outVoid, int* lengthPtr, Datatype* datatype )
EL_NO_EXCEPT
{
    ValueInt<T> a, b;
    auto inData  = static_cast<const byte*>(inVoid);
    auto outData = static_cast<      byte*>(outVoid);
    const int length = *lengthPtr;
    for( int j=0; j<length; ++j )
    {
        inData = Deserialize( 1, inData,  &a );
                 Deserialize( 1, outData, &b );

        if( a.value > b.value || (a.value == b.value && a.index < b.index) )
            outData = Serialize( 1, &a, outData );
        else
            outData += a.value.SerializedSize();
    }
}

template<typename T,typename=EnableIf<IsPacked<T>>>
static void
MinLocFunc( void* inVoid, void* outVoid, int* lengthPtr, Datatype* datatype )
EL_NO_EXCEPT
{           
    auto inData  = static_cast<const ValueInt<T>*>(inVoid);
    auto outData = static_cast<      ValueInt<T>*>(outVoid);
    const int length = *lengthPtr;
    for( int j=0; j<length; ++j )
    {
        const T inVal = inData[j].value;
        const T outVal = outData[j].value;
        const Int inInd = inData[j].index;
        const Int outInd = outData[j].index; 
        if( inVal < outVal || (inVal == outVal && inInd < outInd) )
            outData[j] = inData[j];
    }
}
template<typename T,typename=DisableIf<IsPacked<T>>,typename=void>
static void
MinLocFunc( void* inVoid, void* outVoid, int* lengthPtr, Datatype* datatype )
EL_NO_EXCEPT
{
    ValueInt<T> a, b;
    auto inData  = static_cast<const byte*>(inVoid);
    auto outData = static_cast<      byte*>(outVoid);
    const int length = *lengthPtr;
    for( int j=0; j<length; ++j )
    {
        inData = Deserialize( 1, inData,  &a );
                 Deserialize( 1, outData, &b );

        if( a.value < b.value || (a.value == b.value && a.index < b.index) )
            outData = Serialize( 1, &a, outData );
        else
            outData += a.value.SerializedSize();
    }
}

template<typename T,typename=EnableIf<IsPacked<T>>>
static void
MaxLocPairFunc
( void* inVoid, void* outVoid, int* lengthPtr, Datatype* datatype )
EL_NO_EXCEPT
{
    auto inData  = static_cast<const Entry<T>*>(inVoid);
    auto outData = static_cast<      Entry<T>*>(outVoid);
    const int length = *lengthPtr;
    for( int k=0; k<length; ++k )
    {
        const Entry<T>& in  = inData[k];
              Entry<T>& out = outData[k];
        bool inIndLess = ( in.i < out.i || (in.i == out.i && in.j < out.j) );
        if( in.value > out.value || (in.value == out.value && inIndLess) )
            out = in;
    }
}
template<typename T,typename=DisableIf<IsPacked<T>>,typename=void>
static void
MaxLocPairFunc
( void* inVoid, void* outVoid, int* lengthPtr, Datatype* datatype )
EL_NO_EXCEPT
{
    Entry<T> a, b;
    auto inData  = static_cast<const byte*>(inVoid);
    auto outData = static_cast<      byte*>(outVoid);
    const int length = *lengthPtr;
    for( int j=0; j<length; ++j )
    {
        inData = Deserialize( 1, inData,  &a );
                 Deserialize( 1, outData, &b );

        bool inIndLess = ( a.i < b.i || (a.i == b.i && a.j < b.j) );
        if( a.value > b.value || (a.value == b.value && inIndLess) )
            outData = Serialize( 1, &a, outData );
        else
            outData += a.value.SerializedSize();
    }
}

template<typename T,typename=EnableIf<IsPacked<T>>>
static void
MinLocPairFunc
( void* inVoid, void* outVoid, int* lengthPtr, Datatype* datatype )
EL_NO_EXCEPT
{
    auto inData  = static_cast<const Entry<T>*>(inVoid);
    auto outData = static_cast<      Entry<T>*>(outVoid);
    const int length = *lengthPtr;
    for( int k=0; k<length; ++k )
    {
        const Entry<T>& in  = inData[k];
              Entry<T>& out = outData[k];
        bool inIndLess = ( in.i < out.i || (in.i == out.i && in.j < out.j) );
        if( in.value < out.value || (in.value == out.value && inIndLess) )
            out = in;
    }
}
template<typename T,typename=DisableIf<IsPacked<T>>,typename=void>
static void
MinLocPairFunc
( void* inVoid, void* outVoid, int* lengthPtr, Datatype* datatype )
EL_NO_EXCEPT
{
    Entry<T> a, b;
    auto inData  = static_cast<const byte*>(inVoid);
    auto outData = static_cast<      byte*>(outVoid);
    const int length = *lengthPtr;
    for( int j=0; j<length; ++j )
    {
        inData = Deserialize( 1, inData,  &a );
                 Deserialize( 1, outData, &b );

        bool inIndLess = ( a.i < b.i || (a.i == b.i && a.j < b.j) );
        if( a.value < b.value || (a.value == b.value && inIndLess) )
            outData = Serialize( 1, &a, outData );
        else
            outData += a.value.SerializedSize();
    }
}

template<typename T,typename=EnableIf<IsPacked<T>>>
static void CreateValueIntType() EL_NO_EXCEPT
{
    DEBUG_CSE

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

    Datatype tempType;
    CreateStruct( 2, blockLengths, displs, typeList, tempType );

    MPI_Aint extent = sizeof(v);
    CreateResized( tempType, 0, extent, TypeMap<ValueInt<T>>() );
}

template<typename T,typename=DisableIf<IsPacked<T>>,typename=void>
void CreateValueIntType() EL_NO_EXCEPT
{
    DEBUG_CSE

    Datatype typeList[2];
    typeList[0] = TypeMap<T>();
    typeList[1] = TypeMap<Int>();
    
    int blockLengths[2];
    blockLengths[0] = 1;
    blockLengths[1] = 1; 

    T alpha;
    const size_t packedSize = alpha.SerializedSize();

    MPI_Aint displs[2];
    displs[0] = 0;
    displs[1] = packedSize;

    Datatype tempType;
    CreateStruct( 2, blockLengths, displs, typeList, tempType );

    MPI_Aint extent = packedSize + sizeof(Int);
    CreateResized( tempType, 0, extent, TypeMap<ValueInt<T>>() );
}

template<typename T,typename=EnableIf<IsPacked<T>>>
static void CreateEntryType() EL_NO_EXCEPT
{
    DEBUG_CSE

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

    Datatype tempType;
    CreateStruct( 3, blockLengths, displs, typeList, tempType );

    MPI_Aint extent = sizeof(v);
    CreateResized( tempType, 0, extent, TypeMap<Entry<T>>() );
}

template<typename T,typename=DisableIf<IsPacked<T>>,typename=void>
void CreateEntryType() EL_NO_EXCEPT
{
    DEBUG_CSE

    Datatype typeList[3];
    typeList[0] = TypeMap<Int>();
    typeList[1] = TypeMap<Int>();
    typeList[2] = TypeMap<T>();
    
    int blockLengths[3];
    blockLengths[0] = 1;
    blockLengths[1] = 1; 
    blockLengths[2] = 1; 

    MPI_Aint displs[3];
    displs[0] = 0;
    displs[1] = sizeof(Int);
    displs[2] = 2*sizeof(Int);

    Datatype tempType;
    CreateStruct( 3, blockLengths, displs, typeList, tempType );

    T alpha;
    const auto packedSize = alpha.SerializedSize();

    MPI_Aint extent = 2*sizeof(Int) + packedSize;
    CreateResized( tempType, 0, extent, TypeMap<Entry<T>>() );
}

#ifdef EL_HAVE_MPC
void CreateBigIntFamily()
{
    CreateBigIntType();
    CreateValueIntType<BigInt>();
    CreateEntryType<BigInt>();
}

void CreateBigFloatFamily()
{
    CreateBigFloatType();
    CreateContiguous
    ( 2, Types<BigFloat>::type, Types<Complex<BigFloat>>::type );

    CreateValueIntType<BigFloat>();
    CreateValueIntType<Complex<BigFloat>>();

    CreateEntryType<BigFloat>();
    CreateEntryType<Complex<BigFloat>>();
}

void DestroyBigIntFamily()
{
    Free( Types<Entry<BigInt>>::type );
    Free( Types<ValueInt<BigInt>>::type );
    Free( Types<BigInt>::type );
}

void DestroyBigFloatFamily()
{

    Free( Types<Entry<Complex<BigFloat>>>::type );
    Free( Types<ValueInt<Complex<BigFloat>>>::type );
    Free( Types<Complex<BigFloat>>::type );

    Free( Types<Entry<BigFloat>>::type );
    Free( Types<ValueInt<BigFloat>>::type );
    Free( Types<BigFloat>::type );
}
#endif

template<typename T>
void CreateUserOps()
{
    Create( (UserFunction*)UserReduce<T>, false, UserOp<T>() );
    Create( (UserFunction*)UserReduceComm<T>, true, UserCommOp<T>() );
}

template<typename T>
void CreateSumOp()
{
    Create( (UserFunction*)SumFunc<T>, true, SumOp<T>() );
}
template<typename T>
void CreateMaxOp()
{
    Create( (UserFunction*)MaxFunc<T>, true, MaxOp<T>() );
}
template<typename T>
void CreateMinOp()
{
    Create( (UserFunction*)MinFunc<T>, true, MinOp<T>() );
}

template<typename T>
void CreateMaxLocOp()
{
    Create( (UserFunction*)MaxLocFunc<T>, true, MaxLocOp<T>() );
}
template<typename T>
void CreateMinLocOp()
{
    Create( (UserFunction*)MinLocFunc<T>, true, MinLocOp<T>() );
}

template<typename T>
void CreateMaxLocPairOp()
{
    Create( (UserFunction*)MaxLocPairFunc<T>, true, MaxLocPairOp<T>() );
}
template<typename T>
void CreateMinLocPairOp()
{
    Create( (UserFunction*)MinLocPairFunc<T>, true, MinLocPairOp<T>() );
}

void CreateCustom() EL_NO_RELEASE_EXCEPT
{
    // Create the necessary types
    // ==========================
#ifdef EL_HAVE_QD
    CreateContiguous( 2, MPI_DOUBLE, TypeMap<DoubleDouble>() );
    CreateContiguous( 4, MPI_DOUBLE, TypeMap<QuadDouble>() );
    CreateContiguous
    ( 2, TypeMap<DoubleDouble>(), TypeMap<Complex<DoubleDouble>>() );
    CreateContiguous
    ( 2, TypeMap<QuadDouble>(), TypeMap<Complex<QuadDouble>>() );
#endif
#ifdef EL_HAVE_QUAD
    CreateContiguous( 2, MPI_DOUBLE, TypeMap<Quad>() );
    CreateContiguous( 4, MPI_DOUBLE, TypeMap<Complex<Quad>>() );
#endif
    // NOTE: The BigFloat types are created by mpfr::SetPrecision previously
    //       within El::Initialize

    // A value and an integer
    // ----------------------
    CreateValueIntType<Int>();
#ifdef EL_USE_64BIT_INTS
    CreateValueIntType<float>();
    CreateValueIntType<double>();
#else
    TypeMap<ValueInt<float>>() = MPI_FLOAT_INT;
    TypeMap<ValueInt<double>>() = MPI_DOUBLE_INT;
#endif
    CreateValueIntType<Complex<float>>();
    CreateValueIntType<Complex<double>>();
#ifdef EL_HAVE_QD
    CreateValueIntType<DoubleDouble>();
    CreateValueIntType<QuadDouble>();
    CreateValueIntType<Complex<DoubleDouble>>();
    CreateValueIntType<Complex<QuadDouble>>();
#endif
#ifdef EL_HAVE_QUAD
    CreateValueIntType<Quad>();
    CreateValueIntType<Complex<Quad>>();
#endif

    // A triplet of a value and a pair of integers
    // -------------------------------------------
    CreateEntryType<Int>();
    CreateEntryType<float>();
    CreateEntryType<double>();
    CreateEntryType<Complex<float>>();
    CreateEntryType<Complex<double>>();
#ifdef EL_HAVE_QD
    CreateEntryType<DoubleDouble>();
    CreateEntryType<QuadDouble>();
    CreateEntryType<Complex<DoubleDouble>>();
    CreateEntryType<Complex<QuadDouble>>();
#endif
#ifdef EL_HAVE_QUAD
    CreateEntryType<Quad>();
    CreateEntryType<Complex<Quad>>();
#endif

    // Create the necessary MPI operations
    // ===================================
    // Functions for user-defined ops
    // ------------------------------
    CreateUserOps<Int>();
    CreateUserOps<float>();
    CreateUserOps<double>();
    CreateUserOps<Complex<float>>();
    CreateUserOps<Complex<double>>();
#ifdef EL_HAVE_QD
    CreateUserOps<DoubleDouble>();
    CreateUserOps<QuadDouble>();
    CreateUserOps<Complex<DoubleDouble>>();
    CreateUserOps<Complex<QuadDouble>>();
#endif
#ifdef EL_HAVE_QUAD
    CreateUserOps<Quad>();
    CreateUserOps<Complex<Quad>>();
#endif
#ifdef EL_HAVE_MPC
    CreateUserOps<BigInt>();
    CreateUserOps<BigFloat>();
    CreateUserOps<Complex<BigFloat>>();
#endif
   
    // Functions for scalar types
    // --------------------------
#ifdef EL_HAVE_QD
    CreateMaxOp<DoubleDouble>();
    CreateMinOp<DoubleDouble>();
    CreateSumOp<DoubleDouble>();

    CreateMaxOp<QuadDouble>();
    CreateMinOp<QuadDouble>();
    CreateSumOp<QuadDouble>();

    CreateSumOp<Complex<DoubleDouble>>();
    CreateSumOp<Complex<QuadDouble>>();
#endif
#ifdef EL_HAVE_QUAD
    CreateMaxOp<Quad>();
    CreateMinOp<Quad>();
    CreateSumOp<Quad>();

    CreateSumOp<Complex<Quad>>();
#endif
#ifdef EL_HAVE_MPC
    CreateMaxOp<BigInt>();
    CreateMinOp<BigInt>();
    CreateSumOp<BigInt>();

    CreateMaxOp<BigFloat>();
    CreateMinOp<BigFloat>();
    CreateSumOp<BigFloat>();

    CreateSumOp<Complex<BigFloat>>();
#endif
    // Functions for the value and integer
    // -----------------------------------
    CreateMaxLocOp<Int>();
    CreateMinLocOp<Int>();
#ifdef EL_USE_64BIT_INTS
    CreateMaxLocOp<float>();
    CreateMinLocOp<float>();

    CreateMaxLocOp<double>();
    CreateMinLocOp<double>();
#else
    MaxLocOp<float>() = MAXLOC;
    MinLocOp<float>() = MINLOC;

    MaxLocOp<double>() = MAXLOC;
    MinLocOp<double>() = MINLOC;
#endif
#ifdef EL_HAVE_QD
    CreateMaxLocOp<DoubleDouble>();
    CreateMinLocOp<DoubleDouble>();

    CreateMaxLocOp<QuadDouble>();
    CreateMinLocOp<QuadDouble>();
#endif
#ifdef EL_HAVE_QUAD
    CreateMaxLocOp<Quad>();
    CreateMinLocOp<Quad>();
#endif
#ifdef EL_HAVE_MPC
    CreateMaxLocOp<BigInt>();
    CreateMinLocOp<BigInt>();

    CreateMaxLocOp<BigFloat>();
    CreateMinLocOp<BigFloat>();
#endif

    // Functions for the triplet of a value and a pair of integers
    // -----------------------------------------------------------
    CreateMaxLocPairOp<Int>();
    CreateMinLocPairOp<Int>();

    CreateMaxLocPairOp<float>();
    CreateMinLocPairOp<float>();

    CreateMaxLocPairOp<double>();
    CreateMinLocPairOp<double>();
#ifdef EL_HAVE_QD
    CreateMaxLocPairOp<DoubleDouble>();
    CreateMinLocPairOp<DoubleDouble>();

    CreateMaxLocPairOp<QuadDouble>();
    CreateMinLocPairOp<QuadDouble>();
#endif
#ifdef EL_HAVE_QUAD
    CreateMaxLocPairOp<Quad>();
    CreateMinLocPairOp<Quad>();
#endif
#ifdef EL_HAVE_MPC
    CreateMaxLocPairOp<BigInt>();
    CreateMinLocPairOp<BigInt>();

    CreateMaxLocPairOp<BigFloat>();
    CreateMinLocPairOp<BigFloat>();
#endif
}

template<typename T>
void FreeUserOps()
{
    Free( Types<T>::userOp );
    Free( Types<T>::userCommOp );
}

template<typename T>
void FreeMaxMinOps()
{
    Free( Types<T>::maxOp );
    Free( Types<T>::minOp );
}

template<typename T>
void FreeScalarOps()
{
    if( IsReal<T>::value )
    {
        FreeMaxMinOps<T>();
        FreeMaxMinOps<ValueInt<T>>();
        FreeMaxMinOps<Entry<T>>();
    }
    Free( Types<T>::sumOp );
}

void DestroyCustom() EL_NO_RELEASE_EXCEPT
{
    // Destroy the created operations
    // ==============================

    // User-defined operations
    // -----------------------
    FreeUserOps<Int>();
    FreeUserOps<float>();
    FreeUserOps<double>();
    FreeUserOps<Complex<float>>();
    FreeUserOps<Complex<double>>();
#ifdef EL_HAVE_QD
    FreeUserOps<DoubleDouble>();
    FreeUserOps<QuadDouble>();
    FreeUserOps<Complex<DoubleDouble>>();
    FreeUserOps<Complex<QuadDouble>>();
#endif
#ifdef EL_HAVE_QUAD
    FreeUserOps<Quad>();
    FreeUserOps<Complex<Quad>>();
#endif
#ifdef EL_HAVE_MPC
    FreeUserOps<BigInt>();
    FreeUserOps<BigFloat>();
    FreeUserOps<Complex<BigFloat>>();
#endif

#ifdef EL_HAVE_QD
    FreeScalarOps<DoubleDouble>();
    FreeScalarOps<QuadDouble>();
    FreeScalarOps<Complex<DoubleDouble>>();
    FreeScalarOps<Complex<QuadDouble>>();
#endif
#ifdef EL_HAVE_QUAD
    FreeScalarOps<Quad>();
    FreeScalarOps<Complex<Quad>>();
#endif
#ifdef EL_HAVE_MPC
    FreeScalarOps<BigInt>();
    FreeScalarOps<BigFloat>();
    FreeScalarOps<Complex<BigFloat>>();
#endif

    FreeMaxMinOps<ValueInt<Int>>();
#ifdef EL_USE_64BIT_INTS
    FreeMaxMinOps<ValueInt<float>>();
    FreeMaxMinOps<ValueInt<double>>();
#endif
    FreeMaxMinOps<Entry<Int>>();
    FreeMaxMinOps<Entry<float>>();
    FreeMaxMinOps<Entry<double>>();

    // Destroy the created types
    // =========================
    Free( Types<Entry<Int>>::type );
    Free( Types<Entry<float>>::type );
    Free( Types<Entry<double>>::type );
    Free( Types<Entry<Complex<float>>>::type );
    Free( Types<Entry<Complex<double>>>::type );
#ifdef EL_HAVE_QD
    Free( Types<Entry<DoubleDouble>>::type );
    Free( Types<Entry<QuadDouble>>::type );
    Free( Types<Entry<Complex<DoubleDouble>>>::type );
    Free( Types<Entry<Complex<QuadDouble>>>::type );
#endif
#ifdef EL_HAVE_QUAD
    Free( Types<Entry<Quad>>::type );
    Free( Types<Entry<Complex<Quad>>>::type );
#endif

    Free( Types<ValueInt<Int>>::type );
#ifdef EL_USE_64BIT_INTS
    Free( Types<ValueInt<float>>::type );
    Free( Types<ValueInt<double>>::type );
#endif
    Free( Types<ValueInt<Complex<float>>>::type );
    Free( Types<ValueInt<Complex<double>>>::type );
#ifdef EL_HAVE_QD
    Free( Types<ValueInt<DoubleDouble>>::type );
    Free( Types<ValueInt<QuadDouble>>::type );
    Free( Types<ValueInt<Complex<DoubleDouble>>>::type );
    Free( Types<ValueInt<Complex<QuadDouble>>>::type );
#endif
#ifdef EL_HAVE_QUAD
    Free( Types<ValueInt<Quad>>::type );
    Free( Types<ValueInt<Complex<Quad>>>::type );
#endif

#ifdef EL_HAVE_QD
    Free( Types<Complex<DoubleDouble>>::type );
    Free( Types<Complex<QuadDouble>>::type );
    Free( Types<DoubleDouble>::type );
    Free( Types<QuadDouble>::type );
#endif
#ifdef EL_HAVE_QUAD
    Free( Types<Complex<Quad>>::type );
    Free( Types<Quad>::type );
#endif
#ifdef EL_HAVE_MPC
    DestroyBigIntFamily();
    DestroyBigFloatFamily();
#endif
}

} // namespace mpi
} // namespace El
