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

// Provide default datatype (it should not end up being used)
template<typename T>
Datatype Types<T>::type = MPI_DATATYPE_NULL;

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

template<typename T>
bool Types<T>::createdTypeBeforeResize = false;
template<typename T>
Datatype Types<T>::typeBeforeResize = MPI_DATATYPE_NULL;
template<typename T>
bool Types<T>::createdType = false;

template<typename T>
bool Types<T>::haveSumOp = true;
template<typename T>
bool Types<T>::createdSumOp = false;
template<typename T>
Op Types<T>::sumOp = MPI_SUM;

template<typename T>
bool Types<T>::haveProdOp = true;
template<typename T>
bool Types<T>::createdProdOp = false;
template<typename T>
Op Types<T>::prodOp = MPI_PROD;

template<typename T>
bool Types<T>::haveMaxOp = true;
template<typename T>
bool Types<T>::createdMaxOp = false;
template<typename T>
Op Types<T>::maxOp = MPI_MAX;

template<typename T>
bool Types<T>::haveMinOp = true;
template<typename T>
bool Types<T>::createdMinOp = false;
template<typename T>
Op Types<T>::minOp = MPI_MIN;

template<typename T>
bool Types<T>::haveUserOp = true;
template<typename T>
bool Types<T>::createdUserOp = false;
template<typename T>
Op Types<T>::userOp = MPI_SUM;

template<typename T>
bool Types<T>::haveUserCommOp = true;
template<typename T>
bool Types<T>::createdUserCommOp = false;
template<typename T>
Op Types<T>::userCommOp = MPI_SUM;

template<typename T> function<T(const T&,const T&)> Types<T>::userFunc;
template<typename T> function<T(const T&,const T&)> Types<T>::userCommFunc;

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

// TODO(poulson): ValueInt<Real> user functions and ops
// TODO(poulson): ValueIntPair<Real> user functions and ops

template<typename T>
void Types<T>::Destroy()
{
    if( createdUserCommOp )
    {
        Free( userCommOp );
        haveUserCommOp = false;
        createdUserCommOp = false;
    }
    if( createdUserOp )
    {
        Free( userOp );
        haveUserOp = false;
        createdUserOp = false;
    }
    if( createdMaxOp )
    {
        Free( maxOp );
        haveMaxOp = false;
        createdMaxOp = false;
    }
    if( createdMinOp )
    {
        Free( minOp );
        haveMinOp = false;
        createdMinOp = false;
    }
    if( createdProdOp )
    {
        Free( prodOp );
        haveProdOp = false;
        createdProdOp = false;
    }
    if( createdSumOp )
    {
        Free( sumOp );
        haveSumOp = false;
        createdSumOp = false;
    }
    if( createdType )
    {
        EL_DEBUG_ONLY(
          if( type == MPI_DATATYPE_NULL )
              LogicError
              ("Types<",TypeName<T>(),">::type == MPI_DATATYPE_NULL in Free");
        )
        Free( type );
        createdType = false;
    }
    if( createdTypeBeforeResize )
    {
        EL_DEBUG_ONLY(
          if( typeBeforeResize == MPI_DATATYPE_NULL )
              LogicError
              ("Types<",TypeName<T>(),
               ">::typeBeforeResize == MPI_DATATYPE_NULL in Free");
        )
        Free( typeBeforeResize );
        createdTypeBeforeResize = false;
    }
}

namespace {

// The first stage; to be finalized via 'CreateResized'.
template<typename T>
void CreateStruct
( int numItems,
  int* blockLengths,
  MPI_Aint* displs,
  Datatype* typeList )
{
    int err =
      MPI_Type_create_struct
      ( numItems, blockLengths, displs, typeList, &Types<T>::typeBeforeResize );
    if( err != MPI_SUCCESS )
        RuntimeError("MPI_Type_create_struct returned with err=",err);
    Types<T>::createdTypeBeforeResize = true;
}

template<typename T>
void Commit()
{
    // Ensure that we can communicate using this datatype
    int err = MPI_Type_commit( &Types<T>::type );
    if( err != MPI_SUCCESS )
        RuntimeError("MPI_Type_commit returned with err=",err);
}

template<typename T>
void CreateResized( MPI_Aint lowerBound, MPI_Aint extent )
{
    int err =
      MPI_Type_create_resized
      ( Types<T>::typeBeforeResize, 0, extent, &Types<T>::type );
    if( err != MPI_SUCCESS )
        RuntimeError("MPI_Type_create_resized returned with err=",err);
    Types<T>::createdType = true;
    Commit<T>();
}

template<typename T>
void CreateContiguous( int numItems, Datatype type )
{
    int err = MPI_Type_contiguous( numItems, type, &Types<T>::type );
    if( err != MPI_SUCCESS )
        RuntimeError("MPI_Type_contiguous returned with err=",err);
    Types<T>::createdType = true;
    Commit<T>();
}

template<typename T>
void FreeResized()
{
    if( Types<T>::createdType )
    {
        if( Types<T>::type == MPI_DATATYPE_NULL )
            LogicError
            ("Types<",TypeName<T>(),">::type == MPI_DATATYPE_NULL in Free");
        Free( Types<T>::type );
        Types<T>::createdType = false;
    }
    if( Types<T>::createdTypeBeforeResize )
    {
        if( Types<T>::typeBeforeResize == MPI_DATATYPE_NULL )
            LogicError
            ("Types<",TypeName<T>(),
             ">::typeBeforeResize == MPI_DATATYPE_NULL in Free");
        Free( Types<T>::typeBeforeResize );
        Types<T>::createdTypeBeforeResize = false;
    }
}

template<typename T>
void FreeResizedFamily()
{
    FreeResized<Entry<T>>();
    FreeResized<ValueInt<T>>();
    FreeResized<T>();
}

template<typename Real>
void FreeResizedScalarFamily()
{
    FreeResizedFamily<Complex<Real>>();
    FreeResizedFamily<Real>();
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

    CreateStruct<BigInt>( 2, blockLengths, displs, typeList );
    CreateResized<BigInt>( 0, packedSize );
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

    CreateStruct<BigFloat>( 4, blockLengths, displs, typeList );
    CreateResized<BigFloat>( 0, packedSize );
}
#endif // ifdef EL_HAVE_MPC

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
ProdFunc( void* inVoid, void* outVoid, int* lengthPtr, Datatype* datatype )
EL_NO_EXCEPT
{
    auto inData  = static_cast<const T*>(inVoid);
    auto outData = static_cast<      T*>(outVoid);
    const int length = *lengthPtr;
    for( int j=0; j<length; ++j )
        outData[j] *= inData[j];
}
template<typename T,typename=DisableIf<IsPacked<T>>,typename=void>
static void
ProdFunc( void* inVoid, void* outVoid, int* lengthPtr, Datatype* datatype )
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

        b *= a;
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
    EL_DEBUG_CSE

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

    MPI_Aint extent = sizeof(v);

    CreateStruct<ValueInt<T>>( 2, blockLengths, displs, typeList );
    CreateResized<ValueInt<T>>( 0, extent );
}

template<typename T,typename=DisableIf<IsPacked<T>>,typename=void>
void CreateValueIntType() EL_NO_EXCEPT
{
    EL_DEBUG_CSE

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

    MPI_Aint extent = packedSize + sizeof(Int);

    CreateStruct<ValueInt<T>>( 2, blockLengths, displs, typeList );
    CreateResized<ValueInt<T>>( 0, extent );
}

template<typename T,typename=EnableIf<IsPacked<T>>>
static void CreateEntryType() EL_NO_EXCEPT
{
    EL_DEBUG_CSE

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

    MPI_Aint extent = sizeof(v);

    CreateStruct<Entry<T>>( 3, blockLengths, displs, typeList );
    CreateResized<Entry<T>>( 0, extent );
}

template<typename T,typename=DisableIf<IsPacked<T>>,typename=void>
void CreateEntryType() EL_NO_EXCEPT
{
    EL_DEBUG_CSE

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

    T alpha;
    const auto packedSize = alpha.SerializedSize();
    MPI_Aint extent = 2*sizeof(Int) + packedSize;

    CreateStruct<Entry<T>>( 3, blockLengths, displs, typeList );
    CreateResized<Entry<T>>( 0, extent );
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
    CreateContiguous<Complex<BigFloat>>( 2, Types<BigFloat>::type );

    CreateValueIntType<BigFloat>();
    CreateValueIntType<Complex<BigFloat>>();

    CreateEntryType<BigFloat>();
    CreateEntryType<Complex<BigFloat>>();
}

void DestroyBigIntFamily()
{
    FreeResizedFamily<BigInt>();
}

void DestroyBigFloatFamily()
{
    FreeResizedScalarFamily<BigFloat>();
}
#endif

template<typename T>
void CreateUserOps()
{
    Create( (UserFunction*)UserReduce<T>, false, UserOp<T>() );
    Create( (UserFunction*)UserReduceComm<T>, true, UserCommOp<T>() );
    Types<T>::createdUserOp = true;
    Types<T>::createdUserCommOp = true;
}

template<typename T>
void CreateSumOp()
{
    Create( (UserFunction*)SumFunc<T>, true, SumOp<T>() );
    Types<T>::createdSumOp = true;
}
template<typename T>
void CreateProdOp()
{
    Create( (UserFunction*)ProdFunc<T>, true, ProdOp<T>() );
    Types<T>::createdProdOp = true;
}
template<typename T>
void CreateMaxOp()
{
    Create( (UserFunction*)MaxFunc<T>, true, MaxOp<T>() );
    Types<T>::createdMaxOp = true;
}
template<typename T>
void CreateMinOp()
{
    Create( (UserFunction*)MinFunc<T>, true, MinOp<T>() );
    Types<T>::createdMinOp = true;
}

template<typename T>
void CreateMaxLocOp()
{
    Create( (UserFunction*)MaxLocFunc<T>, true, MaxLocOp<T>() );
    Types<ValueInt<T>>::createdMaxOp = true;
}
template<typename T>
void CreateMinLocOp()
{
    Create( (UserFunction*)MinLocFunc<T>, true, MinLocOp<T>() );
    Types<ValueInt<T>>::createdMinOp = true;
}

template<typename T>
void CreateMaxLocPairOp()
{
    Create( (UserFunction*)MaxLocPairFunc<T>, true, MaxLocPairOp<T>() );
    Types<Entry<T>>::createdMaxOp = true;
}
template<typename T>
void CreateMinLocPairOp()
{
    Create( (UserFunction*)MinLocPairFunc<T>, true, MinLocPairOp<T>() );
    Types<Entry<T>>::createdMinOp = true;
}

void CreateCustom() EL_NO_RELEASE_EXCEPT
{
    // Int
    // ===
    CreateValueIntType<Int>();
    CreateEntryType<Int>();
    CreateUserOps<Int>();
    CreateMaxLocOp<Int>();
    CreateMinLocOp<Int>();
    CreateMaxLocPairOp<Int>();
    CreateMinLocPairOp<Int>();

    // float
    // =====
#ifdef EL_USE_64BIT_INTS
    CreateValueIntType<float>();
#else
    TypeMap<ValueInt<float>>() = MPI_FLOAT_INT;
#endif
    CreateValueIntType<Complex<float>>();
    CreateEntryType<float>();
    CreateEntryType<Complex<float>>();
    CreateUserOps<float>();
    CreateUserOps<Complex<float>>();
#ifdef EL_USE_64BIT_INTS
    CreateMaxLocOp<float>();
    CreateMinLocOp<float>();
#else
    MaxLocOp<float>() = MAXLOC;
    MinLocOp<float>() = MINLOC;
#endif
    CreateMaxLocPairOp<float>();
    CreateMinLocPairOp<float>();

    // double
    // ======
#ifdef EL_USE_64BIT_INTS
    CreateValueIntType<double>();
#else
    TypeMap<ValueInt<double>>() = MPI_DOUBLE_INT;
#endif
    CreateValueIntType<Complex<double>>();
    CreateEntryType<double>();
    CreateEntryType<Complex<double>>();
    CreateUserOps<double>();
    CreateUserOps<Complex<double>>();
#ifdef EL_USE_64BIT_INTS
    CreateMaxLocOp<double>();
    CreateMinLocOp<double>();
#else
    MaxLocOp<double>() = MAXLOC;
    MinLocOp<double>() = MINLOC;
#endif
    CreateMaxLocPairOp<double>();
    CreateMinLocPairOp<double>();

#ifdef EL_HAVE_QD
    // DoubleDouble
    // ============
    CreateContiguous<DoubleDouble>( 2, MPI_DOUBLE );
    CreateContiguous<Complex<DoubleDouble>>( 2, TypeMap<DoubleDouble>() );
    CreateValueIntType<DoubleDouble>();
    CreateValueIntType<Complex<DoubleDouble>>();
    CreateEntryType<DoubleDouble>();
    CreateEntryType<Complex<DoubleDouble>>();
    CreateUserOps<DoubleDouble>();
    CreateUserOps<Complex<DoubleDouble>>();
    CreateMaxOp<DoubleDouble>();
    CreateMinOp<DoubleDouble>();
    CreateSumOp<DoubleDouble>();
    CreateProdOp<DoubleDouble>();
    CreateSumOp<Complex<DoubleDouble>>();
    CreateProdOp<Complex<DoubleDouble>>();
    CreateMaxLocOp<DoubleDouble>();
    CreateMinLocOp<DoubleDouble>();
    CreateMaxLocPairOp<DoubleDouble>();
    CreateMinLocPairOp<DoubleDouble>();

    // QuadDouble
    // ==========
    CreateContiguous<QuadDouble>( 4, MPI_DOUBLE );
    CreateContiguous<Complex<QuadDouble>>( 2, TypeMap<QuadDouble>() );
    CreateValueIntType<QuadDouble>();
    CreateValueIntType<Complex<QuadDouble>>();
    CreateEntryType<QuadDouble>();
    CreateEntryType<Complex<QuadDouble>>();
    CreateUserOps<QuadDouble>();
    CreateUserOps<Complex<QuadDouble>>();
    CreateMaxOp<QuadDouble>();
    CreateMinOp<QuadDouble>();
    CreateSumOp<QuadDouble>();
    CreateProdOp<QuadDouble>();
    CreateSumOp<Complex<QuadDouble>>();
    CreateProdOp<Complex<QuadDouble>>();
    CreateMaxLocOp<QuadDouble>();
    CreateMinLocOp<QuadDouble>();
    CreateMaxLocPairOp<QuadDouble>();
    CreateMinLocPairOp<QuadDouble>();
#endif

#ifdef EL_HAVE_QUAD
    // Quad
    // ====
    CreateContiguous<Quad>( 2, MPI_DOUBLE );
    CreateContiguous<Complex<Quad>>( 2, TypeMap<Quad>() );
    CreateValueIntType<Quad>();
    CreateValueIntType<Complex<Quad>>();
    CreateEntryType<Quad>();
    CreateEntryType<Complex<Quad>>();
    CreateUserOps<Quad>();
    CreateUserOps<Complex<Quad>>();
    CreateMaxOp<Quad>();
    CreateMinOp<Quad>();
    CreateSumOp<Quad>();
    CreateProdOp<Quad>();
    CreateSumOp<Complex<Quad>>();
    CreateProdOp<Complex<Quad>>();
    CreateMaxLocOp<Quad>();
    CreateMinLocOp<Quad>();
    CreateMaxLocPairOp<Quad>();
    CreateMinLocPairOp<Quad>();
#endif

#ifdef EL_HAVE_MPC
    // BigFloat
    // ========
    // NOTE: The BigFloat types are created by mpfr::SetPrecision previously
    //       within El::Initialize
    CreateUserOps<BigFloat>();
    CreateUserOps<Complex<BigFloat>>();
    CreateMaxOp<BigFloat>();
    CreateMinOp<BigFloat>();
    CreateSumOp<BigFloat>();
    CreateProdOp<BigFloat>();
    CreateSumOp<Complex<BigFloat>>();
    CreateProdOp<Complex<BigFloat>>();
    CreateMaxLocOp<BigFloat>();
    CreateMinLocOp<BigFloat>();
    CreateMaxLocPairOp<BigFloat>();
    CreateMinLocPairOp<BigFloat>();

    // BigInt
    // ======
    CreateUserOps<BigInt>();
    CreateMaxOp<BigInt>();
    CreateMinOp<BigInt>();
    CreateSumOp<BigInt>();
    CreateProdOp<BigInt>();
    CreateMaxLocOp<BigInt>();
    CreateMinLocOp<BigInt>();
    CreateMaxLocPairOp<BigInt>();
    CreateMinLocPairOp<BigInt>();
#endif
}

template<typename T>
void DestroyFamily()
{
    Types<Entry<T>>::Destroy();
    Types<ValueInt<T>>::Destroy();
    Types<T>::Destroy();
}

template<typename T>
void DestroyScalarFamily()
{
    DestroyFamily<Complex<T>>();
    DestroyFamily<T>();
}

void DestroyCustom() EL_NO_RELEASE_EXCEPT
{
    DestroyFamily<Int>();
    DestroyScalarFamily<float>();
    DestroyScalarFamily<double>();
#ifdef EL_HAVE_QD
    DestroyScalarFamily<DoubleDouble>();
    DestroyScalarFamily<QuadDouble>();
#endif
#ifdef EL_HAVE_QUAD
    DestroyScalarFamily<Quad>();
#endif
#ifdef EL_HAVE_MPC
    DestroyScalarFamily<BigFloat>();
    DestroyFamily<BigInt>();
#endif
}

} // namespace mpi
} // namespace El
