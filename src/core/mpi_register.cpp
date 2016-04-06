/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace {

using El::Int;
#ifdef EL_HAVE_QD
using El::DoubleDouble;
using El::QuadDouble;
#endif
#ifdef EL_HAVE_QUAD
using El::Quad;
#endif
using El::Complex;
#ifdef EL_HAVE_MPC
using El::BigInt;
using El::BigFloat;
#endif
using std::function;

// Datatypes
// =========

// Scalar datatypes
// ----------------
#ifdef EL_HAVE_QD
El::mpi::Datatype DoubleDoubleType, QuadDoubleType;
#endif
#ifdef EL_HAVE_QUAD
El::mpi::Datatype QuadType, QuadComplexType;
#endif
#ifdef EL_HAVE_MPC
El::mpi::Datatype BigIntType, BigFloatType;
#endif

// (Int,Scalar) datatypes
// ----------------------
El::mpi::Datatype IntIntType, floatIntType, doubleIntType,
                  floatComplexIntType, doubleComplexIntType;
#ifdef EL_HAVE_QD
El::mpi::Datatype DoubleDoubleIntType, QuadDoubleIntType;
#endif
#ifdef EL_HAVE_QUAD
El::mpi::Datatype QuadIntType, QuadComplexIntType;
#endif
#ifdef EL_HAVE_MPC
El::mpi::Datatype BigIntIntType, BigFloatIntType;
#endif

// (Int,Int,Scalar) datatypes
// --------------------------
El::mpi::Datatype IntEntryType, floatEntryType, doubleEntryType,
                  floatComplexEntryType, doubleComplexEntryType;
#ifdef EL_HAVE_QD
El::mpi::Datatype DoubleDoubleEntryType, QuadDoubleEntryType;
#endif
#ifdef EL_HAVE_QUAD
El::mpi::Datatype QuadEntryType, QuadComplexEntryType;
#endif
#ifdef EL_HAVE_MPC
El::mpi::Datatype BigIntEntryType, BigFloatEntryType;
#endif

// Operations
// ==========

// Scalar datatype operations
// --------------------------
#ifdef EL_HAVE_QD
El::mpi::Op minDoubleDoubleOp, maxDoubleDoubleOp, sumDoubleDoubleOp;
El::mpi::Op minQuadDoubleOp, maxQuadDoubleOp, sumQuadDoubleOp;
#endif
#ifdef EL_HAVE_QUAD
El::mpi::Op minQuadOp, maxQuadOp;
El::mpi::Op sumQuadOp, sumQuadComplexOp;
#endif
#ifdef EL_HAVE_MPC
El::mpi::Op minBigIntOp, maxBigIntOp;
El::mpi::Op sumBigIntOp;

El::mpi::Op minBigFloatOp, maxBigFloatOp;
El::mpi::Op sumBigFloatOp;
#endif

// (Int,Scalar) datatype operations
// --------------------------------
El::mpi::Op minLocIntOp,    maxLocIntOp,
            minLocFloatOp,  maxLocFloatOp,
            minLocDoubleOp, maxLocDoubleOp;
#ifdef EL_HAVE_QD
El::mpi::Op minLocDoubleDoubleOp, maxLocDoubleDoubleOp;
El::mpi::Op minLocQuadDoubleOp, maxLocQuadDoubleOp;
#endif
#ifdef EL_HAVE_QUAD
El::mpi::Op minLocQuadOp, maxLocQuadOp;
#endif
#ifdef EL_HAVE_MPC
El::mpi::Op minLocBigIntOp, maxLocBigIntOp;
El::mpi::Op minLocBigFloatOp, maxLocBigFloatOp;
#endif

// (Int,Int,Scalar) datatype operations
// ------------------------------------
El::mpi::Op minLocPairIntOp,    maxLocPairIntOp,
            minLocPairFloatOp,  maxLocPairFloatOp,
            minLocPairDoubleOp, maxLocPairDoubleOp;
#ifdef EL_HAVE_QD
El::mpi::Op minLocPairDoubleDoubleOp, maxLocPairDoubleDoubleOp;
El::mpi::Op minLocPairQuadDoubleOp, maxLocPairQuadDoubleOp;
#endif
#ifdef EL_HAVE_QUAD
El::mpi::Op minLocPairQuadOp, maxLocPairQuadOp;
#endif
#ifdef EL_HAVE_MPC
El::mpi::Op minLocPairBigIntOp, maxLocPairBigIntOp;
El::mpi::Op minLocPairBigFloatOp, maxLocPairBigFloatOp;
#endif

function<Int(const Int&,const Int&)>
  userIntFunc, userIntCommFunc;
El::mpi::Op userIntOp, userIntCommOp;

function<float(const float&,const float&)>
  userFloatFunc, userFloatCommFunc;
El::mpi::Op userFloatOp, userFloatCommOp;

function<double(const double&,const double&)>
  userDoubleFunc, userDoubleCommFunc;
El::mpi::Op userDoubleOp, userDoubleCommOp;

function<Complex<float>(const Complex<float>&,const Complex<float>&)>
  userComplexFloatFunc, userComplexFloatCommFunc;
El::mpi::Op userComplexFloatOp, userComplexFloatCommOp;

function<Complex<double>(const Complex<double>&,const Complex<double>&)>
  userComplexDoubleFunc, userComplexDoubleCommFunc;
El::mpi::Op userComplexDoubleOp, userComplexDoubleCommOp;

#ifdef EL_HAVE_QD
function<DoubleDouble(const DoubleDouble&,const DoubleDouble&)>
  userDoubleDoubleFunc, userDoubleDoubleCommFunc;
El::mpi::Op userDoubleDoubleOp, userDoubleDoubleCommOp;

function<QuadDouble(const QuadDouble&,const QuadDouble&)>
  userQuadDoubleFunc, userQuadDoubleCommFunc;
El::mpi::Op userQuadDoubleOp, userQuadDoubleCommOp;
#endif

#ifdef EL_HAVE_QUAD
function<Quad(const Quad&,const Quad&)>
  userQuadFunc, userQuadCommFunc;
El::mpi::Op userQuadOp, userQuadCommOp;

function<Complex<Quad>(const Complex<Quad>&,const Complex<Quad>&)>
  userComplexQuadFunc, userComplexQuadCommFunc;
El::mpi::Op userComplexQuadOp, userComplexQuadCommOp;
#endif

#ifdef EL_HAVE_MPC
function<BigInt(const BigInt&,const BigInt&)>
  userBigIntFunc, userBigIntCommFunc;
El::mpi::Op userBigIntOp, userBigIntCommOp;

function<BigFloat(const BigFloat&,const BigFloat&)>
  userBigFloatFunc, userBigFloatCommFunc;
El::mpi::Op userBigFloatOp, userBigFloatCommOp;

// TODO: Complex BigFloat functions and ops
#endif

// TODO: ValueInt<Real> user functions and ops
// TODO: ValueIntPair<Real> user functions and ops

} // anonymous namespace   

namespace El {

namespace mpi {

namespace {

// The first stage; to be finalized via 'CreateResized'.
// TODO: Decide whether the intermediate type (i.e., 'newType') needs to be
//       committed and/or freed.
void CreateStruct
( int numItems,
  int* blockLengths,
  MPI_Aint* displs,
  mpi::Datatype* typeList,
  mpi::Datatype& newType )
{
    int err =
      MPI_Type_create_struct
      ( numItems, blockLengths, displs, typeList, &newType );
    if( err != MPI_SUCCESS )
        RuntimeError("MPI_Type_create_struct returned with err=",err);
}

void CreateResized
( mpi::Datatype& type,
  MPI_Aint lowerBound,
  MPI_Aint extent,
  mpi::Datatype& newType )
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
  mpi::Datatype type,
  mpi::Datatype& newType )
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
    const int numLimbs = mpc::NumIntLimbs();

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
    CreateResized( tempType, 0, packedSize, ::BigIntType );
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
    CreateResized( tempType, 0, packedSize, ::BigFloatType );
}
#endif // ifdef EL_HAVE_MPC

// TODO: Expose hooks for these routines so that the user could run
//
//  mpi::AllReduce
//  ( value,
//    []( const T& alpha, const T& beta ) { return Min(alpha,beta); }
//    comm )
//
template<>
void SetUserReduceFunc
( function<Int(const Int&,const Int&)> func, bool commutative )
{
    if( commutative )
        ::userIntCommFunc = func;
    else
        ::userIntFunc = func;
}
static void
UserIntReduce
( void* inVoid, void* outVoid, int* lengthPtr, Datatype* datatype )
EL_NO_EXCEPT
{
    auto inData  = static_cast<const Int*>(inVoid);
    auto outData = static_cast<      Int*>(outVoid);
    const int length = *lengthPtr;
    for( int j=0; j<length; ++j )
        outData[j] = ::userIntFunc(inData[j],outData[j]);
}
static void
UserIntReduceComm
( void* inVoid, void* outVoid, int* lengthPtr, Datatype* datatype )
EL_NO_EXCEPT
{
    auto inData  = static_cast<const Int*>(inVoid);
    auto outData = static_cast<      Int*>(outVoid);
    const int length = *lengthPtr;
    for( int j=0; j<length; ++j )
        outData[j] = ::userIntCommFunc(inData[j],outData[j]);
}

template<>
void SetUserReduceFunc
( function<float(const float&,const float&)> func, bool commutative )
{
    if( commutative )
        ::userFloatCommFunc = func;
    else
        ::userFloatFunc = func;
}
static void
UserFloatReduce
( void* inVoid, void* outVoid, int* lengthPtr, Datatype* datatype )
EL_NO_EXCEPT
{
    auto inData  = static_cast<const float*>(inVoid);
    auto outData = static_cast<      float*>(outVoid);
    const int length = *lengthPtr;
    for( int j=0; j<length; ++j )
        outData[j] = ::userFloatFunc(inData[j],outData[j]);
}
static void
UserFloatReduceComm
( void* inVoid, void* outVoid, int* lengthPtr, Datatype* datatype )
EL_NO_EXCEPT
{
    auto inData  = static_cast<const float*>(inVoid);
    auto outData = static_cast<      float*>(outVoid);
    const int length = *lengthPtr;
    for( int j=0; j<length; ++j )
        outData[j] = ::userFloatCommFunc(inData[j],outData[j]);
}

template<>
void SetUserReduceFunc
( function<double(const double&,const double&)> func, bool commutative )
{
    if( commutative )
        ::userDoubleCommFunc = func;
    else
        ::userDoubleFunc = func;
}
static void
UserDoubleReduce
( void* inVoid, void* outVoid, int* lengthPtr, Datatype* datatype )
EL_NO_EXCEPT
{
    auto inData  = static_cast<const double*>(inVoid);
    auto outData = static_cast<      double*>(outVoid);
    const int length = *lengthPtr;
    for( int j=0; j<length; ++j )
        outData[j] = ::userDoubleFunc(inData[j],outData[j]);
}
static void
UserDoubleReduceComm
( void* inVoid, void* outVoid, int* lengthPtr, Datatype* datatype )
EL_NO_EXCEPT
{
    auto inData  = static_cast<const double*>(inVoid);
    auto outData = static_cast<      double*>(outVoid);
    const int length = *lengthPtr;
    for( int j=0; j<length; ++j )
        outData[j] = ::userDoubleCommFunc(inData[j],outData[j]);
}

template<>
void SetUserReduceFunc
( function<Complex<float>(const Complex<float>&,const Complex<float>&)> func,
  bool commutative )
{
    if( commutative )
        ::userComplexFloatCommFunc = func;
    else
        ::userComplexFloatFunc = func;
}
static void
UserComplexFloatReduce
( void* inVoid, void* outVoid, int* lengthPtr, Datatype* datatype )
EL_NO_EXCEPT
{
    auto inData  = static_cast<const Complex<float>*>(inVoid);
    auto outData = static_cast<      Complex<float>*>(outVoid);
    const int length = *lengthPtr;
    for( int j=0; j<length; ++j )
        outData[j] = ::userComplexFloatFunc(inData[j],outData[j]);
}
static void
UserComplexFloatReduceComm
( void* inVoid, void* outVoid, int* lengthPtr, Datatype* datatype )
EL_NO_EXCEPT
{
    auto inData  = static_cast<const Complex<float>*>(inVoid);
    auto outData = static_cast<      Complex<float>*>(outVoid);
    const int length = *lengthPtr;
    for( int j=0; j<length; ++j )
        outData[j] = ::userComplexFloatCommFunc(inData[j],outData[j]);
}

template<>
void SetUserReduceFunc
( function<Complex<double>(const Complex<double>&,const Complex<double>&)> func,
  bool commutative )
{
    if( commutative )
        ::userComplexDoubleCommFunc = func;
    else
        ::userComplexDoubleFunc = func;
}
static void
UserComplexDoubleReduce
( void* inVoid, void* outVoid, int* lengthPtr, Datatype* datatype )
EL_NO_EXCEPT
{
    auto inData  = static_cast<const Complex<double>*>(inVoid);
    auto outData = static_cast<      Complex<double>*>(outVoid);
    const int length = *lengthPtr;
    for( int j=0; j<length; ++j )
        outData[j] = ::userComplexDoubleFunc(inData[j],outData[j]);
}
static void
UserComplexDoubleReduceComm
( void* inVoid, void* outVoid, int* lengthPtr, Datatype* datatype )
EL_NO_EXCEPT
{
    auto inData  = static_cast<const Complex<double>*>(inVoid);
    auto outData = static_cast<      Complex<double>*>(outVoid);
    const int length = *lengthPtr;
    for( int j=0; j<length; ++j )
        outData[j] = ::userComplexDoubleCommFunc(inData[j],outData[j]);
}

#ifdef EL_HAVE_QD
template<>
void SetUserReduceFunc
( function<DoubleDouble(const DoubleDouble&,const DoubleDouble&)> func,
  bool commutative )
{
    if( commutative )
        ::userDoubleDoubleCommFunc = func;
    else
        ::userDoubleDoubleFunc = func;
}
static void
UserDoubleDoubleReduce
( void* inVoid, void* outVoid, int* lengthPtr, Datatype* datatype )
EL_NO_EXCEPT
{
    auto inData  = static_cast<const DoubleDouble*>(inVoid);
    auto outData = static_cast<      DoubleDouble*>(outVoid);
    const int length = *lengthPtr;
    for( int j=0; j<length; ++j )
        outData[j] = ::userDoubleDoubleFunc(inData[j],outData[j]);
}
static void
UserDoubleDoubleReduceComm
( void* inVoid, void* outVoid, int* lengthPtr, Datatype* datatype )
EL_NO_EXCEPT
{
    auto inData  = static_cast<const DoubleDouble*>(inVoid);
    auto outData = static_cast<      DoubleDouble*>(outVoid);
    const int length = *lengthPtr;
    for( int j=0; j<length; ++j )
        outData[j] = ::userDoubleDoubleCommFunc(inData[j],outData[j]);
}

static void
MaxDoubleDouble
( void* inVoid, void* outVoid, int* lengthPtr, Datatype* datatype )
EL_NO_EXCEPT
{
    auto inData  = static_cast<const DoubleDouble*>(inVoid);
    auto outData = static_cast<      DoubleDouble*>(outVoid);
    const int length = *lengthPtr;
    for( int j=0; j<length; ++j )
    {
        if( inData[j] > outData[j] )
            outData[j] = inData[j];
    }
}

static void
MinDoubleDouble
( void* inVoid, void* outVoid, int* lengthPtr, Datatype* datatype )
EL_NO_EXCEPT
{
    auto inData  = static_cast<const DoubleDouble*>(inVoid);
    auto outData = static_cast<      DoubleDouble*>(outVoid);
    const int length = *lengthPtr;
    for( int j=0; j<length; ++j )
    {
        if( inData[j] < outData[j] )
            outData[j] = inData[j];
    }
}

static void
SumDoubleDouble
( void* inVoid, void* outVoid, int* lengthPtr, Datatype* datatype )
EL_NO_EXCEPT
{
    auto inData  = static_cast<const DoubleDouble*>(inVoid);
    auto outData = static_cast<      DoubleDouble*>(outVoid);
    const int length = *lengthPtr;
    for( int j=0; j<length; ++j )
        outData[j] += inData[j];
}

template<>
void SetUserReduceFunc
( function<QuadDouble(const QuadDouble&,const QuadDouble&)> func,
  bool commutative )
{
    if( commutative )
        ::userQuadDoubleCommFunc = func;
    else
        ::userQuadDoubleFunc = func;
}
static void
UserQuadDoubleReduce
( void* inVoid, void* outVoid, int* lengthPtr, Datatype* datatype )
EL_NO_EXCEPT
{
    auto inData  = static_cast<const QuadDouble*>(inVoid);
    auto outData = static_cast<      QuadDouble*>(outVoid);
    const int length = *lengthPtr;
    for( int j=0; j<length; ++j )
        outData[j] = ::userQuadDoubleFunc(inData[j],outData[j]);
}
static void
UserQuadDoubleReduceComm
( void* inVoid, void* outVoid, int* lengthPtr, Datatype* datatype )
EL_NO_EXCEPT
{
    auto inData  = static_cast<const QuadDouble*>(inVoid);
    auto outData = static_cast<      QuadDouble*>(outVoid);
    const int length = *lengthPtr;
    for( int j=0; j<length; ++j )
        outData[j] = ::userQuadDoubleCommFunc(inData[j],outData[j]);
}

static void
MaxQuadDouble
( void* inVoid, void* outVoid, int* lengthPtr, Datatype* datatype )
EL_NO_EXCEPT
{
    auto inData  = static_cast<const QuadDouble*>(inVoid);
    auto outData = static_cast<      QuadDouble*>(outVoid);
    const int length = *lengthPtr;
    for( int j=0; j<length; ++j )
    {
        if( inData[j] > outData[j] )
            outData[j] = inData[j];
    }
}

static void
MinQuadDouble
( void* inVoid, void* outVoid, int* lengthPtr, Datatype* datatype )
EL_NO_EXCEPT
{
    auto inData  = static_cast<const QuadDouble*>(inVoid);
    auto outData = static_cast<      QuadDouble*>(outVoid);
    const int length = *lengthPtr;
    for( int j=0; j<length; ++j )
    {
        if( inData[j] < outData[j] )
            outData[j] = inData[j];
    }
}

static void
SumQuadDouble
( void* inVoid, void* outVoid, int* lengthPtr, Datatype* datatype )
EL_NO_EXCEPT
{
    auto inData  = static_cast<const QuadDouble*>(inVoid);
    auto outData = static_cast<      QuadDouble*>(outVoid);
    const int length = *lengthPtr;
    for( int j=0; j<length; ++j )
        outData[j] += inData[j];
}
#endif // ifdef EL_HAVE_QD

#ifdef EL_HAVE_QUAD
template<>
void SetUserReduceFunc
( function<Quad(const Quad&,const Quad&)> func, bool commutative )
{
    if( commutative )
        ::userQuadCommFunc = func;
    else
        ::userQuadFunc = func;
}
static void
UserQuadReduce
( void* inVoid, void* outVoid, int* lengthPtr, Datatype* datatype )
EL_NO_EXCEPT
{
    auto inData  = static_cast<const Quad*>(inVoid);
    auto outData = static_cast<      Quad*>(outVoid);
    const int length = *lengthPtr;
    for( int j=0; j<length; ++j )
        outData[j] = ::userQuadFunc(inData[j],outData[j]);
}
static void
UserQuadReduceComm
( void* inVoid, void* outVoid, int* lengthPtr, Datatype* datatype )
EL_NO_EXCEPT
{
    auto inData  = static_cast<const Quad*>(inVoid);
    auto outData = static_cast<      Quad*>(outVoid);
    const int length = *lengthPtr;
    for( int j=0; j<length; ++j )
        outData[j] = ::userQuadCommFunc(inData[j],outData[j]);
}

template<>
void SetUserReduceFunc
( function<Complex<Quad>(const Complex<Quad>&,const Complex<Quad>&)> func,
  bool commutative )
{
    if( commutative )
        ::userComplexQuadCommFunc = func;
    else
        ::userComplexQuadFunc = func;
}
static void
UserComplexQuadReduce
( void* inVoid, void* outVoid, int* lengthPtr, Datatype* datatype )
EL_NO_EXCEPT
{
    auto inData  = static_cast<const Complex<Quad>*>(inVoid);
    auto outData = static_cast<      Complex<Quad>*>(outVoid);
    const int length = *lengthPtr;
    for( int j=0; j<length; ++j )
        outData[j] = ::userComplexQuadFunc(inData[j],outData[j]);
}
static void
UserComplexQuadReduceComm
( void* inVoid, void* outVoid, int* lengthPtr, Datatype* datatype )
EL_NO_EXCEPT
{
    auto inData  = static_cast<const Complex<Quad>*>(inVoid);
    auto outData = static_cast<      Complex<Quad>*>(outVoid);
    const int length = *lengthPtr;
    for( int j=0; j<length; ++j )
        outData[j] = ::userComplexQuadCommFunc(inData[j],outData[j]);
}

static void
MaxQuad( void* inVoid, void* outVoid, int* lengthPtr, Datatype* datatype )
EL_NO_EXCEPT
{
    auto inData  = static_cast<const Quad*>(inVoid);
    auto outData = static_cast<      Quad*>(outVoid);
    const int length = *lengthPtr;
    for( int j=0; j<length; ++j )
    {
        if( inData[j] > outData[j] )
            outData[j] = inData[j];
    }
}

static void
MinQuad( void* inVoid, void* outVoid, int* lengthPtr, Datatype* datatype )
EL_NO_EXCEPT
{
    auto inData  = static_cast<const Quad*>(inVoid);
    auto outData = static_cast<      Quad*>(outVoid);
    const int length = *lengthPtr;
    for( int j=0; j<length; ++j )
    {
        if( inData[j] < outData[j] )
            outData[j] = inData[j];
    }
}

static void
SumQuad( void* inVoid, void* outVoid, int* lengthPtr, Datatype* datatype )
EL_NO_EXCEPT
{
    auto inData  = static_cast<const Quad*>(inVoid);
    auto outData = static_cast<      Quad*>(outVoid);
    const int length = *lengthPtr;
    for( int j=0; j<length; ++j )
        outData[j] += inData[j];
}

static void
SumQuadComplex
( void* inVoid, void* outVoid, int* lengthPtr, Datatype* datatype )
EL_NO_EXCEPT
{
    auto inData  = static_cast<const Complex<Quad>*>(inVoid);
    auto outData = static_cast<      Complex<Quad>*>(outVoid);
    const int length = *lengthPtr;
    for( int j=0; j<length; ++j )
        outData[j] += inData[j];
}
#endif // ifdef EL_HAVE_QUAD

#ifdef EL_HAVE_MPC

template<>
void SetUserReduceFunc
( function<BigInt(const BigInt&,const BigInt&)> func, bool commutative )
{
    if( commutative )
        ::userBigIntCommFunc = func;
    else
        ::userBigIntFunc = func;
}

static void
UserBigIntReduce
( void* inVoid, void* outVoid, int* lengthPtr, Datatype* datatype )
EL_NO_EXCEPT
{
    BigInt a, b;
    auto inData  = static_cast<const byte*>(inVoid);
    auto outData = static_cast<      byte*>(outVoid);
    const int length = *lengthPtr;
    for( int j=0; j<length; ++j )
    {
        inData = a.Deserialize(inData);
        b.Deserialize(outData);

        b = ::userBigIntFunc(a,b);
        outData = b.Serialize(outData); 
    }
}

static void
UserBigIntReduceComm
( void* inVoid, void* outVoid, int* lengthPtr, Datatype* datatype )
EL_NO_EXCEPT
{
    BigInt a, b;
    auto inData  = static_cast<const byte*>(inVoid);
    auto outData = static_cast<      byte*>(outVoid);
    const int length = *lengthPtr;
    for( int j=0; j<length; ++j )
    {
        inData = a.Deserialize(inData);
        b.Deserialize(outData);

        b = ::userBigIntCommFunc(a,b);
        outData = b.Serialize(outData); 
    }
}

static void
MaxBigInt( void* inVoid, void* outVoid, int* lengthPtr, Datatype* datatype )
EL_NO_EXCEPT
{
    BigInt a, b;
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

static void
MinBigInt( void* inVoid, void* outVoid, int* lengthPtr, Datatype* datatype )
EL_NO_EXCEPT
{
    BigInt a, b;
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

static void
SumBigInt( void* inVoid, void* outVoid, int* lengthPtr, Datatype* datatype )
EL_NO_EXCEPT
{
    BigInt a, b;
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

template<>
void SetUserReduceFunc
( function<BigFloat(const BigFloat&,const BigFloat&)> func, bool commutative )
{
    if( commutative )
        ::userBigFloatCommFunc = func;
    else
        ::userBigFloatFunc = func;
}

static void
UserBigFloatReduce
( void* inVoid, void* outVoid, int* lengthPtr, Datatype* datatype )
EL_NO_EXCEPT
{
    BigFloat a, b;
    auto inData  = static_cast<const byte*>(inVoid);
    auto outData = static_cast<      byte*>(outVoid);
    const int length = *lengthPtr;
    for( int j=0; j<length; ++j )
    {
        inData = a.Deserialize(inData);
        b.Deserialize(outData);

        b = ::userBigFloatFunc(a,b);
        outData = b.Serialize(outData); 
    }
}

static void
UserBigFloatReduceComm
( void* inVoid, void* outVoid, int* lengthPtr, Datatype* datatype )
EL_NO_EXCEPT
{
    BigFloat a, b;
    auto inData  = static_cast<const byte*>(inVoid);
    auto outData = static_cast<      byte*>(outVoid);
    const int length = *lengthPtr;
    for( int j=0; j<length; ++j )
    {
        inData = a.Deserialize(inData);
        b.Deserialize(outData);

        b = ::userBigFloatCommFunc(a,b);
        outData = b.Serialize(outData); 
    }
}

static void
MaxBigFloat( void* inVoid, void* outVoid, int* lengthPtr, Datatype* datatype )
EL_NO_EXCEPT
{
    BigFloat a, b;
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

static void
MinBigFloat( void* inVoid, void* outVoid, int* lengthPtr, Datatype* datatype )
EL_NO_EXCEPT
{
    BigFloat a, b;
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

static void
SumBigFloat( void* inVoid, void* outVoid, int* lengthPtr, Datatype* datatype )
EL_NO_EXCEPT
{
    BigFloat a, b;
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
#endif // ifdef EL_HAVE_MPC

template<typename T>
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

#ifdef EL_HAVE_MPC
template<>
void MaxLocFunc<BigInt>
( void* inVoid, void* outVoid, int* lengthPtr, Datatype* datatype )
EL_NO_EXCEPT
{
    ValueInt<BigInt> a, b;
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

template<>
void MaxLocFunc<BigFloat>
( void* inVoid, void* outVoid, int* lengthPtr, Datatype* datatype )
EL_NO_EXCEPT
{
    ValueInt<BigFloat> a, b;
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
#endif

template void
MaxLocFunc<Int>( void* in, void* out, int* length, Datatype* datatype )
EL_NO_EXCEPT;
template void
MaxLocFunc<float>( void* in, void* out, int* length, Datatype* datatype )
EL_NO_EXCEPT;
template void
MaxLocFunc<double>( void* in, void* out, int* length, Datatype* datatype )
EL_NO_EXCEPT;
#ifdef EL_HAVE_QD
template void
MaxLocFunc<DoubleDouble>( void* in, void* out, int* length, Datatype* datatype )
EL_NO_EXCEPT;
template void
MaxLocFunc<QuadDouble>( void* in, void* out, int* length, Datatype* datatype )
EL_NO_EXCEPT;
#endif
#ifdef EL_HAVE_QUAD
template void
MaxLocFunc<Quad>( void* in, void* out, int* length, Datatype* datatype )
EL_NO_EXCEPT;
#endif

template<typename T>
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

#ifdef EL_HAVE_MPC
template<>
void MaxLocPairFunc<BigInt>
( void* inVoid, void* outVoid, int* lengthPtr, Datatype* datatype )
EL_NO_EXCEPT
{
    Entry<BigInt> a, b;
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

template<>
void MaxLocPairFunc<BigFloat>
( void* inVoid, void* outVoid, int* lengthPtr, Datatype* datatype )
EL_NO_EXCEPT
{
    Entry<BigFloat> a, b;
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
#endif

template void
MaxLocPairFunc<Int>( void* in, void* out, int* length, Datatype* datatype )
EL_NO_EXCEPT;
template void
MaxLocPairFunc<float>( void* in, void* out, int* length, Datatype* datatype )
EL_NO_EXCEPT;
template void
MaxLocPairFunc<double>( void* in, void* out, int* length, Datatype* datatype )
EL_NO_EXCEPT;
#ifdef EL_HAVE_QD
template void
MaxLocPairFunc<DoubleDouble>
( void* in, void* out, int* length, Datatype* datatype )
EL_NO_EXCEPT;
template void
MaxLocPairFunc<QuadDouble>
( void* in, void* out, int* length, Datatype* datatype )
EL_NO_EXCEPT;
#endif
#ifdef EL_HAVE_QUAD
template void
MaxLocPairFunc<Quad>( void* in, void* out, int* length, Datatype* datatype )
EL_NO_EXCEPT;
#endif

template<typename T>
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

#ifdef EL_HAVE_MPC
template<>
void MinLocFunc<BigInt>
( void* inVoid, void* outVoid, int* lengthPtr, Datatype* datatype )
EL_NO_EXCEPT
{
    ValueInt<BigInt> a, b;
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

template<>
void MinLocFunc<BigFloat>
( void* inVoid, void* outVoid, int* lengthPtr, Datatype* datatype )
EL_NO_EXCEPT
{
    ValueInt<BigFloat> a, b;
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
#endif

template void
MinLocFunc<Int>( void* in, void* out, int* length, Datatype* datatype )
EL_NO_EXCEPT;
template void
MinLocFunc<float>( void* in, void* out, int* length, Datatype* datatype )
EL_NO_EXCEPT;
template void
MinLocFunc<double>( void* in, void* out, int* length, Datatype* datatype )
EL_NO_EXCEPT;
#ifdef EL_HAVE_QD
template void
MinLocFunc<DoubleDouble>( void* in, void* out, int* length, Datatype* datatype )
EL_NO_EXCEPT;
template void
MinLocFunc<QuadDouble>( void* in, void* out, int* length, Datatype* datatype )
EL_NO_EXCEPT;
#endif
#ifdef EL_HAVE_QUAD
template void
MinLocFunc<Quad>( void* in, void* out, int* length, Datatype* datatype )
EL_NO_EXCEPT;
#endif

template<typename T>
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
        const Entry<T>& in = inData[k];
        Entry<T>& out = outData[k];
        bool inIndLess = ( in.i < out.i || (in.i == out.i && in.j < out.j) );
        if( in.value < out.value || (in.value == out.value && inIndLess) )
            out = in;
    }
}

#ifdef EL_HAVE_MPC
template<>
void MinLocPairFunc<BigInt>
( void* inVoid, void* outVoid, int* lengthPtr, Datatype* datatype )
EL_NO_EXCEPT
{
    Entry<BigInt> a, b;
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

template<>
void MinLocPairFunc<BigFloat>
( void* inVoid, void* outVoid, int* lengthPtr, Datatype* datatype )
EL_NO_EXCEPT
{
    Entry<BigFloat> a, b;
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
#endif

template void
MinLocPairFunc<Int>( void* in, void* out, int* length, Datatype* datatype )
EL_NO_EXCEPT;
template void
MinLocPairFunc<float>( void* in, void* out, int* length, Datatype* datatype )
EL_NO_EXCEPT;
template void
MinLocPairFunc<double>( void* in, void* out, int* length, Datatype* datatype )
EL_NO_EXCEPT;
#ifdef EL_HAVE_QD
template void
MinLocPairFunc<DoubleDouble>
( void* in, void* out, int* length, Datatype* datatype )
EL_NO_EXCEPT;
template void
MinLocPairFunc<QuadDouble>
( void* in, void* out, int* length, Datatype* datatype )
EL_NO_EXCEPT;
#endif
#ifdef EL_HAVE_QUAD
template void
MinLocPairFunc<Quad>( void* in, void* out, int* length, Datatype* datatype )
EL_NO_EXCEPT;
#endif

template<typename R> static Datatype& ValueIntType() EL_NO_EXCEPT;
template<>
Datatype& ValueIntType<Int>() EL_NO_EXCEPT { return ::IntIntType; }
template<>
Datatype& ValueIntType<float>() EL_NO_EXCEPT { return ::floatIntType; }
template<>
Datatype& ValueIntType<Complex<float>>() EL_NO_EXCEPT
{ return ::floatComplexIntType; }
template<>
Datatype& ValueIntType<double>() EL_NO_EXCEPT { return ::doubleIntType; }
template<>
Datatype& ValueIntType<Complex<double>>() EL_NO_EXCEPT
{ return ::doubleComplexIntType; }
#ifdef EL_HAVE_QD
template<>
Datatype& ValueIntType<DoubleDouble>() EL_NO_EXCEPT
{ return ::DoubleDoubleIntType; }
template<>
Datatype& ValueIntType<QuadDouble>() EL_NO_EXCEPT
{ return ::QuadDoubleIntType; }
#endif
#ifdef EL_HAVE_QUAD
template<>
Datatype& ValueIntType<Quad>() EL_NO_EXCEPT { return ::QuadIntType; }
template<>
Datatype& ValueIntType<Complex<Quad>>() EL_NO_EXCEPT
{ return ::QuadComplexIntType; }
#endif
#ifdef EL_HAVE_MPC
template<>
Datatype& ValueIntType<BigInt>() EL_NO_EXCEPT { return ::BigIntIntType; }
template<>
Datatype& ValueIntType<BigFloat>() EL_NO_EXCEPT { return ::BigFloatIntType; }
#endif

template<typename R> static Datatype& EntryType() EL_NO_EXCEPT;
template<>
Datatype& EntryType<Int>() EL_NO_EXCEPT { return ::IntEntryType; }
template<>
Datatype& EntryType<float>() EL_NO_EXCEPT { return ::floatEntryType; }
template<>
Datatype& EntryType<Complex<float>>() EL_NO_EXCEPT
{ return ::floatComplexEntryType; }
template<>
Datatype& EntryType<double>() EL_NO_EXCEPT { return ::doubleEntryType; }
template<>
Datatype& EntryType<Complex<double>>() EL_NO_EXCEPT
{ return ::doubleComplexEntryType; }
#ifdef EL_HAVE_QD
template<>
Datatype& EntryType<DoubleDouble>() EL_NO_EXCEPT
{ return ::DoubleDoubleEntryType; }
template<>
Datatype& EntryType<QuadDouble>() EL_NO_EXCEPT
{ return ::QuadDoubleEntryType; }
#endif
#ifdef EL_HAVE_QUAD
template<>
Datatype& EntryType<Quad>() EL_NO_EXCEPT { return ::QuadEntryType; }
template<>
Datatype& EntryType<Complex<Quad>>() EL_NO_EXCEPT
{ return ::QuadComplexEntryType; }
#endif
#ifdef EL_HAVE_MPC
template<>
Datatype& EntryType<BigInt>() EL_NO_EXCEPT { return ::BigIntEntryType; }
template<>
Datatype& EntryType<BigFloat>() EL_NO_EXCEPT { return ::BigFloatEntryType; }
#endif

template<> Datatype TypeMap<byte>() EL_NO_EXCEPT
{ return MPI_UNSIGNED_CHAR; }
template<> Datatype TypeMap<short>() EL_NO_EXCEPT
{ return MPI_SHORT; }
template<> Datatype TypeMap<int>() EL_NO_EXCEPT
{ return MPI_INT; }
template<> Datatype TypeMap<unsigned>() EL_NO_EXCEPT
{ return MPI_UNSIGNED; }
template<> Datatype TypeMap<long int>() EL_NO_EXCEPT
{ return MPI_LONG_INT; }
template<> Datatype TypeMap<long unsigned>() EL_NO_EXCEPT
{ return MPI_UNSIGNED_LONG; }

#ifdef EL_HAVE_MPI_LONG_LONG
template<> Datatype TypeMap<long long int>() EL_NO_EXCEPT
{ return MPI_LONG_LONG_INT; }
template<> Datatype TypeMap<unsigned long long>() EL_NO_EXCEPT
{ return MPI_UNSIGNED_LONG_LONG; }
#endif

template<> Datatype TypeMap<float>() EL_NO_EXCEPT { return MPI_FLOAT; }
template<> Datatype TypeMap<double>() EL_NO_EXCEPT{ return MPI_DOUBLE; }
#ifdef EL_HAVE_QD
template<> Datatype TypeMap<DoubleDouble>() EL_NO_EXCEPT
{ return ::DoubleDoubleType; }
template<> Datatype TypeMap<QuadDouble>() EL_NO_EXCEPT
{ return ::QuadDoubleType; }
#endif
#ifdef EL_HAVE_QUAD
template<> Datatype TypeMap<Quad>() EL_NO_EXCEPT { return ::QuadType; }
template<> Datatype TypeMap<Complex<Quad>>() EL_NO_EXCEPT 
{ return ::QuadComplexType; }
#endif
#ifdef EL_HAVE_MPC
template<> Datatype TypeMap<BigInt>() EL_NO_EXCEPT { return ::BigIntType; }
template<> Datatype TypeMap<BigFloat>() EL_NO_EXCEPT { return ::BigFloatType; }
// TODO: Complex<BigFloat>?
#endif

/* I'm not sure of whether it is better to manually implement these
   or not. MPI_COMPLEX and MPI_DOUBLE_COMPLEX are dangerous since it 
   appears that recent versions of MPICH leave them as NULL when 
   compiling with Clang. 

   It also appears that certain versions of OpenMPI do not support 
   MPI_C_FLOAT_COMPLEX and MPI_C_DOUBLE_COMPLEX, and so we will, for now,
   use these by default and fall back to MPI_COMPLEX and 
   MPI_DOUBLE_COMPLEX otherwise. */
template<> Datatype TypeMap<Complex<float>>() EL_NO_EXCEPT
{ 
#ifdef EL_HAVE_MPI_C_COMPLEX
    return MPI_C_FLOAT_COMPLEX; 
#else
    return MPI_COMPLEX;
#endif
}
template<> Datatype TypeMap<Complex<double>>() EL_NO_EXCEPT
{ 
#ifdef EL_HAVE_MPI_C_COMPLEX
    return MPI_C_DOUBLE_COMPLEX; 
#else
    return MPI_DOUBLE_COMPLEX;
#endif
}

template<> Datatype TypeMap<ValueInt<Int>>() EL_NO_EXCEPT
{ return ValueIntType<Int>(); }
template<> Datatype TypeMap<ValueInt<float>>() EL_NO_EXCEPT
{ return ValueIntType<float>(); }
template<> Datatype TypeMap<ValueInt<Complex<float>>>() EL_NO_EXCEPT
{ return ValueIntType<Complex<float>>(); }
template<> Datatype TypeMap<ValueInt<double>>() EL_NO_EXCEPT
{ return ValueIntType<double>(); }
template<> Datatype TypeMap<ValueInt<Complex<double>>>() EL_NO_EXCEPT
{ return ValueIntType<Complex<double>>(); }
#ifdef EL_HAVE_QD
template<> Datatype TypeMap<ValueInt<DoubleDouble>>() EL_NO_EXCEPT
{ return ValueIntType<DoubleDouble>(); }
template<> Datatype TypeMap<ValueInt<QuadDouble>>() EL_NO_EXCEPT
{ return ValueIntType<QuadDouble>(); }
#endif
#ifdef EL_HAVE_QUAD
template<> Datatype TypeMap<ValueInt<Quad>>() EL_NO_EXCEPT
{ return ValueIntType<Quad>(); }
template<> Datatype TypeMap<ValueInt<Complex<Quad>>>() EL_NO_EXCEPT
{ return ValueIntType<Complex<Quad>>(); }
#endif
#ifdef EL_HAVE_MPC
template<> Datatype TypeMap<ValueInt<BigInt>>() EL_NO_EXCEPT
{ return ValueIntType<BigInt>(); }
template<> Datatype TypeMap<ValueInt<BigFloat>>() EL_NO_EXCEPT
{ return ValueIntType<BigFloat>(); }
#endif

template<> Datatype TypeMap<Entry<Int>>() EL_NO_EXCEPT
{ return EntryType<Int>(); }
template<> Datatype TypeMap<Entry<float>>() EL_NO_EXCEPT
{ return EntryType<float>(); }
template<> Datatype TypeMap<Entry<Complex<float>>>() EL_NO_EXCEPT
{ return EntryType<Complex<float>>(); }
template<> Datatype TypeMap<Entry<double>>() EL_NO_EXCEPT
{ return EntryType<double>(); }
template<> Datatype TypeMap<Entry<Complex<double>>>() EL_NO_EXCEPT
{ return EntryType<Complex<double>>(); }
#ifdef EL_HAVE_QD
template<> Datatype TypeMap<Entry<DoubleDouble>>() EL_NO_EXCEPT
{ return EntryType<DoubleDouble>(); }
template<> Datatype TypeMap<Entry<QuadDouble>>() EL_NO_EXCEPT
{ return EntryType<QuadDouble>(); }
#endif
#ifdef EL_HAVE_QUAD
template<> Datatype TypeMap<Entry<Quad>>() EL_NO_EXCEPT
{ return EntryType<Quad>(); }
template<> Datatype TypeMap<Entry<Complex<Quad>>>() EL_NO_EXCEPT
{ return EntryType<Complex<Quad>>(); }
#endif
#ifdef EL_HAVE_MPC
template<> Datatype TypeMap<Entry<BigInt>>() EL_NO_EXCEPT
{ return EntryType<BigInt>(); }
template<> Datatype TypeMap<Entry<BigFloat>>() EL_NO_EXCEPT
{ return EntryType<BigFloat>(); }
#endif

template<typename T>
static void CreateValueIntType() EL_NO_EXCEPT
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

    Datatype tempType;
    CreateStruct( 2, blockLengths, displs, typeList, tempType );

    Datatype& type = ValueIntType<T>();
    MPI_Aint extent = sizeof(v);
    CreateResized( tempType, 0, extent, type );
}

#ifdef EL_HAVE_MPC
template<>
void CreateValueIntType<BigInt>() EL_NO_EXCEPT
{
    DEBUG_ONLY(CSE cse("CreateValueIntType [BigInt]"))

    Datatype typeList[2];
    typeList[0] = TypeMap<BigInt>();
    typeList[1] = TypeMap<Int>();
    
    int blockLengths[2];
    blockLengths[0] = 1;
    blockLengths[1] = 1; 

    BigInt alpha;
    const size_t packedSize = alpha.SerializedSize();

    MPI_Aint displs[2];
    displs[0] = 0;
    displs[1] = packedSize;

    Datatype tempType;
    CreateStruct( 2, blockLengths, displs, typeList, tempType );

    Datatype& type = ValueIntType<BigInt>();
    MPI_Aint extent = packedSize + sizeof(Int);
    CreateResized( tempType, 0, extent, type );
}

template<>
void CreateValueIntType<BigFloat>() EL_NO_EXCEPT
{
    DEBUG_ONLY(CSE cse("CreateValueIntType [BigFloat]"))

    Datatype typeList[2];
    typeList[0] = TypeMap<BigFloat>();
    typeList[1] = TypeMap<Int>();
    
    int blockLengths[2];
    blockLengths[0] = 1;
    blockLengths[1] = 1; 

    BigFloat alpha;
    const size_t packedSize = alpha.SerializedSize();

    MPI_Aint displs[2];
    displs[0] = 0;
    displs[1] = packedSize;

    Datatype tempType;
    CreateStruct( 2, blockLengths, displs, typeList, tempType );

    Datatype& type = ValueIntType<BigFloat>();
    MPI_Aint extent = packedSize + sizeof(Int);
    CreateResized( tempType, 0, extent, type );
}
#endif

template void CreateValueIntType<Int>() EL_NO_EXCEPT;
template void CreateValueIntType<float>() EL_NO_EXCEPT;
template void CreateValueIntType<double>() EL_NO_EXCEPT;
template void CreateValueIntType<Complex<float>>() EL_NO_EXCEPT;
template void CreateValueIntType<Complex<double>>() EL_NO_EXCEPT;
#ifdef EL_HAVE_QD
template void CreateValueIntType<DoubleDouble>() EL_NO_EXCEPT;
template void CreateValueIntType<QuadDouble>() EL_NO_EXCEPT;
#endif
#ifdef EL_HAVE_QUAD
template void CreateValueIntType<Quad>() EL_NO_EXCEPT;
template void CreateValueIntType<Complex<Quad>>() EL_NO_EXCEPT;
#endif

template<typename T>
static void CreateEntryType() EL_NO_EXCEPT
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

    Datatype tempType;
    CreateStruct( 3, blockLengths, displs, typeList, tempType );

    Datatype& type = EntryType<T>();
    MPI_Aint extent = sizeof(v);
    CreateResized( tempType, 0, extent, type );
}

#ifdef EL_HAVE_MPC
template<>
void CreateEntryType<BigInt>() EL_NO_EXCEPT
{
    DEBUG_ONLY(CSE cse("CreateEntryType [BigInt]"))

    Datatype typeList[3];
    typeList[0] = TypeMap<Int>();
    typeList[1] = TypeMap<Int>();
    typeList[2] = TypeMap<BigInt>();
    
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

    BigInt alpha;
    const auto packedSize = alpha.SerializedSize();

    Datatype& type = EntryType<BigInt>();
    MPI_Aint extent = 2*sizeof(Int) + packedSize;
    CreateResized( tempType, 0, extent, type );
}

template<>
void CreateEntryType<BigFloat>() EL_NO_EXCEPT
{
    DEBUG_ONLY(CSE cse("CreateEntryType [BigFloat]"))

    Datatype typeList[3];
    typeList[0] = TypeMap<Int>();
    typeList[1] = TypeMap<Int>();
    typeList[2] = TypeMap<BigFloat>();
    
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

    BigFloat alpha;
    const auto packedSize = alpha.SerializedSize();

    Datatype& type = EntryType<BigFloat>();
    MPI_Aint extent = 2*sizeof(Int) + packedSize;
    CreateResized( tempType, 0, extent, type );
}
#endif

template void CreateEntryType<Int>() EL_NO_EXCEPT;
template void CreateEntryType<float>() EL_NO_EXCEPT;
template void CreateEntryType<double>() EL_NO_EXCEPT;
template void CreateEntryType<Complex<float>>() EL_NO_EXCEPT;
template void CreateEntryType<Complex<double>>() EL_NO_EXCEPT;
#ifdef EL_HAVE_QD
template void CreateEntryType<DoubleDouble>() EL_NO_EXCEPT;
template void CreateEntryType<QuadDouble>() EL_NO_EXCEPT;
#endif
#ifdef EL_HAVE_QUAD
template void CreateEntryType<Quad>() EL_NO_EXCEPT;
template void CreateEntryType<Complex<Quad>>() EL_NO_EXCEPT;
#endif

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
    CreateValueIntType<BigFloat>();
    CreateEntryType<BigFloat>();
}

void DestroyBigIntFamily()
{
    Free( EntryType<BigInt>() );
    Free( ValueIntType<BigInt>() );
    Free( ::BigIntType );
}

void DestroyBigFloatFamily()
{
    Free( EntryType<BigFloat>() );
    Free( ValueIntType<BigFloat>() );
    Free( ::BigFloatType );
}
#endif

void CreateCustom() EL_NO_RELEASE_EXCEPT
{
    // Create the necessary types
    // ==========================
#ifdef EL_HAVE_QD
    // Create an MPI type for DoubleDouble
    // -----------------------------------
    CreateContiguous( 2, MPI_DOUBLE, ::DoubleDoubleType );

    // Create an MPI type for QuadDouble
    // ---------------------------------
    CreateContiguous( 4, MPI_DOUBLE, ::QuadDoubleType );

#endif
#ifdef EL_HAVE_QUAD
    // Create an MPI type for Quad
    // ---------------------------
    CreateContiguous( 2, MPI_DOUBLE, ::QuadType );

    // Create an MPI type for Complex<Quad>
    // ------------------------------------
    CreateContiguous( 4, MPI_DOUBLE, ::QuadComplexType );

#endif
    // NOTE: The BigFloat types are created by mpc::SetPrecision previously
    //       within El::Initialize

    // A value and an integer
    // ----------------------
    mpi::CreateValueIntType<Int>();
#ifdef EL_USE_64BIT_INTS
    mpi::CreateValueIntType<float>();
    mpi::CreateValueIntType<double>();
#else
    ::floatIntType = MPI_FLOAT_INT;
    ::doubleIntType = MPI_DOUBLE_INT;
#endif
    mpi::CreateValueIntType<Complex<float>>();
    mpi::CreateValueIntType<Complex<double>>();
#ifdef EL_HAVE_QD
    mpi::CreateValueIntType<DoubleDouble>();
    mpi::CreateValueIntType<QuadDouble>();
#endif
#ifdef EL_HAVE_QUAD
    mpi::CreateValueIntType<Quad>();
    mpi::CreateValueIntType<Complex<Quad>>();
#endif

    // A triplet of a value and a pair of integers
    // -------------------------------------------
    mpi::CreateEntryType<Int>();
    mpi::CreateEntryType<float>();
    mpi::CreateEntryType<double>();
    mpi::CreateEntryType<Complex<float>>();
    mpi::CreateEntryType<Complex<double>>();
#ifdef EL_HAVE_QD
    mpi::CreateEntryType<DoubleDouble>();
    mpi::CreateEntryType<QuadDouble>();
#endif
#ifdef EL_HAVE_QUAD
    mpi::CreateEntryType<Quad>();
    mpi::CreateEntryType<Complex<Quad>>();
#endif

    // Create the necessary MPI operations
    // ===================================
    // Functions for user-defined ops
    // ------------------------------
    Create
    ( (UserFunction*)UserIntReduce, false, ::userIntOp );
    Create
    ( (UserFunction*)UserIntReduceComm, true, ::userIntCommOp );
    Create
    ( (UserFunction*)UserFloatReduce, false, ::userFloatOp );
    Create
    ( (UserFunction*)UserFloatReduceComm, true, ::userFloatCommOp );
    Create
    ( (UserFunction*)UserDoubleReduce, false, ::userDoubleOp );
    Create
    ( (UserFunction*)UserDoubleReduceComm, true, ::userDoubleCommOp );
    Create
    ( (UserFunction*)UserComplexFloatReduce, false,
                   ::userComplexFloatOp );
    Create
    ( (UserFunction*)UserComplexFloatReduceComm, true,
                   ::userComplexFloatCommOp );
    Create
    ( (UserFunction*)UserComplexDoubleReduce, false,
                   ::userComplexDoubleOp );
    Create
    ( (UserFunction*)UserComplexDoubleReduceComm, true,
                   ::userComplexDoubleCommOp );
#ifdef EL_HAVE_QD
    Create
    ( (UserFunction*)UserDoubleDoubleReduce, false,
      ::userDoubleDoubleOp );
    Create
    ( (UserFunction*)UserDoubleDoubleReduceComm, true,
      ::userDoubleDoubleCommOp );
    Create
    ( (UserFunction*)UserQuadDoubleReduce, false,
      ::userQuadDoubleOp );
    Create
    ( (UserFunction*)UserQuadDoubleReduceComm, true,
      ::userQuadDoubleCommOp );
#endif
#ifdef EL_HAVE_QUAD
    Create
    ( (UserFunction*)UserQuadReduce, false, ::userQuadOp );
    Create
    ( (UserFunction*)UserQuadReduceComm, true, ::userQuadCommOp );
    Create
    ( (UserFunction*)UserComplexQuadReduce, false, ::userComplexQuadOp );
    Create
    ( (UserFunction*)UserComplexQuadReduceComm, true, ::userComplexQuadCommOp );
#endif
#ifdef EL_HAVE_MPC
    Create
    ( (UserFunction*)UserBigIntReduce, false, ::userBigIntOp );
    Create
    ( (UserFunction*)UserBigIntReduceComm, true, ::userBigIntCommOp );

    Create
    ( (UserFunction*)UserBigFloatReduce, false, ::userBigFloatOp );
    Create
    ( (UserFunction*)UserBigFloatReduceComm, true, ::userBigFloatCommOp );
    // TODO: Complex versions
#endif
   
    // Functions for scalar types
    // --------------------------
#ifdef EL_HAVE_QD
    Create( (UserFunction*)MaxDoubleDouble, true, ::maxDoubleDoubleOp );
    Create( (UserFunction*)MinDoubleDouble, true, ::minDoubleDoubleOp );
    Create( (UserFunction*)SumDoubleDouble, true, ::sumDoubleDoubleOp );

    Create( (UserFunction*)MaxQuadDouble, true, ::maxQuadDoubleOp );
    Create( (UserFunction*)MinQuadDouble, true, ::minQuadDoubleOp );
    Create( (UserFunction*)SumQuadDouble, true, ::sumQuadDoubleOp );
#endif
#ifdef EL_HAVE_QUAD
    Create( (UserFunction*)MaxQuad, true, ::maxQuadOp );
    Create( (UserFunction*)MinQuad, true, ::minQuadOp );
    Create( (UserFunction*)SumQuad, true, ::sumQuadOp );
    Create( (UserFunction*)SumQuadComplex, true, ::sumQuadComplexOp );
#endif
#ifdef EL_HAVE_MPC
    Create( (UserFunction*)MaxBigInt, true, ::maxBigIntOp );
    Create( (UserFunction*)MinBigInt, true, ::minBigIntOp );
    Create( (UserFunction*)SumBigInt, true, ::sumBigIntOp );

    Create( (UserFunction*)MaxBigFloat, true, ::maxBigFloatOp );
    Create( (UserFunction*)MinBigFloat, true, ::minBigFloatOp );
    Create( (UserFunction*)SumBigFloat, true, ::sumBigFloatOp );
    // TODO: Complex sum
#endif
    // Functions for the value and integer
    // -----------------------------------
    Create( (UserFunction*)MaxLocFunc<Int>,    true, ::maxLocIntOp );
    Create( (UserFunction*)MinLocFunc<Int>,    true, ::minLocIntOp );
#ifdef EL_USE_64BIT_INTS
    Create( (UserFunction*)MaxLocFunc<float>,  true, ::maxLocFloatOp  );
    Create( (UserFunction*)MinLocFunc<float>,  true, ::minLocFloatOp  );
    Create( (UserFunction*)MaxLocFunc<double>, true, ::maxLocDoubleOp );
    Create( (UserFunction*)MinLocFunc<double>, true, ::minLocDoubleOp );
#else
    ::maxLocFloatOp = MAXLOC; 
    ::minLocFloatOp = MINLOC;
    ::maxLocDoubleOp = MAXLOC;
    ::minLocDoubleOp = MINLOC;
#endif
#ifdef EL_HAVE_QD
    Create
    ( (UserFunction*)MaxLocFunc<DoubleDouble>, true, ::maxLocDoubleDoubleOp );
    Create
    ( (UserFunction*)MinLocFunc<DoubleDouble>, true, ::minLocDoubleDoubleOp );
    Create
    ( (UserFunction*)MaxLocFunc<QuadDouble>, true, ::maxLocQuadDoubleOp );
    Create
    ( (UserFunction*)MinLocFunc<QuadDouble>, true, ::minLocQuadDoubleOp );
#endif
#ifdef EL_HAVE_QUAD
    Create( (UserFunction*)MaxLocFunc<Quad>, true, ::maxLocQuadOp );
    Create( (UserFunction*)MinLocFunc<Quad>, true, ::minLocQuadOp );
#endif
#ifdef EL_HAVE_MPC
    Create( (UserFunction*)MaxLocFunc<BigInt>, true, ::maxLocBigIntOp );
    Create( (UserFunction*)MinLocFunc<BigInt>, true, ::minLocBigIntOp );

    Create( (UserFunction*)MaxLocFunc<BigFloat>, true, ::maxLocBigFloatOp );
    Create( (UserFunction*)MinLocFunc<BigFloat>, true, ::minLocBigFloatOp );
#endif

    // Functions for the triplet of a value and a pair of integers
    // -----------------------------------------------------------
    Create( (UserFunction*)MaxLocPairFunc<Int>,    true, ::maxLocPairIntOp    );
    Create( (UserFunction*)MinLocPairFunc<Int>,    true, ::minLocPairIntOp    );
    Create( (UserFunction*)MaxLocPairFunc<float>,  true, ::maxLocPairFloatOp  );
    Create( (UserFunction*)MinLocPairFunc<float>,  true, ::minLocPairFloatOp  );
    Create( (UserFunction*)MaxLocPairFunc<double>, true, ::maxLocPairDoubleOp );
    Create( (UserFunction*)MinLocPairFunc<double>, true, ::minLocPairDoubleOp );
#ifdef EL_HAVE_QD
    Create
    ( (UserFunction*)MaxLocPairFunc<DoubleDouble>, true,
      ::maxLocPairDoubleDoubleOp );
    Create
    ( (UserFunction*)MinLocPairFunc<DoubleDouble>, true,
      ::minLocPairDoubleDoubleOp );
    Create
    ( (UserFunction*)MaxLocPairFunc<QuadDouble>, true,
      ::maxLocPairQuadDoubleOp );
    Create
    ( (UserFunction*)MinLocPairFunc<QuadDouble>, true,
      ::minLocPairQuadDoubleOp );
#endif
#ifdef EL_HAVE_QUAD
    Create( (UserFunction*)MaxLocPairFunc<Quad>,   true, ::maxLocPairQuadOp );
    Create( (UserFunction*)MinLocPairFunc<Quad>,   true, ::minLocPairQuadOp );
#endif
#ifdef EL_HAVE_MPC
    Create
    ( (UserFunction*)MaxLocPairFunc<BigInt>, true, ::maxLocPairBigIntOp );
    Create
    ( (UserFunction*)MinLocPairFunc<BigInt>, true, ::minLocPairBigIntOp );

    Create
    ( (UserFunction*)MaxLocPairFunc<BigFloat>, true, ::maxLocPairBigFloatOp );
    Create
    ( (UserFunction*)MinLocPairFunc<BigFloat>, true, ::minLocPairBigFloatOp );
#endif
}

void DestroyCustom() EL_NO_RELEASE_EXCEPT
{
    // Destroy the created types
    // =========================
    Free( EntryType<Int>() );
    Free( EntryType<float>() );
    Free( EntryType<double>() );
    Free( EntryType<Complex<float>>() );
    Free( EntryType<Complex<double>>() );
#ifdef EL_HAVE_QD
    Free( EntryType<DoubleDouble>() );
    Free( EntryType<QuadDouble>() );
#endif
#ifdef EL_HAVE_QUAD
    Free( EntryType<Quad>() );
    Free( EntryType<Complex<Quad>>() );
#endif

    Free( ValueIntType<Int>() );
#ifdef EL_USE_64BIT_INTS
    Free( ValueIntType<float>() );
    Free( ValueIntType<double>() );
#endif
    Free( ValueIntType<Complex<float>>() );
    Free( ValueIntType<Complex<double>>() );
#ifdef EL_HAVE_QD
    Free( ValueIntType<DoubleDouble>() );
    Free( ValueIntType<QuadDouble>() );
#endif
#ifdef EL_HAVE_QUAD
    Free( ValueIntType<Quad>() );
    Free( ValueIntType<Complex<Quad>>() );
#endif

#ifdef EL_HAVE_QD
    Free( ::DoubleDoubleType );
    Free( ::QuadDoubleType );
#endif
#ifdef EL_HAVE_QUAD
    Free( ::QuadType );
    Free( ::QuadComplexType );
#endif
#ifdef EL_HAVE_MPC
    mpi::DestroyBigIntFamily();
    mpi::DestroyBigFloatFamily();
#endif

    // Destroy the created operations
    // ==============================

    // User-defined operations
    // -----------------------
    Free( ::userIntOp );
    Free( ::userIntCommOp );
    Free( ::userFloatOp );
    Free( ::userFloatCommOp );
    Free( ::userDoubleOp );
    Free( ::userDoubleCommOp );
    Free( ::userComplexFloatOp );
    Free( ::userComplexFloatCommOp );
    Free( ::userComplexDoubleOp );
    Free( ::userComplexDoubleCommOp );
#ifdef EL_HAVE_QD
    Free( ::userDoubleDoubleOp );
    Free( ::userDoubleDoubleCommOp );
    Free( ::userQuadDoubleOp );
    Free( ::userQuadDoubleCommOp );
#endif
#ifdef EL_HAVE_QUAD
    Free( ::userQuadOp );
    Free( ::userQuadCommOp );
    Free( ::userComplexQuadOp );
    Free( ::userComplexQuadCommOp );
#endif
#ifdef EL_HAVE_MPC
    Free( ::userBigIntOp );
    Free( ::userBigIntCommOp );
    Free( ::userBigFloatOp );
    Free( ::userBigFloatCommOp );
#endif

#ifdef EL_HAVE_QD
    Free( ::maxDoubleDoubleOp );
    Free( ::minDoubleDoubleOp );
    Free( ::sumDoubleDoubleOp );
    Free( ::maxQuadDoubleOp );
    Free( ::minQuadDoubleOp );
    Free( ::sumQuadDoubleOp );
#endif
#ifdef EL_HAVE_QUAD
    Free( ::maxQuadOp );
    Free( ::minQuadOp );
    Free( ::sumQuadOp );
    Free( ::sumQuadComplexOp );
#endif
#ifdef EL_HAVE_MPC
    Free( ::maxBigIntOp );
    Free( ::minBigIntOp );
    Free( ::sumBigIntOp );
    Free( ::maxBigFloatOp );
    Free( ::minBigFloatOp );
    Free( ::sumBigFloatOp );
#endif

    Free( ::maxLocIntOp );
    Free( ::minLocIntOp );
#ifdef EL_USE_64BIT_INTS
    Free( ::maxLocFloatOp );
    Free( ::minLocFloatOp );
    Free( ::maxLocDoubleOp );
    Free( ::minLocDoubleOp );
#endif
#ifdef EL_HAVE_QD
    Free( ::maxLocDoubleDoubleOp );
    Free( ::minLocDoubleDoubleOp );
    Free( ::maxLocQuadDoubleOp );
    Free( ::minLocQuadDoubleOp );
#endif
#ifdef EL_HAVE_QUAD
    Free( ::maxLocQuadOp );
    Free( ::minLocQuadOp );
#endif
#ifdef EL_HAVE_MPC
    Free( ::maxLocBigIntOp );
    Free( ::minLocBigIntOp );
    Free( ::maxLocBigFloatOp );
    Free( ::minLocBigFloatOp );
#endif

    Free( ::maxLocPairIntOp );
    Free( ::minLocPairIntOp );
    Free( ::maxLocPairFloatOp );
    Free( ::minLocPairFloatOp );
    Free( ::maxLocPairDoubleOp );
    Free( ::minLocPairDoubleOp );
#ifdef EL_HAVE_QD
    Free( ::maxLocPairDoubleDoubleOp );
    Free( ::minLocPairDoubleDoubleOp );
    Free( ::maxLocPairQuadDoubleOp );
    Free( ::minLocPairQuadDoubleOp );
#endif
#ifdef EL_HAVE_QUAD
    Free( ::maxLocPairQuadOp );
    Free( ::minLocPairQuadOp );
#endif
#ifdef EL_HAVE_MPC
    Free( ::maxLocPairBigIntOp );
    Free( ::minLocPairBigIntOp );
    Free( ::maxLocPairBigFloatOp );
    Free( ::minLocPairBigFloatOp );
#endif
}

template<> Op UserOp<Int>() EL_NO_EXCEPT { return ::userIntOp; }
template<> Op UserCommOp<Int>() EL_NO_EXCEPT { return ::userIntCommOp; }
template<> Op UserOp<float>() EL_NO_EXCEPT { return ::userFloatOp; }
template<> Op UserCommOp<float>() EL_NO_EXCEPT { return ::userFloatCommOp; }
template<> Op UserOp<double>() EL_NO_EXCEPT { return ::userDoubleOp; }
template<> Op UserCommOp<double>() EL_NO_EXCEPT { return ::userDoubleCommOp; }
template<> Op UserOp<Complex<float>>() EL_NO_EXCEPT
{ return ::userComplexFloatOp; }
template<> Op UserCommOp<Complex<float>>() EL_NO_EXCEPT
{ return ::userComplexFloatCommOp; }
template<> Op UserOp<Complex<double>>() EL_NO_EXCEPT
{ return ::userComplexDoubleOp; }
template<> Op UserCommOp<Complex<double>>() EL_NO_EXCEPT
{ return ::userComplexDoubleCommOp; }
#ifdef EL_HAVE_QD
template<> Op UserOp<DoubleDouble>() EL_NO_EXCEPT
{ return ::userDoubleDoubleOp; }
template<> Op UserCommOp<DoubleDouble>() EL_NO_EXCEPT
{ return ::userDoubleDoubleCommOp; }
template<> Op UserOp<QuadDouble>() EL_NO_EXCEPT
{ return ::userQuadDoubleOp; }
template<> Op UserCommOp<QuadDouble>() EL_NO_EXCEPT
{ return ::userQuadDoubleCommOp; }
#endif
#ifdef EL_HAVE_QUAD
template<> Op UserOp<Quad>() EL_NO_EXCEPT
{ return ::userQuadOp; }
template<> Op UserCommOp<Quad>() EL_NO_EXCEPT
{ return ::userQuadCommOp; }
template<> Op UserOp<Complex<Quad>>() EL_NO_EXCEPT
{ return ::userComplexQuadOp; }
template<> Op UserCommOp<Complex<Quad>>() EL_NO_EXCEPT
{ return ::userComplexQuadCommOp; }
#endif
#ifdef EL_HAVE_MPC
template<> Op UserOp<BigInt>() EL_NO_EXCEPT
{ return ::userBigIntOp; }
template<> Op UserCommOp<BigInt>() EL_NO_EXCEPT
{ return ::userBigIntCommOp; }

template<> Op UserOp<BigFloat>() EL_NO_EXCEPT
{ return ::userBigFloatOp; }
template<> Op UserCommOp<BigFloat>() EL_NO_EXCEPT
{ return ::userBigFloatCommOp; }
#endif

#ifdef EL_HAVE_QD
template<> Op MaxOp<DoubleDouble>() EL_NO_EXCEPT { return ::maxDoubleDoubleOp; }
template<> Op MinOp<DoubleDouble>() EL_NO_EXCEPT { return ::minDoubleDoubleOp; }
template<> Op SumOp<DoubleDouble>() EL_NO_EXCEPT { return ::sumDoubleDoubleOp; }

template<> Op MaxOp<QuadDouble>() EL_NO_EXCEPT { return ::maxQuadDoubleOp; }
template<> Op MinOp<QuadDouble>() EL_NO_EXCEPT { return ::minQuadDoubleOp; }
template<> Op SumOp<QuadDouble>() EL_NO_EXCEPT { return ::sumQuadDoubleOp; }
#endif
#ifdef EL_HAVE_QUAD
template<> Op MaxOp<Quad>() EL_NO_EXCEPT { return ::maxQuadOp; }
template<> Op MinOp<Quad>() EL_NO_EXCEPT { return ::minQuadOp; }

template<> Op SumOp<Quad>() EL_NO_EXCEPT { return ::sumQuadOp; }
template<> Op SumOp<Complex<Quad>>() EL_NO_EXCEPT { return ::sumQuadComplexOp; }
#endif
#ifdef EL_HAVE_MPC
template<> Op MaxOp<BigInt>() EL_NO_EXCEPT { return ::maxBigIntOp; }
template<> Op MinOp<BigInt>() EL_NO_EXCEPT { return ::minBigIntOp; }
template<> Op SumOp<BigInt>() EL_NO_EXCEPT { return ::sumBigIntOp; }

template<> Op MaxOp<BigFloat>() EL_NO_EXCEPT { return ::maxBigFloatOp; }
template<> Op MinOp<BigFloat>() EL_NO_EXCEPT { return ::minBigFloatOp; }
template<> Op SumOp<BigFloat>() EL_NO_EXCEPT { return ::sumBigFloatOp; }
#endif

template<> Op MaxLocOp<Int>() EL_NO_EXCEPT { return ::maxLocIntOp; }
template<> Op MinLocOp<Int>() EL_NO_EXCEPT { return ::minLocIntOp; }
template<> Op MaxLocOp<float>() EL_NO_EXCEPT { return ::maxLocFloatOp; }
template<> Op MinLocOp<float>() EL_NO_EXCEPT { return ::minLocFloatOp; }
template<> Op MaxLocOp<double>() EL_NO_EXCEPT { return ::maxLocDoubleOp; }
template<> Op MinLocOp<double>() EL_NO_EXCEPT { return ::minLocDoubleOp; }
#ifdef EL_HAVE_QD
template<> Op MaxLocOp<DoubleDouble>() EL_NO_EXCEPT
{ return ::maxLocDoubleDoubleOp; }
template<> Op MinLocOp<DoubleDouble>() EL_NO_EXCEPT
{ return ::minLocDoubleDoubleOp; }

template<> Op MaxLocOp<QuadDouble>() EL_NO_EXCEPT
{ return ::maxLocQuadDoubleOp; }
template<> Op MinLocOp<QuadDouble>() EL_NO_EXCEPT
{ return ::minLocQuadDoubleOp; }
#endif
#ifdef EL_HAVE_QUAD
template<> Op MaxLocOp<Quad>() EL_NO_EXCEPT { return ::maxLocQuadOp; }
template<> Op MinLocOp<Quad>() EL_NO_EXCEPT { return ::minLocQuadOp; }
#endif
#ifdef EL_HAVE_MPC
template<> Op MaxLocOp<BigInt>() EL_NO_EXCEPT { return ::maxLocBigIntOp; }
template<> Op MinLocOp<BigInt>() EL_NO_EXCEPT { return ::minLocBigIntOp; }

template<> Op MaxLocOp<BigFloat>() EL_NO_EXCEPT { return ::maxLocBigFloatOp; }
template<> Op MinLocOp<BigFloat>() EL_NO_EXCEPT { return ::minLocBigFloatOp; }
#endif

template<> Op MaxLocPairOp<Int>() EL_NO_EXCEPT
{ return ::maxLocPairIntOp; }
template<> Op MinLocPairOp<Int>() EL_NO_EXCEPT
{ return ::minLocPairIntOp; }
template<> Op MaxLocPairOp<float>() EL_NO_EXCEPT
{ return ::maxLocPairFloatOp; }
template<> Op MinLocPairOp<float>() EL_NO_EXCEPT
{ return ::minLocPairFloatOp; }
template<> Op MaxLocPairOp<double>() EL_NO_EXCEPT
{ return ::maxLocPairDoubleOp; }
template<> Op MinLocPairOp<double>() EL_NO_EXCEPT
{ return ::minLocPairDoubleOp; }
#ifdef EL_HAVE_QD
template<> Op MaxLocPairOp<DoubleDouble>() EL_NO_EXCEPT
{ return ::maxLocPairDoubleDoubleOp; }
template<> Op MinLocPairOp<DoubleDouble>() EL_NO_EXCEPT
{ return ::minLocPairDoubleDoubleOp; }

template<> Op MaxLocPairOp<QuadDouble>() EL_NO_EXCEPT
{ return ::maxLocPairQuadDoubleOp; }
template<> Op MinLocPairOp<QuadDouble>() EL_NO_EXCEPT
{ return ::minLocPairQuadDoubleOp; }
#endif
#ifdef EL_HAVE_QUAD
template<> Op MaxLocPairOp<Quad>() EL_NO_EXCEPT
{ return ::maxLocPairQuadOp; }
template<> Op MinLocPairOp<Quad>() EL_NO_EXCEPT
{ return ::minLocPairQuadOp; }
#endif
#ifdef EL_HAVE_MPC
template<> Op MaxLocPairOp<BigInt>() EL_NO_EXCEPT
{ return ::maxLocPairBigIntOp; }
template<> Op MinLocPairOp<BigInt>() EL_NO_EXCEPT
{ return ::minLocPairBigIntOp; }

template<> Op MaxLocPairOp<BigFloat>() EL_NO_EXCEPT
{ return ::maxLocPairBigFloatOp; }
template<> Op MinLocPairOp<BigFloat>() EL_NO_EXCEPT
{ return ::minLocPairBigFloatOp; }
#endif

} // namespace mpi
} // namespace El
