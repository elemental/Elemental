/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"
#include "El-C.h"
using namespace El;

#define RCM_s(AHandle) reinterpret_cast<Matrix<float>*>(AHandle)
#define RCM_d(AHandle) reinterpret_cast<Matrix<double>*>(AHandle)
#define RCM_c(AHandle) reinterpret_cast<Matrix<Complex<float>>*>(AHandle)
#define RCM_z(AHandle) reinterpret_cast<Matrix<Complex<double>>*>(AHandle)

#define RCM_s_const(AHandle) reinterpret_cast<const Matrix<float>*>(AHandle)
#define RCM_d_const(AHandle) reinterpret_cast<const Matrix<double>*>(AHandle)
#define RCM_c_const(AHandle) \
    reinterpret_cast<const Matrix<Complex<float>>*>(AHandle)
#define RCM_z_const(AHandle) \
    reinterpret_cast<const Matrix<Complex<double>>*>(AHandle)

#define RCB_c(buffer) reinterpret_cast<Complex<float>*>(buffer)
#define RCB_z(buffer) reinterpret_cast<Complex<double>*>(buffer)

#define RCB_c_const(buffer) reinterpret_cast<const Complex<float>*>(buffer)
#define RCB_z_const(buffer) reinterpret_cast<const Complex<double>*>(buffer)

#define CATCH catch( std::exception& e ) { ReportException(e); }

extern "C" {

// Matrix::Matrix()
// ----------------

ElMatrix_s* ElMatrixCreate_s()
{
    try { return reinterpret_cast<ElMatrix_s*>( new Matrix<float>() ); }
    CATCH
}

ElMatrix_d* ElMatrixCreate_d()
{
    try { return reinterpret_cast<ElMatrix_d*>( new Matrix<double>() ); }
    CATCH
}

ElMatrix_c* ElMatrixCreate_c()
{
    try 
    { return reinterpret_cast<ElMatrix_c*>( new Matrix<Complex<float>>() ); }
    CATCH
}

ElMatrix_z* ElMatrixCreate_z()
{
    try 
    { return reinterpret_cast<ElMatrix_z*>( new Matrix<Complex<double>>() ); }
    CATCH
}

// Matrix::~Matrix()
// -----------------

void ElMatrixDestroy_s( const ElMatrix_s* AHandle )
{ delete RCM_s_const(AHandle); }

void ElMatrixDestroy_d( const ElMatrix_d* AHandle )
{ delete RCM_d_const(AHandle); }

void ElMatrixDestroy_c( const ElMatrix_c* AHandle )
{ delete RCM_c_const(AHandle); }

void ElMatrixDestroy_z( const ElMatrix_z* AHandle )
{ delete RCM_z_const(AHandle); }

// void Matrix::Empty()
// --------------------

void ElMatrixEmpty_s( ElMatrix_s* AHandle )
{ RCM_s(AHandle)->Empty(); }

void ElMatrixEmpty_d( ElMatrix_d* AHandle )
{ RCM_d(AHandle)->Empty(); }

void ElMatrixEmpty_c( ElMatrix_c* AHandle )
{ RCM_c(AHandle)->Empty(); }

void ElMatrixEmpty_z( ElMatrix_z* AHandle )
{ RCM_z(AHandle)->Empty(); }

// void Matrix::Resize( Int height, Int width )
// --------------------------------------------

void ElMatrixResize_s( ElMatrix_s* AHandle, ElInt height, ElInt width )
{
    try { RCM_s(AHandle)->Resize(height,width); }
    CATCH
}

void ElMatrixResize_d( ElMatrix_d* AHandle, ElInt height, ElInt width )
{
    try { RCM_d(AHandle)->Resize(height,width); }
    CATCH
}

void ElMatrixResize_c( ElMatrix_c* AHandle, ElInt height, ElInt width )
{
    try { RCM_c(AHandle)->Resize(height,width); }
    CATCH
}

void ElMatrixResize_z( ElMatrix_z* AHandle, ElInt height, ElInt width )
{
    try { RCM_z(AHandle)->Resize(height,width); }
    CATCH
}

// void Matrix::Resize( Int height, Int width, Int ldim )
// ------------------------------------------------------

void ElMatrixResizeWithLDim_s
( ElMatrix_s* AHandle, ElInt height, ElInt width, ElInt ldim )
{
    try { RCM_s(AHandle)->Resize(height,width,ldim); }
    CATCH
}

void ElMatrixResizeWithLDim_d
( ElMatrix_d* AHandle, ElInt height, ElInt width, ElInt ldim )
{
    try { RCM_d(AHandle)->Resize(height,width,ldim); }
    CATCH
}

void ElMatrixResizeWithLDim_c
( ElMatrix_c* AHandle, ElInt height, ElInt width, ElInt ldim )
{
    try { RCM_c(AHandle)->Resize(height,width,ldim); }
    CATCH
}

void ElMatrixResizeWithLDim_z
( ElMatrix_z* AHandle, ElInt height, ElInt width, ElInt ldim )
{
    try { RCM_z(AHandle)->Resize(height,width,ldim); }
    CATCH
}

// void Matrix::Attach( Int height, Int width, T* buffer, Int ldim )
// -----------------------------------------------------------------

void ElMatrixAttach_s
( ElMatrix_s* AHandle, ElInt height, ElInt width, float* buffer, ElInt ldim )
{
    try { RCM_s(AHandle)->Attach(height,width,buffer,ldim); }
    CATCH
}

void ElMatrixAttach_d
( ElMatrix_d* AHandle, ElInt height, ElInt width, double* buffer, ElInt ldim )
{
    try { RCM_d(AHandle)->Attach(height,width,buffer,ldim); }
    CATCH
}

void ElMatrixAttach_c
( ElMatrix_c* AHandle, ElInt height, ElInt width, void* buffer, ElInt ldim )
{
    try { RCM_c(AHandle)->Attach(height,width,RCB_c(buffer),ldim); }
    CATCH
}

void ElMatrixAttach_z
( ElMatrix_z* AHandle, ElInt height, ElInt width, void* buffer, ElInt ldim )
{
    try { RCM_z(AHandle)->Attach(height,width,RCB_z(buffer),ldim); }
    CATCH
}

// void Matrix::LockedAttach( Int height, Int width, const T* buffer, Int ldim )
// -----------------------------------------------------------------------------

void ElMatrixLockedAttach_s
( ElMatrix_s* AHandle, 
  ElInt height, ElInt width, const float* buffer, ElInt ldim )
{
    try { RCM_s(AHandle)->LockedAttach(height,width,buffer,ldim); }
    CATCH
}

void ElMatrixLockedAttach_d
( ElMatrix_d* AHandle, 
  ElInt height, ElInt width, const double* buffer, ElInt ldim )
{
    try { RCM_d(AHandle)->LockedAttach(height,width,buffer,ldim); }
    CATCH
}

void ElMatrixLockedAttach_c
( ElMatrix_c* AHandle, 
  ElInt height, ElInt width, const void* buffer, ElInt ldim )
{
    try { RCM_c(AHandle)->LockedAttach(height,width,RCB_c_const(buffer),ldim); }
    CATCH
}

void ElMatrixLockedAttach_z
( ElMatrix_z* AHandle, 
  ElInt height, ElInt width, const void* buffer, ElInt ldim )
{
    try { RCM_z(AHandle)->LockedAttach(height,width,RCB_z_const(buffer),ldim); }
    CATCH
}

// void Matrix::Control( Int height, Int width, T* buffer, Int ldim )
// ------------------------------------------------------------------

void ElMatrixControl_s
( ElMatrix_s* AHandle, ElInt height, ElInt width, float* buffer, ElInt ldim )
{
    try { RCM_s(AHandle)->Control(height,width,buffer,ldim); }
    CATCH
}

void ElMatrixControl_d
( ElMatrix_d* AHandle, ElInt height, ElInt width, double* buffer, ElInt ldim )
{
    try { RCM_d(AHandle)->Control(height,width,buffer,ldim); }
    CATCH
}

void ElMatrixControl_c
( ElMatrix_c* AHandle, ElInt height, ElInt width, void* buffer, ElInt ldim )
{
    try { RCM_c(AHandle)->Control(height,width,RCB_c(buffer),ldim); }
    CATCH
}

void ElMatrixControl_z
( ElMatrix_z* AHandle, ElInt height, ElInt width, void* buffer, ElInt ldim )
{
    try { RCM_z(AHandle)->Control(height,width,RCB_z(buffer),ldim); }
    CATCH
}

// B = A
// -----

void ElMatrixCopy_s( const ElMatrix_s* AHandle, ElMatrix_s* BHandle )
{
    try { *RCM_s(BHandle) = *RCM_s_const(AHandle); }
    CATCH
}

void ElMatrixCopy_d( const ElMatrix_d* AHandle, ElMatrix_d* BHandle )
{
    try { *RCM_d(BHandle) = *RCM_d_const(AHandle); }
    CATCH
}

void ElMatrixCopy_c( const ElMatrix_c* AHandle, ElMatrix_c* BHandle )
{
    try { *RCM_c(BHandle) = *RCM_c_const(AHandle); }
    CATCH
}

void ElMatrixCopy_z( const ElMatrix_z* AHandle, ElMatrix_z* BHandle )
{
    try { *RCM_z(BHandle) = *RCM_z_const(AHandle); }
    CATCH
}

// Int Matrix::Height() const
// --------------------------

ElInt ElMatrixHeight_s( const ElMatrix_s* AHandle )
{ return RCM_s_const(AHandle)->Height(); }

ElInt ElMatrixHeight_d( const ElMatrix_d* AHandle )
{ return RCM_d_const(AHandle)->Height(); }

ElInt ElMatrixHeight_c( const ElMatrix_c* AHandle )
{ return RCM_c_const(AHandle)->Height(); }

ElInt ElMatrixHeight_z( const ElMatrix_z* AHandle )
{ return RCM_z_const(AHandle)->Height(); }

// Int Matrix::Width() const
// -------------------------

ElInt ElMatrixWidth_s( const ElMatrix_s* AHandle )
{ return RCM_s_const(AHandle)->Width(); }

ElInt ElMatrixWidth_d( const ElMatrix_d* AHandle )
{ return RCM_d_const(AHandle)->Width(); }

ElInt ElMatrixWidth_c( const ElMatrix_c* AHandle )
{ return RCM_c_const(AHandle)->Width(); }

ElInt ElMatrixWidth_z( const ElMatrix_z* AHandle )
{ return RCM_z_const(AHandle)->Width(); }

// Int Matrix::LDim() const
// -------------------------

ElInt ElMatrixLDim_s( const ElMatrix_s* AHandle )
{ return RCM_s_const(AHandle)->LDim(); }

ElInt ElMatrixLDim_d( const ElMatrix_d* AHandle )
{ return RCM_d_const(AHandle)->LDim(); }

ElInt ElMatrixLDim_c( const ElMatrix_c* AHandle )
{ return RCM_c_const(AHandle)->LDim(); }

ElInt ElMatrixLDim_z( const ElMatrix_z* AHandle )
{ return RCM_z_const(AHandle)->LDim(); }

// Int Matrix::MemorySize() const
// ------------------------------

ElInt ElMatrixMemorySize_s( const ElMatrix_s* AHandle )
{ return RCM_s_const(AHandle)->MemorySize(); }

ElInt ElMatrixMemorySize_d( const ElMatrix_d* AHandle )
{ return RCM_d_const(AHandle)->MemorySize(); }

ElInt ElMatrixMemorySize_c( const ElMatrix_c* AHandle )
{ return RCM_c_const(AHandle)->MemorySize(); }

ElInt ElMatrixMemorySize_z( const ElMatrix_z* AHandle )
{ return RCM_z_const(AHandle)->MemorySize(); }

// Int Matrix::DiagonalLength( Int offset ) const
// ----------------------------------------------

ElInt ElMatrixDiagonalLength_s( const ElMatrix_s* AHandle, ElInt offset )
{ return RCM_s_const(AHandle)->DiagonalLength(offset); }

ElInt ElMatrixDiagonalLength_d( const ElMatrix_d* AHandle, ElInt offset )
{ return RCM_d_const(AHandle)->DiagonalLength(offset); }

ElInt ElMatrixDiagonalLength_c( const ElMatrix_c* AHandle, ElInt offset )
{ return RCM_c_const(AHandle)->DiagonalLength(offset); }

ElInt ElMatrixDiagonalLength_z( const ElMatrix_z* AHandle, ElInt offset )
{ return RCM_z_const(AHandle)->DiagonalLength(offset); }

// T* Matrix::Buffer()
// -------------------

float* ElMatrixBuffer_s( ElMatrix_s* AHandle )
{ return RCM_s(AHandle)->Buffer(); }

double* ElMatrixBuffer_d( ElMatrix_d* AHandle )
{ return RCM_d(AHandle)->Buffer(); }

void* ElMatrixBuffer_c( ElMatrix_c* AHandle )
{ return RCM_c(AHandle)->Buffer(); }

void* ElMatrixBuffer_z( ElMatrix_z* AHandle )
{ return RCM_z(AHandle)->Buffer(); }

// const T* Matrix::LockedBuffer() const
// -------------------------------------

const float* ElMatrixLockedBuffer_s( const ElMatrix_s* AHandle )
{ return RCM_s_const(AHandle)->LockedBuffer(); }

const double* ElMatrixLockedBuffer_d( const ElMatrix_d* AHandle )
{ return RCM_d_const(AHandle)->LockedBuffer(); }

const void* ElMatrixLockedBuffer_c( const ElMatrix_c* AHandle )
{ return RCM_c_const(AHandle)->LockedBuffer(); }

const void* ElMatrixLockedBuffer_z( const ElMatrix_z* AHandle )
{ return RCM_z_const(AHandle)->LockedBuffer(); }

// bool Matrix::Viewing() const
// ----------------------------

bool ElMatrixViewing_s( const ElMatrix_s* AHandle )
{ return RCM_s_const(AHandle)->Viewing(); }

bool ElMatrixViewing_d( const ElMatrix_d* AHandle )
{ return RCM_d_const(AHandle)->Viewing(); }

bool ElMatrixViewing_c( const ElMatrix_c* AHandle )
{ return RCM_c_const(AHandle)->Viewing(); }

bool ElMatrixViewing_z( const ElMatrix_z* AHandle )
{ return RCM_z_const(AHandle)->Viewing(); }

// bool Matrix::FixedSize() const
// ------------------------------

bool ElMatrixFixedSize_s( const ElMatrix_s* AHandle )
{ return RCM_s_const(AHandle)->FixedSize(); }

bool ElMatrixFixedSize_d( const ElMatrix_d* AHandle )
{ return RCM_d_const(AHandle)->FixedSize(); }

bool ElMatrixFixedSize_c( const ElMatrix_c* AHandle )
{ return RCM_c_const(AHandle)->FixedSize(); }

bool ElMatrixFixedSize_z( const ElMatrix_z* AHandle )
{ return RCM_z_const(AHandle)->FixedSize(); }

// bool Matrix::Locked() const
// ---------------------------

bool ElMatrixLocked_s( const ElMatrix_s* AHandle )
{ return RCM_s_const(AHandle)->Locked(); }

bool ElMatrixLocked_d( const ElMatrix_d* AHandle )
{ return RCM_d_const(AHandle)->Locked(); }

bool ElMatrixLocked_c( const ElMatrix_c* AHandle )
{ return RCM_c_const(AHandle)->Locked(); }

bool ElMatrixLocked_z( const ElMatrix_z* AHandle )
{ return RCM_z_const(AHandle)->Locked(); }

// T Matrix::Get( Int i, Int j ) const
// -----------------------------------

float ElMatrixGet_s( const ElMatrix_s* AHandle, ElInt i, ElInt j )
{
    float alpha;
    try { alpha = RCM_s_const(AHandle)->Get(i,j); }
    CATCH
    return alpha;
}

double ElMatrixGet_d( const ElMatrix_d* AHandle, ElInt i, ElInt j )
{
    float alpha;
    try { alpha = RCM_s_const(AHandle)->Get(i,j); }
    CATCH
    return alpha;
}

void ElMatrixGet_c( const ElMatrix_c* AHandle, ElInt i, ElInt j, void* alpha )
{
    try { *RCB_c(alpha) = RCM_c_const(AHandle)->Get(i,j); }
    CATCH
}

void ElMatrixGet_z( const ElMatrix_z* AHandle, ElInt i, ElInt j, void* alpha )
{
    try { *RCB_z(alpha) = RCM_z_const(AHandle)->Get(i,j); }
    CATCH
}

// Base<T> Matrix::GetRealPart( Int i, Int j ) const
// -------------------------------------------------

float ElMatrixGetRealPart_c( const ElMatrix_c* AHandle, ElInt i, ElInt j )
{
    float alpha;
    try { alpha = RCM_s_const(AHandle)->GetRealPart(i,j); }
    CATCH
    return alpha;
}

double ElMatrixGetRealPart_z( const ElMatrix_z* AHandle, ElInt i, ElInt j )
{
    double alpha;
    try { alpha = RCM_d_const(AHandle)->GetRealPart(i,j); }
    CATCH
    return alpha;
}

// Base<T> Matrix::GetImagPart( Int i, Int j ) const
// -------------------------------------------------

float ElMatrixGetImagPart_c( const ElMatrix_c* AHandle, ElInt i, ElInt j )
{
    float alpha;
    try { alpha = RCM_s_const(AHandle)->GetImagPart(i,j); }
    CATCH
    return alpha;
}

double ElMatrixGetImagPart_z( const ElMatrix_z* AHandle, ElInt i, ElInt j )
{
    double alpha;
    try { alpha = RCM_d_const(AHandle)->GetImagPart(i,j); }
    CATCH
    return alpha;
}

// void Matrix::Set( Int i, Int j, T alpha )
// -----------------------------------------

void ElMatrixSet_s( ElMatrix_s* AHandle, ElInt i, ElInt j, float alpha )
{
    try { RCM_s(AHandle)->Set(i,j,alpha); }
    CATCH
}

void ElMatrixSet_d( ElMatrix_d* AHandle, ElInt i, ElInt j, double alpha )
{
    try { RCM_d(AHandle)->Set(i,j,alpha); }
    CATCH
}

void ElMatrixSet_c( ElMatrix_c* AHandle, ElInt i, ElInt j, void* alpha )
{
    try { RCM_c(AHandle)->Set(i,j,*RCB_c(alpha)); }
    CATCH
}

void ElMatrixSet_z( ElMatrix_z* AHandle, ElInt i, ElInt j, void* alpha )
{
    try { RCM_z(AHandle)->Set(i,j,*RCB_z(alpha)); }
    CATCH
}

// void Matrix::SetRealPart( Int i, Int j, Base<T> alpha )
// -------------------------------------------------------

void ElMatrixSetRealPart_c
( ElMatrix_c* AHandle, ElInt i, ElInt j, float alpha )
{
    try { RCM_c(AHandle)->SetRealPart(i,j,alpha); }
    CATCH
}

void ElMatrixSetRealPart_z
( ElMatrix_z* AHandle, ElInt i, ElInt j, double alpha )
{
    try { RCM_z(AHandle)->SetRealPart(i,j,alpha); }
    CATCH
}

// void Matrix::SetImagPart( Int i, Int j, Base<T> alpha )
// -------------------------------------------------------

void ElMatrixSetImagPart_c
( ElMatrix_c* AHandle, ElInt i, ElInt j, float alpha )
{
    try { RCM_c(AHandle)->SetImagPart(i,j,alpha); }
    CATCH
}

void ElMatrixSetImagPart_z
( ElMatrix_z* AHandle, ElInt i, ElInt j, double alpha )
{
    try { RCM_z(AHandle)->SetImagPart(i,j,alpha); }
    CATCH
}

// void Matrix::Update( Int i, Int j, T alpha )
// --------------------------------------------

void ElMatrixUpdate_s( ElMatrix_s* AHandle, ElInt i, ElInt j, float alpha )
{
    try { RCM_s(AHandle)->Update(i,j,alpha); }
    CATCH
}

void ElMatrixUpdate_d( ElMatrix_d* AHandle, ElInt i, ElInt j, double alpha )
{
    try { RCM_d(AHandle)->Update(i,j,alpha); }
    CATCH
}

void ElMatrixUpdate_c( ElMatrix_c* AHandle, ElInt i, ElInt j, void* alpha )
{
    try { RCM_c(AHandle)->Update(i,j,*RCB_c(alpha)); }
    CATCH
}

void ElMatrixUpdate_z( ElMatrix_z* AHandle, ElInt i, ElInt j, void* alpha )
{
    try { RCM_z(AHandle)->Update(i,j,*RCB_z(alpha)); }
    CATCH
}

// void Matrix::UpdateRealPart( Int i, Int j, Base<T> alpha )
// ----------------------------------------------------------

void ElMatrixUpdateRealPart_c
( ElMatrix_c* AHandle, ElInt i, ElInt j, float alpha )
{
    try { RCM_c(AHandle)->UpdateRealPart(i,j,alpha); }
    CATCH
}

void ElMatrixUpdateRealPart_z
( ElMatrix_z* AHandle, ElInt i, ElInt j, double alpha )
{
    try { RCM_z(AHandle)->UpdateRealPart(i,j,alpha); }
    CATCH
}

// void Matrix::UpdateImagPart( Int i, Int j, Base<T> alpha )
// ----------------------------------------------------------

void ElMatrixUpdateImagPart_c
( ElMatrix_c* AHandle, ElInt i, ElInt j, float alpha )
{
    try { RCM_c(AHandle)->UpdateImagPart(i,j,alpha); }
    CATCH
}

void ElMatrixUpdateImagPart_z
( ElMatrix_z* AHandle, ElInt i, ElInt j, double alpha )
{
    try { RCM_z(AHandle)->UpdateImagPart(i,j,alpha); }
    CATCH
}

// void Matrix::MakeReal( Int i, Int j )
// -------------------------------------

// TODO

} // extern "C"
