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

#define CATCH \
  catch( std::bad_alloc& e ) \
  { ReportException(e); return EL_ALLOC_ERROR; } \
  catch( std::logic_error& e ) \
  { ReportException(e); return EL_LOGIC_ERROR; } \
  catch( std::runtime_error& e ) \
  { ReportException(e); return EL_RUNTIME_ERROR; } \
  catch( std::exception& e ) \
  { ReportException(e); return EL_ERROR; }

extern "C" {

// Matrix<T>::Matrix()
// -------------------
ElError ElMatrixCreate_s( ElMatrix_s* A )
{
    try { *A = (ElMatrix_s)reinterpret_cast<struct ElMatrix_sDummy*>
               ( new Matrix<float>() ); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixCreate_d( ElMatrix_d* A )
{
    try { *A = (ElMatrix_d)reinterpret_cast<struct ElMatrix_dDummy*>
               ( new Matrix<double>() ); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixCreate_c( ElMatrix_c* A )
{
    try 
    { *A = (ElMatrix_c)reinterpret_cast<struct ElMatrix_cDummy*>
           ( new Matrix<Complex<float>>() ); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixCreate_z( ElMatrix_z* A )
{
    try 
    { *A = (ElMatrix_z)reinterpret_cast<struct ElMatrix_zDummy*>
           ( new Matrix<Complex<double>>() ); }
    CATCH
    return EL_SUCCESS;
}

// Matrix<T>::~Matrix()
// --------------------
ElError ElMatrixDestroy_s( ElConstMatrix_s AHandle )
{ 
    try { delete RCM_s_const(AHandle); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixDestroy_d( ElConstMatrix_d AHandle )
{ 
    try { delete RCM_d_const(AHandle); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixDestroy_c( ElConstMatrix_c AHandle )
{ 
    try { delete RCM_c_const(AHandle); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixDestroy_z( ElConstMatrix_z AHandle )
{ 
    try { delete RCM_z_const(AHandle); }
    CATCH
    return EL_SUCCESS;
}

// void Matrix<T>::Empty()
// -----------------------
ElError ElMatrixEmpty_s( ElMatrix_s AHandle )
{ 
    try { RCM_s(AHandle)->Empty(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixEmpty_d( ElMatrix_d AHandle )
{ 
    try { RCM_d(AHandle)->Empty(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixEmpty_c( ElMatrix_c AHandle )
{ 
    try { RCM_c(AHandle)->Empty(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixEmpty_z( ElMatrix_z AHandle )
{ 
    try { RCM_z(AHandle)->Empty(); }
    CATCH
    return EL_SUCCESS;
}

// void Matrix<T>::Resize( Int height, Int width )
// -----------------------------------------------
ElError ElMatrixResize_s( ElMatrix_s AHandle, ElInt height, ElInt width )
{
    try { RCM_s(AHandle)->Resize(height,width); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixResize_d( ElMatrix_d AHandle, ElInt height, ElInt width )
{
    try { RCM_d(AHandle)->Resize(height,width); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixResize_c( ElMatrix_c AHandle, ElInt height, ElInt width )
{
    try { RCM_c(AHandle)->Resize(height,width); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixResize_z( ElMatrix_z AHandle, ElInt height, ElInt width )
{
    try { RCM_z(AHandle)->Resize(height,width); }
    CATCH
    return EL_SUCCESS;
}

// void Matrix<T>::Resize( Int height, Int width, Int ldim )
// ---------------------------------------------------------
ElError ElMatrixResizeWithLDim_s
( ElMatrix_s AHandle, ElInt height, ElInt width, ElInt ldim )
{
    try { RCM_s(AHandle)->Resize(height,width,ldim); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixResizeWithLDim_d
( ElMatrix_d AHandle, ElInt height, ElInt width, ElInt ldim )
{
    try { RCM_d(AHandle)->Resize(height,width,ldim); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixResizeWithLDim_c
( ElMatrix_c AHandle, ElInt height, ElInt width, ElInt ldim )
{
    try { RCM_c(AHandle)->Resize(height,width,ldim); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixResizeWithLDim_z
( ElMatrix_z AHandle, ElInt height, ElInt width, ElInt ldim )
{
    try { RCM_z(AHandle)->Resize(height,width,ldim); }
    CATCH
    return EL_SUCCESS;
}

// void Matrix<T>::Attach( Int height, Int width, T* buffer, Int ldim )
// --------------------------------------------------------------------
ElError ElMatrixAttach_s
( ElMatrix_s AHandle, ElInt height, ElInt width, float* buffer, ElInt ldim )
{
    try { RCM_s(AHandle)->Attach(height,width,buffer,ldim); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixAttach_d
( ElMatrix_d AHandle, ElInt height, ElInt width, double* buffer, ElInt ldim )
{
    try { RCM_d(AHandle)->Attach(height,width,buffer,ldim); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixAttach_c
( ElMatrix_c AHandle, ElInt height, ElInt width, 
  complex_float* buffer, ElInt ldim )
{
    try { RCM_c(AHandle)->Attach(height,width,RCB_c(buffer),ldim); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixAttach_z
( ElMatrix_z AHandle, ElInt height, ElInt width, 
  complex_double* buffer, ElInt ldim )
{
    try { RCM_z(AHandle)->Attach(height,width,RCB_z(buffer),ldim); }
    CATCH
    return EL_SUCCESS;
}

// void Matrix<T>::LockedAttach
// ( Int height, Int width, const T* buffer, Int ldim )
// ----------------------------------------------------
ElError ElMatrixLockedAttach_s
( ElMatrix_s AHandle, 
  ElInt height, ElInt width, const float* buffer, ElInt ldim )
{
    try { RCM_s(AHandle)->LockedAttach(height,width,buffer,ldim); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixLockedAttach_d
( ElMatrix_d AHandle, 
  ElInt height, ElInt width, const double* buffer, ElInt ldim )
{
    try { RCM_d(AHandle)->LockedAttach(height,width,buffer,ldim); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixLockedAttach_c
( ElMatrix_c AHandle, 
  ElInt height, ElInt width, const complex_float* buffer, ElInt ldim )
{
    try { RCM_c(AHandle)->LockedAttach(height,width,RCB_c_const(buffer),ldim); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixLockedAttach_z
( ElMatrix_z AHandle, 
  ElInt height, ElInt width, const complex_double* buffer, ElInt ldim )
{
    try { RCM_z(AHandle)->LockedAttach(height,width,RCB_z_const(buffer),ldim); }
    CATCH
    return EL_SUCCESS;
}

// void Matrix<T>::Control( Int height, Int width, T* buffer, Int ldim )
// ---------------------------------------------------------------------
ElError ElMatrixControl_s
( ElMatrix_s AHandle, ElInt height, ElInt width, float* buffer, ElInt ldim )
{
    try { RCM_s(AHandle)->Control(height,width,buffer,ldim); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixControl_d
( ElMatrix_d AHandle, ElInt height, ElInt width, double* buffer, ElInt ldim )
{
    try { RCM_d(AHandle)->Control(height,width,buffer,ldim); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixControl_c
( ElMatrix_c AHandle, ElInt height, ElInt width, 
  complex_float* buffer, ElInt ldim )
{
    try { RCM_c(AHandle)->Control(height,width,RCB_c(buffer),ldim); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixControl_z
( ElMatrix_z AHandle, ElInt height, ElInt width, 
  complex_double* buffer, ElInt ldim )
{
    try { RCM_z(AHandle)->Control(height,width,RCB_z(buffer),ldim); }
    CATCH
    return EL_SUCCESS;
}

// B = A
// -----
ElError ElMatrixCopy_s( ElConstMatrix_s AHandle, ElMatrix_s BHandle )
{
    try { *RCM_s(BHandle) = *RCM_s_const(AHandle); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixCopy_d( ElConstMatrix_d AHandle, ElMatrix_d BHandle )
{
    try { *RCM_d(BHandle) = *RCM_d_const(AHandle); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixCopy_c( ElConstMatrix_c AHandle, ElMatrix_c BHandle )
{
    try { *RCM_c(BHandle) = *RCM_c_const(AHandle); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixCopy_z( ElConstMatrix_z AHandle, ElMatrix_z BHandle )
{
    try { *RCM_z(BHandle) = *RCM_z_const(AHandle); }
    CATCH
    return EL_SUCCESS;
}

// Int Matrix<T>::Height() const
// -----------------------------
ElError ElMatrixHeight_s( ElConstMatrix_s AHandle, ElInt* height )
{ 
    try { *height = RCM_s_const(AHandle)->Height(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixHeight_d( ElConstMatrix_d AHandle, ElInt* height )
{ 
    try { *height = RCM_d_const(AHandle)->Height(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixHeight_c( ElConstMatrix_c AHandle, ElInt* height )
{ 
    try { *height = RCM_c_const(AHandle)->Height(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixHeight_z( ElConstMatrix_z AHandle, ElInt* height )
{ 
    try { *height = RCM_z_const(AHandle)->Height(); }
    CATCH
    return EL_SUCCESS;
}

// Int Matrix<T>::Width() const
// ----------------------------
ElError ElMatrixWidth_s( ElConstMatrix_s AHandle, ElInt* width )
{ 
    try { *width = RCM_s_const(AHandle)->Width(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixWidth_d( ElConstMatrix_d AHandle, ElInt* width )
{ 
    try { *width = RCM_d_const(AHandle)->Width(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixWidth_c( ElConstMatrix_c AHandle, ElInt* width )
{ 
    try { *width = RCM_c_const(AHandle)->Width(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixWidth_z( ElConstMatrix_z AHandle, ElInt* width )
{ 
    try { *width = RCM_z_const(AHandle)->Width(); }
    CATCH
    return EL_SUCCESS;
}

// Int Matrix<T>::LDim() const
// ---------------------------
ElError ElMatrixLDim_s( ElConstMatrix_s AHandle, ElInt* ldim )
{ 
    try { *ldim = RCM_s_const(AHandle)->LDim(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixLDim_d( ElConstMatrix_d AHandle, ElInt* ldim )
{ 
    try { *ldim = RCM_d_const(AHandle)->LDim(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixLDim_c( ElConstMatrix_c AHandle, ElInt* ldim )
{ 
    try { *ldim = RCM_c_const(AHandle)->LDim(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixLDim_z( ElConstMatrix_z AHandle, ElInt* ldim )
{ 
    try { *ldim = RCM_z_const(AHandle)->LDim(); }
    CATCH
    return EL_SUCCESS;
}

// Int Matrix<T>::MemorySize() const
// ---------------------------------
ElError ElMatrixMemorySize_s( ElConstMatrix_s AHandle, ElInt* mem )
{ 
    try { *mem = RCM_s_const(AHandle)->MemorySize(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixMemorySize_d( ElConstMatrix_d AHandle, ElInt* mem )
{ 
    try { *mem = RCM_d_const(AHandle)->MemorySize(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixMemorySize_c( ElConstMatrix_c AHandle, ElInt* mem )
{ 
    try { *mem = RCM_c_const(AHandle)->MemorySize(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixMemorySize_z( ElConstMatrix_z AHandle, ElInt* mem )
{ 
    try { *mem = RCM_z_const(AHandle)->MemorySize(); }
    CATCH
    return EL_SUCCESS;
}

// Int Matrix<T>::DiagonalLength( Int offset ) const
// -------------------------------------------------
ElError ElMatrixDiagonalLength_s
( ElConstMatrix_s AHandle, ElInt offset, ElInt* length )
{ 
    try { *length = RCM_s_const(AHandle)->DiagonalLength(offset); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixDiagonalLength_d
( ElConstMatrix_d AHandle, ElInt offset, ElInt* length )
{ 
    try { *length = RCM_d_const(AHandle)->DiagonalLength(offset); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixDiagonalLength_c
( ElConstMatrix_c AHandle, ElInt offset, ElInt* length )
{ 
    try { *length = RCM_c_const(AHandle)->DiagonalLength(offset); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixDiagonalLength_z
( ElConstMatrix_z AHandle, ElInt offset, ElInt* length )
{ 
    try { *length = RCM_z_const(AHandle)->DiagonalLength(offset); }
    CATCH
    return EL_SUCCESS;
}

// T* Matrix<T>::Buffer()
// ----------------------
ElError ElMatrixBuffer_s( ElMatrix_s AHandle, float** buffer )
{ 
    try { *buffer = RCM_s(AHandle)->Buffer(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixBuffer_d( ElMatrix_d AHandle, double** buffer )
{ 
    try { *buffer = RCM_d(AHandle)->Buffer(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixBuffer_c( ElMatrix_c AHandle, complex_float** buffer )
{ 
    try { *buffer = (complex_float*)RCM_c(AHandle)->Buffer(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixBuffer_z( ElMatrix_z AHandle, complex_double** buffer )
{ 
    try { *buffer = (complex_double*)RCM_z(AHandle)->Buffer(); }
    CATCH
    return EL_SUCCESS;
}

// const T* Matrix<T>::LockedBuffer() const
// ----------------------------------------
ElError ElMatrixLockedBuffer_s
( ElConstMatrix_s AHandle, const float** buffer )
{ 
    try { *buffer = RCM_s_const(AHandle)->LockedBuffer(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixLockedBuffer_d
( ElConstMatrix_d AHandle, const double** buffer )
{ 
    try { *buffer = RCM_d_const(AHandle)->LockedBuffer(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixLockedBuffer_c
( ElConstMatrix_c AHandle, const complex_float** buffer )
{ 
    try 
    { *buffer = (const complex_float*)RCM_c_const(AHandle)->LockedBuffer(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixLockedBuffer_z
( ElConstMatrix_z AHandle, const complex_double** buffer )
{ 
    try 
    { *buffer = (const complex_double*)RCM_z_const(AHandle)->LockedBuffer(); }
    CATCH
    return EL_SUCCESS;
}

// bool Matrix<T>::Viewing() const
// -------------------------------
ElError ElMatrixViewing_s( ElConstMatrix_s AHandle, bool* viewing )
{ 
    try { *viewing = RCM_s_const(AHandle)->Viewing(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixViewing_d( ElConstMatrix_d AHandle, bool* viewing )
{ 
    try { *viewing = RCM_d_const(AHandle)->Viewing(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixViewing_c( ElConstMatrix_c AHandle, bool* viewing )
{ 
    try { *viewing = RCM_c_const(AHandle)->Viewing(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixViewing_z( ElConstMatrix_z AHandle, bool* viewing )
{ 
    try { *viewing = RCM_z_const(AHandle)->Viewing(); }
    CATCH
    return EL_SUCCESS;
}

// bool Matrix<T>::FixedSize() const
// ---------------------------------
ElError ElMatrixFixedSize_s( ElConstMatrix_s AHandle, bool* fixedSize )
{ 
    try { *fixedSize = RCM_s_const(AHandle)->FixedSize(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixFixedSize_d( ElConstMatrix_d AHandle, bool* fixedSize )
{ 
    try { *fixedSize = RCM_d_const(AHandle)->FixedSize(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixFixedSize_c( ElConstMatrix_c AHandle, bool* fixedSize )
{ 
    try { *fixedSize = RCM_c_const(AHandle)->FixedSize(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixFixedSize_z( ElConstMatrix_z AHandle, bool* fixedSize )
{ 
    try { *fixedSize = RCM_z_const(AHandle)->FixedSize(); }
    CATCH
    return EL_SUCCESS;
}

// bool Matrix<T>::Locked() const
// ------------------------------
ElError ElMatrixLocked_s( ElConstMatrix_s AHandle, bool* locked )
{ 
    try { *locked = RCM_s_const(AHandle)->Locked(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixLocked_d( ElConstMatrix_d AHandle, bool* locked )
{ 
    try { *locked = RCM_d_const(AHandle)->Locked(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixLocked_c( ElConstMatrix_c AHandle, bool* locked )
{ 
    try { *locked = RCM_c_const(AHandle)->Locked(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixLocked_z( ElConstMatrix_z AHandle, bool* locked )
{ 
    try { *locked = RCM_z_const(AHandle)->Locked(); }
    CATCH
    return EL_SUCCESS;
}

// T Matrix<T>::Get( Int i, Int j ) const
// --------------------------------------
ElError ElMatrixGet_s
( ElConstMatrix_s AHandle, ElInt i, ElInt j, float* val )
{
    try { *val = RCM_s_const(AHandle)->Get(i,j); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixGet_d
( ElConstMatrix_d AHandle, ElInt i, ElInt j, double* val )
{
    try { *val = RCM_s_const(AHandle)->Get(i,j); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixGet_c
( ElConstMatrix_c AHandle, ElInt i, ElInt j, complex_float* val )
{
    try 
    { 
        Complex<float> alpha = RCM_c_const(AHandle)->Get(i,j); 
        val->real = alpha.real();
        val->imag = alpha.imag();
    }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixGet_z
( ElConstMatrix_z AHandle, ElInt i, ElInt j, complex_double* val )
{
    try 
    { 
        Complex<double> alpha = RCM_z_const(AHandle)->Get(i,j); 
        val->real = alpha.real();
        val->imag = alpha.imag();
    }
    CATCH
    return EL_SUCCESS;
}

// Base<T> Matrix<T>::GetRealPart( Int i, Int j ) const
// ----------------------------------------------------
ElError ElMatrixGetRealPart_c
( ElConstMatrix_c AHandle, ElInt i, ElInt j, float* val )
{
    try { *val = RCM_s_const(AHandle)->GetRealPart(i,j); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixGetRealPart_z
( ElConstMatrix_z AHandle, ElInt i, ElInt j, double* val )
{
    try { *val = RCM_d_const(AHandle)->GetRealPart(i,j); }
    CATCH
    return EL_SUCCESS;
}

// Base<T> Matrix<T>::GetImagPart( Int i, Int j ) const
// ----------------------------------------------------
ElError ElMatrixGetImagPart_c
( ElConstMatrix_c AHandle, ElInt i, ElInt j, float* val )
{
    try { *val = RCM_s_const(AHandle)->GetImagPart(i,j); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixGetImagPart_z
( ElConstMatrix_z AHandle, ElInt i, ElInt j, double* val )
{
    try { *val = RCM_d_const(AHandle)->GetImagPart(i,j); }
    CATCH
    return EL_SUCCESS;
}

// void Matrix<T>::Set( Int i, Int j, T alpha )
// --------------------------------------------
ElError ElMatrixSet_s( ElMatrix_s AHandle, ElInt i, ElInt j, float alpha )
{
    try { RCM_s(AHandle)->Set(i,j,alpha); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixSet_d( ElMatrix_d AHandle, ElInt i, ElInt j, double alpha )
{
    try { RCM_d(AHandle)->Set(i,j,alpha); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixSet_c
( ElMatrix_c AHandle, ElInt i, ElInt j, complex_float alpha )
{
    try { RCM_c(AHandle)->Set(i,j,Complex<float>(alpha.real,alpha.imag)); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixSet_z
( ElMatrix_z AHandle, ElInt i, ElInt j, complex_double alpha )
{
    try { RCM_z(AHandle)->Set(i,j,Complex<double>(alpha.real,alpha.imag)); }
    CATCH
    return EL_SUCCESS;
}

// void Matrix<T>::SetRealPart( Int i, Int j, Base<T> alpha )
// ----------------------------------------------------------
ElError ElMatrixSetRealPart_c
( ElMatrix_c AHandle, ElInt i, ElInt j, float alpha )
{
    try { RCM_c(AHandle)->SetRealPart(i,j,alpha); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixSetRealPart_z
( ElMatrix_z AHandle, ElInt i, ElInt j, double alpha )
{
    try { RCM_z(AHandle)->SetRealPart(i,j,alpha); }
    CATCH
    return EL_SUCCESS;
}

// void Matrix<T>::SetImagPart( Int i, Int j, Base<T> alpha )
// ----------------------------------------------------------
ElError ElMatrixSetImagPart_c
( ElMatrix_c AHandle, ElInt i, ElInt j, float alpha )
{
    try { RCM_c(AHandle)->SetImagPart(i,j,alpha); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixSetImagPart_z
( ElMatrix_z AHandle, ElInt i, ElInt j, double alpha )
{
    try { RCM_z(AHandle)->SetImagPart(i,j,alpha); }
    CATCH
    return EL_SUCCESS;
}

// void Matrix<T>::Update( Int i, Int j, T alpha )
// -----------------------------------------------
ElError ElMatrixUpdate_s
( ElMatrix_s AHandle, ElInt i, ElInt j, float alpha )
{
    try { RCM_s(AHandle)->Update(i,j,alpha); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixUpdate_d
( ElMatrix_d AHandle, ElInt i, ElInt j, double alpha )
{
    try { RCM_d(AHandle)->Update(i,j,alpha); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixUpdate_c
( ElMatrix_c AHandle, ElInt i, ElInt j, complex_float alpha )
{
    try { RCM_c(AHandle)->Update(i,j,Complex<float>(alpha.real,alpha.imag)); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixUpdate_z
( ElMatrix_z AHandle, ElInt i, ElInt j, complex_double alpha )
{
    try { RCM_z(AHandle)->Update(i,j,Complex<double>(alpha.real,alpha.imag)); }
    CATCH
    return EL_SUCCESS;
}

// void Matrix<T>::UpdateRealPart( Int i, Int j, Base<T> alpha )
// -------------------------------------------------------------
ElError ElMatrixUpdateRealPart_c
( ElMatrix_c AHandle, ElInt i, ElInt j, float alpha )
{
    try { RCM_c(AHandle)->UpdateRealPart(i,j,alpha); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixUpdateRealPart_z
( ElMatrix_z AHandle, ElInt i, ElInt j, double alpha )
{
    try { RCM_z(AHandle)->UpdateRealPart(i,j,alpha); }
    CATCH
    return EL_SUCCESS;
}

// void Matrix<T>::UpdateImagPart( Int i, Int j, Base<T> alpha )
// -------------------------------------------------------------
ElError ElMatrixUpdateImagPart_c
( ElMatrix_c AHandle, ElInt i, ElInt j, float alpha )
{
    try { RCM_c(AHandle)->UpdateImagPart(i,j,alpha); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixUpdateImagPart_z
( ElMatrix_z AHandle, ElInt i, ElInt j, double alpha )
{
    try { RCM_z(AHandle)->UpdateImagPart(i,j,alpha); }
    CATCH
    return EL_SUCCESS;
}

// void Matrix<T>::MakeReal( Int i, Int j )
// ----------------------------------------
ElError ElMatrixMakeReal_c( ElMatrix_c AHandle, ElInt i, ElInt j )
{ 
    try { RCM_c(AHandle)->MakeReal(i,j); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixMakeReal_z( ElMatrix_z AHandle, ElInt i, ElInt j )
{ 
    try { RCM_z(AHandle)->MakeReal(i,j); }
    CATCH
    return EL_SUCCESS;
}

// void Matrix<T>::Conjugate( Int i, Int j )
// -----------------------------------------
ElError ElMatrixConjugate_c( ElMatrix_c AHandle, ElInt i, ElInt j )
{
    try { RCM_c(AHandle)->Conjugate(i,j); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixConjugate_z( ElMatrix_z AHandle, ElInt i, ElInt j )
{
    try { RCM_z(AHandle)->Conjugate(i,j); }
    CATCH
    return EL_SUCCESS;
}

// Matrix<T> Matrix<T>::GetDiagonal( Int offset ) const
// ----------------------------------------------------
ElError ElMatrixGetDiagonal_s
( ElConstMatrix_s AHandle, ElInt offset, ElMatrix_s* dHandle )
{
    try 
    {
        auto d = new Matrix<float>;
        RCM_s_const(AHandle)->GetDiagonal( *d, offset );
        *dHandle = (ElMatrix_s)d;
    }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixGetDiagonal_d
( ElConstMatrix_d AHandle, ElInt offset, ElMatrix_d* dHandle )
{
    try 
    {
        auto d = new Matrix<double>;
        RCM_d_const(AHandle)->GetDiagonal( *d, offset );
        *dHandle = (ElMatrix_d)d;
    }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixGetDiagonal_c
( ElConstMatrix_c AHandle, ElInt offset, ElMatrix_c* dHandle )
{
    try 
    {
        auto d = new Matrix<Complex<float>>;
        RCM_c_const(AHandle)->GetDiagonal( *d, offset );
        *dHandle = (ElMatrix_c)d;
    }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixGetDiagonal_z
( ElConstMatrix_z AHandle, ElInt offset, ElMatrix_z* dHandle )
{
    try 
    {
        auto d = new Matrix<Complex<double>>;
        RCM_z_const(AHandle)->GetDiagonal( *d, offset );
        *dHandle = (ElMatrix_z)d;
    }
    CATCH
    return EL_SUCCESS;
}

// Matrix<Base<T>> Matrix<T>::GetRealPartOfDiagonal( Int offset ) const
// --------------------------------------------------------------------
ElError ElMatrixGetRealPartOfDiagonal_c
( ElConstMatrix_c AHandle, ElInt offset, ElMatrix_s* dHandle )
{
    try 
    {
        auto d = new Matrix<float>;
        RCM_c_const(AHandle)->GetRealPartOfDiagonal( *d, offset );
        *dHandle = (ElMatrix_s)d;
    }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixGetRealPartOfDiagonal_z
( ElConstMatrix_z AHandle, ElInt offset, ElMatrix_d* dHandle )
{
    try 
    {
        auto d = new Matrix<double>;
        RCM_z_const(AHandle)->GetRealPartOfDiagonal( *d, offset );
        *dHandle = (ElMatrix_d)d;
    }
    CATCH
    return EL_SUCCESS;
}

// Matrix<Base<T>> Matrix<T>::GetImagPartOfDiagonal( Int offset ) const
// --------------------------------------------------------------------
ElError ElMatrixGetImagPartOfDiagonal_c
( ElConstMatrix_c AHandle, ElInt offset, ElMatrix_s* dHandle )
{
    try 
    {
        auto d = new Matrix<float>;
        RCM_c_const(AHandle)->GetImagPartOfDiagonal( *d, offset );
        *dHandle = (ElMatrix_s)d;
    }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixGetImagPartOfDiagonal_z
( ElConstMatrix_z AHandle, ElInt offset, ElMatrix_d* dHandle )
{
    try 
    {
        auto d = new Matrix<double>;
        RCM_z_const(AHandle)->GetImagPartOfDiagonal( *d, offset );
        *dHandle = (ElMatrix_d)d;
    }
    CATCH
    return EL_SUCCESS;
}

// void Matrix<T>::SetDiagonal( const Matrix<T>& d, Int offset )
// -------------------------------------------------------------
ElError ElMatrixSetDiagonal_s
( ElMatrix_s AHandle, ElConstMatrix_s dHandle, ElInt offset )
{
    try { RCM_s(AHandle)->SetDiagonal( RCM_s_const(dHandle), offset ); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixSetDiagonal_d
( ElMatrix_d AHandle, ElConstMatrix_d dHandle, ElInt offset )
{
    try { RCM_d(AHandle)->SetDiagonal( RCM_d_const(dHandle), offset ); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixSetDiagonal_c
( ElMatrix_c AHandle, ElConstMatrix_c dHandle, ElInt offset )
{
    try { RCM_c(AHandle)->SetDiagonal( RCM_c_const(dHandle), offset ); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixSetDiagonal_z
( ElMatrix_z AHandle, ElConstMatrix_z dHandle, ElInt offset )
{
    try { RCM_z(AHandle)->SetDiagonal( RCM_z_const(dHandle), offset ); }
    CATCH
    return EL_SUCCESS;
}

// void Matrix<T>::SetRealPartOfDiagonal( const Matrix<Base<T>>& d, Int offset )
// -----------------------------------------------------------------------------
ElError ElMatrixSetRealPartOfDiagonal_c
( ElMatrix_c AHandle, ElConstMatrix_s dHandle, ElInt offset )
{
    try { RCM_c(AHandle)->SetRealPartOfDiagonal
          ( RCM_s_const(dHandle), offset ); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixSetRealPartOfDiagonal_z
( ElMatrix_z AHandle, ElConstMatrix_d dHandle, ElInt offset )
{
    try { RCM_z(AHandle)->SetRealPartOfDiagonal
          ( RCM_d_const(dHandle), offset ); }
    CATCH
    return EL_SUCCESS;
}

// void Matrix<T>::SetImagPartOfDiagonal( const Matrix<Base<T>>& d, Int offset )
// -----------------------------------------------------------------------------
ElError ElMatrixSetImagPartOfDiagonal_c
( ElMatrix_c AHandle, ElConstMatrix_s dHandle, ElInt offset )
{
    try { RCM_c(AHandle)->SetImagPartOfDiagonal
          ( RCM_s_const(dHandle), offset ); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixSetImagPartOfDiagonal_z
( ElMatrix_z AHandle, ElConstMatrix_d dHandle, ElInt offset )
{
    try { RCM_z(AHandle)->SetImagPartOfDiagonal
          ( RCM_d_const(dHandle), offset ); }
    CATCH
    return EL_SUCCESS;
}

// void Matrix<T>::UpdateDiagonal( const Matrix<T>& d, Int offset )
// ----------------------------------------------------------------
ElError ElMatrixUpdateDiagonal_s
( ElMatrix_s AHandle, ElConstMatrix_s dHandle, ElInt offset )
{
    try { RCM_s(AHandle)->UpdateDiagonal( RCM_s_const(dHandle), offset ); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixUpdateDiagonal_d
( ElMatrix_d AHandle, ElConstMatrix_d dHandle, ElInt offset )
{
    try { RCM_d(AHandle)->UpdateDiagonal( RCM_d_const(dHandle), offset ); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixUpdateDiagonal_c
( ElMatrix_c AHandle, ElConstMatrix_c dHandle, ElInt offset )
{
    try { RCM_c(AHandle)->UpdateDiagonal( RCM_c_const(dHandle), offset ); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixUpdateDiagonal_z
( ElMatrix_z AHandle, ElConstMatrix_z dHandle, ElInt offset )
{
    try { RCM_z(AHandle)->UpdateDiagonal( RCM_z_const(dHandle), offset ); }
    CATCH
    return EL_SUCCESS;
}

// void Matrix<T>::UpdateRealPartOfDiagonal
// ( const Matrix<Base<T>>& d, Int offset )
// ----------------------------------------
ElError ElMatrixUpdateRealPartOfDiagonal_c
( ElMatrix_c AHandle, ElConstMatrix_s dHandle, ElInt offset )
{
    try { RCM_c(AHandle)->UpdateRealPartOfDiagonal
          ( RCM_s_const(dHandle), offset ); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixUpdateRealPartOfDiagonal_z
( ElMatrix_z AHandle, ElConstMatrix_d dHandle, ElInt offset )
{
    try { RCM_z(AHandle)->UpdateRealPartOfDiagonal
          ( RCM_d_const(dHandle), offset ); }
    CATCH
    return EL_SUCCESS;
}

// void Matrix<T>::UpdateImagPartOfDiagonal
// ( const Matrix<Base<T>>& d, Int offset )
// ----------------------------------------
ElError ElMatrixUpdateImagPartOfDiagonal_c
( ElMatrix_c AHandle, ElConstMatrix_s dHandle, ElInt offset )
{
    try { RCM_c(AHandle)->UpdateImagPartOfDiagonal
          ( RCM_s_const(dHandle), offset ); }
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixUpdateImagPartOfDiagonal_z
( ElMatrix_z AHandle, ElConstMatrix_d dHandle, ElInt offset )
{
    try { RCM_z(AHandle)->UpdateImagPartOfDiagonal
          ( RCM_d_const(dHandle), offset ); }
    CATCH
    return EL_SUCCESS;
}

// void Matrix<T>::MakeDiaogonalReal( Int offset )
// -----------------------------------------------
ElError ElMatrixMakeDiagonalReal_c( ElMatrix_c AHandle, ElInt offset )
{ RCM_c(AHandle)->MakeDiagonalReal(offset); return EL_SUCCESS; }

ElError ElMatrixMakeDiagonalReal_z( ElMatrix_z AHandle, ElInt offset )
{ RCM_z(AHandle)->MakeDiagonalReal(offset); return EL_SUCCESS; }

// void Matrix<T>::ConjugateDiagonal Int offset )
// ----------------------------------------------
ElError ElMatrixConjugateDiagonal_c( ElMatrix_c AHandle, ElInt offset )
{ RCM_c(AHandle)->ConjugateDiagonal(offset); return EL_SUCCESS; }

ElError ElMatrixConjugateDiagonal_z( ElMatrix_z AHandle, ElInt offset )
{ RCM_z(AHandle)->ConjugateDiagonal(offset); return EL_SUCCESS; }

// Matrix<T> Matrix<T>::GetSubmatrix
// ( const std::vector<Int>& rowInds, const std::vector<Int>& colInds ) const
// --------------------------------------------------------------------------
ElError ElMatrixGetSubmatrix_s
( ElConstMatrix_s AHandle,
  ElInt numRowInds, const ElInt* rowInds, 
  ElInt numColInds, const ElInt* colInds, ElMatrix_s* ASubHandle )
{
    try 
    { 
        std::vector<Int> rowIndVec( rowInds, rowInds+numRowInds ),
                         colIndVec( colInds, colInds+numColInds );
        auto ASub = new Matrix<float>;
        RCM_s_const(AHandle)->GetSubmatrix( rowIndVec, colIndVec, *ASub );
        *ASubHandle = (ElMatrix_s)ASub;
    } 
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixGetSubmatrix_d
( ElConstMatrix_d AHandle,
  ElInt numRowInds, const ElInt* rowInds, 
  ElInt numColInds, const ElInt* colInds, ElMatrix_d* ASubHandle )
{
    try 
    { 
        std::vector<Int> rowIndVec( rowInds, rowInds+numRowInds ),
                         colIndVec( colInds, colInds+numColInds );
        auto ASub = new Matrix<double>;
        RCM_d_const(AHandle)->GetSubmatrix( rowIndVec, colIndVec, *ASub );
        *ASubHandle = (ElMatrix_d)ASub;
    } 
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixGetSubmatrix_c
( ElConstMatrix_c AHandle,
  ElInt numRowInds, const ElInt* rowInds, 
  ElInt numColInds, const ElInt* colInds, ElMatrix_c* ASubHandle )
{
    try 
    { 
        std::vector<Int> rowIndVec( rowInds, rowInds+numRowInds ),
                         colIndVec( colInds, colInds+numColInds );
        auto ASub = new Matrix<Complex<float>>;
        RCM_c_const(AHandle)->GetSubmatrix( rowIndVec, colIndVec, *ASub );
        *ASubHandle = (ElMatrix_c)ASub;
    } 
    CATCH
    return EL_SUCCESS;
}

ElError ElMatrixGetSubmatrix_z
( ElConstMatrix_z AHandle,
  ElInt numRowInds, const ElInt* rowInds, 
  ElInt numColInds, const ElInt* colInds, ElMatrix_z* ASubHandle )
{
    try 
    { 
        std::vector<Int> rowIndVec( rowInds, rowInds+numRowInds ),
                         colIndVec( colInds, colInds+numColInds );
        auto ASub = new Matrix<Complex<double>>;
        RCM_z_const(AHandle)->GetSubmatrix( rowIndVec, colIndVec, *ASub );
        *ASubHandle = (ElMatrix_z)ASub;
    } 
    CATCH
    return EL_SUCCESS;
}

} // extern "C"
