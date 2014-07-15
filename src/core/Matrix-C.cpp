/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
#include "El-C.h"
using namespace El;

#define CATCH \
  catch( std::bad_alloc& e ) \
  { ReportException(e); return EL_ALLOC_ERROR; } \
  catch( std::logic_error& e ) \
  { ReportException(e); return EL_LOGIC_ERROR; } \
  catch( std::runtime_error& e ) \
  { ReportException(e); return EL_RUNTIME_ERROR; } \
  catch( std::exception& e ) \
  { ReportException(e); return EL_ERROR; }

#define EL_TRY(payload) \
  try { payload; } CATCH \
  return EL_SUCCESS;

extern "C" {

#define MATRIX_CONSTRUCT(SIG,T) \
  /* Matrix<T>::Matrix() */ \
  ElError ElMatrixCreate_ ## SIG ( ElMatrix_ ## SIG * A ) \
  { EL_TRY( *A = Reinterpret( new Matrix<T> ) ) } \
  /* Matrix<T>::~Matrix() */ \
  ElError ElMatrixDestroy_ ## SIG ( ElConstMatrix_ ## SIG AHandle ) \
  { EL_TRY( delete Reinterpret(AHandle) ) }

#define MATRIX_RECONFIG(SIG,T) \
  /* void Matrix<T>::Empty() */ \
  ElError ElMatrixEmpty_ ## SIG ( ElMatrix_ ## SIG AHandle ) \
  { EL_TRY( Reinterpret(AHandle)->Empty() ) } \
  /* void Matrix<T>::Resize( Int height, Int width ) */ \
  ElError ElMatrixResize_ ## SIG \
  ( ElMatrix_ ## SIG AHandle, ElInt height, ElInt width ) \
  { EL_TRY( Reinterpret(AHandle)->Resize(height,width) ) } \
  /* void Matrix<T>::Resize( Int height, Int width, Int ldim ) */ \
  ElError ElMatrixResizeWithLDim_ ## SIG \
  ( ElMatrix_ ## SIG AHandle, ElInt height, ElInt width, ElInt ldim ) \
  { EL_TRY( Reinterpret(AHandle)->Resize(height,width,ldim) ) } \
  /* void Matrix<T>::Attach( Int height, Int width, T* buffer, Int ldim ) */ \
  ElError ElMatrixAttach_ ## SIG \
  ( ElMatrix_ ## SIG AHandle, \
    ElInt height, ElInt width, CREFLECT(T)* buffer, ElInt ldim ) \
  { EL_TRY\
    ( Reinterpret(AHandle)->Attach(height,width,Reinterpret(buffer),ldim) ) } \
  /* void Matrix<T>::LockedAttach
     ( Int height, Int width, const T* buffer, Int ldim ) */ \
  ElError ElMatrixLockedAttach_ ## SIG \
  ( ElMatrix_ ## SIG AHandle, \
    ElInt height, ElInt width, const CREFLECT(T)* buffer, ElInt ldim ) \
  { EL_TRY \
    ( Reinterpret(AHandle)->LockedAttach \
      (height,width,Reinterpret(buffer),ldim) ) } \
  /* void Matrix<T>::Control( Int height, Int width, T* buffer, Int ldim ) */ \
  ElError ElMatrixControl_ ## SIG \
  ( ElMatrix_ ## SIG AHandle, ElInt height, ElInt width, CREFLECT(T)* buffer, \
    ElInt ldim ) \
  { EL_TRY \
    ( Reinterpret(AHandle)->Control(height,width,Reinterpret(buffer),ldim) ) } \
  /* B = A */ \
  ElError ElMatrixCopy_ ## SIG \
  ( ElConstMatrix_ ## SIG AHandle, ElMatrix_ ## SIG BHandle ) \
  { EL_TRY( *Reinterpret(BHandle) = *Reinterpret(AHandle) ) }

#define MATRIX_BASIC(SIG,T) \
  /* Int Matrix<T>::Height() const */ \
  ElError ElMatrixHeight_ ## SIG \
  ( ElConstMatrix_ ## SIG AHandle, ElInt* height ) \
  { EL_TRY( *height = Reinterpret(AHandle)->Height() ) } \
  /* Int Matrix<T>::Width() const */ \
  ElError ElMatrixWidth_ ## SIG \
  ( ElConstMatrix_ ## SIG AHandle, ElInt* width ) \
  { EL_TRY( *width = Reinterpret(AHandle)->Width() ) } \
  /* Int Matrix<T>::LDim() const */ \
  ElError ElMatrixLDim_ ## SIG \
  ( ElConstMatrix_ ## SIG AHandle, ElInt* ldim ) \
  { EL_TRY( *ldim = Reinterpret(AHandle)->LDim() ) } \
  /* Int Matrix<T>::MemorySize() const */ \
  ElError ElMatrixMemorySize_ ## SIG \
  ( ElConstMatrix_ ## SIG AHandle, ElInt* mem ) \
  { EL_TRY( *mem = Reinterpret(AHandle)->MemorySize() ) } \
  /* Int Matrix<T>::DiagonalLength( Int offset ) const */ \
  ElError ElMatrixDiagonalLength_ ## SIG \
  ( ElConstMatrix_ ## SIG AHandle, ElInt offset, ElInt* length ) \
  { EL_TRY( *length = Reinterpret(AHandle)->DiagonalLength(offset) ) } \
  /* T* Matrix<T>::Buffer() */ \
  ElError ElMatrixBuffer_ ## SIG \
  ( ElMatrix_ ## SIG AHandle, CREFLECT(T)** buffer ) \
  { EL_TRY( *buffer = Reinterpret(Reinterpret(AHandle)->Buffer()) ) } \
  /* const T* Matrix<T>::LockedBuffer() const */ \
  ElError ElMatrixLockedBuffer_ ## SIG \
  ( ElConstMatrix_ ## SIG AHandle, const CREFLECT(T)** buffer ) \
  { EL_TRY( *buffer = Reinterpret(Reinterpret(AHandle)->LockedBuffer()) ) } \
  /* bool Matrix<T>::Viewing() const */ \
  ElError ElMatrixViewing_ ## SIG \
  ( ElConstMatrix_ ## SIG AHandle, bool* viewing ) \
  { EL_TRY( *viewing = Reinterpret(AHandle)->Viewing() ) } \
  /* bool Matrix<T>::FixedSize() const */ \
  ElError ElMatrixFixedSize_ ## SIG \
  ( ElConstMatrix_ ## SIG AHandle, bool* fixedSize ) \
  { EL_TRY( *fixedSize = Reinterpret(AHandle)->FixedSize() ) } \
  /* bool Matrix<T>::Locked() const */ \
  ElError ElMatrixLocked_ ## SIG \
  ( ElConstMatrix_ ## SIG AHandle, bool* locked ) \
  { EL_TRY( *locked = Reinterpret(AHandle)->Locked() ) }

#define MATRIX_SINGLEENTRY(SIG,T) \
  /* T Matrix<T>::Get( Int i, Int j ) const */ \
  ElError ElMatrixGet_ ## SIG \
  ( ElConstMatrix_ ## SIG AHandle, ElInt i, ElInt j, CREFLECT(T)* val ) \
  { EL_TRY( *val = Reinterpret(Reinterpret(AHandle)->Get(i,j)) ) } \
  /* void Matrix<T>::Set( Int i, Int j, T alpha ) */ \
  ElError ElMatrixSet_ ## SIG \
  ( ElMatrix_ ## SIG AHandle, ElInt i, ElInt j, CREFLECT(T) alpha ) \
  { EL_TRY( Reinterpret(AHandle)->Set(i,j,Reinterpret(alpha)) ) } \
  /* void Matrix<T>::Update( Int i, Int j, T alpha ) */ \
  ElError ElMatrixUpdate_ ## SIG \
  ( ElMatrix_ ## SIG AHandle, ElInt i, ElInt j, CREFLECT(T) alpha ) \
  { EL_TRY( Reinterpret(AHandle)->Update(i,j,Reinterpret(alpha)) ) }

#define MATRIX_SINGLEENTRY_COMPLEX(SIG,SIGBASE,T) \
  /* Base<T> Matrix<T>::GetRealPart( Int i, Int j ) const */ \
  ElError ElMatrixGetRealPart_ ## SIG \
  ( ElConstMatrix_ ## SIG AHandle, ElInt i, ElInt j, Base<T>* val ) \
  { EL_TRY( *val = Reinterpret(AHandle)->GetRealPart(i,j) ) } \
  /* Base<T> Matrix<T>::GetImagPart( Int i, Int j ) const */ \
  ElError ElMatrixGetImagPart_ ## SIG \
  ( ElConstMatrix_ ## SIG AHandle, ElInt i, ElInt j, Base<T>* val ) \
  { EL_TRY( *val = Reinterpret(AHandle)->GetImagPart(i,j) ) } \
  /* void Matrix<T>::SetRealPart( Int i, Int j, Base<T> alpha ) */ \
  ElError ElMatrixSetRealPart_ ## SIG \
  ( ElMatrix_ ## SIG AHandle, ElInt i, ElInt j, Base<T> alpha ) \
  { EL_TRY( Reinterpret(AHandle)->SetRealPart(i,j,alpha) ) } \
  /* void Matrix<T>::SetImagPart( Int i, Int j, Base<T> alpha ) */ \
  ElError ElMatrixSetImagPart_ ## SIG \
  ( ElMatrix_ ## SIG AHandle, ElInt i, ElInt j, Base<T> alpha ) \
  { EL_TRY( Reinterpret(AHandle)->SetImagPart(i,j,alpha) ) } \
  /* void Matrix<T>::UpdateRealPart( Int i, Int j, Base<T> alpha ) */ \
  ElError ElMatrixUpdateRealPart_ ## SIG \
  ( ElMatrix_ ## SIG AHandle, ElInt i, ElInt j, Base<T> alpha ) \
  { EL_TRY( Reinterpret(AHandle)->UpdateRealPart(i,j,alpha) ) } \
  /* void Matrix<T>::UpdateImagPart( Int i, Int j, Base<T> alpha ) */ \
  ElError ElMatrixUpdateImagPart_ ## SIG \
  ( ElMatrix_ ## SIG AHandle, ElInt i, ElInt j, Base<T> alpha ) \
  { EL_TRY( Reinterpret(AHandle)->UpdateImagPart(i,j,alpha) ) } \
  /* void Matrix<T>::MakeReal( Int i, Int j ) */ \
  ElError ElMatrixMakeReal_ ## SIG \
  ( ElMatrix_ ## SIG AHandle, ElInt i, ElInt j ) \
  { EL_TRY( Reinterpret(AHandle)->MakeReal(i,j) ) } \
  /* void Matrix<T>::Conjugate( Int i, Int j ) */ \
  ElError ElMatrixConjugate_ ## SIG \
  ( ElMatrix_ ## SIG AHandle, ElInt i, ElInt j ) \
  { EL_TRY( Reinterpret(AHandle)->Conjugate(i,j) ) }

#define MATRIX_DIAGONAL(SIG,T) \
  /* Matrix<T> Matrix<T>::GetDiagonal( Int offset ) const */ \
  ElError ElMatrixGetDiagonal_ ## SIG \
  ( ElConstMatrix_ ## SIG AHandle, ElInt offset, ElMatrix_ ## SIG *dHandle ) \
  { EL_TRY( auto d = new Matrix<T>; \
            Reinterpret(AHandle)->GetDiagonal( *d, offset ); \
            *dHandle = Reinterpret(d) ) } \
  /* void Matrix<T>::SetDiagonal( const Matrix<T>& d, Int offset ) */ \
  ElError ElMatrixSetDiagonal_ ## SIG \
  ( ElMatrix_ ## SIG AHandle, ElConstMatrix_ ## SIG dHandle, ElInt offset ) \
  { EL_TRY \
    ( Reinterpret(AHandle)->SetDiagonal( Reinterpret(dHandle), offset ) ) } \
  /* void Matrix<T>::UpdateDiagonal( const Matrix<T>& d, Int offset ) */ \
  ElError ElMatrixUpdateDiagonal_ ## SIG \
  ( ElMatrix_ ## SIG AHandle, ElConstMatrix_ ## SIG dHandle, ElInt offset ) \
  { EL_TRY \
    ( Reinterpret(AHandle)->UpdateDiagonal( Reinterpret(dHandle), offset ) ) }

#define MATRIX_DIAGONAL_COMPLEX(SIG,SIGBASE,T) \
  /* Matrix<Base<T>> Matrix<T>::GetRealPartOfDiagonal( Int offset ) const */ \
  ElError ElMatrixGetRealPartOfDiagonal_ ## SIG \
  ( ElConstMatrix_ ## SIG AHandle, ElInt offset, \
    ElMatrix_ ## SIGBASE *dHandle ) \
  { EL_TRY( auto d = new Matrix<Base<T>>; \
            Reinterpret(AHandle)->GetRealPartOfDiagonal( *d, offset ); \
            *dHandle = Reinterpret(d) ) } \
  /* Matrix<Base<T>> Matrix<T>::GetImagPartOfDiagonal( Int offset ) const */ \
  ElError ElMatrixGetImagPartOfDiagonal_ ## SIG \
  ( ElConstMatrix_ ## SIG AHandle, ElInt offset, \
    ElMatrix_ ## SIGBASE *dHandle ) \
  { EL_TRY( auto d = new Matrix<Base<T>>; \
            Reinterpret(AHandle)->GetImagPartOfDiagonal( *d, offset ); \
            *dHandle = Reinterpret(d) ) } \
  /* void Matrix<T>::SetRealPartOfDiagonal \
     ( const Matrix<Base<T>>& d, Int offset ) */ \
  ElError ElMatrixSetRealPartOfDiagonal_ ## SIG \
  ( ElMatrix_ ## SIG AHandle, \
    ElConstMatrix_ ## SIGBASE dHandle, ElInt offset ) \
  { EL_TRY( Reinterpret(AHandle)->SetRealPartOfDiagonal \
            (Reinterpret(dHandle),offset) ) } \
  /* void Matrix<T>::SetImagPartOfDiagonal \
     ( const Matrix<Base<T>>& d, Int offset ) */ \
  ElError ElMatrixSetImagPartOfDiagonal_ ## SIG \
  ( ElMatrix_ ## SIG AHandle, \
    ElConstMatrix_ ## SIGBASE dHandle, ElInt offset ) \
  { EL_TRY( Reinterpret(AHandle)->SetImagPartOfDiagonal \
            (Reinterpret(dHandle),offset) ) } \
  /* void Matrix<T>::UpdateRealPartOfDiagonal \
     ( const Matrix<Base<T>>& d, Int offset ) */ \
  ElError ElMatrixUpdateRealPartOfDiagonal_ ## SIG \
  ( ElMatrix_ ## SIG AHandle, \
    ElConstMatrix_ ## SIGBASE dHandle, ElInt offset ) \
  { EL_TRY( Reinterpret(AHandle)->UpdateRealPartOfDiagonal \
            (Reinterpret(dHandle),offset) ) } \
  /* void Matrix<T>::UpdateImagPartOfDiagonal \
     ( const Matrix<Base<T>>& d, Int offset ) */ \
  ElError ElMatrixUpdateImagPartOfDiagonal_ ## SIG \
  ( ElMatrix_ ## SIG AHandle, \
    ElConstMatrix_ ## SIGBASE dHandle, ElInt offset ) \
  { EL_TRY( Reinterpret(AHandle)->UpdateImagPartOfDiagonal \
            (Reinterpret(dHandle),offset) ) } \
  /* void Matrix<T>::MakeDiagonalReal( Int offset ) */ \
  ElError ElMatrixMakeDiagonalReal_ ## SIG \
  ( ElMatrix_ ## SIG AHandle, ElInt offset ) \
  { EL_TRY( Reinterpret(AHandle)->MakeDiagonalReal(offset) ) } \
  /* void Matrix<T>::ConjugateDiagonal( Int offset ) */ \
  ElError ElMatrixConjugateDiagonal_ ## SIG \
  ( ElMatrix_ ## SIG AHandle, ElInt offset ) \
  { EL_TRY( Reinterpret(AHandle)->ConjugateDiagonal(offset) ) }

#define MATRIX_SUBMATRIX(SIG,T) \
  /* Matrix<T> Matrix<T>::GetSubmatrix
     ( const std::vector<Int>& rowInds, \
       const std::vector<Int>& colInds ) const */ \
  ElError ElMatrixGetSubmatrix_ ## SIG \
  ( ElConstMatrix_ ## SIG AHandle, \
    ElInt numRowInds, const ElInt* rowInds, \
    ElInt numColInds, const ElInt* colInds, ElMatrix_ ## SIG *ASubHandle ) \
  { EL_TRY( std::vector<Int> rowIndVec( rowInds, rowInds+numRowInds ); \
            std::vector<Int> colIndVec( colInds, colInds+numColInds ); \
            auto ASub = new Matrix<T>; \
            Reinterpret(AHandle)->GetSubmatrix( rowIndVec, colIndVec, *ASub ); \
            *ASubHandle = Reinterpret(ASub) ) }

#define C_PROTO(SIG,T) \
  MATRIX_CONSTRUCT(SIG,T) \
  MATRIX_RECONFIG(SIG,T) \
  MATRIX_BASIC(SIG,T) \
  MATRIX_SINGLEENTRY(SIG,T) \
  MATRIX_DIAGONAL(SIG,T) \
  MATRIX_SUBMATRIX(SIG,T)

#define C_PROTO_COMPLEX(SIG,SIGBASE,T) \
  C_PROTO(SIG,T) \
  MATRIX_SINGLEENTRY_COMPLEX(SIG,SIGBASE,T) \
  MATRIX_DIAGONAL_COMPLEX(SIG,SIGBASE,T)

#include "El/macros/CInstantiate.h"

} // extern "C"
