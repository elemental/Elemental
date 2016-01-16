/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
#include "El.h"
using namespace El;

extern "C" {

ElError ElQtImageFormat( ElFileFormat format, const char** formatStr )
{ EL_TRY( *formatStr = QtImageFormat(CReflect(format)) ) }

ElError ElFileExtension( ElFileFormat format, const char** fileExt )
{ EL_TRY( *fileExt = CReflect(FileExtension(CReflect(format))) ) }

ElError ElFormatFromExtension( const char* ext, ElFileFormat* format )
{ EL_TRY( *format = CReflect(FormatFromExtension(CReflect(ext))) ) }

ElError ElDetectFormat( const char* filename, ElFileFormat* format )
{ EL_TRY( *format = CReflect(DetectFormat(CReflect(filename))) ) }

/* Color maps
   ========== */
ElError ElSetColorMap( ElColorMap map )
{ EL_TRY( SetColorMap(CReflect(map)) ) }

ElError ElGetColorMap( ElColorMap* map )
{ EL_TRY( *map = CReflect(GetColorMap()) ) }

ElError ElSetNumDiscreteColors( ElInt numColors )
{ EL_TRY( SetNumDiscreteColors(numColors) ) }

ElError ElNumDiscreteColors( ElInt* numColors )
{ EL_TRY( *numColors = NumDiscreteColors() ) }

/* Display
   ======= */
ElError ElProcessEvents( int numMsecs )
{ EL_TRY( ProcessEvents(numMsecs) ) }

ElError ElDisplayGraph( ElConstGraph graph, const char* title )
{ EL_TRY( Display( *CReflect(graph), string(title) ) ) }
ElError ElDisplayDistGraph( ElConstDistGraph graph, const char* title )
{ EL_TRY( Display( *CReflect(graph), string(title) ) ) }

/* Print
   ===== */
ElError ElPrintGraph( ElConstGraph graph, const char* title )
{ EL_TRY( Print( *CReflect(graph), string(title) ) ) }
ElError ElPrintDistGraph( ElConstDistGraph graph, const char* title )
{ EL_TRY( Print( *CReflect(graph), string(title) ) ) }

#define C_PROTO(SIG,SIGBASE,T) \
  /* Display */ \
  ElError ElDisplay_ ## SIG \
  ( ElConstMatrix_ ## SIG A, const char* title ) \
  { EL_TRY( Display( *CReflect(A), string(title) ) ) } \
  ElError ElDisplayDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, const char* title ) \
  { EL_TRY( Display( *CReflect(A), string(title) ) ) } \
  ElError ElDisplayDistMultiVec_ ## SIG \
  ( ElConstDistMultiVec_ ## SIG X, const char* title ) \
  { EL_TRY( Display( *CReflect(X), string(title) ) ) } \
  ElError ElDisplaySparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, const char* title ) \
  { EL_TRY( Display( *CReflect(A), string(title) ) ) } \
  ElError ElDisplayDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, const char* title ) \
  { EL_TRY( Display( *CReflect(A), string(title) ) ) } \
  /* Print */ \
  ElError ElPrint_ ## SIG \
  ( ElConstMatrix_ ## SIG A, const char* title ) \
  { EL_TRY( Print( *CReflect(A), string(title) ) ) } \
  ElError ElPrintDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, const char* title ) \
  { EL_TRY( Print( *CReflect(A), string(title) ) ) } \
  ElError ElPrintDistMultiVec_ ## SIG \
  ( ElConstDistMultiVec_ ## SIG X, const char* title ) \
  { EL_TRY( Print( *CReflect(X), string(title) ) ) } \
  ElError ElPrintSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, const char* title ) \
  { EL_TRY( Print( *CReflect(A), string(title) ) ) } \
  ElError ElPrintDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, const char* title ) \
  { EL_TRY( Print( *CReflect(A), string(title) ) ) } \
  /* Read */ \
  ElError ElRead_ ## SIG \
  ( ElMatrix_ ## SIG A, const char* filename, ElFileFormat format ) \
  { EL_TRY( Read( *CReflect(A), CReflect(filename), CReflect(format) ) ) } \
  ElError ElReadDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, const char* filename, ElFileFormat format, \
    bool sequential ) \
  { EL_TRY( Read( *CReflect(A), CReflect(filename), CReflect(format), \
                  sequential ) ) } \
  /* Spy */ \
  ElError ElSpy_ ## SIG \
  ( ElConstMatrix_ ## SIG A, const char* title, Base<T> tol ) \
  { EL_TRY( Spy( *CReflect(A), string(title), tol ) ) } \
  ElError ElSpyDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, const char* title, Base<T> tol ) \
  { EL_TRY( Spy( *CReflect(A), string(title), tol ) ) } \
  /* Write */ \
  ElError ElWrite_ ## SIG \
  ( ElConstMatrix_ ## SIG A, const char* basename, ElFileFormat format, \
    const char* title ) \
  { EL_TRY( Write( \
      *CReflect(A), CReflect(basename), CReflect(format), \
      CReflect(title) ) ) } \
  ElError ElWriteDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, const char* basename, ElFileFormat format, \
    const char* title ) \
  { EL_TRY( Write( \
      *CReflect(A), CReflect(basename), CReflect(format), \
      CReflect(title) ) ) }

#include "El/macros/CInstantiate.h"

} // extern "C"
