/*
   Copyright (c) 2009-2014, Jack Poulson
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
{ EL_TRY( Display( *CReflect(graph), std::string(title) ) ) }

/* Print
   ===== */
ElError ElPrintGraph( ElConstGraph graph, const char* title )
{ EL_TRY( Print( *CReflect(graph), std::string(title) ) ) }

#define C_PROTO(SIG,SIGBASE,T) \
  /* Display */ \
  ElError ElDisplay_ ## SIG \
  ( ElConstMatrix_ ## SIG A, const char* title ) \
  { EL_TRY( Display( *CReflect(A), std::string(title) ) ) } \
  ElError ElDisplayDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, const char* title ) \
  { EL_TRY( Display( *CReflect(A), std::string(title) ) ) } \
  /* Print */ \
  ElError ElPrint_ ## SIG \
  ( ElConstMatrix_ ## SIG A, const char* title ) \
  { EL_TRY( Print( *CReflect(A), std::string(title) ) ) } \
  ElError ElPrintDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, const char* title ) \
  { EL_TRY( Print( *CReflect(A), std::string(title) ) ) } \
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
  { EL_TRY( Spy( *CReflect(A), std::string(title), tol ) ) } \
  ElError ElSpyDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, const char* title, Base<T> tol ) \
  { EL_TRY( Spy( *CReflect(A), std::string(title), tol ) ) } \
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
