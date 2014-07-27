/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_IO_C_H
#define EL_IO_C_H

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
  EL_AUTO,
  EL_ASCII,
  EL_ASCII_MATLAB,
  EL_BINARY,
  EL_BINARY_FLAT,
  EL_BMP,
  EL_JPG,
  EL_JPEG,
  EL_MATRIX_MARKET,
  EL_PNG,
  EL_PPM,
  EL_XBM,
  EL_XPM,
  EL_FileFormat_MAX
} ElFileFormat;

ElError ElQtImageFormat( ElFileFormat format, const char** formatStr );
ElError ElFileExtension( ElFileFormat format, const char** fileExt );
ElError ElFormatFromExtension( const char* ext, ElFileFormat* format );
ElError ElDetectFormat( const char* filename, ElFileFormat* format );

/* TODO: FileSize wrapper. What should the input/output types be? */

typedef enum {
  EL_GRAYSCALE,
  EL_GRAYSCALE_DISCRETE,
  EL_RED_BLACK_GREEN,
  EL_BLUE_RED
} ElColorMap;

/* Color maps
   ========== */
ElError ElSetColorMap( ElColorMap map );
ElError ElGetColorMap( ElColorMap* map );
ElError ElSetNumDiscreteColors( ElInt numColors );
ElError ElNumDiscreteColors( ElInt* numColors );

/* Display
   ======= */
ElError ElDisplay_i( ElConstMatrix_i A, const char* title );
ElError ElDisplay_s( ElConstMatrix_s A, const char* title );
ElError ElDisplay_d( ElConstMatrix_d A, const char* title );
ElError ElDisplay_c( ElConstMatrix_c A, const char* title );
ElError ElDisplay_z( ElConstMatrix_z A, const char* title );

ElError ElDisplayDist_i( ElConstDistMatrix_i A, const char* title );
ElError ElDisplayDist_s( ElConstDistMatrix_s A, const char* title );
ElError ElDisplayDist_d( ElConstDistMatrix_d A, const char* title );
ElError ElDisplayDist_c( ElConstDistMatrix_c A, const char* title );
ElError ElDisplayDist_z( ElConstDistMatrix_z A, const char* title );

/* Print 
   ===== */
/* 
  TODO: Extend to support more general output streams. This is a bit 
        problematic due to difficulties in converting between FILE* and
        std::ostream 
*/
ElError ElPrint_i( ElConstMatrix_i A, const char* title );
ElError ElPrint_s( ElConstMatrix_s A, const char* title );
ElError ElPrint_d( ElConstMatrix_d A, const char* title );
ElError ElPrint_c( ElConstMatrix_c A, const char* title );
ElError ElPrint_z( ElConstMatrix_z A, const char* title );

ElError ElPrintDist_i( ElConstDistMatrix_i A, const char* title );
ElError ElPrintDist_s( ElConstDistMatrix_s A, const char* title );
ElError ElPrintDist_d( ElConstDistMatrix_d A, const char* title );
ElError ElPrintDist_c( ElConstDistMatrix_c A, const char* title );
ElError ElPrintDist_z( ElConstDistMatrix_z A, const char* title );

/* Read
   ==== */
ElError ElRead_i( ElMatrix_i A, const char* filename, ElFileFormat format );
ElError ElRead_s( ElMatrix_s A, const char* filename, ElFileFormat format );
ElError ElRead_d( ElMatrix_d A, const char* filename, ElFileFormat format );
ElError ElRead_c( ElMatrix_c A, const char* filename, ElFileFormat format );
ElError ElRead_z( ElMatrix_z A, const char* filename, ElFileFormat format );

ElError ElReadDist_i
( ElDistMatrix_i A, const char* filename, ElFileFormat format, 
  bool sequential );
ElError ElReadDist_s
( ElDistMatrix_s A, const char* filename, ElFileFormat format, 
  bool sequential );
ElError ElReadDist_d
( ElDistMatrix_d A, const char* filename, ElFileFormat format, 
  bool sequential );
ElError ElReadDist_c
( ElDistMatrix_c A, const char* filename, ElFileFormat format, 
  bool sequential );
ElError ElReadDist_z
( ElDistMatrix_z A, const char* filename, ElFileFormat format, 
  bool sequential );

/* Spy
   === */
ElError ElSpy_i( ElConstMatrix_i A, const char* title, ElInt tol );
ElError ElSpy_s( ElConstMatrix_s A, const char* title, float tol );
ElError ElSpy_d( ElConstMatrix_d A, const char* title, double tol );
ElError ElSpy_c( ElConstMatrix_c A, const char* title, float tol );
ElError ElSpy_z( ElConstMatrix_z A, const char* title, double tol );

ElError ElSpyDist_i( ElConstDistMatrix_i A, const char* title, ElInt tol );
ElError ElSpyDist_s( ElConstDistMatrix_s A, const char* title, float tol );
ElError ElSpyDist_d( ElConstDistMatrix_d A, const char* title, double tol );
ElError ElSpyDist_c( ElConstDistMatrix_c A, const char* title, float tol );
ElError ElSpyDist_z( ElConstDistMatrix_z A, const char* title, double tol );

/* Write
   ===== */
ElError ElWrite_i
( ElConstMatrix_i A, const char* basename, ElFileFormat format, 
  const char* title );
ElError ElWrite_s
( ElConstMatrix_s A, const char* basename, ElFileFormat format, 
  const char* title );
ElError ElWrite_d
( ElConstMatrix_d A, const char* basename, ElFileFormat format, 
  const char* title );
ElError ElWrite_c
( ElConstMatrix_c A, const char* basename, ElFileFormat format, 
  const char* title );
ElError ElWrite_z
( ElConstMatrix_z A, const char* basename, ElFileFormat format, 
  const char* title );

ElError ElWriteDist_i
( ElConstDistMatrix_i A, const char* basename, ElFileFormat format, 
  const char* title );
ElError ElWriteDist_s
( ElConstDistMatrix_s A, const char* basename, ElFileFormat format, 
  const char* title );
ElError ElWriteDist_d
( ElConstDistMatrix_d A, const char* basename, ElFileFormat format, 
  const char* title );
ElError ElWriteDist_c
( ElConstDistMatrix_c A, const char* basename, ElFileFormat format, 
  const char* title );
ElError ElWriteDist_z
( ElConstDistMatrix_z A, const char* basename, ElFileFormat format, 
  const char* title );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_IO_C_H */
