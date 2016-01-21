/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
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

EL_EXPORT ElError ElQtImageFormat
( ElFileFormat format, const char** formatStr );
EL_EXPORT ElError ElFileExtension
( ElFileFormat format, const char** fileExt );
EL_EXPORT ElError ElFormatFromExtension
( const char* ext, ElFileFormat* format );
EL_EXPORT ElError ElDetectFormat
( const char* filename, ElFileFormat* format );

/* TODO: FileSize wrapper. What should the input/output types be? */

typedef enum {
  EL_GRAYSCALE,
  EL_GRAYSCALE_DISCRETE,
  EL_RED_BLACK_GREEN,
  EL_BLUE_RED
} ElColorMap;

/* Color maps
   ========== */
EL_EXPORT ElError ElSetColorMap( ElColorMap map );
EL_EXPORT ElError ElGetColorMap( ElColorMap* map );
EL_EXPORT ElError ElSetNumDiscreteColors( ElInt numColors );
EL_EXPORT ElError ElNumDiscreteColors( ElInt* numColors );

/* Display
   ======= */
EL_EXPORT ElError ElProcessEvents( int numMsecs );

EL_EXPORT ElError ElDisplay_i( ElConstMatrix_i A, const char* title );
EL_EXPORT ElError ElDisplay_s( ElConstMatrix_s A, const char* title );
EL_EXPORT ElError ElDisplay_d( ElConstMatrix_d A, const char* title );
EL_EXPORT ElError ElDisplay_c( ElConstMatrix_c A, const char* title );
EL_EXPORT ElError ElDisplay_z( ElConstMatrix_z A, const char* title );

EL_EXPORT ElError ElDisplayDist_i( ElConstDistMatrix_i A, const char* title );
EL_EXPORT ElError ElDisplayDist_s( ElConstDistMatrix_s A, const char* title );
EL_EXPORT ElError ElDisplayDist_d( ElConstDistMatrix_d A, const char* title );
EL_EXPORT ElError ElDisplayDist_c( ElConstDistMatrix_c A, const char* title );
EL_EXPORT ElError ElDisplayDist_z( ElConstDistMatrix_z A, const char* title );

EL_EXPORT ElError ElDisplayDistMultiVec_i
( ElConstDistMultiVec_i X, const char* title );
EL_EXPORT ElError ElDisplayDistMultiVec_s
( ElConstDistMultiVec_s X, const char* title );
EL_EXPORT ElError ElDisplayDistMultiVec_d
( ElConstDistMultiVec_d X, const char* title );
EL_EXPORT ElError ElDisplayDistMultiVec_c
( ElConstDistMultiVec_c X, const char* title );
EL_EXPORT ElError ElDisplayDistMultiVec_z
( ElConstDistMultiVec_z X, const char* title );

EL_EXPORT ElError ElDisplayGraph( ElConstGraph graph, const char* title );
EL_EXPORT ElError ElDisplayDistGraph
( ElConstDistGraph graph, const char* title );

EL_EXPORT ElError ElDisplaySparse_i
( ElConstSparseMatrix_i A, const char* title );
EL_EXPORT ElError ElDisplaySparse_s
( ElConstSparseMatrix_s A, const char* title );
EL_EXPORT ElError ElDisplaySparse_d
( ElConstSparseMatrix_d A, const char* title );
EL_EXPORT ElError ElDisplaySparse_c
( ElConstSparseMatrix_c A, const char* title );
EL_EXPORT ElError ElDisplaySparse_z
( ElConstSparseMatrix_z A, const char* title );

EL_EXPORT ElError ElDisplayDistSparse_i
( ElConstDistSparseMatrix_i A, const char* title );
EL_EXPORT ElError ElDisplayDistSparse_s
( ElConstDistSparseMatrix_s A, const char* title );
EL_EXPORT ElError ElDisplayDistSparse_d
( ElConstDistSparseMatrix_d A, const char* title );
EL_EXPORT ElError ElDisplayDistSparse_c
( ElConstDistSparseMatrix_c A, const char* title );
EL_EXPORT ElError ElDisplayDistSparse_z
( ElConstDistSparseMatrix_z A, const char* title );

/* Print 
   ===== */
/* 
  TODO: Extend to support more general output streams. This is a bit 
        problematic due to difficulties in converting between FILE* and
        std::ostream 
*/
EL_EXPORT ElError ElPrint_i( ElConstMatrix_i A, const char* title );
EL_EXPORT ElError ElPrint_s( ElConstMatrix_s A, const char* title );
EL_EXPORT ElError ElPrint_d( ElConstMatrix_d A, const char* title );
EL_EXPORT ElError ElPrint_c( ElConstMatrix_c A, const char* title );
EL_EXPORT ElError ElPrint_z( ElConstMatrix_z A, const char* title );

EL_EXPORT ElError ElPrintDist_i( ElConstDistMatrix_i A, const char* title );
EL_EXPORT ElError ElPrintDist_s( ElConstDistMatrix_s A, const char* title );
EL_EXPORT ElError ElPrintDist_d( ElConstDistMatrix_d A, const char* title );
EL_EXPORT ElError ElPrintDist_c( ElConstDistMatrix_c A, const char* title );
EL_EXPORT ElError ElPrintDist_z( ElConstDistMatrix_z A, const char* title );

EL_EXPORT ElError ElPrintDistMultiVec_i
( ElConstDistMultiVec_i X, const char* title );
EL_EXPORT ElError ElPrintDistMultiVec_s
( ElConstDistMultiVec_s X, const char* title );
EL_EXPORT ElError ElPrintDistMultiVec_d
( ElConstDistMultiVec_d X, const char* title );
EL_EXPORT ElError ElPrintDistMultiVec_c
( ElConstDistMultiVec_c X, const char* title );
EL_EXPORT ElError ElPrintDistMultiVec_z
( ElConstDistMultiVec_z X, const char* title );

EL_EXPORT ElError ElPrintGraph( ElConstGraph graph, const char* title );
EL_EXPORT ElError ElPrintDistGraph( ElConstDistGraph graph, const char* title );

EL_EXPORT ElError ElPrintSparse_i( ElConstSparseMatrix_i A, const char* title );
EL_EXPORT ElError ElPrintSparse_s( ElConstSparseMatrix_s A, const char* title );
EL_EXPORT ElError ElPrintSparse_d( ElConstSparseMatrix_d A, const char* title );
EL_EXPORT ElError ElPrintSparse_c( ElConstSparseMatrix_c A, const char* title );
EL_EXPORT ElError ElPrintSparse_z( ElConstSparseMatrix_z A, const char* title );

EL_EXPORT ElError ElPrintDistSparse_i
( ElConstDistSparseMatrix_i A, const char* title );
EL_EXPORT ElError ElPrintDistSparse_s
( ElConstDistSparseMatrix_s A, const char* title );
EL_EXPORT ElError ElPrintDistSparse_d
( ElConstDistSparseMatrix_d A, const char* title );
EL_EXPORT ElError ElPrintDistSparse_c
( ElConstDistSparseMatrix_c A, const char* title );
EL_EXPORT ElError ElPrintDistSparse_z
( ElConstDistSparseMatrix_z A, const char* title );

/* Read
   ==== */
EL_EXPORT ElError ElRead_i
( ElMatrix_i A, const char* filename, ElFileFormat format );
EL_EXPORT ElError ElRead_s
( ElMatrix_s A, const char* filename, ElFileFormat format );
EL_EXPORT ElError ElRead_d
( ElMatrix_d A, const char* filename, ElFileFormat format );
EL_EXPORT ElError ElRead_c
( ElMatrix_c A, const char* filename, ElFileFormat format );
EL_EXPORT ElError ElRead_z
( ElMatrix_z A, const char* filename, ElFileFormat format );

EL_EXPORT ElError ElReadDist_i
( ElDistMatrix_i A, const char* filename, ElFileFormat format, 
  bool sequential );
EL_EXPORT ElError ElReadDist_s
( ElDistMatrix_s A, const char* filename, ElFileFormat format, 
  bool sequential );
EL_EXPORT ElError ElReadDist_d
( ElDistMatrix_d A, const char* filename, ElFileFormat format, 
  bool sequential );
EL_EXPORT ElError ElReadDist_c
( ElDistMatrix_c A, const char* filename, ElFileFormat format, 
  bool sequential );
EL_EXPORT ElError ElReadDist_z
( ElDistMatrix_z A, const char* filename, ElFileFormat format, 
  bool sequential );

/* Spy
   === */
EL_EXPORT ElError ElSpy_i( ElConstMatrix_i A, const char* title, ElInt tol );
EL_EXPORT ElError ElSpy_s( ElConstMatrix_s A, const char* title, float tol );
EL_EXPORT ElError ElSpy_d( ElConstMatrix_d A, const char* title, double tol );
EL_EXPORT ElError ElSpy_c( ElConstMatrix_c A, const char* title, float tol );
EL_EXPORT ElError ElSpy_z( ElConstMatrix_z A, const char* title, double tol );

EL_EXPORT ElError ElSpyDist_i
( ElConstDistMatrix_i A, const char* title, ElInt tol );
EL_EXPORT ElError ElSpyDist_s
( ElConstDistMatrix_s A, const char* title, float tol );
EL_EXPORT ElError ElSpyDist_d
( ElConstDistMatrix_d A, const char* title, double tol );
EL_EXPORT ElError ElSpyDist_c
( ElConstDistMatrix_c A, const char* title, float tol );
EL_EXPORT ElError ElSpyDist_z
( ElConstDistMatrix_z A, const char* title, double tol );

/* Write
   ===== */
EL_EXPORT ElError ElWrite_i
( ElConstMatrix_i A, const char* basename, ElFileFormat format, 
  const char* title );
EL_EXPORT ElError ElWrite_s
( ElConstMatrix_s A, const char* basename, ElFileFormat format, 
  const char* title );
EL_EXPORT ElError ElWrite_d
( ElConstMatrix_d A, const char* basename, ElFileFormat format, 
  const char* title );
EL_EXPORT ElError ElWrite_c
( ElConstMatrix_c A, const char* basename, ElFileFormat format, 
  const char* title );
EL_EXPORT ElError ElWrite_z
( ElConstMatrix_z A, const char* basename, ElFileFormat format, 
  const char* title );

EL_EXPORT ElError ElWriteDist_i
( ElConstDistMatrix_i A, const char* basename, ElFileFormat format, 
  const char* title );
EL_EXPORT ElError ElWriteDist_s
( ElConstDistMatrix_s A, const char* basename, ElFileFormat format, 
  const char* title );
EL_EXPORT ElError ElWriteDist_d
( ElConstDistMatrix_d A, const char* basename, ElFileFormat format, 
  const char* title );
EL_EXPORT ElError ElWriteDist_c
( ElConstDistMatrix_c A, const char* basename, ElFileFormat format, 
  const char* title );
EL_EXPORT ElError ElWriteDist_z
( ElConstDistMatrix_z A, const char* basename, ElFileFormat format, 
  const char* title );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_IO_C_H */
