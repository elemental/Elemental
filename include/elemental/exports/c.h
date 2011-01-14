/*
   Copyright (c) 2009-2011, Jack Poulson
   All rights reserved.

   This file is part of Elemental.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

    - Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

    - Neither the name of the owner nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/
#ifndef ELEMENTAL_EXPORTS_C_H
#define ELEMENTAL_EXPORTS_C_H 1

#include "mpi.h"

/* Make sure that all of our configuration definitions are pulled in */
#include "elemental/config.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifndef WITHOUT_COMPLEX
/* We should not assume C99 support and use <complex.h> */
typedef struct { float real; float imag; } SComplex;
typedef struct { double real; double imag; } DComplex;
#endif /* WITHOUT_COMPLEX */

/* Grid A */
typedef int Grid;

/* Handles for real distributed matrices */
typedef int MC_MR_Single;
typedef int MC_MR_Double;
typedef int MC_Star_Single;
typedef int MC_Star_Double;
typedef int MD_Star_Single;
typedef int MD_Star_Double;
typedef int MR_Star_Single;
typedef int MR_Star_Double;
typedef int Star_MC_Single;
typedef int Star_MC_Double;
typedef int Star_MD_Single;
typedef int Star_MD_Double;
typedef int Star_MR_Single;
typedef int Star_MR_Double;
typedef int Star_Star_Single;
typedef int Star_Star_Double;
typedef int Star_VC_Single;
typedef int Star_VC_Double;
typedef int Star_VR_Single;
typedef int Star_VR_Double;
typedef int VC_Star_Single;
typedef int VC_Star_Double;
typedef int VR_Star_Single;
typedef int VR_Star_Double;

#ifndef WITHOUT_COMPLEX
/* Handles for complex distributed matrices */
typedef int MC_MR_SComplex;
typedef int MC_MR_DComplex;
typedef int MC_Star_SComplex;
typedef int MC_Star_DComplex;
typedef int MD_Star_SComplex;
typedef int MD_Star_DComplex;
typedef int MR_Star_SComplex;
typedef int MR_Star_DComplex;
typedef int Star_MC_SComplex;
typedef int Star_MC_DComplex;
typedef int Star_MD_SComplex;
typedef int Star_MD_DComplex;
typedef int Star_MR_SComplex;
typedef int Star_MR_DComplex;
typedef int Star_Star_SComplex;
typedef int Star_Star_DComplex;
typedef int Star_VC_SComplex;
typedef int Star_VC_DComplex;
typedef int Star_VR_SComplex;
typedef int Star_VR_DComplex;
typedef int VC_Star_SComplex;
typedef int VC_Star_DComplex;
typedef int VR_Star_SComplex;
typedef int VR_Star_DComplex;
#endif /* WITHOUT_COMPLEX */

/* Elemental's environment */
void ElementalInit( int* argc, char** argv[] );
void ElementalFinalize();
int ElementalBlocksize();
void ElementalSetBlocksize( int blocksize );
void ElementalPushBlocksizeStack( int blocksize );
void ElementalPopBlocksizeStack();

/* Grid manipulation */
Grid ElementalDefaultGrid( MPI_Comm comm );
Grid ElementalGrid( MPI_Comm comm, int r, int c );
Grid ElementalGridHeight( Grid g );
Grid ElementalGridWidth( Grid g );
Grid ElementalGridSize( Grid g );
Grid ElementalInGrid( Grid g );
Grid ElementalGridVCRank( Grid g );
Grid ElementalGridVRRank( Grid g );
Grid ElementalGridMCRank( Grid g );
Grid ElementalGridMRRank( Grid g );
MPI_Comm ElementalGridVCComm( Grid g );
MPI_Comm ElementalGridVRComm( Grid g );
MPI_Comm ElementalGridMCComm( Grid g );
MPI_Comm ElementalGridMRComm( Grid g );

/* Clean up */
void ElementalClearGrids();
void ElementalClearDistMatrices();

/* Real double-precision distributed matrices */
MC_MR_Single ElementalCreateEmpty_MC_MR_Single( Grid g );
MC_MR_Double ElementalCreateEmpty_MC_MR_Double( Grid g );
Star_VR_Single ElementalCreateEmpty_Star_VR_Single( Grid g );
Star_VR_Double ElementalCreateEmpty_Star_VR_Double( Grid g );

MC_MR_Single ElementalRegister_MC_MR_Single
( int height, int width, int colAlignment, int rowAlignment,
  float* buffer, int ldim, Grid g );
MC_MR_Double ElementalRegister_MC_MR_Double
( int height, int width, int colAlignment, int rowAlignment,
  double* buffer, int ldim, Grid g );
Star_VR_Single ElementalRegister_Star_VR_Single
( int height, int width, int rowAlignment, 
  float* buffer, int ldim, Grid g);
Star_VR_Double ElementalRegister_Star_VR_Double
( int height, int width, int rowAlignment, 
  double* buffer, int ldim, Grid g );

void ElementalPrint_MC_MR_Single( char* msg, MC_MR_Single A );
void ElementalPrint_MC_MR_Double( char* msg, MC_MR_Double A );
void ElementalPrint_MC_Star_Single( char* msg, MC_Star_Single A );
void ElementalPrint_MC_Star_Double( char* msg, MC_Star_Double A );
void ElementalPrint_MD_Star_Single( char* msg, MD_Star_Single A );
void ElementalPrint_MD_Star_Double( char* msg, MD_Star_Double A );
void ElementalPrint_MR_Star_Single( char* msg, MR_Star_Single A );
void ElementalPrint_MR_Star_Double( char* msg, MR_Star_Double A );
void ElementalPrint_Star_MC_Single( char* msg, Star_MC_Single A );
void ElementalPrint_Star_MC_Double( char* msg, Star_MC_Double A );
void ElementalPrint_Star_MR_Single( char* msg, Star_MR_Single A );
void ElementalPrint_Star_MR_Double( char* msg, Star_MR_Double A );
void ElementalPrint_Star_Star_Single( char* msg, Star_Star_Single A );
void ElementalPrint_Star_Star_Double( char* msg, Star_Star_Double A );
void ElementalPrint_Star_VC_Single( char* msg, Star_VC_Single A );
void ElementalPrint_Star_VC_Double( char* msg, Star_VC_Double A );
void ElementalPrint_Star_VR_Single( char* msg, Star_VR_Single A );
void ElementalPrint_Star_VR_Double( char* msg, Star_VR_Double A );
void ElementalPrint_VC_Star_Single( char* msg, VC_Star_Single A );
void ElementalPrint_VC_Star_Double( char* msg, VC_Star_Double A );
void ElementalPrint_VR_Star_Single( char* msg, VR_Star_Single A );
void ElementalPrint_VR_Star_Double( char* msg, VR_Star_Double A );

#ifndef WITHOUT_COMPLEX
/* Complex double-precision distributed matrices */
MC_MR_SComplex ElementalCreateEmpty_MC_MR_SComplex( Grid g );
MC_MR_DComplex ElementalCreateEmpty_MC_MR_DComplex( Grid g );
Star_VR_SComplex ElementalCreateEmpty_Star_VR_SComplex( Grid g );
Star_VR_DComplex ElementalCreateEmpty_Star_VR_DComplex( Grid g );

MC_MR_SComplex ElementalRegister_MC_MR_SComplex
( int height, int width, int colAlignment, int rowAlignment,
  SComplex* buffer, int ldim, Grid g );
MC_MR_DComplex ElementalRegister_MC_MR_DComplex
( int height, int width, int colAlignment, int rowAlignment,
  DComplex* buffer, int ldim, Grid g );  

Star_VR_SComplex ElementalRegister_Star_VR_SComplex
( int height, int width, int rowAlignment,
  SComplex* buffer, int ldim, Grid g );
Star_VR_DComplex ElementalRegister_Star_VR_DComplex
( int height, int width, int rowAlignment, 
  DComplex* buffer, int ldim, Grid g );

void ElementalPrint_MC_MR_SComplex( char* msg, MC_MR_SComplex A );
void ElementalPrint_MC_MR_DComplex( char* msg, MC_MR_DComplex A );
void ElementalPrint_MC_Star_SComplex( char* msg, MC_Star_SComplex A );
void ElementalPrint_MC_Star_DComplex( char* msg, MC_Star_DComplex A );
void ElementalPrint_MD_Star_SComplex( char* msg, MD_Star_SComplex A );
void ElementalPrint_MD_Star_DComplex( char* msg, MD_Star_DComplex A );
void ElementalPrint_MR_Star_SComplex( char* msg, MR_Star_SComplex A );
void ElementalPrint_MR_Star_DComplex( char* msg, MR_Star_DComplex A );
void ElementalPrint_Star_MC_SComplex( char* msg, Star_MC_SComplex A );
void ElementalPrint_Star_MC_DComplex( char* msg, Star_MC_DComplex A );
void ElementalPrint_Star_MR_SComplex( char* msg, Star_MR_SComplex A );
void ElementalPrint_Star_MR_DComplex( char* msg, Star_MR_DComplex A );
void ElementalPrint_Star_Star_SComplex( char* msg, Star_Star_SComplex A );
void ElementalPrint_Star_Star_DComplex( char* msg, Star_Star_DComplex A );
void ElementalPrint_Star_VC_SComplex( char* msg, Star_VC_SComplex A );
void ElementalPrint_Star_VC_DComplex( char* msg, Star_VC_DComplex A );
void ElementalPrint_Star_VR_SComplex( char* msg, Star_VR_SComplex A );
void ElementalPrint_Star_VR_DComplex( char* msg, Star_VR_DComplex A );
void ElementalPrint_VC_Star_SComplex( char* msg, VC_Star_SComplex A );
void ElementalPrint_VC_Star_DComplex( char* msg, VC_Star_DComplex A );
void ElementalPrint_VR_Star_SComplex( char* msg, VR_Star_SComplex A );
void ElementalPrint_VR_Star_DComplex( char* msg, VR_Star_DComplex A );
#endif /* WITHOUT_COMPLEX */

/* Utilities */
int ElementalLocalLength
( int globalLength, int myIndex, int alignment, int modulus );

/* LAPACK-level interface */
#ifndef WITHOUT_PMRRR
void
ElementalHermitianEigDouble
( char uplo, MC_MR_Double A, Star_VR_Double w, MC_MR_Double Z,
  int tryForHighAccuracy );
void
ElementalHermitianEigDComplex
( char uplo, MC_MR_DComplex A, Star_VR_Double w, MC_MR_DComplex Z, 
  int tryForHighAccuracy );
#endif /* WITHOUT_PMRRR */

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* ELEMENTAL_EXPORTS_C_H */

