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

typedef MPI_Comm ElementalComm;
const ElementalComm ELEMENTAL_COMM_WORLD = MPI_COMM_WORLD;

/* Grid A */
typedef int Grid;

/* Handles for real distributed matrices */
typedef int MC_MR_Single;
typedef int MC_MR_Double;
typedef int MC_Star_Single;
typedef int MC_Star_Double;
typedef int MD_Star_Single;
typedef int MD_Star_Double;
typedef int MR_MC_Single;
typedef int MR_MC_Double;
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
typedef int MR_MC_SComplex;
typedef int MR_MC_DComplex;
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
int Blocksize();
void SetBlocksize( int blocksize );
void PushBlocksizeStack( int blocksize );
void PopBlocksizeStack();

/* Grid manipulation */
Grid CreateDefaultGrid( ElementalComm comm );
Grid CreateGrid( ElementalComm comm, int r, int c );
void DestroyGrid( Grid g );
Grid GridHeight( Grid g );
Grid GridWidth( Grid g );
Grid GridSize( Grid g );
Grid InGrid( Grid g );
Grid GridVCRank( Grid g );
Grid GridVRRank( Grid g );
Grid GridMCRank( Grid g );
Grid GridMRRank( Grid g );
ElementalComm GridVCComm( Grid g );
ElementalComm GridVRComm( Grid g );
ElementalComm GridMCComm( Grid g );
ElementalComm GridMRComm( Grid g );

/* Clean up */
void ClearGrids();
void ClearDistMatrices();

/*******************************************************************************
 * [MC,MR] manipulation routines                                               *
 ******************************************************************************/
/* Create empty */
MC_MR_Single CreateEmpty_MC_MR_Single( Grid g );
MC_MR_Double CreateEmpty_MC_MR_Double( Grid g );
#ifndef WITHOUT_COMPLEX
MC_MR_SComplex CreateEmpty_MC_MR_SComplex( Grid g );
MC_MR_DComplex CreateEmpty_MC_MR_DComplex( Grid g );
#endif /* WITHOUT_COMPLEX */
/* Create from exiting buffer */
MC_MR_Single Register_MC_MR_Single
( int height, int width, int colAlignment, int rowAlignment, 
  float* buffer, int ldim, Grid g );
MC_MR_Double Register_MC_MR_Double
( int height, int width, int colAlignment, int rowAlignment,
  double* buffer, int ldim, Grid g );
#ifndef WITHOUT_COMPLEX
MC_MR_SComplex Register_MC_MR_SComplex
( int height, int width, int colAlignment, int rowAlignment, 
  SComplex* buffer, int ldim, Grid g );
MC_MR_DComplex Register_MC_MR_DComplex
( int height, int width, int colAlignment, int rowAlignment,
  DComplex* buffer, int ldim, Grid g );
#endif /* WITHOUT_COMPLEX */
/* Print */
void Print_MC_MR_Single( char* msg, MC_MR_Single A );
void Print_MC_MR_Double( char* msg, MC_MR_Double A );
#ifndef WITHOUT_COMPLEX
void Print_MC_MR_SComplex( char* msg, MC_MR_SComplex A );
void Print_MC_MR_DComplex( char* msg, MC_MR_DComplex A );
#endif /* WITHOUT_COMPLEX */

/*******************************************************************************
 * [MC,* ] manipulation routines                                               *
 ******************************************************************************/
/* Create empty */
MC_Star_Single CreateEmpty_MC_Star_Single( Grid g );
MC_Star_Double CreateEmpty_MC_Star_Double( Grid g );
#ifndef WITHOUT_COMPLEX
MC_Star_SComplex CreateEmpty_MC_Star_SComplex( Grid g );
MC_Star_DComplex CreateEmpty_MC_Star_DComplex( Grid g );
#endif /* WITHOUT_COMPLEX */
/* Create from exiting buffer */
MC_Star_Single Register_MC_Star_Single
( int height, int width, int colAlignment, float* buffer, int ldim, Grid g );
MC_Star_Double Register_MC_Star_Double
( int height, int width, int colAlignment, double* buffer, int ldim, Grid g );
#ifndef WITHOUT_COMPLEX
MC_Star_SComplex Register_MC_Star_SComplex
( int height, int width, int colAlignment, SComplex* buffer, int ldim, Grid g );
MC_Star_DComplex Register_MC_Star_DComplex
( int height, int width, int colAlignment, DComplex* buffer, int ldim, Grid g );
#endif /* WITHOUT_COMPLEX */
/* Print */
void Print_MC_Star_Single( char* msg, MC_Star_Single A );
void Print_MC_Star_Double( char* msg, MC_Star_Double A );
#ifndef WITHOUT_COMPLEX
void Print_MC_Star_SComplex( char* msg, MC_Star_SComplex A );
void Print_MC_Star_DComplex( char* msg, MC_Star_DComplex A );
#endif /* WITHOUT_COMPLEX */

/*******************************************************************************
 * [MD,* ] manipulation routines                                               *
 ******************************************************************************/
/* Create empty */
MD_Star_Single CreateEmpty_MD_Star_Single( Grid g );
MD_Star_Double CreateEmpty_MD_Star_Double( Grid g );
#ifndef WITHOUT_COMPLEX
MD_Star_SComplex CreateEmpty_MD_Star_SComplex( Grid g );
MD_Star_DComplex CreateEmpty_MD_Star_DComplex( Grid g );
#endif /* WITHOUT_COMPLEX */
/* Create from exiting buffer */
MD_Star_Single Register_MD_Star_Single
( int height, int width, int colAlignment, float* buffer, int ldim, Grid g );
MD_Star_Double Register_MD_Star_Double
( int height, int width, int colAlignment, double* buffer, int ldim, Grid g );
#ifndef WITHOUT_COMPLEX
MD_Star_SComplex Register_MD_Star_SComplex
( int height, int width, int colAlignment, SComplex* buffer, int ldim, Grid g );
MD_Star_DComplex Register_MD_Star_DComplex
( int height, int width, int colAlignment, DComplex* buffer, int ldim, Grid g );
#endif /* WITHOUT_COMPLEX */
/* Print */
void Print_MD_Star_Single( char* msg, MD_Star_Single A );
void Print_MD_Star_Double( char* msg, MD_Star_Double A );
#ifndef WITHOUT_COMPLEX
void Print_MD_Star_SComplex( char* msg, MD_Star_SComplex A );
void Print_MD_Star_DComplex( char* msg, MD_Star_DComplex A );
#endif /* WITHOUT_COMPLEX */

/*******************************************************************************
 * [MR,MC] manipulation routines                                               *
 ******************************************************************************/
/* Create empty */
MR_MC_Single CreateEmpty_MR_MC_Single( Grid g );
MR_MC_Double CreateEmpty_MR_MC_Double( Grid g );
#ifndef WITHOUT_COMPLEX
MR_MC_SComplex CreateEmpty_MR_MC_SComplex( Grid g );
MR_MC_DComplex CreateEmpty_MR_MC_DComplex( Grid g );
#endif /* WITHOUT_COMPLEX */
/* Create from exiting buffer */
MR_MC_Single Register_MR_MC_Single
( int height, int width, int colAlignment, int rowAlignment, 
  float* buffer, int ldim, Grid g );
MR_MC_Double Register_MR_MC_Double
( int height, int width, int colAlignment, int rowAlignment,
  double* buffer, int ldim, Grid g );
#ifndef WITHOUT_COMPLEX
MR_MC_SComplex Register_MR_MC_SComplex
( int height, int width, int colAlignment, int rowAlignment, 
  SComplex* buffer, int ldim, Grid g );
MR_MC_DComplex Register_MR_MC_DComplex
( int height, int width, int colAlignment, int rowAlignment,
  DComplex* buffer, int ldim, Grid g );
#endif /* WITHOUT_COMPLEX */
/* Print */
void Print_MR_MC_Single( char* msg, MR_MC_Single A );
void Print_MR_MC_Double( char* msg, MR_MC_Double A );
#ifndef WITHOUT_COMPLEX
void Print_MR_MC_SComplex( char* msg, MR_MC_SComplex A );
void Print_MR_MC_DComplex( char* msg, MR_MC_DComplex A );
#endif /* WITHOUT_COMPLEX */

/*******************************************************************************
 * [MR,* ] manipulation routines                                               *
 ******************************************************************************/
/* Create empty */
MR_Star_Single CreateEmpty_MR_Star_Single( Grid g );
MR_Star_Double CreateEmpty_MR_Star_Double( Grid g );
#ifndef WITHOUT_COMPLEX
MR_Star_SComplex CreateEmpty_MR_Star_SComplex( Grid g );
MR_Star_DComplex CreateEmpty_MR_Star_DComplex( Grid g );
#endif /* WITHOUT_COMPLEX */
/* Create from exiting buffer */
MR_Star_Single Register_MR_Star_Single
( int height, int width, int colAlignment, float* buffer, int ldim, Grid g );
MR_Star_Double Register_MR_Star_Double
( int height, int width, int colAlignment, double* buffer, int ldim, Grid g );
#ifndef WITHOUT_COMPLEX
MR_Star_SComplex Register_MR_Star_SComplex
( int height, int width, int colAlignment, SComplex* buffer, int ldim, Grid g );
MR_Star_DComplex Register_MR_Star_DComplex
( int height, int width, int colAlignment, DComplex* buffer, int ldim, Grid g );
#endif /* WITHOUT_COMPLEX */
/* Print */
void Print_MR_Star_Single( char* msg, MR_Star_Single A );
void Print_MR_Star_Double( char* msg, MR_Star_Double A );
#ifndef WITHOUT_COMPLEX
void Print_MR_Star_SComplex( char* msg, MR_Star_SComplex A );
void Print_MR_Star_DComplex( char* msg, MR_Star_DComplex A );
#endif /* WITHOUT_COMPLEX */

/*******************************************************************************
 * [* ,MC] manipulation routines                                               *
 ******************************************************************************/
/* Create empty */
Star_MC_Single CreateEmpty_Star_MC_Single( Grid g );
Star_MC_Double CreateEmpty_Star_MC_Double( Grid g );
#ifndef WITHOUT_COMPLEX
Star_MC_SComplex CreateEmpty_Star_MC_SComplex( Grid g );
Star_MC_DComplex CreateEmpty_Star_MC_DComplex( Grid g );
#endif /* WITHOUT_COMPLEX */
/* Create from exiting buffer */
Star_MC_Single Register_Star_MC_Single
( int height, int width, int rowAlignment, float* buffer, int ldim, Grid g );
Star_MC_Double Register_Star_MC_Double
( int height, int width, int rowAlignment, double* buffer, int ldim, Grid g );
#ifndef WITHOUT_COMPLEX
Star_MC_SComplex Register_Star_MC_SComplex
( int height, int width, int rowAlignment, SComplex* buffer, int ldim, Grid g );
Star_MC_DComplex Register_Star_MC_DComplex
( int height, int width, int rowAlignment, DComplex* buffer, int ldim, Grid g );
#endif /* WITHOUT_COMPLEX */
/* Print */
void Print_Star_MC_Single( char* msg, Star_MC_Single A );
void Print_Star_MC_Double( char* msg, Star_MC_Double A );
#ifndef WITHOUT_COMPLEX
void Print_Star_MC_SComplex( char* msg, Star_MC_SComplex A );
void Print_Star_MC_DComplex( char* msg, Star_MC_DComplex A );
#endif /* WITHOUT_COMPLEX */

/*******************************************************************************
 * [* ,MD] manipulation routines                                               *
 ******************************************************************************/
/* Create empty */
Star_MD_Single CreateEmpty_Star_MD_Single( Grid g );
Star_MD_Double CreateEmpty_Star_MD_Double( Grid g );
#ifndef WITHOUT_COMPLEX
Star_MD_SComplex CreateEmpty_Star_MD_SComplex( Grid g );
Star_MD_DComplex CreateEmpty_Star_MD_DComplex( Grid g );
#endif /* WITHOUT_COMPLEX */
/* Create from exiting buffer */
Star_MD_Single Register_Star_MD_Single
( int height, int width, int rowAlignment, float* buffer, int ldim, Grid g );
Star_MD_Double Register_Star_MD_Double
( int height, int width, int rowAlignment, double* buffer, int ldim, Grid g );
#ifndef WITHOUT_COMPLEX
Star_MD_SComplex Register_Star_MD_SComplex
( int height, int width, int rowAlignment, SComplex* buffer, int ldim, Grid g );
Star_MD_DComplex Register_Star_MD_DComplex
( int height, int width, int rowAlignment, DComplex* buffer, int ldim, Grid g );
#endif /* WITHOUT_COMPLEX */
/* Print */
void Print_Star_MD_Single( char* msg, Star_MD_Single A );
void Print_Star_MD_Double( char* msg, Star_MD_Double A );
#ifndef WITHOUT_COMPLEX
void Print_Star_MD_SComplex( char* msg, Star_MD_SComplex A );
void Print_Star_MD_DComplex( char* msg, Star_MD_DComplex A );
#endif /* WITHOUT_COMPLEX */

/*******************************************************************************
 * [* ,MR] manipulation routines                                               *
 ******************************************************************************/
/* Create empty */
Star_MR_Single CreateEmpty_Star_MR_Single( Grid g );
Star_MR_Double CreateEmpty_Star_MR_Double( Grid g );
#ifndef WITHOUT_COMPLEX
Star_MR_SComplex CreateEmpty_Star_MR_SComplex( Grid g );
Star_MR_DComplex CreateEmpty_Star_MR_DComplex( Grid g );
#endif /* WITHOUT_COMPLEX */
/* Create from exiting buffer */
Star_MR_Single Register_Star_MR_Single
( int height, int width, int rowAlignment, float* buffer, int ldim, Grid g );
Star_MR_Double Register_Star_MR_Double
( int height, int width, int rowAlignment, double* buffer, int ldim, Grid g );
#ifndef WITHOUT_COMPLEX
Star_MR_SComplex Register_Star_MR_SComplex
( int height, int width, int rowAlignment, SComplex* buffer, int ldim, Grid g );
Star_MR_DComplex Register_Star_MR_DComplex
( int height, int width, int rowAlignment, DComplex* buffer, int ldim, Grid g );
#endif /* WITHOUT_COMPLEX */
/* Print */
void Print_Star_MR_Single( char* msg, Star_MR_Single A );
void Print_Star_MR_Double( char* msg, Star_MR_Double A );
#ifndef WITHOUT_COMPLEX
void Print_Star_MR_SComplex( char* msg, Star_MR_SComplex A );
void Print_Star_MR_DComplex( char* msg, Star_MR_DComplex A );
#endif /* WITHOUT_COMPLEX */

/*******************************************************************************
 * [* ,* ] manipulation routines                                               *
 ******************************************************************************/
/* Create empty */
Star_Star_Single CreateEmpty_Star_Star_Single( Grid g );
Star_Star_Double CreateEmpty_Star_Star_Double( Grid g );
#ifndef WITHOUT_COMPLEX
Star_Star_SComplex CreateEmpty_Star_Star_SComplex( Grid g );
Star_Star_DComplex CreateEmpty_Star_Star_DComplex( Grid g );
#endif /* WITHOUT_COMPLEX */
/* Create from exiting buffer */
Star_Star_Single Register_Star_Star_Single
( int height, int width, float* buffer, int ldim, Grid g );
Star_Star_Double Register_Star_Star_Double
( int height, int width, double* buffer, int ldim, Grid g );
#ifndef WITHOUT_COMPLEX
Star_Star_SComplex Register_Star_Star_SComplex
( int height, int width, SComplex* buffer, int ldim, Grid g );
Star_Star_DComplex Register_Star_Star_DComplex
( int height, int width, DComplex* buffer, int ldim, Grid g );
#endif /* WITHOUT_COMPLEX */
/* Print */
void Print_Star_Star_Single( char* msg, Star_Star_Single A );
void Print_Star_Star_Double( char* msg, Star_Star_Double A );
#ifndef WITHOUT_COMPLEX
void Print_Star_Star_SComplex( char* msg, Star_Star_SComplex A );
void Print_Star_Star_DComplex( char* msg, Star_Star_DComplex A );
#endif /* WITHOUT_COMPLEX */

/*******************************************************************************
 * [* ,VC] manipulation routines                                               *
 ******************************************************************************/
/* Create empty */
Star_VC_Single CreateEmpty_Star_VC_Single( Grid g );
Star_VC_Double CreateEmpty_Star_VC_Double( Grid g );
#ifndef WITHOUT_COMPLEX
Star_VC_SComplex CreateEmpty_Star_VC_SComplex( Grid g );
Star_VC_DComplex CreateEmpty_Star_VC_DComplex( Grid g );
#endif /* WITHOUT_COMPLEX */
/* Create from exiting buffer */
Star_VC_Single Register_Star_VC_Single
( int height, int width, int rowAlignment, float* buffer, int ldim, Grid g );
Star_VC_Double Register_Star_VC_Double
( int height, int width, int rowAlignment, double* buffer, int ldim, Grid g );
#ifndef WITHOUT_COMPLEX
Star_VC_SComplex Register_Star_VC_SComplex
( int height, int width, int rowAlignment, SComplex* buffer, int ldim, Grid g );
Star_VC_DComplex Register_Star_VC_DComplex
( int height, int width, int rowAlignment, DComplex* buffer, int ldim, Grid g );
#endif /* WITHOUT_COMPLEX */
/* Print */
void Print_Star_VC_Single( char* msg, Star_VC_Single A );
void Print_Star_VC_Double( char* msg, Star_VC_Double A );
#ifndef WITHOUT_COMPLEX
void Print_Star_VC_SComplex( char* msg, Star_VC_SComplex A );
void Print_Star_VC_DComplex( char* msg, Star_VC_DComplex A );
#endif /* WITHOUT_COMPLEX */

/*******************************************************************************
 * [* ,VR] manipulation routines                                               *
 ******************************************************************************/
/* Create empty */
Star_VR_Single CreateEmpty_Star_VR_Single( Grid g );
Star_VR_Double CreateEmpty_Star_VR_Double( Grid g );
#ifndef WITHOUT_COMPLEX
Star_VR_SComplex CreateEmpty_Star_VR_SComplex( Grid g );
Star_VR_DComplex CreateEmpty_Star_VR_DComplex( Grid g );
#endif /* WITHOUT_COMPLEX */
/* Create from exiting buffer */
Star_VR_Single Register_Star_VR_Single
( int height, int width, int rowAlignment, float* buffer, int ldim, Grid g );
Star_VR_Double Register_Star_VR_Double
( int height, int width, int rowAlignment, double* buffer, int ldim, Grid g );
#ifndef WITHOUT_COMPLEX
Star_VR_SComplex Register_Star_VR_SComplex
( int height, int width, int rowAlignment, SComplex* buffer, int ldim, Grid g );
Star_VR_DComplex Register_Star_VR_DComplex
( int height, int width, int rowAlignment, DComplex* buffer, int ldim, Grid g );
#endif /* WITHOUT_COMPLEX */
/* Print */
void Print_Star_VR_Single( char* msg, Star_VR_Single A );
void Print_Star_VR_Double( char* msg, Star_VR_Double A );
#ifndef WITHOUT_COMPLEX
void Print_Star_VR_SComplex( char* msg, Star_VR_SComplex A );
void Print_Star_VR_DComplex( char* msg, Star_VR_DComplex A );
#endif /* WITHOUT_COMPLEX */

/*******************************************************************************
 * [VC,* ] manipulation routines                                               *
 ******************************************************************************/
/* Create empty */
VC_Star_Single CreateEmpty_VC_Star_Single( Grid g );
VC_Star_Double CreateEmpty_VC_Star_Double( Grid g );
#ifndef WITHOUT_COMPLEX
VC_Star_SComplex CreateEmpty_VC_Star_SComplex( Grid g );
VC_Star_DComplex CreateEmpty_VC_Star_DComplex( Grid g );
#endif /* WITHOUT_COMPLEX */
/* Create from exiting buffer */
VC_Star_Single Register_VC_Star_Single
( int height, int width, int colAlignment, float* buffer, int ldim, Grid g );
VC_Star_Double Register_VC_Star_Double
( int height, int width, int colAlignment, double* buffer, int ldim, Grid g );
#ifndef WITHOUT_COMPLEX
VC_Star_SComplex Register_VC_Star_SComplex
( int height, int width, int colAlignment, SComplex* buffer, int ldim, Grid g );
VC_Star_DComplex Register_VC_Star_DComplex
( int height, int width, int colAlignment, DComplex* buffer, int ldim, Grid g );
#endif /* WITHOUT_COMPLEX */
/* Print */
void Print_VC_Star_Single( char* msg, VC_Star_Single A );
void Print_VC_Star_Double( char* msg, VC_Star_Double A );
#ifndef WITHOUT_COMPLEX
void Print_VC_Star_SComplex( char* msg, VC_Star_SComplex A );
void Print_VC_Star_DComplex( char* msg, VC_Star_DComplex A );
#endif /* WITHOUT_COMPLEX */

/*******************************************************************************
 * [VR,* ] manipulation routines                                               *
 ******************************************************************************/
/* Create empty */
VR_Star_Single CreateEmpty_VR_Star_Single( Grid g );
VR_Star_Double CreateEmpty_VR_Star_Double( Grid g );
#ifndef WITHOUT_COMPLEX
VR_Star_SComplex CreateEmpty_VR_Star_SComplex( Grid g );
VR_Star_DComplex CreateEmpty_VR_Star_DComplex( Grid g );
#endif /* WITHOUT_COMPLEX */
/* Create from exiting buffer */
VR_Star_Single Register_VR_Star_Single
( int height, int width, int colAlignment, float* buffer, int ldim, Grid g );
VR_Star_Double Register_VR_Star_Double
( int height, int width, int colAlignment, double* buffer, int ldim, Grid g );
#ifndef WITHOUT_COMPLEX
VR_Star_SComplex Register_VR_Star_SComplex
( int height, int width, int colAlignment, SComplex* buffer, int ldim, Grid g );
VR_Star_DComplex Register_VR_Star_DComplex
( int height, int width, int colAlignment, DComplex* buffer, int ldim, Grid g );
#endif /* WITHOUT_COMPLEX */
/* Print */
void Print_VR_Star_Single( char* msg, VR_Star_Single A );
void Print_VR_Star_Double( char* msg, VR_Star_Double A );
#ifndef WITHOUT_COMPLEX
void Print_VR_Star_SComplex( char* msg, VR_Star_SComplex A );
void Print_VR_Star_DComplex( char* msg, VR_Star_DComplex A );
#endif /* WITHOUT_COMPLEX */

/* Utilities */
int LocalLength
( int globalLength, int myIndex, int alignment, int modulus );

/* LAPACK-level interface */
void CholSingle( char uplo, MC_MR_Single A );
void CholDouble( char uplo, MC_MR_Double A );
#ifndef WITHOUT_COMPLEX
void CholSComplex( char uplo, MC_MR_SComplex A );
void CholDComplex( char uplo, MC_MR_DComplex A );
#endif

#ifndef WITHOUT_PMRRR
void HermitianEigDouble
( char uplo, int tryForHighAccuracy, 
  char job, int a, int b, double u, double v,
  MC_MR_Double A, Star_VR_Double w, MC_MR_Double Z );

void GeneralizedHermitianEigDouble
( int genEigInt, char uplo, int tryForHighAccuracy, 
  char job, int a, int b, double u, double v,
  MC_MR_Double A, MC_MR_Double B, 
  Star_VR_Double w, MC_MR_Double Z );

void HermitianEigDouble_OnlyEigvals
( char uplo, int tryForHighAccuracy, 
  char job, int a, int b, double u, double v,
  MC_MR_Double A, Star_VR_Double w );

void GeneralizedHermitianEigDouble_OnlyEigvals
( int genEigInt, char uplo, int tryForHighAccuracy,
  char job, int a, int b, double u, double v,
  MC_MR_Double A, MC_MR_Double B, Star_VR_Double w );

#ifndef WITHOUT_COMPLEX
void HermitianEigDComplex
( char uplo, int tryForHighAccuracy,
  char job, int a, int b, double u, double v,
  MC_MR_DComplex A, Star_VR_Double w, MC_MR_DComplex Z );

void GeneralizedHermitianEigDComplex
( int genEigInt, char uplo, int tryForHighAccuracy,
  char job, int a, int b, double u, double v,
  MC_MR_DComplex A, MC_MR_DComplex B, Star_VR_Double w, MC_MR_DComplex Z );

void HermitianEigDComplex_OnlyEigvals
( char uplo, int tryForHighAccuracy,
  char job, int a, int b, double u, double v,
  MC_MR_DComplex A, Star_VR_Double w );

void GeneralizedHermitianEigDComplex_OnlyEigvals
( int genEigInt, char uplo, int tryForHighAccuracy,
  char job, int a, int b, double u, double v,
  MC_MR_DComplex A, MC_MR_DComplex B, Star_VR_Double w );
#endif /* WITHOUT_COMPLEX */
#endif /* WITHOUT_PMRRR */

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* ELEMENTAL_EXPORTS_C_H */

