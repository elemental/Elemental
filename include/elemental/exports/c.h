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

// Make sure that all of our configuration definitions are pulled in
#include "elemental/config.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifndef WITHOUT_COMPLEX
// We should not assume C99 support and use <complex.h>
typedef struct { float real; float imag; } ElementalSComplex;
typedef struct { double real; double imag; } ElementalDComplex;
#endif // WITHOUT_COMPLEX

void ElementalInit( int* argc, char** argv[] );
void ElementalFinalize();
int ElementalBlocksize();
void ElementalSetBlocksize( int blocksize );
void ElementalPushBlocksizeStack( int blocksize );
void ElementalPopBlocksizeStack();

void ElementalClearGrids();
int ElementalDefaultGrid( MPI_Comm comm );
int ElementalGrid( MPI_Comm comm, int r, int c );

int ElementalGridHeight( int gridHandle );
int ElementalGridWidth( int gridHandle );
int ElementalGridSize( int gridHandle );
int ElementalInGrid( int gridHandle );
int ElementalGridVCRank( int gridHandle );
int ElementalGridVRRank( int gridHandle );
int ElementalGridMCRank( int gridHandle );
int ElementalGridMRRank( int gridHandle );
MPI_Comm ElementalGridVCComm( int gridHandle );
MPI_Comm ElementalGridVRComm( int gridHandle );
MPI_Comm ElementalGridMCComm( int gridHandle );
MPI_Comm ElementalGridMRComm( int gridHandle );

void ElementalClearDistMatrices();

int ElementalDistMatrixDouble // Default to [MC,MR]
( int height, int width, int colAlignment, int rowAlignment,
  double* buffer, int ldim, int gridHandle );
int ElementalDistMatrix_MC_MR_Double
( int height, int width, int colAlignment, int rowAlignment,
  double* buffer, int ldim, int gridHandle );
void ElementalDistMatrixDoublePrint( char* msg, int distMatrixDoubleHandle );

#ifndef WITHOUT_COMPLEX
int ElementalDistMatrixDComplex // Default to [MC,MR]
( int height, int width, int colAlignment, int rowAlignment,
  ElementalDComplex* buffer, int ldim, int gridHandle );  
int ElementalDistMatrix_MC_MR_DComplex
( int height, int width, int colAlignment, int rowAlignment,
  ElementalDComplex* buffer, int ldim, int gridHandle );  
void ElementalDistMatrixDComplexPrint( char* msg, int distMatrixDComplexHandle );
#endif // WITHOUT_COMPLEX

int ElementalLocalLength
( int globalLength, int myIndex, int alignment, int modulus );

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* ELEMENTAL_EXPORTS_C_H */

