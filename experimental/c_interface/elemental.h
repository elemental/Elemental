/*
   Copyright (c) 2011-2012, Jack Poulson
   All rights reserved.

   This file is a part of a prototype interface to a few generalized eigensolver
   routines of Elemental.

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
#include "mpi.h"

typedef unsigned GridHandle;
typedef unsigned RealDistMatHandle;
typedef unsigned ComplexDistMatHandle;
typedef unsigned RealDistColVecHandle;

/* Environment controls */
void Initialize( int* argc, char** argv[] );
void Finalize();
void SetBlocksize( int blocksize );
int Blocksize();
/* TODO: other tuning parameters */

/* Process grid management */
GridHandle CreateGrid( MPI_Comm comm );
void FreeGrid( GridHandle handle );
int GridHeight( GridHandle gridHandle );
int GridWidth( GridHandle gridHandle );
int GridSize( GridHandle gridHandle );
int GridRow( GridHandle gridHandle );
int GridCol( GridHandle gridHandle );
int GridRank( GridHandle gridHandle );

/* Distributed matrix management */
RealDistMatHandle CreateEmptyRealDistMat( GridHandle gridHandle );
ComplexDistMatHandle CreateEmptyComplexDistMat( GridHandle gridHandle );
RealDistMatHandle RegisterRealDistMat
( int height, int width, int colAlignment, int rowAlignment, 
  double* buffer, int ldim, GridHandle gridHandle );
ComplexDistMatHandle RegisterComplexDistMat
( int height, int width, int colAlignment, int rowAlignment, 
  void* voidBuffer, int ldim, GridHandle gridHandle );
void FreeRealDistMat( RealDistMatHandle handle );
void FreeComplexDistMat( ComplexDistMatHandle handle );
void PrintRealDistMat( RealDistMatHandle AHandle );
void PrintComplexDistMat( ComplexDistMatHandle AHandle );

/* Distributed column vector management */
RealDistColVecHandle CreateEmptyRealDistColVec( GridHandle gridHandle );
void FreeRealDistColVec( RealDistColVecHandle handle );
void PrintRealDistColVec( RealDistColVecHandle AHandle );

/* Generalized Hermitian-definite eigensolvers */
void SymmetricAxBx
( RealDistMatHandle AHandle, RealDistMatHandle BHandle,
  RealDistColVecHandle* wHandle, RealDistMatHandle* XHandle );
void SymmetricAxBxPartialRange
( RealDistMatHandle AHandle, RealDistMatHandle BHandle,
  RealDistColVecHandle* wHandle, RealDistMatHandle* XHandle,
  double a, double b );
void SymmetricAxBxPartialIndices
( RealDistMatHandle AHandle, RealDistMatHandle BHandle,
  RealDistColVecHandle* wHandle, RealDistMatHandle* XHandle,
  int a, int b );
void HermitianAxBx
( ComplexDistMatHandle AHandle, ComplexDistMatHandle BHandle,
  RealDistColVecHandle* wHandle, ComplexDistMatHandle* XHandle );
void HermitianAxBxPartialRange
( ComplexDistMatHandle AHandle, ComplexDistMatHandle BHandle,
  RealDistColVecHandle* wHandle, ComplexDistMatHandle* XHandle,
  double a, double b );
void HermitianAxBxPartialIndices
( ComplexDistMatHandle AHandle, ComplexDistMatHandle BHandle,
  RealDistColVecHandle* wHandle, ComplexDistMatHandle* XHandle,
  int a, int b );

/* Utility functions */
int LocalLength( int n, int shift, int modulus );
