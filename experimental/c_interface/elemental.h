/*
   Copyright (c) 2011-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
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
