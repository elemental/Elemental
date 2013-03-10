/*
   Copyright (c) 2011-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "mpi.h"

typedef unsigned ElemGrid;

typedef unsigned ElemMat;
typedef unsigned ElemCpxMat;

typedef unsigned ElemDistMat;
typedef unsigned ElemCpxDistMat;

typedef unsigned ElemDistMat_VC_STAR;
typedef unsigned ElemDistMat_VR_STAR;

/* Environment controls */
void ElemInitialize( int* argc, char** argv[] );
void ElemFinalize();
void ElemSetBlocksize( int blocksize );
int ElemBlocksize();
/* TODO: other tuning parameters */

/* Process grid management */
ElemGrid ElemDefaultGrid();
ElemGrid ElemCreateGrid( MPI_Comm comm );
void ElemFreeGrid( ElemGrid grid );
int ElemGridHeight( ElemGrid grid );
int ElemGridWidth( ElemGrid grid );
int ElemGridSize( ElemGrid grid );
int ElemGridRow( ElemGrid grid );
int ElemGridCol( ElemGrid grid );
int ElemGridRank( ElemGrid grid );

/* [MC,MR] management */
ElemDistMat ElemCreateDistMat( ElemGrid grid );
ElemCpxDistMat ElemCreateCpxDistMat( ElemGrid grid );
ElemDistMat ElemRegisterDistMat
( int height, int width, int colAlignment, int rowAlignment, 
  double* buffer, int ldim, ElemGrid grid );
ElemCpxDistMat ElemRegisterCpxDistMat
( int height, int width, int colAlignment, int rowAlignment, 
  void* voidBuffer, int ldim, ElemGrid grid );
/* TODO: Resize */
void ElemUniformDistMat( ElemDistMat A, int height, int width );
void ElemUniformCpxDistMat( ElemCpxDistMat A, int height, int width );
void ElemFreeDistMat( ElemDistMat A );
void ElemFreeCpxDistMat( ElemCpxDistMat A );
void ElemPrintDistMat( ElemDistMat A );
void ElemPrintCpxDistMat( ElemCpxDistMat A );

/* [VC,*] management */
ElemDistMat_VC_STAR ElemCreateDistMat_VC_STAR( ElemGrid grid );
void ElemFreeDistMat_VC_STAR( ElemDistMat_VC_STAR A );
void ElemPrintDistMat_VC_STAR( ElemDistMat_VC_STAR A );

/* [VR,*] management */
ElemDistMat_VR_STAR ElemCreateDistMat_VR_STAR( ElemGrid grid );
void ElemFreeDistMat_VR_STAR( ElemDistMat_VR_STAR A );
void ElemPrintDistMat_VR_STAR( ElemDistMat_VR_STAR A );

/* Generalized Hermitian-definite eigensolvers */
void ElemSymmetricAxBx
( ElemDistMat A, ElemDistMat B,
  ElemDistMat_VR_STAR* w, ElemDistMat* X );
void ElemSymmetricAxBxRange
( ElemDistMat A, ElemDistMat B,
  ElemDistMat_VR_STAR* w, ElemDistMat* X, 
  double a, double b );
void ElemSymmetricAxBxIndices
( ElemDistMat A, ElemDistMat B,
  ElemDistMat_VR_STAR* w, ElemDistMat* X, 
  int a, int b );
void HermitianAxBx
( ElemCpxDistMat A, ElemCpxDistMat B,
  ElemDistMat_VR_STAR* w, ElemCpxDistMat* X );
void HermitianAxBxRange
( ElemCpxDistMat A, ElemCpxDistMat B,
  ElemDistMat_VR_STAR* w, ElemCpxDistMat* X,
  double a, double b );
void HermitianAxBxIndices
( ElemCpxDistMat A, ElemCpxDistMat B,
  ElemDistMat_VR_STAR* w, ElemCpxDistMat* X,
  int a, int b );

/* Utility functions */
int ElemLength( int n, int shift, int modulus );
/* TODO: ElemShift */
