/*
   Copyright (c) 2011-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "mpi.h"

/* Environment controls */
void Initialize( int* argc, char** argv[] );
void Finalize();
void SetBlocksize( int blocksize );
int Blocksize();
/* TODO: other tuning parameters */

/* Process grid management */
const Grid* DefaultGrid();
Grid* CreateGrid( MPI_Comm comm );
void FreeGrid( Grid** grid );
int GridHeight( const Grid* grid );
int GridWidth( const Grid* grid );
int GridSize( const Grid* grid );
int GridRow( const Grid* grid );
int GridCol( const Grid* grid );
int GridRank( const Grid* grid );

/* [MC,MR] management */
DistMatrix<double>* CreateDistMat( const Grid* grid );
DistMatrix<Complex<double> >* CreateCpxDistMat( const Grid* grid );
void AttachToDistMat
( DistMatrix<double>* A,
  int height, int width, int colAlignment, int rowAlignment, 
  void* buffer, int ldim, const Grid* grid );
void AttachToCpxDistMat
( DistMatrix<Complex<double> >* A,
  int height, int width, int colAlignment, int rowAlignment, 
  void* buffer, int ldim, const Grid* grid );
/* TODO: Resize */
void UniformDistMat( DistMatrix<double>* A, int height, int width );
void UniformCpxDistMat
( DistMatrix<Complex<double> >* A, int height, int width );
void FreeDistMat( DistMatrix<double>** A );
void FreeCpxDistMat( DistMatrix<Complex<double> >** A );
void PrintDistMat( const DistMatrix<double>* A );
void PrintCpxDistMat( const DistMatrix<Complex<double> >* A );

/* [VC,*] management */
DistMatrix<double,VC,STAR>* CreateDistMat_VC_STAR( const Grid* grid );
DistMatrix<Complex<double>,VC,STAR>* 
CreateCpxDistMat_VC_STAR( const Grid* grid );
void FreeDistMat_VC_STAR( DistMatrix<double,VC,STAR>** A );
void FreeCpxDistMat_VC_STAR( DistMatrix<Complex<double>,VC,STAR>** A );
void PrintDistMat_VC_STAR( const DistMatrix<double,VC,STAR>* A );
void PrintCpxDistMat_VC_STAR( const DistMatrix<Complex<double>,VC,STAR>* A );

/* [VR,*] management */
DistMatrix<double,VR,STAR>* CreateDistMat_VR_STAR( const Grid* grid );
DistMatrix<Complex<double>,VR,STAR>* 
CreateCpxDistMat_VR_STAR( const Grid* grid );
void FreeDistMat_VR_STAR( DistMatrix<double,VR,STAR>** A );
void FreeCpxDistMat_VR_STAR( DistMatrix<Complex<double>,VR,STAR>** A );
void PrintDistMat_VR_STAR( const DistMatrix<double,VR,STAR>* A );
void PrintCpxDistMat_VR_STAR( const DistMatrix<Complex<double>,VR,STAR>* A );

/* Singular Value Decomposition */
void SVD
( DistMatrix<double>* A, 
  DistMatrix<double,VR,STAR>* s, 
  DistMatrix<double>* V );
void CpxSVD
( DistMatrix<Complex<double> >* A, 
  DistMatrix<double,VR,STAR>* s, 
  DistMatrix<Complex<double> >* V );

/* Generalized Hermitian-definite eigensolvers */
void SymmetricAxBx
( DistMatrix<double>* A, DistMatrix<double>* B,
  DistMatrix<double,VR,STAR>* w, DistMatrix<double>* X );
void ElemSymmetricAxBxRange
( DistMatrix<double>* A, DistMatrix<double>* B,
  DistMatrix<double,VR,STAR>* w, DistMatrix<double>* X,
  double a, double b );
void ElemSymmetricAxBxIndices
( DistMatrix<double>* A, DistMatrix<double>* B,
  DistMatrix<double,VR,STAR>* w, DistMatrix<double>* X,
  int a, int b );
void HermitianAxBx
( DistMatrix<Complex<double> >* A, DistMatrix<Complex<double> >* B,
  DistMatrix<double,VR,STAR>* w, DistMatrix<Complex<double> >* X );
void HermitianAxBxRange
( DistMatrix<Complex<double> >* A, DistMatrix<Complex<double> >* B,
  DistMatrix<double,VR,STAR>* w, DistMatrix<Complex<double> >* X,
  double a, double b );
void HermitianAxBxIndices
( DistMatrix<Complex<double> >* A, DistMatrix<Complex<double> >* B,
  DistMatrix<double,VR,STAR>* w, DistMatrix<Complex<double> >* X,
  int a, int b );

/* Utility functions */
int ElemLength( int n, int shift, int modulus );
/* TODO: ElemShift */
