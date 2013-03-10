/*
   Copyright (c) 2011-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "mpi.h"

typedef Complex<double> Cpx;

typedef Matrix<double> Mat;
typedef DistMatrix<double> DistMat;
typedef DistMatrix<double,VC,STAR> DistMat_VC_STAR;
typedef DistMatrix<double,VR,STAR> DistMat_VR_STAR;

typedef Matrix<Cpx> CpxMat;
typedef DistMatrix<Cpx> CpxDistMat;
typedef DistMatrix<Cpx,VC,STAR> CpxDistMat_VC_STAR;
typedef DistMatrix<Cpx,VR,STAR> CpxDistMat_VR_STAR;

extern "C" {

/* Environment controls */
void Initialize( int* argc, char** argv[] );
void Finalize();
void SetBlocksize( int blocksize );
int Blocksize();
/* TODO: other tuning parameters */

/* Process grid management */
const Grid* DefaultGrid();
Grid* CreateGrid( MPI_Comm comm );
void FreeGrid( Grid* grid );
int GridHeight( const Grid* grid );
int GridWidth( const Grid* grid );
int GridSize( const Grid* grid );
int GridRow( const Grid* grid );
int GridCol( const Grid* grid );
int GridRank( const Grid* grid );

/* Matrix management */
Mat* CreateMat();
CpxMat* CreateCpxMat();
void ResizeMat( Mat* A, int height, int width );
void ResizeCpxMat( CpxMat* A, int height, int width );
void AttachToMat( Mat* A, int height, int width, void* buffer, int ldim );
void AttachToCpxMat( CpxMat* A, int height, int width, void* buffer, int ldim );
void UniformMat( Mat* A, int height, int width );
void UniformCpxMat( CpxMat* A, int height, int width );
void FreeMat( Mat* A );
void FreeCpxMat( CpxMat* A );
void PrintMat( const Mat* A );
void PrintCpxMat( const CpxMat* A );
int MatHeight( const Mat* A );
int MatWidth( const Mat* A );
int MatLDim( const Mat* A );
int CpxMatHeight( const CpxMat* A );
int CpxMatWidth( const CpxMat* A );
int CpxMatLDim( const CpxMat* A );
double GetMatEntry( const Mat* A, int i, int j );
void SetMatEntry( Mat* A, int i, int j, double alpha );
/* How to handle passing complex numbers? */
/* 
Cpx GetCpxMatEntry( const CpxMat* A, int i, int j );
void SetCpxMatEntry( CpxMat* A, int i, int j, Cpx alpha );
*/

/* [MC,MR] management */
DistMat* CreateDistMat( const Grid* grid );
CpxDistMat* CreateCpxDistMat( const Grid* grid );
void ResizeDistMat( DistMat* A, int height, int width );
void ResizeCpxDistMat( CpxDistMat* A, int height, int width );
void AttachToDistMat
( DistMat* A,
  int height, int width, int colAlignment, int rowAlignment, 
  void* buffer, int ldim, const Grid* grid );
void AttachToCpxDistMat
( CpxDistMat* A,
  int height, int width, int colAlignment, int rowAlignment, 
  void* buffer, int ldim, const Grid* grid );
void UniformDistMat( DistMat* A, int height, int width );
void UniformCpxDistMat( CpxDistMat* A, int height, int width );
void FreeDistMat( DistMat* A );
void FreeCpxDistMat( CpxDistMat* A );
void PrintDistMat( const DistMat* A );
void PrintCpxDistMat( const CpxDistMat* A );
int DistMatHeight( const DistMat* A );
int CpxDistMatHeight( const CpxDistMat* A );
int DistMatWidth( const DistMat* A );
int CpxDistMatWidth( const CpxDistMat* A );
int DistMatLocalHeight( const DistMat* A );
int CpxDistMatLocalHeight( const CpxDistMat* A );
int DistMatLocalWidth( const DistMat* A );
int CpxDistMatLocalWidth( const CpxDistMat* A );
int DistMatLDim( const DistMat* A );
int CpxDistMatLDim( const CpxDistMat* A );
int DistMatColShift( const DistMat* A );
int CpxDistMatColShift( const CpxDistMat* A );
int DistMatRowShift( const DistMat* A );
int CpxDistMatRowShift( const CpxDistMat* A );
int DistMatColStride( const DistMat* A );
int CpxDistMatColStride( const CpxDistMat* A );
int DistMatRowStride( const DistMat* A );
int CpxDistMatRowStride( const CpxDistMat* A );
double GetDistMatEntry( const DistMat* A, int i, int j );
void SetDistMatEntry( DistMat* A, int i, int j, double alpha );
double GetLocalDistMatEntry( const DistMat* A, int iLocal, int jLocal );
void SetLocalDistMatEntry( DistMat* A, int iLocal, int jLocal, double alpha );
/* How to handle passing complex numbers? */

/* [VC,*] management */
DistMat_VC_STAR* CreateDistMat_VC_STAR( const Grid* grid );
CpxDistMat_VC_STAR* CreateCpxDistMat_VC_STAR( const Grid* grid );
void FreeDistMat_VC_STAR( DistMat_VC_STAR* A );
void FreeCpxDistMat_VC_STAR( CpxDistMat_VC_STAR* A );
void PrintDistMat_VC_STAR( const DistMat_VC_STAR* A );
void PrintCpxDistMat_VC_STAR( const CpxDistMat_VC_STAR* A );
int DistMatHeight_VC_STAR( const DistMat_VC_STAR* A );
int CpxDistMatHeight_VC_STAR( const CpxDistMat_VC_STAR* A );
int DistMatWidth_VC_STAR( const DistMat_VC_STAR* A );
int CpxDistMatWidth_VC_STAR( const CpxDistMat_VC_STAR* A );
int DistMatLocalHeight_VC_STAR( const DistMat_VC_STAR* A );
int CpxDistMatLocalHeight_VC_STAR( const CpxDistMat_VC_STAR* A );
int DistMatLocalWidth_VC_STAR( const DistMat_VC_STAR* A );
int CpxDistMatLocalWidth_VC_STAR( const CpxDistMat_VC_STAR* A );
int DistMatLDim_VC_STAR( const DistMat_VC_STAR* A );
int CpxDistMatLDim_VC_STAR( const CpxDistMat_VC_STAR* A );
int DistMatColShift_VC_STAR( const DistMat_VC_STAR* A );
int CpxDistMatColShift_VC_STAR( const CpxDistMat_VC_STAR* A );
int DistMatRowShift_VC_STAR( const DistMat_VC_STAR* A );
int CpxDistMatRowShift_VC_STAR( const CpxDistMat_VC_STAR* A );
int DistMatColStride_VC_STAR( const DistMat_VC_STAR* A );
int CpxDistMatColStride_VC_STAR( const CpxDistMat_VC_STAR* A );
int DistMatRowStride_VC_STAR( const DistMat_VC_STAR* A );
int CpxDistMatRowStride_VC_STAR( const CpxDistMat_VC_STAR* A );
double GetDistMatEntry_VC_STAR( const DistMat_VC_STAR* A, int i, int j );
void SetDistMatEntry_VC_STAR( DistMat_VC_STAR* A, int i, int j, double alpha );
double GetLocalDistMatEntry_VC_STAR
( const DistMat_VC_STAR* A, int iLocal, int jLocal );
void SetLocalDistMatEntry_VC_STAR
( DistMat_VC_STAR* A, int iLocal, int jLocal, double alpha );
/* How to handle passing complex numbers? */

/* [VR,*] management */
DistMat_VR_STAR* CreateDistMat_VR_STAR( const Grid* grid );
CpxDistMat_VR_STAR* CreateCpxDistMat_VR_STAR( const Grid* grid );
void FreeDistMat_VR_STAR( DistMat_VR_STAR* A );
void FreeCpxDistMat_VR_STAR( CpxDistMat_VR_STAR* A );
void PrintDistMat_VR_STAR( const DistMat_VR_STAR* A );
void PrintCpxDistMat_VR_STAR( const CpxDistMat_VR_STAR* A );
int DistMatHeight_VR_STAR( const DistMat_VR_STAR* A );
int CpxDistMatHeight_VR_STAR( const CpxDistMat_VR_STAR* A );
int DistMatWidth_VR_STAR( const DistMat_VR_STAR* A );
int CpxDistMatWidth_VR_STAR( const CpxDistMat_VR_STAR* A );
int DistMatLocalHeight_VR_STAR( const DistMat_VR_STAR* A );
int CpxDistMatLocalHeight_VR_STAR( const CpxDistMat_VR_STAR* A );
int DistMatLocalWidth_VR_STAR( const DistMat_VR_STAR* A );
int CpxDistMatLocalWidth_VR_STAR( const CpxDistMat_VR_STAR* A );
int DistMatLDim_VR_STAR( const DistMat_VR_STAR* A );
int CpxDistMatLDim_VR_STAR( const CpxDistMat_VR_STAR* A );
int DistMatColShift_VR_STAR( const DistMat_VR_STAR* A );
int CpxDistMatColShift_VR_STAR( const CpxDistMat_VR_STAR* A );
int DistMatRowShift_VR_STAR( const DistMat_VR_STAR* A );
int CpxDistMatRowShift_VR_STAR( const CpxDistMat_VR_STAR* A );
int DistMatColStride_VR_STAR( const DistMat_VR_STAR* A );
int CpxDistMatColStride_VR_STAR( const CpxDistMat_VR_STAR* A );
int DistMatRowStride_VR_STAR( const DistMat_VR_STAR* A );
int CpxDistMatRowStride_VR_STAR( const CpxDistMat_VR_STAR* A );
double GetDistMatEntry_VR_STAR( const DistMat_VR_STAR* A, int i, int j );
void SetDistMatEntry_VR_STAR( DistMat_VR_STAR* A, int i, int j, double alpha );
double GetLocalDistMatEntry_VR_STAR
( const DistMat_VR_STAR* A, int iLocal, int jLocal );
void SetLocalDistMatEntry_VR_STAR
( DistMat_VR_STAR* A, int iLocal, int jLocal, double alpha );

/* QR factorization */
void ExplicitQR( DistMat* A, DistMat* R );
void CpxExplicitQR( CpxDistMat* A, CpxDistMat* R );

/* Singular Value Decomposition */
void SVD( DistMat* A, DistMat_VR_STAR* s, DistMat* V );
void CpxSVD( CpxDistMat* A, DistMat_VR_STAR* s, CpxDistMat* V );
void SingularValues( DistMat* A, DistMat_VR_STAR* s );
void CpxSingularValues( CpxDistMat* A, DistMat_VR_STAR* s );

/* Generalized Hermitian-definite eigensolvers */
void SymmetricAxBx( DistMat* A, DistMat* B, DistMat_VR_STAR* w, DistMat* X );
void SymmetricAxBxRange
( DistMat* A, DistMat* B, DistMat_VR_STAR* w, DistMat* X, 
  double a, double b );
void SymmetricAxBxIndices
( DistMat* A, DistMat* B, DistMat_VR_STAR* w, DistMat* X,
  int a, int b );
void HermitianAxBx
( CpxDistMat* A, CpxDistMat* B, DistMat_VR_STAR* w, CpxDistMat* X );
void HermitianAxBxRange
( CpxDistMat* A, CpxDistMat* B, DistMat_VR_STAR* w, CpxDistMat* X,
  double a, double b );
void HermitianAxBxIndices
( CpxDistMat* A, CpxDistMat* B, DistMat_VR_STAR* w, CpxDistMat* X,
  int a, int b );

/* Utility functions */
int Length( int n, int shift, int modulus );
/* TODO: Shift */

} // extern "C"
