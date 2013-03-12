/*
   Copyright (c) 2011-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "mpi.h"

/* TODO: Test these macros in Windows */
#ifdef ELEMPY_DLL
# ifdef ELEMPY_EXPORTS
#  define ELEMPY_API __declspec(dllexport)
# else
#  define ELEMPY_API __declspec(dllimport)
# endif
#else
# define ELEMPY_API extern
#endif

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
ELEMPY_API void Initialize( int* argc, char** argv[] );
ELEMPY_API void Finalize();
ELEMPY_API void SetBlocksize( int blocksize );
ELEMPY_API int Blocksize();
/* TODO: other tuning parameters */

/* Process grid management */
ELEMPY_API const Grid* DefaultGrid();
ELEMPY_API Grid* CreateGrid( MPI_Comm comm );
ELEMPY_API void FreeGrid( Grid* grid );
ELEMPY_API int GridHeight( const Grid* grid );
ELEMPY_API int GridWidth( const Grid* grid );
ELEMPY_API int GridSize( const Grid* grid );
ELEMPY_API int GridRow( const Grid* grid );
ELEMPY_API int GridCol( const Grid* grid );
ELEMPY_API int GridRank( const Grid* grid );

/* Matrix management */
ELEMPY_API Mat* CreateMat();
ELEMPY_API CpxMat* CreateCpxMat();
ELEMPY_API void ResizeMat( Mat* A, int height, int width );
ELEMPY_API void ResizeCpxMat( CpxMat* A, int height, int width );
ELEMPY_API void AttachToMat
( Mat* A, int height, int width, void* buffer, int ldim );
ELEMPY_API void AttachToCpxMat
( CpxMat* A, int height, int width, void* buffer, int ldim );
ELEMPY_API void UniformMat( Mat* A, int height, int width );
ELEMPY_API void UniformCpxMat( CpxMat* A, int height, int width );
ELEMPY_API void FreeMat( Mat* A );
ELEMPY_API void FreeCpxMat( CpxMat* A );
ELEMPY_API void PrintMat( const Mat* A );
ELEMPY_API void PrintCpxMat( const CpxMat* A );
ELEMPY_API int MatHeight( const Mat* A );
ELEMPY_API int MatWidth( const Mat* A );
ELEMPY_API int MatLDim( const Mat* A );
ELEMPY_API int CpxMatHeight( const CpxMat* A );
ELEMPY_API int CpxMatWidth( const CpxMat* A );
ELEMPY_API int CpxMatLDim( const CpxMat* A );
ELEMPY_API double GetMatEntry( const Mat* A, int i, int j );
ELEMPY_API void SetMatEntry( Mat* A, int i, int j, double alpha );
/* How to handle passing complex numbers? */
/* 
ELEMPY_API Cpx GetCpxMatEntry( const CpxMat* A, int i, int j );
ELEMPY_API void SetCpxMatEntry( CpxMat* A, int i, int j, Cpx alpha );
*/

/* [MC,MR] management */
ELEMPY_API DistMat* CreateDistMat( const Grid* grid );
ELEMPY_API CpxDistMat* CreateCpxDistMat( const Grid* grid );
ELEMPY_API void ResizeDistMat( DistMat* A, int height, int width );
ELEMPY_API void ResizeCpxDistMat( CpxDistMat* A, int height, int width );
ELEMPY_API void AttachToDistMat
( DistMat* A,
  int height, int width, int colAlignment, int rowAlignment, 
  void* buffer, int ldim, const Grid* grid );
ELEMPY_API void AttachToCpxDistMat
( CpxDistMat* A,
  int height, int width, int colAlignment, int rowAlignment, 
  void* buffer, int ldim, const Grid* grid );
ELEMPY_API void UniformDistMat( DistMat* A, int height, int width );
ELEMPY_API void UniformCpxDistMat( CpxDistMat* A, int height, int width );
ELEMPY_API void FreeDistMat( DistMat* A );
ELEMPY_API void FreeCpxDistMat( CpxDistMat* A );
ELEMPY_API void PrintDistMat( const DistMat* A );
ELEMPY_API void PrintCpxDistMat( const CpxDistMat* A );
ELEMPY_API int DistMatHeight( const DistMat* A );
ELEMPY_API int CpxDistMatHeight( const CpxDistMat* A );
ELEMPY_API int DistMatWidth( const DistMat* A );
ELEMPY_API int CpxDistMatWidth( const CpxDistMat* A );
ELEMPY_API int DistMatLocalHeight( const DistMat* A );
ELEMPY_API int CpxDistMatLocalHeight( const CpxDistMat* A );
ELEMPY_API int DistMatLocalWidth( const DistMat* A );
ELEMPY_API int CpxDistMatLocalWidth( const CpxDistMat* A );
ELEMPY_API int DistMatLDim( const DistMat* A );
ELEMPY_API int CpxDistMatLDim( const CpxDistMat* A );
ELEMPY_API int DistMatColShift( const DistMat* A );
ELEMPY_API int CpxDistMatColShift( const CpxDistMat* A );
ELEMPY_API int DistMatRowShift( const DistMat* A );
ELEMPY_API int CpxDistMatRowShift( const CpxDistMat* A );
ELEMPY_API int DistMatColStride( const DistMat* A );
ELEMPY_API int CpxDistMatColStride( const CpxDistMat* A );
ELEMPY_API int DistMatRowStride( const DistMat* A );
ELEMPY_API int CpxDistMatRowStride( const CpxDistMat* A );
ELEMPY_API double GetDistMatEntry( const DistMat* A, int i, int j );
ELEMPY_API void SetDistMatEntry( DistMat* A, int i, int j, double alpha );
ELEMPY_API double GetLocalDistMatEntry
( const DistMat* A, int iLocal, int jLocal );
ELEMPY_API void SetLocalDistMatEntry
( DistMat* A, int iLocal, int jLocal, double alpha );
/* How to handle passing complex numbers? */

/* [VC,*] management */
ELEMPY_API DistMat_VC_STAR* CreateDistMat_VC_STAR( const Grid* grid );
ELEMPY_API CpxDistMat_VC_STAR* CreateCpxDistMat_VC_STAR( const Grid* grid );
ELEMPY_API void FreeDistMat_VC_STAR( DistMat_VC_STAR* A );
ELEMPY_API void FreeCpxDistMat_VC_STAR( CpxDistMat_VC_STAR* A );
ELEMPY_API void PrintDistMat_VC_STAR( const DistMat_VC_STAR* A );
ELEMPY_API void PrintCpxDistMat_VC_STAR( const CpxDistMat_VC_STAR* A );
ELEMPY_API int DistMatHeight_VC_STAR( const DistMat_VC_STAR* A );
ELEMPY_API int CpxDistMatHeight_VC_STAR( const CpxDistMat_VC_STAR* A );
ELEMPY_API int DistMatWidth_VC_STAR( const DistMat_VC_STAR* A );
ELEMPY_API int CpxDistMatWidth_VC_STAR( const CpxDistMat_VC_STAR* A );
ELEMPY_API int DistMatLocalHeight_VC_STAR( const DistMat_VC_STAR* A );
ELEMPY_API int CpxDistMatLocalHeight_VC_STAR( const CpxDistMat_VC_STAR* A );
ELEMPY_API int DistMatLocalWidth_VC_STAR( const DistMat_VC_STAR* A );
ELEMPY_API int CpxDistMatLocalWidth_VC_STAR( const CpxDistMat_VC_STAR* A );
ELEMPY_API int DistMatLDim_VC_STAR( const DistMat_VC_STAR* A );
ELEMPY_API int CpxDistMatLDim_VC_STAR( const CpxDistMat_VC_STAR* A );
ELEMPY_API int DistMatColShift_VC_STAR( const DistMat_VC_STAR* A );
ELEMPY_API int CpxDistMatColShift_VC_STAR( const CpxDistMat_VC_STAR* A );
ELEMPY_API int DistMatRowShift_VC_STAR( const DistMat_VC_STAR* A );
ELEMPY_API int CpxDistMatRowShift_VC_STAR( const CpxDistMat_VC_STAR* A );
ELEMPY_API int DistMatColStride_VC_STAR( const DistMat_VC_STAR* A );
ELEMPY_API int CpxDistMatColStride_VC_STAR( const CpxDistMat_VC_STAR* A );
ELEMPY_API int DistMatRowStride_VC_STAR( const DistMat_VC_STAR* A );
ELEMPY_API int CpxDistMatRowStride_VC_STAR( const CpxDistMat_VC_STAR* A );
ELEMPY_API double GetDistMatEntry_VC_STAR
( const DistMat_VC_STAR* A, int i, int j );
ELEMPY_API void SetDistMatEntry_VC_STAR
( DistMat_VC_STAR* A, int i, int j, double alpha );
ELEMPY_API double GetLocalDistMatEntry_VC_STAR
( const DistMat_VC_STAR* A, int iLocal, int jLocal );
ELEMPY_API void SetLocalDistMatEntry_VC_STAR
( DistMat_VC_STAR* A, int iLocal, int jLocal, double alpha );
/* How to handle passing complex numbers? */

/* [VR,*] management */
ELEMPY_API DistMat_VR_STAR* CreateDistMat_VR_STAR( const Grid* grid );
ELEMPY_API CpxDistMat_VR_STAR* CreateCpxDistMat_VR_STAR( const Grid* grid );
ELEMPY_API void FreeDistMat_VR_STAR( DistMat_VR_STAR* A );
ELEMPY_API void FreeCpxDistMat_VR_STAR( CpxDistMat_VR_STAR* A );
ELEMPY_API void PrintDistMat_VR_STAR( const DistMat_VR_STAR* A );
ELEMPY_API void PrintCpxDistMat_VR_STAR( const CpxDistMat_VR_STAR* A );
ELEMPY_API int DistMatHeight_VR_STAR( const DistMat_VR_STAR* A );
ELEMPY_API int CpxDistMatHeight_VR_STAR( const CpxDistMat_VR_STAR* A );
ELEMPY_API int DistMatWidth_VR_STAR( const DistMat_VR_STAR* A );
ELEMPY_API int CpxDistMatWidth_VR_STAR( const CpxDistMat_VR_STAR* A );
ELEMPY_API int DistMatLocalHeight_VR_STAR( const DistMat_VR_STAR* A );
ELEMPY_API int CpxDistMatLocalHeight_VR_STAR( const CpxDistMat_VR_STAR* A );
ELEMPY_API int DistMatLocalWidth_VR_STAR( const DistMat_VR_STAR* A );
ELEMPY_API int CpxDistMatLocalWidth_VR_STAR( const CpxDistMat_VR_STAR* A );
ELEMPY_API int DistMatLDim_VR_STAR( const DistMat_VR_STAR* A );
ELEMPY_API int CpxDistMatLDim_VR_STAR( const CpxDistMat_VR_STAR* A );
ELEMPY_API int DistMatColShift_VR_STAR( const DistMat_VR_STAR* A );
ELEMPY_API int CpxDistMatColShift_VR_STAR( const CpxDistMat_VR_STAR* A );
ELEMPY_API int DistMatRowShift_VR_STAR( const DistMat_VR_STAR* A );
ELEMPY_API int CpxDistMatRowShift_VR_STAR( const CpxDistMat_VR_STAR* A );
ELEMPY_API int DistMatColStride_VR_STAR( const DistMat_VR_STAR* A );
ELEMPY_API int CpxDistMatColStride_VR_STAR( const CpxDistMat_VR_STAR* A );
ELEMPY_API int DistMatRowStride_VR_STAR( const DistMat_VR_STAR* A );
ELEMPY_API int CpxDistMatRowStride_VR_STAR( const CpxDistMat_VR_STAR* A );
ELEMPY_API double GetDistMatEntry_VR_STAR
( const DistMat_VR_STAR* A, int i, int j );
ELEMPY_API void SetDistMatEntry_VR_STAR
( DistMat_VR_STAR* A, int i, int j, double alpha );
ELEMPY_API double GetLocalDistMatEntry_VR_STAR
( const DistMat_VR_STAR* A, int iLocal, int jLocal );
ELEMPY_API void SetLocalDistMatEntry_VR_STAR
( DistMat_VR_STAR* A, int iLocal, int jLocal, double alpha );

/* QR factorization */
ELEMPY_API void ExplicitQR( DistMat* A, DistMat* R );
ELEMPY_API void CpxExplicitQR( CpxDistMat* A, CpxDistMat* R );

/* Singular Value Decomposition */
ELEMPY_API void SVD( DistMat* A, DistMat_VR_STAR* s, DistMat* V );
ELEMPY_API void CpxSVD( CpxDistMat* A, DistMat_VR_STAR* s, CpxDistMat* V );
ELEMPY_API void SingularValues( DistMat* A, DistMat_VR_STAR* s );
ELEMPY_API void CpxSingularValues( CpxDistMat* A, DistMat_VR_STAR* s );

/* Generalized Hermitian-definite eigensolvers */
ELEMPY_API void SymmetricAxBx
( DistMat* A, DistMat* B, DistMat_VR_STAR* w, DistMat* X );
ELEMPY_API void SymmetricAxBxRange
( DistMat* A, DistMat* B, DistMat_VR_STAR* w, DistMat* X, 
  double a, double b );
ELEMPY_API void SymmetricAxBxIndices
( DistMat* A, DistMat* B, DistMat_VR_STAR* w, DistMat* X,
  int a, int b );
ELEMPY_API void HermitianAxBx
( CpxDistMat* A, CpxDistMat* B, DistMat_VR_STAR* w, CpxDistMat* X );
ELEMPY_API void HermitianAxBxRange
( CpxDistMat* A, CpxDistMat* B, DistMat_VR_STAR* w, CpxDistMat* X,
  double a, double b );
ELEMPY_API void HermitianAxBxIndices
( CpxDistMat* A, CpxDistMat* B, DistMat_VR_STAR* w, CpxDistMat* X,
  int a, int b );

/* Utility functions */
ELEMPY_API int Length( int n, int shift, int modulus );
/* TODO: Shift */

} // extern "C"
