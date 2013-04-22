/*
   Copyright (c) 2011-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <limits>
#include "elemental.hpp"
// TODO: Switch to elemental-lite.hpp
using namespace elem;

#include "./elempy.hpp"

extern "C" {

//
// Environment controls
//

ELEMPY_API void Initialize( int* argc, char** argv[] )
{ elem::Initialize( *argc, *argv ); }

ELEMPY_API void Finalize()
{ elem::Finalize(); }

ELEMPY_API void SetBlocksize( int blocksize )
{ elem::SetBlocksize( blocksize ); }

ELEMPY_API int Blocksize()
{ return elem::Blocksize(); }

//
// Process grid management
//

ELEMPY_API const Grid* DefaultGrid()
{
    const Grid& grid = elem::DefaultGrid();
    return &grid;
}

ELEMPY_API Grid* CreateGrid( MPI_Comm comm )
{ return new Grid( comm ); }

ELEMPY_API int GridHeight( const Grid* grid )
{ return grid->Height(); }

ELEMPY_API int GridWidth( const Grid* grid )
{ return grid->Width(); }

ELEMPY_API int GridSize( const Grid* grid )
{ return grid->Size(); }

ELEMPY_API int GridRow( const Grid* grid )
{ return grid->MCRank(); }

ELEMPY_API int GridCol( const Grid* grid )
{ return grid->MRRank(); }

ELEMPY_API int GridRank( const Grid* grid )
{ return grid->Rank(); }

ELEMPY_API void FreeGrid( Grid* grid )
{ delete grid; }

//
// Matrix management
//

ELEMPY_API Mat* CreateMat()
{ return new Mat; }

ELEMPY_API CpxMat* CreateCpxMat()
{ return new CpxMat; }

ELEMPY_API void ResizeMat( Mat* A, int height, int width )
{ A->ResizeTo( height, width ); }

ELEMPY_API void ResizeCpxMat( CpxMat* A, int height, int width )
{ A->ResizeTo( height, width ); }

ELEMPY_API void 
AttachToMat( Mat* A, int height, int width, void* buffer, int ldim )
{ A->Attach( height, width, reinterpret_cast<double*>(buffer), ldim ); }

ELEMPY_API void 
AttachToCpxMat( CpxMat* A, int height, int width, void* buffer, int ldim )
{ A->Attach( height, width, reinterpret_cast<Cpx*>(buffer), ldim ); }

ELEMPY_API void FreeMat( Mat* A )
{ delete A; }

ELEMPY_API void FreeCpxMat( CpxMat* A )
{ delete A; }

ELEMPY_API void PrintMat( const Mat* A )
{ A->Print(); }

ELEMPY_API void PrintCpxMat( const CpxMat* A )
{ A->Print(); }

ELEMPY_API int MatHeight( const Mat* A )
{ return A->Height(); }

ELEMPY_API int CpxMatHeight( const CpxMat* A )
{ return A->Height(); }

ELEMPY_API int MatWidth( const Mat* A )
{ return A->Width(); }

ELEMPY_API int CpxMatWidth( const CpxMat* A )
{ return A->Width(); }

ELEMPY_API int MatLDim( const Mat* A )
{ return A->LDim(); }

ELEMPY_API int CpxMatLDim( const CpxMat* A )
{ return A->LDim(); }

ELEMPY_API double GetMatEntry( const Mat* A, int i, int j )
{ return A->Get( i, j ); }

ELEMPY_API void SetMatEntry( Mat* A, int i, int j, double alpha )
{ A->Set( i, j, alpha ); }

ELEMPY_API void GetCpxMatEntry
( const CpxMat* A, int i, int j, double* real, double* imag )
{ 
    const Complex<double>& alpha = A->Get( i, j );
    *real = RealPart(alpha);
    *imag = ImagPart(alpha);
}
    
ELEMPY_API void SetCpxMatEntry
( CpxMat* A, int i, int j, double real, double imag )
{ A->Set( i, j, Complex<double>(real,imag) ); }

ELEMPY_API void* MatBuffer( Mat* A )
{ return A->Buffer(); }

ELEMPY_API void* CpxMatBuffer( CpxMat* A )
{ return A->Buffer(); }

//
// Distributed matrix management
//

ELEMPY_API DistMat* CreateDistMat( const Grid* grid )
{ return new DistMat( *grid ); }

ELEMPY_API CpxDistMat* CreateCpxDistMat( const Grid* grid )
{ return new CpxDistMat( *grid ); }

ELEMPY_API void ResizeDistMat( DistMat* A, int height, int width )
{ A->ResizeTo( height, width ); }

ELEMPY_API void ResizeCpxDistMat( CpxDistMat* A, int height, int width )
{ A->ResizeTo( height, width ); }

ELEMPY_API void AttachToDistMat
( DistMat* A, 
  int height, int width, int colAlignment, int rowAlignment, 
  void* buffer, int ldim, const Grid* grid )
{ A->Attach
  ( height, width, colAlignment, rowAlignment, 
    reinterpret_cast<double*>(buffer), ldim, *grid ); }

ELEMPY_API void AttachToCpxDistMat
( CpxDistMat* A, 
  int height, int width, int colAlignment, int rowAlignment, 
  void* buffer, int ldim, const Grid* grid )
{ A->Attach
  ( height, width, colAlignment, rowAlignment, 
    reinterpret_cast<Cpx*>(buffer), ldim, *grid ); }

ELEMPY_API void FreeDistMat( DistMat* A )
{ delete A; }

ELEMPY_API void FreeCpxDistMat( CpxDistMat* A )
{ delete A; }

ELEMPY_API void PrintDistMat( const DistMat* A )
{ A->Print(); }

ELEMPY_API void PrintCpxDistMat( const CpxDistMat* A )
{ A->Print(); }

ELEMPY_API int DistMatHeight( const DistMat* A )
{ return A->Height(); }

ELEMPY_API int CpxDistMatHeight( const CpxDistMat* A )
{ return A->Height(); }

ELEMPY_API int DistMatWidth( const DistMat* A )
{ return A->Width(); }

ELEMPY_API int CpxDistMatWidth( const CpxDistMat* A )
{ return A->Width(); }

ELEMPY_API int DistMatLocalHeight( const DistMat* A )
{ return A->LocalHeight(); }

ELEMPY_API int CpxDistMatLocalHeight( const CpxDistMat* A )
{ return A->LocalHeight(); }

ELEMPY_API int DistMatLocalWidth( const DistMat* A )
{ return A->LocalWidth(); }

ELEMPY_API int CpxDistMatLocalWidth( const CpxDistMat* A )
{ return A->LocalWidth(); }

ELEMPY_API int DistMatLDim( const DistMat* A )
{ return A->LDim(); }

ELEMPY_API int CpxDistMatLDim( const CpxDistMat* A )
{ return A->LDim(); }

ELEMPY_API int DistMatColShift( const DistMat* A )
{ return A->ColShift(); }

ELEMPY_API int CpxDistMatColShift( const CpxDistMat* A )
{ return A->ColShift(); }

ELEMPY_API int DistMatRowShift( const DistMat* A )
{ return A->RowShift(); }

ELEMPY_API int CpxDistMatRowShift( const CpxDistMat* A )
{ return A->RowShift(); }

ELEMPY_API int DistMatColStride( const DistMat* A )
{ return A->ColStride(); }

ELEMPY_API int CpxDistMatColStride( const CpxDistMat* A )
{ return A->ColStride(); }

ELEMPY_API int DistMatRowStride( const DistMat* A )
{ return A->RowStride(); }

ELEMPY_API int CpxDistMatRowStride( const CpxDistMat* A )
{ return A->RowStride(); }

ELEMPY_API double GetDistMatEntry( const DistMat* A, int i, int j )
{ return A->Get( i, j ); }

ELEMPY_API void SetDistMatEntry( DistMat* A, int i, int j, double alpha )
{ A->Set( i, j, alpha ); }

ELEMPY_API double 
GetLocalDistMatEntry( const DistMat* A, int iLocal, int jLocal )
{ return A->GetLocal( iLocal, jLocal ); }

ELEMPY_API void 
SetLocalDistMatEntry( DistMat* A, int iLocal, int jLocal, double alpha )
{ A->SetLocal( iLocal, jLocal, alpha ); }

ELEMPY_API void GetCpxDistMatEntry
( const CpxDistMat* A, int i, int j, double* real, double* imag )
{ 
    const Complex<double>& alpha = A->Get( i, j );
    *real = RealPart(alpha);
    *imag = ImagPart(alpha);
}
    
ELEMPY_API void SetCpxDistMatEntry
( CpxDistMat* A, int i, int j, double real, double imag )
{ A->Set( i, j, Complex<double>(real,imag) ); }

ELEMPY_API void GetLocalCpxDistMatEntry
( const CpxDistMat* A, int iLocal, int jLocal, double* real, double* imag )
{ 
    const Complex<double>& alpha = A->GetLocal( iLocal, jLocal );
    *real = RealPart(alpha);
    *imag = ImagPart(alpha);
}
    
ELEMPY_API void SetLocalCpxDistMatEntry
( CpxDistMat* A, int iLocal, int jLocal, double real, double imag )
{ A->SetLocal( iLocal, jLocal, Complex<double>(real,imag) ); }

ELEMPY_API void* DistMatBuffer( DistMat* A )
{ return A->Buffer(); }

ELEMPY_API void* CpxDistMatBuffer( CpxDistMat* A )
{ return A->Buffer(); }

//
// [VC,* ] management
//

ELEMPY_API DistMat_VC_STAR* CreateDistMat_VC_STAR( const Grid* grid )
{ return new DistMat_VC_STAR( *grid ); }

ELEMPY_API CpxDistMat_VC_STAR* CreateCpxDistMat_VC_STAR( const Grid* grid )
{ return new CpxDistMat_VC_STAR( *grid ); }

ELEMPY_API void 
ResizeDistMat_VC_STAR( DistMat_VC_STAR* A, int height, int width )
{ A->ResizeTo( height, width ); }

ELEMPY_API void 
ResizeCpxDistMat_VC_STAR( CpxDistMat_VC_STAR* A, int height, int width )
{ A->ResizeTo( height, width ); }

ELEMPY_API void FreeDistMat_VC_STAR( DistMat_VC_STAR* A )
{ delete A; }

ELEMPY_API void FreeCpxDistMat_VC_STAR( CpxDistMat_VC_STAR* A )
{ delete A; }

ELEMPY_API void PrintDistMat_VC_STAR( const DistMat_VC_STAR* A )
{ A->Print(); }

ELEMPY_API void PrintCpxDistMat_VC_STAR( const CpxDistMat_VC_STAR* A )
{ A->Print(); }

ELEMPY_API int DistMatHeight_VC_STAR( const DistMat_VC_STAR* A )
{ return A->Height(); }

ELEMPY_API int CpxDistMatHeight_VC_STAR( const CpxDistMat_VC_STAR* A )
{ return A->Height(); }

ELEMPY_API int DistMatWidth_VC_STAR( const DistMat_VC_STAR* A )
{ return A->Width(); }

ELEMPY_API int CpxDistMatWidth_VC_STAR( const CpxDistMat_VC_STAR* A )
{ return A->Width(); }

ELEMPY_API int DistMatLocalHeight_VC_STAR( const DistMat_VC_STAR* A )
{ return A->LocalHeight(); }

ELEMPY_API int CpxDistMatLocalHeight_VC_STAR( const CpxDistMat_VC_STAR* A )
{ return A->LocalHeight(); }

ELEMPY_API int DistMatLocalWidth_VC_STAR( const DistMat_VC_STAR* A )
{ return A->LocalWidth(); }

ELEMPY_API int CpxDistMatLocalWidth_VC_STAR( const CpxDistMat_VC_STAR* A )
{ return A->LocalWidth(); }

ELEMPY_API int DistMatLDim_VC_STAR( const DistMat_VC_STAR* A )
{ return A->LDim(); }

ELEMPY_API int CpxDistMatLDim_VC_STAR( const CpxDistMat_VC_STAR* A )
{ return A->LDim(); }

ELEMPY_API int DistMatColShift_VC_STAR( const DistMat_VC_STAR* A )
{ return A->ColShift(); }

ELEMPY_API int CpxDistMatColShift_VC_STAR( const CpxDistMat_VC_STAR* A )
{ return A->ColShift(); }

ELEMPY_API int DistMatRowShift_VC_STAR( const DistMat_VC_STAR* A )
{ return A->RowShift(); }

ELEMPY_API int CpxDistMatRowShift_VC_STAR( const CpxDistMat_VC_STAR* A )
{ return A->RowShift(); }

ELEMPY_API int DistMatColStride_VC_STAR( const DistMat_VC_STAR* A )
{ return A->ColStride(); }

ELEMPY_API int CpxDistMatColStride_VC_STAR( const CpxDistMat_VC_STAR* A )
{ return A->ColStride(); }

ELEMPY_API int DistMatRowStride_VC_STAR( const DistMat_VC_STAR* A )
{ return A->RowStride(); }

ELEMPY_API int CpxDistMatRowStride_VC_STAR( const CpxDistMat_VC_STAR* A )
{ return A->RowStride(); }

ELEMPY_API double 
GetDistMatEntry_VC_STAR( const DistMat_VC_STAR* A, int i, int j )
{ return A->Get( i, j ); }

ELEMPY_API void 
SetDistMatEntry_VC_STAR( DistMat_VC_STAR* A, int i, int j, double alpha )
{ A->Set( i, j, alpha ); }

ELEMPY_API double GetLocalDistMatEntry_VC_STAR
( const DistMat_VC_STAR* A, int iLocal, int jLocal )
{ return A->GetLocal( iLocal, jLocal ); }

ELEMPY_API void SetLocalDistMatEntry_VC_STAR
( DistMat_VC_STAR* A, int iLocal, int jLocal, double alpha )
{ A->SetLocal( iLocal, jLocal, alpha ); }

ELEMPY_API void GetCpxDistMatEntry_VC_STAR
( const CpxDistMat_VC_STAR* A, int i, int j, double* real, double* imag )
{ 
    const Complex<double>& alpha = A->Get( i, j );
    *real = RealPart(alpha);
    *imag = ImagPart(alpha);
}
    
ELEMPY_API void SetCpxDistMatEntry_VC_STAR
( CpxDistMat_VC_STAR* A, int i, int j, double real, double imag )
{ A->Set( i, j, Complex<double>(real,imag) ); }

ELEMPY_API void GetLocalCpxDistMatEntry_VC_STAR
( const CpxDistMat_VC_STAR* A, 
  int iLocal, int jLocal, double* real, double* imag )
{ 
    const Complex<double>& alpha = A->GetLocal( iLocal, jLocal );
    *real = RealPart(alpha);
    *imag = ImagPart(alpha);
}
    
ELEMPY_API void SetLocalCpxDistMatEntry_VC_STAR
( CpxDistMat_VC_STAR* A, int iLocal, int jLocal, double real, double imag )
{ A->SetLocal( iLocal, jLocal, Complex<double>(real,imag) ); }

ELEMPY_API void* DistMatBuffer_VC_STAR( DistMat_VC_STAR* A )
{ return A->Buffer(); }

ELEMPY_API void* CpxDistMatBuffer_VC_STAR( CpxDistMat_VC_STAR* A )
{ return A->Buffer(); }

//
// [VR,* ] management
//

ELEMPY_API DistMat_VR_STAR* CreateDistMat_VR_STAR( const Grid* grid )
{ return new DistMat_VR_STAR( *grid ); }

ELEMPY_API CpxDistMat_VR_STAR* CreateCpxDistMat_VR_STAR( const Grid* grid )
{ return new CpxDistMat_VR_STAR( *grid ); }

ELEMPY_API void 
ResizeDistMat_VR_STAR( DistMat_VR_STAR* A, int height, int width )
{ A->ResizeTo( height, width ); }

ELEMPY_API void 
ResizeCpxDistMat_VR_STAR( CpxDistMat_VR_STAR* A, int height, int width )
{ A->ResizeTo( height, width ); }

ELEMPY_API void FreeDistMat_VR_STAR( DistMat_VR_STAR* A )
{ delete A; }

ELEMPY_API void FreeCpxDistMat_VR_STAR( CpxDistMat_VR_STAR* A )
{ delete A; }

ELEMPY_API void PrintDistMat_VR_STAR( const DistMat_VR_STAR* A )
{ A->Print(); }

ELEMPY_API void PrintCpxDistMat_VR_STAR( const CpxDistMat_VR_STAR* A )
{ A->Print(); }

ELEMPY_API int DistMatHeight_VR_STAR( const DistMat_VR_STAR* A )
{ return A->Height(); }

ELEMPY_API int CpxDistMatHeight_VR_STAR( const CpxDistMat_VR_STAR* A )
{ return A->Height(); }

ELEMPY_API int DistMatWidth_VR_STAR( const DistMat_VR_STAR* A )
{ return A->Width(); }

ELEMPY_API int CpxDistMatWidth_VR_STAR( const CpxDistMat_VR_STAR* A )
{ return A->Width(); }

ELEMPY_API int DistMatLocalHeight_VR_STAR( const DistMat_VR_STAR* A )
{ return A->LocalHeight(); }

ELEMPY_API int CpxDistMatLocalHeight_VR_STAR( const CpxDistMat_VR_STAR* A )
{ return A->LocalHeight(); }

ELEMPY_API int DistMatLocalWidth_VR_STAR( const DistMat_VR_STAR* A )
{ return A->LocalWidth(); }

ELEMPY_API int CpxDistMatLocalWidth_VR_STAR( const CpxDistMat_VR_STAR* A )
{ return A->LocalWidth(); }

ELEMPY_API int DistMatLDim_VR_STAR( const DistMat_VR_STAR* A )
{ return A->LDim(); }

ELEMPY_API int CpxDistMatLDim_VR_STAR( const CpxDistMat_VR_STAR* A )
{ return A->LDim(); }

ELEMPY_API int DistMatColShift_VR_STAR( const DistMat_VR_STAR* A )
{ return A->ColShift(); }

ELEMPY_API int CpxDistMatColShift_VR_STAR( const CpxDistMat_VR_STAR* A )
{ return A->ColShift(); }

ELEMPY_API int DistMatRowShift_VR_STAR( const DistMat_VR_STAR* A )
{ return A->RowShift(); }

ELEMPY_API int CpxDistMatRowShift_VR_STAR( const CpxDistMat_VR_STAR* A )
{ return A->RowShift(); }

ELEMPY_API int DistMatColStride_VR_STAR( const DistMat_VR_STAR* A )
{ return A->ColStride(); }

ELEMPY_API int CpxDistMatColStride_VR_STAR( const CpxDistMat_VR_STAR* A )
{ return A->ColStride(); }

ELEMPY_API int DistMatRowStride_VR_STAR( const DistMat_VR_STAR* A )
{ return A->RowStride(); }

ELEMPY_API int CpxDistMatRowStride_VR_STAR( const CpxDistMat_VR_STAR* A )
{ return A->RowStride(); }

ELEMPY_API double 
GetDistMatEntry_VR_STAR( const DistMat_VR_STAR* A, int i, int j )
{ return A->Get( i, j ); }

ELEMPY_API void 
SetDistMatEntry_VR_STAR( DistMat_VR_STAR* A, int i, int j, double alpha )
{ A->Set( i, j, alpha ); }

ELEMPY_API double GetLocalDistMatEntry_VR_STAR
( const DistMat_VR_STAR* A, int iLocal, int jLocal )
{ return A->GetLocal( iLocal, jLocal ); }

ELEMPY_API void SetLocalDistMatEntry_VR_STAR
( DistMat_VR_STAR* A, int iLocal, int jLocal, double alpha )
{ A->SetLocal( iLocal, jLocal, alpha ); }

ELEMPY_API void GetCpxDistMatEntry_VR_STAR
( const CpxDistMat_VR_STAR* A, int i, int j, double* real, double* imag )
{ 
    const Complex<double>& alpha = A->Get( i, j );
    *real = RealPart(alpha);
    *imag = ImagPart(alpha);
}
    
ELEMPY_API void SetCpxDistMatEntry_VR_STAR
( CpxDistMat_VR_STAR* A, int i, int j, double real, double imag )
{ A->Set( i, j, Complex<double>(real,imag) ); }

ELEMPY_API void GetLocalCpxDistMatEntry_VR_STAR
( const CpxDistMat_VR_STAR* A, 
  int iLocal, int jLocal, double* real, double* imag )
{ 
    const Complex<double>& alpha = A->GetLocal( iLocal, jLocal );
    *real = RealPart(alpha);
    *imag = ImagPart(alpha);
}
    
ELEMPY_API void SetLocalCpxDistMatEntry_VR_STAR
( CpxDistMat_VR_STAR* A, int iLocal, int jLocal, double real, double imag )
{ A->SetLocal( iLocal, jLocal, Complex<double>(real,imag) ); }

ELEMPY_API void* DistMatBuffer_VR_STAR( DistMat_VR_STAR* A )
{ return A->Buffer(); }

ELEMPY_API void* CpxDistMatBuffer_VR_STAR( CpxDistMat_VR_STAR* A )
{ return A->Buffer(); }

//
// QR factorization
//

ELEMPY_API void ExplicitQR( DistMat* A, DistMat* R )
{ elem::qr::Explicit( *A, *R ); }

ELEMPY_API void CpxExplicitQR( CpxDistMat* A, CpxDistMat* R )
{ elem::qr::Explicit( *A, *R ); }

//
// Singular Value Decomposition
//

ELEMPY_API void SVD( DistMat* A, DistMat_VR_STAR* s, DistMat* V )
{ elem::SVD( *A, *s, *V ); }

ELEMPY_API void CpxSVD( CpxDistMat* A, DistMat_VR_STAR* s, CpxDistMat* V )
{ elem::SVD( *A, *s, *V ); }

ELEMPY_API void SingularValues( DistMat* A, DistMat_VR_STAR* s );
ELEMPY_API void CpxSingularValues( CpxDistMat* A, DistMat_VR_STAR* s );

//
// Generalized Hermitian-definite eigensolvers for A X = B X \Lambda
//

ELEMPY_API void SymmetricAxBx
( DistMat* A, DistMat* B, DistMat_VR_STAR* w, DistMat* X )
{ HermitianGenDefiniteEig( AXBX, LOWER, *A, *B, *w, *X ); }

ELEMPY_API void SymmetricAxBxRange
( DistMat* A, DistMat* B, DistMat_VR_STAR* w, DistMat* X,
  double a, double b )
{ HermitianGenDefiniteEig( AXBX, LOWER, *A, *B, *w, *X, a, b ); }

ELEMPY_API void SymmetricAxBxIndices
( DistMat* A, DistMat* B, DistMat_VR_STAR* w, DistMat* X,
  int a, int b )
{ HermitianGenDefiniteEig( AXBX, LOWER, *A, *B, *w, *X, a, b ); }

ELEMPY_API void HermitianAxBx
( CpxDistMat* A, CpxDistMat* B, DistMat_VR_STAR* w, CpxDistMat* X )
{ HermitianGenDefiniteEig( AXBX, LOWER, *A, *B, *w, *X ); }

ELEMPY_API void HermitianAxBxRange
( CpxDistMat* A, CpxDistMat* B, DistMat_VR_STAR* w, CpxDistMat* X,
  double a, double b )
{ HermitianGenDefiniteEig( AXBX, LOWER, *A, *B, *w, *X, a, b ); }

ELEMPY_API void HermitianAxBxIndices
( CpxDistMat* A, CpxDistMat* B, DistMat_VR_STAR* w, CpxDistMat* X,
  int a, int b )
{ HermitianGenDefiniteEig( AXBX, LOWER, *A, *B, *w, *X, a, b ); }

//
// Special matrices
//

ELEMPY_API void UniformMat( Mat* A, int height, int width )
{ Uniform( *A, height, width ); }

ELEMPY_API void UniformCpxMat( CpxMat* A, int height, int width )
{ Uniform( *A, height, width ); }

ELEMPY_API void UniformDistMat( DistMat* A, int height, int width )
{ Uniform( *A, height, width ); }

ELEMPY_API void UniformCpxDistMat( CpxDistMat* A, int height, int width )
{ Uniform( *A, height, width ); }

//
// Utilities
//

ELEMPY_API int Length( int n, int shift, int modulus )
{ return elem::Length<int>(n,shift,modulus); }

} // extern "C"
