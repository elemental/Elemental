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

void Initialize( int* argc, char** argv[] )
{ elem::Initialize( *argc, *argv ); }

void Finalize()
{ elem::Finalize(); }

void SetBlocksize( int blocksize )
{ elem::SetBlocksize( blocksize ); }

int Blocksize()
{ return elem::Blocksize(); }

//
// Process grid management
//

const Grid* DefaultGrid()
{
    const Grid& grid = elem::DefaultGrid();
    return &grid;
}

Grid* CreateGrid( MPI_Comm comm )
{ return new Grid( comm ); }

int GridHeight( const Grid* grid )
{ return grid->Height(); }

int GridWidth( const Grid* grid )
{ return grid->Width(); }

int GridSize( const Grid* grid )
{ return grid->Size(); }

int GridRow( const Grid* grid )
{ return grid->MCRank(); }

int GridCol( const Grid* grid )
{ return grid->MRRank(); }

int GridRank( const Grid* grid )
{ return grid->Rank(); }

void FreeGrid( Grid* grid )
{ delete grid; }

//
// Matrix management
//

Mat* CreateMat()
{ return new Mat; }

CpxMat* CreateCpxMat()
{ return new CpxMat; }

void ResizeMat( Mat* A, int height, int width )
{ A->ResizeTo( height, width ); }

void ResizeCpxMat( CpxMat* A, int height, int width )
{ A->ResizeTo( height, width ); }

void AttachToMat( Mat* A, int height, int width, void* buffer, int ldim )
{ A->Attach( height, width, reinterpret_cast<double*>(buffer), ldim ); }

void AttachToCpxMat( CpxMat* A, int height, int width, void* buffer, int ldim )
{ A->Attach( height, width, reinterpret_cast<Cpx*>(buffer), ldim ); }

void FreeMat( Mat* A )
{ delete A; }

void FreeCpxMat( CpxMat* A )
{ delete A; }

void PrintMat( const Mat* A )
{ A->Print(); }

void PrintCpxMat( const CpxMat* A )
{ A->Print(); }

int MatHeight( const Mat* A )
{ return A->Height(); }

int CpxMatHeight( const CpxMat* A )
{ return A->Height(); }

int MatWidth( const Mat* A )
{ return A->Width(); }

int CpxMatWidth( const CpxMat* A )
{ return A->Width(); }

int MatLDim( const Mat* A )
{ return A->LDim(); }

int CpxMatLDim( const CpxMat* A )
{ return A->LDim(); }

double GetMatEntry( const Mat* A, int i, int j )
{ return A->Get( i, j ); }

void SetMatEntry( Mat* A, int i, int j, double alpha )
{ A->Set( i, j, alpha ); }

//
// Distributed matrix management
//

DistMat* CreateDistMat( const Grid* grid )
{ return new DistMat( *grid ); }

CpxDistMat* CreateCpxDistMat( const Grid* grid )
{ return new CpxDistMat( *grid ); }

void ResizeDistMat( DistMat* A, int height, int width )
{ A->ResizeTo( height, width ); }

void ResizeCpxDistMat( CpxDistMat* A, int height, int width )
{ A->ResizeTo( height, width ); }

void AttachToDistMat
( DistMat* A, 
  int height, int width, int colAlignment, int rowAlignment, 
  void* buffer, int ldim, const Grid* grid )
{ A->Attach
  ( height, width, colAlignment, rowAlignment, 
    reinterpret_cast<double*>(buffer), ldim, *grid ); }

void AttachToCpxDistMat
( CpxDistMat* A, 
  int height, int width, int colAlignment, int rowAlignment, 
  void* buffer, int ldim, const Grid* grid )
{ A->Attach
  ( height, width, colAlignment, rowAlignment, 
    reinterpret_cast<Cpx*>(buffer), ldim, *grid ); }

void FreeDistMat( DistMat* A )
{ delete A; }

void FreeCpxDistMat( CpxDistMat* A )
{ delete A; }

void PrintDistMat( const DistMat* A )
{ A->Print(); }

void PrintCpxDistMat( const CpxDistMat* A )
{ A->Print(); }

int DistMatHeight( const DistMat* A )
{ return A->Height(); }

int CpxDistMatHeight( const CpxDistMat* A )
{ return A->Height(); }

int DistMatWidth( const DistMat* A )
{ return A->Width(); }

int CpxDistMatWidth( const CpxDistMat* A )
{ return A->Width(); }

int DistMatLocalHeight( const DistMat* A )
{ return A->LocalHeight(); }

int CpxDistMatLocalHeight( const CpxDistMat* A )
{ return A->LocalHeight(); }

int DistMatLocalWidth( const DistMat* A )
{ return A->LocalWidth(); }

int CpxDistMatLocalWidth( const CpxDistMat* A )
{ return A->LocalWidth(); }

int DistMatLDim( const DistMat* A )
{ return A->LDim(); }

int CpxDistMatLDim( const CpxDistMat* A )
{ return A->LDim(); }

int DistMatColShift( const DistMat* A )
{ return A->ColShift(); }

int CpxDistMatColShift( const CpxDistMat* A )
{ return A->ColShift(); }

int DistMatRowShift( const DistMat* A )
{ return A->RowShift(); }

int CpxDistMatRowShift( const CpxDistMat* A )
{ return A->RowShift(); }

int DistMatColStride( const DistMat* A )
{ return A->ColStride(); }

int CpxDistMatColStride( const CpxDistMat* A )
{ return A->ColStride(); }

int DistMatRowStride( const DistMat* A )
{ return A->RowStride(); }

int CpxDistMatRowStride( const CpxDistMat* A )
{ return A->RowStride(); }

double GetDistMatEntry( const DistMat* A, int i, int j )
{ return A->Get( i, j ); }

void SetDistMatEntry( DistMat* A, int i, int j, double alpha )
{ A->Set( i, j, alpha ); }

double GetLocalDistMatEntry( const DistMat* A, int iLocal, int jLocal )
{ return A->GetLocal( iLocal, jLocal ); }

void SetLocalDistMatEntry( DistMat* A, int iLocal, int jLocal, double alpha )
{ A->SetLocal( iLocal, jLocal, alpha ); }

//
// [VC,* ] management
//

DistMat_VC_STAR* CreateDistMat_VC_STAR( const Grid* grid )
{ return new DistMat_VC_STAR( *grid ); }

CpxDistMat_VC_STAR* CreateCpxDistMat_VC_STAR( const Grid* grid )
{ return new CpxDistMat_VC_STAR( *grid ); }

void ResizeDistMat_VC_STAR( DistMat_VC_STAR* A, int height, int width )
{ A->ResizeTo( height, width ); }

void ResizeCpxDistMat_VC_STAR( CpxDistMat_VC_STAR* A, int height, int width )
{ A->ResizeTo( height, width ); }

void FreeDistMat_VC_STAR( DistMat_VC_STAR* A )
{ delete A; }

void FreeCpxDistMat_VC_STAR( CpxDistMat_VC_STAR* A )
{ delete A; }

void PrintDistMat_VC_STAR( const DistMat_VC_STAR* A )
{ A->Print(); }

void PrintCpxDistMat_VC_STAR( const CpxDistMat_VC_STAR* A )
{ A->Print(); }

int DistMatHeight_VC_STAR( const DistMat_VC_STAR* A )
{ return A->Height(); }

int CpxDistMatHeight_VC_STAR( const CpxDistMat_VC_STAR* A )
{ return A->Height(); }

int DistMatWidth_VC_STAR( const DistMat_VC_STAR* A )
{ return A->Width(); }

int CpxDistMatWidth_VC_STAR( const CpxDistMat_VC_STAR* A )
{ return A->Width(); }

int DistMatLocalHeight_VC_STAR( const DistMat_VC_STAR* A )
{ return A->LocalHeight(); }

int CpxDistMatLocalHeight_VC_STAR( const CpxDistMat_VC_STAR* A )
{ return A->LocalHeight(); }

int DistMatLocalWidth_VC_STAR( const DistMat_VC_STAR* A )
{ return A->LocalWidth(); }

int CpxDistMatLocalWidth_VC_STAR( const CpxDistMat_VC_STAR* A )
{ return A->LocalWidth(); }

int DistMatLDim_VC_STAR( const DistMat_VC_STAR* A )
{ return A->LDim(); }

int CpxDistMatLDim_VC_STAR( const CpxDistMat_VC_STAR* A )
{ return A->LDim(); }

int DistMatColShift_VC_STAR( const DistMat_VC_STAR* A )
{ return A->ColShift(); }

int CpxDistMatColShift_VC_STAR( const CpxDistMat_VC_STAR* A )
{ return A->ColShift(); }

int DistMatRowShift_VC_STAR( const DistMat_VC_STAR* A )
{ return A->RowShift(); }

int CpxDistMatRowShift_VC_STAR( const CpxDistMat_VC_STAR* A )
{ return A->RowShift(); }

int DistMatColStride_VC_STAR( const DistMat_VC_STAR* A )
{ return A->ColStride(); }

int CpxDistMatColStride_VC_STAR( const CpxDistMat_VC_STAR* A )
{ return A->ColStride(); }

int DistMatRowStride_VC_STAR( const DistMat_VC_STAR* A )
{ return A->RowStride(); }

int CpxDistMatRowStride_VC_STAR( const CpxDistMat_VC_STAR* A )
{ return A->RowStride(); }

double GetDistMatEntry_VC_STAR( const DistMat_VC_STAR* A, int i, int j )
{ return A->Get( i, j ); }

void SetDistMatEntry_VC_STAR( DistMat_VC_STAR* A, int i, int j, double alpha )
{ A->Set( i, j, alpha ); }

double GetLocalDistMatEntry_VC_STAR
( const DistMat_VC_STAR* A, int iLocal, int jLocal )
{ return A->GetLocal( iLocal, jLocal ); }

void SetLocalDistMatEntry_VC_STAR
( DistMat_VC_STAR* A, int iLocal, int jLocal, double alpha )
{ A->SetLocal( iLocal, jLocal, alpha ); }

//
// [VR,* ] management
//

DistMat_VR_STAR* CreateDistMat_VR_STAR( const Grid* grid )
{ return new DistMat_VR_STAR( *grid ); }

CpxDistMat_VR_STAR* CreateCpxDistMat_VR_STAR( const Grid* grid )
{ return new CpxDistMat_VR_STAR( *grid ); }

void ResizeDistMat_VR_STAR( DistMat_VR_STAR* A, int height, int width )
{ A->ResizeTo( height, width ); }

void ResizeCpxDistMat_VR_STAR( CpxDistMat_VR_STAR* A, int height, int width )
{ A->ResizeTo( height, width ); }

void FreeDistMat_VR_STAR( DistMat_VR_STAR* A )
{ delete A; }

void FreeCpxDistMat_VR_STAR( CpxDistMat_VR_STAR* A )
{ delete A; }

void PrintDistMat_VR_STAR( const DistMat_VR_STAR* A )
{ A->Print(); }

void PrintCpxDistMat_VR_STAR( const CpxDistMat_VR_STAR* A )
{ A->Print(); }

int DistMatHeight_VR_STAR( const DistMat_VR_STAR* A )
{ return A->Height(); }

int CpxDistMatHeight_VR_STAR( const CpxDistMat_VR_STAR* A )
{ return A->Height(); }

int DistMatWidth_VR_STAR( const DistMat_VR_STAR* A )
{ return A->Width(); }

int CpxDistMatWidth_VR_STAR( const CpxDistMat_VR_STAR* A )
{ return A->Width(); }

int DistMatLocalHeight_VR_STAR( const DistMat_VR_STAR* A )
{ return A->LocalHeight(); }

int CpxDistMatLocalHeight_VR_STAR( const CpxDistMat_VR_STAR* A )
{ return A->LocalHeight(); }

int DistMatLocalWidth_VR_STAR( const DistMat_VR_STAR* A )
{ return A->LocalWidth(); }

int CpxDistMatLocalWidth_VR_STAR( const CpxDistMat_VR_STAR* A )
{ return A->LocalWidth(); }

int DistMatLDim_VR_STAR( const DistMat_VR_STAR* A )
{ return A->LDim(); }

int CpxDistMatLDim_VR_STAR( const CpxDistMat_VR_STAR* A )
{ return A->LDim(); }

int DistMatColShift_VR_STAR( const DistMat_VR_STAR* A )
{ return A->ColShift(); }

int CpxDistMatColShift_VR_STAR( const CpxDistMat_VR_STAR* A )
{ return A->ColShift(); }

int DistMatRowShift_VR_STAR( const DistMat_VR_STAR* A )
{ return A->RowShift(); }

int CpxDistMatRowShift_VR_STAR( const CpxDistMat_VR_STAR* A )
{ return A->RowShift(); }

int DistMatColStride_VR_STAR( const DistMat_VR_STAR* A )
{ return A->ColStride(); }

int CpxDistMatColStride_VR_STAR( const CpxDistMat_VR_STAR* A )
{ return A->ColStride(); }

int DistMatRowStride_VR_STAR( const DistMat_VR_STAR* A )
{ return A->RowStride(); }

int CpxDistMatRowStride_VR_STAR( const CpxDistMat_VR_STAR* A )
{ return A->RowStride(); }

double GetDistMatEntry_VR_STAR( const DistMat_VR_STAR* A, int i, int j )
{ return A->Get( i, j ); }

void SetDistMatEntry_VR_STAR( DistMat_VR_STAR* A, int i, int j, double alpha )
{ A->Set( i, j, alpha ); }

double GetLocalDistMatEntry_VR_STAR
( const DistMat_VR_STAR* A, int iLocal, int jLocal )
{ return A->GetLocal( iLocal, jLocal ); }

void SetLocalDistMatEntry_VR_STAR
( DistMat_VR_STAR* A, int iLocal, int jLocal, double alpha )
{ A->SetLocal( iLocal, jLocal, alpha ); }

//
// QR factorization
//

void ExplicitQR( DistMat* A, DistMat* R )
{ elem::ExplicitQR( *A, *R ); }

void CpxExplicitQR( CpxDistMat* A, CpxDistMat* R )
{ elem::ExplicitQR( *A, *R ); }

//
// Singular Value Decomposition
//

void SVD( DistMat* A, DistMat_VR_STAR* s, DistMat* V )
{ elem::SVD( *A, *s, *V ); }

void CpxSVD( CpxDistMat* A, DistMat_VR_STAR* s, CpxDistMat* V )
{ elem::SVD( *A, *s, *V ); }

void SingularValues( DistMat* A, DistMat_VR_STAR* s );
void CpxSingularValues( CpxDistMat* A, DistMat_VR_STAR* s );

//
// Generalized Hermitian-definite eigensolvers for A X = B X \Lambda
//

void SymmetricAxBx
( DistMat* A, DistMat* B, DistMat_VR_STAR* w, DistMat* X )
{ HermitianGenDefiniteEig( AXBX, LOWER, *A, *B, *w, *X ); }

void SymmetricAxBxRange
( DistMat* A, DistMat* B, DistMat_VR_STAR* w, DistMat* X,
  double a, double b )
{ HermitianGenDefiniteEig( AXBX, LOWER, *A, *B, *w, *X, a, b ); }

void SymmetricAxBxIndices
( DistMat* A, DistMat* B, DistMat_VR_STAR* w, DistMat* X,
  int a, int b )
{ HermitianGenDefiniteEig( AXBX, LOWER, *A, *B, *w, *X, a, b ); }

void HermitianAxBx
( CpxDistMat* A, CpxDistMat* B, DistMat_VR_STAR* w, CpxDistMat* X )
{ HermitianGenDefiniteEig( AXBX, LOWER, *A, *B, *w, *X ); }

void HermitianAxBxRange
( CpxDistMat* A, CpxDistMat* B, DistMat_VR_STAR* w, CpxDistMat* X,
  double a, double b )
{ HermitianGenDefiniteEig( AXBX, LOWER, *A, *B, *w, *X, a, b ); }

void HermitianAxBxIndices
( CpxDistMat* A, CpxDistMat* B, DistMat_VR_STAR* w, CpxDistMat* X,
  int a, int b )
{ HermitianGenDefiniteEig( AXBX, LOWER, *A, *B, *w, *X, a, b ); }

//
// Special matrices
//

void UniformMat( Mat* A, int height, int width )
{ Uniform( height, width, *A ); }

void UniformCpxMat( CpxMat* A, int height, int width )
{ Uniform( height, width, *A ); }

void UniformDistMat( DistMat* A, int height, int width )
{ Uniform( height, width, *A ); }

void UniformCpxDistMat( CpxDistMat* A, int height, int width )
{ Uniform( height, width, *A ); }

//
// Utilities
//

int Length( int n, int shift, int modulus )
{ return elem::Length<int>(n,shift,modulus); }

} // extern "C"
