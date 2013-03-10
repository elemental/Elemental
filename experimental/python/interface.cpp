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

extern "C" {
#include "./elemental.h"
}

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

void FreeGrid( Grid** grid )
{
    delete *grid;
    *grid = 0; 
}

//
// Distributed matrix management
//

DistMatrix<double>* CreateDistMat( const Grid* grid )
{ return new DistMatrix<double>( *grid ); }

DistMatrix<Complex<double> >* CreateCpxDistMat( const Grid* grid )
{ return new DistMatrix<Complex<double> >( *grid ); }

void AttachToDistMat
( DistMatrix<double>* A, 
  int height, int width, int colAlignment, int rowAlignment, 
  void* buffer, int ldim, const Grid* grid )
{ A->Attach
  ( height, width, colAlignment, rowAlignment, 
    reinterpret_cast<double*>(buffer), ldim, *grid ); }

void AttachToCpxDistMat
( DistMatrix<Complex<double> >* A, 
  int height, int width, int colAlignment, int rowAlignment, 
  void* buffer, int ldim, const Grid* grid )
{ A->Attach
  ( height, width, colAlignment, rowAlignment, 
    reinterpret_cast<Complex<double>*>(buffer), ldim, *grid ); }

void UniformDistMat( DistMatrix<double>* A, int height, int width )
{ Uniform( height, width, *A ); }

void UniformCpxDistMat( DistMatrix<Complex<double> >* A, int height, int width )
{ Uniform( height, width, *A ); }

void FreeDistMat( DistMatrix<double>** A )
{
    delete *A;
    *A = 0;
}

void FreeCpxDistMat( DistMatrix<Complex<double> >** A )
{
    delete *A;
    *A = 0;
}

void PrintDistMat( const DistMatrix<double>* A )
{ A->Print(); }

void PrintCpxDistMat( const DistMatrix<Complex<double> >* A )
{ A->Print(); }

//
// [VC,* ] management
//

DistMatrix<double,VC,STAR>* CreateDistMat_VC_STAR( const Grid* grid )
{ return new DistMatrix<double,VC,STAR>( *grid ); }

DistMatrix<Complex<double>,VC,STAR>* 
CreateCpxDistMat_VC_STAR( const Grid* grid )
{ return new DistMatrix<Complex<double>,VC,STAR>( *grid ); }

void FreeDistMat_VC_STAR( DistMatrix<double,VC,STAR>** A )
{
    delete *A;
    *A = 0;
}

void FreeCpxDistMat_VC_STAR( DistMatrix<Complex<double>,VC,STAR>** A )
{
    delete *A;
    *A = 0;
}

void PrintDistMat_VC_STAR( const DistMatrix<double,VC,STAR>* A )
{ A->Print(); }

void PrintCpxDistMat_VC_STAR( const DistMatrix<Complex<double>,VC,STAR>* A )
{ A->Print(); }

//
// [VR,* ] management
//

DistMatrix<double,VR,STAR>* CreateDistMat_VR_STAR( const Grid* grid )
{ return new DistMatrix<double,VR,STAR>( *grid ); }

DistMatrix<Complex<double>,VR,STAR>* 
CreateCpxDistMat_VR_STAR( const Grid* grid )
{ return new DistMatrix<Complex<double>,VR,STAR>( *grid ); }

void FreeDistMat_VR_STAR( DistMatrix<double,VR,STAR>** A )
{
    delete *A;
    *A = 0;
}

void FreeCpxDistMat_VR_STAR( DistMatrix<Complex<double>,VR,STAR>** A )
{
    delete *A;
    *A = 0;
}

void PrintDistMat_VR_STAR( const DistMatrix<double,VR,STAR>* A )
{ A->Print(); }

void PrintCpxDistMat_VR_STAR( const DistMatrix<Complex<double>,VR,STAR>* A )
{ A->Print(); }

//
// Singular Value Decomposition
//

void SVD
( DistMatrix<double>* A, DistMatrix<double,VR,STAR>* s, DistMatrix<double>* V )
{
    s->SetGrid( A->Grid() );
    V->SetGrid( A->Grid() );
    elem::SVD( *A, *s, *V );
}

void CpxSVD
( DistMatrix<Complex<double> >* A, 
  DistMatrix<double,VR,STAR>* s, 
  DistMatrix<Complex<double> >* V )
{
    s->SetGrid( A->Grid() );
    V->SetGrid( A->Grid() );
    elem::SVD( *A, *s, *V );
}

//
// Generalized Hermitian-definite eigensolvers for A X = B X \Lambda
//

void SymmetricAxBx
( DistMatrix<double>* A, DistMatrix<double>* B,
  DistMatrix<double,VR,STAR>* w, DistMatrix<double>* X )
{
    w->SetGrid( A->Grid() );
    X->SetGrid( A->Grid() );
    HermitianGenDefiniteEig( AXBX, LOWER, *A, *B, *w, *X );
}

void SymmetricAxBxRange
( DistMatrix<double>* A, DistMatrix<double>* B,
  DistMatrix<double,VR,STAR>* w, DistMatrix<double>* X,
  double a, double b )
{
    w->SetGrid( A->Grid() );
    X->SetGrid( A->Grid() );
    HermitianGenDefiniteEig( AXBX, LOWER, *A, *B, *w, *X, a, b );
}

void SymmetricAxBxIndices
( DistMatrix<double>* A, DistMatrix<double>* B,
  DistMatrix<double,VR,STAR>* w, DistMatrix<double>* X,
  int a, int b )
{
    w->SetGrid( A->Grid() );
    X->SetGrid( A->Grid() );
    HermitianGenDefiniteEig( AXBX, LOWER, *A, *B, *w, *X, a, b );
}

void HermitianAxBx
( DistMatrix<Complex<double> >* A, DistMatrix<Complex<double> >* B,
  DistMatrix<double,VR,STAR>* w, DistMatrix<Complex<double> >* X )
{
    w->SetGrid( A->Grid() );
    X->SetGrid( A->Grid() );
    HermitianGenDefiniteEig( AXBX, LOWER, *A, *B, *w, *X );
}

void HermitianAxBxRange
( DistMatrix<Complex<double> >* A, DistMatrix<Complex<double> >* B,
  DistMatrix<double,VR,STAR>* w, DistMatrix<Complex<double> >* X,
  double a, double b )
{
    w->SetGrid( A->Grid() );
    X->SetGrid( A->Grid() );
    HermitianGenDefiniteEig( AXBX, LOWER, *A, *B, *w, *X, a, b );
}

void HermitianAxBxIndices
( DistMatrix<Complex<double> >* A, DistMatrix<Complex<double> >* B,
  DistMatrix<double,VR,STAR>* w, DistMatrix<Complex<double> >* X,
  int a, int b )
{
    w->SetGrid( A->Grid() );
    X->SetGrid( A->Grid() );
    HermitianGenDefiniteEig( AXBX, LOWER, *A, *B, *w, *X, a, b );
}

//
// Utilities
//

int Length( int n, int shift, int modulus )
{ return elem::Length<int>(n,shift,modulus); }

} // extern "C"
