/*
   Copyright (c) 2011-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "elemental.hpp"
using namespace elem;

typedef int GridHandle;
typedef int RealDistMatHandle;
typedef int ComplexDistMatHandle;
typedef int RealDistColVecHandle;
typedef int ComplexDistColVecHandle;

namespace {

std::vector<Grid*> gridList;

std::vector<DistMatrix<double,MC,MR>*> realDistMatList;
std::vector<DistMatrix<Complex<double>,MC,MR>*> complexDistMatList;

std::vector<DistMatrix<double,VR,STAR>*> realDistColVecList;

Grid& TranslateGridHandle( GridHandle handle )
{ return *gridList[handle]; }

DistMatrix<double,MC,MR>& TranslateRealDistMatHandle
( RealDistMatHandle handle )
{ return *realDistMatList[handle]; }

DistMatrix<Complex<double>,MC,MR>& TranslateComplexDistMatHandle
( ComplexDistMatHandle handle )
{ return *complexDistMatList[handle]; }

DistMatrix<double,VR,STAR>& TranslateRealDistColVecHandle
( RealDistColVecHandle handle )
{ return *realDistColVecList[handle]; }

template<typename T>
unsigned GetOpenIndex( std::vector<T*>& list )
{
    unsigned index;
    for( index=0; index<list.size(); ++index )
        if( list[index] == 0 )
            break;
    if( index == list.size() )
        list.push_back( 0 );
    return index;
}

RealDistMatHandle CreateEmptyRealDistMat( const Grid& grid )
{
    const unsigned index = GetOpenIndex( realDistMatList );
    realDistMatList[index] = new DistMatrix<double,MC,MR>( grid );
    return index;
}

ComplexDistMatHandle CreateEmptyComplexDistMat( const Grid& grid )
{
    const unsigned index = GetOpenIndex( complexDistMatList );
    complexDistMatList[index] = new DistMatrix<Complex<double>,MC,MR>( grid );
    return index;
}

RealDistColVecHandle CreateEmptyRealDistColVec( const Grid& grid )
{
    const unsigned index = GetOpenIndex( realDistColVecList );
    realDistColVecList[index] = new DistMatrix<double,VR,STAR>( grid );
    return index;
}

} // anonymous namespace

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

GridHandle CreateGrid( MPI_Comm comm )
{
    const unsigned index = GetOpenIndex( gridList );
    gridList[index] = new Grid( comm );
    return index;
}

int GridHeight( GridHandle handle )
{ return gridList[handle]->Height(); }

int GridWidth( GridHandle handle )
{ return gridList[handle]->Width(); }

int GridSize( GridHandle handle )
{ return gridList[handle]->Size(); }

int GridRow( GridHandle handle )
{ return gridList[handle]->MCRank(); }

int GridCol( GridHandle handle )
{ return gridList[handle]->MRRank(); }

int GridRank( GridHandle handle )
{ return gridList[handle]->Rank(); }

void FreeGrid( GridHandle handle )
{
    if( gridList[handle] != 0 )
    {
        delete gridList[handle];
        gridList[handle] = 0;
    }
}

//
// Distributed matrix management
//

RealDistMatHandle RegisterRealDistMat
( int height, int width, int colAlignment, int rowAlignment, 
  double* buffer, int ldim, GridHandle gridHandle )
{
    const Grid& grid = TranslateGridHandle( gridHandle );
    const unsigned index = GetOpenIndex( realDistMatList );
    realDistMatList[index] = 
        new DistMatrix<double,MC,MR>
        (height,width,colAlignment,rowAlignment,buffer,ldim,grid);
    return index;
}

RealDistMatHandle CreateEmptyRealDistMat( GridHandle gridHandle )
{
    const Grid& grid = TranslateGridHandle( gridHandle );
    const unsigned index = GetOpenIndex( realDistMatList );
    realDistMatList[index] = new DistMatrix<double,MC,MR>( grid );
    return index;
}

ComplexDistMatHandle RegisterComplexDistMat
( int height, int width, int colAlignment, int rowAlignment, 
  void* voidBuffer, int ldim, GridHandle gridHandle )
{
    typedef Complex<double> C;
    C* buffer = static_cast<C*>(voidBuffer);

    const Grid& grid = TranslateGridHandle( gridHandle );
    const unsigned index = GetOpenIndex( complexDistMatList );
    complexDistMatList[index] = 
        new DistMatrix<Complex<double>,MC,MR>
        (height,width,colAlignment,rowAlignment,buffer,ldim,grid);
    return index;
}

ComplexDistMatHandle CreateEmptyComplexDistMat( GridHandle gridHandle )
{
    const Grid& grid = TranslateGridHandle( gridHandle );
    const unsigned index = GetOpenIndex( complexDistMatList );
    complexDistMatList[index] = new DistMatrix<Complex<double>,MC,MR>( grid );
    return index;
}

void FreeRealDistMat( RealDistMatHandle handle )
{
    if( realDistMatList[handle] != 0 )
    {
        delete realDistMatList[handle];
        realDistMatList[handle] = 0;
    }
}

void FreeComplexDistMat( ComplexDistMatHandle handle )
{
    if( complexDistMatList[handle] != 0 )
    {
        delete complexDistMatList[handle];
        complexDistMatList[handle] = 0;
    }
}

void PrintRealDistMat( RealDistMatHandle AHandle )
{
    const DistMatrix<double,MC,MR>& A = TranslateRealDistMatHandle( AHandle );
    A.Print();
}

void PrintComplexDistMat( ComplexDistMatHandle AHandle )
{
    typedef Complex<double> C;
    const DistMatrix<C,MC,MR>& A = TranslateComplexDistMatHandle( AHandle );
    A.Print();
}

//
// Distributed column vector management
//

RealDistColVecHandle CreateEmptyRealDistColVec( GridHandle gridHandle )
{
    const Grid& grid = TranslateGridHandle( gridHandle );
    const unsigned index = GetOpenIndex( realDistColVecList );
    realDistColVecList[index] = new DistMatrix<double,VR,STAR>( grid );
    return index;
}

void FreeRealDistColVec( RealDistColVecHandle handle )
{
    if( realDistColVecList[handle] != 0 )
    {
        delete realDistColVecList[handle];
        realDistColVecList[handle] = 0;
    }
}

void PrintRealDistColVec( RealDistColVecHandle AHandle )
{
    const DistMatrix<double,VR,STAR>& A = 
        TranslateRealDistColVecHandle( AHandle );
    A.Print();
}

//
// Generalized Hermitian-definite eigensolvers for A X = B X \Lambda
//

void SymmetricAxBx
( RealDistMatHandle AHandle, RealDistMatHandle BHandle,
  RealDistColVecHandle* wHandle, RealDistMatHandle* XHandle )
{
    DistMatrix<double,MC,MR>& A   = TranslateRealDistMatHandle( AHandle );
    DistMatrix<double,MC,MR>& B   = TranslateRealDistMatHandle( BHandle );

    *wHandle = CreateEmptyRealDistColVec( A.Grid() ); 
    *XHandle = CreateEmptyRealDistMat( A.Grid() );
    DistMatrix<double,VR,STAR>& w = TranslateRealDistColVecHandle( *wHandle );
    DistMatrix<double,MC,MR>& X   = TranslateRealDistMatHandle( *XHandle );
    
    HermitianGenDefiniteEig( AXBX, LOWER, A, B, w, X );
}

void SymmetricAxBxPartialRange
( RealDistMatHandle AHandle, RealDistMatHandle BHandle,
  RealDistColVecHandle* wHandle, RealDistMatHandle* XHandle,
  double a, double b )
{
    DistMatrix<double,MC,MR>& A   = TranslateRealDistMatHandle( AHandle );
    DistMatrix<double,MC,MR>& B   = TranslateRealDistMatHandle( BHandle );
    
    *wHandle = CreateEmptyRealDistColVec( A.Grid() ); 
    *XHandle = CreateEmptyRealDistMat( A.Grid() );
    DistMatrix<double,VR,STAR>& w = TranslateRealDistColVecHandle( *wHandle );
    DistMatrix<double,MC,MR>& X   = TranslateRealDistMatHandle( *XHandle );
    
    HermitianGenDefiniteEig( AXBX, LOWER, A, B, w, X, a, b );
}

void SymmetricAxBxPartialIndices
( RealDistMatHandle AHandle, RealDistMatHandle BHandle,
  RealDistColVecHandle* wHandle, RealDistMatHandle* XHandle,
  int a, int b )
{
    DistMatrix<double,MC,MR>& A   = TranslateRealDistMatHandle( AHandle );
    DistMatrix<double,MC,MR>& B   = TranslateRealDistMatHandle( BHandle );

    *wHandle = CreateEmptyRealDistColVec( A.Grid() ); 
    *XHandle = CreateEmptyRealDistMat( A.Grid() );
    DistMatrix<double,VR,STAR>& w = TranslateRealDistColVecHandle( *wHandle );
    DistMatrix<double,MC,MR>& X   = TranslateRealDistMatHandle( *XHandle );
    
    HermitianGenDefiniteEig( AXBX, LOWER, A, B, w, X, a, b );
}

void HermitianAxBx
( ComplexDistMatHandle AHandle, ComplexDistMatHandle BHandle,
  RealDistColVecHandle* wHandle, ComplexDistMatHandle* XHandle )
{
    typedef Complex<double> C;

    DistMatrix<C,MC,MR>& A        = TranslateComplexDistMatHandle( AHandle );
    DistMatrix<C,MC,MR>& B        = TranslateComplexDistMatHandle( BHandle );
    
    *wHandle = CreateEmptyRealDistColVec( A.Grid() ); 
    *XHandle = CreateEmptyComplexDistMat( A.Grid() );
    DistMatrix<double,VR,STAR>& w = TranslateRealDistColVecHandle( *wHandle );
    DistMatrix<C,MC,MR>& X        = TranslateComplexDistMatHandle( *XHandle );
    
    HermitianGenDefiniteEig( AXBX, LOWER, A, B, w, X );
}

void HermitianAxBxPartialRange
( ComplexDistMatHandle AHandle, ComplexDistMatHandle BHandle,
  RealDistColVecHandle* wHandle, ComplexDistMatHandle* XHandle,
  double a, double b )
{
    typedef Complex<double> C;

    DistMatrix<C,MC,MR>& A        = TranslateComplexDistMatHandle( AHandle );
    DistMatrix<C,MC,MR>& B        = TranslateComplexDistMatHandle( BHandle );

    *wHandle = CreateEmptyRealDistColVec( A.Grid() ); 
    *XHandle = CreateEmptyComplexDistMat( A.Grid() );
    DistMatrix<double,VR,STAR>& w = TranslateRealDistColVecHandle( *wHandle );
    DistMatrix<C,MC,MR>& X        = TranslateComplexDistMatHandle( *XHandle );
    
    HermitianGenDefiniteEig( AXBX, LOWER, A, B, w, X, a, b );
}

void HermitianAxBxPartialIndices
( ComplexDistMatHandle AHandle, ComplexDistMatHandle BHandle,
  RealDistColVecHandle* wHandle, ComplexDistMatHandle* XHandle,
  int a, int b )
{
    typedef Complex<double> C;

    DistMatrix<C,MC,MR>& A        = TranslateComplexDistMatHandle( AHandle );
    DistMatrix<C,MC,MR>& B        = TranslateComplexDistMatHandle( BHandle );

    *wHandle = CreateEmptyRealDistColVec( A.Grid() ); 
    *XHandle = CreateEmptyComplexDistMat( A.Grid() );
    DistMatrix<double,VR,STAR>& w = TranslateRealDistColVecHandle( *wHandle );
    DistMatrix<C,MC,MR>& X        = TranslateComplexDistMatHandle( *XHandle );
    
    HermitianGenDefiniteEig( AXBX, LOWER, A, B, w, X, a, b );
}

//
// Utilities
//

int LocalLength( int n, int shift, int modulus )
{ return elem::LocalLength<int>(n,shift,modulus); }

} // extern "C"
