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
#include "elemental.hpp"
using namespace elemental;

typedef int GridHandle;
typedef int RealDistMatHandle;
typedef int ComplexDistMatHandle;
typedef int RealDistColVecHandle;
typedef int ComplexDistColVecHandle;

namespace {

std::vector<Grid*> gridList;

std::vector<DistMatrix<double,MC,MR>*> realDistMatList;
std::vector<DistMatrix<std::complex<double>,MC,MR>*> complexDistMatList;

std::vector<DistMatrix<double,VR,STAR>*> realDistColVecList;

Grid& TranslateGridHandle( GridHandle handle )
{ return *gridList[handle]; }

DistMatrix<double,MC,MR>& TranslateRealDistMatHandle
( RealDistMatHandle realDistMatHandle )
{ return *realDistMatList[realDistMatHandle]; }

DistMatrix<std::complex<double>,MC,MR>& TranslateComplexDistMatHandle
( ComplexDistMatHandle complexDistMatHandle )
{ return *complexDistMatList[complexDistMatHandle]; }

DistMatrix<double,VR,STAR>& TranslateRealDistColVecHandle
( RealDistColVecHandle realDistColVecHandle )
{ return *realDistColVecList[realDistColVecHandle]; }

RealDistMatHandle CreateEmptyRealDistMat( const Grid& grid )
{
    int index;
    for( index=0; index<realDistMatList.size(); ++index )
    {
        if( realDistMatList[index] == 0 )
        {
            realDistMatList[index] = new DistMatrix<double,MC,MR>(grid);
            break;
        }
    }

    if( index == realDistMatList.size() )
        realDistMatList.push_back( new DistMatrix<double,MC,MR>(grid) );
    return index;
}

ComplexDistMatHandle CreateEmptyComplexDistMat( const Grid& grid )
{
    int index;
    for( index=0; index<complexDistMatList.size(); ++index )
    {
        if( complexDistMatList[index] == 0 )
        {
            complexDistMatList[index] =
                new DistMatrix<std::complex<double>,MC,MR>(grid);
            break;
        }
    }

    if( index == complexDistMatList.size() )
        complexDistMatList.push_back(
            new DistMatrix<std::complex<double>,MC,MR>(grid)
        );
    return index;
}

RealDistColVecHandle CreateEmptyRealDistColVec( const Grid& grid )
{
    int index;
    for( index=0; index<realDistColVecList.size(); ++index )
    {
        if( realDistColVecList[index] == 0 )
        {
            realDistColVecList[index] = new DistMatrix<double,VR,STAR>(grid);
            break;
        }
    }

    if( index == realDistColVecList.size() )
        realDistColVecList.push_back( new DistMatrix<double,VR,STAR>(grid) );
    return index;
}

}

extern "C" {

void Initialize( int* argc, char** argv[] )
{
    elemental::Initialize( *argc, *argv );
}

void Finalize()
{
    elemental::Finalize();
}

int LocalLength( int n, int shift, int modulus )
{
    return elemental::LocalLength<int>(n,shift,modulus);
}

GridHandle CreateGrid( MPI_Comm comm )
{
    // TODO: Switch from O(n) algorithm to O(1) algorithm by keeping static 
    //       data which tracks the last free index (if it exists)
    unsigned index;
    for( index=0; index<gridList.size(); ++index )
    {
        if( gridList[index] == 0 )
        {
            gridList[index] = new Grid( comm );        
            break;
        }
    }

    if( index == gridList.size() )
        gridList.push_back( new Grid( comm ) );
    return index;
}

int GridHeight( GridHandle handle )
{
    return gridList[handle]->Height();
}

int GridWidth( GridHandle handle )
{
    return gridList[handle]->Width();
}

int GridSize( GridHandle handle )
{
    return gridList[handle]->Size();
}

int GridRow( GridHandle handle )
{
    return gridList[handle]->MCRank();
}

int GridCol( GridHandle handle )
{
    return gridList[handle]->MRRank();
}

int GridRank( GridHandle handle )
{
    return gridList[handle]->Rank();
}

RealDistMatHandle RegisterRealDistMat
( int height, int width, int colAlignment, int rowAlignment, 
  double* buffer, int ldim, GridHandle gridHandle )
{
    const Grid& grid = TranslateGridHandle( gridHandle );

    unsigned index;
    for( index=0; index<realDistMatList.size(); ++index )
    {
        if( realDistMatList[index] == 0 )
        {
            realDistMatList[index] = 
                new DistMatrix<double,MC,MR>
                (height,width,colAlignment,rowAlignment,buffer,ldim,grid);
            break;
        }
    }

    if( index == realDistMatList.size() )
    {
        realDistMatList.push_back( 
            new DistMatrix<double,MC,MR>
            (height,width,colAlignment,rowAlignment,buffer,ldim,grid)
        );
    }
    return index;
}

RealDistMatHandle CreateEmptyRealDistMat( GridHandle gridHandle )
{
    const Grid& grid = TranslateGridHandle( gridHandle );

    unsigned index;
    for( index=0; index<realDistMatList.size(); ++index )
    {
        if( realDistMatList[index] == 0 )
        {
            realDistMatList[index] = new DistMatrix<double,MC,MR>(grid);
            break;
        }
    }

    if( index == realDistMatList.size() )
        realDistMatList.push_back( new DistMatrix<double,MC,MR>(grid) );
    return index;
}

ComplexDistMatHandle RegisterComplexDistMat
( int height, int width, int colAlignment, int rowAlignment, 
  void* voidBuffer, int ldim, GridHandle gridHandle )
{
    typedef std::complex<double> C;
    const Grid& grid = TranslateGridHandle( gridHandle );
    C* buffer = static_cast<C*>(voidBuffer);

    unsigned index;
    for( index=0; index<complexDistMatList.size(); ++index )
    {
        if( complexDistMatList[index] == 0 )
        {
            complexDistMatList[index] = 
                new DistMatrix<std::complex<double>,MC,MR>
                (height,width,colAlignment,rowAlignment,buffer,ldim,grid);
            break;
        }
    }

    if( index == complexDistMatList.size() )
    {
        complexDistMatList.push_back( 
            new DistMatrix<std::complex<double>,MC,MR>
            (height,width,colAlignment,rowAlignment,buffer,ldim,grid)
        );
    }
    return index;
}

ComplexDistMatHandle CreateEmptyComplexDistMat( GridHandle gridHandle )
{
    const Grid& grid = TranslateGridHandle( gridHandle );

    unsigned index;
    for( index=0; index<complexDistMatList.size(); ++index )
    {
        if( complexDistMatList[index] == 0 )
        {
            complexDistMatList[index] = 
                new DistMatrix<std::complex<double>,MC,MR>(grid);
            break;
        }
    }

    if( index == complexDistMatList.size() )
    {
        complexDistMatList.push_back( 
                new DistMatrix<std::complex<double>,MC,MR>(grid) 
        );
    }
    return index;
}

RealDistColVecHandle CreateEmptyRealDistColVec( GridHandle gridHandle )
{
    const Grid& grid = TranslateGridHandle( gridHandle );

    unsigned index;
    for( index=0; index<realDistColVecList.size(); ++index )
    {
        if( realDistColVecList[index] == 0 )
        {
            realDistColVecList[index] = new DistMatrix<double,VR,STAR>(grid);
            break;
        }
    }

    if( index == realDistColVecList.size() )
        realDistColVecList.push_back( new DistMatrix<double,VR,STAR>(grid) );
    return index;
}

void FreeGrid( GridHandle handle )
{
    if( gridList[handle] != 0 )
    {
        delete gridList[handle];
        gridList[handle] = 0;
    }
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

void FreeRealDistColVec( RealDistColVecHandle handle )
{
    if( realDistColVecList[handle] != 0 )
    {
        delete realDistColVecList[handle];
        realDistColVecList[handle] = 0;
    }
}

void PrintRealDistMat( RealDistMatHandle AHandle )
{
    const DistMatrix<double,MC,MR>& A = TranslateRealDistMatHandle( AHandle );
    A.Print();
}

void PrintComplexDistMat( ComplexDistMatHandle AHandle )
{
    typedef std::complex<double> C;
    const DistMatrix<C,MC,MR>& A = TranslateComplexDistMatHandle( AHandle );
    A.Print();
}

void PrintRealDistColVec( RealDistColVecHandle AHandle )
{
    const DistMatrix<double,VR,STAR>& A = 
        TranslateRealDistColVecHandle( AHandle );
    A.Print();
}

void SymmetricAxBx
( RealDistMatHandle AHandle, RealDistMatHandle BHandle,
  RealDistColVecHandle* wHandle, RealDistMatHandle* XHandle )
{
    DistMatrix<double,MC,MR>& A   = TranslateRealDistMatHandle( AHandle );
    DistMatrix<double,MC,MR>& B   = TranslateRealDistMatHandle( BHandle );

    *wHandle = CreateEmptyRealDistColVec( A.Grid() ); 
    *XHandle = CreateEmptyRealDistMat( A.Grid() );
    DistMatrix<double,VR,STAR>& w = TranslateRealDistColVecHandle( wHandle );
    DistMatrix<double,MC,MR>& X   = TranslateRealDistMatHandle( XHandle );
    
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
    DistMatrix<double,VR,STAR>& w = TranslateRealDistColVecHandle( wHandle );
    DistMatrix<double,MC,MR>& X   = TranslateRealDistMatHandle( XHandle );
    
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
    DistMatrix<double,VR,STAR>& w = TranslateRealDistColVecHandle( wHandle );
    DistMatrix<double,MC,MR>& X   = TranslateRealDistMatHandle( XHandle );
    
    HermitianGenDefiniteEig( AXBX, LOWER, A, B, w, X, a, b );
}

void HermitianAxBx
( ComplexDistMatHandle AHandle, ComplexDistMatHandle BHandle,
  RealDistColVecHandle* wHandle, ComplexDistMatHandle* XHandle )
{
    typedef std::complex<double> C;

    DistMatrix<C,MC,MR>& A        = TranslateComplexDistMatHandle( AHandle );
    DistMatrix<C,MC,MR>& B        = TranslateComplexDistMatHandle( BHandle );
    
    *wHandle = CreateEmptyRealDistColVec( A.Grid() ); 
    *XHandle = CreateEmptyComplexDistMat( A.Grid() );
    DistMatrix<double,VR,STAR>& w = TranslateRealDistColVecHandle( wHandle );
    DistMatrix<C,MC,MR>& X        = TranslateComplexDistMatHandle( XHandle );
    
    HermitianGenDefiniteEig( AXBX, LOWER, A, B, w, X );
}

void HermitianAxBxPartialRange
( ComplexDistMatHandle AHandle, ComplexDistMatHandle BHandle,
  RealDistColVecHandle* wHandle, ComplexDistMatHandle* XHandle,
  double a, double b )
{
    typedef std::complex<double> C;

    DistMatrix<C,MC,MR>& A        = TranslateComplexDistMatHandle( AHandle );
    DistMatrix<C,MC,MR>& B        = TranslateComplexDistMatHandle( BHandle );

    *wHandle = CreateEmptyRealDistColVec( A.Grid() ); 
    *XHandle = CreateEmptyComplexDistMat( A.Grid() );
    DistMatrix<double,VR,STAR>& w = TranslateRealDistColVecHandle( wHandle );
    DistMatrix<C,MC,MR>& X        = TranslateComplexDistMatHandle( XHandle );
    
    HermitianGenDefiniteEig( AXBX, LOWER, A, B, w, X, a, b );
}

void HermitianAxBxPartialIndices
( ComplexDistMatHandle AHandle, ComplexDistMatHandle BHandle,
  RealDistColVecHandle* wHandle, ComplexDistMatHandle* XHandle,
  int a, int b )
{
    typedef std::complex<double> C;

    DistMatrix<C,MC,MR>& A        = TranslateComplexDistMatHandle( AHandle );
    DistMatrix<C,MC,MR>& B        = TranslateComplexDistMatHandle( BHandle );

    *wHandle = CreateEmptyRealDistColVec( A.Grid() ); 
    *XHandle = CreateEmptyComplexDistMat( A.Grid() );
    DistMatrix<double,VR,STAR>& w = TranslateRealDistColVecHandle( wHandle );
    DistMatrix<C,MC,MR>& X        = TranslateComplexDistMatHandle( XHandle );
    
    HermitianGenDefiniteEig( AXBX, LOWER, A, B, w, X, a, b );
}

} // extern "C"
