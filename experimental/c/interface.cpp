/*
   Copyright (c) 2011-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <limits>
#include "El.hpp"
// TODO: Switch to El-lite.hpp
using namespace El;

extern "C" {
#include "./elemental.h"
}

namespace {

typedef double Real;
typedef Complex<Real> Cpx;

const unsigned maxUnsigned = std::numeric_limits<unsigned>::max();

std::vector<Grid*> gridList;

std::vector<DistMatrix<Real>*> distMatList;
std::vector<DistMatrix<Cpx >*> cpxDistMatList;

std::vector<DistMatrix<Real,VC,STAR>*> distMatList_VC_STAR;
std::vector<DistMatrix<Real,VR,STAR>*> distMatList_VR_STAR;

const Grid& GetGrid( ElemGrid grid )
{ 
    if( grid == maxUnsigned )
        return DefaultGrid();
    else
        return *gridList[grid]; 
}

DistMatrix<Real>& GetDistMat( ElemDistMat A )
{ return *distMatList[A]; }

DistMatrix<Cpx>& GetCpxDistMat( ElemCpxDistMat A )
{ return *cpxDistMatList[A]; }

DistMatrix<Real,VC,STAR>& GetDistMat_VC_STAR( ElemDistMat_VC_STAR A )
{ return *distMatList_VC_STAR[A]; }

DistMatrix<Real,VR,STAR>& GetDistMat_VR_STAR( ElemDistMat_VR_STAR A )
{ return *distMatList_VR_STAR[A]; }

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

ElemDistMat CreateDistMat( const Grid& grid )
{
    const unsigned index = GetOpenIndex( distMatList );
    distMatList[index] = new DistMatrix<Real>( grid );
    return index;
}

ElemCpxDistMat CreateCpxDistMat( const Grid& grid )
{
    const unsigned index = GetOpenIndex( cpxDistMatList );
    cpxDistMatList[index] = new DistMatrix<Cpx>( grid );
    return index;
}

ElemDistMat_VC_STAR CreateDistMat_VC_STAR( const Grid& grid )
{
    const unsigned index = GetOpenIndex( distMatList_VC_STAR );
    distMatList_VC_STAR[index] = new DistMatrix<Real,VC,STAR>( grid );
    return index;
}

ElemDistMat_VR_STAR CreateDistMat_VR_STAR( const Grid& grid )
{
    const unsigned index = GetOpenIndex( distMatList_VR_STAR );
    distMatList_VR_STAR[index] = new DistMatrix<Real,VR,STAR>( grid );
    return index;
}

} // anonymous namespace

extern "C" {

//
// Environment controls
//

void ElemInitialize( int* argc, char** argv[] )
{ 
    for( unsigned j=0; j<gridList.size(); ++j )
        ElemFreeGrid( ElemGrid(j) );
    gridList.clear();
    Initialize( *argc, *argv ); 
}

void ElemFinalize()
{ 
    for( unsigned j=0; j<gridList.size(); ++j )
        ElemFreeGrid( ElemGrid(j) );
    gridList.clear();
    Finalize(); 
}

void ElemSetBlocksize( int blocksize )
{ SetBlocksize( blocksize ); }

int ElemBlocksize()
{ return Blocksize(); }

//
// Process grid management
//

ElemGrid ElemDefaultGrid()
{ return maxUnsigned; }

ElemGrid ElemCreateGrid( MPI_Comm comm )
{
    const unsigned index = GetOpenIndex( gridList );
    gridList[index] = new Grid( comm );
    return index;
}

int ElemGridHeight( ElemGrid grid )
{ return gridList[grid]->Height(); }

int ElemGridWidth( ElemGrid grid )
{ return gridList[grid]->Width(); }

int ElemGridSize( ElemGrid grid )
{ return gridList[grid]->Size(); }

int ElemGridRow( ElemGrid grid )
{ return gridList[grid]->MCRank(); }

int ElemGridCol( ElemGrid grid )
{ return gridList[grid]->MRRank(); }

int ElemGridRank( ElemGrid grid )
{ return gridList[grid]->Rank(); }

void ElemFreeGrid( ElemGrid grid )
{
    if( grid != maxUnsigned && gridList[grid] != 0 )
    {
        delete gridList[grid];
        gridList[grid] = 0;
    }
}

//
// Distributed matrix management
//

ElemDistMat ElemRegisterDistMat
( int height, int width, int colAlign, int rowAlign, 
  Real* buffer, int ldim, ElemGrid gridHandle )
{
    const Grid& grid = GetGrid( gridHandle );
    const unsigned index = GetOpenIndex( distMatList );
    distMatList[index] = new DistMatrix<Real>(grid);
    distMatList[index]->Attach(height,width,grid,colAlign,rowAlign,buffer,ldim);
    return index;
}

ElemDistMat ElemCreateDistMat( ElemGrid gridHandle )
{
    const Grid& grid = GetGrid( gridHandle );
    const unsigned index = GetOpenIndex( distMatList );
    distMatList[index] = new DistMatrix<Real>( grid );
    return index;
}

ElemCpxDistMat ElemRegisterCpxDistMat
( int height, int width, int colAlign, int rowAlign, 
  void* voidBuffer, int ldim, ElemGrid gridHandle )
{
    Cpx* buffer = static_cast<Cpx*>(voidBuffer);

    const Grid& grid = GetGrid( gridHandle );
    const unsigned index = GetOpenIndex( cpxDistMatList );
    cpxDistMatList[index] = new DistMatrix<Cpx>(grid);
    cpxDistMatList[index]->Attach
    (height,width,grid,colAlign,rowAlign,buffer,ldim);
    return index;
}

ElemCpxDistMat ElemCreateCpxDistMat( ElemGrid gridHandle )
{
    const Grid& grid = GetGrid( gridHandle );
    const unsigned index = GetOpenIndex( cpxDistMatList );
    cpxDistMatList[index] = new DistMatrix<Cpx>( grid );
    return index;
}

void ElemUniformDistMat( ElemDistMat AHandle, int height, int width )
{
    auto& A = GetDistMat( AHandle );
    Uniform( A, height, width );
}

void ElemUniformCpxDistMat( ElemCpxDistMat AHandle, int height, int width )
{
    auto& A = GetCpxDistMat( AHandle );
    Uniform( A, height, width );
}

void ElemFreeDistMat( ElemDistMat A )
{
    if( distMatList[A] != 0 )
    {
        delete distMatList[A];
        distMatList[A] = 0;
    }
}

void ElemFreeCpxDistMat( ElemCpxDistMat A )
{
    if( cpxDistMatList[A] != 0 )
    {
        delete cpxDistMatList[A];
        cpxDistMatList[A] = 0;
    }
}

void ElemPrintDistMat( ElemDistMat AHandle )
{
    const auto& A = GetDistMat( AHandle );
    Print( A );
}

void ElemPrintCpxDistMat( ElemCpxDistMat AHandle )
{
    const auto& A = GetCpxDistMat( AHandle );
    Print( A );
}

//
// [VC,* ] management
//

ElemDistMat_VC_STAR ElemCreateDistMat_VC_STAR( ElemGrid gridHandle )
{
    const Grid& grid = GetGrid( gridHandle );
    const unsigned index = GetOpenIndex( distMatList_VC_STAR );
    distMatList_VC_STAR[index] = new DistMatrix<Real,VC,STAR>( grid );
    return index;
}

void ElemFreeDistMat_VC_STAR( ElemDistMat_VC_STAR A )
{
    if( distMatList_VC_STAR[A] != 0 )
    {
        delete distMatList_VC_STAR[A];
        distMatList_VC_STAR[A] = 0;
    }
}

void ElemPrintDistMat_VC_STAR( ElemDistMat_VC_STAR AHandle )
{
    const auto& A = GetDistMat_VC_STAR( AHandle );
    Print( A );
}

//
// [VR,* ] management
//

ElemDistMat_VR_STAR ElemCreateDistMat_VR_STAR( ElemGrid gridHandle )
{
    const Grid& grid = GetGrid( gridHandle );
    const unsigned index = GetOpenIndex( distMatList_VR_STAR );
    distMatList_VR_STAR[index] = new DistMatrix<Real,VR,STAR>( grid );
    return index;
}

void ElemFreeDistMat_VR_STAR( ElemDistMat_VR_STAR A )
{
    if( distMatList_VR_STAR[A] != 0 )
    {
        delete distMatList_VR_STAR[A];
        distMatList_VR_STAR[A] = 0;
    }
}

void ElemPrintDistMat_VR_STAR( ElemDistMat_VR_STAR AHandle )
{
    const auto& A = GetDistMat_VR_STAR( AHandle );
    Print( A );
}

//
// Generalized Hermitian-definite eigensolvers for A X = B X \Lambda
//

void ElemSymmetricAxBx
( ElemDistMat AHandle, ElemDistMat BHandle,
  ElemDistMat_VR_STAR* wHandle, ElemDistMat* XHandle )
{
    auto& A = GetDistMat( AHandle );
    auto& B = GetDistMat( BHandle );

    *wHandle = CreateDistMat_VR_STAR( A.Grid() ); 
    *XHandle = CreateDistMat( A.Grid() );
    auto& w = GetDistMat_VR_STAR( *wHandle );
    auto& X = GetDistMat( *XHandle );
    
    HermitianGenDefiniteEig( AXBX, LOWER, A, B, w, X );
}

void ElemSymmetricAxBxRange
( ElemDistMat AHandle, ElemDistMat BHandle,
  ElemDistMat_VR_STAR* wHandle, ElemDistMat* XHandle,
  Real a, Real b )
{
    auto& A = GetDistMat( AHandle );
    auto& B = GetDistMat( BHandle );
    
    *wHandle = CreateDistMat_VR_STAR( A.Grid() ); 
    *XHandle = CreateDistMat( A.Grid() );
    auto& w = GetDistMat_VR_STAR( *wHandle );
    auto& X = GetDistMat( *XHandle );
    
    HermitianGenDefiniteEig( AXBX, LOWER, A, B, w, X, a, b );
}

void ElemSymmetricAxBxIndices
( ElemDistMat AHandle, ElemDistMat BHandle,
  ElemDistMat_VR_STAR* wHandle, ElemDistMat* XHandle,
  int a, int b )
{
    auto& A = GetDistMat( AHandle );
    auto& B = GetDistMat( BHandle );

    *wHandle = CreateDistMat_VR_STAR( A.Grid() ); 
    *XHandle = CreateDistMat( A.Grid() );
    auto& w = GetDistMat_VR_STAR( *wHandle );
    auto& X = GetDistMat( *XHandle );
    
    HermitianGenDefiniteEig( AXBX, LOWER, A, B, w, X, a, b );
}

void ElemHermitianAxBx
( ElemCpxDistMat AHandle, ElemCpxDistMat BHandle,
  ElemDistMat_VR_STAR* wHandle, ElemCpxDistMat* XHandle )
{
    auto& A = GetCpxDistMat( AHandle );
    auto& B = GetCpxDistMat( BHandle );
    
    *wHandle = CreateDistMat_VR_STAR( A.Grid() ); 
    *XHandle = CreateCpxDistMat( A.Grid() );
    auto& w = GetDistMat_VR_STAR( *wHandle );
    auto& X = GetCpxDistMat( *XHandle );
    
    HermitianGenDefiniteEig( AXBX, LOWER, A, B, w, X );
}

void ElemHermitianAxBxRange
( ElemCpxDistMat AHandle, ElemCpxDistMat BHandle,
  ElemDistMat_VR_STAR* wHandle, ElemCpxDistMat* XHandle,
  Real a, Real b )
{
    auto& A = GetCpxDistMat( AHandle );
    auto& B = GetCpxDistMat( BHandle );

    *wHandle = CreateDistMat_VR_STAR( A.Grid() ); 
    *XHandle = CreateCpxDistMat( A.Grid() );
    auto& w = GetDistMat_VR_STAR( *wHandle );
    auto& X = GetCpxDistMat( *XHandle );
    
    HermitianGenDefiniteEig( AXBX, LOWER, A, B, w, X, a, b );
}

void ElemHermitianAxBxIndices
( ElemCpxDistMat AHandle, ElemCpxDistMat BHandle,
  ElemDistMat_VR_STAR* wHandle, ElemCpxDistMat* XHandle,
  int a, int b )
{
    auto& A = GetCpxDistMat( AHandle );
    auto& B = GetCpxDistMat( BHandle );

    *wHandle = CreateDistMat_VR_STAR( A.Grid() ); 
    *XHandle = CreateCpxDistMat( A.Grid() );
    auto& w = GetDistMat_VR_STAR( *wHandle );
    auto& X = GetCpxDistMat( *XHandle );
    
    HermitianGenDefiniteEig( AXBX, LOWER, A, B, w, X, a, b );
}

//
// Utilities
//

int ElemLength( int n, int shift, int modulus )
{ return Length(n,shift,modulus); }

} // extern "C"
