/*
   Copyright (c) 2011-2013, Jack Poulson
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
int GetOpenIndex( std::vector<T*>& list )
{
    int index;
    for( index=0; index<list.size(); ++index )
        if( list[index] == 0 )
            break;
    if( index == list.size() )
        list.push_back( 0 );
    return index;
}

RealDistMatHandle CreateEmptyRealDistMat( const Grid& grid )
{
    const int index = GetOpenIndex( realDistMatList );
    realDistMatList[index] = new DistMatrix<double,MC,MR>( grid );
    return index;
}

ComplexDistMatHandle CreateEmptyComplexDistMat( const Grid& grid )
{
    const int index = GetOpenIndex( complexDistMatList );
    complexDistMatList[index] = new DistMatrix<Complex<double>,MC,MR>( grid );
    return index;
}

RealDistColVecHandle CreateEmptyRealDistColVec( const Grid& grid )
{
    const int index = GetOpenIndex( realDistColVecList );
    realDistColVecList[index] = new DistMatrix<double,VR,STAR>( grid );
    return index;
}

} // anonymous namespace

extern "C" {

//
// Environment controls
//

void FC_GLOBAL(initialize,NAME)()
{
    int argc = 0;
    char** argv = 0;
    elem::Initialize( argc, argv );
}

void FC_GLOBAL(finalize,NAME)()
{ elem::Finalize(); }

void FC_GLOBAL_(set_blocksize,NAME)( int blocksize )
{ elem::SetBlocksize( blocksize ); }

void FC_GLOBAL(blocksize,NAME)( int* blocksize )
{ *blocksize = elem::Blocksize(); }

void FC_GLOBAL_(set_normal_tridiag_approach,NAME)()
{ elem::SetHermitianTridiagApproach( HERMITIAN_TRIDIAG_NORMAL ); }

void FC_GLOBAL_(set_square_tridiag_approach,NAME)()
{ elem::SetHermitianTridiagApproach( HERMITIAN_TRIDIAG_SQUARE ); }

void FC_GLOBAL_(set_default_tridiag_approach,NAME)()
{ elem::SetHermitianTridiagApproach( HERMITIAN_TRIDIAG_DEFAULT ); }

void FC_GLOBAL_(set_row_major_tridiag_subgrid,NAME)()
{ elem::SetHermitianTridiagGridOrder( ROW_MAJOR ); }

void FC_GLOBAL_(set_col_major_tridiag_subgrid,NAME)()
{ elem::SetHermitianTridiagGridOrder( COLUMN_MAJOR ); }

//
// Process grid management
//

void FC_GLOBAL_(create_grid,NAME)( MPI_Fint* fComm, int* gridHandle )
{
    MPI_Comm comm = MPI_Comm_f2c( *fComm );
    const int index = GetOpenIndex( gridList );
    gridList[index] = new Grid( comm );
    *gridHandle = index;
}

void FC_GLOBAL_(grid_height,NAME)( int* handle, int* height )
{ *height = gridList[*handle]->Height(); }

void FC_GLOBAL_(grid_width,NAME)( int* handle, int* width )
{ *width = gridList[*handle]->Width(); }

void FC_GLOBAL_(grid_size,NAME)( int* handle, int* size )
{ *size = gridList[*handle]->Size(); }

void FC_GLOBAL_(grid_row,NAME)( int* handle, int* row )
{ *row = gridList[*handle]->MCRank(); }

void FC_GLOBAL_(grid_col,NAME)( int* handle, int* col )
{ *col = gridList[*handle]->MRRank(); }

void FC_GLOBAL_(grid_rank,NAME)( int* handle, int* rank )
{ *rank = gridList[*handle]->Rank(); }

void FC_GLOBAL_(free_grid,NAME)( int* handle )
{
    if( gridList[*handle] != 0 )
    {
        delete gridList[*handle];
        gridList[*handle] = 0;
    }
}

//
// Distributed matrix management
//

void FC_GLOBAL_(register_real_dist_mat,NAME)
( int* height, int* width, int* colAlignment, int* rowAlignment, 
  double* buffer, int* ldim, int* gridHandle, int* matHandle )
{
    const Grid& grid = TranslateGridHandle( *gridHandle );
    const int index = GetOpenIndex( realDistMatList );
    realDistMatList[index] = 
        new DistMatrix<double,MC,MR>
        (*height,*width,*colAlignment,*rowAlignment,buffer,*ldim,grid);
    *matHandle = index;
}

void FC_GLOBAL_(create_empty_real_dist_mat,NAME)
( int* gridHandle, int* matHandle )
{
    const Grid& grid = TranslateGridHandle( *gridHandle );
    const int index = GetOpenIndex( realDistMatList );
    realDistMatList[index] = new DistMatrix<double,MC,MR>( grid );
    *matHandle = index;
}

void FC_GLOBAL_(register_complex_dist_mat,NAME)
( int* height, int* width, int* colAlignment, int* rowAlignment, 
  void* voidBuffer, int* ldim, int* gridHandle, int* matHandle )
{
    typedef Complex<double> C;
    C* buffer = static_cast<C*>(voidBuffer);

    const Grid& grid = TranslateGridHandle( *gridHandle );
    const int index = GetOpenIndex( complexDistMatList );
    complexDistMatList[index] = 
        new DistMatrix<C,MC,MR>
        (*height,*width,*colAlignment,*rowAlignment,buffer,*ldim,grid);
    *matHandle = index;
}

void FC_GLOBAL_(create_empty_complex_dist_mat,NAME)
( int* gridHandle, int* matHandle )
{
    const Grid& grid = TranslateGridHandle( *gridHandle );
    const int index = GetOpenIndex( complexDistMatList );
    complexDistMatList[index] = new DistMatrix<Complex<double>,MC,MR>( grid );
    *matHandle = index;
}

void FC_GLOBAL_(free_real_dist_mat,NAME)( int* handle )
{
    if( realDistMatList[*handle] != 0 )
    {
        delete realDistMatList[*handle];
        realDistMatList[*handle] = 0;
    }
}

void FC_GLOBAL_(free_complex_dist_mat,NAME)( int* handle )
{
    if( complexDistMatList[*handle] != 0 )
    {
        delete complexDistMatList[*handle];
        complexDistMatList[*handle] = 0;
    }
}

void FC_GLOBAL_(print_real_dist_mat,NAME)( int* AHandle )
{
    const DistMatrix<double,MC,MR>& A = TranslateRealDistMatHandle( *AHandle );
    A.Print();
}

void FC_GLOBAL_(print_complex_dist_mat,NAME)( int* AHandle )
{
    typedef Complex<double> C;
    const DistMatrix<C,MC,MR>& A = TranslateComplexDistMatHandle( *AHandle );
    A.Print();
}

//
// Distributed column vector management
//

void FC_GLOBAL_(create_empty_real_dist_col_vec,NAME)
( int* gridHandle, int* vecHandle )
{
    const Grid& grid = TranslateGridHandle( *gridHandle );
    const int index = GetOpenIndex( realDistColVecList );
    realDistColVecList[index] = new DistMatrix<double,VR,STAR>( grid );
    *vecHandle = index;
}

void FC_GLOBAL_(free_real_dist_col_vec,NAME)( int* handle )
{
    if( realDistColVecList[*handle] != 0 )
    {
        delete realDistColVecList[*handle];
        realDistColVecList[*handle] = 0;
    }
}

void FC_GLOBAL_(print_real_dist_col_vec,NAME)( int* AHandle )
{
    const DistMatrix<double,VR,STAR>& A = 
        TranslateRealDistColVecHandle( *AHandle );
    A.Print();
}

//
// Generalized Hermitian-definite eigensolvers for A X = B X \Lambda
//

void FC_GLOBAL_(symmetric_axbx,NAME)
( int* AHandle, int* BHandle,
  int* wHandle, int* XHandle )
{
    DistMatrix<double,MC,MR>& A   = TranslateRealDistMatHandle( *AHandle );
    DistMatrix<double,MC,MR>& B   = TranslateRealDistMatHandle( *BHandle );

    *wHandle = CreateEmptyRealDistColVec( A.Grid() );
    *XHandle = CreateEmptyRealDistMat( A.Grid() );
    DistMatrix<double,VR,STAR>& w = TranslateRealDistColVecHandle( *wHandle );
    DistMatrix<double,MC,MR>& X   = TranslateRealDistMatHandle( *XHandle );
    
    HermitianGenDefiniteEig( AXBX, LOWER, A, B, w, X );
}

void FC_GLOBAL_(symmetric_axbx_partial_range,NAME)
( int* AHandle, int* BHandle,
  int* wHandle, int* XHandle,
  double* a, double* b )
{
    DistMatrix<double,MC,MR>& A   = TranslateRealDistMatHandle( *AHandle );
    DistMatrix<double,MC,MR>& B   = TranslateRealDistMatHandle( *BHandle );

    *wHandle = CreateEmptyRealDistColVec( A.Grid() );
    *XHandle = CreateEmptyRealDistMat( A.Grid() );
    DistMatrix<double,VR,STAR>& w = TranslateRealDistColVecHandle( *wHandle );
    DistMatrix<double,MC,MR>& X   = TranslateRealDistMatHandle( *XHandle );
    
    HermitianGenDefiniteEig( AXBX, LOWER, A, B, w, X, *a, *b );
}

void FC_GLOBAL_(symmetric_axbx_partial_indices,NAME)
( int* AHandle, int* BHandle,
  int* wHandle, int* XHandle,
  int* a, int* b ) 
{
    DistMatrix<double,MC,MR>& A   = TranslateRealDistMatHandle( *AHandle );
    DistMatrix<double,MC,MR>& B   = TranslateRealDistMatHandle( *BHandle );

    *wHandle = CreateEmptyRealDistColVec( A.Grid() );
    *XHandle = CreateEmptyRealDistMat( A.Grid() );
    DistMatrix<double,VR,STAR>& w = TranslateRealDistColVecHandle( *wHandle );
    DistMatrix<double,MC,MR>& X   = TranslateRealDistMatHandle( *XHandle );

    // Convert from Fortran to C indexing
    int aC = *a-1;
    int bC = *b-1;
    
    HermitianGenDefiniteEig( AXBX, LOWER, A, B, w, X, aC, bC );
}

void FC_GLOBAL_(hermitian_axbx,NAME)
( int* AHandle, int* BHandle,
  int* wHandle, int* XHandle )
{
    typedef Complex<double> C;

    DistMatrix<C,MC,MR>& A        = TranslateComplexDistMatHandle( *AHandle );
    DistMatrix<C,MC,MR>& B        = TranslateComplexDistMatHandle( *BHandle );

    *wHandle = CreateEmptyRealDistColVec( A.Grid() );
    *XHandle = CreateEmptyComplexDistMat( A.Grid() );
    DistMatrix<double,VR,STAR>& w = TranslateRealDistColVecHandle( *wHandle );
    DistMatrix<C,MC,MR>& X        = TranslateComplexDistMatHandle( *XHandle );
    
    HermitianGenDefiniteEig( AXBX, LOWER, A, B, w, X );
}

void FC_GLOBAL_(hermitian_axbx_partial_range,NAME)
( int* AHandle, int* BHandle,
  int* wHandle, int* XHandle,
  double* a, double* b )
{
    typedef Complex<double> C;

    DistMatrix<C,MC,MR>& A        = TranslateComplexDistMatHandle( *AHandle );
    DistMatrix<C,MC,MR>& B        = TranslateComplexDistMatHandle( *BHandle );
    
    *wHandle = CreateEmptyRealDistColVec( A.Grid() );
    *XHandle = CreateEmptyComplexDistMat( A.Grid() );
    DistMatrix<double,VR,STAR>& w = TranslateRealDistColVecHandle( *wHandle );
    DistMatrix<C,MC,MR>& X        = TranslateComplexDistMatHandle( *XHandle );
    
    HermitianGenDefiniteEig( AXBX, LOWER, A, B, w, X, *a, *b );
}

void FC_GLOBAL_(hermitian_axbx_partial_indices,NAME)
( int* AHandle, int* BHandle,
  int* wHandle, int* XHandle,
  int* a, int* b )
{
    typedef Complex<double> C;

    DistMatrix<C,MC,MR>& A        = TranslateComplexDistMatHandle( *AHandle );
    DistMatrix<C,MC,MR>& B        = TranslateComplexDistMatHandle( *BHandle );
    
    *wHandle = CreateEmptyRealDistColVec( A.Grid() );
    *XHandle = CreateEmptyComplexDistMat( A.Grid() );
    DistMatrix<double,VR,STAR>& w = TranslateRealDistColVecHandle( *wHandle );
    DistMatrix<C,MC,MR>& X        = TranslateComplexDistMatHandle( *XHandle );
    
    // Convert from Fortran to C indexing
    int aC = *a-1;
    int bC = *b-1;

    HermitianGenDefiniteEig( AXBX, LOWER, A, B, w, X, aC, bC );
}

//
// Utilities
//

void FC_GLOBAL_(local_length,NAME)
( int* n, int* shift, int* modulus, int* localLength )
{ *localLength = elem::LocalLength<int>(*n,*shift,*modulus); }

} // extern "C"
