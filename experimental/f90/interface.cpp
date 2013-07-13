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

typedef int ElemGrid;

typedef int ElemDistMat;
typedef int ElemCpxDistMat;

typedef int ElemDistMat_VC_STAR;
typedef int ElemCpxDistMat_VC_STAR;

typedef int ElemDistMat_VR_STAR;
typedef int ElemCpxDistMat_VR_STAR;

typedef int ElemDistMat_STAR_STAR;
typedef int ElemCpxDistMat_STAR_STAR;

namespace {

const int maxInt = std::numeric_limits<int>::max();

std::vector<Grid*> gridList;

std::vector<DistMatrix<double>*> distMatList;
std::vector<DistMatrix<Complex<double> >*> cpxDistMatList;

std::vector<DistMatrix<double,VC,  STAR>*> distMatList_VC_STAR;
std::vector<DistMatrix<double,VR,  STAR>*> distMatList_VR_STAR;
std::vector<DistMatrix<double,STAR,STAR>*> distMatList_STAR_STAR;

const Grid& GetGrid( ElemGrid grid )
{
    if( grid == maxInt )
        return DefaultGrid();
    else
        return *gridList[grid];
}

DistMatrix<double>& GetDistMat( ElemDistMat A )
{ return *distMatList[A]; }

DistMatrix<Complex<double> >& GetCpxDistMat( ElemCpxDistMat A )
{ return *cpxDistMatList[A]; }

DistMatrix<double,VC,STAR>& GetDistMat_VC_STAR( ElemDistMat_VC_STAR A )
{ return *distMatList_VC_STAR[A]; }

DistMatrix<double,VR,STAR>& GetDistMat_VR_STAR( ElemDistMat_VR_STAR A )
{ return *distMatList_VR_STAR[A]; }

DistMatrix<double,STAR,STAR>& GetDistMat_STAR_STAR( ElemDistMat_STAR_STAR A )
{ return *distMatList_STAR_STAR[A]; }

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

ElemDistMat CreateDistMat( const Grid& grid )
{
    const int index = GetOpenIndex( distMatList );
    distMatList[index] = new DistMatrix<double>( grid );
    return index;
}

ElemCpxDistMat CreateCpxDistMat( const Grid& grid )
{
    const int index = GetOpenIndex( cpxDistMatList );
    cpxDistMatList[index] = new DistMatrix<Complex<double> >( grid );
    return index;
}

ElemDistMat_VC_STAR CreateDistMat_VC_STAR( const Grid& grid )
{
    const int index = GetOpenIndex( distMatList_VC_STAR );
    distMatList_VC_STAR[index] = new DistMatrix<double,VC,STAR>( grid );
    return index;
}

ElemDistMat_VR_STAR CreateDistMat_VR_STAR( const Grid& grid )
{
    const int index = GetOpenIndex( distMatList_VR_STAR );
    distMatList_VR_STAR[index] = new DistMatrix<double,VR,STAR>( grid );
    return index;
}

ElemDistMat_STAR_STAR CreateDistMat_STAR_STAR( const Grid& grid )
{
    const int index = GetOpenIndex( distMatList_STAR_STAR );
    distMatList_STAR_STAR[index] = new DistMatrix<double,STAR,STAR>( grid );
    return index;
}

} // anonymous namespace

extern "C" {

//
// Environment controls
//

void FC_GLOBAL(elem_initialize,NAME)()
{
    int argc = 0;
    char** argv = 0;
    Initialize( argc, argv );
}

void FC_GLOBAL(elem_finalize,NAME)()
{ Finalize(); }

void FC_GLOBAL_(elem_set_blocksize,NAME)( int blocksize )
{ SetBlocksize( blocksize ); }

void FC_GLOBAL(elem_blocksize,NAME)( int* blocksize )
{ *blocksize = Blocksize(); }

void FC_GLOBAL_(elem_set_normal_tridiag_approach,NAME)()
{ SetHermitianTridiagApproach( HERMITIAN_TRIDIAG_NORMAL ); }

void FC_GLOBAL_(elem_set_square_tridiag_approach,NAME)()
{ SetHermitianTridiagApproach( HERMITIAN_TRIDIAG_SQUARE ); }

void FC_GLOBAL_(elem_set_default_tridiag_approach,NAME)()
{ SetHermitianTridiagApproach( HERMITIAN_TRIDIAG_DEFAULT ); }

void FC_GLOBAL_(elem_set_row_major_tridiag_subgrid,NAME)()
{ SetHermitianTridiagGridOrder( ROW_MAJOR ); }

void FC_GLOBAL_(elem_set_col_major_tridiag_subgrid,NAME)()
{ SetHermitianTridiagGridOrder( COLUMN_MAJOR ); }

//
// Process grid management
//

void FC_GLOBAL_(elem_default_grid,NAME)( int* gridHandle )
{ *gridHandle = maxInt; }

void FC_GLOBAL_(elem_create_grid,NAME)( int* gridHandle, MPI_Fint* fComm )
{
    MPI_Comm comm = MPI_Comm_f2c( *fComm );
    const int index = GetOpenIndex( gridList );
    gridList[index] = new Grid( comm );
    *gridHandle = index;
}

void FC_GLOBAL_(elem_create_specific_grid,NAME)
( int* gridHandle, MPI_Fint* fComm, int* height )
{
    MPI_Comm comm = MPI_Comm_f2c( *fComm );
    const int index = GetOpenIndex( gridList );
    gridList[index] = new Grid( comm, *height );
    *gridHandle = index;
}

void FC_GLOBAL_(elem_grid_height,NAME)( int* handle, int* height )
{ *height = gridList[*handle]->Height(); }

void FC_GLOBAL_(elem_grid_width,NAME)( int* handle, int* width )
{ *width = gridList[*handle]->Width(); }

void FC_GLOBAL_(elem_grid_size,NAME)( int* handle, int* size )
{ *size = gridList[*handle]->Size(); }

void FC_GLOBAL_(elem_grid_row,NAME)( int* handle, int* row )
{ *row = gridList[*handle]->MCRank(); }

void FC_GLOBAL_(elem_grid_col,NAME)( int* handle, int* col )
{ *col = gridList[*handle]->MRRank(); }

void FC_GLOBAL_(elem_grid_rank,NAME)( int* handle, int* rank )
{ *rank = gridList[*handle]->Rank(); }

void FC_GLOBAL_(elem_free_grid,NAME)( int* handle )
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

void FC_GLOBAL_(elem_register_dist_mat,NAME)
( int* handle, 
  int* height, int* width, int* colAlignment, int* rowAlignment, 
  double* buffer, int* ldim, int* gridHandle )
{
    const Grid& grid = GetGrid( *gridHandle );
    const int index = GetOpenIndex( distMatList );
    try 
    {
        distMatList[index] = 
            new DistMatrix<double>
            (*height,*width,*colAlignment,*rowAlignment,buffer,*ldim,grid);
    } catch( std::exception& e ) { ReportException(e); }
    *handle = index;
}

void FC_GLOBAL_(elem_register_cpx_dist_mat,NAME)
( int* handle, 
  int* height, int* width, int* colAlignment, int* rowAlignment, 
  void* voidBuffer, int* ldim, int* gridHandle )
{
    typedef Complex<double> C;
    C* buffer = static_cast<C*>(voidBuffer);

    const Grid& grid = GetGrid( *gridHandle );
    const int index = GetOpenIndex( cpxDistMatList );
    try
    {
        cpxDistMatList[index] = 
            new DistMatrix<C>
            (*height,*width,*colAlignment,*rowAlignment,buffer,*ldim,grid);
    } catch( std::exception& e ) { ReportException(e); }
    *handle = index;
}

void FC_GLOBAL_(elem_create_dist_mat,NAME)
( int* handle, int* gridHandle )
{
    const Grid& grid = GetGrid( *gridHandle );
    const int index = GetOpenIndex( distMatList );
    try
    {
        distMatList[index] = new DistMatrix<double>( grid );
    } catch( std::exception& e ) { ReportException(e); }
    *handle = index;
}

void FC_GLOBAL_(elem_create_cpx_dist_mat,NAME)
( int* handle, int* gridHandle )
{
    const Grid& grid = GetGrid( *gridHandle );
    const int index = GetOpenIndex( cpxDistMatList );
    try
    {
        cpxDistMatList[index] = new DistMatrix<Complex<double> >( grid );
    } catch( std::exception& e ) { ReportException(e); }
    *handle = index;
}

void FC_GLOBAL_(elem_free_dist_mat,NAME)( int* handle )
{
    if( distMatList[*handle] != 0 )
    {
        delete distMatList[*handle];
        distMatList[*handle] = 0;
    }
}

void FC_GLOBAL_(elem_free_cpx_dist_mat,NAME)( int* handle )
{
    if( cpxDistMatList[*handle] != 0 )
    {
        delete cpxDistMatList[*handle];
        cpxDistMatList[*handle] = 0;
    }
}

void FC_GLOBAL_(elem_print_dist_mat,NAME)( int* AHandle )
{
    const DistMatrix<double>& A = GetDistMat( *AHandle );
    try
    {
        Print( A );
    } catch( std::exception& e ) { ReportException(e); }
}

void FC_GLOBAL_(elem_print_cpx_dist_mat,NAME)( int* AHandle )
{
    typedef Complex<double> C;
    const DistMatrix<C>& A = GetCpxDistMat( *AHandle );
    try
    {
        Print( A );
    } catch( std::exception& e ) { ReportException(e); }
}

//
// [VC,* ] management
//

void FC_GLOBAL_(elem_register_dist_mat_vc_star,NAME)
( int* handle, 
  int* height, int* width, int* colAlignment, 
  double* buffer, int* ldim, int* gridHandle )
{
    const Grid& grid = GetGrid( *gridHandle );
    const int index = GetOpenIndex( distMatList_VC_STAR );
    try
    {
        distMatList_VC_STAR[index] = 
            new DistMatrix<double,VC,STAR>
            (*height,*width,*colAlignment,buffer,*ldim,grid);
    } catch( std::exception& e ) { ReportException(e); }
    *handle = index;
}

void FC_GLOBAL_(elem_create_dist_mat_vc_star,NAME)
( int* handle, int* gridHandle )
{
    const Grid& grid = GetGrid( *gridHandle );
    const int index = GetOpenIndex( distMatList_VC_STAR );
    try
    {
        distMatList_VC_STAR[index] = new DistMatrix<double,VC,STAR>( grid );
    } catch( std::exception& e ) { ReportException(e); }
    *handle = index;
}

void FC_GLOBAL_(elem_free_dist_mat_vc_star,NAME)( int* handle )
{
    if( distMatList_VC_STAR[*handle] != 0 )
    {
        delete distMatList_VC_STAR[*handle];
        distMatList_VC_STAR[*handle] = 0;
    }
}

void FC_GLOBAL_(elem_print_dist_mat_vc_star,NAME)( int* AHandle )
{
    const DistMatrix<double,VC,STAR>& A = GetDistMat_VC_STAR( *AHandle );
    try
    {
        Print( A );
    } catch( std::exception& e ) { ReportException(e); }
}

//
// [VR,* ] management
//

void FC_GLOBAL_(elem_register_dist_mat_vr_star,NAME)
( int* handle, 
  int* height, int* width, int* colAlignment, 
  double* buffer, int* ldim, int* gridHandle )
{
    const Grid& grid = GetGrid( *gridHandle );
    const int index = GetOpenIndex( distMatList_VR_STAR );
    try
    {
        distMatList_VR_STAR[index] = 
            new DistMatrix<double,VR,STAR>
            (*height,*width,*colAlignment,buffer,*ldim,grid);
    } catch( std::exception& e ) { ReportException(e); }
    *handle = index;
}

void FC_GLOBAL_(elem_create_dist_mat_vr_star,NAME)
( int* handle, int* gridHandle )
{
    const Grid& grid = GetGrid( *gridHandle );
    const int index = GetOpenIndex( distMatList_VR_STAR );
    try
    {
        distMatList_VR_STAR[index] = new DistMatrix<double,VR,STAR>( grid );
    } catch( std::exception& e ) { ReportException(e); }
    *handle = index;
}

void FC_GLOBAL_(elem_free_dist_mat_vr_star,NAME)( int* handle )
{
    if( distMatList_VR_STAR[*handle] != 0 )
    {
        delete distMatList_VR_STAR[*handle];
        distMatList_VR_STAR[*handle] = 0;
    }
}

void FC_GLOBAL_(elem_print_dist_mat_vr_star,NAME)( int* AHandle )
{
    const DistMatrix<double,VR,STAR>& A = GetDistMat_VR_STAR( *AHandle );
    try
    {
        Print( A );
    } catch( std::exception& e ) { ReportException(e); }
}

//
// [* ,* ] management
//

void FC_GLOBAL_(elem_register_dist_mat_star_star,NAME)
( int* handle, 
  int* height, int* width, double* buffer, int* ldim, int* gridHandle )
{
    const Grid& grid = GetGrid( *gridHandle );
    const int index = GetOpenIndex( distMatList_STAR_STAR );
    try
    {
        distMatList_STAR_STAR[index] = 
            new DistMatrix<double,STAR,STAR>(*height,*width,buffer,*ldim,grid);
    } catch( std::exception& e ) { ReportException(e); }
    *handle = index;
}

void FC_GLOBAL_(elem_create_dist_mat_star_star,NAME)
( int* handle, int* gridHandle )
{
    const Grid& grid = GetGrid( *gridHandle );
    const int index = GetOpenIndex( distMatList_STAR_STAR );
    try
    {
        distMatList_STAR_STAR[index] = new DistMatrix<double,STAR,STAR>( grid );
    } catch( std::exception& e ) { ReportException(e); }
    *handle = index;
}

void FC_GLOBAL_(elem_free_dist_mat_star_star,NAME)( int* handle )
{
    if( distMatList_STAR_STAR[*handle] != 0 )
    {
        delete distMatList_STAR_STAR[*handle];
        distMatList_STAR_STAR[*handle] = 0;
    }
}

void FC_GLOBAL_(elem_print_dist_mat_star_star,NAME)( int* AHandle )
{
    const DistMatrix<double,STAR,STAR>& A = GetDistMat_STAR_STAR( *AHandle );
    try
    {
        Print( A );
    } catch( std::exception& e ) { ReportException(e); }
}

//
// Generalized Hermitian-definite eigensolvers for A X = B X \Lambda
//

void FC_GLOBAL_(elem_symmetric_eig,NAME)
( int* AHandle, int* wHandle, int* XHandle )
{
    DistMatrix<double>& A = GetDistMat( *AHandle );
    DistMatrix<double>& X = GetDistMat( *XHandle );
    DistMatrix<double,STAR,STAR>& w = GetDistMat_STAR_STAR( *wHandle );
    
    try
    {
        DistMatrix<double,VR,STAR> w_VR_STAR( A.Grid() );
        HermitianEig( LOWER, A, w_VR_STAR, X );
        w = w_VR_STAR;
    } catch( std::exception& e ) { ReportException(e); }
}

void FC_GLOBAL_(elem_hermitian_eig,NAME)
( int* AHandle, int* wHandle, int* XHandle )
{
    typedef Complex<double> C;
    DistMatrix<C>& A = GetCpxDistMat( *AHandle );
    DistMatrix<C>& X = GetCpxDistMat( *XHandle );
    DistMatrix<double,STAR,STAR>& w = GetDistMat_STAR_STAR( *wHandle );
    
    try
    {
        DistMatrix<double,VR,STAR> w_VR_STAR( A.Grid() );
        HermitianEig( LOWER, A, w_VR_STAR, X );
        w = w_VR_STAR;
    } catch( std::exception& e ) { ReportException(e); }
}

void FC_GLOBAL_(elem_symmetric_axbx_reduce,NAME)
( int* AHandle, int* BHandle )
{
    DistMatrix<double>& A = GetDistMat( *AHandle );
    DistMatrix<double>& B = GetDistMat( *BHandle );

    try
    {
        Cholesky( LOWER, B );
        TwoSidedTrsm( LOWER, NON_UNIT, A, B );
    } catch( std::exception& e ) { ReportException(e); }
}

void FC_GLOBAL_(elem_hermitian_axbx_reduce,NAME)
( int* AHandle, int* BHandle )
{
    DistMatrix<Complex<double> >& A = GetCpxDistMat( *AHandle );
    DistMatrix<Complex<double> >& B = GetCpxDistMat( *BHandle );

    try
    {
        Cholesky( LOWER, B );
        TwoSidedTrsm( LOWER, NON_UNIT, A, B );
    } catch( std::exception& e ) { ReportException(e); }
}

void FC_GLOBAL_(elem_symmetric_axbx_expand,NAME)
( int* AHandle, int* BHandle, int* XHandle )
{
    DistMatrix<double>& A = GetDistMat( *AHandle );
    DistMatrix<double>& B = GetDistMat( *BHandle );
    DistMatrix<double>& X = GetDistMat( *XHandle );
    
    try
    {
        Trsm( LEFT, LOWER, ADJOINT, NON_UNIT, 1., B, X );
    } catch( std::exception& e ) { ReportException(e); }
}

void FC_GLOBAL_(elem_hermitian_axbx_expand,NAME)
( int* AHandle, int* BHandle, int* XHandle )
{
    typedef Complex<double> C;
    DistMatrix<C>& A = GetCpxDistMat( *AHandle );
    DistMatrix<C>& B = GetCpxDistMat( *BHandle );
    DistMatrix<C>& X = GetCpxDistMat( *XHandle );
    
    try
    {
        Trsm( LEFT, LOWER, ADJOINT, NON_UNIT, C(1.), B, X );
    } catch( std::exception& e ) { ReportException(e); }
}

void FC_GLOBAL_(elem_symmetric_axbx,NAME)
( int* AHandle, int* BHandle,
  int* wHandle, int* XHandle )
{
    DistMatrix<double>& A = GetDistMat( *AHandle );
    DistMatrix<double>& B = GetDistMat( *BHandle );
    DistMatrix<double>& X = GetDistMat( *XHandle );
    DistMatrix<double,STAR,STAR>& w = GetDistMat_STAR_STAR( *wHandle );
    
    try
    {
        DistMatrix<double,VR,STAR> w_VR_STAR( A.Grid() );
        HermitianGenDefiniteEig( AXBX, LOWER, A, B, w_VR_STAR, X );
        w = w_VR_STAR;
    } catch( std::exception& e ) { ReportException(e); }
}

void FC_GLOBAL_(elem_symmetric_axbx_range,NAME)
( int* AHandle, int* BHandle,
  int* wHandle, int* XHandle,
  double* a, double* b )
{
    DistMatrix<double>& A = GetDistMat( *AHandle );
    DistMatrix<double>& B = GetDistMat( *BHandle );
    DistMatrix<double>& X = GetDistMat( *XHandle );
    DistMatrix<double,STAR,STAR>& w = GetDistMat_STAR_STAR( *wHandle );
    
    try
    {
        DistMatrix<double,VR,STAR> w_VR_STAR( A.Grid() );
        HermitianGenDefiniteEig( AXBX, LOWER, A, B, w_VR_STAR, X, *a, *b );
        w = w_VR_STAR; 
    } catch( std::exception& e ) { ReportException(e); }
}

void FC_GLOBAL_(elem_symmetric_axbx_indices,NAME)
( int* AHandle, int* BHandle,
  int* wHandle, int* XHandle,
  int* a, int* b ) 
{
    DistMatrix<double>& A = GetDistMat( *AHandle );
    DistMatrix<double>& B = GetDistMat( *BHandle );
    DistMatrix<double>& X = GetDistMat( *XHandle );
    DistMatrix<double,STAR,STAR>& w = GetDistMat_STAR_STAR( *wHandle );

    // Convert from Fortran to C indexing
    int aC = *a-1;
    int bC = *b-1;
    
    try
    {
        DistMatrix<double,VR,STAR> w_VR_STAR( A.Grid() );
        HermitianGenDefiniteEig( AXBX, LOWER, A, B, w_VR_STAR, X, aC, bC );
        w = w_VR_STAR;
    } catch( std::exception& e ) { ReportException(e); }
}

void FC_GLOBAL_(elem_hermitian_axbx,NAME)
( int* AHandle, int* BHandle,
  int* wHandle, int* XHandle )
{
    typedef Complex<double> C;
    DistMatrix<C>& A = GetCpxDistMat( *AHandle );
    DistMatrix<C>& B = GetCpxDistMat( *BHandle );
    DistMatrix<C>& X = GetCpxDistMat( *XHandle );
    DistMatrix<double,STAR,STAR>& w = GetDistMat_STAR_STAR( *wHandle );
    
    try
    {
        DistMatrix<double,VR,STAR> w_VR_STAR( A.Grid() );
        HermitianGenDefiniteEig( AXBX, LOWER, A, B, w_VR_STAR, X );
        w = w_VR_STAR;
    } catch( std::exception& e ) { ReportException(e); }
}

void FC_GLOBAL_(elem_hermitian_axbx_range,NAME)
( int* AHandle, int* BHandle,
  int* wHandle, int* XHandle,
  double* a, double* b )
{
    typedef Complex<double> C;
    DistMatrix<C>& A = GetCpxDistMat( *AHandle );
    DistMatrix<C>& B = GetCpxDistMat( *BHandle );
    DistMatrix<C>& X = GetCpxDistMat( *XHandle );
    DistMatrix<double,STAR,STAR>& w = GetDistMat_STAR_STAR( *wHandle );
    
    try
    {
        DistMatrix<double,VR,STAR> w_VR_STAR( A.Grid() );
        HermitianGenDefiniteEig( AXBX, LOWER, A, B, w_VR_STAR, X, *a, *b );
        w = w_VR_STAR;
    } catch( std::exception& e ) { ReportException(e); }
}

void FC_GLOBAL_(elem_hermitian_axbx_indices,NAME)
( int* AHandle, int* BHandle,
  int* wHandle, int* XHandle,
  int* a, int* b )
{
    typedef Complex<double> C;
    DistMatrix<C>& A = GetCpxDistMat( *AHandle );
    DistMatrix<C>& B = GetCpxDistMat( *BHandle );
    DistMatrix<C>& X = GetCpxDistMat( *XHandle );
    DistMatrix<double,STAR,STAR>& w = GetDistMat_STAR_STAR( *wHandle );
    
    // Convert from Fortran to C indexing
    int aC = *a-1;
    int bC = *b-1;

    try
    {
        DistMatrix<double,VR,STAR> w_VR_STAR( A.Grid() );
        HermitianGenDefiniteEig( AXBX, LOWER, A, B, w_VR_STAR, X, aC, bC );
        w = w_VR_STAR;
    } catch( std::exception& e ) { ReportException(e); }
}

//
// Utilities
//

void FC_GLOBAL_(elem_length,NAME)
( int* n, int* shift, int* modulus, int* length )
{ *length = Length<int>(*n,*shift,*modulus); }

// Returns the global and local sizes of zero-aligned n x k matrix
void FC_GLOBAL_(elem_padded_eigvec_size,NAME)
( int* n, int* k, int* gridHandle, int* nPad, int* kPad, int* nLoc, int* kLoc )
{
    const Grid& grid = GetGrid( *gridHandle );    
    const int r = grid.Height();
    const int c = grid.Width();
    const int p = grid.Size();
    const int row = grid.Row();
    const int col = grid.Col();
    *nPad = MaxLength(*n,r)*r;
    *kPad = MaxLength(*n,p)*p;
    *nLoc = Length(*nPad,row,r);
    *kLoc = Length(*kPad,col,c);
}

} // extern "C"
