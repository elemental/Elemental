/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

#ifdef EL_HAVE_SCALAPACK

extern "C" {

int Csys2blacs_handle( MPI_Comm comm );
void Cblacs_gridinit
( int* context, const char* order, int gridHeight, int gridWidth );
void Cblacs_gridinfo
( int  context, int* gridHeight, int* gridWidth, int* gridRow, int* gridCol );
void Cfree_blacs_system_handle( int bhandle );
void Cblacs_gridexit( int context );
void Cblacs_exit( int notDone );

// TODO: Decide if the context pointers can be const
void EL_SCALAPACK(psgemr2d)
( const int* m, const int* n,
  const float* A, const int* iA, const int* jA, const int* descA,
        float* B, const int* iB, const int* jB, const int* descB,
  int* context );
void EL_SCALAPACK(pdgemr2d)
( const int* m, const int* n,
  const double* A, const int* iA, const int* jA, const int* descA,
        double* B, const int* iB, const int* jB, const int* descB,
  int* context );
void EL_SCALAPACK(pcgemr2d)
( const int* m, const int* n,
  const El::scomplex* A, const int* iA, const int* jA, const int* descA,
        El::scomplex* B, const int* iB, const int* jB, const int* descB,
  int* context );
void EL_SCALAPACK(pzgemr2d)
( const int* m, const int* n,
  const El::dcomplex* A, const int* iA, const int* jA, const int* descA,
        El::dcomplex* B, const int* iB, const int* jB, const int* descB,
  int* context );

} // extern "C"

namespace El {
namespace blacs {

int Handle( MPI_Comm comm )
{ return Csys2blacs_handle( comm ); }

int GridInit( int bhandle, bool colMajor, int gridHeight, int gridWidth )
{ 
    int context = bhandle;
    const char* order = ( colMajor ? "Col" : "Row" );
    Cblacs_gridinit( &context, order, gridHeight, gridWidth ); 
    return context;
}

int GridHeight( int context )
{
    int gridHeight, gridWidth, gridRow, gridCol;
    Cblacs_gridinfo( context, &gridHeight, &gridWidth, &gridRow, &gridCol );
    return gridHeight;
}

int GridWidth( int context )
{
    int gridHeight, gridWidth, gridRow, gridCol;
    Cblacs_gridinfo( context, &gridHeight, &gridWidth, &gridRow, &gridCol );
    return gridWidth;
}

int GridRow( int context )
{
    int gridHeight, gridWidth, gridRow, gridCol;
    Cblacs_gridinfo( context, &gridHeight, &gridWidth, &gridRow, &gridCol );
    return gridRow;
}

int GridCol( int context )
{
    int gridHeight, gridWidth, gridRow, gridCol;
    Cblacs_gridinfo( context, &gridHeight, &gridWidth, &gridRow, &gridCol );
    return gridCol;
}

void FreeHandle( int bhandle )
{ Cfree_blacs_system_handle( bhandle ); }

void FreeGrid( int context )
{ Cblacs_gridexit( context ); }

void Exit( bool finished )
{ 
    int notDone = ( finished ? 0 : 1 );
    Cblacs_exit( notDone ); 
}

void Redistribute
( int m, int n,
  const float* A, const int* descA, 
        float* B, const int* descB, int context )
{
    int iA=1, jA=1, iB=1, jB=1;
    EL_SCALAPACK(psgemr2d)
    ( &m, &n,
      A, &iA, &jA, descA,
      B, &iB, &jB, descB, &context );
}

void Redistribute
( int m, int n,
  const double* A, const int* descA, 
        double* B, const int* descB, int context )
{
    int iA=1, jA=1, iB=1, jB=1;
    EL_SCALAPACK(pdgemr2d)
    ( &m, &n,
      A, &iA, &jA, descA,
      B, &iB, &jB, descB, &context );
}

void Redistribute
( int m, int n,
  const Complex<float>* A, const int* descA, 
        Complex<float>* B, const int* descB, int context )
{
    int iA=1, jA=1, iB=1, jB=1;
    EL_SCALAPACK(pcgemr2d)
    ( &m, &n,
      A, &iA, &jA, descA,
      B, &iB, &jB, descB, &context );
}

void Redistribute
( int m, int n,
  const Complex<double>* A, const int* descA, 
        Complex<double>* B, const int* descB, int context )
{
    int iA=1, jA=1, iB=1, jB=1;
    EL_SCALAPACK(pzgemr2d)
    ( &m, &n,
      A, &iA, &jA, descA,
      B, &iB, &jB, descB, &context );
}

} // namespace blacs
} // namespace El

#endif // ifdef EL_HAVE_SCALAPACK
