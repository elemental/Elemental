/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "elemental-lite.hpp"

#ifdef ELEM_HAVE_SCALAPACK

using elem::scomplex;
using elem::dcomplex;

extern "C" {

// Basic Linear Algebra Communication Subprograms
// ==============================================
int Csys2blacs_handle( MPI_Comm comm );
void Cblacs_gridinit
( int* context, const char* order, int gridHeight, int gridWidth );
void Cblacs_gridinfo
( int  context, int* gridHeight, int* gridWidth, int* gridRow, int* gridCol );
void Cfree_blacs_system_handle( int bhandle );
void Cblacs_gridexit( int context );
void Cblacs_exit( int notDone );

// Hessenberg QR algorithm
// =======================

// Aggressive Early Deflation
// --------------------------
// NOTE: ScaLAPACK currently only supports AED for real matrices
void ELEM_SCALAPACK(pshseqr)
( const char* job, const char* compz, 
  const int* n, const int* ilo, const int* ihi, 
  float* H, const int* desch, float* wr, float* wi, 
  float* Z, const int* descz, float* work, const int* lwork, 
  int* iwork, const int* liwork, int* info );
void ELEM_SCALAPACK(pdhseqr)
( const char* job, const char* compz, 
  const int* n, const int* ilo, const int* ihi, 
  double* H, const int* desch, double* wr, double* wi, 
  double* Z, const int* descz, double* work, const int* lwork, 
  int* iwork, const int* liwork, int* info );

// Pipelined QR algorithm without AED
// ----------------------------------
void ELEM_SCALAPACK(pclahqr)
( const ELEM_FORT_LOGICAL* wantt, const ELEM_FORT_LOGICAL* wantz, const int* n,
  const int* ilo, const int* ihi, scomplex* H, const int* desch,
  scomplex* w, const int* iloz, const int* ihiz, scomplex* Z, const int* descz,
  scomplex* work, const int* lwork, int* iwork, const int* liwork, int* info );
void ELEM_SCALAPACK(pzlahqr)
( const ELEM_FORT_LOGICAL* wantt, const ELEM_FORT_LOGICAL* wantz, const int* n,
  const int* ilo, const int* ihi, dcomplex* H, const int* desch,
  dcomplex* w, const int* iloz, const int* ihiz, dcomplex* Z, const int* descz,
  dcomplex* work, const int* lwork, int* iwork, const int* liwork, int* info );

// Pipelined QR algorithm with AED for big matrices
// ------------------------------------------------
void ELEM_SCALAPACK(pslaqr0)
( const ELEM_FORT_LOGICAL* wantt, const ELEM_FORT_LOGICAL* wantz, const int* n,
  const int* ilo, const int* ihi, float* H, const int* desch, 
  float* wr, float* wi, const int* iloz, const int* ihiz, 
  float* Z, const int* descz, float* work, const int* lwork, 
  int* iwork, const int* liwork, int* info, const int* reclevel );
void ELEM_SCALAPACK(pdlaqr0)
( const ELEM_FORT_LOGICAL* wantt, const ELEM_FORT_LOGICAL* wantz, const int* n,
  const int* ilo, const int* ihi, double* H, const int* desch, 
  double* wr, double* wi, const int* iloz, const int* ihiz,
  double* Z, const int* descz, double* work, const int* lwork, 
  int* iwork, const int* liwork, int* info, const int* reclevel );

// Pipelined QR algorithm with AED for small matrices
// --------------------------------------------------
void ELEM_SCALAPACK(pslaqr1)
( const ELEM_FORT_LOGICAL* wantt, const ELEM_FORT_LOGICAL* wantz, const int* n,
  const int* ilo, const int* ihi, float* H, const int* desch, 
  float* wr, float* wi, const int* iloz, const int* ihiz, 
  float* Z, const int* descz, float* work, const int* lwork, 
  int* iwork, const int* liwork, int* info );
void ELEM_SCALAPACK(pdlaqr1)
( const ELEM_FORT_LOGICAL* wantt, const ELEM_FORT_LOGICAL* wantz, const int* n,
  const int* ilo, const int* ihi, double* H, const int* desch, 
  double* wr, double* wi, const int* iloz, const int* ihiz, 
  double* Z, const int* descz, double* work, const int* lwork, 
  int* iwork, const int* liwork, int* info );

} // extern "C"

namespace elem {

namespace blacs {

// BLACS
// =====

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

} // namespace blacs

namespace scalapack {
// ScaLAPACK
// =========

// Hessenberg QR algorithm
// -----------------------

void HessenbergSchur( int n, float* H, const int* desch, scomplex* w ) 
{
    DEBUG_ONLY(CallStackEntry cse("scalapack::HessenbergSchur"))
    const char job='E', compz='N';
    const int ilo=1, ihi=n;

    // Query the workspace sizes
    int fakeLDim=1, lwork=-1, dummyIWork, liwork=-1, info;
    int descz[9] = 
        { 1, desch[1], 0, 0, desch[4], desch[5], desch[6], desch[7], 1 };
    float dummyWork;
    std::vector<float> wr(n), wi(n);
    ELEM_SCALAPACK(pshseqr)
    ( &job, &compz, &n, &ilo, &ihi, H, desch, wr.data(), wi.data(), 0, descz,
      &dummyWork, &lwork, &dummyIWork, &liwork, &info );

    // Compute the eigenvalues in parallel
    lwork = dummyWork;
    liwork = dummyIWork;
    std::vector<float> work(lwork);
    std::vector<int> iwork(liwork);
    ELEM_SCALAPACK(pshseqr)
    ( &job, &compz, &n, &ilo, &ihi, H, desch, wr.data(), wi.data(), 0, descz,
      work.data(), &lwork, iwork.data(), &liwork, &info );

    // Combine the real and imaginary components of the eigenvalues
    for( int j=0; j<n; ++j )
        w[j] = scomplex(wr[j],wi[j]);
}

void HessenbergSchur( int n, double* H, const int* desch, dcomplex* w ) 
{
    DEBUG_ONLY(CallStackEntry cse("scalapack::HessenbergSchur"))
    const char job = 'E', compz='N';
    const int ilo=1, ihi=n;

    // Query the workspace sizes
    int fakeLDim=1, lwork=-1, dummyIWork, liwork=-1, info;
    int descz[9] = 
        { 1, desch[1], 0, 0, desch[4], desch[5], desch[6], desch[7], 1 };
    double dummyWork;
    std::vector<double> wr(n), wi(n);
    ELEM_SCALAPACK(pdhseqr)
    ( &job, &compz, &n, &ilo, &ihi, H, desch, wr.data(), wi.data(), 0, descz,
      &dummyWork, &lwork, &dummyIWork, &liwork, &info );

    // Compute the eigenvalues in parallel
    lwork = dummyWork;
    liwork = dummyIWork;
    std::vector<double> work(lwork);
    std::vector<int> iwork(liwork);
    ELEM_SCALAPACK(pdhseqr)
    ( &job, &compz, &n, &ilo, &ihi, H, desch, wr.data(), wi.data(), 0, descz,
      work.data(), &lwork, iwork.data(), &liwork, &info );

    // Combine the real and imaginary components of the eigenvalues
    for( int j=0; j<n; ++j )
        w[j] = dcomplex(wr[j],wi[j]);
}

void HessenbergSchur( int n, scomplex* H, const int* desch, scomplex* w ) 
{
    DEBUG_ONLY(CallStackEntry cse("scalapack::HessenbergSchur"))
    ELEM_FORT_LOGICAL wantt=ELEM_FORT_FALSE, wantz=ELEM_FORT_FALSE;
    const int ilo=1, ihi=n;

    // Query the workspace sizes
    int fakeLDim=1, lwork=-1, dummyIWork, liwork=-1, info;
    int descz[9] = 
        { 1, desch[1], 0, 0, desch[4], desch[5], desch[6], desch[7], 1 };
    scomplex dummyWork;
    ELEM_SCALAPACK(pclahqr)
    ( &wantt, &wantz, &n, &ilo, &ihi, H, desch, w, &ilo, &ihi, 0, descz,
      &dummyWork, &lwork, &dummyIWork, &liwork, &info );

    // Compute the eigenvalues in parallel
    lwork = dummyWork.real();
    liwork = dummyIWork;
    std::vector<scomplex> work(lwork);
    std::vector<int> iwork(liwork);
    ELEM_SCALAPACK(pclahqr)
    ( &wantt, &wantz, &n, &ilo, &ihi, H, desch, w, &ilo, &ihi, 0, descz,
      work.data(), &lwork, iwork.data(), &liwork, &info );
}

void HessenbergSchur( int n, dcomplex* H, const int* desch, dcomplex* w ) 
{
    DEBUG_ONLY(CallStackEntry cse("scalapack::HessenbergSchur"))
    ELEM_FORT_LOGICAL wantt=ELEM_FORT_FALSE, wantz=ELEM_FORT_FALSE;
    const int ilo=1, ihi=n;

    // Query the workspace sizes
    int fakeLDim=1, lwork=-1, dummyIWork, liwork=-1, info;
    int descz[9] = 
        { 1, desch[1], 0, 0, desch[4], desch[5], desch[6], desch[7], 1 };
    dcomplex dummyWork;
    ELEM_SCALAPACK(pzlahqr)
    ( &wantt, &wantz, &n, &ilo, &ihi, H, desch, w, &ilo, &ihi, 0, descz,
      &dummyWork, &lwork, &dummyIWork, &liwork, &info );

    // Compute the eigenvalues in parallel
    lwork = dummyWork.real();
    liwork = dummyIWork;
    std::vector<dcomplex> work(lwork);
    std::vector<int> iwork(liwork);
    ELEM_SCALAPACK(pzlahqr)
    ( &wantt, &wantz, &n, &ilo, &ihi, H, desch, w, &ilo, &ihi, 0, descz,
      work.data(), &lwork, iwork.data(), &liwork, &info );
}

void HessenbergSchur
( int n, float* H, const int* desch, scomplex* w, float* U, const int* descu ) 
{
    DEBUG_ONLY(CallStackEntry cse("scalapack::HessenbergSchur"))
    LogicError("This routine not yet written");
}

void HessenbergSchur
( int n, double* H, const int* desch, dcomplex* w, 
  double* U, const int* descu ) 
{
    DEBUG_ONLY(CallStackEntry cse("scalapack::HessenbergSchur"))
    LogicError("This routine not yet written");
}

void HessenbergSchur
( int n, scomplex* H, const int* desch, scomplex* w, 
  scomplex* U, const int* descu ) 
{
    DEBUG_ONLY(CallStackEntry cse("scalapack::HessenbergSchur"))
    LogicError("This routine not yet written");
}

void HessenbergSchur
( int n, dcomplex* H, const int* desch, dcomplex* w, 
  dcomplex* U, const int* descu ) 
{
    DEBUG_ONLY(CallStackEntry cse("scalapack::HessenbergSchur"))
    LogicError("This routine not yet written");
}

} // namespace scalapack
} // namespace elem

#endif // ifdef ELEM_HAVE_SCALAPACK
