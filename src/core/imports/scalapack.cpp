/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El-lite.hpp>

#ifdef EL_HAVE_SCALAPACK

using El::FortranLogical;
using El::scomplex;
using El::dcomplex;

extern "C" {

// Factorizations
// ==============

// Cholesky
// --------
void EL_SCALAPACK(pspotrf)
( const char* uplo,
  const int* n,
  float* A, const int* iA, const int* jA, const int* descA,
  int* info );
void EL_SCALAPACK(pdpotrf)
( const char* uplo,
  const int* n,
  double* A, const int* iA, const int* jA, const int* descA,
  int* info );
void EL_SCALAPACK(pcpotrf)
( const char* uplo,
  const int* n,
  scomplex* A, const int* iA, const int* jA, const int* descA,
  int* info );
void EL_SCALAPACK(pzpotrf)
( const char* uplo,
  const int* n,
  dcomplex* A, const int* iA, const int* jA, const int* descA,
  int* info );

// QR
// --
void EL_SCALAPACK(psgeqrf)
( const int* m, const int* n,
  float* A, const int* iA, const int* jA, const int* descA,
  float* tau,
  float* work, const int* lwork,
  int* info );
void EL_SCALAPACK(pdgeqrf)
( const int* m, const int* n,
  double* A, const int* iA, const int* jA, const int* descA,
  double* tau,
  double* work, const int* lwork,
  int* info );
void EL_SCALAPACK(pcgeqrf)
( const int* m, const int* n,
  scomplex* A, const int* iA, const int* jA, const int* descA,
  scomplex* tau,
  scomplex* work, const int* lwork,
  int* info );
void EL_SCALAPACK(pzgeqrf)
( const int* m, const int* n,
  dcomplex* A, const int* iA, const int* jA, const int* descA,
  dcomplex* tau,
  dcomplex* work, const int* lwork,
  int* info );

// Solvers
// =======

// General linear solver
// ---------------------
void EL_SCALAPACK(psgesv)
( const int* n, const int* numRhs,
  float* A, const int* iA, const int* jA, const int* descA,
  int* ipiv,
  float* B, const int* iB, const int* jB, const int* descB,
  int* info );
void EL_SCALAPACK(pdgesv)
( const int* n, const int* numRhs,
  double* A, const int* iA, const int* jA, const int* descA,
  int* ipiv,
  double* B, const int* iB, const int* jB, const int* descB,
  int* info );
void EL_SCALAPACK(pcgesv)
( const int* n, const int* numRhs,
  scomplex* A, const int* iA, const int* jA, const int* descA,
  int* ipiv,
  scomplex* B, const int* iB, const int* jB, const int* descB,
  int* info );
void EL_SCALAPACK(pzgesv)
( const int* n, const int* numRhs,
  dcomplex* A, const int* iA, const int* jA, const int* descA,
  int* ipiv,
  dcomplex* B, const int* iB, const int* jB, const int* descB,
  int* info );

// Spectral analysis
// =================

// Hermitian Eig
// -------------

// MRRR
// ^^^^
void EL_SCALAPACK(pssyevr)
( const char* jobZ,
  const char* range,
  const char* uplo,
  const int* n,
  float* A, const int* iA, const int* jA, const int* descA,
  const float* vl, const float* vu,
  const int* il, const int* iu,
  int* numCompEigVals,
  int* numCompEigVecs,
  float* w,
  float* Z, const int* iZ, const int* jZ, const int* descZ,
  float* work, const int* lwork,
  int* iWork, const int* liWork,
  int* info );
void EL_SCALAPACK(pdsyevr)
( const char* jobZ,
  const char* range,
  const char* uplo,
  const int* n,
  double* A, const int* iA, const int* jA, const int* descA,
  const double* vl, const double* vu,
  const int* il, const int* iu,
  int* numCompEigVals,
  int* numCompEigVecs,
  double* w,
  double* Z, const int* iZ, const int* jZ, const int* descZ,
  double* work, const int* lwork,
  int* iWork, const int* liWork,
  int* info );
void EL_SCALAPACK(pcheevr)
( const char* jobZ,
  const char* range,
  const char* uplo,
  const int* n,
  scomplex* A, const int* iA, const int* jA, const int* descA,
  const float* vl, const float* vu,
  const int* il, const int* iu,
  int* numCompEigVals,
  int* numCompEigVecs,
  float* w,
  scomplex* Z, const int* iZ, const int* jZ, const int* descZ,
  scomplex* work, const int* lwork,
  float* realWork, const int* lrWork,
  int* iWork, const int* liWork,
  int* info );
void EL_SCALAPACK(pzheevr)
( const char* jobZ,
  const char* range,
  const char* uplo,
  const int* n,
  dcomplex* A, const int* iA, const int* jA, const int* descA,
  const double* vl, const double* vu,
  const int* il, const int* iu,
  int* numCompEigVals,
  int* numCompEigVecs,
  double* w,
  dcomplex* Z, const int* iZ, const int* jZ, const int* descZ,
  dcomplex* work, const int* lwork,
  double* realWork, const int* lrWork,
  int* iWork, const int* liWork,
  int* info );

// SVD
// ---
void EL_SCALAPACK(psgesvd)
( const char* jobU,
  const char* jobVH,
  const int* m, const int* n,
  float* A, const int* iA, const int* jA, const int* descA,
  float* S,
  float* U, const int* iU, const int* jU, const int* descU,
  float* VH, const int* iVH, const int* jVH, const int* descVH,
  float* work, const int* lwork,
  int* info );
void EL_SCALAPACK(pdgesvd)
( const char* jobU,
  const char* jobVH,
  const int* m, const int* n,
  double* A, const int* iA, const int* jA, const int* descA,
  double* S,
  double* U, const int* iU, const int* jU, const int* descU,
  double* VH, const int* iVH, const int* jVH, const int* descVH,
  double* work, const int* lwork,
  int* info );
void EL_SCALAPACK(pcgesvd)
( const char* jobU,
  const char* jobVH,
  const int* m, const int* n,
  scomplex* A, const int* iA, const int* jA, const int* descA,
  float* S,
  scomplex* U, const int* iU, const int* jU, const int* descU,
  scomplex* VH, const int* iVH, const int* jVH, const int* descVH,
  scomplex* work, const int* lwork,
  float* rwork,
  int* info );
void EL_SCALAPACK(pzgesvd)
( const char* jobU,
  const char* jobVH,
  const int* m, const int* n,
  dcomplex* A, const int* iA, const int* jA, const int* descA,
  double* S,
  dcomplex* U, const int* iU, const int* jU, const int* descU,
  dcomplex* VH, const int* iVH, const int* jVH, const int* descVH,
  dcomplex* work, const int* lwork,
  double* rwork,
  int* info );

// Reduction of a Hermitian positive-definite EVP to standard form
// ---------------------------------------------------------------
void EL_SCALAPACK(pssyngst)
( const int* typeB,
  const char* uplo,
  const int* n, 
        float* A, const int* iA, const int* jA, const int* descA,
  const float* B, const int* iB, const int* jB, const int* descB,
        float* scale,
        float* work, const int* workSize,
        int* info );
void EL_SCALAPACK(pdsyngst)
( const int* typeB,
  const char* uplo,
  const int* n, 
        double* A, const int* iA, const int* jA, const int* descA,
  const double* B, const int* iB, const int* jB, const int* descB,
        double* scale,
        double* work, const int* workSize,
        int* info );
void EL_SCALAPACK(pchengst)
( const int* typeB, const char* uplo, const int* n, 
        scomplex* A, const int* iA, const int* jA, const int* descA,
  const scomplex* B, const int* iB, const int* jB, const int* descB,
        float* scale,
        scomplex* work, const int* workSize,
        int* info );
void EL_SCALAPACK(pzhengst)
( const int* typeB, const char* uplo, const int* n, 
        dcomplex* A, const int* iA, const int* jA, const int* descA,
  const dcomplex* B, const int* iB, const int* jB, const int* descB,
        double* scale,
        dcomplex* work, const int* workSize,
        int* info );

// Hessenberg QR algorithm
// -----------------------

// Aggressive Early Deflation
// ^^^^^^^^^^^^^^^^^^^^^^^^^^
// NOTE: ScaLAPACK currently only supports AED for real matrices
void EL_SCALAPACK(pshseqr)
( const char* job,
  const char* compz, 
  const int* n,
  const int* ilo, const int* ihi, 
  float* H, const int* descH,
  float* wr, float* wi, 
  float* Q, const int* descQ,
  float* work, const int* workSize, 
  int* iWork, const int* iWorkSize,
  int* info );
void EL_SCALAPACK(pdhseqr)
( const char* job,
  const char* compz, 
  const int* n,
  const int* ilo, const int* ihi, 
  double* H, const int* descH,
  double* wr, double* wi, 
  double* Q, const int* descQ,
  double* work, const int* workSize, 
  int* iWork, const int* iWorkSize,
  int* info );

// Pipelined without AED
// ^^^^^^^^^^^^^^^^^^^^^
void EL_SCALAPACK(pslahqr)
( const FortranLogical* wantt,
  const FortranLogical* wantz,
  const int* n,
  const int* ilo, const int* ihi,
  float* H, const int* descH,
  float* wr, float* wi,
  const int* iloQ, const int* ihiQ, 
  float* Q, const int* descQ,
  float* work, const int* workSize,
  int* iWork, const int* iWorkSize, 
  int* info );
void EL_SCALAPACK(pdlahqr)
( const FortranLogical* wantt,
  const FortranLogical* wantz,
  const int* n,
  const int* ilo, const int* ihi,
  double* H, const int* descH,
  double* wr, double* wi,
  const int* iloQ, const int* ihiQ, 
  double* Q, const int* descQ,
  double* work, const int* workSize,
  int* iWork, const int* iWorkSize, 
  int* info );
void EL_SCALAPACK(pclahqr)
( const FortranLogical* wantt,
  const FortranLogical* wantz,
  const int* n,
  const int* ilo, const int* ihi,
  scomplex* H, const int* descH,
  scomplex* w,
  const int* iloQ, const int* ihiQ,
  scomplex* Q, const int* descQ,
  scomplex* work, const int* workSize,
  int* iWork, const int* iWorkSize, 
  int* info );
void EL_SCALAPACK(pzlahqr)
( const FortranLogical* wantt,
  const FortranLogical* wantz,
  const int* n,
  const int* ilo, const int* ihi,
  dcomplex* H, const int* descH,
  dcomplex* w,
  const int* iloQ, const int* ihiQ,
  dcomplex* Q, const int* descQ,
  dcomplex* work, const int* workSize,
  int* iWork, const int* iWorkSize, 
  int* info );

// Pipelined with AED for big matrices
// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
void EL_SCALAPACK(pslaqr0)
( const FortranLogical* wantt,
  const FortranLogical* wantz,
  const int* n,
  const int* ilo, const int* ihi,
  float* H, const int* descH, 
  float* wr, float* wi,
  const int* iloQ, const int* ihiQ, 
  float* Q, const int* descQ,
  float* work, const int* workSize, 
  int* iWork, const int* iWorkSize,
  int* info,
  const int* reclevel );
void EL_SCALAPACK(pdlaqr0)
( const FortranLogical* wantt,
  const FortranLogical* wantz,
  const int* n,
  const int* ilo, const int* ihi,
  double* H, const int* descH, 
  double* wr, double* wi,
  const int* iloQ, const int* ihiQ,
  double* Q, const int* descQ,
  double* work, const int* workSize, 
  int* iWork, const int* iWorkSize,
  int* info,
  const int* reclevel );

// Pipelined with AED for small matrices
// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
void EL_SCALAPACK(pslaqr1)
( const FortranLogical* wantt,
  const FortranLogical* wantz,
  const int* n,
  const int* ilo, const int* ihi,
  float* H, const int* descH, 
  float* wr, float* wi,
  const int* iloQ, const int* ihiQ, 
  float* Q, const int* descQ,
  float* work, const int* workSize, 
  int* iWork, const int* iWorkSize,
  int* info );
void EL_SCALAPACK(pdlaqr1)
( const FortranLogical* wantt,
  const FortranLogical* wantz,
  const int* n,
  const int* ilo, const int* ihi,
  double* H, const int* descH, 
  double* wr, double* wi,
  const int* iloQ, const int* ihiQ, 
  double* Q, const int* descQ,
  double* work, const int* workSize, 
  int* iWork, const int* iWorkSize,
  int* info );

} // extern "C"

namespace El {
namespace scalapack {

// Factorization
// =============

// Cholesky
// --------
void Cholesky( char uplo, int n, float* A, const int* descA )
{
    EL_DEBUG_CSE
    int iA=1,jA=1,info;
    EL_SCALAPACK(pspotrf)( &uplo, &n, A, &iA, &jA, descA, &info );
    if( info != 0 )
        RuntimeError("pspotrf returned with info=",info);
}

void Cholesky( char uplo, int n, double* A, const int* descA )
{
    EL_DEBUG_CSE
    int iA=1,jA=1,info;
    EL_SCALAPACK(pdpotrf)( &uplo, &n, A, &iA, &jA, descA, &info );
    if( info != 0 )
        RuntimeError("pdpotrf returned with info=",info);
}

void Cholesky( char uplo, int n, scomplex* A, const int* descA )
{
    EL_DEBUG_CSE
    int iA=1,jA=1,info;
    EL_SCALAPACK(pcpotrf)( &uplo, &n, A, &iA, &jA, descA, &info );
    if( info != 0 )
        RuntimeError("pcpotrf returned with info=",info);
}

void Cholesky( char uplo, int n, dcomplex* A, const int* descA )
{
    EL_DEBUG_CSE
    int iA=1,jA=1,info;
    EL_SCALAPACK(pzpotrf)( &uplo, &n, A, &iA, &jA, descA, &info );
    if( info != 0 )
        RuntimeError("pzpotrf returned with info=",info);
}

// QR
// --
void QR( int m, int n, float* A, const int* descA, float* tau )
{
    EL_DEBUG_CSE
    int iA=1, jA=1, info;

    int lwork=-1;
    float dummyWork;
    EL_SCALAPACK(psgeqrf)
    ( &m, &n, A, &iA, &jA, descA, tau, &dummyWork, &lwork, &info );
    
    lwork = dummyWork;
    vector<float> work(lwork);
    EL_SCALAPACK(psgeqrf)
    ( &m, &n, A, &iA, &jA, descA, tau, work.data(), &lwork, &info );
    if( info != 0 )
        RuntimeError("psgeqrf returned with info=",info);
}

void QR( int m, int n, double* A, const int* descA, double* tau )
{
    EL_DEBUG_CSE
    int iA=1, jA=1, info;

    int lwork=-1;
    double dummyWork;
    EL_SCALAPACK(pdgeqrf)
    ( &m, &n, A, &iA, &jA, descA, tau, &dummyWork, &lwork, &info );
    
    lwork = dummyWork;
    vector<double> work(lwork);
    EL_SCALAPACK(pdgeqrf)
    ( &m, &n, A, &iA, &jA, descA, tau, work.data(), &lwork, &info );
    if( info != 0 )
        RuntimeError("pdgeqrf returned with info=",info);
}

void QR( int m, int n, scomplex* A, const int* descA, scomplex* tau )
{
    EL_DEBUG_CSE
    int iA=1, jA=1, info;

    int lwork=-1;
    scomplex dummyWork;
    EL_SCALAPACK(pcgeqrf)
    ( &m, &n, A, &iA, &jA, descA, tau, &dummyWork, &lwork, &info );
    
    lwork = dummyWork.real();
    vector<scomplex> work(lwork);
    EL_SCALAPACK(pcgeqrf)
    ( &m, &n, A, &iA, &jA, descA, tau, work.data(), &lwork, &info );
    if( info != 0 )
        RuntimeError("pcgeqrf returned with info=",info);
}

void QR( int m, int n, dcomplex* A, const int* descA, dcomplex* tau )
{
    EL_DEBUG_CSE
    int iA=1, jA=1, info;

    int lwork=-1;
    dcomplex dummyWork;
    EL_SCALAPACK(pzgeqrf)
    ( &m, &n, A, &iA, &jA, descA, tau, &dummyWork, &lwork, &info );
    
    lwork = dummyWork.real();
    vector<dcomplex> work(lwork);
    EL_SCALAPACK(pzgeqrf)
    ( &m, &n, A, &iA, &jA, descA, tau, work.data(), &lwork, &info );
    if( info != 0 )
        RuntimeError("pzgeqrf returned with info=",info);
}

// Solvers
// =======

// General linear solver
// ---------------------

void LinearSolve
( int n, int numRhs,
  float* A, const int* descA,
  int* ipiv,
  float* B, const int* descB )
{
    EL_DEBUG_CSE
    int iA=1, jA=1, iB=1, jB=1, info;
    EL_SCALAPACK(psgesv)
    ( &n, &numRhs, A, &iA, &jA, descA, ipiv, B, &iB, &jB, descB, &info );
    if( info != 0 )
        RuntimeError("psgesv returned with info=",info);
}

void LinearSolve
( int n, int numRhs,
  double* A, const int* descA,
  int* ipiv,
  double* B, const int* descB )
{
    EL_DEBUG_CSE
    int iA=1, jA=1, iB=1, jB=1, info;
    EL_SCALAPACK(pdgesv)
    ( &n, &numRhs, A, &iA, &jA, descA, ipiv, B, &iB, &jB, descB, &info );
    if( info != 0 )
        RuntimeError("psgesv returned with info=",info);
}

void LinearSolve
( int n, int numRhs,
  scomplex* A, const int* descA,
  int* ipiv,
  scomplex* B, const int* descB )
{
    EL_DEBUG_CSE
    int iA=1, jA=1, iB=1, jB=1, info;
    EL_SCALAPACK(pcgesv)
    ( &n, &numRhs, A, &iA, &jA, descA, ipiv, B, &iB, &jB, descB, &info );
    if( info != 0 )
        RuntimeError("psgesv returned with info=",info);
}

void LinearSolve
( int n, int numRhs,
  dcomplex* A, const int* descA,
  int* ipiv,
  dcomplex* B, const int* descB )
{
    EL_DEBUG_CSE
    int iA=1, jA=1, iB=1, jB=1, info;
    EL_SCALAPACK(pzgesv)
    ( &n, &numRhs, A, &iA, &jA, descA, ipiv, B, &iB, &jB, descB, &info );
    if( info != 0 )
        RuntimeError("psgesv returned with info=",info);
}

// Spectral analysis
// =================

// Hermitian Eig
// -------------

// All eigenvalues
// ^^^^^^^^^^^^^^^
void HermitianEig
( char uplo, int n, float* A, const int* descA, float* w )
{
    EL_DEBUG_CSE
    char jobZ='N', range='A';
    int iA=1, jA=1, iZ=1, jZ=1,
        iL, iU,
        numCompEigVals, numCompEigVecs,
        lWork, liWork,
        info;
    float vL, vU;
    float* Z=nullptr;

    // Workspace query
    int dummyIWork;
    float dummyWork;
    lWork = liWork = -1;
    EL_SCALAPACK(pssyevr)
    ( &jobZ, &range, &uplo, &n,
      A, &iA, &jA, descA,
      &vL, &vU, &iL, &iU,
      &numCompEigVals, &numCompEigVecs,
      w,
      Z, &iZ, &jZ, descA,
      &dummyWork, &lWork,
      &dummyIWork, &liWork,
      &info );

    // Actual call
    lWork = dummyWork;
    liWork = dummyIWork;
    vector<float> work(lWork);
    vector<int> iWork(liWork);
    EL_SCALAPACK(pssyevr)
    ( &jobZ, &range, &uplo, &n,
      A, &iA, &jA, descA,
      &vL, &vU, &iL, &iU,
      &numCompEigVals, &numCompEigVecs,
      w,
      Z, &iZ, &jZ, descA,
      work.data(), &lWork,
      iWork.data(), &liWork,
      &info );
    if( info != 0 )
        RuntimeError("pssyevr exited with info=",info);
}

void HermitianEig
( char uplo, int n, double* A, const int* descA, double* w )
{
    EL_DEBUG_CSE
    char jobZ='N', range='A';
    int iA=1, jA=1, iZ=1, jZ=1,
        iL, iU,
        numCompEigVals, numCompEigVecs,
        lWork, liWork,
        info;
    double vL, vU;
    double* Z=nullptr;

    // Workspace query
    int dummyIWork;
    double dummyWork;
    lWork = liWork = -1;
    EL_SCALAPACK(pdsyevr)
    ( &jobZ, &range, &uplo, &n,
      A, &iA, &jA, descA,
      &vL, &vU, &iL, &iU,
      &numCompEigVals, &numCompEigVecs,
      w,
      Z, &iZ, &jZ, descA,
      &dummyWork, &lWork,
      &dummyIWork, &liWork,
      &info );

    // Actual call
    lWork = dummyWork;
    liWork = dummyIWork;
    vector<double> work(lWork);
    vector<int> iWork(liWork);
    EL_SCALAPACK(pdsyevr)    
    ( &jobZ, &range, &uplo, &n,
      A, &iA, &jA, descA,
      &vL, &vU, &iL, &iU,
      &numCompEigVals, &numCompEigVecs,
      w,
      Z, &iZ, &jZ, descA,
      work.data(), &lWork,
      iWork.data(), &liWork,
      &info );
    if( info != 0 )
        RuntimeError("pdsyevr exited with info=",info);
}

void HermitianEig
( char uplo, int n, scomplex* A, const int* descA, float* w )
{
    EL_DEBUG_CSE
    char jobZ='N', range='A';
    int iA=1, jA=1, iZ=1, jZ=1,
        iL, iU,
        numCompEigVals, numCompEigVecs,
        lWork, lrWork, liWork,
        info;
    float vL, vU;
    scomplex* Z=nullptr;

    // Workspace query
    int dummyIWork;
    scomplex dummyWork;
    float dummyRWork;
    lWork = lrWork = liWork = -1;
    EL_SCALAPACK(pcheevr)
    ( &jobZ, &range, &uplo, &n,
      A, &iA, &jA, descA,
      &vL, &vU, &iL, &iU,
      &numCompEigVals, &numCompEigVecs,
      w,
      Z, &iZ, &jZ, descA,
      &dummyWork, &lWork,
      &dummyRWork, &lrWork,
      &dummyIWork, &liWork,
      &info );

    // Actual call
    lWork = dummyWork.real();
    lrWork = dummyRWork;
    liWork = dummyIWork;
    vector<scomplex> work(lWork);
    vector<float> rWork(lrWork);
    vector<int> iWork(liWork);
    EL_SCALAPACK(pcheevr)
    ( &jobZ, &range, &uplo, &n,
      A, &iA, &jA, descA,
      &vL, &vU, &iL, &iU,
      &numCompEigVals, &numCompEigVecs,
      w,
      Z, &iZ, &jZ, descA,
      work.data(), &lWork,
      rWork.data(), &lrWork,
      iWork.data(), &liWork,
      &info );
    if( info != 0 )
        RuntimeError("pcheevr exited with info=",info);
}

void HermitianEig
( char uplo, int n, dcomplex* A, const int* descA, double* w )
{
    EL_DEBUG_CSE
    char jobZ='N', range='A';
    int iA=1, jA=1, iZ=1, jZ=1,
        iL, iU,
        numCompEigVals, numCompEigVecs,
        lWork, lrWork, liWork,
        info;
    double vL, vU;
    dcomplex* Z=nullptr;

    // Workspace query
    int dummyIWork;
    dcomplex dummyWork;
    double dummyRWork;
    lWork = lrWork = liWork = -1;
    EL_SCALAPACK(pzheevr)
    ( &jobZ, &range, &uplo, &n,
      A, &iA, &jA, descA,
      &vL, &vU, &iL, &iU,
      &numCompEigVals, &numCompEigVecs,
      w,
      Z, &iZ, &jZ, descA,
      &dummyWork, &lWork,
      &dummyRWork, &lrWork,
      &dummyIWork, &liWork,
      &info );

    // Actual call
    lWork = dummyWork.real();
    lrWork = dummyRWork;
    liWork = dummyIWork;
    vector<dcomplex> work(lWork);
    vector<double> rWork(lrWork);
    vector<int> iWork(liWork);
    EL_SCALAPACK(pzheevr)
    ( &jobZ, &range, &uplo, &n,
      A, &iA, &jA, descA,
      &vL, &vU, &iL, &iU,
      &numCompEigVals, &numCompEigVecs,
      w,
      Z, &iZ, &jZ, descA,
      work.data(), &lWork,
      rWork.data(), &lrWork,
      iWork.data(), &liWork,
      &info );
    if( info != 0 )
        RuntimeError("pzheevr exited with info=",info);
}

// All eigenpairs
// ^^^^^^^^^^^^^^
void HermitianEig
( char uplo, int n,
  float* A, const int* descA,
  float* w,
  float* Z, const int* descZ )
{
    EL_DEBUG_CSE
    char jobZ='V', range='A';
    int iA=1, jA=1, iZ=1, jZ=1,
        iL, iU,
        numCompEigVals, numCompEigVecs,
        lWork, liWork,
        info;
    float vL, vU;

    // Workspace query
    int dummyIWork;
    float dummyWork;
    lWork = liWork = -1;
    EL_SCALAPACK(pssyevr)
    ( &jobZ, &range, &uplo, &n,
      A, &iA, &jA, descA,
      &vL, &vU, &iL, &iU,
      &numCompEigVals, &numCompEigVecs,
      w,
      Z, &iZ, &jZ, descZ,
      &dummyWork, &lWork,
      &dummyIWork, &liWork,
      &info );

    // Actual call
    lWork = dummyWork;
    liWork = dummyIWork;
    vector<float> work(lWork);
    vector<int> iWork(liWork);
    EL_SCALAPACK(pssyevr)
    ( &jobZ, &range, &uplo, &n,
      A, &iA, &jA, descA,
      &vL, &vU, &iL, &iU,
      &numCompEigVals, &numCompEigVecs,
      w,
      Z, &iZ, &jZ, descZ,
      work.data(), &lWork,
      iWork.data(), &liWork,
      &info );
    if( info != 0 )
        RuntimeError("pssyevr exited with info=",info);
}

void HermitianEig
( char uplo, int n,
  double* A, const int* descA,
  double* w,
  double* Z, const int* descZ )
{
    EL_DEBUG_CSE
    char jobZ='V', range='A';
    int iA=1, jA=1, iZ=1, jZ=1,
        iL, iU,
        numCompEigVals, numCompEigVecs,
        lWork, liWork,
        info;
    double vL, vU;

    // Workspace query
    int dummyIWork;
    double dummyWork;
    lWork = liWork = -1;
    EL_SCALAPACK(pdsyevr)
    ( &jobZ, &range, &uplo, &n,
      A, &iA, &jA, descA,
      &vL, &vU, &iL, &iU,
      &numCompEigVals, &numCompEigVecs,
      w,
      Z, &iZ, &jZ, descZ,
      &dummyWork, &lWork,
      &dummyIWork, &liWork,
      &info );

    // Actual call
    lWork = dummyWork;
    liWork = dummyIWork;
    vector<double> work(lWork);
    vector<int> iWork(liWork);
    EL_SCALAPACK(pdsyevr)    
    ( &jobZ, &range, &uplo, &n,
      A, &iA, &jA, descA,
      &vL, &vU, &iL, &iU,
      &numCompEigVals, &numCompEigVecs,
      w,
      Z, &iZ, &jZ, descZ,
      work.data(), &lWork,
      iWork.data(), &liWork,
      &info );
    if( info != 0 )
        RuntimeError("pdsyevr exited with info=",info);
}

void HermitianEig
( char uplo, int n,
  scomplex* A, const int* descA,
  float* w,
  scomplex* Z, const int* descZ )
{
    EL_DEBUG_CSE
    char jobZ='V', range='A';
    int iA=1, jA=1, iZ=1, jZ=1,
        iL, iU,
        numCompEigVals, numCompEigVecs,
        lWork, lrWork, liWork,
        info;
    float vL, vU;

    // Workspace query
    int dummyIWork;
    scomplex dummyWork;
    float dummyRWork;
    lWork = lrWork = liWork = -1;
    EL_SCALAPACK(pcheevr)
    ( &jobZ, &range, &uplo, &n,
      A, &iA, &jA, descA,
      &vL, &vU, &iL, &iU,
      &numCompEigVals, &numCompEigVecs,
      w,
      Z, &iZ, &jZ, descZ,
      &dummyWork, &lWork,
      &dummyRWork, &lrWork,
      &dummyIWork, &liWork,
      &info );

    // Actual call
    lWork = dummyWork.real();
    lrWork = dummyRWork;
    liWork = dummyIWork;
    vector<scomplex> work(lWork);
    vector<float> rWork(lrWork);
    vector<int> iWork(liWork);
    EL_SCALAPACK(pcheevr)
    ( &jobZ, &range, &uplo, &n,
      A, &iA, &jA, descA,
      &vL, &vU, &iL, &iU,
      &numCompEigVals, &numCompEigVecs,
      w,
      Z, &iZ, &jZ, descZ,
      work.data(), &lWork,
      rWork.data(), &lrWork,
      iWork.data(), &liWork,
      &info );
    if( info != 0 )
        RuntimeError("pcheevr exited with info=",info);
}

void HermitianEig
( char uplo, int n,
  dcomplex* A, const int* descA,
  double* w,
  dcomplex* Z, const int* descZ )
{
    EL_DEBUG_CSE
    char jobZ='V', range='A';
    int iA=1, jA=1, iZ=1, jZ=1,
        iL, iU,
        numCompEigVals, numCompEigVecs,
        lWork, lrWork, liWork,
        info;
    double vL, vU;

    // Workspace query
    int dummyIWork;
    dcomplex dummyWork;
    double dummyRWork;
    lWork = lrWork = liWork = -1;
    EL_SCALAPACK(pzheevr)
    ( &jobZ, &range, &uplo, &n,
      A, &iA, &jA, descA,
      &vL, &vU, &iL, &iU,
      &numCompEigVals, &numCompEigVecs,
      w,
      Z, &iZ, &jZ, descZ,
      &dummyWork, &lWork,
      &dummyRWork, &lrWork,
      &dummyIWork, &liWork,
      &info );

    // Actual call
    lWork = dummyWork.real();
    lrWork = dummyRWork;
    liWork = dummyIWork;
    vector<dcomplex> work(lWork);
    vector<double> rWork(lrWork);
    vector<int> iWork(liWork);
    EL_SCALAPACK(pzheevr)
    ( &jobZ, &range, &uplo, &n,
      A, &iA, &jA, descA,
      &vL, &vU, &iL, &iU,
      &numCompEigVals, &numCompEigVecs,
      w,
      Z, &iZ, &jZ, descZ,
      work.data(), &lWork,
      rWork.data(), &lrWork,
      iWork.data(), &liWork,
      &info );
    if( info != 0 )
        RuntimeError("pzheevr exited with info=",info);
}

// SVD
// ---

// Singular values
// ^^^^^^^^^^^^^^^
void SingularValues( int m, int n, float* A, const int* descA, float* s )
{
    EL_DEBUG_CSE
    const char jobU='N', jobVH='N';
    int iA=1, jA=1, iU=1, jU=1, iVH=1, jVH=1, info;
    float* U=nullptr;
    float* VH=nullptr;

    // Workspace query
    float dummyWork;
    int lwork=-1;
    EL_SCALAPACK(psgesvd)
    ( &jobU, &jobVH, &m, &n,
      A, &iA, &jA, descA,
      s,
      U, &iU, &jU, descA,
      VH, &iVH, &jVH, descA,
      &dummyWork, &lwork,
      &info );

    // Actual call
    lwork = dummyWork;
    vector<float> work(lwork);
    EL_SCALAPACK(psgesvd)
    ( &jobU, &jobVH, &m, &n,
      A, &iA, &jA, descA,
      s,
      U, &iU, &jU, descA,
      VH, &iVH, &jVH, descA,
      work.data(), &lwork,
      &info );
    if( info != 0 )
        RuntimeError("psgesvd exited with info=",info);
}

void SingularValues( int m, int n, double* A, const int* descA, double* s )
{
    EL_DEBUG_CSE
    const char jobU='N', jobVH='N';
    int iA=1, jA=1, iU=1, jU=1, iVH=1, jVH=1, info;
    double* U=nullptr;
    double* VH=nullptr;

    // Workspace query
    double dummyWork;
    int lwork=-1;
    EL_SCALAPACK(pdgesvd)
    ( &jobU, &jobVH, &m, &n,
      A, &iA, &jA, descA,
      s,
      U, &iU, &jU, descA,
      VH, &iVH, &jVH, descA,
      &dummyWork, &lwork,
      &info );

    // Actual call
    lwork = dummyWork;
    vector<double> work(lwork);
    EL_SCALAPACK(pdgesvd)
    ( &jobU, &jobVH, &m, &n,
      A, &iA, &jA, descA,
      s,
      U, &iU, &jU, descA,
      VH, &iVH, &jVH, descA,
      work.data(), &lwork,
      &info );
    if( info != 0 )
        RuntimeError("pdgesvd exited with info=",info);
}

void SingularValues( int m, int n, scomplex* A, const int* descA, float* s )
{
    EL_DEBUG_CSE
    const char jobU='N', jobVH='N';
    int iA=1, jA=1, iU=1, jU=1, iVH=1, jVH=1, info;
    scomplex* U=nullptr;
    scomplex* VH=nullptr;

    const int maxDim = Max(m,n);
    vector<float> rwork(1+4*maxDim);

    // Workspace query
    scomplex dummyWork;
    int lwork=-1;
    EL_SCALAPACK(pcgesvd)
    ( &jobU, &jobVH, &m, &n,
      A, &iA, &jA, descA,
      s,
      U, &iU, &jU, descA,
      VH, &iVH, &jVH, descA,
      &dummyWork, &lwork,
      rwork.data(),
      &info );
    if( int(rwork[0]) > int(rwork.size()) )
    {
        if( mpi::Rank() == 0 )
            Output("WARNING: Resized rwork from ",rwork.size()," to ",rwork[0]);
        rwork.resize( int(rwork[0]) );
    }

    // Actual call
    lwork = dummyWork.real();
    vector<scomplex> work(lwork);
    EL_SCALAPACK(pcgesvd)
    ( &jobU, &jobVH, &m, &n,
      A, &iA, &jA, descA,
      s, U, &iU, &jU, descA,
      VH, &iVH, &jVH, descA,
      work.data(), &lwork,
      rwork.data(),
      &info );
    if( info != 0 )
        RuntimeError("pcgesvd exited with info=",info);
}

void SingularValues( int m, int n, dcomplex* A, const int* descA, double* s )
{
    EL_DEBUG_CSE
    const char jobU='N', jobVH='N';
    int iA=1, jA=1, iU=1, jU=1, iVH=1, jVH=1, info;
    dcomplex* U=nullptr;
    dcomplex* VH=nullptr;

    const int maxDim = Max(m,n);
    vector<double> rwork(1+4*maxDim);

    // Workspace query
    dcomplex dummyWork;
    int lwork=-1;
    EL_SCALAPACK(pzgesvd)
    ( &jobU, &jobVH, &m, &n,
      A, &iA, &jA, descA,
      s,
      U, &iU, &jU, descA,
      VH, &iVH, &jVH, descA,
      &dummyWork, &lwork,
      rwork.data(),
      &info );
    if( int(rwork[0]) > int(rwork.size()) )
    {
        if( mpi::Rank() == 0 )
            Output("WARNING: Resized rwork from ",rwork.size()," to ",rwork[0]);
        rwork.resize( int(rwork[0]) );
    }

    // Actual call
    lwork = dummyWork.real();
    vector<dcomplex> work(lwork);
    EL_SCALAPACK(pzgesvd)
    ( &jobU, &jobVH, &m, &n,
      A, &iA, &jA, descA,
      s,
      U, &iU, &jU, descA,
      VH, &iVH, &jVH, descA,
      work.data(), &lwork,
      rwork.data(),
      &info );
    if( info != 0 )
        RuntimeError("pzgesvd exited with info=",info);
}

// Values and vectors
// ^^^^^^^^^^^^^^^^^^
void SVD
( int m, int n,
  float* A, const int* descA,
  float* s,
  float* U, const int* descU,
  float* VH, const int* descVH )
{
    EL_DEBUG_CSE
    int iA=1, jA=1, iU=1, jU=1, iVH=1, jVH=1, info;

    // Workspace query
    const char jobU='V', jobVH='V';
    float dummyWork;
    int lwork=-1;
    EL_SCALAPACK(psgesvd)
    ( &jobU, &jobVH, &m, &n,
      A, &iA, &jA, descA,
      s,
      U, &iU, &jU, descU,
      VH, &iVH, &jVH, descVH,
      &dummyWork, &lwork,
      &info );

    // Actual call
    lwork = dummyWork;
    vector<float> work(lwork);
    EL_SCALAPACK(psgesvd)
    ( &jobU, &jobVH, &m, &n,
      A, &iA, &jA, descA,
      s,
      U, &iU, &jU, descU,
      VH, &iVH, &jVH, descVH,
      work.data(), &lwork,
      &info );
    if( info != 0 )
        RuntimeError("psgesvd exited with info=",info);
}

void SVD
( int m, int n,
  double* A, const int* descA,
  double* s,
  double* U, const int* descU,
  double* VH, const int* descVH )
{
    EL_DEBUG_CSE
    int iA=1, jA=1, iU=1, jU=1, iVH=1, jVH=1, info;

    // Workspace query
    const char jobU='V', jobVH='V';
    double dummyWork;
    int lwork=-1;
    EL_SCALAPACK(pdgesvd)
    ( &jobU, &jobVH, &m, &n,
      A, &iA, &jA, descA,
      s,
      U, &iU, &jU, descU,
      VH, &iVH, &jVH, descVH,
      &dummyWork, &lwork,
      &info );

    // Actual call
    lwork = dummyWork;
    vector<double> work(lwork);
    EL_SCALAPACK(pdgesvd)
    ( &jobU, &jobVH, &m, &n,
      A, &iA, &jA, descA,
      s,
      U, &iU, &jU, descU,
      VH, &iVH, &jVH, descVH,
      work.data(), &lwork,
      &info );
    if( info != 0 )
        RuntimeError("pdgesvd exited with info=",info);
}

void SVD
( int m, int n,
  scomplex* A, const int* descA,
  float* s,
  scomplex* U, const int* descU,
  scomplex* VH, const int* descVH )
{
    EL_DEBUG_CSE
    int iA=1, jA=1, iU=1, jU=1, iVH=1, jVH=1, info;

    const int maxDim = Max(m,n);
    vector<float> rwork(1+4*maxDim);

    const char jobU='V', jobVH='V';
    scomplex dummyWork;
    int lwork=-1;
    EL_SCALAPACK(pcgesvd)
    ( &jobU, &jobVH, &m, &n,
      A, &iA, &jA, descA,
      s,
      U, &iU, &jU, descU,
      VH, &iVH, &jVH, descVH,
      &dummyWork, &lwork, rwork.data(), &info );
    if( int(rwork[0]) > int(rwork.size()) )
    {
        if( mpi::Rank() == 0 )
            Output("WARNING: Resized rwork from ",rwork.size()," to ",rwork[0]);
        rwork.resize( int(rwork[0]) );
    }

    lwork = dummyWork.real();
    vector<scomplex> work(lwork);
    EL_SCALAPACK(pcgesvd)
    ( &jobU, &jobVH, &m, &n,
      A, &iA, &jA, descA,
      s,
      U, &iU, &jU, descU,
      VH, &iVH, &jVH, descVH,
      work.data(), &lwork, rwork.data(), &info );
}

void SVD
( int m, int n,
  dcomplex* A, const int* descA,
  double* s,
  dcomplex* U, const int* descU,
  dcomplex* VH, const int* descVH )
{
    EL_DEBUG_CSE
    int iA=1, jA=1, iU=1, jU=1, iVH=1, jVH=1, info;

    const int maxDim = Max(m,n);
    vector<double> rwork(1+4*maxDim);

    const char jobU='V', jobVH='V';
    dcomplex dummyWork;
    int lwork=-1;
    EL_SCALAPACK(pzgesvd)
    ( &jobU, &jobVH, &m, &n,
      A, &iA, &jA, descA,
      s,
      U, &iU, &jU, descU,
      VH, &iVH, &jVH, descVH,
      &dummyWork, &lwork, rwork.data(), &info );
    if( int(rwork[0]) > int(rwork.size()) )
    {
        if( mpi::Rank() == 0 )
            Output("WARNING: Resized rwork from ",rwork.size()," to ",rwork[0]);
        rwork.resize( int(rwork[0]) );
    }

    lwork = dummyWork.real();
    vector<dcomplex> work(lwork);
    EL_SCALAPACK(pzgesvd)
    ( &jobU, &jobVH, &m, &n,
      A, &iA, &jA, descA,
      s,
      U, &iU, &jU, descU,
      VH, &iVH, &jVH, descVH,
      work.data(), &lwork, rwork.data(), &info );
}

// Hermitian eigenvalue problems
// -----------------------------

// Compute eigenvalues
// ^^^^^^^^^^^^^^^^^^^

// All eigenvalues
// """""""""""""""
// TODO

// Floating-point range
// """"""""""""""""""""
// TODO

// Index range
// """""""""""
// TODO

// Compute eigenpairs
// ^^^^^^^^^^^^^^^^^^

// All eigenpairs
// """"""""""""""
// TODO

// Floating-point range
// """"""""""""""""""""
// TODO

// Index range
// """""""""""
// TODO

// Reduction of a generalized Hermitian positive-definite EVP to standard form
// ---------------------------------------------------------------------------
// NOTE: It is required that B have a positive diagonal

// Two-sided Trsm
// ^^^^^^^^^^^^^^
void TwoSidedTrsm
( char uplo, int n, 
        float* A, const int* descA, 
  const float* B, const int* descB )
{
    EL_DEBUG_CSE
    int typeB=1,iA=1,jA=1,iB=1,jB=1,workSize=-1,info;
    float scale, dummyWork;
    EL_SCALAPACK(pssyngst)
    ( &typeB, &uplo, &n,
      A, &iA, &jA, descA,
      B, &iB, &jB, descB,
      &scale, &dummyWork, &workSize, &info );

    workSize = dummyWork;
    vector<float> work( workSize );
    EL_SCALAPACK(pssyngst)
    ( &typeB, &uplo, &n,
      A, &iA, &jA, descA,
      B, &iB, &jB, descB,
      &scale, work.data(), &workSize, &info );
    if( info != 0 )
        RuntimeError("pssyngst exited with info=",info);
}

void TwoSidedTrsm
( char uplo, int n, 
        double* A, const int* descA, 
  const double* B, const int* descB )
{
    EL_DEBUG_CSE
    int typeB=1,iA=1,jA=1,iB=1,jB=1,workSize=-1,info;
    double scale, dummyWork;
    EL_SCALAPACK(pdsyngst)
    ( &typeB, &uplo, &n,
      A, &iA, &jA, descA,
      B, &iB, &jB, descB,
      &scale, &dummyWork, &workSize, &info );

    workSize = dummyWork;
    vector<double> work( workSize );
    EL_SCALAPACK(pdsyngst)
    ( &typeB, &uplo, &n,
      A, &iA, &jA, descA,
      B, &iB, &jB, descB,
      &scale, work.data(), &workSize, &info );
    if( info != 0 )
        RuntimeError("pdsyngst exited with info=",info);
}

void TwoSidedTrsm
( char uplo, int n, 
        scomplex* A, const int* descA, 
  const scomplex* B, const int* descB )
{
    EL_DEBUG_CSE
    int typeB=1,iA=1,jA=1,iB=1,jB=1,workSize=-1,info;
    float scale;
    scomplex dummyWork;
    EL_SCALAPACK(pchengst)
    ( &typeB, &uplo, &n,
      A, &iA, &jA, descA,
      B, &iB, &jB, descB,
      &scale, &dummyWork, &workSize, &info );

    workSize = dummyWork.real();
    vector<scomplex> work( workSize );
    EL_SCALAPACK(pchengst)
    ( &typeB, &uplo, &n,
      A, &iA, &jA, descA,
      B, &iB, &jB, descB,
      &scale, work.data(), &workSize, &info );
    if( info != 0 )
        RuntimeError("pchengst exited with info=",info);
}

void TwoSidedTrsm
( char uplo, int n, 
        dcomplex* A, const int* descA, 
  const dcomplex* B, const int* descB )
{
    EL_DEBUG_CSE
    int typeB=1,iA=1,jA=1,iB=1,jB=1,workSize=-1,info;
    double scale;
    dcomplex dummyWork;
    EL_SCALAPACK(pzhengst)
    ( &typeB, &uplo, &n,
      A, &iA, &jA, descA,
      B, &iB, &jB, descB,
      &scale, &dummyWork, &workSize, &info );

    workSize = dummyWork.real();
    vector<dcomplex> work( workSize );
    EL_SCALAPACK(pzhengst)
    ( &typeB, &uplo, &n,
      A, &iA, &jA, descA,
      B, &iB, &jB, descB,
      &scale, work.data(), &workSize, &info );
    if( info != 0 )
        RuntimeError("pzhengst exited with info=",info);
}

// Two-sided Trmm
// ^^^^^^^^^^^^^^
void TwoSidedTrmm
( char uplo, int n, 
        float* A, const int* descA, 
  const float* B, const int* descB )
{
    EL_DEBUG_CSE
    int typeB=2,iA=1,jA=1,iB=1,jB=1,workSize=-1,info;
    float scale, dummyWork;
    EL_SCALAPACK(pssyngst)
    ( &typeB, &uplo, &n,
      A, &iA, &jA, descA,
      B, &iB, &jB, descB,
      &scale, &dummyWork, &workSize, &info );

    workSize = dummyWork;
    vector<float> work( workSize );
    EL_SCALAPACK(pssyngst)
    ( &typeB, &uplo, &n,
      A, &iA, &jA, descA,
      B, &iB, &jB, descB,
      &scale, work.data(), &workSize, &info );
    if( info != 0 )
        RuntimeError("pssyngst exited with info=",info);
}

void TwoSidedTrmm
( char uplo, int n, 
        double* A, const int* descA, 
  const double* B, const int* descB )
{
    EL_DEBUG_CSE
    int typeB=2,iA=1,jA=1,iB=1,jB=1,workSize=-1,info;
    double scale, dummyWork;
    EL_SCALAPACK(pdsyngst)
    ( &typeB, &uplo, &n,
      A, &iA, &jA, descA,
      B, &iB, &jB, descB,
      &scale, &dummyWork, &workSize, &info );

    workSize = dummyWork;
    vector<double> work( workSize );
    EL_SCALAPACK(pdsyngst)
    ( &typeB, &uplo, &n,
      A, &iA, &jA, descA,
      B, &iB, &jB, descB,
      &scale, work.data(), &workSize, &info );
    if( info != 0 )
        RuntimeError("pdsyngst exited with info=",info);
}

void TwoSidedTrmm
( char uplo, int n, 
        scomplex* A, const int* descA, 
  const scomplex* B, const int* descB )
{
    EL_DEBUG_CSE
    int typeB=2,iA=1,jA=1,iB=1,jB=1,workSize=-1,info;
    float scale;
    scomplex dummyWork;
    EL_SCALAPACK(pchengst)
    ( &typeB, &uplo, &n,
      A, &iA, &jA, descA,
      B, &iB, &jB, descB,
      &scale, &dummyWork, &workSize, &info );

    workSize = dummyWork.real();
    vector<scomplex> work( workSize );
    EL_SCALAPACK(pchengst)
    ( &typeB, &uplo, &n,
      A, &iA, &jA, descA,
      B, &iB, &jB, descB,
      &scale, work.data(), &workSize, &info );
    if( info != 0 )
        RuntimeError("pchengst exited with info=",info);
}

void TwoSidedTrmm
( char uplo, int n, 
        dcomplex* A, const int* descA, 
  const dcomplex* B, const int* descB )
{
    EL_DEBUG_CSE
    int typeB=2,iA=1,jA=1,iB=1,jB=1,workSize=-1,info;
    double scale;
    dcomplex dummyWork;
    EL_SCALAPACK(pzhengst)
    ( &typeB, &uplo, &n,
      A, &iA, &jA, descA,
      B, &iB, &jB, descB,
      &scale, &dummyWork, &workSize, &info );

    workSize = dummyWork.real();
    vector<dcomplex> work( workSize );
    EL_SCALAPACK(pzhengst)
    ( &typeB, &uplo, &n,
      A, &iA, &jA, descA,
      B, &iB, &jB, descB,
      &scale, work.data(), &workSize, &info );
    if( info != 0 )
        RuntimeError("pzhengst exited with info=",info);
}

// Hessenberg Schur decomposition via the QR algorithm
// ---------------------------------------------------
void HessenbergSchur
( int n, float* H, const int* descH, scomplex* w, bool fullTriangle, bool aed ) 
{
    EL_DEBUG_CSE
    const int ilo=1, ihi=n;
    vector<float> wr(n), wi(n);
    int descQ[9] = 
        { 1, descH[1], 0, 0, descH[4], descH[5], descH[6], descH[7], 1 };
    int info;
    if( aed )
    {
        cerr << 
          "WARNING: PSHSEQR seems to have a bug in its eigenvalue reordering" 
          << endl;
        const char job=(fullTriangle ? 'S' : 'E'), compz='N';

        // Query the workspace sizes
        int workSize=-1, dummyIWork, iWorkSize=-1;
        float dummyWork;

        EL_SCALAPACK(pshseqr)
        ( &job, &compz, &n, &ilo, &ihi, H, descH, wr.data(), wi.data(), 
          0, descQ, &dummyWork, &workSize, &dummyIWork, &iWorkSize, &info );

        // Compute the eigenvalues in parallel
        workSize = dummyWork;
        iWorkSize = dummyIWork;
        vector<float> work(workSize);
        vector<int> iWork(iWorkSize);
        EL_SCALAPACK(pshseqr)
        ( &job, &compz, &n, &ilo, &ihi, H, descH, wr.data(), wi.data(), 
          0, descQ, work.data(), &workSize, iWork.data(), &iWorkSize, &info );
        if( info != 0 )
            RuntimeError("pshseqr exited with info=",info);
    }
    else
    {
        FortranLogical wantt=(fullTriangle ? FORTRAN_TRUE : FORTRAN_FALSE), 
                       wantz=FORTRAN_FALSE;

        // PSLAHQR does not support a workspace query and instead assumes
        // that the workspace is at least 
        //     3*N + max(2*max(ldH,ldQ)+2*localWidth,numRowBlocks)

        const int context  = descH[1];
        const int mb       = descH[4];
        const int nb       = descH[5];
        const int rowAlign = descH[7];
        const int ldH      = descH[8];
        const int ldQ      = descQ[8];

        const int rowCut = 0;
        const int rowRank = blacs::GridCol( context );
        const int rowStride = blacs::GridWidth( context );

        const Int localWidth =
            BlockedLength(n,rowRank,rowAlign,nb,rowCut,rowStride);
        const Int numRowBlocks = n/mb + 1;
        int workSize= 3*n + Max(2*Max(ldH,ldQ)+2*localWidth,numRowBlocks),
            iWorkSize=0;
        vector<float> work(workSize);
        vector<int> iWork(iWorkSize);

        // Compute the eigenvalues in parallel
        EL_SCALAPACK(pslahqr)
        ( &wantt, &wantz, &n, &ilo, &ihi, H, descH, wr.data(), wi.data(), 
          &ilo, &ihi, 0, descQ, work.data(), &workSize, 
          iWork.data(), &iWorkSize, &info );
        if( info != 0 )
            RuntimeError("pslahqr exited with info=",info);
    }
    // Combine the real and imaginary components of the eigenvalues
    for( int j=0; j<n; ++j )
        w[j] = std::complex<float>(wr[j],wi[j]);
}

void HessenbergSchur
( int n, double* H, const int* descH, dcomplex* w, bool fullTriangle, bool aed )
{
    EL_DEBUG_CSE
    const int ilo=1, ihi=n;
    vector<double> wr(n), wi(n);
    int descQ[9] = 
        { 1, descH[1], 0, 0, descH[4], descH[5], descH[6], descH[7], 1 };
    int info;
    if( aed )
    {
        cerr << 
          "WARNING: PDHSEQR seems to have a bug in its eigenvalue reordering" 
          << endl;
        const char job=(fullTriangle ? 'S' : 'E'), compz='N';

        // Query the workspace sizes
        int workSize=-1, dummyIWork, iWorkSize=-1;
        double dummyWork;
        EL_SCALAPACK(pdhseqr)
        ( &job, &compz, &n, &ilo, &ihi, H, descH, wr.data(), wi.data(), 
          0, descQ, &dummyWork, &workSize, &dummyIWork, &iWorkSize, &info );

        // Compute the eigenvalues in parallel
        workSize = dummyWork;
        iWorkSize = dummyIWork;
        vector<double> work(workSize);
        vector<int> iWork(iWorkSize);
        EL_SCALAPACK(pdhseqr)
        ( &job, &compz, &n, &ilo, &ihi, H, descH, wr.data(), wi.data(), 
          0, descQ, work.data(), &workSize, iWork.data(), &iWorkSize, &info );
        if( info != 0 )
            RuntimeError("pdhseqr exited with info=",info);
    }
    else
    {
        FortranLogical wantt=(fullTriangle ? FORTRAN_TRUE : FORTRAN_FALSE), 
                       wantz=FORTRAN_FALSE;

        // PDLAHQR does not support a workspace query and instead assumes
        // that the workspace is at least 
        //     3*N + max(2*max(ldH,ldQ)+2*localWidth,numRowBlocks)

        const int context  = descH[1];
        const int mb       = descH[4];
        const int nb       = descH[5];
        const int rowAlign = descH[7];
        const int ldH      = descH[8];
        const int ldQ      = descQ[8];

        const int rowCut = 0;
        const int rowRank = blacs::GridCol( context );
        const int rowStride = blacs::GridWidth( context );

        const Int localWidth = 
            BlockedLength(n,rowRank,rowAlign,nb,rowCut,rowStride);
        const Int numRowBlocks = n/mb + 1;
        int workSize= 3*n + Max(2*Max(ldH,ldQ)+2*localWidth,numRowBlocks),
            iWorkSize=0;
        vector<double> work(workSize);
        vector<int> iWork(iWorkSize);

        // Compute the eigenvalues in parallel
        EL_SCALAPACK(pdlahqr)
        ( &wantt, &wantz, &n, &ilo, &ihi, H, descH, wr.data(), wi.data(), 
          &ilo, &ihi, 0, descQ, work.data(), &workSize, 
          iWork.data(), &iWorkSize, &info );
        if( info != 0 )
            RuntimeError("pdlahqr exited with info=",info);
    }
    // Combine the real and imaginary components of the eigenvalues
    for( int j=0; j<n; ++j )
        w[j] = std::complex<double>(wr[j],wi[j]);
}

void HessenbergSchur
( int n, scomplex* H, const int* descH, scomplex* w, bool fullTriangle, 
  bool aed ) 
{
    EL_DEBUG_CSE
    FortranLogical wantt=(fullTriangle ? FORTRAN_TRUE : FORTRAN_FALSE), 
                   wantz=FORTRAN_FALSE;
    if( aed )
        LogicError("AED is not supported for complex matrices");
    const int ilo=1, ihi=n;

    // Query the workspace sizes
    // NOTE: dummyIWork will currently be left unmodified (hence the
    //       zero initialization)!
    int workSize=-1, dummyIWork=0, iWorkSize=-1, info;
    int descQ[9] = 
        { 1, descH[1], 0, 0, descH[4], descH[5], descH[6], descH[7], 1 };
    scomplex dummyWork;
    EL_SCALAPACK(pclahqr)
    ( &wantt, &wantz, &n, &ilo, &ihi, H, descH, w, &ilo, &ihi, 0, descQ,
      &dummyWork, &workSize, &dummyIWork, &iWorkSize, &info );

    // Compute the eigenvalues in parallel
    workSize = dummyWork.real();
    iWorkSize = dummyIWork;
    vector<scomplex> work(workSize);
    vector<int> iWork(iWorkSize);
    EL_SCALAPACK(pclahqr)
    ( &wantt, &wantz, &n, &ilo, &ihi, H, descH, w, &ilo, &ihi, 0, descQ,
      work.data(), &workSize, iWork.data(), &iWorkSize, &info );
    if( info != 0 )
        RuntimeError("pclahqr exited with info=",info);
}

void HessenbergSchur
( int n, dcomplex* H, const int* descH, dcomplex* w, bool fullTriangle,
  bool aed ) 
{
    EL_DEBUG_CSE
    FortranLogical wantt=(fullTriangle ? FORTRAN_TRUE : FORTRAN_FALSE), 
                   wantz=FORTRAN_FALSE;
    if( aed )
        LogicError("AED is not supported for complex matrices");
    const int ilo=1, ihi=n;

    // Query the workspace sizes
    // NOTE: dummyIWork will currently be left unmodified (hence the
    //       zero initialization)!
    int workSize=-1, dummyIWork=0, iWorkSize=-1, info;
    int descQ[9] = 
        { 1, descH[1], 0, 0, descH[4], descH[5], descH[6], descH[7], 1 };
    dcomplex dummyWork;
    EL_SCALAPACK(pzlahqr)
    ( &wantt, &wantz, &n, &ilo, &ihi, H, descH, w, &ilo, &ihi, 0, descQ,
      &dummyWork, &workSize, &dummyIWork, &iWorkSize, &info );

    // Compute the eigenvalues in parallel
    workSize = dummyWork.real();
    iWorkSize = dummyIWork;
    vector<dcomplex> work(workSize);
    vector<int> iWork(iWorkSize);
    EL_SCALAPACK(pzlahqr)
    ( &wantt, &wantz, &n, &ilo, &ihi, H, descH, w, &ilo, &ihi, 0, descQ,
      work.data(), &workSize, iWork.data(), &iWorkSize, &info );
    if( info != 0 )
        RuntimeError("pzlahqr exited with info=",info);
}

void HessenbergSchur
( int n, float* H, const int* descH, scomplex* w, float* Q, const int* descQ, 
  bool fullTriangle, bool multiplyQ, bool aed ) 
{
    EL_DEBUG_CSE
    const int ilo=1, ihi=n;
    vector<float> wr(n), wi(n);
    int info;
    if( aed )
    {
        cerr << 
          "WARNING: PSHSEQR seems to have a bug in its eigenvalue reordering" 
          << endl;
        const char job=(fullTriangle ? 'S' : 'E'),
                   compz=(multiplyQ ? 'V' : 'I');

        // Query the workspace sizes. Due to a bug in p{s,d}hseqr's workspace
        // querying, which is located in p{s,d}laqr1, 
        //    https://github.com/poulson/scalapack/commits/master 
        // we must be a bit more careful.
        int workSize=-1, dummyIWork=3, iWorkSize=-1;
        float dummyWork;
        EL_SCALAPACK(pshseqr)
        ( &job, &compz, &n, &ilo, &ihi, H, descH, wr.data(), wi.data(), 
          Q, descQ, &dummyWork, &workSize, &dummyIWork, &iWorkSize, &info );

        // Compute the eigenvalues in parallel
        workSize = dummyWork;
        iWorkSize = dummyIWork;
        vector<float> work(workSize);
        vector<int> iWork(iWorkSize);
        EL_SCALAPACK(pshseqr)
        ( &job, &compz, &n, &ilo, &ihi, H, descH, wr.data(), wi.data(), 
          Q, descQ, work.data(), &workSize, iWork.data(), &iWorkSize, &info );
        if( info != 0 )
            RuntimeError("pshseqr exited with info=",info);
    }
    else
    {
        if( multiplyQ == false )
            LogicError("Forcing the matrix to identity is not yet supported");
        FortranLogical wantt=(fullTriangle ? FORTRAN_TRUE : FORTRAN_FALSE), 
                       wantz=FORTRAN_TRUE;

        // PSLAHQR does not support a workspace query and instead assumes
        // that the workspace is at least 
        //     3*N + max(2*max(ldH,ldQ)+2*localWidth,numRowBlocks)

        const int context  = descH[1];
        const int mb       = descH[4];
        const int nb       = descH[5];
        const int rowAlign = descH[7];
        const int ldH      = descH[8];
        const int ldQ      = descQ[8];

        const int rowCut = 0;
        const int rowRank = blacs::GridCol( context );
        const int rowStride = blacs::GridWidth( context );

        const Int localWidth =
            BlockedLength(n,rowRank,rowAlign,nb,rowCut,rowStride);
        const Int numRowBlocks = n/mb + 1;
        int workSize= 3*n + Max(2*Max(ldH,ldQ)+2*localWidth,numRowBlocks),
            iWorkSize=0;
        vector<float> work(workSize);
        vector<int> iWork(iWorkSize);

        // Compute the eigenvalues in parallel
        EL_SCALAPACK(pslahqr)
        ( &wantt, &wantz, &n, &ilo, &ihi, H, descH, wr.data(), wi.data(), 
          &ilo, &ihi, Q, descQ, work.data(), &workSize, 
          iWork.data(), &iWorkSize, &info );
        if( info != 0 )
            RuntimeError("pslahqr exited with info=",info);
    }
    // Combine the real and imaginary components of the eigenvalues
    for( int j=0; j<n; ++j )
        w[j] = std::complex<float>(wr[j],wi[j]);
}

void HessenbergSchur
( int n,
  double* H, const int* descH,
  dcomplex* w, 
  double* Q, const int* descQ,
  bool fullTriangle,
  bool multiplyQ,
  bool aed ) 
{
    EL_DEBUG_CSE
    const int ilo=1, ihi=n;
    vector<double> wr(n), wi(n);
    int info;
    if( aed )
    {
        cerr << 
          "WARNING: PDHSEQR seems to have a bug in its eigenvalue reordering" 
          << endl;
        const char job=(fullTriangle ? 'S' : 'E'),
                   compz=(multiplyQ ? 'V' : 'I');

        // Query the workspace sizes. Due to a bug in p{s,d}hseqr's workspace
        // querying, which is located in p{s,d}laqr1, 
        //    https://github.com/poulson/scalapack/commits/master 
        // we must be a bit more careful.
        int workSize=-1, dummyIWork=3, iWorkSize=-1;
        double dummyWork;
        EL_SCALAPACK(pdhseqr)
        ( &job, &compz, &n, &ilo, &ihi, H, descH, wr.data(), wi.data(), 
          Q, descQ, &dummyWork, &workSize, &dummyIWork, &iWorkSize, &info );

        // Compute the eigenvalues in parallel
        workSize = dummyWork;
        iWorkSize = dummyIWork;
        vector<double> work(workSize);
        vector<int> iWork(iWorkSize);
        EL_SCALAPACK(pdhseqr)
        ( &job, &compz, &n, &ilo, &ihi, H, descH, wr.data(), wi.data(), 
          Q, descQ, work.data(), &workSize, iWork.data(), &iWorkSize, &info );
        if( info != 0 )
            RuntimeError("pdhseqr exited with info=",info);
    }
    else
    {
        if( multiplyQ == false )
            LogicError("Forcing the matrix to identity is not yet supported");
        FortranLogical wantt=(fullTriangle ? FORTRAN_TRUE : FORTRAN_FALSE), 
                       wantz=FORTRAN_TRUE;

        // PDLAHQR does not support a workspace query and instead assumes
        // that the workspace is at least 
        //     3*N + max(2*max(ldH,ldQ)+2*localWidth,numRowBlocks)

        const int context  = descH[1];
        const int mb       = descH[4];
        const int nb       = descH[5];
        const int rowAlign = descH[7];
        const int ldH      = descH[8];
        const int ldQ      = descQ[8];

        const int rowCut = 0;
        const int rowRank = blacs::GridCol( context );
        const int rowStride = blacs::GridWidth( context );

        const Int localWidth =
            BlockedLength(n,rowRank,rowAlign,nb,rowCut,rowStride);
        const Int numRowBlocks = n/mb + 1;
        int workSize= 3*n + Max(2*Max(ldH,ldQ)+2*localWidth,numRowBlocks),
            iWorkSize=0;
        vector<double> work(workSize);
        vector<int> iWork(iWorkSize);

        // Compute the eigenvalues in parallel
        EL_SCALAPACK(pdlahqr)
        ( &wantt, &wantz, &n, &ilo, &ihi, H, descH, wr.data(), wi.data(), 
          &ilo, &ihi, Q, descQ, work.data(), &workSize, 
          iWork.data(), &iWorkSize, &info );
        if( info != 0 )
            RuntimeError("pdlahqr exited with info=",info);
    }
    // Combine the real and imaginary components of the eigenvalues
    for( int j=0; j<n; ++j )
        w[j] = std::complex<double>(wr[j],wi[j]);
}

void HessenbergSchur
( int n,
  scomplex* H, const int* descH,
  scomplex* w, 
  scomplex* Q, const int* descQ,
  bool fullTriangle,
  bool multiplyQ,
  bool aed ) 
{
    EL_DEBUG_CSE
    if( !multiplyQ )
        LogicError("Forcing the matrix to identity is not yet supported");
    if( aed )
        LogicError("AED is not supported for complex matrices");
    FortranLogical wantt=(fullTriangle ? FORTRAN_TRUE : FORTRAN_FALSE), 
                   wantz=FORTRAN_TRUE;
    const int ilo=1, ihi=n;

    // Query the workspace sizes
    // NOTE: dummyIWork will currently be left unmodified (hence the
    //       zero initialization)!
    int workSize=-1, dummyIWork=0, iWorkSize=-1, info;
    scomplex dummyWork;
    EL_SCALAPACK(pclahqr)
    ( &wantt, &wantz, &n, &ilo, &ihi, H, descH, w, &ilo, &ihi, Q, descQ,
      &dummyWork, &workSize, &dummyIWork, &iWorkSize, &info );

    // Compute the eigenvalues in parallel
    workSize = dummyWork.real();
    iWorkSize = dummyIWork;
    vector<scomplex> work(workSize);
    vector<int> iWork(iWorkSize);
    EL_SCALAPACK(pclahqr)
    ( &wantt, &wantz, &n, &ilo, &ihi, H, descH, w, &ilo, &ihi, Q, descQ,
      work.data(), &workSize, iWork.data(), &iWorkSize, &info );
    if( info != 0 )
        RuntimeError("pclahqr exited with info=",info);
}

void HessenbergSchur
( int n,
  dcomplex* H, const int* descH,
  dcomplex* w, 
  dcomplex* Q, const int* descQ,
  bool fullTriangle,
  bool multiplyQ,
  bool aed ) 
{
    EL_DEBUG_CSE
    if( !multiplyQ )
        LogicError("Forcing the matrix to identity is not yet supported");
    if( aed )
        LogicError("AED is not supported for complex matrices");
    FortranLogical wantt=(fullTriangle ? FORTRAN_TRUE : FORTRAN_FALSE), 
                   wantz=FORTRAN_TRUE;
    const int ilo=1, ihi=n;

    // Query the workspace sizes
    // NOTE: dummyIWork will currently be left unmodified (hence the
    //       zero initialization)!
    int workSize=-1, dummyIWork=0, iWorkSize=-1, info;
    dcomplex dummyWork;
    EL_SCALAPACK(pzlahqr)
    ( &wantt, &wantz, &n, &ilo, &ihi, H, descH, w, &ilo, &ihi, Q, descQ,
      &dummyWork, &workSize, &dummyIWork, &iWorkSize, &info );

    // Compute the eigenvalues in parallel
    workSize = dummyWork.real();
    iWorkSize = dummyIWork;
    vector<dcomplex> work(workSize);
    vector<int> iWork(iWorkSize);
    EL_SCALAPACK(pzlahqr)
    ( &wantt, &wantz, &n, &ilo, &ihi, H, descH, w, &ilo, &ihi, Q, descQ,
      work.data(), &workSize, iWork.data(), &iWorkSize, &info );
    if( info != 0 )
        RuntimeError("pzlahqr exited with info=",info);
}

// Hessenberg eigenvalues/pairs via the QR algorithm
// -------------------------------------------------
void HessenbergEig( int n, float* H, const int* descH, scomplex* w ) 
{
    EL_DEBUG_CSE
    HessenbergSchur( n, H, descH, w, false );
}

void HessenbergEig( int n, double* H, const int* descH, dcomplex* w ) 
{
    EL_DEBUG_CSE
    HessenbergSchur( n, H, descH, w, false );
}

void HessenbergEig( int n, scomplex* H, const int* descH, scomplex* w ) 
{
    EL_DEBUG_CSE
    HessenbergSchur( n, H, descH, w, false );
}

void HessenbergEig( int n, dcomplex* H, const int* descH, dcomplex* w ) 
{
    EL_DEBUG_CSE
    HessenbergSchur( n, H, descH, w, false );
}

} // namespace scalapack
} // namespace El

#endif // ifdef EL_HAVE_SCALAPACK
