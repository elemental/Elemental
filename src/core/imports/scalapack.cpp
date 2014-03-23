/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "elemental-lite.hpp"

using elem::scomplex;
using elem::dcomplex;

extern "C" {

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
  const int* ilo, const int* ihi, scomplex* A, const int* desca,
  scomplex* w, const int* iloz, const int* ihiz, scomplex* z, const int* descz,
  scomplex* work, const int* lwork, int* iwork, const int* ilwork, int* info );
void ELEM_SCALAPACK(pzlahqr)
( const ELEM_FORT_LOGICAL* wantt, const ELEM_FORT_LOGICAL* wantz, const int* n,
  const int* ilo, const int* ihi, dcomplex* A, const int* desca,
  dcomplex* w, const int* iloz, const int* ihiz, dcomplex* z, const int* descz,
  dcomplex* work, const int* lwork, int* iwork, const int* ilwork, int* info );

// Pipelined QR algorithm with AED for big matrices
// ------------------------------------------------
void ELEM_SCALAPACK(pslaqr0)
( const ELEM_FORT_LOGICAL* wantt, const ELEM_FORT_LOGICAL* wantz, const int* n,
  const int* ilo, const int* ihi, float* H, const int* desch, 
  float* wr, float* wi, const int* iloz, const int* ihiz, const int* descz, 
  float* work, const int* lwork, int* iwork, const int* liwork, int* info,
  const int* reclevel );
void ELEM_SCALAPACK(pdlaqr0)
( const ELEM_FORT_LOGICAL* wantt, const ELEM_FORT_LOGICAL* wantz, const int* n,
  const int* ilo, const int* ihi, double* H, const int* desch, 
  double* wr, double* wi, const int* iloz, const int* ihiz, const int* descz, 
  double* work, const int* lwork, int* iwork, const int* liwork, int* info,
  const int* reclevel );

// Pipelined QR algorithm with AED for small matrices
// --------------------------------------------------
void ELEM_SCALAPACK(pslaqr1)
( const ELEM_FORT_LOGICAL* wantt, const ELEM_FORT_LOGICAL* wantz, const int* n,
  const int* ilo, const int* ihi, float* H, const int* desch, 
  float* wr, float* wi, const int* iloz, const int* ihiz, const int* descz, 
  float* work, const int* lwork, int* iwork, const int* liwork, int* info );
void ELEM_SCALAPACK(pdlaqr1)
( const ELEM_FORT_LOGICAL* wantt, const ELEM_FORT_LOGICAL* wantz, const int* n,
  const int* ilo, const int* ihi, double* H, const int* desch, 
  double* wr, double* wi, const int* iloz, const int* ihiz, const int* descz, 
  double* work, const int* lwork, int* iwork, const int* liwork, int* info );

} // extern "C"

namespace elem {
namespace scalapack {

void HessenbergSchur( int n, float* H, const int* desch, scomplex* w ) 
{
    DEBUG_ONLY(CallStackEntry cse("scalapack::HessenbergSchur"))
    LogicError("This routine not yet written");
}

void HessenbergSchur( int n, double* H, const int* desch, dcomplex* w ) 
{
    DEBUG_ONLY(CallStackEntry cse("scalapack::HessenbergSchur"))
    LogicError("This routine not yet written");
}

void HessenbergSchur( int n, scomplex* H, const int* desch, scomplex* w ) 
{
    DEBUG_ONLY(CallStackEntry cse("scalapack::HessenbergSchur"))
    LogicError("This routine not yet written");
}

void HessenbergSchur( int n, dcomplex* H, const int* desch, dcomplex* w ) 
{
    DEBUG_ONLY(CallStackEntry cse("scalapack::HessenbergSchur"))
    LogicError("This routine not yet written");
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
