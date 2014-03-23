/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "elemental-lite.hpp"
#ifdef ELEM_HAVE_FLA_BSVD

extern "C" {

typedef int FLA_Error;
FLA_Error FLA_Bsvd_v_opd_var1
( int       k,
  int       mU,
  int       mV,
  int       nGH,
  int       nIterMax,
  double*   d, int dInc,
  double*   e, int eInc,
  elem::dcomplex* G, int rsG, int csG,
  elem::dcomplex* H, int rsH, int csH,
  double*   U, int rsU, int csU,
  double*   V, int rsV, int csV,
  int       nb );

FLA_Error FLA_Bsvd_v_opz_var1
( int       k,
  int       mU,
  int       mV,
  int       nGH,
  int       nIterMax,
  double*   d, int dInc,
  double*   e, int eInc,
  elem::dcomplex* G, int rsG, int csG,
  elem::dcomplex* H, int rsH, int csH,
  elem::dcomplex* U, int rsU, int csU,
  elem::dcomplex* V, int rsV, int csV,
  int       nb );

} // extern "C"

namespace elem {

void FlaBidiagSVD
( int k, int mU, int mV, double* d, double* e, 
  double* U, int ldu, double* V, int ldv, 
  int numAccum, int maxNumIts, int bAlg )
{
    std::vector<Complex<double>> G( (k-1)*numAccum ), H( (k-1)*numAccum ); 
    FLA_Bsvd_v_opd_var1
    ( k, mU, mV, numAccum, maxNumIts, d, 1, e, 1, 
      G.data(), 1, k-1, H.data(), 1, k-1, U, 1, ldu, V, 1, ldv, bAlg );
}

void FlaBidiagSVD
( int k, int mU, int mV, double* d, double* e, 
  Complex<double>* U, int ldu, Complex<double>* V, int ldv, 
  int numAccum, int maxNumIts, int bAlg )
{
    std::vector<Complex<double>> G( (k-1)*numAccum ), H( (k-1)*numAccum ); 
    FLA_Bsvd_v_opz_var1
    ( k, mU, mV, numAccum, maxNumIts, d, 1, e, 1, 
      G.data(), 1, k-1, H.data(), 1, k-1, U, 1, ldu, V, 1, ldv, bAlg );
}

} // namespace elem

#endif // ifdef ELEM_HAVE_FLA_BSVD
