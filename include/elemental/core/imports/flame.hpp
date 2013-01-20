/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef CORE_FLAME_HPP
#define CORE_FLAME_HPP

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

#endif // ifndef CORE_FLAME_HPP
