/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_LATTICE_C_H
#define EL_LATTICE_C_H

#ifdef __cplusplus
extern "C" {
#endif

EL_EXPORT ElError ElLLL_s
( ElMatrix_s B, ElMatrix_s QR, float delta, float innerTol,
  bool presort, bool smallestFirst, bool progress, ElInt* numBacktrack );
EL_EXPORT ElError ElLLL_d
( ElMatrix_d B, ElMatrix_d QR, double delta, double innerTol,
  bool presort, bool smallestFirst, bool progress, ElInt* numBacktrack );
EL_EXPORT ElError ElLLL_c
( ElMatrix_c B, ElMatrix_c QR, float delta, float innerTol,
  bool presort, bool smallestFirst, bool progress, ElInt* numBacktrack );
EL_EXPORT ElError ElLLL_z
( ElMatrix_z B, ElMatrix_z QR, double delta, double innerTol,
  bool presort, bool smallestFirst, bool progress, ElInt* numBacktrack );

EL_EXPORT ElError ElLLLDelta_s( ElConstMatrix_s QR, float* delta );
EL_EXPORT ElError ElLLLDelta_d( ElConstMatrix_d QR, double* delta );
EL_EXPORT ElError ElLLLDelta_c( ElConstMatrix_c QR, float* delta );
EL_EXPORT ElError ElLLLDelta_z( ElConstMatrix_z QR, double* delta );

EL_EXPORT ElError ElLatticeGramSchmidt_s
( ElConstMatrix_s B, ElMatrix_s G, ElMatrix_s M );
EL_EXPORT ElError ElLatticeGramSchmidt_d
( ElConstMatrix_d B, ElMatrix_d G, ElMatrix_d M );
EL_EXPORT ElError ElLatticeGramSchmidt_c
( ElConstMatrix_c B, ElMatrix_c G, ElMatrix_c M );
EL_EXPORT ElError ElLatticeGramSchmidt_z
( ElConstMatrix_z B, ElMatrix_z G, ElMatrix_z M );

#ifdef __cplusplus
}
#endif

#endif // ifndef EL_LATTICE_C_H
