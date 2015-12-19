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

typedef struct
{
    ElInt nullity;    
    ElInt numSwaps;
} ElLLLInfo;

typedef struct
{
    float delta;
    bool weak;
    bool presort;
    bool smallestFirst;
    float reorthogTol;
    float zeroTol;
    bool progress;
    bool time;
} ElLLLCtrl_s;
EL_EXPORT ElError ElLLLCtrlDefault_s( ElLLLCtrl_s* ctrl );

typedef struct
{
    double delta;
    bool weak;
    bool presort;
    bool smallestFirst;
    double reorthogTol;
    double zeroTol;
    bool progress;
    bool time;
} ElLLLCtrl_d;
EL_EXPORT ElError ElLLLCtrlDefault_d( ElLLLCtrl_d* ctrl );

EL_EXPORT ElError ElLLL_s
( ElMatrix_s B, ElMatrix_s QR, ElLLLCtrl_s ctrl, ElLLLInfo* info );
EL_EXPORT ElError ElLLL_d
( ElMatrix_d B, ElMatrix_d QR, ElLLLCtrl_d ctrl, ElLLLInfo* info );
EL_EXPORT ElError ElLLL_c
( ElMatrix_c B, ElMatrix_c QR, ElLLLCtrl_s ctrl, ElLLLInfo* info );
EL_EXPORT ElError ElLLL_z
( ElMatrix_z B, ElMatrix_z QR, ElLLLCtrl_d ctrl, ElLLLInfo* info );

EL_EXPORT ElError ElLLLFull_s
( ElMatrix_s B, ElMatrix_s U, ElMatrix_s UInv, ElMatrix_s QR,
  ElLLLCtrl_s ctrl, ElLLLInfo* info );
EL_EXPORT ElError ElLLLFull_d
( ElMatrix_d B, ElMatrix_d U, ElMatrix_d UInv, ElMatrix_d QR,
  ElLLLCtrl_d ctrl, ElLLLInfo* info );
EL_EXPORT ElError ElLLLFull_c
( ElMatrix_c B, ElMatrix_c U, ElMatrix_c UInv, ElMatrix_c QR,
  ElLLLCtrl_s ctrl, ElLLLInfo* info );
EL_EXPORT ElError ElLLLFull_z
( ElMatrix_z B, ElMatrix_z U, ElMatrix_z UInv, ElMatrix_z QR,
  ElLLLCtrl_d ctrl, ElLLLInfo* info );

EL_EXPORT ElError ElLLLDelta_s
( ElConstMatrix_s QR, ElLLLCtrl_s ctrl, float* delta );
EL_EXPORT ElError ElLLLDelta_d
( ElConstMatrix_d QR, ElLLLCtrl_d ctrl, double* delta );
EL_EXPORT ElError ElLLLDelta_c
( ElConstMatrix_c QR, ElLLLCtrl_s ctrl, float* delta );
EL_EXPORT ElError ElLLLDelta_z
( ElConstMatrix_z QR, ElLLLCtrl_d ctrl, double* delta );

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
