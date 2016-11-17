/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   Copyright (c) 2016, Ron Estrin
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_LATTICE_C_H
#define EL_LATTICE_C_H

#ifdef __cplusplus
extern "C" {
#endif

/* Lenstra-Lenstra-Lovasz lattice reduction
   ======================================== */

typedef struct
{
    float delta;
    float eta;
    ElInt rank;
    ElInt nullity;
    ElInt numSwaps;
    ElInt firstSwap;
    float logVol;
} ElLLLInfo_s;

typedef struct
{
    double delta;
    double eta;
    ElInt rank;
    ElInt nullity;
    ElInt numSwaps;
    ElInt firstSwap;
    double logVol;
} ElLLLInfo_d;

typedef enum {
  EL_LLL_WEAK,
  EL_LLL_NORMAL,
  EL_LLL_DEEP,
  EL_LLL_DEEP_REDUCE
} ElLLLVariant;

typedef struct
{
    float delta;
    float eta;
    ElLLLVariant variant;
    bool recursive;
    ElInt cutoff;
    float precisionFudge;
    ElInt minColThresh;
    bool unsafeSizeReduct;
    bool presort;
    bool smallestFirst;
    float reorthogTol;
    ElInt numOrthog;
    float zeroTol;
    float blockingThresh;
    bool progress;
    bool time;
    bool jumpstart;
    ElInt startCol;
} ElLLLCtrl_s;
EL_EXPORT ElError ElLLLCtrlDefault_s( ElLLLCtrl_s* ctrl );

typedef struct
{
    double delta;
    double eta;
    ElLLLVariant variant;
    bool recursive;
    ElInt cutoff;
    double precisionFudge;
    ElInt minColThresh;
    bool unsafeSizeReduct;
    bool presort;
    bool smallestFirst;
    double reorthogTol;
    ElInt numOrthog;
    double zeroTol;
    float blockingThresh;
    bool progress;
    bool time;
    bool jumpstart;
    ElInt startCol;
} ElLLLCtrl_d;
EL_EXPORT ElError ElLLLCtrlDefault_d( ElLLLCtrl_d* ctrl );

EL_EXPORT ElError ElLLL_s( ElMatrix_s B, ElLLLCtrl_s ctrl, ElLLLInfo_s* info );
EL_EXPORT ElError ElLLL_d( ElMatrix_d B, ElLLLCtrl_d ctrl, ElLLLInfo_d* info );
EL_EXPORT ElError ElLLL_c( ElMatrix_c B, ElLLLCtrl_s ctrl, ElLLLInfo_s* info );
EL_EXPORT ElError ElLLL_z( ElMatrix_z B, ElLLLCtrl_d ctrl, ElLLLInfo_d* info );

EL_EXPORT ElError ElLLLFormR_s
( ElMatrix_s B, ElMatrix_s R, ElLLLCtrl_s ctrl, ElLLLInfo_s* info );
EL_EXPORT ElError ElLLLFormR_d
( ElMatrix_d B, ElMatrix_d R, ElLLLCtrl_d ctrl, ElLLLInfo_d* info );
EL_EXPORT ElError ElLLLFormR_c
( ElMatrix_c B, ElMatrix_c R, ElLLLCtrl_s ctrl, ElLLLInfo_s* info );
EL_EXPORT ElError ElLLLFormR_z
( ElMatrix_z B, ElMatrix_z R, ElLLLCtrl_d ctrl, ElLLLInfo_d* info );

EL_EXPORT ElError ElLLLFull_s
( ElMatrix_s B, ElMatrix_s U, ElMatrix_s R,
  ElLLLCtrl_s ctrl, ElLLLInfo_s* info );
EL_EXPORT ElError ElLLLFull_d
( ElMatrix_d B, ElMatrix_d U, ElMatrix_d R,
  ElLLLCtrl_d ctrl, ElLLLInfo_d* info );
EL_EXPORT ElError ElLLLFull_c
( ElMatrix_c B, ElMatrix_c U, ElMatrix_c R,
  ElLLLCtrl_s ctrl, ElLLLInfo_s* info );
EL_EXPORT ElError ElLLLFull_z
( ElMatrix_z B, ElMatrix_z U, ElMatrix_z R,
  ElLLLCtrl_d ctrl, ElLLLInfo_d* info );

/* Lattice image and kernel
   ======================== */
EL_EXPORT ElError ElLatticeImageAndKernel_s
( ElConstMatrix_s B, ElMatrix_s M, ElMatrix_s K, ElLLLCtrl_s ctrl );
EL_EXPORT ElError ElLatticeImageAndKernel_d
( ElConstMatrix_d B, ElMatrix_d M, ElMatrix_d K, ElLLLCtrl_d ctrl );
EL_EXPORT ElError ElLatticeImageAndKernel_c
( ElConstMatrix_c B, ElMatrix_c M, ElMatrix_c K, ElLLLCtrl_s ctrl );
EL_EXPORT ElError ElLatticeImageAndKernel_z
( ElConstMatrix_z B, ElMatrix_z M, ElMatrix_z K, ElLLLCtrl_d ctrl );

EL_EXPORT ElError ElLatticeImage_s
( ElConstMatrix_s B, ElMatrix_s M, ElLLLCtrl_s ctrl );
EL_EXPORT ElError ElLatticeImage_d
( ElConstMatrix_d B, ElMatrix_d M, ElLLLCtrl_d ctrl );
EL_EXPORT ElError ElLatticeImage_c
( ElConstMatrix_c B, ElMatrix_c M, ElLLLCtrl_s ctrl );
EL_EXPORT ElError ElLatticeImage_z
( ElConstMatrix_z B, ElMatrix_z M, ElLLLCtrl_d ctrl );

EL_EXPORT ElError ElLatticeKernel_s
( ElConstMatrix_s B, ElMatrix_s K, ElLLLCtrl_s ctrl );
EL_EXPORT ElError ElLatticeKernel_d
( ElConstMatrix_d B, ElMatrix_d K, ElLLLCtrl_d ctrl );
EL_EXPORT ElError ElLatticeKernel_c
( ElConstMatrix_c B, ElMatrix_c K, ElLLLCtrl_s ctrl );
EL_EXPORT ElError ElLatticeKernel_z
( ElConstMatrix_z B, ElMatrix_z K, ElLLLCtrl_d ctrl );

/* Search for Z-dependence
   ======================= */
EL_EXPORT ElError ElZDependenceSearch_s
( ElConstMatrix_s z, float NSqrt,
  ElMatrix_s B, ElMatrix_s U, ElLLLCtrl_s ctrl,
  ElInt* numExact );
EL_EXPORT ElError ElZDependenceSearch_d
( ElConstMatrix_d z, double NSqrt,
  ElMatrix_d B, ElMatrix_d U, ElLLLCtrl_d ctrl,
  ElInt* numExact );
EL_EXPORT ElError ElZDependenceSearch_c
( ElConstMatrix_c z, float NSqrt,
  ElMatrix_c B, ElMatrix_c U, ElLLLCtrl_s ctrl,
  ElInt* numExact );
EL_EXPORT ElError ElZDependenceSearch_z
( ElConstMatrix_z z, double NSqrt,
  ElMatrix_z B, ElMatrix_z U, ElLLLCtrl_d ctrl,
  ElInt* numExact );

/* Search for an algebraic relation
   ================================ */
EL_EXPORT ElError ElAlgebraicRelationSearch_s
( float alpha, ElInt n, float NSqrt,
  ElMatrix_s B, ElMatrix_s U, ElLLLCtrl_s ctrl,
  ElInt* numExact );
EL_EXPORT ElError ElAlgebraicRelationSearch_d
( double alpha, ElInt n, double NSqrt,
  ElMatrix_d B, ElMatrix_d U, ElLLLCtrl_d ctrl,
  ElInt* numExact );
EL_EXPORT ElError ElAlgebraicRelationSearch_c
( complex_float alpha, ElInt n, float NSqrt,
  ElMatrix_c B, ElMatrix_c U, ElLLLCtrl_s ctrl,
  ElInt* numExact );
EL_EXPORT ElError ElAlgebraicRelationSearch_z
( complex_double alpha, ElInt n, double NSqrt,
  ElMatrix_z B, ElMatrix_z U, ElLLLCtrl_d ctrl,
  ElInt* numExact );

#ifdef __cplusplus
}
#endif

#endif // ifndef EL_LATTICE_C_H
