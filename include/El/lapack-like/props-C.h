/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_LAPACK_PROPS_C_H
#define EL_LAPACK_PROPS_C_H

#include "El/core/DistMatrix-C.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
  EL_PS_TWO_NORM,
  EL_PS_ONE_NORM
} ElPseudospecNorm;

typedef struct {
  ElInt realSize, imagSize;

  ElInt imgSaveFreq, numSaveFreq, imgDispFreq;
  ElInt imgSaveCount, numSaveCount, imgDispCount;
  const char *imgBase, *numBase;
  ElFileFormat imgFormat, numFormat; 
  bool itCounts;
} ElSnapshotCtrl;
ElError ElSnapshotCtrlDefault( ElSnapshotCtrl* ctrl );
/* NOTE: Since conversion from SnapshotCtrl involves deep copies of char* */
ElError ElSnapshotCtrlDestroy( const ElSnapshotCtrl* ctrl );

typedef struct {
  ElPseudospecNorm norm;
  ElInt blockWidth;

  bool schur;
  bool forceComplexSchur;
  bool forceComplexPs;
  ElSchurCtrl_s schurCtrl;

  ElInt maxIts;
  float tol;
  bool deflate;
  
  bool arnoldi;
  ElInt basisSize;
  bool reorthog;

  bool progress;
} ElPseudospecCtrl_s;
ElError ElPseudospecCtrlDefault_s( ElPseudospecCtrl_s* ctrl );
/* NOTE: Since conversion from SnapshotCtrl involves deep copies of char* */
ElError ElPseudospecCtrlDestroy_s( const ElPseudospecCtrl_s* ctrl );

typedef struct {
  ElPseudospecNorm norm;
  ElInt blockWidth;

  bool schur;
  bool forceComplexSchur;
  bool forceComplexPs;
  ElSchurCtrl_d schurCtrl;

  ElInt maxIts;
  double tol;
  bool deflate;
  
  bool arnoldi;
  ElInt basisSize;
  bool reorthog;

  bool progress;
} ElPseudospecCtrl_d;
ElError ElPseudospecCtrlDefault_d( ElPseudospecCtrl_d* ctrl );
/* NOTE: Since conversion from SnapshotCtrl involves deep copies of char* */
ElError ElPseudospecCtrlDestroy_d( const ElPseudospecCtrl_d* ctrl );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_LAPACK_PROPS_C_H */
