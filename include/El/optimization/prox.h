/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_OPTIMIZATION_PROX_C_H
#define EL_OPTIMIZATION_PROX_C_H

#ifdef __cplusplus
extern "C" {
#endif

/* Clipping
   -------- */
EL_EXPORT ElError ElLowerClip_s( ElMatrix_s X, float lowerBound );
EL_EXPORT ElError ElLowerClip_d( ElMatrix_d X, double lowerBound );

EL_EXPORT ElError ElLowerClipDist_s( ElDistMatrix_s X, float lowerBound );
EL_EXPORT ElError ElLowerClipDist_d( ElDistMatrix_d X, double lowerBound );

EL_EXPORT ElError ElUpperClip_s( ElMatrix_s X, float upperBound );
EL_EXPORT ElError ElUpperClip_d( ElMatrix_d X, double upperBound );

EL_EXPORT ElError ElUpperClipDist_s( ElDistMatrix_s X, float upperBound );
EL_EXPORT ElError ElUpperClipDist_d( ElDistMatrix_d X, double upperBound );

EL_EXPORT ElError ElClip_s
( ElMatrix_s X, float lowerBound, float upperBound );
EL_EXPORT ElError ElClip_d
( ElMatrix_d X, double lowerBound, double upperBound );

EL_EXPORT ElError ElClipDist_s
( ElDistMatrix_s X, float lowerBound, float upperBound );
EL_EXPORT ElError ElClipDist_d
( ElDistMatrix_d X, double lowerBound, double upperBound );

/* Frobenius-norm proximal map
   --------------------------- */
/* arg min || A ||_F + rho/2 || A - A0 ||_F^2 
      A                                       */
EL_EXPORT ElError ElFrobeniusProx_s( ElMatrix_s A, float rho );
EL_EXPORT ElError ElFrobeniusProx_d( ElMatrix_d A, double rho );
EL_EXPORT ElError ElFrobeniusProx_c( ElMatrix_c A, float rho );
EL_EXPORT ElError ElFrobeniusProx_z( ElMatrix_z A, double rho );

EL_EXPORT ElError ElFrobeniusProxDist_s( ElDistMatrix_s A, float rho );
EL_EXPORT ElError ElFrobeniusProxDist_d( ElDistMatrix_d A, double rho );
EL_EXPORT ElError ElFrobeniusProxDist_c( ElDistMatrix_c A, float rho );
EL_EXPORT ElError ElFrobeniusProxDist_z( ElDistMatrix_z A, double rho );

/* Hinge-loss proximal map
   ----------------------- */
EL_EXPORT ElError ElHingeLossProx_s( ElMatrix_s A, float rho );
EL_EXPORT ElError ElHingeLossProx_d( ElMatrix_d A, double rho );

EL_EXPORT ElError ElHingeLossProxDist_s( ElDistMatrix_s A, float rho );
EL_EXPORT ElError ElHingeLossProxDist_d( ElDistMatrix_d A, double rho );

/* Logistic proximal map
   --------------------- */
EL_EXPORT ElError ElLogisticProx_s( ElMatrix_s A, float rho );
EL_EXPORT ElError ElLogisticProx_d( ElMatrix_d A, double rho );

EL_EXPORT ElError ElLogisticProxDist_s( ElDistMatrix_s A, float rho );
EL_EXPORT ElError ElLogisticProxDist_d( ElDistMatrix_d A, double rho );

/* Singular-value soft thresholding
   -------------------------------- */
EL_EXPORT ElError ElSVT_s( ElMatrix_s A, float rho, bool relative );
EL_EXPORT ElError ElSVT_d( ElMatrix_d A, double rho, bool relative );
EL_EXPORT ElError ElSVT_c( ElMatrix_c A, float rho, bool relative );
EL_EXPORT ElError ElSVT_z( ElMatrix_z A, double rho, bool relative );

EL_EXPORT ElError ElSVTDist_s( ElDistMatrix_s A, float rho, bool relative );
EL_EXPORT ElError ElSVTDist_d( ElDistMatrix_d A, double rho, bool relative );
EL_EXPORT ElError ElSVTDist_c( ElDistMatrix_c A, float rho, bool relative );
EL_EXPORT ElError ElSVTDist_z( ElDistMatrix_z A, double rho, bool relative );

/* TODO: Expert versions */

/* Soft-thresholding
   ----------------- */
EL_EXPORT ElError ElSoftThreshold_s( ElMatrix_s A, float rho, bool relative );
EL_EXPORT ElError ElSoftThreshold_d( ElMatrix_d A, double rho, bool relative );
EL_EXPORT ElError ElSoftThreshold_c( ElMatrix_c A, float rho, bool relative );
EL_EXPORT ElError ElSoftThreshold_z( ElMatrix_z A, double rho, bool relative );

EL_EXPORT ElError ElSoftThresholdDist_s
( ElDistMatrix_s A, float rho, bool relative );
EL_EXPORT ElError ElSoftThresholdDist_d
( ElDistMatrix_d A, double rho, bool relative );
EL_EXPORT ElError ElSoftThresholdDist_c
( ElDistMatrix_c A, float rho, bool relative );
EL_EXPORT ElError ElSoftThresholdDist_z
( ElDistMatrix_z A, double rho, bool relative );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_OPTIMIZATION_PROX_C_H */
