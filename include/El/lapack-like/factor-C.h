/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_LAPACK_FACTOR_C_H
#define EL_LAPACK_FACTOR_C_H

#include "El/core/DistMatrix-C.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Cholesky
   ======== */
ElError ElCholesky_s( ElUpperOrLower uplo, ElMatrix_s A );
ElError ElCholesky_d( ElUpperOrLower uplo, ElMatrix_d A );
ElError ElCholesky_c( ElUpperOrLower uplo, ElMatrix_c A );
ElError ElCholesky_z( ElUpperOrLower uplo, ElMatrix_z A );

/* NOTE: 'A' must be in a [MC,MR] distribution */
ElError ElCholeskyDist_s( ElUpperOrLower uplo, ElDistMatrix_s A );
ElError ElCholeskyDist_d( ElUpperOrLower uplo, ElDistMatrix_d A );
ElError ElCholeskyDist_c( ElUpperOrLower uplo, ElDistMatrix_c A );
ElError ElCholeskyDist_z( ElUpperOrLower uplo, ElDistMatrix_z A );

ElError ElReverseCholesky_s( ElUpperOrLower uplo, ElMatrix_s A );
ElError ElReverseCholesky_d( ElUpperOrLower uplo, ElMatrix_d A );
ElError ElReverseCholesky_c( ElUpperOrLower uplo, ElMatrix_c A );
ElError ElReverseCholesky_z( ElUpperOrLower uplo, ElMatrix_z A );

/* NOTE: 'A' must be in a [MC,MR] distribution */
ElError ElReverseCholeskyDist_s( ElUpperOrLower uplo, ElDistMatrix_s A );
ElError ElReverseCholeskyDist_d( ElUpperOrLower uplo, ElDistMatrix_d A );
ElError ElReverseCholeskyDist_c( ElUpperOrLower uplo, ElDistMatrix_c A );
ElError ElReverseCholeskyDist_z( ElUpperOrLower uplo, ElDistMatrix_z A );

ElError ElCholeskyPiv_s( ElUpperOrLower uplo, ElMatrix_s A, ElMatrix_i p );
ElError ElCholeskyPiv_d( ElUpperOrLower uplo, ElMatrix_d A, ElMatrix_i p );
ElError ElCholeskyPiv_c( ElUpperOrLower uplo, ElMatrix_c A, ElMatrix_i p );
ElError ElCholeskyPiv_z( ElUpperOrLower uplo, ElMatrix_z A, ElMatrix_i p );

/* NOTE: 'A' must be in a [MC,MR] distribution, and
         'p' must be in a [VC,STAR] distribution */
ElError ElCholeskyPivDist_s
( ElUpperOrLower uplo, ElDistMatrix_s A, ElDistMatrix_i p );
ElError ElCholeskyPivDist_d
( ElUpperOrLower uplo, ElDistMatrix_d A, ElDistMatrix_i p );
ElError ElCholeskyPivDist_c
( ElUpperOrLower uplo, ElDistMatrix_c A, ElDistMatrix_i p );
ElError ElCholeskyPivDist_z
( ElUpperOrLower uplo, ElDistMatrix_z A, ElDistMatrix_i p );

ElError ElCholeskyMod_s
( ElUpperOrLower uplo, ElMatrix_s T, float alpha, ElMatrix_s V );
ElError ElCholeskyMod_d
( ElUpperOrLower uplo, ElMatrix_d T, double alpha, ElMatrix_d V );
ElError ElCholeskyMod_c
( ElUpperOrLower uplo, ElMatrix_c T, float alpha, ElMatrix_c V );
ElError ElCholeskyMod_z
( ElUpperOrLower uplo, ElMatrix_z T, double alpha, ElMatrix_z V );

/* NOTE: 'T' and 'V' must be in a [MC,MR] distribution */
ElError ElCholeskyModDist_s
( ElUpperOrLower uplo, ElDistMatrix_s T, float alpha, ElDistMatrix_s V );
ElError ElCholeskyModDist_d
( ElUpperOrLower uplo, ElDistMatrix_d T, double alpha, ElDistMatrix_d V );
ElError ElCholeskyModDist_c
( ElUpperOrLower uplo, ElDistMatrix_c T, float alpha, ElDistMatrix_c V );
ElError ElCholeskyModDist_z
( ElUpperOrLower uplo, ElDistMatrix_z T, double alpha, ElDistMatrix_z V );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_LAPACK_FACTOR_C_H */
