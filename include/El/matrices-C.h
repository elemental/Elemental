/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_MATRICES_C_H
#define EL_MATRICES_C_H

#ifdef __cplusplus
extern "C" {
#endif

/* Bull's head
   =========== */
ElError ElBullsHead_c( ElMatrix_c A, ElInt n );
ElError ElBullsHead_z( ElMatrix_z A, ElInt n );

ElError ElBullsHeadDist_c( ElDistMatrix_c A, ElInt n );
ElError ElBullsHeadDist_z( ElDistMatrix_z A, ElInt n );

/* Cauchy
   ====== */
ElError ElCauchy_s
( ElMatrix_s A, ElInt xSize, float* xBuf,
                ElInt ySize, float* yBuf ); 
ElError ElCauchy_d
( ElMatrix_d A, ElInt xSize, double* xBuf,
                ElInt ySize, double* yBuf ); 
ElError ElCauchy_c
( ElMatrix_c A, ElInt xSize, complex_float* xBuf,
                ElInt ySize, complex_float* yBuf ); 
ElError ElCauchy_z
( ElMatrix_z A, ElInt xSize, complex_double* xBuf,
                ElInt ySize, complex_double* yBuf ); 

ElError ElCauchyDist_s
( ElDistMatrix_s A, ElInt xSize, float* xBuf,
                    ElInt ySize, float* yBuf ); 
ElError ElCauchyDist_d
( ElDistMatrix_d A, ElInt xSize, double* xBuf,
                    ElInt ySize, double* yBuf ); 
ElError ElCauchyDist_c
( ElDistMatrix_c A, ElInt xSize, complex_float* xBuf,
                    ElInt ySize, complex_float* yBuf ); 
ElError ElCauchyDist_z
( ElDistMatrix_z A, ElInt xSize, complex_double* xBuf,
                    ElInt ySize, complex_double* yBuf ); 

/* Cauchy-like
   =========== */
ElError ElCauchyLike_s
( ElMatrix_s A, ElInt rSize, float* rBuf,
                ElInt sSize, float* sBuf,
                ElInt xSize, float* xBuf,
                ElInt ySize, float* yBuf ); 
ElError ElCauchyLike_d
( ElMatrix_d A, ElInt rSize, double* rBuf,
                ElInt sSize, double* sBuf,
                ElInt xSize, double* xBuf,
                ElInt ySize, double* yBuf ); 
ElError ElCauchyLike_c
( ElMatrix_c A, ElInt rSize, complex_float* rBuf,
                ElInt sSize, complex_float* sBuf,
                ElInt xSize, complex_float* xBuf,
                ElInt ySize, complex_float* yBuf ); 
ElError ElCauchyLike_z
( ElMatrix_z A, ElInt rSize, complex_double* rBuf,
                ElInt sSize, complex_double* sBuf,
                ElInt xSize, complex_double* xBuf,
                ElInt ySize, complex_double* yBuf ); 

ElError ElCauchyLikeDist_s
( ElDistMatrix_s A, ElInt rSize, float* rBuf,
                    ElInt sSize, float* sBuf,
                    ElInt xSize, float* xBuf,
                    ElInt ySize, float* yBuf ); 
ElError ElCauchyLikeDist_d
( ElDistMatrix_d A, ElInt rSize, double* rBuf,
                    ElInt sSize, double* sBuf,
                    ElInt xSize, double* xBuf,
                    ElInt ySize, double* yBuf ); 
ElError ElCauchyLikeDist_c
( ElDistMatrix_c A, ElInt rSize, complex_float* rBuf,
                    ElInt sSize, complex_float* sBuf,
                    ElInt xSize, complex_float* xBuf,
                    ElInt ySize, complex_float* yBuf ); 
ElError ElCauchyLikeDist_z
( ElDistMatrix_z A, ElInt rSize, complex_double* rBuf,
                    ElInt sSize, complex_double* sBuf,
                    ElInt xSize, complex_double* xBuf,
                    ElInt ySize, complex_double* yBuf ); 

/* Circulant
   ========= */
ElError ElCirculant_i( ElMatrix_i A, ElInt aSize, ElInt* aBuf );
ElError ElCirculant_s( ElMatrix_s A, ElInt aSize, float* aBuf );
ElError ElCirculant_d( ElMatrix_d A, ElInt aSize, double* aBuf );
ElError ElCirculant_c( ElMatrix_c A, ElInt aSize, complex_float* aBuf );
ElError ElCirculant_z( ElMatrix_z A, ElInt aSize, complex_double* aBuf );

ElError ElCirculantDist_i
( ElDistMatrix_i A, ElInt aSize, ElInt* aBuf );
ElError ElCirculantDist_s
( ElDistMatrix_s A, ElInt aSize, float* aBuf );
ElError ElCirculantDist_d
( ElDistMatrix_d A, ElInt aSize, double* aBuf );
ElError ElCirculantDist_c
( ElDistMatrix_c A, ElInt aSize, complex_float* aBuf );
ElError ElCirculantDist_z
( ElDistMatrix_z A, ElInt aSize, complex_double* aBuf );

/* Demmel
   ====== */
ElError ElDemmel_s( ElMatrix_s A, ElInt n );
ElError ElDemmel_d( ElMatrix_d A, ElInt n );
ElError ElDemmel_c( ElMatrix_c A, ElInt n );
ElError ElDemmel_z( ElMatrix_z A, ElInt n );

ElError ElDemmelDist_s( ElDistMatrix_s A, ElInt n );
ElError ElDemmelDist_d( ElDistMatrix_d A, ElInt n );
ElError ElDemmelDist_c( ElDistMatrix_c A, ElInt n );
ElError ElDemmelDist_z( ElDistMatrix_z A, ElInt n );

/* Diagonal
   ======== */
ElError ElDiagonal_i( ElMatrix_i A, ElInt dSize, ElInt* dBuf );
ElError ElDiagonal_s( ElMatrix_s A, ElInt dSize, float* dBuf );
ElError ElDiagonal_d( ElMatrix_d A, ElInt dSize, double* dBuf );
ElError ElDiagonal_c( ElMatrix_c A, ElInt dSize, complex_float* dBuf );
ElError ElDiagonal_z( ElMatrix_z A, ElInt dSize, complex_double* dBuf );

ElError ElDiagonalDist_i( ElDistMatrix_i A, ElInt dSize, ElInt* dBuf );
ElError ElDiagonalDist_s( ElDistMatrix_s A, ElInt dSize, float* dBuf );
ElError ElDiagonalDist_d( ElDistMatrix_d A, ElInt dSize, double* dBuf );
ElError ElDiagonalDist_c( ElDistMatrix_c A, ElInt dSize, complex_float* dBuf );
ElError ElDiagonalDist_z( ElDistMatrix_z A, ElInt dSize, complex_double* dBuf );

/* Egorov
   ====== */
ElError ElEgorov_c( ElMatrix_c A, float (*phase)( ElInt, ElInt ), ElInt n );
ElError ElEgorov_z( ElMatrix_z A, double (*phase)( ElInt, ElInt ), ElInt n );

ElError ElEgorovDist_c
( ElDistMatrix_c A, float (*phase)( ElInt, ElInt ), ElInt n );
ElError ElEgorovDist_z
( ElDistMatrix_z A, double (*phase)( ElInt, ElInt ), ElInt n );

/* Ehrenfest
   ========= */
ElError ElEhrenfest_s( ElMatrix_s P, ElInt n );
ElError ElEhrenfest_d( ElMatrix_d P, ElInt n );
ElError ElEhrenfest_c( ElMatrix_c P, ElInt n );
ElError ElEhrenfest_z( ElMatrix_z P, ElInt n );

ElError ElEhrenfestDist_s( ElDistMatrix_s P, ElInt n );
ElError ElEhrenfestDist_d( ElDistMatrix_d P, ElInt n );
ElError ElEhrenfestDist_c( ElDistMatrix_c P, ElInt n );
ElError ElEhrenfestDist_z( ElDistMatrix_z P, ElInt n );

ElError ElEhrenfestStationary_s( ElMatrix_s PInf, ElInt n );
ElError ElEhrenfestStationary_d( ElMatrix_d PInf, ElInt n );
ElError ElEhrenfestStationary_c( ElMatrix_c PInf, ElInt n );
ElError ElEhrenfestStationary_z( ElMatrix_z PInf, ElInt n );

ElError ElEhrenfestStationaryDist_s( ElDistMatrix_s PInf, ElInt n );
ElError ElEhrenfestStationaryDist_d( ElDistMatrix_d PInf, ElInt n );
ElError ElEhrenfestStationaryDist_c( ElDistMatrix_c PInf, ElInt n );
ElError ElEhrenfestStationaryDist_z( ElDistMatrix_z PInf, ElInt n );

ElError ElEhrenfestDecay_s( ElMatrix_s A, ElInt n );
ElError ElEhrenfestDecay_d( ElMatrix_d A, ElInt n );
ElError ElEhrenfestDecay_c( ElMatrix_c A, ElInt n );
ElError ElEhrenfestDecay_z( ElMatrix_z A, ElInt n );

/* TODO: Distributed EhrenfestDecay */

/* ExtendedKahan
   ============= */
ElError ElExtendedKahan_s( ElMatrix_s A, ElInt k, float phi, float mu );
ElError ElExtendedKahan_d( ElMatrix_d A, ElInt k, double phi, double mu );
ElError ElExtendedKahan_c( ElMatrix_c A, ElInt k, float phi, float mu );
ElError ElExtendedKahan_z( ElMatrix_z A, ElInt k, double phi, double mu );

/* TODO: Distributed ExtendedKahan */

/* Fiedler
   ======= */
ElError ElFiedler_s( ElMatrix_s A, ElInt cSize, float* cBuf );
ElError ElFiedler_d( ElMatrix_d A, ElInt cSize, double* cBuf );
ElError ElFiedler_c( ElMatrix_c A, ElInt cSize, complex_float* cBuf );
ElError ElFiedler_z( ElMatrix_z A, ElInt cSize, complex_double* cBuf );

ElError ElFiedlerDist_s( ElDistMatrix_s A, ElInt cSize, float* cBuf );
ElError ElFiedlerDist_d( ElDistMatrix_d A, ElInt cSize, double* cBuf );
ElError ElFiedlerDist_c( ElDistMatrix_c A, ElInt cSize, complex_float* cBuf );
ElError ElFiedlerDist_z( ElDistMatrix_z A, ElInt cSize, complex_double* cBuf );

/* Forsythe
   ======== */
ElError ElForsythe_i
( ElMatrix_i J, ElInt n, ElInt alpha, ElInt lambda );
ElError ElForsythe_s
( ElMatrix_s J, ElInt n, float alpha, float lambda );
ElError ElForsythe_d
( ElMatrix_d J, ElInt n, double alpha, double lambda );
ElError ElForsythe_c
( ElMatrix_c J, ElInt n, complex_float alpha, complex_float lambda );
ElError ElForsythe_z
( ElMatrix_z J, ElInt n, complex_double alpha, complex_double lambda );

ElError ElForsytheDist_i
( ElDistMatrix_i J, ElInt n, ElInt alpha, ElInt lambda );
ElError ElForsytheDist_s
( ElDistMatrix_s J, ElInt n, float alpha, float lambda );
ElError ElForsytheDist_d
( ElDistMatrix_d J, ElInt n, double alpha, double lambda );
ElError ElForsytheDist_c
( ElDistMatrix_c J, ElInt n, complex_float alpha, complex_float lambda );
ElError ElForsytheDist_z
( ElDistMatrix_z J, ElInt n, complex_double alpha, complex_double lambda );

/* Ones
   ==== */
ElError ElOnes_i( ElMatrix_i A, ElInt m, ElInt n );
ElError ElOnes_s( ElMatrix_s A, ElInt m, ElInt n );
ElError ElOnes_d( ElMatrix_d A, ElInt m, ElInt n );
ElError ElOnes_c( ElMatrix_c A, ElInt m, ElInt n );
ElError ElOnes_z( ElMatrix_z A, ElInt m, ElInt n );

ElError ElOnesDist_i( ElDistMatrix_i A, ElInt m, ElInt n );
ElError ElOnesDist_s( ElDistMatrix_s A, ElInt m, ElInt n );
ElError ElOnesDist_d( ElDistMatrix_d A, ElInt m, ElInt n );
ElError ElOnesDist_c( ElDistMatrix_c A, ElInt m, ElInt n );
ElError ElOnesDist_z( ElDistMatrix_z A, ElInt m, ElInt n );

/* Uniform
   ======= */
ElError ElUniform_i
( ElMatrix_i A, ElInt m, ElInt n, ElInt center, ElInt radius );
ElError ElUniform_s
( ElMatrix_s A, ElInt m, ElInt n, float center, float radius );
ElError ElUniform_d
( ElMatrix_d A, ElInt m, ElInt n, double center, double radius );
ElError ElUniform_c
( ElMatrix_c A, ElInt m, ElInt n, complex_float center, float radius );
ElError ElUniform_z
( ElMatrix_z A, ElInt m, ElInt n, complex_double center, double radius );

ElError ElUniformDist_i
( ElDistMatrix_i A, ElInt m, ElInt n, ElInt center, ElInt radius );
ElError ElUniformDist_s
( ElDistMatrix_s A, ElInt m, ElInt n, float center, float radius );
ElError ElUniformDist_d
( ElDistMatrix_d A, ElInt m, ElInt n, double center, double radius );
ElError ElUniformDist_c
( ElDistMatrix_c A, ElInt m, ElInt n, complex_float center, float radius );
ElError ElUniformDist_z
( ElDistMatrix_z A, ElInt m, ElInt n, complex_double center, double radius );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_MATRICES_C_H */
