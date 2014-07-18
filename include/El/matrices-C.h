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

/* Fox-Li
   ====== */
ElError ElFoxLi_c( ElMatrix_c A, ElInt n, float omega );
ElError ElFoxLi_z( ElMatrix_z A, ElInt n, double omega );

/* TODO: Distributed Fox-Li */

/* Fourier
   ======= */
ElError ElFourier_c( ElMatrix_c A, ElInt n );
ElError ElFourier_z( ElMatrix_z A, ElInt n );

ElError ElFourierDist_c( ElDistMatrix_c A, ElInt n );
ElError ElFourierDist_z( ElDistMatrix_z A, ElInt n );

/* GCD matrix
   ========== */
ElError ElGCDMatrix_i( ElMatrix_i G, ElInt m, ElInt n );
ElError ElGCDMatrix_s( ElMatrix_s G, ElInt m, ElInt n );
ElError ElGCDMatrix_d( ElMatrix_d G, ElInt m, ElInt n );
ElError ElGCDMatrix_c( ElMatrix_c G, ElInt m, ElInt n );
ElError ElGCDMatrix_z( ElMatrix_z G, ElInt m, ElInt n );

ElError ElGCDMatrixDist_i( ElDistMatrix_i G, ElInt m, ElInt n );
ElError ElGCDMatrixDist_s( ElDistMatrix_s G, ElInt m, ElInt n );
ElError ElGCDMatrixDist_d( ElDistMatrix_d G, ElInt m, ElInt n );
ElError ElGCDMatrixDist_c( ElDistMatrix_c G, ElInt m, ElInt n );
ElError ElGCDMatrixDist_z( ElDistMatrix_z G, ElInt m, ElInt n );

/* Gear matrix
   =========== */
ElError ElGear_i( ElMatrix_i G, ElInt n, ElInt s, ElInt t );
ElError ElGear_s( ElMatrix_s G, ElInt n, ElInt s, ElInt t );
ElError ElGear_d( ElMatrix_d G, ElInt n, ElInt s, ElInt t );
ElError ElGear_c( ElMatrix_c G, ElInt n, ElInt s, ElInt t );
ElError ElGear_z( ElMatrix_z G, ElInt n, ElInt s, ElInt t );

ElError ElGearDist_i( ElDistMatrix_i G, ElInt n, ElInt s, ElInt t );
ElError ElGearDist_s( ElDistMatrix_s G, ElInt n, ElInt s, ElInt t );
ElError ElGearDist_d( ElDistMatrix_d G, ElInt n, ElInt s, ElInt t );
ElError ElGearDist_c( ElDistMatrix_c G, ElInt n, ElInt s, ElInt t );
ElError ElGearDist_z( ElDistMatrix_z G, ElInt n, ElInt s, ElInt t );

/* Golub Klema Stewart
   =================== */
ElError ElGKS_s( ElMatrix_s A, ElInt n );
ElError ElGKS_d( ElMatrix_d A, ElInt n );
ElError ElGKS_c( ElMatrix_c A, ElInt n );
ElError ElGKS_z( ElMatrix_z A, ElInt n );

ElError ElGKSDist_s( ElDistMatrix_s A, ElInt n );
ElError ElGKSDist_d( ElDistMatrix_d A, ElInt n );
ElError ElGKSDist_c( ElDistMatrix_c A, ElInt n );
ElError ElGKSDist_z( ElDistMatrix_z A, ElInt n );

/* Grcar
   ===== */
ElError ElGrcar_i( ElMatrix_i A, ElInt n, ElInt k );
ElError ElGrcar_s( ElMatrix_s A, ElInt n, ElInt k );
ElError ElGrcar_d( ElMatrix_d A, ElInt n, ElInt k );
ElError ElGrcar_c( ElMatrix_c A, ElInt n, ElInt k );
ElError ElGrcar_z( ElMatrix_z A, ElInt n, ElInt k );

ElError ElGrcarDist_i( ElDistMatrix_i A, ElInt n, ElInt k );
ElError ElGrcarDist_s( ElDistMatrix_s A, ElInt n, ElInt k );
ElError ElGrcarDist_d( ElDistMatrix_d A, ElInt n, ElInt k );
ElError ElGrcarDist_c( ElDistMatrix_c A, ElInt n, ElInt k );
ElError ElGrcarDist_z( ElDistMatrix_z A, ElInt n, ElInt k );

/* Hankel
   ====== */
ElError ElHankel_i
( ElMatrix_i A, ElInt m, ElInt n, ElInt aSize, ElInt* aBuf );
ElError ElHankel_s
( ElMatrix_s A, ElInt m, ElInt n, ElInt aSize, float* aBuf );
ElError ElHankel_d
( ElMatrix_d A, ElInt m, ElInt n, ElInt aSize, double* aBuf );
ElError ElHankel_c
( ElMatrix_c A, ElInt m, ElInt n, ElInt aSize, complex_float* aBuf );
ElError ElHankel_z
( ElMatrix_z A, ElInt m, ElInt n, ElInt aSize, complex_double* aBuf );

ElError ElHankelDist_i
( ElDistMatrix_i A, ElInt m, ElInt n, ElInt aSize, ElInt* aBuf );
ElError ElHankelDist_s
( ElDistMatrix_s A, ElInt m, ElInt n, ElInt aSize, float* aBuf );
ElError ElHankelDist_d
( ElDistMatrix_d A, ElInt m, ElInt n, ElInt aSize, double* aBuf );
ElError ElHankelDist_c
( ElDistMatrix_c A, ElInt m, ElInt n, ElInt aSize, complex_float* aBuf );
ElError ElHankelDist_z
( ElDistMatrix_z A, ElInt m, ElInt n, ElInt aSize, complex_double* aBuf );

/* Hanowa
   ====== */
ElError ElHanowa_i( ElMatrix_i A, ElInt n, ElInt mu );
ElError ElHanowa_s( ElMatrix_s A, ElInt n, float mu );
ElError ElHanowa_d( ElMatrix_d A, ElInt n, double mu );
ElError ElHanowa_c( ElMatrix_c A, ElInt n, complex_float mu );
ElError ElHanowa_z( ElMatrix_z A, ElInt n, complex_double mu );

/* TODO: Distributed Hanowa */

/* Hatano-Nelson
   ============= */
ElError ElHatanoNelson_s
( ElMatrix_s A, ElInt n, float center, float radius, float g, 
  bool periodic );
ElError ElHatanoNelson_d
( ElMatrix_d A, ElInt n, double center, double radius, double g, 
  bool periodic );
ElError ElHatanoNelson_c
( ElMatrix_c A, ElInt n, complex_float center, float radius, complex_float g, 
  bool periodic );
ElError ElHatanoNelson_z
( ElMatrix_z A, ElInt n, complex_double center, double radius, complex_double g,
  bool periodic );

/* TODO: Distributed Hatano-Nelson */

/* Helmholtz
   ========= */
ElError ElHelmholtz1D_s( ElMatrix_s H, ElInt nx, float shift );
ElError ElHelmholtz1D_d( ElMatrix_d H, ElInt nx, double shift );
ElError ElHelmholtz1D_c( ElMatrix_c H, ElInt nx, complex_float shift );
ElError ElHelmholtz1D_z( ElMatrix_z H, ElInt nx, complex_double shift );

ElError ElHelmholtz1DDist_s( ElDistMatrix_s H, ElInt nx, float shift );
ElError ElHelmholtz1DDist_d( ElDistMatrix_d H, ElInt nx, double shift );
ElError ElHelmholtz1DDist_c( ElDistMatrix_c H, ElInt nx, complex_float shift );
ElError ElHelmholtz1DDist_z( ElDistMatrix_z H, ElInt nx, complex_double shift );

ElError ElHelmholtz2D_s
( ElMatrix_s H, ElInt nx, ElInt ny, float shift );
ElError ElHelmholtz2D_d
( ElMatrix_d H, ElInt nx, ElInt ny, double shift );
ElError ElHelmholtz2D_c
( ElMatrix_c H, ElInt nx, ElInt ny, complex_float shift );
ElError ElHelmholtz2D_z
( ElMatrix_z H, ElInt nx, ElInt ny, complex_double shift );

ElError ElHelmholtz2DDist_s
( ElDistMatrix_s H, ElInt nx, ElInt ny, float shift );
ElError ElHelmholtz2DDist_d
( ElDistMatrix_d H, ElInt nx, ElInt ny, double shift );
ElError ElHelmholtz2DDist_c
( ElDistMatrix_c H, ElInt nx, ElInt ny, complex_float shift );
ElError ElHelmholtz2DDist_z
( ElDistMatrix_z H, ElInt nx, ElInt ny, complex_double shift );

ElError ElHelmholtz3D_s
( ElMatrix_s H, ElInt nx, ElInt ny, ElInt nz, float shift );
ElError ElHelmholtz3D_d
( ElMatrix_d H, ElInt nx, ElInt ny, ElInt nz, double shift );
ElError ElHelmholtz3D_c
( ElMatrix_c H, ElInt nx, ElInt ny, ElInt nz, complex_float shift );
ElError ElHelmholtz3D_z
( ElMatrix_z H, ElInt nx, ElInt ny, ElInt nz, complex_double shift );

ElError ElHelmholtz3DDist_s
( ElDistMatrix_s H, ElInt nx, ElInt ny, ElInt nz, float shift );
ElError ElHelmholtz3DDist_d
( ElDistMatrix_d H, ElInt nx, ElInt ny, ElInt nz, double shift );
ElError ElHelmholtz3DDist_c
( ElDistMatrix_c H, ElInt nx, ElInt ny, ElInt nz, complex_float shift );
ElError ElHelmholtz3DDist_z
( ElDistMatrix_z H, ElInt nx, ElInt ny, ElInt nz, complex_double shift );

/* Helmholtz with PML
   ================== */
ElError ElHelmholtzPML1D_c
( ElMatrix_c H, ElInt nx, complex_float omega, 
  ElInt numPmlPoints, float sigma, float pmlExp );
ElError ElHelmholtzPML1D_z
( ElMatrix_z H, ElInt nx, complex_double omega, 
  ElInt numPmlPoints, double sigma, double pmlExp );

ElError ElHelmholtzPML1DDist_c
( ElDistMatrix_c H, ElInt nx, complex_float omega, 
  ElInt numPmlPoints, float sigma, float pmlExp );
ElError ElHelmholtzPML1DDist_z
( ElDistMatrix_z H, ElInt nx, complex_double omega, 
  ElInt numPmlPoints, double sigma, double pmlExp );

ElError ElHelmholtzPML2D_c
( ElMatrix_c H, ElInt nx, ElInt ny, complex_float omega, 
  ElInt numPmlPoints, float sigma, float pmlExp );
ElError ElHelmholtzPML2D_z
( ElMatrix_z H, ElInt nx, ElInt ny, complex_double omega, 
  ElInt numPmlPoints, double sigma, double pmlExp );

ElError ElHelmholtzPML2DDist_c
( ElDistMatrix_c H, ElInt nx, ElInt ny, complex_float omega, 
  ElInt numPmlPoints, float sigma, float pmlExp );
ElError ElHelmholtzPML2DDist_z
( ElDistMatrix_z H, ElInt nx, ElInt ny, complex_double omega, 
  ElInt numPmlPoints, double sigma, double pmlExp );

ElError ElHelmholtzPML3D_c
( ElMatrix_c H, ElInt nx, ElInt ny, ElInt nz, complex_float omega, 
  ElInt numPmlPoints, float sigma, float pmlExp );
ElError ElHelmholtzPML3D_z
( ElMatrix_z H, ElInt nx, ElInt ny, ElInt nz, complex_double omega, 
  ElInt numPmlPoints, double sigma, double pmlExp );

ElError ElHelmholtzPML3DDist_c
( ElDistMatrix_c H, ElInt nx, ElInt ny, ElInt nz, complex_float omega, 
  ElInt numPmlPoints, float sigma, float pmlExp );
ElError ElHelmholtzPML3DDist_z
( ElDistMatrix_z H, ElInt nx, ElInt ny, ElInt nz, complex_double omega, 
  ElInt numPmlPoints, double sigma, double pmlExp );

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
