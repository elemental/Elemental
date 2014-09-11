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

ElError ElEhrenfestDecayDist_s( ElDistMatrix_s A, ElInt n );
ElError ElEhrenfestDecayDist_d( ElDistMatrix_d A, ElInt n );
ElError ElEhrenfestDecayDist_c( ElDistMatrix_c A, ElInt n );
ElError ElEhrenfestDecayDist_z( ElDistMatrix_z A, ElInt n );

/* ExtendedKahan
   ============= */
ElError ElExtendedKahan_s( ElMatrix_s A, ElInt k, float phi, float mu );
ElError ElExtendedKahan_d( ElMatrix_d A, ElInt k, double phi, double mu );
ElError ElExtendedKahan_c( ElMatrix_c A, ElInt k, float phi, float mu );
ElError ElExtendedKahan_z( ElMatrix_z A, ElInt k, double phi, double mu );

ElError ElExtendedKahanDist_s
( ElDistMatrix_s A, ElInt k, float phi, float mu );
ElError ElExtendedKahanDist_d
( ElDistMatrix_d A, ElInt k, double phi, double mu );
ElError ElExtendedKahanDist_c
( ElDistMatrix_c A, ElInt k, float phi, float mu );
ElError ElExtendedKahanDist_z
( ElDistMatrix_z A, ElInt k, double phi, double mu );

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

ElError ElFoxLiDist_c( ElDistMatrix_c A, ElInt n, float omega );
ElError ElFoxLiDist_z( ElDistMatrix_z A, ElInt n, double omega );

/* Fourier
   ======= */
ElError ElFourier_c( ElMatrix_c A, ElInt n );
ElError ElFourier_z( ElMatrix_z A, ElInt n );

ElError ElFourierDist_c( ElDistMatrix_c A, ElInt n );
ElError ElFourierDist_z( ElDistMatrix_z A, ElInt n );

/* Gaussian
   ======== */
ElError ElGaussian_s
( ElMatrix_s A, ElInt m, ElInt n, float mean, float stddev );
ElError ElGaussian_d
( ElMatrix_d A, ElInt m, ElInt n, double mean, double stddev );
ElError ElGaussian_c
( ElMatrix_c A, ElInt m, ElInt n, complex_float mean, float stddev );
ElError ElGaussian_z
( ElMatrix_z A, ElInt m, ElInt n, complex_double mean, double stddev );

ElError ElGaussianDist_s
( ElDistMatrix_s A, ElInt m, ElInt n, float mean, float stddev );
ElError ElGaussianDist_d
( ElDistMatrix_d A, ElInt m, ElInt n, double mean, double stddev );
ElError ElGaussianDist_c
( ElDistMatrix_c A, ElInt m, ElInt n, complex_float mean, float stddev );
ElError ElGaussianDist_z
( ElDistMatrix_z A, ElInt m, ElInt n, complex_double mean, double stddev );

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

/* Haar 
   ==== */
ElError ElHaar_s( ElMatrix_s A, ElInt n );
ElError ElHaar_d( ElMatrix_d A, ElInt n );
ElError ElHaar_c( ElMatrix_c A, ElInt n );
ElError ElHaar_z( ElMatrix_z A, ElInt n );

/* TODO: Distributed Haar */

ElError ElImplicitHaar_s( ElMatrix_s A, ElMatrix_s t, ElMatrix_s d, ElInt n );
ElError ElImplicitHaar_d( ElMatrix_d A, ElMatrix_d t, ElMatrix_d d, ElInt n );
ElError ElImplicitHaar_c( ElMatrix_c A, ElMatrix_c t, ElMatrix_s d, ElInt n );
ElError ElImplicitHaar_z( ElMatrix_z A, ElMatrix_z t, ElMatrix_d d, ElInt n );

/* TODO: Distributed implicit Haar */

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

ElError ElHanowaDist_i( ElDistMatrix_i A, ElInt n, ElInt mu );
ElError ElHanowaDist_s( ElDistMatrix_s A, ElInt n, float mu );
ElError ElHanowaDist_d( ElDistMatrix_d A, ElInt n, double mu );
ElError ElHanowaDist_c( ElDistMatrix_c A, ElInt n, complex_float mu );
ElError ElHanowaDist_z( ElDistMatrix_z A, ElInt n, complex_double mu );

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

/* Hermitian from EVD
   ================== */
ElError ElHermitianFromEVD_s
( ElUpperOrLower uplo, ElMatrix_s A, 
  ElConstMatrix_s w, ElConstMatrix_s Z );
ElError ElHermitianFromEVD_d
( ElUpperOrLower uplo, ElMatrix_d A, 
  ElConstMatrix_d w, ElConstMatrix_d Z );
ElError ElHermitianFromEVD_c
( ElUpperOrLower uplo, ElMatrix_c A, 
  ElConstMatrix_s w, ElConstMatrix_c Z );
ElError ElHermitianFromEVD_z
( ElUpperOrLower uplo, ElMatrix_z A, 
  ElConstMatrix_d w, ElConstMatrix_z Z );

/* TODO: Distributed HermitianFromEVD */

/* Hermitian uniform spectrum
   ========================== */
ElError ElHermitianUniformSpectrum_s
( ElMatrix_s A, ElInt n, float lower, float upper );
ElError ElHermitianUniformSpectrum_d
( ElMatrix_d A, ElInt n, double lower, double upper );
ElError ElHermitianUniformSpectrum_c
( ElMatrix_c A, ElInt n, float lower, float upper );
ElError ElHermitianUniformSpectrum_z
( ElMatrix_z A, ElInt n, double lower, double upper );

ElError ElHermitianUniformSpectrumDist_s
( ElDistMatrix_s A, ElInt n, float lower, float upper );
ElError ElHermitianUniformSpectrumDist_d
( ElDistMatrix_d A, ElInt n, double lower, double upper );
ElError ElHermitianUniformSpectrumDist_c
( ElDistMatrix_c A, ElInt n, float lower, float upper );
ElError ElHermitianUniformSpectrumDist_z
( ElDistMatrix_z A, ElInt n, double lower, double upper );

/* Hilbert
   ======= */
ElError ElHilbert_s( ElMatrix_s A, ElInt n );
ElError ElHilbert_d( ElMatrix_d A, ElInt n );
ElError ElHilbert_c( ElMatrix_c A, ElInt n );
ElError ElHilbert_z( ElMatrix_z A, ElInt n );

ElError ElHilbertDist_s( ElDistMatrix_s A, ElInt n );
ElError ElHilbertDist_d( ElDistMatrix_d A, ElInt n );
ElError ElHilbertDist_c( ElDistMatrix_c A, ElInt n );
ElError ElHilbertDist_z( ElDistMatrix_z A, ElInt n );

/* Identity
   ======== */
ElError ElIdentity_i( ElMatrix_i A, ElInt m, ElInt n );
ElError ElIdentity_s( ElMatrix_s A, ElInt m, ElInt n );
ElError ElIdentity_d( ElMatrix_d A, ElInt m, ElInt n );
ElError ElIdentity_c( ElMatrix_c A, ElInt m, ElInt n );
ElError ElIdentity_z( ElMatrix_z A, ElInt m, ElInt n );

ElError ElIdentityDist_i( ElDistMatrix_i A, ElInt m, ElInt n );
ElError ElIdentityDist_s( ElDistMatrix_s A, ElInt m, ElInt n );
ElError ElIdentityDist_d( ElDistMatrix_d A, ElInt m, ElInt n );
ElError ElIdentityDist_c( ElDistMatrix_c A, ElInt m, ElInt n );
ElError ElIdentityDist_z( ElDistMatrix_z A, ElInt m, ElInt n );

/* Jordan
   ====== */
ElError ElJordan_i( ElMatrix_i J, ElInt n, ElInt lambda );
ElError ElJordan_s( ElMatrix_s J, ElInt n, float lambda );
ElError ElJordan_d( ElMatrix_d J, ElInt n, double lambda );
ElError ElJordan_c( ElMatrix_c J, ElInt n, complex_float lambda );
ElError ElJordan_z( ElMatrix_z J, ElInt n, complex_double lambda );

ElError ElJordanDist_i( ElDistMatrix_i J, ElInt n, ElInt lambda );
ElError ElJordanDist_s( ElDistMatrix_s J, ElInt n, float lambda );
ElError ElJordanDist_d( ElDistMatrix_d J, ElInt n, double lambda );
ElError ElJordanDist_c( ElDistMatrix_c J, ElInt n, complex_float lambda );
ElError ElJordanDist_z( ElDistMatrix_z J, ElInt n, complex_double lambda );

/* Kahan
   ===== */
ElError ElKahan_s( ElMatrix_s A, ElInt n, float phi );
ElError ElKahan_d( ElMatrix_d A, ElInt n, double phi );
ElError ElKahan_c( ElMatrix_c A, ElInt n, complex_float phi );
ElError ElKahan_z( ElMatrix_z A, ElInt n, complex_double phi );

ElError ElKahanDist_s( ElDistMatrix_s A, ElInt n, float phi );
ElError ElKahanDist_d( ElDistMatrix_d A, ElInt n, double phi );
ElError ElKahanDist_c( ElDistMatrix_c A, ElInt n, complex_float phi );
ElError ElKahanDist_z( ElDistMatrix_z A, ElInt n, complex_double phi );

/* KMS
   === */
ElError ElKMS_i( ElMatrix_i K, ElInt n, ElInt rho );
ElError ElKMS_s( ElMatrix_s K, ElInt n, float rho );
ElError ElKMS_d( ElMatrix_d K, ElInt n, double rho );
ElError ElKMS_c( ElMatrix_c K, ElInt n, complex_float rho );
ElError ElKMS_z( ElMatrix_z K, ElInt n, complex_double rho );

ElError ElKMSDist_i( ElDistMatrix_i K, ElInt n, ElInt rho );
ElError ElKMSDist_s( ElDistMatrix_s K, ElInt n, float rho );
ElError ElKMSDist_d( ElDistMatrix_d K, ElInt n, double rho );
ElError ElKMSDist_c( ElDistMatrix_c K, ElInt n, complex_float rho );
ElError ElKMSDist_z( ElDistMatrix_z K, ElInt n, complex_double rho );

/* Laplacian
   ========= */
ElError ElLaplacian1D_s( ElMatrix_s L, ElInt nx );
ElError ElLaplacian1D_d( ElMatrix_d L, ElInt nx );
ElError ElLaplacian1D_c( ElMatrix_c L, ElInt nx );
ElError ElLaplacian1D_z( ElMatrix_z L, ElInt nx );

ElError ElLaplacian1DDist_s( ElDistMatrix_s L, ElInt nx );
ElError ElLaplacian1DDist_d( ElDistMatrix_d L, ElInt nx );
ElError ElLaplacian1DDist_c( ElDistMatrix_c L, ElInt nx );
ElError ElLaplacian1DDist_z( ElDistMatrix_z L, ElInt nx );

ElError ElLaplacian2D_s( ElMatrix_s L, ElInt nx, ElInt ny );
ElError ElLaplacian2D_d( ElMatrix_d L, ElInt nx, ElInt ny );
ElError ElLaplacian2D_c( ElMatrix_c L, ElInt nx, ElInt ny );
ElError ElLaplacian2D_z( ElMatrix_z L, ElInt nx, ElInt ny );

ElError ElLaplacian2DDist_s( ElDistMatrix_s L, ElInt nx, ElInt ny );
ElError ElLaplacian2DDist_d( ElDistMatrix_d L, ElInt nx, ElInt ny );
ElError ElLaplacian2DDist_c( ElDistMatrix_c L, ElInt nx, ElInt ny );
ElError ElLaplacian2DDist_z( ElDistMatrix_z L, ElInt nx, ElInt ny );

ElError ElLaplacian3D_s( ElMatrix_s L, ElInt nx, ElInt ny, ElInt nz );
ElError ElLaplacian3D_d( ElMatrix_d L, ElInt nx, ElInt ny, ElInt nz );
ElError ElLaplacian3D_c( ElMatrix_c L, ElInt nx, ElInt ny, ElInt nz );
ElError ElLaplacian3D_z( ElMatrix_z L, ElInt nx, ElInt ny, ElInt nz );

ElError ElLaplacian3DDist_s( ElDistMatrix_s L, ElInt nx, ElInt ny, ElInt nz );
ElError ElLaplacian3DDist_d( ElDistMatrix_d L, ElInt nx, ElInt ny, ElInt nz );
ElError ElLaplacian3DDist_c( ElDistMatrix_c L, ElInt nx, ElInt ny, ElInt nz );
ElError ElLaplacian3DDist_z( ElDistMatrix_z L, ElInt nx, ElInt ny, ElInt nz );

/* Lauchli
   ======= */
ElError ElLauchli_i( ElMatrix_i A, ElInt n, ElInt mu );
ElError ElLauchli_s( ElMatrix_s A, ElInt n, float mu );
ElError ElLauchli_d( ElMatrix_d A, ElInt n, double mu );
ElError ElLauchli_c( ElMatrix_c A, ElInt n, complex_float mu );
ElError ElLauchli_z( ElMatrix_z A, ElInt n, complex_double mu );

ElError ElLauchliDist_i( ElDistMatrix_i A, ElInt n, ElInt mu );
ElError ElLauchliDist_s( ElDistMatrix_s A, ElInt n, float mu );
ElError ElLauchliDist_d( ElDistMatrix_d A, ElInt n, double mu );
ElError ElLauchliDist_c( ElDistMatrix_c A, ElInt n, complex_float mu );
ElError ElLauchliDist_z( ElDistMatrix_z A, ElInt n, complex_double mu );

/* Legendre
   ======== */
ElError ElLegendre_s( ElMatrix_s A, ElInt n );
ElError ElLegendre_d( ElMatrix_d A, ElInt n );
ElError ElLegendre_c( ElMatrix_c A, ElInt n );
ElError ElLegendre_z( ElMatrix_z A, ElInt n );

ElError ElLegendreDist_s( ElDistMatrix_s A, ElInt n );
ElError ElLegendreDist_d( ElDistMatrix_d A, ElInt n );
ElError ElLegendreDist_c( ElDistMatrix_c A, ElInt n );
ElError ElLegendreDist_z( ElDistMatrix_z A, ElInt n );

/* Lehmer
   ====== */
ElError ElLehmer_s( ElMatrix_s L, ElInt n );
ElError ElLehmer_d( ElMatrix_d L, ElInt n );
ElError ElLehmer_c( ElMatrix_c L, ElInt n );
ElError ElLehmer_z( ElMatrix_z L, ElInt n );

ElError ElLehmerDist_s( ElDistMatrix_s L, ElInt n );
ElError ElLehmerDist_d( ElDistMatrix_d L, ElInt n );
ElError ElLehmerDist_c( ElDistMatrix_c L, ElInt n );
ElError ElLehmerDist_z( ElDistMatrix_z L, ElInt n );

/* Lotkin
   ====== */
ElError ElLotkin_s( ElMatrix_s A, ElInt n );
ElError ElLotkin_d( ElMatrix_d A, ElInt n );
ElError ElLotkin_c( ElMatrix_c A, ElInt n );
ElError ElLotkin_z( ElMatrix_z A, ElInt n );

ElError ElLotkinDist_s( ElDistMatrix_s A, ElInt n );
ElError ElLotkinDist_d( ElDistMatrix_d A, ElInt n );
ElError ElLotkinDist_c( ElDistMatrix_c A, ElInt n );
ElError ElLotkinDist_z( ElDistMatrix_z A, ElInt n );

/* MinIJ
   ===== */
ElError ElMinIJ_i( ElMatrix_i A, ElInt n );
ElError ElMinIJ_s( ElMatrix_s A, ElInt n );
ElError ElMinIJ_d( ElMatrix_d A, ElInt n );
ElError ElMinIJ_c( ElMatrix_c A, ElInt n );
ElError ElMinIJ_z( ElMatrix_z A, ElInt n );

ElError ElMinIJDist_i( ElDistMatrix_i A, ElInt n );
ElError ElMinIJDist_s( ElDistMatrix_s A, ElInt n );
ElError ElMinIJDist_d( ElDistMatrix_d A, ElInt n );
ElError ElMinIJDist_c( ElDistMatrix_c A, ElInt n );
ElError ElMinIJDist_z( ElDistMatrix_z A, ElInt n );

/* NormalFromEVD
   ============= */
ElError ElNormalFromEVD_c( ElMatrix_c A, ElConstMatrix_c w, ElConstMatrix_c Z );
ElError ElNormalFromEVD_z( ElMatrix_z A, ElConstMatrix_z w, ElConstMatrix_z Z );

/* TODO: Distributed NormalFromEVD */

/* Normal uniform spectrum
   ======================= */
ElError ElNormalUniformSpectrum_c
( ElMatrix_c A, ElInt n, complex_float center, float radius );
ElError ElNormalUniformSpectrum_z
( ElMatrix_z A, ElInt n, complex_double center, double radius );

/* TODO: Distributed NormalUniformSpectrum */

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

/* 1-2-1
   ===== */
ElError ElOneTwoOne_i( ElMatrix_i A, ElInt n );
ElError ElOneTwoOne_s( ElMatrix_s A, ElInt n );
ElError ElOneTwoOne_d( ElMatrix_d A, ElInt n );
ElError ElOneTwoOne_c( ElMatrix_c A, ElInt n );
ElError ElOneTwoOne_z( ElMatrix_z A, ElInt n );

ElError ElOneTwoOneDist_i( ElDistMatrix_i A, ElInt n );
ElError ElOneTwoOneDist_s( ElDistMatrix_s A, ElInt n );
ElError ElOneTwoOneDist_d( ElDistMatrix_d A, ElInt n );
ElError ElOneTwoOneDist_c( ElDistMatrix_c A, ElInt n );
ElError ElOneTwoOneDist_z( ElDistMatrix_z A, ElInt n );

/* Parter
   ====== */
ElError ElParter_s( ElMatrix_s A, ElInt n );
ElError ElParter_d( ElMatrix_d A, ElInt n );
ElError ElParter_c( ElMatrix_c A, ElInt n );
ElError ElParter_z( ElMatrix_z A, ElInt n );

ElError ElParterDist_s( ElDistMatrix_s A, ElInt n );
ElError ElParterDist_d( ElDistMatrix_d A, ElInt n );
ElError ElParterDist_c( ElDistMatrix_c A, ElInt n );
ElError ElParterDist_z( ElDistMatrix_z A, ElInt n );

/* Pei
   === */
ElError ElPei_s( ElMatrix_s A, ElInt n, float alpha );
ElError ElPei_d( ElMatrix_d A, ElInt n, double alpha );
ElError ElPei_c( ElMatrix_c A, ElInt n, complex_float alpha );
ElError ElPei_z( ElMatrix_z A, ElInt n, complex_double alpha );

ElError ElPeiDist_s( ElDistMatrix_s A, ElInt n, float alpha );
ElError ElPeiDist_d( ElDistMatrix_d A, ElInt n, double alpha );
ElError ElPeiDist_c( ElDistMatrix_c A, ElInt n, complex_float alpha );
ElError ElPeiDist_z( ElDistMatrix_z A, ElInt n, complex_double alpha );

/* Redheffer
   ========= */
ElError ElRedheffer_i( ElMatrix_i A, ElInt n );
ElError ElRedheffer_s( ElMatrix_s A, ElInt n );
ElError ElRedheffer_d( ElMatrix_d A, ElInt n );
ElError ElRedheffer_c( ElMatrix_c A, ElInt n );
ElError ElRedheffer_z( ElMatrix_z A, ElInt n );

ElError ElRedhefferDist_i( ElDistMatrix_i A, ElInt n );
ElError ElRedhefferDist_s( ElDistMatrix_s A, ElInt n );
ElError ElRedhefferDist_d( ElDistMatrix_d A, ElInt n );
ElError ElRedhefferDist_c( ElDistMatrix_c A, ElInt n );
ElError ElRedhefferDist_z( ElDistMatrix_z A, ElInt n );

/* Riemann
   ======= */
ElError ElRiemann_i( ElMatrix_i A, ElInt n );
ElError ElRiemann_s( ElMatrix_s A, ElInt n );
ElError ElRiemann_d( ElMatrix_d A, ElInt n );
ElError ElRiemann_c( ElMatrix_c A, ElInt n );
ElError ElRiemann_z( ElMatrix_z A, ElInt n );

ElError ElRiemannDist_i( ElDistMatrix_i A, ElInt n );
ElError ElRiemannDist_s( ElDistMatrix_s A, ElInt n );
ElError ElRiemannDist_d( ElDistMatrix_d A, ElInt n );
ElError ElRiemannDist_c( ElDistMatrix_c A, ElInt n );
ElError ElRiemannDist_z( ElDistMatrix_z A, ElInt n );

/* Riffle
   ====== */
ElError ElRiffle_s( ElMatrix_s P, ElInt n );
ElError ElRiffle_d( ElMatrix_d P, ElInt n );
ElError ElRiffle_c( ElMatrix_c P, ElInt n );
ElError ElRiffle_z( ElMatrix_z P, ElInt n );

ElError ElRiffleDist_s( ElDistMatrix_s P, ElInt n );
ElError ElRiffleDist_d( ElDistMatrix_d P, ElInt n );
ElError ElRiffleDist_c( ElDistMatrix_c P, ElInt n );
ElError ElRiffleDist_z( ElDistMatrix_z P, ElInt n );

ElError ElRiffleStationary_s( ElMatrix_s PInf, ElInt n );
ElError ElRiffleStationary_d( ElMatrix_d PInf, ElInt n );
ElError ElRiffleStationary_c( ElMatrix_c PInf, ElInt n );
ElError ElRiffleStationary_z( ElMatrix_z PInf, ElInt n );

ElError ElRiffleStationaryDist_s( ElDistMatrix_s PInf, ElInt n );
ElError ElRiffleStationaryDist_d( ElDistMatrix_d PInf, ElInt n );
ElError ElRiffleStationaryDist_c( ElDistMatrix_c PInf, ElInt n );
ElError ElRiffleStationaryDist_z( ElDistMatrix_z PInf, ElInt n );

ElError ElRiffleDecay_s( ElMatrix_s A, ElInt n );
ElError ElRiffleDecay_d( ElMatrix_d A, ElInt n );
ElError ElRiffleDecay_c( ElMatrix_c A, ElInt n );
ElError ElRiffleDecay_z( ElMatrix_z A, ElInt n );

ElError ElRiffleDecayDist_s( ElDistMatrix_s A, ElInt n );
ElError ElRiffleDecayDist_d( ElDistMatrix_d A, ElInt n );
ElError ElRiffleDecayDist_c( ElDistMatrix_c A, ElInt n );
ElError ElRiffleDecayDist_z( ElDistMatrix_z A, ElInt n );

/* Ris
   === */
ElError ElRis_s( ElMatrix_s A, ElInt n );
ElError ElRis_d( ElMatrix_d A, ElInt n );
ElError ElRis_c( ElMatrix_c A, ElInt n );
ElError ElRis_z( ElMatrix_z A, ElInt n );

ElError ElRisDist_s( ElDistMatrix_s A, ElInt n );
ElError ElRisDist_d( ElDistMatrix_d A, ElInt n );
ElError ElRisDist_c( ElDistMatrix_c A, ElInt n );
ElError ElRisDist_z( ElDistMatrix_z A, ElInt n );

/* Toeplitz
   ======== */
ElError ElToeplitz_i
( ElMatrix_i A, ElInt m, ElInt n, ElInt aSize, ElInt* aBuf );
ElError ElToeplitz_s
( ElMatrix_s A, ElInt m, ElInt n, ElInt aSize, float* aBuf );
ElError ElToeplitz_d
( ElMatrix_d A, ElInt m, ElInt n, ElInt aSize, double* aBuf );
ElError ElToeplitz_c
( ElMatrix_c A, ElInt m, ElInt n, ElInt aSize, complex_float* aBuf );
ElError ElToeplitz_z
( ElMatrix_z A, ElInt m, ElInt n, ElInt aSize, complex_double* aBuf );

ElError ElToeplitzDist_i
( ElDistMatrix_i A, ElInt m, ElInt n, ElInt aSize, ElInt* aBuf );
ElError ElToeplitzDist_s
( ElDistMatrix_s A, ElInt m, ElInt n, ElInt aSize, float* aBuf );
ElError ElToeplitzDist_d
( ElDistMatrix_d A, ElInt m, ElInt n, ElInt aSize, double* aBuf );
ElError ElToeplitzDist_c
( ElDistMatrix_c A, ElInt m, ElInt n, ElInt aSize, complex_float* aBuf );
ElError ElToeplitzDist_z
( ElDistMatrix_z A, ElInt m, ElInt n, ElInt aSize, complex_double* aBuf );

/* Trefethen
   ========= */
ElError ElTrefethen_c( ElMatrix_c A, ElInt n );
ElError ElTrefethen_z( ElMatrix_z A, ElInt n );

ElError ElTrefethenDist_c( ElDistMatrix_c A, ElInt n );
ElError ElTrefethenDist_z( ElDistMatrix_z A, ElInt n );

/* Triangle
   ======== */
ElError ElTriangle_s( ElMatrix_s A, ElInt n );
ElError ElTriangle_d( ElMatrix_d A, ElInt n );
ElError ElTriangle_c( ElMatrix_c A, ElInt n );
ElError ElTriangle_z( ElMatrix_z A, ElInt n );

ElError ElTriangleDist_s( ElDistMatrix_s A, ElInt n );
ElError ElTriangleDist_d( ElDistMatrix_d A, ElInt n );
ElError ElTriangleDist_c( ElDistMatrix_c A, ElInt n );
ElError ElTriangleDist_z( ElDistMatrix_z A, ElInt n );

/* TriW
   ==== */
ElError ElTriW_i
( ElMatrix_i A, ElInt m, ElInt n, ElInt alpha, ElInt k );
ElError ElTriW_s
( ElMatrix_s A, ElInt m, ElInt n, float alpha, ElInt k );
ElError ElTriW_d
( ElMatrix_d A, ElInt m, ElInt n, double alpha, ElInt k );
ElError ElTriW_c
( ElMatrix_c A, ElInt m, ElInt n, complex_float alpha, ElInt k );
ElError ElTriW_z
( ElMatrix_z A, ElInt m, ElInt n, complex_double alpha, ElInt k );

ElError ElTriWDist_i
( ElDistMatrix_i A, ElInt m, ElInt n, ElInt alpha, ElInt k );
ElError ElTriWDist_s
( ElDistMatrix_s A, ElInt m, ElInt n, float alpha, ElInt k );
ElError ElTriWDist_d
( ElDistMatrix_d A, ElInt m, ElInt n, double alpha, ElInt k );
ElError ElTriWDist_c
( ElDistMatrix_c A, ElInt m, ElInt n, complex_float alpha, ElInt k );
ElError ElTriWDist_z
( ElDistMatrix_z A, ElInt m, ElInt n, complex_double alpha, ElInt k );

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

/* Uniform Helmholtz Green's
   ========================= */
ElError ElUniformHelmholtzGreens_c( ElMatrix_c A, ElInt n, float lambda );
ElError ElUniformHelmholtzGreens_z( ElMatrix_z A, ElInt n, double lambda );

ElError ElUniformHelmholtzGreensDist_c
( ElDistMatrix_c A, ElInt n, float lambda );
ElError ElUniformHelmholtzGreensDist_z
( ElDistMatrix_z A, ElInt n, double lambda );

/* Walsh
   ===== */
ElError ElWalsh_i( ElMatrix_i A, ElInt k, bool binary );
ElError ElWalsh_s( ElMatrix_s A, ElInt k, bool binary );
ElError ElWalsh_d( ElMatrix_d A, ElInt k, bool binary );
ElError ElWalsh_c( ElMatrix_c A, ElInt k, bool binary );
ElError ElWalsh_z( ElMatrix_z A, ElInt k, bool binary );

ElError ElWalshDist_i( ElDistMatrix_i A, ElInt k, bool binary );
ElError ElWalshDist_s( ElDistMatrix_s A, ElInt k, bool binary );
ElError ElWalshDist_d( ElDistMatrix_d A, ElInt k, bool binary );
ElError ElWalshDist_c( ElDistMatrix_c A, ElInt k, bool binary );
ElError ElWalshDist_z( ElDistMatrix_z A, ElInt k, bool binary );

/* Whale
   ===== */
ElError ElWhale_c( ElMatrix_c A, ElInt n );
ElError ElWhale_z( ElMatrix_z A, ElInt n );

ElError ElWhaleDist_c( ElDistMatrix_c A, ElInt n );
ElError ElWhaleDist_z( ElDistMatrix_z A, ElInt n );

/* Wigner
   ====== */
ElError ElWigner_s
( ElMatrix_s A, ElInt n, float mean, float stddev );
ElError ElWigner_d
( ElMatrix_d A, ElInt n, double mean, double stddev );
ElError ElWigner_c
( ElMatrix_c A, ElInt n, complex_float mean, float stddev );
ElError ElWigner_z
( ElMatrix_z A, ElInt n, complex_double mean, double stddev );

ElError ElWignerDist_s
( ElDistMatrix_s A, ElInt n, float mean, float stddev );
ElError ElWignerDist_d
( ElDistMatrix_d A, ElInt n, double mean, double stddev );
ElError ElWignerDist_c
( ElDistMatrix_c A, ElInt n, complex_float mean, float stddev );
ElError ElWignerDist_z
( ElDistMatrix_z A, ElInt n, complex_double mean, double stddev );

/* Wilkinson
   ========= */
ElError ElWilkinson_i( ElMatrix_i A, ElInt k );
ElError ElWilkinson_s( ElMatrix_s A, ElInt k );
ElError ElWilkinson_d( ElMatrix_d A, ElInt k );
ElError ElWilkinson_c( ElMatrix_c A, ElInt k );
ElError ElWilkinson_z( ElMatrix_z A, ElInt k );

ElError ElWilkinsonDist_i( ElDistMatrix_i A, ElInt k );
ElError ElWilkinsonDist_s( ElDistMatrix_s A, ElInt k );
ElError ElWilkinsonDist_d( ElDistMatrix_d A, ElInt k );
ElError ElWilkinsonDist_c( ElDistMatrix_c A, ElInt k );
ElError ElWilkinsonDist_z( ElDistMatrix_z A, ElInt k );

/* Zeros
   ===== */
ElError ElZeros_i( ElMatrix_i A, ElInt m, ElInt n );
ElError ElZeros_s( ElMatrix_s A, ElInt m, ElInt n );
ElError ElZeros_d( ElMatrix_d A, ElInt m, ElInt n );
ElError ElZeros_c( ElMatrix_c A, ElInt m, ElInt n );
ElError ElZeros_z( ElMatrix_z A, ElInt m, ElInt n );

ElError ElZerosDist_i( ElDistMatrix_i A, ElInt m, ElInt n );
ElError ElZerosDist_s( ElDistMatrix_s A, ElInt m, ElInt n );
ElError ElZerosDist_d( ElDistMatrix_d A, ElInt m, ElInt n );
ElError ElZerosDist_c( ElDistMatrix_c A, ElInt m, ElInt n );
ElError ElZerosDist_z( ElDistMatrix_z A, ElInt m, ElInt n );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_MATRICES_C_H */
