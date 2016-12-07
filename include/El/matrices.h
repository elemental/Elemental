/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_MATRICES_C_H
#define EL_MATRICES_C_H

#ifdef __cplusplus
extern "C" {
#endif

/* Deterministic
   ############# */

/* Bull's head
   =========== */
EL_EXPORT ElError ElBullsHead_c( ElMatrix_c A, ElInt n );
EL_EXPORT ElError ElBullsHead_z( ElMatrix_z A, ElInt n );

EL_EXPORT ElError ElBullsHeadDist_c( ElDistMatrix_c A, ElInt n );
EL_EXPORT ElError ElBullsHeadDist_z( ElDistMatrix_z A, ElInt n );

/* Cauchy
   ====== */
EL_EXPORT ElError ElCauchy_s
( ElMatrix_s A, ElInt xSize, float* xBuf,
                ElInt ySize, float* yBuf ); 
EL_EXPORT ElError ElCauchy_d
( ElMatrix_d A, ElInt xSize, double* xBuf,
                ElInt ySize, double* yBuf ); 
EL_EXPORT ElError ElCauchy_c
( ElMatrix_c A, ElInt xSize, complex_float* xBuf,
                ElInt ySize, complex_float* yBuf ); 
EL_EXPORT ElError ElCauchy_z
( ElMatrix_z A, ElInt xSize, complex_double* xBuf,
                ElInt ySize, complex_double* yBuf ); 

EL_EXPORT ElError ElCauchyDist_s
( ElDistMatrix_s A, ElInt xSize, float* xBuf,
                    ElInt ySize, float* yBuf ); 
EL_EXPORT ElError ElCauchyDist_d
( ElDistMatrix_d A, ElInt xSize, double* xBuf,
                    ElInt ySize, double* yBuf ); 
EL_EXPORT ElError ElCauchyDist_c
( ElDistMatrix_c A, ElInt xSize, complex_float* xBuf,
                    ElInt ySize, complex_float* yBuf ); 
EL_EXPORT ElError ElCauchyDist_z
( ElDistMatrix_z A, ElInt xSize, complex_double* xBuf,
                    ElInt ySize, complex_double* yBuf ); 

/* Cauchy-like
   =========== */
EL_EXPORT ElError ElCauchyLike_s
( ElMatrix_s A, ElInt rSize, float* rBuf,
                ElInt sSize, float* sBuf,
                ElInt xSize, float* xBuf,
                ElInt ySize, float* yBuf ); 
EL_EXPORT ElError ElCauchyLike_d
( ElMatrix_d A, ElInt rSize, double* rBuf,
                ElInt sSize, double* sBuf,
                ElInt xSize, double* xBuf,
                ElInt ySize, double* yBuf ); 
EL_EXPORT ElError ElCauchyLike_c
( ElMatrix_c A, ElInt rSize, complex_float* rBuf,
                ElInt sSize, complex_float* sBuf,
                ElInt xSize, complex_float* xBuf,
                ElInt ySize, complex_float* yBuf ); 
EL_EXPORT ElError ElCauchyLike_z
( ElMatrix_z A, ElInt rSize, complex_double* rBuf,
                ElInt sSize, complex_double* sBuf,
                ElInt xSize, complex_double* xBuf,
                ElInt ySize, complex_double* yBuf ); 

EL_EXPORT ElError ElCauchyLikeDist_s
( ElDistMatrix_s A, ElInt rSize, float* rBuf,
                    ElInt sSize, float* sBuf,
                    ElInt xSize, float* xBuf,
                    ElInt ySize, float* yBuf ); 
EL_EXPORT ElError ElCauchyLikeDist_d
( ElDistMatrix_d A, ElInt rSize, double* rBuf,
                    ElInt sSize, double* sBuf,
                    ElInt xSize, double* xBuf,
                    ElInt ySize, double* yBuf ); 
EL_EXPORT ElError ElCauchyLikeDist_c
( ElDistMatrix_c A, ElInt rSize, complex_float* rBuf,
                    ElInt sSize, complex_float* sBuf,
                    ElInt xSize, complex_float* xBuf,
                    ElInt ySize, complex_float* yBuf ); 
EL_EXPORT ElError ElCauchyLikeDist_z
( ElDistMatrix_z A, ElInt rSize, complex_double* rBuf,
                    ElInt sSize, complex_double* sBuf,
                    ElInt xSize, complex_double* xBuf,
                    ElInt ySize, complex_double* yBuf ); 

/* Circulant
   ========= */
EL_EXPORT ElError ElCirculant_i
( ElMatrix_i A, ElInt aSize, ElInt* aBuf );
EL_EXPORT ElError ElCirculant_s
( ElMatrix_s A, ElInt aSize, float* aBuf );
EL_EXPORT ElError ElCirculant_d
( ElMatrix_d A, ElInt aSize, double* aBuf );
EL_EXPORT ElError ElCirculant_c
( ElMatrix_c A, ElInt aSize, complex_float* aBuf );
EL_EXPORT ElError ElCirculant_z
( ElMatrix_z A, ElInt aSize, complex_double* aBuf );

EL_EXPORT ElError ElCirculantDist_i
( ElDistMatrix_i A, ElInt aSize, ElInt* aBuf );
EL_EXPORT ElError ElCirculantDist_s
( ElDistMatrix_s A, ElInt aSize, float* aBuf );
EL_EXPORT ElError ElCirculantDist_d
( ElDistMatrix_d A, ElInt aSize, double* aBuf );
EL_EXPORT ElError ElCirculantDist_c
( ElDistMatrix_c A, ElInt aSize, complex_float* aBuf );
EL_EXPORT ElError ElCirculantDist_z
( ElDistMatrix_z A, ElInt aSize, complex_double* aBuf );

/* Demmel
   ====== */
EL_EXPORT ElError ElDemmel_s( ElMatrix_s A, ElInt n );
EL_EXPORT ElError ElDemmel_d( ElMatrix_d A, ElInt n );
EL_EXPORT ElError ElDemmel_c( ElMatrix_c A, ElInt n );
EL_EXPORT ElError ElDemmel_z( ElMatrix_z A, ElInt n );

EL_EXPORT ElError ElDemmelDist_s( ElDistMatrix_s A, ElInt n );
EL_EXPORT ElError ElDemmelDist_d( ElDistMatrix_d A, ElInt n );
EL_EXPORT ElError ElDemmelDist_c( ElDistMatrix_c A, ElInt n );
EL_EXPORT ElError ElDemmelDist_z( ElDistMatrix_z A, ElInt n );

/* Diagonal
   ======== */
EL_EXPORT ElError ElDiagonal_i( ElMatrix_i A, ElMatrix_i d);
EL_EXPORT ElError ElDiagonal_s( ElMatrix_s A, ElMatrix_s d );
EL_EXPORT ElError ElDiagonal_d( ElMatrix_d A, ElMatrix_d d );
EL_EXPORT ElError ElDiagonal_c( ElMatrix_c A, ElMatrix_c d );
EL_EXPORT ElError ElDiagonal_z( ElMatrix_z A, ElMatrix_z d );

EL_EXPORT ElError ElDiagonalDist_i( ElDistMatrix_i A, ElDistMatrix_i d );
EL_EXPORT ElError ElDiagonalDist_s( ElDistMatrix_s A, ElDistMatrix_s d );
EL_EXPORT ElError ElDiagonalDist_d( ElDistMatrix_d A, ElDistMatrix_d d );
EL_EXPORT ElError ElDiagonalDist_c( ElDistMatrix_c A, ElDistMatrix_c d );
EL_EXPORT ElError ElDiagonalDist_z( ElDistMatrix_z A, ElDistMatrix_z d );

EL_EXPORT ElError ElDiagonalSparse_i( ElSparseMatrix_i A, ElMatrix_i d);
EL_EXPORT ElError ElDiagonalSparse_s( ElSparseMatrix_s A, ElMatrix_s d );
EL_EXPORT ElError ElDiagonalSparse_d( ElSparseMatrix_d A, ElMatrix_d d );
EL_EXPORT ElError ElDiagonalSparse_c( ElSparseMatrix_c A, ElMatrix_c d );
EL_EXPORT ElError ElDiagonalSparse_z( ElSparseMatrix_z A, ElMatrix_z d );

EL_EXPORT ElError ElDiagonalDistSparse_i
( ElDistSparseMatrix_i A, ElDistMultiVec_i d);
EL_EXPORT ElError ElDiagonalDistSparse_s
( ElDistSparseMatrix_s A, ElDistMultiVec_s d );
EL_EXPORT ElError ElDiagonalDistSparse_d
( ElDistSparseMatrix_d A, ElDistMultiVec_d d );
EL_EXPORT ElError ElDiagonalDistSparse_c
( ElDistSparseMatrix_c A, ElDistMultiVec_c d );
EL_EXPORT ElError ElDiagonalDistSparse_z
( ElDistSparseMatrix_z A, ElDistMultiVec_z d );

/* Druinsky-Toledo
   =============== */
EL_EXPORT ElError ElDruinskyToledo_s( ElMatrix_s A, ElInt n );
EL_EXPORT ElError ElDruinskyToledo_d( ElMatrix_d A, ElInt n );
EL_EXPORT ElError ElDruinskyToledo_c( ElMatrix_c A, ElInt n );
EL_EXPORT ElError ElDruinskyToledo_z( ElMatrix_z A, ElInt n );

EL_EXPORT ElError ElDruinskyToledoDist_s( ElDistMatrix_s A, ElInt n );
EL_EXPORT ElError ElDruinskyToledoDist_d( ElDistMatrix_d A, ElInt n );
EL_EXPORT ElError ElDruinskyToledoDist_c( ElDistMatrix_c A, ElInt n );
EL_EXPORT ElError ElDruinskyToledoDist_z( ElDistMatrix_z A, ElInt n );

/* Dynamic regularization counter-example
   ====================================== */
EL_EXPORT ElError ElDynamicRegCounter_s( ElMatrix_s A, ElInt n );
EL_EXPORT ElError ElDynamicRegCounter_d( ElMatrix_d A, ElInt n );
EL_EXPORT ElError ElDynamicRegCounter_c( ElMatrix_c A, ElInt n );
EL_EXPORT ElError ElDynamicRegCounter_z( ElMatrix_z A, ElInt n );

EL_EXPORT ElError ElDynamicRegCounterDist_s( ElDistMatrix_s A, ElInt n );
EL_EXPORT ElError ElDynamicRegCounterDist_d( ElDistMatrix_d A, ElInt n );
EL_EXPORT ElError ElDynamicRegCounterDist_c( ElDistMatrix_c A, ElInt n );
EL_EXPORT ElError ElDynamicRegCounterDist_z( ElDistMatrix_z A, ElInt n );

EL_EXPORT ElError ElDynamicRegCounterSparse_s( ElSparseMatrix_s A, ElInt n );
EL_EXPORT ElError ElDynamicRegCounterSparse_d( ElSparseMatrix_d A, ElInt n );
EL_EXPORT ElError ElDynamicRegCounterSparse_c( ElSparseMatrix_c A, ElInt n );
EL_EXPORT ElError ElDynamicRegCounterSparse_z( ElSparseMatrix_z A, ElInt n );

EL_EXPORT ElError ElDynamicRegCounterDistSparse_s
( ElDistSparseMatrix_s A, ElInt n );
EL_EXPORT ElError ElDynamicRegCounterDistSparse_d
( ElDistSparseMatrix_d A, ElInt n );
EL_EXPORT ElError ElDynamicRegCounterDistSparse_c
( ElDistSparseMatrix_c A, ElInt n );
EL_EXPORT ElError ElDynamicRegCounterDistSparse_z
( ElDistSparseMatrix_z A, ElInt n );

/* Egorov
   ====== */
EL_EXPORT ElError ElEgorov_c
( ElMatrix_c A, float (*phase)( ElInt, ElInt ), ElInt n );
EL_EXPORT ElError ElEgorov_z
( ElMatrix_z A, double (*phase)( ElInt, ElInt ), ElInt n );

EL_EXPORT ElError ElEgorovDist_c
( ElDistMatrix_c A, float (*phase)( ElInt, ElInt ), ElInt n );
EL_EXPORT ElError ElEgorovDist_z
( ElDistMatrix_z A, double (*phase)( ElInt, ElInt ), ElInt n );

/* Ehrenfest
   ========= */
EL_EXPORT ElError ElEhrenfest_s( ElMatrix_s P, ElInt n );
EL_EXPORT ElError ElEhrenfest_d( ElMatrix_d P, ElInt n );
EL_EXPORT ElError ElEhrenfest_c( ElMatrix_c P, ElInt n );
EL_EXPORT ElError ElEhrenfest_z( ElMatrix_z P, ElInt n );

EL_EXPORT ElError ElEhrenfestDist_s( ElDistMatrix_s P, ElInt n );
EL_EXPORT ElError ElEhrenfestDist_d( ElDistMatrix_d P, ElInt n );
EL_EXPORT ElError ElEhrenfestDist_c( ElDistMatrix_c P, ElInt n );
EL_EXPORT ElError ElEhrenfestDist_z( ElDistMatrix_z P, ElInt n );

EL_EXPORT ElError ElEhrenfestStationary_s( ElMatrix_s PInf, ElInt n );
EL_EXPORT ElError ElEhrenfestStationary_d( ElMatrix_d PInf, ElInt n );
EL_EXPORT ElError ElEhrenfestStationary_c( ElMatrix_c PInf, ElInt n );
EL_EXPORT ElError ElEhrenfestStationary_z( ElMatrix_z PInf, ElInt n );

EL_EXPORT ElError ElEhrenfestStationaryDist_s( ElDistMatrix_s PInf, ElInt n );
EL_EXPORT ElError ElEhrenfestStationaryDist_d( ElDistMatrix_d PInf, ElInt n );
EL_EXPORT ElError ElEhrenfestStationaryDist_c( ElDistMatrix_c PInf, ElInt n );
EL_EXPORT ElError ElEhrenfestStationaryDist_z( ElDistMatrix_z PInf, ElInt n );

EL_EXPORT ElError ElEhrenfestDecay_s( ElMatrix_s A, ElInt n );
EL_EXPORT ElError ElEhrenfestDecay_d( ElMatrix_d A, ElInt n );
EL_EXPORT ElError ElEhrenfestDecay_c( ElMatrix_c A, ElInt n );
EL_EXPORT ElError ElEhrenfestDecay_z( ElMatrix_z A, ElInt n );

EL_EXPORT ElError ElEhrenfestDecayDist_s( ElDistMatrix_s A, ElInt n );
EL_EXPORT ElError ElEhrenfestDecayDist_d( ElDistMatrix_d A, ElInt n );
EL_EXPORT ElError ElEhrenfestDecayDist_c( ElDistMatrix_c A, ElInt n );
EL_EXPORT ElError ElEhrenfestDecayDist_z( ElDistMatrix_z A, ElInt n );

/* ExtendedKahan
   ============= */
EL_EXPORT ElError ElExtendedKahan_s
( ElMatrix_s A, ElInt k, float phi, float mu );
EL_EXPORT ElError ElExtendedKahan_d
( ElMatrix_d A, ElInt k, double phi, double mu );
EL_EXPORT ElError ElExtendedKahan_c
( ElMatrix_c A, ElInt k, float phi, float mu );
EL_EXPORT ElError ElExtendedKahan_z
( ElMatrix_z A, ElInt k, double phi, double mu );

EL_EXPORT ElError ElExtendedKahanDist_s
( ElDistMatrix_s A, ElInt k, float phi, float mu );
EL_EXPORT ElError ElExtendedKahanDist_d
( ElDistMatrix_d A, ElInt k, double phi, double mu );
EL_EXPORT ElError ElExtendedKahanDist_c
( ElDistMatrix_c A, ElInt k, float phi, float mu );
EL_EXPORT ElError ElExtendedKahanDist_z
( ElDistMatrix_z A, ElInt k, double phi, double mu );

/* Fiedler
   ======= */
EL_EXPORT ElError ElFiedler_s
( ElMatrix_s A, ElInt cSize, float* cBuf );
EL_EXPORT ElError ElFiedler_d
( ElMatrix_d A, ElInt cSize, double* cBuf );
EL_EXPORT ElError ElFiedler_c
( ElMatrix_c A, ElInt cSize, complex_float* cBuf );
EL_EXPORT ElError ElFiedler_z
( ElMatrix_z A, ElInt cSize, complex_double* cBuf );

EL_EXPORT ElError ElFiedlerDist_s
( ElDistMatrix_s A, ElInt cSize, float* cBuf );
EL_EXPORT ElError ElFiedlerDist_d
( ElDistMatrix_d A, ElInt cSize, double* cBuf );
EL_EXPORT ElError ElFiedlerDist_c
( ElDistMatrix_c A, ElInt cSize, complex_float* cBuf );
EL_EXPORT ElError ElFiedlerDist_z
( ElDistMatrix_z A, ElInt cSize, complex_double* cBuf );

/* Forsythe
   ======== */
EL_EXPORT ElError ElForsythe_i
( ElMatrix_i J, ElInt n, ElInt alpha, ElInt lambda );
EL_EXPORT ElError ElForsythe_s
( ElMatrix_s J, ElInt n, float alpha, float lambda );
EL_EXPORT ElError ElForsythe_d
( ElMatrix_d J, ElInt n, double alpha, double lambda );
EL_EXPORT ElError ElForsythe_c
( ElMatrix_c J, ElInt n, complex_float alpha, complex_float lambda );
EL_EXPORT ElError ElForsythe_z
( ElMatrix_z J, ElInt n, complex_double alpha, complex_double lambda );

EL_EXPORT ElError ElForsytheDist_i
( ElDistMatrix_i J, ElInt n, ElInt alpha, ElInt lambda );
EL_EXPORT ElError ElForsytheDist_s
( ElDistMatrix_s J, ElInt n, float alpha, float lambda );
EL_EXPORT ElError ElForsytheDist_d
( ElDistMatrix_d J, ElInt n, double alpha, double lambda );
EL_EXPORT ElError ElForsytheDist_c
( ElDistMatrix_c J, ElInt n, complex_float alpha, complex_float lambda );
EL_EXPORT ElError ElForsytheDist_z
( ElDistMatrix_z J, ElInt n, complex_double alpha, complex_double lambda );

/* Fox-Li
   ====== */
EL_EXPORT ElError ElFoxLi_c( ElMatrix_c A, ElInt n, float omega );
EL_EXPORT ElError ElFoxLi_z( ElMatrix_z A, ElInt n, double omega );

EL_EXPORT ElError ElFoxLiDist_c( ElDistMatrix_c A, ElInt n, float omega );
EL_EXPORT ElError ElFoxLiDist_z( ElDistMatrix_z A, ElInt n, double omega );

/* Fourier
   ======= */
EL_EXPORT ElError ElFourier_c( ElMatrix_c A, ElInt n );
EL_EXPORT ElError ElFourier_z( ElMatrix_z A, ElInt n );

EL_EXPORT ElError ElFourierDist_c( ElDistMatrix_c A, ElInt n );
EL_EXPORT ElError ElFourierDist_z( ElDistMatrix_z A, ElInt n );

/* GCD matrix
   ========== */
EL_EXPORT ElError ElGCDMatrix_i( ElMatrix_i G, ElInt m, ElInt n );
EL_EXPORT ElError ElGCDMatrix_s( ElMatrix_s G, ElInt m, ElInt n );
EL_EXPORT ElError ElGCDMatrix_d( ElMatrix_d G, ElInt m, ElInt n );
EL_EXPORT ElError ElGCDMatrix_c( ElMatrix_c G, ElInt m, ElInt n );
EL_EXPORT ElError ElGCDMatrix_z( ElMatrix_z G, ElInt m, ElInt n );

EL_EXPORT ElError ElGCDMatrixDist_i( ElDistMatrix_i G, ElInt m, ElInt n );
EL_EXPORT ElError ElGCDMatrixDist_s( ElDistMatrix_s G, ElInt m, ElInt n );
EL_EXPORT ElError ElGCDMatrixDist_d( ElDistMatrix_d G, ElInt m, ElInt n );
EL_EXPORT ElError ElGCDMatrixDist_c( ElDistMatrix_c G, ElInt m, ElInt n );
EL_EXPORT ElError ElGCDMatrixDist_z( ElDistMatrix_z G, ElInt m, ElInt n );

/* Gear matrix
   =========== */
EL_EXPORT ElError ElGear_i( ElMatrix_i G, ElInt n, ElInt s, ElInt t );
EL_EXPORT ElError ElGear_s( ElMatrix_s G, ElInt n, ElInt s, ElInt t );
EL_EXPORT ElError ElGear_d( ElMatrix_d G, ElInt n, ElInt s, ElInt t );
EL_EXPORT ElError ElGear_c( ElMatrix_c G, ElInt n, ElInt s, ElInt t );
EL_EXPORT ElError ElGear_z( ElMatrix_z G, ElInt n, ElInt s, ElInt t );

EL_EXPORT ElError ElGearDist_i( ElDistMatrix_i G, ElInt n, ElInt s, ElInt t );
EL_EXPORT ElError ElGearDist_s( ElDistMatrix_s G, ElInt n, ElInt s, ElInt t );
EL_EXPORT ElError ElGearDist_d( ElDistMatrix_d G, ElInt n, ElInt s, ElInt t );
EL_EXPORT ElError ElGearDist_c( ElDistMatrix_c G, ElInt n, ElInt s, ElInt t );
EL_EXPORT ElError ElGearDist_z( ElDistMatrix_z G, ElInt n, ElInt s, ElInt t );

/* GEPP Growth
   =========== */
EL_EXPORT ElError ElGEPPGrowth_s( ElMatrix_s A, ElInt n );
EL_EXPORT ElError ElGEPPGrowth_d( ElMatrix_d A, ElInt n );
EL_EXPORT ElError ElGEPPGrowth_c( ElMatrix_c A, ElInt n );
EL_EXPORT ElError ElGEPPGrowth_z( ElMatrix_z A, ElInt n );

EL_EXPORT ElError ElGEPPGrowthDist_s( ElDistMatrix_s A, ElInt n );
EL_EXPORT ElError ElGEPPGrowthDist_d( ElDistMatrix_d A, ElInt n );
EL_EXPORT ElError ElGEPPGrowthDist_c( ElDistMatrix_c A, ElInt n );
EL_EXPORT ElError ElGEPPGrowthDist_z( ElDistMatrix_z A, ElInt n );

/* Golub Klema Stewart
   =================== */
EL_EXPORT ElError ElGKS_s( ElMatrix_s A, ElInt n );
EL_EXPORT ElError ElGKS_d( ElMatrix_d A, ElInt n );
EL_EXPORT ElError ElGKS_c( ElMatrix_c A, ElInt n );
EL_EXPORT ElError ElGKS_z( ElMatrix_z A, ElInt n );

EL_EXPORT ElError ElGKSDist_s( ElDistMatrix_s A, ElInt n );
EL_EXPORT ElError ElGKSDist_d( ElDistMatrix_d A, ElInt n );
EL_EXPORT ElError ElGKSDist_c( ElDistMatrix_c A, ElInt n );
EL_EXPORT ElError ElGKSDist_z( ElDistMatrix_z A, ElInt n );

/* Grcar
   ===== */
EL_EXPORT ElError ElGrcar_i( ElMatrix_i A, ElInt n, ElInt k );
EL_EXPORT ElError ElGrcar_s( ElMatrix_s A, ElInt n, ElInt k );
EL_EXPORT ElError ElGrcar_d( ElMatrix_d A, ElInt n, ElInt k );
EL_EXPORT ElError ElGrcar_c( ElMatrix_c A, ElInt n, ElInt k );
EL_EXPORT ElError ElGrcar_z( ElMatrix_z A, ElInt n, ElInt k );

EL_EXPORT ElError ElGrcarDist_i( ElDistMatrix_i A, ElInt n, ElInt k );
EL_EXPORT ElError ElGrcarDist_s( ElDistMatrix_s A, ElInt n, ElInt k );
EL_EXPORT ElError ElGrcarDist_d( ElDistMatrix_d A, ElInt n, ElInt k );
EL_EXPORT ElError ElGrcarDist_c( ElDistMatrix_c A, ElInt n, ElInt k );
EL_EXPORT ElError ElGrcarDist_z( ElDistMatrix_z A, ElInt n, ElInt k );

/* Hankel
   ====== */
EL_EXPORT ElError ElHankel_i
( ElMatrix_i A, ElInt m, ElInt n, ElInt aSize, ElInt* aBuf );
EL_EXPORT ElError ElHankel_s
( ElMatrix_s A, ElInt m, ElInt n, ElInt aSize, float* aBuf );
EL_EXPORT ElError ElHankel_d
( ElMatrix_d A, ElInt m, ElInt n, ElInt aSize, double* aBuf );
EL_EXPORT ElError ElHankel_c
( ElMatrix_c A, ElInt m, ElInt n, ElInt aSize, complex_float* aBuf );
EL_EXPORT ElError ElHankel_z
( ElMatrix_z A, ElInt m, ElInt n, ElInt aSize, complex_double* aBuf );

EL_EXPORT ElError ElHankelDist_i
( ElDistMatrix_i A, ElInt m, ElInt n, ElInt aSize, ElInt* aBuf );
EL_EXPORT ElError ElHankelDist_s
( ElDistMatrix_s A, ElInt m, ElInt n, ElInt aSize, float* aBuf );
EL_EXPORT ElError ElHankelDist_d
( ElDistMatrix_d A, ElInt m, ElInt n, ElInt aSize, double* aBuf );
EL_EXPORT ElError ElHankelDist_c
( ElDistMatrix_c A, ElInt m, ElInt n, ElInt aSize, complex_float* aBuf );
EL_EXPORT ElError ElHankelDist_z
( ElDistMatrix_z A, ElInt m, ElInt n, ElInt aSize, complex_double* aBuf );

/* Hanowa
   ====== */
EL_EXPORT ElError ElHanowa_i( ElMatrix_i A, ElInt n, ElInt mu );
EL_EXPORT ElError ElHanowa_s( ElMatrix_s A, ElInt n, float mu );
EL_EXPORT ElError ElHanowa_d( ElMatrix_d A, ElInt n, double mu );
EL_EXPORT ElError ElHanowa_c( ElMatrix_c A, ElInt n, complex_float mu );
EL_EXPORT ElError ElHanowa_z( ElMatrix_z A, ElInt n, complex_double mu );

EL_EXPORT ElError ElHanowaDist_i( ElDistMatrix_i A, ElInt n, ElInt mu );
EL_EXPORT ElError ElHanowaDist_s( ElDistMatrix_s A, ElInt n, float mu );
EL_EXPORT ElError ElHanowaDist_d( ElDistMatrix_d A, ElInt n, double mu );
EL_EXPORT ElError ElHanowaDist_c( ElDistMatrix_c A, ElInt n, complex_float mu );
EL_EXPORT ElError ElHanowaDist_z( ElDistMatrix_z A, ElInt n, complex_double mu );
/* Helmholtz
   ========= */
EL_EXPORT ElError ElHelmholtz1D_s
( ElMatrix_s H, ElInt nx, float shift );
EL_EXPORT ElError ElHelmholtz1D_d
( ElMatrix_d H, ElInt nx, double shift );
EL_EXPORT ElError ElHelmholtz1D_c
( ElMatrix_c H, ElInt nx, complex_float shift );
EL_EXPORT ElError ElHelmholtz1D_z
( ElMatrix_z H, ElInt nx, complex_double shift );

EL_EXPORT ElError ElHelmholtz1DDist_s
( ElDistMatrix_s H, ElInt nx, float shift );
EL_EXPORT ElError ElHelmholtz1DDist_d
( ElDistMatrix_d H, ElInt nx, double shift );
EL_EXPORT ElError ElHelmholtz1DDist_c
( ElDistMatrix_c H, ElInt nx, complex_float shift );
EL_EXPORT ElError ElHelmholtz1DDist_z
( ElDistMatrix_z H, ElInt nx, complex_double shift );

EL_EXPORT ElError ElHelmholtz1DSparse_s
( ElSparseMatrix_s H, ElInt nx, float shift );
EL_EXPORT ElError ElHelmholtz1DSparse_d
( ElSparseMatrix_d H, ElInt nx, double shift );
EL_EXPORT ElError ElHelmholtz1DSparse_c
( ElSparseMatrix_c H, ElInt nx, complex_float shift );
EL_EXPORT ElError ElHelmholtz1DSparse_z
( ElSparseMatrix_z H, ElInt nx, complex_double shift );

EL_EXPORT ElError ElHelmholtz1DDistSparse_s
( ElDistSparseMatrix_s H, ElInt nx, float shift );
EL_EXPORT ElError ElHelmholtz1DDistSparse_d
( ElDistSparseMatrix_d H, ElInt nx, double shift );
EL_EXPORT ElError ElHelmholtz1DDistSparse_c
( ElDistSparseMatrix_c H, ElInt nx, complex_float shift );
EL_EXPORT ElError ElHelmholtz1DDistSparse_z
( ElDistSparseMatrix_z H, ElInt nx, complex_double shift );

EL_EXPORT ElError ElHelmholtz2D_s
( ElMatrix_s H, ElInt nx, ElInt ny, float shift );
EL_EXPORT ElError ElHelmholtz2D_d
( ElMatrix_d H, ElInt nx, ElInt ny, double shift );
EL_EXPORT ElError ElHelmholtz2D_c
( ElMatrix_c H, ElInt nx, ElInt ny, complex_float shift );
EL_EXPORT ElError ElHelmholtz2D_z
( ElMatrix_z H, ElInt nx, ElInt ny, complex_double shift );

EL_EXPORT ElError ElHelmholtz2DDist_s
( ElDistMatrix_s H, ElInt nx, ElInt ny, float shift );
EL_EXPORT ElError ElHelmholtz2DDist_d
( ElDistMatrix_d H, ElInt nx, ElInt ny, double shift );
EL_EXPORT ElError ElHelmholtz2DDist_c
( ElDistMatrix_c H, ElInt nx, ElInt ny, complex_float shift );
EL_EXPORT ElError ElHelmholtz2DDist_z
( ElDistMatrix_z H, ElInt nx, ElInt ny, complex_double shift );

EL_EXPORT ElError ElHelmholtz2DSparse_s
( ElSparseMatrix_s H, ElInt nx, ElInt ny, float shift );
EL_EXPORT ElError ElHelmholtz2DSparse_d
( ElSparseMatrix_d H, ElInt nx, ElInt ny, double shift );
EL_EXPORT ElError ElHelmholtz2DSparse_c
( ElSparseMatrix_c H, ElInt nx, ElInt ny, complex_float shift );
EL_EXPORT ElError ElHelmholtz2DSparse_z
( ElSparseMatrix_z H, ElInt nx, ElInt ny, complex_double shift );

EL_EXPORT ElError ElHelmholtz2DDistSparse_s
( ElDistSparseMatrix_s H, ElInt nx, ElInt ny, float shift );
EL_EXPORT ElError ElHelmholtz2DDistSparse_d
( ElDistSparseMatrix_d H, ElInt nx, ElInt ny, double shift );
EL_EXPORT ElError ElHelmholtz2DDistSparse_c
( ElDistSparseMatrix_c H, ElInt nx, ElInt ny, complex_float shift );
EL_EXPORT ElError ElHelmholtz2DDistSparse_z
( ElDistSparseMatrix_z H, ElInt nx, ElInt ny, complex_double shift );

EL_EXPORT ElError ElHelmholtz3D_s
( ElMatrix_s H, ElInt nx, ElInt ny, ElInt nz, float shift );
EL_EXPORT ElError ElHelmholtz3D_d
( ElMatrix_d H, ElInt nx, ElInt ny, ElInt nz, double shift );
EL_EXPORT ElError ElHelmholtz3D_c
( ElMatrix_c H, ElInt nx, ElInt ny, ElInt nz, complex_float shift );
EL_EXPORT ElError ElHelmholtz3D_z
( ElMatrix_z H, ElInt nx, ElInt ny, ElInt nz, complex_double shift );

EL_EXPORT ElError ElHelmholtz3DDist_s
( ElDistMatrix_s H, ElInt nx, ElInt ny, ElInt nz, float shift );
EL_EXPORT ElError ElHelmholtz3DDist_d
( ElDistMatrix_d H, ElInt nx, ElInt ny, ElInt nz, double shift );
EL_EXPORT ElError ElHelmholtz3DDist_c
( ElDistMatrix_c H, ElInt nx, ElInt ny, ElInt nz, complex_float shift );
EL_EXPORT ElError ElHelmholtz3DDist_z
( ElDistMatrix_z H, ElInt nx, ElInt ny, ElInt nz, complex_double shift );

EL_EXPORT ElError ElHelmholtz3DSparse_s
( ElSparseMatrix_s H, ElInt nx, ElInt ny, ElInt nz, float shift );
EL_EXPORT ElError ElHelmholtz3DSparse_d
( ElSparseMatrix_d H, ElInt nx, ElInt ny, ElInt nz, double shift );
EL_EXPORT ElError ElHelmholtz3DSparse_c
( ElSparseMatrix_c H, ElInt nx, ElInt ny, ElInt nz, complex_float shift );
EL_EXPORT ElError ElHelmholtz3DSparse_z
( ElSparseMatrix_z H, ElInt nx, ElInt ny, ElInt nz, complex_double shift );

EL_EXPORT ElError ElHelmholtz3DDistSparse_s
( ElDistSparseMatrix_s H, ElInt nx, ElInt ny, ElInt nz, float shift );
EL_EXPORT ElError ElHelmholtz3DDistSparse_d
( ElDistSparseMatrix_d H, ElInt nx, ElInt ny, ElInt nz, double shift );
EL_EXPORT ElError ElHelmholtz3DDistSparse_c
( ElDistSparseMatrix_c H, ElInt nx, ElInt ny, ElInt nz, complex_float shift );
EL_EXPORT ElError ElHelmholtz3DDistSparse_z
( ElDistSparseMatrix_z H, ElInt nx, ElInt ny, ElInt nz, complex_double shift );

/* Helmholtz with PML
   ================== */
EL_EXPORT ElError ElHelmholtzPML1D_c
( ElMatrix_c H, ElInt nx, complex_float omega, 
  ElInt numPmlPoints, float sigma, float pmlExp );
EL_EXPORT ElError ElHelmholtzPML1D_z
( ElMatrix_z H, ElInt nx, complex_double omega, 
  ElInt numPmlPoints, double sigma, double pmlExp );

EL_EXPORT ElError ElHelmholtzPML1DDist_c
( ElDistMatrix_c H, ElInt nx, complex_float omega, 
  ElInt numPmlPoints, float sigma, float pmlExp );
EL_EXPORT ElError ElHelmholtzPML1DDist_z
( ElDistMatrix_z H, ElInt nx, complex_double omega, 
  ElInt numPmlPoints, double sigma, double pmlExp );

EL_EXPORT ElError ElHelmholtzPML1DSparse_c
( ElSparseMatrix_c H, ElInt nx, complex_float omega, 
  ElInt numPmlPoints, float sigma, float pmlExp );
EL_EXPORT ElError ElHelmholtzPML1DSparse_z
( ElSparseMatrix_z H, ElInt nx, complex_double omega, 
  ElInt numPmlPoints, double sigma, double pmlExp );

EL_EXPORT ElError ElHelmholtzPML1DDistSparse_c
( ElDistSparseMatrix_c H, ElInt nx, complex_float omega, 
  ElInt numPmlPoints, float sigma, float pmlExp );
EL_EXPORT ElError ElHelmholtzPML1DDistSparse_z
( ElDistSparseMatrix_z H, ElInt nx, complex_double omega, 
  ElInt numPmlPoints, double sigma, double pmlExp );

EL_EXPORT ElError ElHelmholtzPML2D_c
( ElMatrix_c H, ElInt nx, ElInt ny, complex_float omega, 
  ElInt numPmlPoints, float sigma, float pmlExp );
EL_EXPORT ElError ElHelmholtzPML2D_z
( ElMatrix_z H, ElInt nx, ElInt ny, complex_double omega, 
  ElInt numPmlPoints, double sigma, double pmlExp );

EL_EXPORT ElError ElHelmholtzPML2DDist_c
( ElDistMatrix_c H, ElInt nx, ElInt ny, complex_float omega, 
  ElInt numPmlPoints, float sigma, float pmlExp );
EL_EXPORT ElError ElHelmholtzPML2DDist_z
( ElDistMatrix_z H, ElInt nx, ElInt ny, complex_double omega, 
  ElInt numPmlPoints, double sigma, double pmlExp );

EL_EXPORT ElError ElHelmholtzPML2DSparse_c
( ElSparseMatrix_c H, ElInt nx, ElInt ny, complex_float omega, 
  ElInt numPmlPoints, float sigma, float pmlExp );
EL_EXPORT ElError ElHelmholtzPML2DSparse_z
( ElSparseMatrix_z H, ElInt nx, ElInt ny, complex_double omega, 
  ElInt numPmlPoints, double sigma, double pmlExp );

EL_EXPORT ElError ElHelmholtzPML2DDistSparse_c
( ElDistSparseMatrix_c H, ElInt nx, ElInt ny, complex_float omega, 
  ElInt numPmlPoints, float sigma, float pmlExp );
EL_EXPORT ElError ElHelmholtzPML2DDistSparse_z
( ElDistSparseMatrix_z H, ElInt nx, ElInt ny, complex_double omega, 
  ElInt numPmlPoints, double sigma, double pmlExp );

EL_EXPORT ElError ElHelmholtzPML3D_c
( ElMatrix_c H, ElInt nx, ElInt ny, ElInt nz, complex_float omega, 
  ElInt numPmlPoints, float sigma, float pmlExp );
EL_EXPORT ElError ElHelmholtzPML3D_z
( ElMatrix_z H, ElInt nx, ElInt ny, ElInt nz, complex_double omega, 
  ElInt numPmlPoints, double sigma, double pmlExp );

EL_EXPORT ElError ElHelmholtzPML3DDist_c
( ElDistMatrix_c H, ElInt nx, ElInt ny, ElInt nz, complex_float omega, 
  ElInt numPmlPoints, float sigma, float pmlExp );
EL_EXPORT ElError ElHelmholtzPML3DDist_z
( ElDistMatrix_z H, ElInt nx, ElInt ny, ElInt nz, complex_double omega, 
  ElInt numPmlPoints, double sigma, double pmlExp );

EL_EXPORT ElError ElHelmholtzPML3DSparse_c
( ElSparseMatrix_c H, ElInt nx, ElInt ny, ElInt nz, complex_float omega, 
  ElInt numPmlPoints, float sigma, float pmlExp );
EL_EXPORT ElError ElHelmholtzPML3DSparse_z
( ElSparseMatrix_z H, ElInt nx, ElInt ny, ElInt nz, complex_double omega, 
  ElInt numPmlPoints, double sigma, double pmlExp );

EL_EXPORT ElError ElHelmholtzPML3DDistSparse_c
( ElDistSparseMatrix_c H, ElInt nx, ElInt ny, ElInt nz, complex_float omega, 
  ElInt numPmlPoints, float sigma, float pmlExp );
EL_EXPORT ElError ElHelmholtzPML3DDistSparse_z
( ElDistSparseMatrix_z H, ElInt nx, ElInt ny, ElInt nz, complex_double omega, 
  ElInt numPmlPoints, double sigma, double pmlExp );

/* Hilbert
   ======= */
EL_EXPORT ElError ElHilbert_s( ElMatrix_s A, ElInt n );
EL_EXPORT ElError ElHilbert_d( ElMatrix_d A, ElInt n );
EL_EXPORT ElError ElHilbert_c( ElMatrix_c A, ElInt n );
EL_EXPORT ElError ElHilbert_z( ElMatrix_z A, ElInt n );

EL_EXPORT ElError ElHilbertDist_s( ElDistMatrix_s A, ElInt n );
EL_EXPORT ElError ElHilbertDist_d( ElDistMatrix_d A, ElInt n );
EL_EXPORT ElError ElHilbertDist_c( ElDistMatrix_c A, ElInt n );
EL_EXPORT ElError ElHilbertDist_z( ElDistMatrix_z A, ElInt n );

/* Identity
   ======== */
EL_EXPORT ElError ElIdentity_i( ElMatrix_i A, ElInt m, ElInt n );
EL_EXPORT ElError ElIdentity_s( ElMatrix_s A, ElInt m, ElInt n );
EL_EXPORT ElError ElIdentity_d( ElMatrix_d A, ElInt m, ElInt n );
EL_EXPORT ElError ElIdentity_c( ElMatrix_c A, ElInt m, ElInt n );
EL_EXPORT ElError ElIdentity_z( ElMatrix_z A, ElInt m, ElInt n );

EL_EXPORT ElError ElIdentityDist_i( ElDistMatrix_i A, ElInt m, ElInt n );
EL_EXPORT ElError ElIdentityDist_s( ElDistMatrix_s A, ElInt m, ElInt n );
EL_EXPORT ElError ElIdentityDist_d( ElDistMatrix_d A, ElInt m, ElInt n );
EL_EXPORT ElError ElIdentityDist_c( ElDistMatrix_c A, ElInt m, ElInt n );
EL_EXPORT ElError ElIdentityDist_z( ElDistMatrix_z A, ElInt m, ElInt n );

EL_EXPORT ElError ElIdentitySparse_i( ElSparseMatrix_i A, ElInt m, ElInt n );
EL_EXPORT ElError ElIdentitySparse_s( ElSparseMatrix_s A, ElInt m, ElInt n );
EL_EXPORT ElError ElIdentitySparse_d( ElSparseMatrix_d A, ElInt m, ElInt n );
EL_EXPORT ElError ElIdentitySparse_c( ElSparseMatrix_c A, ElInt m, ElInt n );
EL_EXPORT ElError ElIdentitySparse_z( ElSparseMatrix_z A, ElInt m, ElInt n );

EL_EXPORT ElError ElIdentityDistSparse_i
( ElDistSparseMatrix_i A, ElInt m, ElInt n );
EL_EXPORT ElError ElIdentityDistSparse_s
( ElDistSparseMatrix_s A, ElInt m, ElInt n );
EL_EXPORT ElError ElIdentityDistSparse_d
( ElDistSparseMatrix_d A, ElInt m, ElInt n );
EL_EXPORT ElError ElIdentityDistSparse_c
( ElDistSparseMatrix_c A, ElInt m, ElInt n );
EL_EXPORT ElError ElIdentityDistSparse_z
( ElDistSparseMatrix_z A, ElInt m, ElInt n );

/* Jordan
   ====== */
EL_EXPORT ElError ElJordan_i( ElMatrix_i J, ElInt n, ElInt lambda );
EL_EXPORT ElError ElJordan_s( ElMatrix_s J, ElInt n, float lambda );
EL_EXPORT ElError ElJordan_d( ElMatrix_d J, ElInt n, double lambda );
EL_EXPORT ElError ElJordan_c( ElMatrix_c J, ElInt n, complex_float lambda );
EL_EXPORT ElError ElJordan_z( ElMatrix_z J, ElInt n, complex_double lambda );

EL_EXPORT ElError ElJordanDist_i
( ElDistMatrix_i J, ElInt n, ElInt lambda );
EL_EXPORT ElError ElJordanDist_s
( ElDistMatrix_s J, ElInt n, float lambda );
EL_EXPORT ElError ElJordanDist_d
( ElDistMatrix_d J, ElInt n, double lambda );
EL_EXPORT ElError ElJordanDist_c
( ElDistMatrix_c J, ElInt n, complex_float lambda );
EL_EXPORT ElError ElJordanDist_z
( ElDistMatrix_z J, ElInt n, complex_double lambda );

/* Jordan-Cholesky (a matrix whose Cholesky factor is 2 J_{1/2}(n))
   ================================================================ */
EL_EXPORT ElError ElJordanCholesky_s( ElMatrix_s A, ElInt n );
EL_EXPORT ElError ElJordanCholesky_d( ElMatrix_d A, ElInt n );
EL_EXPORT ElError ElJordanCholesky_c( ElMatrix_c A, ElInt n );
EL_EXPORT ElError ElJordanCholesky_z( ElMatrix_z A, ElInt n );

EL_EXPORT ElError ElJordanCholeskyDist_s( ElDistMatrix_s A, ElInt n );
EL_EXPORT ElError ElJordanCholeskyDist_d( ElDistMatrix_d A, ElInt n );
EL_EXPORT ElError ElJordanCholeskyDist_c( ElDistMatrix_c A, ElInt n );
EL_EXPORT ElError ElJordanCholeskyDist_z( ElDistMatrix_z A, ElInt n );

EL_EXPORT ElError ElJordanCholeskySparse_s( ElSparseMatrix_s A, ElInt n );
EL_EXPORT ElError ElJordanCholeskySparse_d( ElSparseMatrix_d A, ElInt n );
EL_EXPORT ElError ElJordanCholeskySparse_c( ElSparseMatrix_c A, ElInt n );
EL_EXPORT ElError ElJordanCholeskySparse_z( ElSparseMatrix_z A, ElInt n );

EL_EXPORT ElError ElJordanCholeskyDistSparse_s
( ElDistSparseMatrix_s A, ElInt n );
EL_EXPORT ElError ElJordanCholeskyDistSparse_d
( ElDistSparseMatrix_d A, ElInt n );
EL_EXPORT ElError ElJordanCholeskyDistSparse_c
( ElDistSparseMatrix_c A, ElInt n );
EL_EXPORT ElError ElJordanCholeskyDistSparse_z
( ElDistSparseMatrix_z A, ElInt n );

/* Kahan
   ===== */
EL_EXPORT ElError ElKahan_s( ElMatrix_s A, ElInt n, float phi );
EL_EXPORT ElError ElKahan_d( ElMatrix_d A, ElInt n, double phi );
EL_EXPORT ElError ElKahan_c( ElMatrix_c A, ElInt n, complex_float phi );
EL_EXPORT ElError ElKahan_z( ElMatrix_z A, ElInt n, complex_double phi );

EL_EXPORT ElError ElKahanDist_s
( ElDistMatrix_s A, ElInt n, float phi );
EL_EXPORT ElError ElKahanDist_d
( ElDistMatrix_d A, ElInt n, double phi );
EL_EXPORT ElError ElKahanDist_c
( ElDistMatrix_c A, ElInt n, complex_float phi );
EL_EXPORT ElError ElKahanDist_z
( ElDistMatrix_z A, ElInt n, complex_double phi );

/* KMS
   === */
EL_EXPORT ElError ElKMS_i( ElMatrix_i K, ElInt n, ElInt rho );
EL_EXPORT ElError ElKMS_s( ElMatrix_s K, ElInt n, float rho );
EL_EXPORT ElError ElKMS_d( ElMatrix_d K, ElInt n, double rho );
EL_EXPORT ElError ElKMS_c( ElMatrix_c K, ElInt n, complex_float rho );
EL_EXPORT ElError ElKMS_z( ElMatrix_z K, ElInt n, complex_double rho );

EL_EXPORT ElError ElKMSDist_i( ElDistMatrix_i K, ElInt n, ElInt rho );
EL_EXPORT ElError ElKMSDist_s( ElDistMatrix_s K, ElInt n, float rho );
EL_EXPORT ElError ElKMSDist_d( ElDistMatrix_d K, ElInt n, double rho );
EL_EXPORT ElError ElKMSDist_c( ElDistMatrix_c K, ElInt n, complex_float rho );
EL_EXPORT ElError ElKMSDist_z( ElDistMatrix_z K, ElInt n, complex_double rho );

/* Laplacian
   ========= */
EL_EXPORT ElError ElLaplacian1D_s( ElMatrix_s L, ElInt nx );
EL_EXPORT ElError ElLaplacian1D_d( ElMatrix_d L, ElInt nx );
EL_EXPORT ElError ElLaplacian1D_c( ElMatrix_c L, ElInt nx );
EL_EXPORT ElError ElLaplacian1D_z( ElMatrix_z L, ElInt nx );

EL_EXPORT ElError ElLaplacian1DDist_s( ElDistMatrix_s L, ElInt nx );
EL_EXPORT ElError ElLaplacian1DDist_d( ElDistMatrix_d L, ElInt nx );
EL_EXPORT ElError ElLaplacian1DDist_c( ElDistMatrix_c L, ElInt nx );
EL_EXPORT ElError ElLaplacian1DDist_z( ElDistMatrix_z L, ElInt nx );

EL_EXPORT ElError ElLaplacian1DSparse_s( ElSparseMatrix_s L, ElInt nx );
EL_EXPORT ElError ElLaplacian1DSparse_d( ElSparseMatrix_d L, ElInt nx );
EL_EXPORT ElError ElLaplacian1DSparse_c( ElSparseMatrix_c L, ElInt nx );
EL_EXPORT ElError ElLaplacian1DSparse_z( ElSparseMatrix_z L, ElInt nx );

EL_EXPORT ElError ElLaplacian1DDistSparse_s( ElDistSparseMatrix_s L, ElInt nx );
EL_EXPORT ElError ElLaplacian1DDistSparse_d( ElDistSparseMatrix_d L, ElInt nx );
EL_EXPORT ElError ElLaplacian1DDistSparse_c( ElDistSparseMatrix_c L, ElInt nx );
EL_EXPORT ElError ElLaplacian1DDistSparse_z( ElDistSparseMatrix_z L, ElInt nx );

EL_EXPORT ElError ElLaplacian2D_s( ElMatrix_s L, ElInt nx, ElInt ny );
EL_EXPORT ElError ElLaplacian2D_d( ElMatrix_d L, ElInt nx, ElInt ny );
EL_EXPORT ElError ElLaplacian2D_c( ElMatrix_c L, ElInt nx, ElInt ny );
EL_EXPORT ElError ElLaplacian2D_z( ElMatrix_z L, ElInt nx, ElInt ny );

EL_EXPORT ElError ElLaplacian2DDist_s( ElDistMatrix_s L, ElInt nx, ElInt ny );
EL_EXPORT ElError ElLaplacian2DDist_d( ElDistMatrix_d L, ElInt nx, ElInt ny );
EL_EXPORT ElError ElLaplacian2DDist_c( ElDistMatrix_c L, ElInt nx, ElInt ny );
EL_EXPORT ElError ElLaplacian2DDist_z( ElDistMatrix_z L, ElInt nx, ElInt ny );

EL_EXPORT ElError ElLaplacian2DSparse_s
( ElSparseMatrix_s L, ElInt nx, ElInt ny );
EL_EXPORT ElError ElLaplacian2DSparse_d
( ElSparseMatrix_d L, ElInt nx, ElInt ny );
EL_EXPORT ElError ElLaplacian2DSparse_c
( ElSparseMatrix_c L, ElInt nx, ElInt ny );
EL_EXPORT ElError ElLaplacian2DSparse_z
( ElSparseMatrix_z L, ElInt nx, ElInt ny );

EL_EXPORT ElError ElLaplacian2DDistSparse_s
( ElDistSparseMatrix_s L, ElInt nx, ElInt ny );
EL_EXPORT ElError ElLaplacian2DDistSparse_d
( ElDistSparseMatrix_d L, ElInt nx, ElInt ny );
EL_EXPORT ElError ElLaplacian2DDistSparse_c
( ElDistSparseMatrix_c L, ElInt nx, ElInt ny );
EL_EXPORT ElError ElLaplacian2DDistSparse_z
( ElDistSparseMatrix_z L, ElInt nx, ElInt ny );

EL_EXPORT ElError ElLaplacian3D_s( ElMatrix_s L, ElInt nx, ElInt ny, ElInt nz );
EL_EXPORT ElError ElLaplacian3D_d( ElMatrix_d L, ElInt nx, ElInt ny, ElInt nz );
EL_EXPORT ElError ElLaplacian3D_c( ElMatrix_c L, ElInt nx, ElInt ny, ElInt nz );
EL_EXPORT ElError ElLaplacian3D_z( ElMatrix_z L, ElInt nx, ElInt ny, ElInt nz );

EL_EXPORT ElError ElLaplacian3DDist_s
( ElDistMatrix_s L, ElInt nx, ElInt ny, ElInt nz );
EL_EXPORT ElError ElLaplacian3DDist_d
( ElDistMatrix_d L, ElInt nx, ElInt ny, ElInt nz );
EL_EXPORT ElError ElLaplacian3DDist_c
( ElDistMatrix_c L, ElInt nx, ElInt ny, ElInt nz );
EL_EXPORT ElError ElLaplacian3DDist_z
( ElDistMatrix_z L, ElInt nx, ElInt ny, ElInt nz );

EL_EXPORT ElError ElLaplacian3DDistSparse_s
( ElDistSparseMatrix_s L, ElInt nx, ElInt ny, ElInt nz );
EL_EXPORT ElError ElLaplacian3DDistSparse_d
( ElDistSparseMatrix_d L, ElInt nx, ElInt ny, ElInt nz );
EL_EXPORT ElError ElLaplacian3DDistSparse_c
( ElDistSparseMatrix_c L, ElInt nx, ElInt ny, ElInt nz );
EL_EXPORT ElError ElLaplacian3DDistSparse_z
( ElDistSparseMatrix_z L, ElInt nx, ElInt ny, ElInt nz );

EL_EXPORT ElError ElLaplacian3DDistSparse_s
( ElDistSparseMatrix_s L, ElInt nx, ElInt ny, ElInt nz );
EL_EXPORT ElError ElLaplacian3DDistSparse_d
( ElDistSparseMatrix_d L, ElInt nx, ElInt ny, ElInt nz );
EL_EXPORT ElError ElLaplacian3DDistSparse_c
( ElDistSparseMatrix_c L, ElInt nx, ElInt ny, ElInt nz );
EL_EXPORT ElError ElLaplacian3DDistSparse_z
( ElDistSparseMatrix_z L, ElInt nx, ElInt ny, ElInt nz );

/* Lauchli
   ======= */
EL_EXPORT ElError ElLauchli_i( ElMatrix_i A, ElInt n, ElInt mu );
EL_EXPORT ElError ElLauchli_s( ElMatrix_s A, ElInt n, float mu );
EL_EXPORT ElError ElLauchli_d( ElMatrix_d A, ElInt n, double mu );
EL_EXPORT ElError ElLauchli_c( ElMatrix_c A, ElInt n, complex_float mu );
EL_EXPORT ElError ElLauchli_z( ElMatrix_z A, ElInt n, complex_double mu );

EL_EXPORT ElError ElLauchliDist_i
( ElDistMatrix_i A, ElInt n, ElInt mu );
EL_EXPORT ElError ElLauchliDist_s
( ElDistMatrix_s A, ElInt n, float mu );
EL_EXPORT ElError ElLauchliDist_d
( ElDistMatrix_d A, ElInt n, double mu );
EL_EXPORT ElError ElLauchliDist_c
( ElDistMatrix_c A, ElInt n, complex_float mu );
EL_EXPORT ElError ElLauchliDist_z
( ElDistMatrix_z A, ElInt n, complex_double mu );

/* Legendre
   ======== */
EL_EXPORT ElError ElLegendre_s( ElMatrix_s A, ElInt n );
EL_EXPORT ElError ElLegendre_d( ElMatrix_d A, ElInt n );
EL_EXPORT ElError ElLegendre_c( ElMatrix_c A, ElInt n );
EL_EXPORT ElError ElLegendre_z( ElMatrix_z A, ElInt n );

EL_EXPORT ElError ElLegendreDist_s( ElDistMatrix_s A, ElInt n );
EL_EXPORT ElError ElLegendreDist_d( ElDistMatrix_d A, ElInt n );
EL_EXPORT ElError ElLegendreDist_c( ElDistMatrix_c A, ElInt n );
EL_EXPORT ElError ElLegendreDist_z( ElDistMatrix_z A, ElInt n );

/* Lehmer
   ====== */
EL_EXPORT ElError ElLehmer_s( ElMatrix_s L, ElInt n );
EL_EXPORT ElError ElLehmer_d( ElMatrix_d L, ElInt n );
EL_EXPORT ElError ElLehmer_c( ElMatrix_c L, ElInt n );
EL_EXPORT ElError ElLehmer_z( ElMatrix_z L, ElInt n );

EL_EXPORT ElError ElLehmerDist_s( ElDistMatrix_s L, ElInt n );
EL_EXPORT ElError ElLehmerDist_d( ElDistMatrix_d L, ElInt n );
EL_EXPORT ElError ElLehmerDist_c( ElDistMatrix_c L, ElInt n );
EL_EXPORT ElError ElLehmerDist_z( ElDistMatrix_z L, ElInt n );

/* Lotkin
   ====== */
EL_EXPORT ElError ElLotkin_s( ElMatrix_s A, ElInt n );
EL_EXPORT ElError ElLotkin_d( ElMatrix_d A, ElInt n );
EL_EXPORT ElError ElLotkin_c( ElMatrix_c A, ElInt n );
EL_EXPORT ElError ElLotkin_z( ElMatrix_z A, ElInt n );

EL_EXPORT ElError ElLotkinDist_s( ElDistMatrix_s A, ElInt n );
EL_EXPORT ElError ElLotkinDist_d( ElDistMatrix_d A, ElInt n );
EL_EXPORT ElError ElLotkinDist_c( ElDistMatrix_c A, ElInt n );
EL_EXPORT ElError ElLotkinDist_z( ElDistMatrix_z A, ElInt n );

/* MinIJ
   ===== */
EL_EXPORT ElError ElMinIJ_i( ElMatrix_i A, ElInt n );
EL_EXPORT ElError ElMinIJ_s( ElMatrix_s A, ElInt n );
EL_EXPORT ElError ElMinIJ_d( ElMatrix_d A, ElInt n );
EL_EXPORT ElError ElMinIJ_c( ElMatrix_c A, ElInt n );
EL_EXPORT ElError ElMinIJ_z( ElMatrix_z A, ElInt n );

EL_EXPORT ElError ElMinIJDist_i( ElDistMatrix_i A, ElInt n );
EL_EXPORT ElError ElMinIJDist_s( ElDistMatrix_s A, ElInt n );
EL_EXPORT ElError ElMinIJDist_d( ElDistMatrix_d A, ElInt n );
EL_EXPORT ElError ElMinIJDist_c( ElDistMatrix_c A, ElInt n );
EL_EXPORT ElError ElMinIJDist_z( ElDistMatrix_z A, ElInt n );

/* Ones
   ==== */
EL_EXPORT ElError ElOnes_i( ElMatrix_i A, ElInt m, ElInt n );
EL_EXPORT ElError ElOnes_s( ElMatrix_s A, ElInt m, ElInt n );
EL_EXPORT ElError ElOnes_d( ElMatrix_d A, ElInt m, ElInt n );
EL_EXPORT ElError ElOnes_c( ElMatrix_c A, ElInt m, ElInt n );
EL_EXPORT ElError ElOnes_z( ElMatrix_z A, ElInt m, ElInt n );

EL_EXPORT ElError ElOnesDist_i( ElDistMatrix_i A, ElInt m, ElInt n );
EL_EXPORT ElError ElOnesDist_s( ElDistMatrix_s A, ElInt m, ElInt n );
EL_EXPORT ElError ElOnesDist_d( ElDistMatrix_d A, ElInt m, ElInt n );
EL_EXPORT ElError ElOnesDist_c( ElDistMatrix_c A, ElInt m, ElInt n );
EL_EXPORT ElError ElOnesDist_z( ElDistMatrix_z A, ElInt m, ElInt n );

EL_EXPORT ElError ElOnesDistMultiVec_i( ElDistMultiVec_i A, ElInt m, ElInt n );
EL_EXPORT ElError ElOnesDistMultiVec_s( ElDistMultiVec_s A, ElInt m, ElInt n );
EL_EXPORT ElError ElOnesDistMultiVec_d( ElDistMultiVec_d A, ElInt m, ElInt n );
EL_EXPORT ElError ElOnesDistMultiVec_c( ElDistMultiVec_c A, ElInt m, ElInt n );
EL_EXPORT ElError ElOnesDistMultiVec_z( ElDistMultiVec_z A, ElInt m, ElInt n );

EL_EXPORT ElError ElOnesSparse_i( ElSparseMatrix_i A, ElInt m, ElInt n );
EL_EXPORT ElError ElOnesSparse_s( ElSparseMatrix_s A, ElInt m, ElInt n );
EL_EXPORT ElError ElOnesSparse_d( ElSparseMatrix_d A, ElInt m, ElInt n );
EL_EXPORT ElError ElOnesSparse_c( ElSparseMatrix_c A, ElInt m, ElInt n );
EL_EXPORT ElError ElOnesSparse_z( ElSparseMatrix_z A, ElInt m, ElInt n );

EL_EXPORT ElError ElOnesDistSparse_i
( ElDistSparseMatrix_i A, ElInt m, ElInt n );
EL_EXPORT ElError ElOnesDistSparse_s
( ElDistSparseMatrix_s A, ElInt m, ElInt n );
EL_EXPORT ElError ElOnesDistSparse_d
( ElDistSparseMatrix_d A, ElInt m, ElInt n );
EL_EXPORT ElError ElOnesDistSparse_c
( ElDistSparseMatrix_c A, ElInt m, ElInt n );
EL_EXPORT ElError ElOnesDistSparse_z
( ElDistSparseMatrix_z A, ElInt m, ElInt n );

/* 1-2-1
   ===== */
EL_EXPORT ElError ElOneTwoOne_i( ElMatrix_i A, ElInt n );
EL_EXPORT ElError ElOneTwoOne_s( ElMatrix_s A, ElInt n );
EL_EXPORT ElError ElOneTwoOne_d( ElMatrix_d A, ElInt n );
EL_EXPORT ElError ElOneTwoOne_c( ElMatrix_c A, ElInt n );
EL_EXPORT ElError ElOneTwoOne_z( ElMatrix_z A, ElInt n );

EL_EXPORT ElError ElOneTwoOneDist_i( ElDistMatrix_i A, ElInt n );
EL_EXPORT ElError ElOneTwoOneDist_s( ElDistMatrix_s A, ElInt n );
EL_EXPORT ElError ElOneTwoOneDist_d( ElDistMatrix_d A, ElInt n );
EL_EXPORT ElError ElOneTwoOneDist_c( ElDistMatrix_c A, ElInt n );
EL_EXPORT ElError ElOneTwoOneDist_z( ElDistMatrix_z A, ElInt n );

/* Parter
   ====== */
EL_EXPORT ElError ElParter_s( ElMatrix_s A, ElInt n );
EL_EXPORT ElError ElParter_d( ElMatrix_d A, ElInt n );
EL_EXPORT ElError ElParter_c( ElMatrix_c A, ElInt n );
EL_EXPORT ElError ElParter_z( ElMatrix_z A, ElInt n );

EL_EXPORT ElError ElParterDist_s( ElDistMatrix_s A, ElInt n );
EL_EXPORT ElError ElParterDist_d( ElDistMatrix_d A, ElInt n );
EL_EXPORT ElError ElParterDist_c( ElDistMatrix_c A, ElInt n );
EL_EXPORT ElError ElParterDist_z( ElDistMatrix_z A, ElInt n );

/* Pei
   === */
EL_EXPORT ElError ElPei_s( ElMatrix_s A, ElInt n, float alpha );
EL_EXPORT ElError ElPei_d( ElMatrix_d A, ElInt n, double alpha );
EL_EXPORT ElError ElPei_c( ElMatrix_c A, ElInt n, complex_float alpha );
EL_EXPORT ElError ElPei_z( ElMatrix_z A, ElInt n, complex_double alpha );

EL_EXPORT ElError ElPeiDist_s
( ElDistMatrix_s A, ElInt n, float alpha );
EL_EXPORT ElError ElPeiDist_d
( ElDistMatrix_d A, ElInt n, double alpha );
EL_EXPORT ElError ElPeiDist_c
( ElDistMatrix_c A, ElInt n, complex_float alpha );
EL_EXPORT ElError ElPeiDist_z
( ElDistMatrix_z A, ElInt n, complex_double alpha );

/* Redheffer
   ========= */
EL_EXPORT ElError ElRedheffer_i( ElMatrix_i A, ElInt n );
EL_EXPORT ElError ElRedheffer_s( ElMatrix_s A, ElInt n );
EL_EXPORT ElError ElRedheffer_d( ElMatrix_d A, ElInt n );
EL_EXPORT ElError ElRedheffer_c( ElMatrix_c A, ElInt n );
EL_EXPORT ElError ElRedheffer_z( ElMatrix_z A, ElInt n );

EL_EXPORT ElError ElRedhefferDist_i( ElDistMatrix_i A, ElInt n );
EL_EXPORT ElError ElRedhefferDist_s( ElDistMatrix_s A, ElInt n );
EL_EXPORT ElError ElRedhefferDist_d( ElDistMatrix_d A, ElInt n );
EL_EXPORT ElError ElRedhefferDist_c( ElDistMatrix_c A, ElInt n );
EL_EXPORT ElError ElRedhefferDist_z( ElDistMatrix_z A, ElInt n );

/* Riffle
   ====== */
EL_EXPORT ElError ElRiffle_s( ElMatrix_s P, ElInt n );
EL_EXPORT ElError ElRiffle_d( ElMatrix_d P, ElInt n );
EL_EXPORT ElError ElRiffle_c( ElMatrix_c P, ElInt n );
EL_EXPORT ElError ElRiffle_z( ElMatrix_z P, ElInt n );

EL_EXPORT ElError ElRiffleDist_s( ElDistMatrix_s P, ElInt n );
EL_EXPORT ElError ElRiffleDist_d( ElDistMatrix_d P, ElInt n );
EL_EXPORT ElError ElRiffleDist_c( ElDistMatrix_c P, ElInt n );
EL_EXPORT ElError ElRiffleDist_z( ElDistMatrix_z P, ElInt n );

EL_EXPORT ElError ElRiffleStationary_s( ElMatrix_s PInf, ElInt n );
EL_EXPORT ElError ElRiffleStationary_d( ElMatrix_d PInf, ElInt n );
EL_EXPORT ElError ElRiffleStationary_c( ElMatrix_c PInf, ElInt n );
EL_EXPORT ElError ElRiffleStationary_z( ElMatrix_z PInf, ElInt n );

EL_EXPORT ElError ElRiffleStationaryDist_s( ElDistMatrix_s PInf, ElInt n );
EL_EXPORT ElError ElRiffleStationaryDist_d( ElDistMatrix_d PInf, ElInt n );
EL_EXPORT ElError ElRiffleStationaryDist_c( ElDistMatrix_c PInf, ElInt n );
EL_EXPORT ElError ElRiffleStationaryDist_z( ElDistMatrix_z PInf, ElInt n );

EL_EXPORT ElError ElRiffleDecay_s( ElMatrix_s A, ElInt n );
EL_EXPORT ElError ElRiffleDecay_d( ElMatrix_d A, ElInt n );
EL_EXPORT ElError ElRiffleDecay_c( ElMatrix_c A, ElInt n );
EL_EXPORT ElError ElRiffleDecay_z( ElMatrix_z A, ElInt n );

EL_EXPORT ElError ElRiffleDecayDist_s( ElDistMatrix_s A, ElInt n );
EL_EXPORT ElError ElRiffleDecayDist_d( ElDistMatrix_d A, ElInt n );
EL_EXPORT ElError ElRiffleDecayDist_c( ElDistMatrix_c A, ElInt n );
EL_EXPORT ElError ElRiffleDecayDist_z( ElDistMatrix_z A, ElInt n );

/* Ris
   === */
EL_EXPORT ElError ElRis_s( ElMatrix_s A, ElInt n );
EL_EXPORT ElError ElRis_d( ElMatrix_d A, ElInt n );
EL_EXPORT ElError ElRis_c( ElMatrix_c A, ElInt n );
EL_EXPORT ElError ElRis_z( ElMatrix_z A, ElInt n );

EL_EXPORT ElError ElRisDist_s( ElDistMatrix_s A, ElInt n );
EL_EXPORT ElError ElRisDist_d( ElDistMatrix_d A, ElInt n );
EL_EXPORT ElError ElRisDist_c( ElDistMatrix_c A, ElInt n );
EL_EXPORT ElError ElRisDist_z( ElDistMatrix_z A, ElInt n );

/* Toeplitz
   ======== */
EL_EXPORT ElError ElToeplitz_i
( ElMatrix_i A, ElInt m, ElInt n, ElInt aSize, ElInt* aBuf );
EL_EXPORT ElError ElToeplitz_s
( ElMatrix_s A, ElInt m, ElInt n, ElInt aSize, float* aBuf );
EL_EXPORT ElError ElToeplitz_d
( ElMatrix_d A, ElInt m, ElInt n, ElInt aSize, double* aBuf );
EL_EXPORT ElError ElToeplitz_c
( ElMatrix_c A, ElInt m, ElInt n, ElInt aSize, complex_float* aBuf );
EL_EXPORT ElError ElToeplitz_z
( ElMatrix_z A, ElInt m, ElInt n, ElInt aSize, complex_double* aBuf );

EL_EXPORT ElError ElToeplitzDist_i
( ElDistMatrix_i A, ElInt m, ElInt n, ElInt aSize, ElInt* aBuf );
EL_EXPORT ElError ElToeplitzDist_s
( ElDistMatrix_s A, ElInt m, ElInt n, ElInt aSize, float* aBuf );
EL_EXPORT ElError ElToeplitzDist_d
( ElDistMatrix_d A, ElInt m, ElInt n, ElInt aSize, double* aBuf );
EL_EXPORT ElError ElToeplitzDist_c
( ElDistMatrix_c A, ElInt m, ElInt n, ElInt aSize, complex_float* aBuf );
EL_EXPORT ElError ElToeplitzDist_z
( ElDistMatrix_z A, ElInt m, ElInt n, ElInt aSize, complex_double* aBuf );

/* Trefethen-Embree
   ================ */
EL_EXPORT ElError ElTrefethenEmbree_c( ElMatrix_c A, ElInt n );
EL_EXPORT ElError ElTrefethenEmbree_z( ElMatrix_z A, ElInt n );

EL_EXPORT ElError ElTrefethenEmbreeDist_c( ElDistMatrix_c A, ElInt n );
EL_EXPORT ElError ElTrefethenEmbreeDist_z( ElDistMatrix_z A, ElInt n );

/* Triangle
   ======== */
EL_EXPORT ElError ElTriangle_s( ElMatrix_s A, ElInt n );
EL_EXPORT ElError ElTriangle_d( ElMatrix_d A, ElInt n );
EL_EXPORT ElError ElTriangle_c( ElMatrix_c A, ElInt n );
EL_EXPORT ElError ElTriangle_z( ElMatrix_z A, ElInt n );

EL_EXPORT ElError ElTriangleDist_s( ElDistMatrix_s A, ElInt n );
EL_EXPORT ElError ElTriangleDist_d( ElDistMatrix_d A, ElInt n );
EL_EXPORT ElError ElTriangleDist_c( ElDistMatrix_c A, ElInt n );
EL_EXPORT ElError ElTriangleDist_z( ElDistMatrix_z A, ElInt n );

/* TriW
   ==== */
EL_EXPORT ElError ElTriW_i
( ElMatrix_i A, ElInt n, ElInt alpha, ElInt k );
EL_EXPORT ElError ElTriW_s
( ElMatrix_s A, ElInt n, float alpha, ElInt k );
EL_EXPORT ElError ElTriW_d
( ElMatrix_d A, ElInt n, double alpha, ElInt k );
EL_EXPORT ElError ElTriW_c
( ElMatrix_c A, ElInt n, complex_float alpha, ElInt k );
EL_EXPORT ElError ElTriW_z
( ElMatrix_z A, ElInt n, complex_double alpha, ElInt k );

EL_EXPORT ElError ElTriWDist_i
( ElDistMatrix_i A, ElInt n, ElInt alpha, ElInt k );
EL_EXPORT ElError ElTriWDist_s
( ElDistMatrix_s A, ElInt n, float alpha, ElInt k );
EL_EXPORT ElError ElTriWDist_d
( ElDistMatrix_d A, ElInt n, double alpha, ElInt k );
EL_EXPORT ElError ElTriWDist_c
( ElDistMatrix_c A, ElInt n, complex_float alpha, ElInt k );
EL_EXPORT ElError ElTriWDist_z
( ElDistMatrix_z A, ElInt n, complex_double alpha, ElInt k );

/* Walsh
   ===== */
EL_EXPORT ElError ElWalsh_i( ElMatrix_i A, ElInt k, bool binary );
EL_EXPORT ElError ElWalsh_s( ElMatrix_s A, ElInt k, bool binary );
EL_EXPORT ElError ElWalsh_d( ElMatrix_d A, ElInt k, bool binary );
EL_EXPORT ElError ElWalsh_c( ElMatrix_c A, ElInt k, bool binary );
EL_EXPORT ElError ElWalsh_z( ElMatrix_z A, ElInt k, bool binary );

EL_EXPORT ElError ElWalshDist_i( ElDistMatrix_i A, ElInt k, bool binary );
EL_EXPORT ElError ElWalshDist_s( ElDistMatrix_s A, ElInt k, bool binary );
EL_EXPORT ElError ElWalshDist_d( ElDistMatrix_d A, ElInt k, bool binary );
EL_EXPORT ElError ElWalshDist_c( ElDistMatrix_c A, ElInt k, bool binary );
EL_EXPORT ElError ElWalshDist_z( ElDistMatrix_z A, ElInt k, bool binary );

/* Whale
   ===== */
EL_EXPORT ElError ElWhale_c( ElMatrix_c A, ElInt n );
EL_EXPORT ElError ElWhale_z( ElMatrix_z A, ElInt n );

EL_EXPORT ElError ElWhaleDist_c( ElDistMatrix_c A, ElInt n );
EL_EXPORT ElError ElWhaleDist_z( ElDistMatrix_z A, ElInt n );

/* Wilkinson
   ========= */
EL_EXPORT ElError ElWilkinson_i( ElMatrix_i A, ElInt k );
EL_EXPORT ElError ElWilkinson_s( ElMatrix_s A, ElInt k );
EL_EXPORT ElError ElWilkinson_d( ElMatrix_d A, ElInt k );
EL_EXPORT ElError ElWilkinson_c( ElMatrix_c A, ElInt k );
EL_EXPORT ElError ElWilkinson_z( ElMatrix_z A, ElInt k );

EL_EXPORT ElError ElWilkinsonDist_i( ElDistMatrix_i A, ElInt k );
EL_EXPORT ElError ElWilkinsonDist_s( ElDistMatrix_s A, ElInt k );
EL_EXPORT ElError ElWilkinsonDist_d( ElDistMatrix_d A, ElInt k );
EL_EXPORT ElError ElWilkinsonDist_c( ElDistMatrix_c A, ElInt k );
EL_EXPORT ElError ElWilkinsonDist_z( ElDistMatrix_z A, ElInt k );

/* Zeros
   ===== */
EL_EXPORT ElError ElZeros_i( ElMatrix_i A, ElInt m, ElInt n );
EL_EXPORT ElError ElZeros_s( ElMatrix_s A, ElInt m, ElInt n );
EL_EXPORT ElError ElZeros_d( ElMatrix_d A, ElInt m, ElInt n );
EL_EXPORT ElError ElZeros_c( ElMatrix_c A, ElInt m, ElInt n );
EL_EXPORT ElError ElZeros_z( ElMatrix_z A, ElInt m, ElInt n );

EL_EXPORT ElError ElZerosDist_i( ElDistMatrix_i A, ElInt m, ElInt n );
EL_EXPORT ElError ElZerosDist_s( ElDistMatrix_s A, ElInt m, ElInt n );
EL_EXPORT ElError ElZerosDist_d( ElDistMatrix_d A, ElInt m, ElInt n );
EL_EXPORT ElError ElZerosDist_c( ElDistMatrix_c A, ElInt m, ElInt n );
EL_EXPORT ElError ElZerosDist_z( ElDistMatrix_z A, ElInt m, ElInt n );

EL_EXPORT ElError ElZerosSparse_i( ElSparseMatrix_i A, ElInt m, ElInt n );
EL_EXPORT ElError ElZerosSparse_s( ElSparseMatrix_s A, ElInt m, ElInt n );
EL_EXPORT ElError ElZerosSparse_d( ElSparseMatrix_d A, ElInt m, ElInt n );
EL_EXPORT ElError ElZerosSparse_c( ElSparseMatrix_c A, ElInt m, ElInt n );
EL_EXPORT ElError ElZerosSparse_z( ElSparseMatrix_z A, ElInt m, ElInt n );

EL_EXPORT ElError ElZerosDistSparse_i
( ElDistSparseMatrix_i A, ElInt m, ElInt n );
EL_EXPORT ElError ElZerosDistSparse_s
( ElDistSparseMatrix_s A, ElInt m, ElInt n );
EL_EXPORT ElError ElZerosDistSparse_d
( ElDistSparseMatrix_d A, ElInt m, ElInt n );
EL_EXPORT ElError ElZerosDistSparse_c
( ElDistSparseMatrix_c A, ElInt m, ElInt n );
EL_EXPORT ElError ElZerosDistSparse_z
( ElDistSparseMatrix_z A, ElInt m, ElInt n );

EL_EXPORT ElError ElZerosDistMultiVec_i( ElDistMultiVec_i A, ElInt m, ElInt n );
EL_EXPORT ElError ElZerosDistMultiVec_s( ElDistMultiVec_s A, ElInt m, ElInt n );
EL_EXPORT ElError ElZerosDistMultiVec_d( ElDistMultiVec_d A, ElInt m, ElInt n );
EL_EXPORT ElError ElZerosDistMultiVec_c( ElDistMultiVec_c A, ElInt m, ElInt n );
EL_EXPORT ElError ElZerosDistMultiVec_z( ElDistMultiVec_z A, ElInt m, ElInt n );

/* Random
   ###### */

/* Bernoulli
   ========= */
EL_EXPORT ElError ElBernoulli_i( ElMatrix_i A, ElInt m, ElInt n, double p );
EL_EXPORT ElError ElBernoulli_s( ElMatrix_s A, ElInt m, ElInt n, double p );
EL_EXPORT ElError ElBernoulli_d( ElMatrix_d A, ElInt m, ElInt n, double p );
EL_EXPORT ElError ElBernoulli_c( ElMatrix_c A, ElInt m, ElInt n, double p );
EL_EXPORT ElError ElBernoulli_z( ElMatrix_z A, ElInt m, ElInt n, double p );

EL_EXPORT ElError ElBernoulliDist_i
( ElDistMatrix_i A, ElInt m, ElInt n, double p );
EL_EXPORT ElError ElBernoulliDist_s
( ElDistMatrix_s A, ElInt m, ElInt n, double p );
EL_EXPORT ElError ElBernoulliDist_d
( ElDistMatrix_d A, ElInt m, ElInt n, double p );
EL_EXPORT ElError ElBernoulliDist_c
( ElDistMatrix_c A, ElInt m, ElInt n, double p );
EL_EXPORT ElError ElBernoulliDist_z
( ElDistMatrix_z A, ElInt m, ElInt n, double p );

/* Gaussian
   ======== */
EL_EXPORT ElError ElGaussian_s
( ElMatrix_s A, ElInt m, ElInt n, float mean, float stddev );
EL_EXPORT ElError ElGaussian_d
( ElMatrix_d A, ElInt m, ElInt n, double mean, double stddev );
EL_EXPORT ElError ElGaussian_c
( ElMatrix_c A, ElInt m, ElInt n, complex_float mean, float stddev );
EL_EXPORT ElError ElGaussian_z
( ElMatrix_z A, ElInt m, ElInt n, complex_double mean, double stddev );

EL_EXPORT ElError ElGaussianDist_s
( ElDistMatrix_s A, ElInt m, ElInt n, float mean, float stddev );
EL_EXPORT ElError ElGaussianDist_d
( ElDistMatrix_d A, ElInt m, ElInt n, double mean, double stddev );
EL_EXPORT ElError ElGaussianDist_c
( ElDistMatrix_c A, ElInt m, ElInt n, complex_float mean, float stddev );
EL_EXPORT ElError ElGaussianDist_z
( ElDistMatrix_z A, ElInt m, ElInt n, complex_double mean, double stddev );

EL_EXPORT ElError ElGaussianDistMultiVec_s
( ElDistMultiVec_s A, ElInt m, ElInt n, float mean, float stddev );
EL_EXPORT ElError ElGaussianDistMultiVec_d
( ElDistMultiVec_d A, ElInt m, ElInt n, double mean, double stddev );
EL_EXPORT ElError ElGaussianDistMultiVec_c
( ElDistMultiVec_c A, ElInt m, ElInt n, complex_float mean, float stddev );
EL_EXPORT ElError ElGaussianDistMultiVec_z
( ElDistMultiVec_z A, ElInt m, ElInt n, complex_double mean, double stddev );

/* Haar 
   ==== */
EL_EXPORT ElError ElHaar_s( ElMatrix_s A, ElInt n );
EL_EXPORT ElError ElHaar_d( ElMatrix_d A, ElInt n );
EL_EXPORT ElError ElHaar_c( ElMatrix_c A, ElInt n );
EL_EXPORT ElError ElHaar_z( ElMatrix_z A, ElInt n );

EL_EXPORT ElError ElHaarDist_s( ElDistMatrix_s A, ElInt n );
EL_EXPORT ElError ElHaarDist_d( ElDistMatrix_d A, ElInt n );
EL_EXPORT ElError ElHaarDist_c( ElDistMatrix_c A, ElInt n );
EL_EXPORT ElError ElHaarDist_z( ElDistMatrix_z A, ElInt n );

EL_EXPORT ElError ElImplicitHaar_s
( ElMatrix_s A, ElMatrix_s t, ElMatrix_s d, ElInt n );
EL_EXPORT ElError ElImplicitHaar_d
( ElMatrix_d A, ElMatrix_d t, ElMatrix_d d, ElInt n );
EL_EXPORT ElError ElImplicitHaar_c
( ElMatrix_c A, ElMatrix_c t, ElMatrix_s d, ElInt n );
EL_EXPORT ElError ElImplicitHaar_z
( ElMatrix_z A, ElMatrix_z t, ElMatrix_d d, ElInt n );

EL_EXPORT ElError ElImplicitHaarDist_s
( ElDistMatrix_s A, ElDistMatrix_s t, ElDistMatrix_s d, ElInt n );
EL_EXPORT ElError ElImplicitHaarDist_d
( ElDistMatrix_d A, ElDistMatrix_d t, ElDistMatrix_d d, ElInt n );
EL_EXPORT ElError ElImplicitHaarDist_c
( ElDistMatrix_c A, ElDistMatrix_c t, ElDistMatrix_s d, ElInt n );
EL_EXPORT ElError ElImplicitHaarDist_z
( ElDistMatrix_z A, ElDistMatrix_z t, ElDistMatrix_d d, ElInt n );

/* Hatano-Nelson
   ============= */
EL_EXPORT ElError ElHatanoNelson_s
( ElMatrix_s A, ElInt n, float center, float radius, float g, 
  bool periodic );
EL_EXPORT ElError ElHatanoNelson_d
( ElMatrix_d A, ElInt n, double center, double radius, double g, 
  bool periodic );
EL_EXPORT ElError ElHatanoNelson_c
( ElMatrix_c A, ElInt n, complex_float center, float radius, complex_float g, 
  bool periodic );
EL_EXPORT ElError ElHatanoNelson_z
( ElMatrix_z A, ElInt n, complex_double center, double radius, complex_double g,
  bool periodic );

EL_EXPORT ElError ElHatanoNelsonDist_s
( ElDistMatrix_s A, ElInt n, float center, float radius, 
  float g, bool periodic );
EL_EXPORT ElError ElHatanoNelsonDist_d
( ElDistMatrix_d A, ElInt n, double center, double radius, 
  double g, bool periodic );
EL_EXPORT ElError ElHatanoNelsonDist_c
( ElDistMatrix_c A, ElInt n, complex_float center, float radius, 
  complex_float g, bool periodic );
EL_EXPORT ElError ElHatanoNelsonDist_z
( ElDistMatrix_z A, ElInt n, complex_double center, double radius, 
  complex_double g, bool periodic );

/* Hermitian uniform spectrum
   ========================== */
EL_EXPORT ElError ElHermitianUniformSpectrum_s
( ElMatrix_s A, ElInt n, float lower, float upper );
EL_EXPORT ElError ElHermitianUniformSpectrum_d
( ElMatrix_d A, ElInt n, double lower, double upper );
EL_EXPORT ElError ElHermitianUniformSpectrum_c
( ElMatrix_c A, ElInt n, float lower, float upper );
EL_EXPORT ElError ElHermitianUniformSpectrum_z
( ElMatrix_z A, ElInt n, double lower, double upper );

EL_EXPORT ElError ElHermitianUniformSpectrumDist_s
( ElDistMatrix_s A, ElInt n, float lower, float upper );
EL_EXPORT ElError ElHermitianUniformSpectrumDist_d
( ElDistMatrix_d A, ElInt n, double lower, double upper );
EL_EXPORT ElError ElHermitianUniformSpectrumDist_c
( ElDistMatrix_c A, ElInt n, float lower, float upper );
EL_EXPORT ElError ElHermitianUniformSpectrumDist_z
( ElDistMatrix_z A, ElInt n, double lower, double upper );

/* Normal uniform spectrum
   ======================= */
EL_EXPORT ElError ElNormalUniformSpectrum_c
( ElMatrix_c A, ElInt n, complex_float center, float radius );
EL_EXPORT ElError ElNormalUniformSpectrum_z
( ElMatrix_z A, ElInt n, complex_double center, double radius );

EL_EXPORT ElError ElNormalUniformSpectrumDist_c
( ElDistMatrix_c A, ElInt n, complex_float center, float radius );
EL_EXPORT ElError ElNormalUniformSpectrumDist_z
( ElDistMatrix_z A, ElInt n, complex_double center, double radius );

/* Rademacher 
   ========== */
EL_EXPORT ElError ElRademacher_i( ElMatrix_i A, ElInt m, ElInt n );
EL_EXPORT ElError ElRademacher_s( ElMatrix_s A, ElInt m, ElInt n );
EL_EXPORT ElError ElRademacher_d( ElMatrix_d A, ElInt m, ElInt n );
EL_EXPORT ElError ElRademacher_c( ElMatrix_c A, ElInt m, ElInt n );
EL_EXPORT ElError ElRademacher_z( ElMatrix_z A, ElInt m, ElInt n );

EL_EXPORT ElError ElRademacherDist_i( ElDistMatrix_i A, ElInt m, ElInt n );
EL_EXPORT ElError ElRademacherDist_s( ElDistMatrix_s A, ElInt m, ElInt n );
EL_EXPORT ElError ElRademacherDist_d( ElDistMatrix_d A, ElInt m, ElInt n );
EL_EXPORT ElError ElRademacherDist_c( ElDistMatrix_c A, ElInt m, ElInt n );
EL_EXPORT ElError ElRademacherDist_z( ElDistMatrix_z A, ElInt m, ElInt n );

/* Three-valued
   ============ */
EL_EXPORT ElError ElThreeValued_i( ElMatrix_i A, ElInt m, ElInt n, double p );
EL_EXPORT ElError ElThreeValued_s( ElMatrix_s A, ElInt m, ElInt n, double p );
EL_EXPORT ElError ElThreeValued_d( ElMatrix_d A, ElInt m, ElInt n, double p );
EL_EXPORT ElError ElThreeValued_c( ElMatrix_c A, ElInt m, ElInt n, double p );
EL_EXPORT ElError ElThreeValued_z( ElMatrix_z A, ElInt m, ElInt n, double p );

EL_EXPORT ElError ElThreeValuedDist_i
( ElDistMatrix_i A, ElInt m, ElInt n, double p );
EL_EXPORT ElError ElThreeValuedDist_s
( ElDistMatrix_s A, ElInt m, ElInt n, double p );
EL_EXPORT ElError ElThreeValuedDist_d
( ElDistMatrix_d A, ElInt m, ElInt n, double p );
EL_EXPORT ElError ElThreeValuedDist_c
( ElDistMatrix_c A, ElInt m, ElInt n, double p );
EL_EXPORT ElError ElThreeValuedDist_z
( ElDistMatrix_z A, ElInt m, ElInt n, double p );

/* Uniform
   ======= */
EL_EXPORT ElError ElUniform_i
( ElMatrix_i A, ElInt m, ElInt n, ElInt center, ElInt radius );
EL_EXPORT ElError ElUniform_s
( ElMatrix_s A, ElInt m, ElInt n, float center, float radius );
EL_EXPORT ElError ElUniform_d
( ElMatrix_d A, ElInt m, ElInt n, double center, double radius );
EL_EXPORT ElError ElUniform_c
( ElMatrix_c A, ElInt m, ElInt n, complex_float center, float radius );
EL_EXPORT ElError ElUniform_z
( ElMatrix_z A, ElInt m, ElInt n, complex_double center, double radius );

EL_EXPORT ElError ElUniformDist_i
( ElDistMatrix_i A, ElInt m, ElInt n, ElInt center, ElInt radius );
EL_EXPORT ElError ElUniformDist_s
( ElDistMatrix_s A, ElInt m, ElInt n, float center, float radius );
EL_EXPORT ElError ElUniformDist_d
( ElDistMatrix_d A, ElInt m, ElInt n, double center, double radius );
EL_EXPORT ElError ElUniformDist_c
( ElDistMatrix_c A, ElInt m, ElInt n, complex_float center, float radius );
EL_EXPORT ElError ElUniformDist_z
( ElDistMatrix_z A, ElInt m, ElInt n, complex_double center, double radius );

EL_EXPORT ElError ElUniformDistMultiVec_i
( ElDistMultiVec_i A, ElInt m, ElInt n, ElInt center, ElInt radius );
EL_EXPORT ElError ElUniformDistMultiVec_s
( ElDistMultiVec_s A, ElInt m, ElInt n, float center, float radius );
EL_EXPORT ElError ElUniformDistMultiVec_d
( ElDistMultiVec_d A, ElInt m, ElInt n, double center, double radius );
EL_EXPORT ElError ElUniformDistMultiVec_c
( ElDistMultiVec_c A, ElInt m, ElInt n, complex_float center, float radius );
EL_EXPORT ElError ElUniformDistMultiVec_z
( ElDistMultiVec_z A, ElInt m, ElInt n, complex_double center, double radius );

/* Uniform Helmholtz Green's
   ========================= */
EL_EXPORT ElError ElUniformHelmholtzGreens_c
( ElMatrix_c A, ElInt n, float lambda );
EL_EXPORT ElError ElUniformHelmholtzGreens_z
( ElMatrix_z A, ElInt n, double lambda );

EL_EXPORT ElError ElUniformHelmholtzGreensDist_c
( ElDistMatrix_c A, ElInt n, float lambda );
EL_EXPORT ElError ElUniformHelmholtzGreensDist_z
( ElDistMatrix_z A, ElInt n, double lambda );

/* Wigner
   ====== */
EL_EXPORT ElError ElWigner_s
( ElMatrix_s A, ElInt n, float mean, float stddev );
EL_EXPORT ElError ElWigner_d
( ElMatrix_d A, ElInt n, double mean, double stddev );
EL_EXPORT ElError ElWigner_c
( ElMatrix_c A, ElInt n, complex_float mean, float stddev );
EL_EXPORT ElError ElWigner_z
( ElMatrix_z A, ElInt n, complex_double mean, double stddev );

EL_EXPORT ElError ElWignerDist_s
( ElDistMatrix_s A, ElInt n, float mean, float stddev );
EL_EXPORT ElError ElWignerDist_d
( ElDistMatrix_d A, ElInt n, double mean, double stddev );
EL_EXPORT ElError ElWignerDist_c
( ElDistMatrix_c A, ElInt n, complex_float mean, float stddev );
EL_EXPORT ElError ElWignerDist_z
( ElDistMatrix_z A, ElInt n, complex_double mean, double stddev );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_MATRICES_C_H */
