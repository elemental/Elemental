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
