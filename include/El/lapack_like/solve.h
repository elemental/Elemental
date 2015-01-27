/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_LAPACK_SOLVE_C_H
#define EL_LAPACK_SOLVE_C_H

#include "El/core/DistMatrix.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Linear solvers
   ============== */

/* Linear solve
   ------------ */
EL_EXPORT ElError ElLinearSolve_s( ElMatrix_s A, ElMatrix_s B );
EL_EXPORT ElError ElLinearSolve_d( ElMatrix_d A, ElMatrix_d B );
EL_EXPORT ElError ElLinearSolve_c( ElMatrix_c A, ElMatrix_c B );
EL_EXPORT ElError ElLinearSolve_z( ElMatrix_z A, ElMatrix_z B );

EL_EXPORT ElError ElLinearSolveDist_s( ElDistMatrix_s A, ElDistMatrix_s B );
EL_EXPORT ElError ElLinearSolveDist_d( ElDistMatrix_d A, ElDistMatrix_d B );
EL_EXPORT ElError ElLinearSolveDist_c( ElDistMatrix_c A, ElDistMatrix_c B );
EL_EXPORT ElError ElLinearSolveDist_z( ElDistMatrix_z A, ElDistMatrix_z B );

EL_EXPORT ElError ElLinearSolveSparse_s( ElSparseMatrix_s A, ElMatrix_s B );
EL_EXPORT ElError ElLinearSolveSparse_d( ElSparseMatrix_d A, ElMatrix_d B );
EL_EXPORT ElError ElLinearSolveSparse_c( ElSparseMatrix_c A, ElMatrix_c B );
EL_EXPORT ElError ElLinearSolveSparse_z( ElSparseMatrix_z A, ElMatrix_z B );

EL_EXPORT ElError ElLinearSolveDistSparse_s
( ElDistSparseMatrix_s A, ElDistMultiVec_s B );
EL_EXPORT ElError ElLinearSolveDistSparse_d
( ElDistSparseMatrix_d A, ElDistMultiVec_d B );
EL_EXPORT ElError ElLinearSolveDistSparse_c
( ElDistSparseMatrix_c A, ElDistMultiVec_c B );
EL_EXPORT ElError ElLinearSolveDistSparse_z
( ElDistSparseMatrix_z A, ElDistMultiVec_z B );

/* Symmetric solve
   --------------- */
EL_EXPORT ElError ElSymmetricSolve_s
( ElUpperOrLower uplo, ElOrientation orientation, ElMatrix_s A, ElMatrix_s B );
EL_EXPORT ElError ElSymmetricSolve_d
( ElUpperOrLower uplo, ElOrientation orientation, ElMatrix_d A, ElMatrix_d B );
EL_EXPORT ElError ElSymmetricSolve_c
( ElUpperOrLower uplo, ElOrientation orientation, ElMatrix_c A, ElMatrix_c B );
EL_EXPORT ElError ElSymmetricSolve_z
( ElUpperOrLower uplo, ElOrientation orientation, ElMatrix_z A, ElMatrix_z B );

EL_EXPORT ElError ElSymmetricSolveDist_s
( ElUpperOrLower uplo, ElOrientation orientation, 
  ElDistMatrix_s A, ElDistMatrix_s B );
EL_EXPORT ElError ElSymmetricSolveDist_d
( ElUpperOrLower uplo, ElOrientation orientation, 
  ElDistMatrix_d A, ElDistMatrix_d B );
EL_EXPORT ElError ElSymmetricSolveDist_c
( ElUpperOrLower uplo, ElOrientation orientation, 
  ElDistMatrix_c A, ElDistMatrix_c B );
EL_EXPORT ElError ElSymmetricSolveDist_z
( ElUpperOrLower uplo, ElOrientation orientation, 
  ElDistMatrix_z A, ElDistMatrix_z B );

/* TODO: Expert version for choosing pivot strategy */

EL_EXPORT ElError ElSymmetricSolveDistSparse_s
( ElConstDistSparseMatrix_s A, ElDistMultiVec_s X );
EL_EXPORT ElError ElSymmetricSolveDistSparse_d
( ElConstDistSparseMatrix_d A, ElDistMultiVec_d X );
EL_EXPORT ElError ElSymmetricSolveDistSparse_c
( ElConstDistSparseMatrix_c A, ElDistMultiVec_c X );
EL_EXPORT ElError ElSymmetricSolveDistSparse_z
( ElConstDistSparseMatrix_z A, ElDistMultiVec_z X );

/* Hermitian solve
   --------------- */
EL_EXPORT ElError ElHermitianSolve_c
( ElUpperOrLower uplo, ElOrientation orientation, ElMatrix_c A, ElMatrix_c B );
EL_EXPORT ElError ElHermitianSolve_z
( ElUpperOrLower uplo, ElOrientation orientation, ElMatrix_z A, ElMatrix_z B );

EL_EXPORT ElError ElHermitianSolveDist_c
( ElUpperOrLower uplo, ElOrientation orientation, 
  ElDistMatrix_c A, ElDistMatrix_c B );
EL_EXPORT ElError ElHermitianSolveDist_z
( ElUpperOrLower uplo, ElOrientation orientation, 
  ElDistMatrix_z A, ElDistMatrix_z B );

/* TODO: Expert version for choosing pivot strategy */

EL_EXPORT ElError ElHermitianSolveDistSparse_c
( ElConstDistSparseMatrix_c A, ElDistMultiVec_c X );
EL_EXPORT ElError ElHermitianSolveDistSparse_z
( ElConstDistSparseMatrix_z A, ElDistMultiVec_z X );

/* HPD solve
   --------- */
EL_EXPORT ElError ElHPDSolve_s
( ElUpperOrLower uplo, ElOrientation orientation, ElMatrix_s A, ElMatrix_s B );
EL_EXPORT ElError ElHPDSolve_d
( ElUpperOrLower uplo, ElOrientation orientation, ElMatrix_d A, ElMatrix_d B );
EL_EXPORT ElError ElHPDSolve_c
( ElUpperOrLower uplo, ElOrientation orientation, ElMatrix_c A, ElMatrix_c B );
EL_EXPORT ElError ElHPDSolve_z
( ElUpperOrLower uplo, ElOrientation orientation, ElMatrix_z A, ElMatrix_z B );

EL_EXPORT ElError ElHPDSolveDist_s
( ElUpperOrLower uplo, ElOrientation orientation, 
  ElDistMatrix_s A, ElDistMatrix_s B );
EL_EXPORT ElError ElHPDSolveDist_d
( ElUpperOrLower uplo, ElOrientation orientation, 
  ElDistMatrix_d A, ElDistMatrix_d B );
EL_EXPORT ElError ElHPDSolveDist_c
( ElUpperOrLower uplo, ElOrientation orientation, 
  ElDistMatrix_c A, ElDistMatrix_c B );
EL_EXPORT ElError ElHPDSolveDist_z
( ElUpperOrLower uplo, ElOrientation orientation, 
  ElDistMatrix_z A, ElDistMatrix_z B );

/* TODO: Hessenberg solve */

/* Multi-shift Hessenberg solve
   ---------------------------- */
EL_EXPORT ElError ElMultiShiftHessSolve_s
( ElUpperOrLower uplo, ElOrientation orientation, float alpha,
  ElConstMatrix_s H, ElConstMatrix_s shifts, ElMatrix_s X );
EL_EXPORT ElError ElMultiShiftHessSolve_d
( ElUpperOrLower uplo, ElOrientation orientation, double alpha,
  ElConstMatrix_d H, ElConstMatrix_d shifts, ElMatrix_d X );
EL_EXPORT ElError ElMultiShiftHessSolve_c
( ElUpperOrLower uplo, ElOrientation orientation, complex_float alpha,
  ElConstMatrix_c H, ElConstMatrix_c shifts, ElMatrix_c X );
EL_EXPORT ElError ElMultiShiftHessSolve_z
( ElUpperOrLower uplo, ElOrientation orientation, complex_double alpha,
  ElConstMatrix_z H, ElConstMatrix_z shifts, ElMatrix_z X );

EL_EXPORT ElError ElMultiShiftHessSolveDist_s
( ElUpperOrLower uplo, ElOrientation orientation, float alpha,
  ElConstDistMatrix_s H, ElConstDistMatrix_s shifts, ElDistMatrix_s X );
EL_EXPORT ElError ElMultiShiftHessSolveDist_d
( ElUpperOrLower uplo, ElOrientation orientation, double alpha,
  ElConstDistMatrix_d H, ElConstDistMatrix_d shifts, ElDistMatrix_d X );
EL_EXPORT ElError ElMultiShiftHessSolveDist_c
( ElUpperOrLower uplo, ElOrientation orientation, complex_float alpha,
  ElConstDistMatrix_c H, ElConstDistMatrix_c shifts, ElDistMatrix_c X );
EL_EXPORT ElError ElMultiShiftHessSolveDist_z
( ElUpperOrLower uplo, ElOrientation orientation, complex_double alpha,
  ElConstDistMatrix_z H, ElConstDistMatrix_z shifts, ElDistMatrix_z X );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_LAPACK_SOLVE_C_H */
