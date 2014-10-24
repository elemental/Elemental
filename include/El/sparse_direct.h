/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_SPARSEDIRECT_C_H
#define EL_SPARSEDIRECT_C_H

#ifdef __cplusplus
extern "C" {
#endif

/* Symmetric/Hermitian solves
   ========================== */
EL_EXPORT ElError ElSymmetricSolveSparseDist_s
( ElConstDistSparseMatrix_s A, ElDistMultiVec_s X );
EL_EXPORT ElError ElSymmetricSolveSparseDist_d
( ElConstDistSparseMatrix_d A, ElDistMultiVec_d X );
EL_EXPORT ElError ElSymmetricSolveSparseDist_c
( ElConstDistSparseMatrix_c A, ElDistMultiVec_c X );
EL_EXPORT ElError ElSymmetricSolveSparseDist_z
( ElConstDistSparseMatrix_z A, ElDistMultiVec_z X );

EL_EXPORT ElError ElHermitianSolveSparseDist_c
( ElConstDistSparseMatrix_c A, ElDistMultiVec_c X );
EL_EXPORT ElError ElHermitianSolveSparseDist_z
( ElConstDistSparseMatrix_z A, ElDistMultiVec_z X );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_SPARSEDIRECT_C_H */
