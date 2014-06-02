/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_BLAS_LEVEL1_C_H
#define EL_BLAS_LEVEL1_C_H

#include "El/core/DistMatrix-C.h"

#ifdef __cplusplus
extern "C" {
#endif

/* B = A 
   ----- */
ElError ElCopyDistMatrix_s( ElConstDistMatrix_s A, ElDistMatrix_s B );
ElError ElCopyDistMatrix_d( ElConstDistMatrix_d A, ElDistMatrix_d B );
ElError ElCopyDistMatrix_c( ElConstDistMatrix_c A, ElDistMatrix_c B );
ElError ElCopyDistMatrix_z( ElConstDistMatrix_z A, ElDistMatrix_z B );

/* B = A^T 
   ------- */
ElError ElTransposeDistMatrix_s( ElConstDistMatrix_s A, ElDistMatrix_s B );
ElError ElTransposeDistMatrix_d( ElConstDistMatrix_d A, ElDistMatrix_d B );
ElError ElTransposeDistMatrix_c( ElConstDistMatrix_c A, ElDistMatrix_c B );
ElError ElTransposeDistMatrix_z( ElConstDistMatrix_z A, ElDistMatrix_z B );

/* B = A^H 
   ------- */
ElError ElAdjointDistMatrix_s( ElConstDistMatrix_s A, ElDistMatrix_s B );
ElError ElAdjointDistMatrix_d( ElConstDistMatrix_d A, ElDistMatrix_d B );
ElError ElAdjointDistMatrix_c( ElConstDistMatrix_c A, ElDistMatrix_c B );
ElError ElAdjointDistMatrix_z( ElConstDistMatrix_z A, ElDistMatrix_z B );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_BLAS_LEVEL1_C_H */
