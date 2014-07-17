/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_UNIFORM_C_H
#define EL_UNIFORM_C_H

#ifdef __cplusplus
extern "C" {
#endif

ElError ElUniform_s
( ElMatrix_s A, ElInt m, ElInt n, float center, float radius );
ElError ElUniform_d
( ElMatrix_d A, ElInt m, ElInt n, double center, double radius );
ElError ElUniform_c
( ElMatrix_c A, ElInt m, ElInt n, complex_float center, float radius );
ElError ElUniform_z
( ElMatrix_z A, ElInt m, ElInt n, complex_double center, double radius );

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

#endif /* ifndef EL_UNIFORM_C_H */
