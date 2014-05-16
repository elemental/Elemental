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

ElError ElUniform_s( ElMatrix_s A, Int m, Int n );
ElError ElUniform_d( ElMatrix_d A, Int m, Int n );
ElError ElUniform_c( ElMatrix_c A, Int m, Int n );
ElError ElUniform_z( ElMatrix_z A, Int m, Int n );

ElError ElUniform_s( ElDistMatrix_s A, Int m, Int n );
ElError ElUniform_d( ElDistMatrix_d A, Int m, Int n );
ElError ElUniform_c( ElDistMatrix_c A, Int m, Int n );
ElError ElUniform_z( ElDistMatrix_z A, Int m, Int n );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_UNIFORM_C_H */
