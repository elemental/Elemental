/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_ONES_C_H
#define EL_ONES_C_H

#ifdef __cplusplus
extern "C" {
#endif

ElError ElOnesMatrix_s( ElMatrix_s A, ElInt m, ElInt n );
ElError ElOnesMatrix_d( ElMatrix_d A, ElInt m, ElInt n );
ElError ElOnesMatrix_c( ElMatrix_c A, ElInt m, ElInt n );
ElError ElOnesMatrix_z( ElMatrix_z A, ElInt m, ElInt n );

ElError ElOnesDistMatrix_s( ElDistMatrix_s A, ElInt m, ElInt n );
ElError ElOnesDistMatrix_d( ElDistMatrix_d A, ElInt m, ElInt n );
ElError ElOnesDistMatrix_c( ElDistMatrix_c A, ElInt m, ElInt n );
ElError ElOnesDistMatrix_z( ElDistMatrix_z A, ElInt m, ElInt n );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_ONES_C_H */
