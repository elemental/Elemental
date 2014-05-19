/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_SPY_C_H
#define EL_SPY_C_H

#ifdef __cplusplus
extern "C" {
#endif

ElError ElSpyMatrix_s( ElConstMatrix_s A, const char* title, float tol );
ElError ElSpyMatrix_d( ElConstMatrix_d A, const char* title, double tol );
ElError ElSpyMatrix_c( ElConstMatrix_c A, const char* title, float tol );
ElError ElSpyMatrix_z( ElConstMatrix_z A, const char* title, double tol );

ElError ElSpyDistMatrix_s
( ElConstDistMatrix_s A, const char* title, float tol );
ElError ElSpyDistMatrix_d
( ElConstDistMatrix_d A, const char* title, double tol );
ElError ElSpyDistMatrix_c
( ElConstDistMatrix_c A, const char* title, float tol );
ElError ElSpyDistMatrix_z
( ElConstDistMatrix_z A, const char* title, double tol );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_SPY_C_H */
