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

ElError ElSpy_s( ElConstMatrix_s A, const char* title, float tol );
ElError ElSpy_d( ElConstMatrix_d A, const char* title, double tol );
ElError ElSpy_c( ElConstMatrix_c A, const char* title, float tol );
ElError ElSpy_z( ElConstMatrix_z A, const char* title, double tol );

ElError ElSpyDist_s( ElConstDistMatrix_s A, const char* title, float tol );
ElError ElSpyDist_d( ElConstDistMatrix_d A, const char* title, double tol );
ElError ElSpyDist_c( ElConstDistMatrix_c A, const char* title, float tol );
ElError ElSpyDist_z( ElConstDistMatrix_z A, const char* title, double tol );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_SPY_C_H */
