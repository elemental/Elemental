/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_DISPLAY_C_H
#define EL_DISPLAY_C_H

#ifdef __cplusplus
extern "C" {
#endif

void ElDisplayMatrix_s( const ElMatrix_s* A, const char* title );
void ElDisplayMatrix_d( const ElMatrix_d* A, const char* title );
void ElDisplayMatrix_c( const ElMatrix_c* A, const char* title );
void ElDisplayMatrix_z( const ElMatrix_z* A, const char* title );

/* TODO: Extend to DynamicDistMatrix */

#ifdef __cplusplus
} // extern "C"
#endif

#endif // ifndef EL_DISPLAY_C_H
