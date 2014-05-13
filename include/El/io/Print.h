/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_PRINT_CINT_H
#define EL_PRINT_CINT_H

#ifdef __cplusplus
extern "C" {
#endif

/* 
  TODO: Extend to support more general output streams. This is a bit 
        problematic due to difficulties in converting between FILE* and
        std::ostream 
*/
void ElPrintMatrix_s( const ElMatrix_s* A, const char* title );
void ElPrintMatrix_d( const ElMatrix_d* A, const char* title );
void ElPrintMatrix_c( const ElMatrix_c* A, const char* title );
void ElPrintMatrix_z( const ElMatrix_z* A, const char* title );

/* TODO: Extend to DynamicDistMatrix */

#ifdef __cplusplus
} // extern "C"
#endif

#endif // ifndef EL_PRINT_CINT_H
