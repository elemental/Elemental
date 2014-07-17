/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_PRINT_C_H
#define EL_PRINT_C_H

#ifdef __cplusplus
extern "C" {
#endif

/* 
  TODO: Extend to support more general output streams. This is a bit 
        problematic due to difficulties in converting between FILE* and
        std::ostream 
*/
ElError ElPrint_s( ElConstMatrix_s A, const char* title );
ElError ElPrint_d( ElConstMatrix_d A, const char* title );
ElError ElPrint_c( ElConstMatrix_c A, const char* title );
ElError ElPrint_z( ElConstMatrix_z A, const char* title );

ElError ElPrintDist_s( ElConstDistMatrix_s A, const char* title );
ElError ElPrintDist_d( ElConstDistMatrix_d A, const char* title );
ElError ElPrintDist_c( ElConstDistMatrix_c A, const char* title );
ElError ElPrintDist_z( ElConstDistMatrix_z A, const char* title );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_PRINT_C_H */
