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

ElError ElDisplayMatrix_s( ElConstMatrix_s A, const char* title );
ElError ElDisplayMatrix_d( ElConstMatrix_d A, const char* title );
ElError ElDisplayMatrix_c( ElConstMatrix_c A, const char* title );
ElError ElDisplayMatrix_z( ElConstMatrix_z A, const char* title );

ElError ElDisplayDistMatrix_s( ElConstDistMatrix_s A, const char* title );
ElError ElDisplayDistMatrix_d( ElConstDistMatrix_d A, const char* title );
ElError ElDisplayDistMatrix_c( ElConstDistMatrix_c A, const char* title );
ElError ElDisplayDistMatrix_z( ElConstDistMatrix_z A, const char* title );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_DISPLAY_C_H */
