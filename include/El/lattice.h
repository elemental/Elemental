/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_LATTICE_C_H
#define EL_LATTICE_C_H

#ifdef __cplusplus
extern "C" {
#endif

EL_EXPORT ElError ElLLL_s
( ElMatrix_s B, float delta, float eta, float theta, float innerTol );
EL_EXPORT ElError ElLLL_d
( ElMatrix_d B, double delta, double eta, double theta, double innerTol );

#ifdef __cplusplus
}
#endif

#endif // ifndef EL_LATTICE_C_H
