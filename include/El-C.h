/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_C_H
#define EL_C_H

#include "mpi.h"

/* TODO: A better include structure for the C interface */
#include "El/config.h"
#include "El/core/types-C.h"
#include "El/core/environment-C.h"
#include "El/core/Grid-C.h"
#include "El/core/Matrix-C.h"
#include "El/core/DistMatrix-C.h"

#include "El/io-C.h"

#include "El/blas-like/level1-C.h"
#include "El/blas-like/level2-C.h"
#include "El/blas-like/level3-C.h"
#include "El/lapack-like/condense-C.h"
#include "El/lapack-like/factor-C.h"
#include "El/lapack-like/solve-C.h"
#include "El/lapack-like/funcs-C.h"
#include "El/lapack-like/decomp-C.h"
#include "El/lapack-like/props-C.h"

#include "El/matrices-C.h"

#ifdef __cplusplus
#include "El/CReflect.hpp"
#endif

#endif /* ifndef EL_C_H */
