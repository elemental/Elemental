/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_LITE_HPP
#define EL_LITE_HPP

#include "El/include-paths.hpp"

#include "El/config.h"
#ifdef EL_HAVE_F90_INTERFACE
# include "El/FCMangle.h"
#endif

#include "El/core.hpp"
#include "El/blas-like.hpp"

#include "El/io.hpp"

#include "El/lapack-like.hpp"
#include "El/optimization.hpp"

#endif // ifndef EL_LITE_HPP
