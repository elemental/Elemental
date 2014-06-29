/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_CONFIGINTERNAL_H
#define EL_CONFIGINTERNAL_H

#include "El/config.h"

#if defined(EL_HAVE_VALGRIND)
# include "valgrind.h"
# define EL_RUNNING_ON_VALGRIND RUNNING_ON_VALGRIND
#else
# define EL_RUNNING_ON_VALGRIND 0
#endif

#endif /* EL_CONFIGINTERNAL_H */
