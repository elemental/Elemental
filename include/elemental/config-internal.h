/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef ELEM_CONFIGINTERNAL_H
#define ELEM_CONFIGINTERNAL_H

#include "elemental/config.h"

#if defined(ELEM_HAVE_VALGRIND)
# include "valgrind.h"
# define ELEM_RUNNING_ON_VALGRIND RUNNING_ON_VALGRIND
#else
# define ELEM_RUNNING_ON_VALGRIND 0
#endif

#endif /* ELEM_CONFIGINTERNAL_H */
