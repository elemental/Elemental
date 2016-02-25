/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_IMPORTS_VALGRIND_HPP
#define EL_IMPORTS_VALGRIND_HPP

#ifdef EL_HAVE_VALGRIND

# include "valgrind/valgrind.h"
# define EL_RUNNING_ON_VALGRIND RUNNING_ON_VALGRIND

#else

# define EL_RUNNING_ON_VALGRIND 0

#endif // ifdef EL_HAVE_VALGRIND

#endif // ifndef EL_IMPORTS_VALGRIND_HPP
