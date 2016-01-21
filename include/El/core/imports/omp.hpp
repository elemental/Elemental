/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_IMPORTS_OMP_HPP
#define EL_IMPORTS_OMP_HPP

#ifdef EL_HYBRID
# include <omp.h>
# define EL_PARALLEL_FOR _Pragma("omp parallel for")
# ifdef EL_HAVE_OMP_COLLAPSE
#  define EL_PARALLEL_FOR_COLLAPSE2 _Pragma("omp parallel for collapse(2)")
# else
#  define EL_PARALLEL_FOR_COLLAPSE2 EL_PARALLEL_FOR
# endif
#else
# define EL_PARALLEL_FOR 
# define EL_PARALLEL_FOR_COLLAPSE2
#endif

#ifdef EL_AVOID_OMP_FMA
# define EL_FMA_PARALLEL_FOR 
#else
# define EL_FMA_PARALLEL_FOR EL_PARALLEL_FOR
#endif
#ifdef EL_PARALLELIZE_INNER_LOOPS
# define EL_INNER_PARALLEL_FOR           EL_PARALLEL_FOR
# define EL_INNER_PARALLEL_FOR_COLLAPSE2 EL_PARALLEL_FOR_COLLAPSE2
# define EL_OUTER_PARALLEL_FOR 
# define EL_OUTER_PARALLEL_FOR_COLLAPSE2
#else
# define EL_INNER_PARALLEL_FOR
# define EL_INNER_PARALLEL_FOR_COLLAPSE2
# define EL_OUTER_PARALLEL_FOR           EL_PARALLEL_FOR
# define EL_OUTER_PARALLEL_FOR_COLLAPSE2 EL_PARALLEL_FOR_COLLAPSE2
#endif

#endif // ifndef EL_IMPORTS_OMP_HPP
