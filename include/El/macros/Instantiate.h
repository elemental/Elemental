/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef PROTO_INT
# define PROTO_INT(T) PROTO(T)
#endif

#ifndef PROTO_REAL 
# define PROTO_REAL(T) PROTO(T)
#endif
#ifndef PROTO_FLOAT
# define PROTO_FLOAT PROTO_REAL(float)
#endif
#ifndef PROTO_DOUBLE
# define PROTO_DOUBLE PROTO_REAL(double)
#endif

#if defined(EL_HAVE_QUAD) && defined(EL_ENABLE_QUAD)
#ifndef PROTO_QUAD
# define PROTO_QUAD PROTO_REAL(Quad)
#endif
#endif

#ifndef PROTO_COMPLEX
# define PROTO_COMPLEX(T) PROTO(T)
#endif
#ifndef PROTO_COMPLEX_FLOAT
# define PROTO_COMPLEX_FLOAT PROTO_COMPLEX(Complex<float>)
#endif
#ifndef PROTO_COMPLEX_DOUBLE
# define PROTO_COMPLEX_DOUBLE PROTO_COMPLEX(Complex<double>)
#endif

#if defined(EL_HAVE_QUAD) && defined(EL_ENABLE_QUAD)
#ifndef PROTO_COMPLEX_QUAD
# define PROTO_COMPLEX_QUAD PROTO_COMPLEX(Complex<Quad>)
#endif
#endif

#ifndef EL_NO_INT_PROTO
PROTO_INT(Int)
#endif

#ifndef EL_NO_REAL_PROTO
# if !defined(EL_NO_FLOAT_PROTO)
PROTO_FLOAT
# endif
# if !defined(EL_NO_DOUBLE_PROTO)
PROTO_DOUBLE
# endif
#if defined(EL_ENABLE_QUAD) && defined(EL_HAVE_QUAD)
PROTO_QUAD
#endif
#endif

#if !defined(EL_NO_COMPLEX_PROTO)
# if !defined(EL_NO_COMPLEX_FLOAT_PROTO)
PROTO_COMPLEX_FLOAT
# endif
# if !defined(EL_NO_COMPLEX_DOUBLE_PROTO)
PROTO_COMPLEX_DOUBLE
# endif
#if defined(EL_ENABLE_QUAD) && defined(EL_HAVE_QUAD)
PROTO_COMPLEX_QUAD
#endif
#endif
