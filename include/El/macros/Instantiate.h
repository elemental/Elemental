/*
   Copyright (c) 2009-2016, Jack Poulson
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

#if defined(EL_HAVE_QD) && defined(EL_ENABLE_DOUBLEDOUBLE)
#ifndef PROTO_DOUBLEDOUBLE
# define PROTO_DOUBLEDOUBLE PROTO_REAL(DoubleDouble)
#endif
#endif

#if defined(EL_HAVE_QD) && defined(EL_ENABLE_QUADDOUBLE)
#ifndef PROTO_QUADDOUBLE
# define PROTO_QUADDOUBLE PROTO_REAL(QuadDouble)
#endif
#endif

#if defined(EL_HAVE_QUAD) && defined(EL_ENABLE_QUAD)
#ifndef PROTO_QUAD
# define PROTO_QUAD PROTO_REAL(Quad)
#endif
#endif

#if defined(EL_HAVE_MPC) && defined(EL_ENABLE_BIGINT)
#ifndef PROTO_BIGINT
# define PROTO_BIGINT PROTO_INT(BigInt)
#endif
#endif

#if defined(EL_HAVE_MPC) && defined(EL_ENABLE_BIGFLOAT)
#ifndef PROTO_BIGFLOAT
# define PROTO_BIGFLOAT PROTO_REAL(BigFloat)
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
#if defined(EL_ENABLE_BIGINT) && defined(EL_HAVE_MPC)
PROTO_BIGINT
#endif
#endif

#ifndef EL_NO_REAL_PROTO
# if !defined(EL_NO_FLOAT_PROTO)
PROTO_FLOAT
# endif
# if !defined(EL_NO_DOUBLE_PROTO)
PROTO_DOUBLE
# endif
#if defined(EL_ENABLE_DOUBLEDOUBLE) && defined(EL_HAVE_QD)
PROTO_DOUBLEDOUBLE
#endif
#if defined(EL_ENABLE_QUADDOUBLE) && defined(EL_HAVE_QD)
PROTO_QUADDOUBLE
#endif
#if defined(EL_ENABLE_QUAD) && defined(EL_HAVE_QUAD)
PROTO_QUAD
#endif
#if defined(EL_ENABLE_BIGFLOAT) && defined(EL_HAVE_MPC)
PROTO_BIGFLOAT
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

#undef PROTO
#undef PROTO_INT
#undef PROTO_BIGINT

#undef PROTO_REAL
#undef PROTO_FLOAT
#undef PROTO_DOUBLE
#undef PROTO_DOUBLEDOUBLE
#undef PROTO_QUADDOUBLE
#undef PROTO_QUAD
#undef PROTO_BIGFLOAT

#undef PROTO_COMPLEX
#undef PROTO_COMPLEX_FLOAT
#undef PROTO_COMPLEX_DOUBLE
#undef PROTO_COMPLEX_QUAD

#undef EL_ENABLE_DOUBLEDOUBLE
#undef EL_ENABLE_QUADDOUBLE
#undef EL_ENABLE_QUAD
#undef EL_ENABLE_BIGINT
#undef EL_ENABLE_BIGFLOAT

#undef EL_NO_INT_PROTO
#undef EL_NO_REAL_PROTO
#undef EL_NO_FLOAT_PROTO
#undef EL_NO_DOUBLE_PROTO
#undef EL_NO_COMPLEX_PROTO
#undef EL_NO_COMPLEX_FLOAT_PROTO
#undef EL_NO_COMPLEX_DOUBLE_PROTO
