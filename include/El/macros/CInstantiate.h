/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

#ifndef C_PROTO_INT
# define C_PROTO_INT(SIG,T) C_PROTO(SIG,T)
#endif

#ifndef C_PROTO_REAL 
# define C_PROTO_REAL(SIG,T) C_PROTO(SIG,T)
#endif
#ifndef C_PROTO_FLOAT
# define C_PROTO_FLOAT C_PROTO_REAL(s,float)
#endif
#ifndef C_PROTO_DOUBLE
# define C_PROTO_DOUBLE C_PROTO_REAL(d,double)
#endif

#ifndef C_PROTO_COMPLEX
# define C_PROTO_COMPLEX(SIG,T) C_PROTO(SIG,T)
#endif
#ifndef C_PROTO_COMPLEX_FLOAT
# define C_PROTO_COMPLEX_FLOAT C_PROTO_COMPLEX(c,Complex<float>)
#endif
#ifndef C_PROTO_COMPLEX_DOUBLE
# define C_PROTO_COMPLEX_DOUBLE C_PROTO_COMPLEX(z,Complex<double>)
#endif

#ifndef EL_NO_INT_PROTO
C_PROTO_INT(i,Int)
#endif

#ifndef EL_NO_REAL_PROTO
# if !defined(EL_NO_FLOAT_PROTO) && !defined(EL_DISABLE_FLOAT)
C_PROTO_FLOAT
# endif
# if !defined(EL_NO_DOUBLE_PROTO) && !defined(EL_DISABLE_DOUBLE)
C_PROTO_DOUBLE
# endif
#endif

#if !defined(EL_NO_COMPLEX_PROTO) && !defined(EL_DISABLE_COMPLEX)
# if !defined(EL_NO_COMPLEX_FLOAT_PROTO) && !defined(EL_DISABLE_FLOAT)
C_PROTO_COMPLEX_FLOAT
# endif
# if !defined(EL_NO_COMPLEX_DOUBLE_PROTO) && !defined(EL_DISABLE_DOUBLE)
C_PROTO_COMPLEX_DOUBLE
# endif
#endif
