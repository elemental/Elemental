/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
#include "El-C.h"
using namespace El;

#define CATCH \
  catch( std::bad_alloc& e ) \
  { ReportException(e); return EL_ALLOC_ERROR; } \
  catch( ArgException& e ) \
  { ReportException(e); return EL_ARG_ERROR; } \
  catch( std::logic_error& e ) \
  { ReportException(e); return EL_LOGIC_ERROR; } \
  catch( std::runtime_error& e ) \
  { ReportException(e); return EL_RUNTIME_ERROR; } \
  catch( std::exception& e ) \
  { ReportException(e); return EL_ERROR; }

extern "C" {

#define EL_UNIFORM_WRAPPER(SIG,T,TBASE) \
  ElError ElUniform ## SIG \
  ( El ## SIG A, ElInt m, ElInt n, T center, TBASE radius ) \
  { \
      try { Uniform( *Reinterpret(A), m, n, Reinterpret(center), radius ); } \
      CATCH \
      return EL_SUCCESS; \
  }
EL_UNIFORM_WRAPPER(Matrix_s,float,float)
EL_UNIFORM_WRAPPER(Matrix_d,double,double)
EL_UNIFORM_WRAPPER(Matrix_c,complex_float,float)
EL_UNIFORM_WRAPPER(Matrix_z,complex_double,double)
EL_UNIFORM_WRAPPER(DistMatrix_s,float,float)
EL_UNIFORM_WRAPPER(DistMatrix_d,double,double)
EL_UNIFORM_WRAPPER(DistMatrix_c,complex_float,float)
EL_UNIFORM_WRAPPER(DistMatrix_z,complex_double,double)
#undef EL_UNIFORM_WRAPPER

} // extern "C"
