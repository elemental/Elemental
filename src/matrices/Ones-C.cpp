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

#define EL_ONES_WRAPPER(SIG) \
  ElError ElOnes ## SIG ( El ## SIG A, ElInt m, ElInt n ) \
  { \
      try { Ones( *Reinterpret(A), m, n ); } \
      CATCH \
      return EL_SUCCESS; \
  }
EL_ONES_WRAPPER(Matrix_s)
EL_ONES_WRAPPER(Matrix_d)
EL_ONES_WRAPPER(Matrix_c)
EL_ONES_WRAPPER(Matrix_z)
EL_ONES_WRAPPER(DistMatrix_s)
EL_ONES_WRAPPER(DistMatrix_d)
EL_ONES_WRAPPER(DistMatrix_c)
EL_ONES_WRAPPER(DistMatrix_z)
#undef EL_ONES_WRAPPER

} // extern "C"
