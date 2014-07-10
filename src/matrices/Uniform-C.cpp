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

#define EL_TRY(payload) \
  try { payload; } CATCH \
  return EL_SUCCESS;

#define CREFLECT(T) typename CReflect<T>::type

extern "C" {

#define C_PROTO(SIG,T) \
  ElError ElUniformMatrix_ ## SIG \
  ( ElMatrix_ ## SIG A, \
    ElInt m, ElInt n, CREFLECT(T) center, CREFLECT(Base<T>) radius ) \
  { EL_TRY( Uniform( *Reinterpret(A), m, n, Reinterpret(center), radius ) ) }
#include "El/macros/CInstantiate.h"
#undef C_PROTO

#define C_PROTO(SIG,T) \
  ElError ElUniformDistMatrix_ ## SIG \
  ( ElDistMatrix_ ## SIG A, \
    ElInt m, ElInt n, CREFLECT(T) center, CREFLECT(Base<T>) radius ) \
  { EL_TRY( Uniform( *Reinterpret(A), m, n, Reinterpret(center), radius ) ) }
#include "El/macros/CInstantiate.h"
#undef C_PROTO

} // extern "C"
