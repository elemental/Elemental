/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"
#include "El/matrices/Uniform.hpp"
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

ElError ElUniformMatrix_s
( ElMatrix_s A, ElInt m, ElInt n, float center, float radius )
{
    try { Uniform( *Reinterpret(A), m, n, center, radius ); }
    CATCH
    return EL_SUCCESS;
}

ElError ElUniformMatrix_d
( ElMatrix_d A, ElInt m, ElInt n, double center, double radius )
{
    try { Uniform( *Reinterpret(A), m, n, center, radius ); }
    CATCH
    return EL_SUCCESS;
}

ElError ElUniformMatrix_c
( ElMatrix_c A, ElInt m, ElInt n, complex_float center, float radius )
{
    try { Uniform( *Reinterpret(A), m, n, 
                   Complex<float>(center.real,center.imag), radius ); }
    CATCH
    return EL_SUCCESS;
}

ElError ElUniformMatrix_z
( ElMatrix_z A, ElInt m, ElInt n, complex_double center, double radius )
{
    try { Uniform( *Reinterpret(A), m, n, 
                   Complex<double>(center.real,center.imag), radius ); }
    CATCH
    return EL_SUCCESS;
}

ElError ElUniformDistMatrix_s
( ElDistMatrix_s A, ElInt m, ElInt n, float center, float radius )
{
    try { Uniform( *Reinterpret(A), m, n, center, radius ); }
    CATCH
    return EL_SUCCESS;
}

ElError ElUniformDistMatrix_d
( ElDistMatrix_d A, ElInt m, ElInt n, double center, double radius )
{
    try { Uniform( *Reinterpret(A), m, n, center, radius ); }
    CATCH
    return EL_SUCCESS;
}

ElError ElUniformDistMatrix_c
( ElDistMatrix_c A, ElInt m, ElInt n, complex_float center, float radius )
{
    try { Uniform( *Reinterpret(A), m, n, 
                   Complex<float>(center.real,center.imag), radius ); }
    CATCH
    return EL_SUCCESS;
}

ElError ElUniformDistMatrix_z
( ElDistMatrix_z A, ElInt m, ElInt n, complex_double center, double radius )
{
    try { Uniform( *Reinterpret(A), m, n, 
                   Complex<double>(center.real,center.imag), radius ); }
    CATCH
    return EL_SUCCESS;
}

} // extern "C"
