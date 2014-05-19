/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"
#include "El/matrices/Ones.hpp"
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

ElError ElOnesMatrix_s( ElMatrix_s A, ElInt m, ElInt n )
{
    try { Ones( *Reinterpret(A), m, n ); }
    CATCH
    return EL_SUCCESS;
}

ElError ElOnesMatrix_d( ElMatrix_d A, ElInt m, ElInt n )
{
    try { Ones( *Reinterpret(A), m, n ); }
    CATCH
    return EL_SUCCESS;
}

ElError ElOnesMatrix_c( ElMatrix_c A, ElInt m, ElInt n )
{
    try { Ones( *Reinterpret(A), m, n ); }
    CATCH
    return EL_SUCCESS;
}

ElError ElOnesMatrix_z( ElMatrix_z A, ElInt m, ElInt n )
{
    try { Ones( *Reinterpret(A), m, n ); }
    CATCH
    return EL_SUCCESS;
}

ElError ElOnesDistMatrix_s( ElDistMatrix_s A, ElInt m, ElInt n )
{
    try { Ones( *Reinterpret(A), m, n ); }
    CATCH
    return EL_SUCCESS;
}

ElError ElOnesDistMatrix_d( ElDistMatrix_d A, ElInt m, ElInt n )
{
    try { Ones( *Reinterpret(A), m, n ); }
    CATCH
    return EL_SUCCESS;
}

ElError ElOnesDistMatrix_c( ElDistMatrix_c A, ElInt m, ElInt n )
{
    try { Ones( *Reinterpret(A), m, n ); }
    CATCH
    return EL_SUCCESS;
}

ElError ElOnesDistMatrix_z( ElDistMatrix_z A, ElInt m, ElInt n )
{
    try { Ones( *Reinterpret(A), m, n ); }
    CATCH
    return EL_SUCCESS;
}

} // extern "C"
