/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"
#include "El-C.h"
using namespace El;

#define CATCH \
  catch( std::bad_alloc& e ) \
  { ReportException(e); return EL_ALLOC_ERROR; } \
  catch( std::logic_error& e ) \
  { ReportException(e); return EL_LOGIC_ERROR; } \
  catch( std::runtime_error& e ) \
  { ReportException(e); return EL_RUNTIME_ERROR; } \
  catch( std::exception& e ) \
  { ReportException(e); return EL_ERROR; }

extern "C" {

// Matrix
// ======

ElError ElSpyMatrix_s( ElConstMatrix_s AHandle, const char* title, float tol )
{
    try { Spy( *Reinterpret(AHandle), std::string(title), tol ); }
    CATCH
    return EL_SUCCESS;
}

ElError ElSpyMatrix_d( ElConstMatrix_d AHandle, const char* title, double tol )
{
    try { Spy( *Reinterpret(AHandle), std::string(title), tol ); }
    CATCH
    return EL_SUCCESS;
}

ElError ElSpyMatrix_c( ElConstMatrix_c AHandle, const char* title, float tol )
{
    try { Spy( *Reinterpret(AHandle), std::string(title), tol ); }
    CATCH
    return EL_SUCCESS;
}

ElError ElSpyMatrix_z( ElConstMatrix_z AHandle, const char* title, double tol )
{
    try { Spy( *Reinterpret(AHandle), std::string(title), tol ); }
    CATCH
    return EL_SUCCESS;
}

// AbstractDistMatrix
// ==================

ElError ElSpyDistMatrix_s
( ElConstDistMatrix_s AHandle, const char* title, float tol )
{
    try { Spy( *Reinterpret(AHandle), std::string(title), tol ); }
    CATCH
    return EL_SUCCESS;
}

ElError ElSpyDistMatrix_d
( ElConstDistMatrix_d AHandle, const char* title, double tol )
{
    try { Spy( *Reinterpret(AHandle), std::string(title), tol ); }
    CATCH
    return EL_SUCCESS;
}

ElError ElSpyDistMatrix_c
( ElConstDistMatrix_c AHandle, const char* title, float tol )
{
    try { Spy( *Reinterpret(AHandle), std::string(title), tol ); }
    CATCH
    return EL_SUCCESS;
}

ElError ElSpyDistMatrix_z
( ElConstDistMatrix_z AHandle, const char* title, double tol )
{
    try { Spy( *Reinterpret(AHandle), std::string(title), tol ); }
    CATCH
    return EL_SUCCESS;
}

} // extern "C"
