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

// B = A
// -----
ElError ElCopyDistMatrix_s
( ElConstDistMatrix_s AHandle, ElDistMatrix_s BHandle )
{
    try { Copy( *Reinterpret(AHandle), *Reinterpret(BHandle) ); }
    CATCH
    return EL_SUCCESS;
}

ElError ElCopyDistMatrix_d
( ElConstDistMatrix_d AHandle, ElDistMatrix_d BHandle )
{
    try { Copy( *Reinterpret(AHandle), *Reinterpret(BHandle) ); }
    CATCH
    return EL_SUCCESS;
}

ElError ElCopyDistMatrix_c
( ElConstDistMatrix_c AHandle, ElDistMatrix_c BHandle )
{
    try { Copy( *Reinterpret(AHandle), *Reinterpret(BHandle) ); }
    CATCH
    return EL_SUCCESS;
}

ElError ElCopyDistMatrix_z
( ElConstDistMatrix_z AHandle, ElDistMatrix_z BHandle )
{
    try { Copy( *Reinterpret(AHandle), *Reinterpret(BHandle) ); }
    CATCH
    return EL_SUCCESS;
}

// B = A^T
// -------
ElError ElTransposeDistMatrix_s
( ElConstDistMatrix_s AHandle, ElDistMatrix_s BHandle )
{
    try { Transpose( *Reinterpret(AHandle), *Reinterpret(BHandle), false ); }
    CATCH
    return EL_SUCCESS;
}

ElError ElTransposeDistMatrix_d
( ElConstDistMatrix_d AHandle, ElDistMatrix_d BHandle )
{
    try { Transpose( *Reinterpret(AHandle), *Reinterpret(BHandle), false ); }
    CATCH
    return EL_SUCCESS;
}

ElError ElTransposeDistMatrix_c
( ElConstDistMatrix_c AHandle, ElDistMatrix_c BHandle )
{
    try { Transpose( *Reinterpret(AHandle), *Reinterpret(BHandle), false ); }
    CATCH
    return EL_SUCCESS;
}

ElError ElTransposeDistMatrix_z
( ElConstDistMatrix_z AHandle, ElDistMatrix_z BHandle )
{
    try { Transpose( *Reinterpret(AHandle), *Reinterpret(BHandle), false ); }
    CATCH
    return EL_SUCCESS;
}

// B = A^H
// -------
ElError ElAdjointDistMatrix_s
( ElConstDistMatrix_s AHandle, ElDistMatrix_s BHandle )
{
    try { Adjoint( *Reinterpret(AHandle), *Reinterpret(BHandle) ); }
    CATCH
    return EL_SUCCESS;
}

ElError ElAdjointDistMatrix_d
( ElConstDistMatrix_d AHandle, ElDistMatrix_d BHandle )
{
    try { Adjoint( *Reinterpret(AHandle), *Reinterpret(BHandle) ); }
    CATCH
    return EL_SUCCESS;
}

ElError ElAdjointDistMatrix_c
( ElConstDistMatrix_c AHandle, ElDistMatrix_c BHandle )
{
    try { Adjoint( *Reinterpret(AHandle), *Reinterpret(BHandle) ); }
    CATCH
    return EL_SUCCESS;
}

ElError ElAdjointDistMatrix_z
( ElConstDistMatrix_z AHandle, ElDistMatrix_z BHandle )
{
    try { Adjoint( *Reinterpret(AHandle), *Reinterpret(BHandle) ); }
    CATCH
    return EL_SUCCESS;
}

} // extern "C"
