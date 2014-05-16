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

#define RC(TYPE,INPUT) reinterpret_cast<TYPE>(INPUT)

#define RCM_s_const(A) RC(const Matrix<float          >*,A)
#define RCM_d_const(A) RC(const Matrix<double         >*,A)
#define RCM_c_const(A) RC(const Matrix<Complex<float >>*,A)
#define RCM_z_const(A) RC(const Matrix<Complex<double>>*,A)

#define RCADM_s_const(A) RC(const AbstractDistMatrix<float          >*,A)
#define RCADM_d_const(A) RC(const AbstractDistMatrix<double         >*,A)
#define RCADM_c_const(A) RC(const AbstractDistMatrix<Complex<float >>*,A)
#define RCADM_z_const(A) RC(const AbstractDistMatrix<Complex<double>>*,A)

#define CATCH catch( std::exception& e ) { ReportException(e); }

extern "C" {

// Matrix
// ======

void ElPrintMatrix_s( ElConstMatrix_s AHandle, const char* title )
{
    try { Print( *RCM_s_const(AHandle), std::string(title) ); }
    CATCH
}

void ElPrintMatrix_d( ElConstMatrix_d AHandle, const char* title )
{
    try { Print( *RCM_d_const(AHandle), std::string(title) ); }
    CATCH
}

void ElPrintMatrix_c( ElConstMatrix_c AHandle, const char* title )
{
    try { Print( *RCM_c_const(AHandle), std::string(title) ); }
    CATCH
}

void ElPrintMatrix_z( ElConstMatrix_z AHandle, const char* title )
{
    try { Print( *RCM_z_const(AHandle), std::string(title) ); }
    CATCH
}

// AbstractDistMatrix
// ==================

void ElPrintDistMatrix_s( ElConstDistMatrix_s AHandle, const char* title )
{
    try { Print( *RCADM_s_const(AHandle), std::string(title) ); }
    CATCH
}

void ElPrintDistMatrix_d( ElConstDistMatrix_d AHandle, const char* title )
{
    try { Print( *RCADM_d_const(AHandle), std::string(title) ); }
    CATCH
}

void ElPrintDistMatrix_c( ElConstDistMatrix_c AHandle, const char* title )
{
    try { Print( *RCADM_c_const(AHandle), std::string(title) ); }
    CATCH
}

void ElPrintDistMatrix_z( ElConstDistMatrix_z AHandle, const char* title )
{
    try { Print( *RCADM_z_const(AHandle), std::string(title) ); }
    CATCH
}

} // extern "C"
