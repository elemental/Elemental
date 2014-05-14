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

#define RCDDM_s_const(A) RC(const DynamicDistMatrix<float          >*,A)
#define RCDDM_d_const(A) RC(const DynamicDistMatrix<double         >*,A)
#define RCDDM_c_const(A) RC(const DynamicDistMatrix<Complex<float >>*,A)
#define RCDDM_z_const(A) RC(const DynamicDistMatrix<Complex<double>>*,A)

#define CATCH catch( std::exception& e ) { ReportException(e); }

extern "C" {

// Matrix
// ======

void ElDisplayMatrix_s( const ElMatrix_s* AHandle, const char* title )
{
    try { Display( *RCM_s_const(AHandle), std::string(title) ); }
    CATCH
}

void ElDisplayMatrix_d( const ElMatrix_d* AHandle, const char* title )
{
    try { Display( *RCM_d_const(AHandle), std::string(title) ); }
    CATCH
}

void ElDisplayMatrix_c( const ElMatrix_c* AHandle, const char* title )
{
    try { Display( *RCM_c_const(AHandle), std::string(title) ); }
    CATCH
}

void ElDisplayMatrix_z( const ElMatrix_z* AHandle, const char* title )
{
    try { Display( *RCM_z_const(AHandle), std::string(title) ); }
    CATCH
}

// DynamicDistMatrix
// =================

void ElDisplayDistMatrix_s( const ElDistMatrix_s* AHandle, const char* title )
{
    try { Display( *RCDDM_s_const(AHandle), std::string(title) ); }
    CATCH
}

void ElDisplayDistMatrix_d( const ElDistMatrix_d* AHandle, const char* title )
{
    try { Display( *RCDDM_d_const(AHandle), std::string(title) ); }
    CATCH
}

void ElDisplayDistMatrix_c( const ElDistMatrix_c* AHandle, const char* title )
{
    try { Display( *RCDDM_c_const(AHandle), std::string(title) ); }
    CATCH
}

void ElDisplayDistMatrix_z( const ElDistMatrix_z* AHandle, const char* title )
{
    try { Display( *RCDDM_z_const(AHandle), std::string(title) ); }
    CATCH
}

} // extern "C"
