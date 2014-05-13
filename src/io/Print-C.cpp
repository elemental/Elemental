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

extern "C" {

void ElPrintMatrix_s( const ElMatrix_s* A, const char* title )
{
    try 
    { 
        const auto* AMat = reinterpret_cast<const Matrix<float>*>(A);
        Print( *AMat, std::string(title) );
    } catch( std::exception& e ) { ReportException(e); }
}

void ElPrintMatrix_d( const ElMatrix_d* A, const char* title )
{
    try 
    { 
        const auto* AMat = reinterpret_cast<const Matrix<double>*>(A);
        Print( *AMat, std::string(title) );
    } catch( std::exception& e ) { ReportException(e); }
}

void ElPrintMatrix_c( const ElMatrix_c* A, const char* title )
{
    try 
    { 
        const auto* AMat = reinterpret_cast<const Matrix<Complex<float>>*>(A);
        Print( *AMat, std::string(title) );
    } catch( std::exception& e ) { ReportException(e); }
}

void ElPrintMatrix_z( const ElMatrix_z* A, const char* title )
{
    try 
    { 
        const auto* AMat = reinterpret_cast<const Matrix<Complex<double>>*>(A);
        Print( *AMat, std::string(title) );
    } catch( std::exception& e ) { ReportException(e); }
}

} // extern "C"
