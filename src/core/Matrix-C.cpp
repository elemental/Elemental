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

ElMatrix_s* ElMatrixCreate_s()
{
    try { return reinterpret_cast<ElMatrix_s*>( new Matrix<float>() ); }
    catch( std::exception& e ) { ReportException(e); }
}

ElMatrix_d* ElMatrixCreate_d()
{
    try { return reinterpret_cast<ElMatrix_d*>( new Matrix<double>() ); }
    catch( std::exception& e ) { ReportException(e); }
}

ElMatrix_c* ElMatrixCreate_c()
{
    try { return reinterpret_cast<ElMatrix_c*>( new Matrix<Complex<float>>() ); }
    catch( std::exception& e ) { ReportException(e); }
}

ElMatrix_z* ElMatrixCreate_z()
{
    try 
    { return reinterpret_cast<ElMatrix_z*>( new Matrix<Complex<double>>() ); }
    catch( std::exception& e ) { ReportException(e); }
}

void ElMatrixResize_s( ElMatrix_s* A, ElInt height, ElInt width )
{
    try 
    {
        auto* AMat = reinterpret_cast<Matrix<float>*>( A );
        AMat->Resize( height, width );
    } catch( std::exception& e ) { ReportException(e); }
}

void ElMatrixResize_d( ElMatrix_d* A, ElInt height, ElInt width )
{
    try 
    {
        auto* AMat = reinterpret_cast<Matrix<double>*>( A );
        AMat->Resize( height, width );
    } catch( std::exception& e ) { ReportException(e); }
}

void ElMatrixResize_c( ElMatrix_c* A, ElInt height, ElInt width )
{
    try 
    {
        auto* AMat = reinterpret_cast<Matrix<Complex<float>>*>( A );
        AMat->Resize( height, width );
    } catch( std::exception& e ) { ReportException(e); }
}

void ElMatrixResize_z( ElMatrix_z* A, ElInt height, ElInt width )
{
    try 
    {
        auto* AMat = reinterpret_cast<Matrix<Complex<double>>*>( A );
        AMat->Resize( height, width );
    } catch( std::exception& e ) { ReportException(e); }
}

void ElMatrixSet_s( ElMatrix_s* A, ElInt i, ElInt j, float alpha )
{
    try 
    {
        auto* AMat = reinterpret_cast<Matrix<float>*>( A );
        AMat->Set( i, j, alpha );
    } catch( std::exception& e ) { ReportException(e); }
}

void ElMatrixSet_d( ElMatrix_d* A, ElInt i, ElInt j, double alpha )
{
    try 
    {
        auto* AMat = reinterpret_cast<Matrix<double>*>( A );
        AMat->Set( i, j, alpha );
    } catch( std::exception& e ) { ReportException(e); }
}

void ElMatrixSet_c( ElMatrix_c* A, ElInt i, ElInt j, void* alpha )
{
    try 
    {
        auto* AMat = reinterpret_cast<Matrix<Complex<float>>*>( A );
        auto* alphaCpx = reinterpret_cast<Complex<float>*>( alpha );
        AMat->Set( i, j, *alphaCpx );
    } catch( std::exception& e ) { ReportException(e); }
}

void ElMatrixSet_z( ElMatrix_z* A, ElInt i, ElInt j, void* alpha )
{
    try 
    {
        auto* AMat = reinterpret_cast<Matrix<Complex<double>>*>( A );
        auto* alphaCpx = reinterpret_cast<Complex<double>*>( alpha );
        AMat->Set( i, j, *alphaCpx );
    } catch( std::exception& e ) { ReportException(e); }
}

} // extern "C"
