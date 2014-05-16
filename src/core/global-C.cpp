/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"
#include "El-C.h"

#define CATCH \
  catch( std::bad_alloc& e ) \
  { El::ReportException(e); return EL_ALLOC_ERROR; } \
  catch( El::ArgException& e ) \
  { El::ReportException(e); return EL_ARG_ERROR; } \
  catch( std::logic_error& e ) \
  { El::ReportException(e); return EL_LOGIC_ERROR; } \
  catch( std::runtime_error& e ) \
  { El::ReportException(e); return EL_RUNTIME_ERROR; } \
  catch( std::exception& e ) \
  { El::ReportException(e); return EL_ERROR; }

extern "C" {

ElError ElPrintVersion( FILE* stream )
{
    // There does not seem to be a portable means of converting C-style
    // filehandles to C++ filestreams, so we will simply reproduce the 
    // functionality of El::PrintVersion.    
    fprintf( stream, "Elemental version information:\n"
                     " Git revision: %s\n"
                     " Version:      %s.%s\n"
                     " Build type:   %s\n\n",
             EL_GIT_SHA1, EL_VERSION_MAJOR, EL_VERSION_MINOR, 
             EL_CMAKE_BUILD_TYPE );
    return EL_SUCCESS;
}

const char* ElErrorString( ElError error )
{
    if( error == EL_SUCCESS )
    {
        const char* errString = "EL_SUCCESS";
        return errString;
    }
    else if( error == EL_ALLOC_ERROR )
    {
        const char* errString = "EL_ALLOC_ERROR";
        return errString;
    }
    else if( error == EL_OUT_OF_BOUNDS_ERROR )
    {
        const char* errString = "EL_OUT_OF_BOUNDS_ERROR";
        return errString;
    }
    else if( error == EL_ARG_ERROR )
    {
        const char* errString = "EL_ARG_ERROR";
        return errString;
    }
    else if( error == EL_LOGIC_ERROR )
    {
        const char* errString = "EL_LOGIC_ERROR";
        return errString;
    }
    else if( error == EL_RUNTIME_ERROR )
    {
        const char* errString = "EL_RUNTIME_ERROR";
        return errString;
    }
    else
    {
        const char* errString = "EL_ERROR";
        return errString;
    }
}

// TODO: ElPrintConfig
// TODO: ElPrintCCompilerInfo
// TODO: ElPrintCxxCompilerInfo

ElError ElInitialize( int* argc, char*** argv )
{
    try { El::Initialize( *argc, *argv ); }
    CATCH
    return EL_SUCCESS;
}

ElError ElFinalize()
{
    try { El::Finalize(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElInitialized( bool* initialized )
{ 
    *initialized = El::Initialized(); 
    return EL_SUCCESS;
}

ElError ElInput_b
( const char* name, const char* desc, bool defaultVal, bool* val )
{ 
    try { *val = El::Input(name,desc,defaultVal); }
    CATCH
    return EL_SUCCESS;
}

ElError ElInput_i
( const char* name, const char* desc, int defaultVal, int* val )
{ 
    try { *val = El::Input(name,desc,defaultVal); }
    CATCH
    return EL_SUCCESS;
}

ElError ElInput_I
( const char* name, const char* desc, ElInt defaultVal, ElInt* val )
{ 
    try { *val = El::Input(name,desc,defaultVal); }
    CATCH
    return EL_SUCCESS;
}

ElError ElInput_s
( const char* name, const char* desc, float defaultVal, float* val )
{ 
    try { *val = El::Input(name,desc,defaultVal); }
    CATCH
    return EL_SUCCESS;
}

ElError ElInput_d
( const char* name, const char* desc, double defaultVal, double* val )
{ 
    try { *val = El::Input(name,desc,defaultVal); }
    CATCH
    return EL_SUCCESS;
}

ElError ElInput_cstr
( const char* name, const char* desc, const char* defaultVal, const char** val )
{ 
    try { *val = El::Input(name,desc,defaultVal); }
    CATCH
    return EL_SUCCESS;
}

ElError ElProcessInput()
{
    try { El::ProcessInput(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElPrintInputReport()
{ 
    El::PrintInputReport(); 
    return EL_SUCCESS;
}

ElError ElBlocksize( ElInt* blocksize )
{ 
    try { *blocksize = El::Blocksize(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElSetBlocksize( ElInt blocksize )
{
    try { El::SetBlocksize(blocksize); }
    CATCH
    return EL_SUCCESS;
}

ElError ElPushBlocksizeStack( ElInt blocksize )
{
    try { El::PushBlocksizeStack(blocksize); }
    CATCH
    return EL_SUCCESS;
}

ElError ElPopBlocksizeStack()
{
    try { El::PopBlocksizeStack(); }
    CATCH
    return EL_SUCCESS;
}

} // extern "C"
