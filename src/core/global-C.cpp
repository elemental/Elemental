/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
#include "El-C.h"

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
{ EL_TRY( El::Initialize( *argc, *argv ) ) }

ElError ElFinalize()
{ EL_TRY( El::Finalize() ) }

ElError ElInitialized( bool* initialized )
{ EL_TRY( *initialized = El::Initialized() ) } 

ElError ElInput_b
( const char* name, const char* desc, bool defaultVal, bool* val )
{ EL_TRY( *val = El::Input(name,desc,defaultVal) ) } 

ElError ElInput_i
( const char* name, const char* desc, int defaultVal, int* val )
{ EL_TRY( *val = El::Input(name,desc,defaultVal) ) } 

ElError ElInput_I
( const char* name, const char* desc, ElInt defaultVal, ElInt* val )
{ EL_TRY( *val = El::Input(name,desc,defaultVal) ) } 

ElError ElInput_s
( const char* name, const char* desc, float defaultVal, float* val )
{ EL_TRY( *val = El::Input(name,desc,defaultVal) ) } 

ElError ElInput_d
( const char* name, const char* desc, double defaultVal, double* val )
{ EL_TRY( *val = El::Input(name,desc,defaultVal) ) } 

ElError ElInput_cstr
( const char* name, const char* desc, const char* defaultVal, const char** val )
{ EL_TRY( *val = El::Input(name,desc,defaultVal) ) } 

ElError ElProcessInput() 
{ EL_TRY( El::ProcessInput() ) }

ElError ElPrintInputReport() 
{ EL_TRY( El::PrintInputReport() ) } 

ElError ElBlocksize( ElInt* blocksize )
{ EL_TRY( *blocksize = El::Blocksize() ) }

ElError ElSetBlocksize( ElInt blocksize )
{ EL_TRY( El::SetBlocksize(blocksize) ) }

ElError ElPushBlocksizeStack( ElInt blocksize )
{ EL_TRY( El::PushBlocksizeStack(blocksize) ) }

ElError ElPopBlocksizeStack()
{ EL_TRY( El::PopBlocksizeStack() ) }

} // extern "C"
