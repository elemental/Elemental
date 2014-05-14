/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"
#include "El-C.h"

#define CATCH catch( std::exception& e ) { El::ReportException(e); }

extern "C" {

void ElPrintVersion( FILE* stream )
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
}

// TODO: ElPrintConfig
// TODO: ElPrintCCompilerInfo
// TODO: ElPrintCxxCompilerInfo

void ElInitialize( int* argc, char*** argv )
{
    try { El::Initialize( *argc, *argv ); }
    CATCH
}

void ElFinalize()
{
    try { El::Finalize(); }
    CATCH
}

bool ElInitialized()
{ return El::Initialized(); }

bool ElInput_b( const char* name, const char* desc, bool defaultVal )
{ return El::Input(name,desc,defaultVal); }

int ElInput_i( const char* name, const char* desc, int defaultVal )
{ return El::Input(name,desc,defaultVal); }

ElInt ElInput_I( const char* name, const char* desc, ElInt defaultVal )
{ return El::Input(name,desc,defaultVal); }

float ElInput_s( const char* name, const char* desc, float defaultVal )
{ return El::Input(name,desc,defaultVal); }

double ElInput_d( const char* name, const char* desc, double defaultVal )
{ return El::Input(name,desc,defaultVal); }

const char* ElInput_cstr
( const char* name, const char* desc, const char* defaultVal )
{ return El::Input(name,desc,defaultVal); }

void ElProcessInput()
{
    try { El::ProcessInput(); }
    CATCH
}

void ElPrintInputReport()
{ El::PrintInputReport(); }

ElInt ElBlocksize()
{ return El::Blocksize(); }

void ElSetBlocksize( ElInt blocksize )
{
    try { El::SetBlocksize(blocksize); }
    CATCH
}

void ElPushBlocksizeStack( ElInt blocksize )
{
    try { El::PushBlocksizeStack(blocksize); }
    CATCH
}

void ElPopBlocksizeStack()
{
    try { El::PopBlocksizeStack(); }
    CATCH
}

} // extern "C"
