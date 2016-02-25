/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
#include "El.h"

extern "C" {

ElError ElPrintVersion( FILE* stream )
{
    // There does not seem to be a portable means of converting C-style
    // filehandles to C++ filestreams, so we will simply reproduce the 
    // functionality of El::PrintVersion.    
    fprintf
    ( stream, 
      "Elemental version information:\n"
      "  Git revision: %s\n"
      "  Version:      %s.%s\n"
      "  Build type:   %s\n\n",
      EL_GIT_SHA1, EL_VERSION_MAJOR, EL_VERSION_MINOR, EL_CMAKE_BUILD_TYPE );
    return EL_SUCCESS;
}

const char* ElErrorString( ElError error )
{
    if( error == EL_SUCCESS )
    {
        static const char* successString = "EL_SUCCESS";
        return successString;
    }
    else if( error == EL_ALLOC_ERROR )
    {
        static const char* allocString = "EL_ALLOC_ERROR";
        return allocString;
    }
    else if( error == EL_OUT_OF_BOUNDS_ERROR )
    {
        static const char* oobString = "EL_OUT_OF_BOUNDS_ERROR";
        return oobString;
    }
    else if( error == EL_ARG_ERROR )
    {
        static const char* argString = "EL_ARG_ERROR";
        return argString;
    }
    else if( error == EL_LOGIC_ERROR )
    {
        static const char* logicString = "EL_LOGIC_ERROR";
        return logicString;
    }
    else if( error == EL_RUNTIME_ERROR )
    {
        static const char* runtimeString = "EL_RUNTIME_ERROR";
        return runtimeString;
    }
    else
    {
        static const char* errString = "EL_ERROR";
        return errString;
    }
}

ElError ElPrintConfig( FILE* stream )
{
    // There does not seem to be a portable means of converting C-style
    // filehandles to C++ filestreams, so we will simply reproduce the 
    // functionality of El::PrintConfig.    
    fprintf
    ( stream, 
      "Elemental configuration:\n"
      "  Math libraries:               %s\n"
#ifdef EL_HAVE_FLA_BSVD
      "  Have FLAME bidiagonal SVD:    YES\n"
#else
      "  Have FLAME bidiagonal SVD:    NO\n"
#endif
#ifdef EL_HYBRID
      "  Hybrid mode:                  YES\n"
#else
      "  Hybrid mode:                  NO\n"
#endif
#ifdef EL_HAVE_QT5
      "  Have Qt5:                     YES\n"
#else
      "  Have Qt5:                     NO\n"
#endif
#ifdef EL_AVOID_COMPLEX_MPI
      "  Avoiding complex MPI:         YES\n"
#else
      "  Avoiding complex MPI:         NO\n"
#endif
#ifdef EL_HAVE_MPI_REDUCE_SCATTER_BLOCK
      "  Have MPI_Reducescatter_block: YES\n"
#else
      "  Have MPI_Reducescatter_block: NO\n"
#endif
#ifdef EL_REDUCE_SCATTER_BLOCK_VIA_ALLREDUCE
      "  AllReduce ReduceScatterBlock: YES\n"
#else
      "  AllReduce ReduceScatterBlock: NO\n"
#endif
#ifdef EL_USE_BYTE_ALLGATHERS
      "  Use byte AllGathers:          YES\n",
#else
      "  Use byte AllGathers:          NO\n",
#endif
      EL_MATH_LIBS );
    return EL_SUCCESS;
}

ElError ElPrintCCompilerInfo( FILE* stream )
{
    // There does not seem to be a portable means of converting C-style
    // filehandles to C++ filestreams, so we will simply reproduce the 
    // functionality of El::PrintCCompilerInfo.    
    fprintf
    ( stream, 
      "Elemental's C compiler info:\n"
      "  EL_CMAKE_C_COMPILER:    %s\n"
      "  EL_MPI_C_COMPILER:      %s\n"
      "  EL_MPI_C_INCLUDE_PATHS: %s\n"
      "  EL_MPI_C_COMPILE_FLAGS: %s\n"
      "  EL_MPI_LINK_FLAGS:      %s\n"
      "  EL_MPI_C_LIBRARIES:     %s\n",
      EL_CMAKE_C_COMPILER, EL_MPI_C_COMPILER, EL_MPI_C_INCLUDE_PATH,
      EL_MPI_C_COMPILE_FLAGS, EL_MPI_LINK_FLAGS, EL_MPI_C_LIBRARIES );
    return EL_SUCCESS;
}

ElError ElPrintCxxCompilerInfo( FILE* stream )
{
    // There does not seem to be a portable means of converting C-style
    // filehandles to C++ filestreams, so we will simply reproduce the 
    // functionality of El::PrintCxxCompilerInfo.    
    fprintf
    ( stream, 
      "Elemental's C compiler info:\n"
      "  EL_CMAKE_CXX_COMPILER:    %s\n"
      "  EL_MPI_CXX_COMPILER:      %s\n"
      "  EL_MPI_CXX_INCLUDE_PATHS: %s\n"
      "  EL_MPI_CXX_COMPILE_FLAGS: %s\n"
      "  EL_MPI_LINK_FLAGS:        %s\n"
      "  EL_MPI_CXX_LIBRARIES:     %s\n",
      EL_CMAKE_CXX_COMPILER, EL_MPI_CXX_COMPILER, EL_MPI_CXX_INCLUDE_PATH,
      EL_MPI_CXX_COMPILE_FLAGS, EL_MPI_LINK_FLAGS, EL_MPI_CXX_LIBRARIES );
    return EL_SUCCESS;
}

ElError ElUsing64BitInt( bool* using64 )
{ EL_TRY( *using64 = El::Using64BitInt() ) }

ElError ElUsing64BitBlasInt( bool* using64 )
{ EL_TRY( *using64 = El::Using64BitBlasInt() ) }

ElError ElSizeOfBool( unsigned* boolSize )
{ EL_TRY( *boolSize = sizeof(bool); ) }

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
