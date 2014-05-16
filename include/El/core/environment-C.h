/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_ENVIRONMENT_C_H
#define EL_ENVIRONMENT_C_H

#ifdef __cplusplus
extern "C" {
#else
#include <stdlib.h>
#endif

ElError ElPrintVersion( FILE* stream );
ElError ElPrintConfig( FILE* stream );
ElError ElPrintCCompilerInfo( FILE* stream );
ElError ElPrintCxxCompilerInfo( FILE* stream );

ElError ElInitialize( int* argc, char*** argv );
ElError ElFinalize();
ElError ElInitialized( bool* initialized );

ElError ElInput_b
( const char* name, const char* desc, bool defaultVal, bool* val );
ElError ElInput_i
( const char* name, const char* desc, int defaultVal, int* val );
ElError ElInput_I
( const char* name, const char* desc, ElInt defaultVal, ElInt* val );
ElError ElInput_s
( const char* name, const char* desc, float defaultVal, float* val );
ElError ElInput_d
( const char* name, const char* desc, double defaultVal, double* val );
ElError ElInput_cstr
( const char* name, const char* desc, const char* defaultVal, 
  const char** val );

ElError ElProcessInput();
ElError ElPrintInputReport();

ElError ElBlocksize( ElInt* blocksize );
ElError ElSetBlocksize( ElInt blocksize );

ElError ElPushBlocksizeStack( ElInt blocksize );
ElError ElPopBlocksizeStack();

#define EL_ABORT_ON_ERROR(error) \
  do \
  { \
      if( error != EL_SUCCESS ) \
      { \
          int commRank_; MPI_Comm_rank( MPI_COMM_WORLD, &commRank_ ); \
          const char* errString = ElErrorString(error); \
          fprintf( stderr, \
            "Rank %d aborting: %s was returned from line %d of file %s\n", \
            commRank_, errString, __LINE__, __FILE__ ); \
          MPI_Abort( MPI_COMM_WORLD, 1 ); \
      } \
  } while( 0 )

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_ENVIRONMENT_C_H */
