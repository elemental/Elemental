/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_ENVIRONMENT_C_H
#define EL_ENVIRONMENT_C_H

#ifdef __cplusplus
extern "C" {
#else
#include <stdlib.h>
#endif

EL_EXPORT ElError ElPrintVersion( FILE* stream );
EL_EXPORT ElError ElPrintConfig( FILE* stream );
EL_EXPORT ElError ElPrintCCompilerInfo( FILE* stream );
EL_EXPORT ElError ElPrintCxxCompilerInfo( FILE* stream );
EL_EXPORT ElError ElUsing64BitInt( bool* using64 );
EL_EXPORT ElError ElUsing64BitBlasInt( bool* using64 );
EL_EXPORT ElError ElSizeOfBool( unsigned* boolSize );
EL_EXPORT ElError ElSizeOfCBool( unsigned* boolSize );

EL_EXPORT ElError ElInitialize( int* argc, char*** argv );
EL_EXPORT ElError ElFinalize();
EL_EXPORT ElError ElInitialized( bool* initialized );

EL_EXPORT ElError ElInput_b
( const char* name, const char* desc, bool defaultVal, bool* val );
EL_EXPORT ElError ElInput_i
( const char* name, const char* desc, int defaultVal, int* val );
EL_EXPORT ElError ElInput_I
( const char* name, const char* desc, ElInt defaultVal, ElInt* val );
EL_EXPORT ElError ElInput_s
( const char* name, const char* desc, float defaultVal, float* val );
EL_EXPORT ElError ElInput_d
( const char* name, const char* desc, double defaultVal, double* val );
EL_EXPORT ElError ElInput_cstr
( const char* name, const char* desc, const char* defaultVal, 
  const char** val );

EL_EXPORT ElError ElProcessInput();
EL_EXPORT ElError ElPrintInputReport();

EL_EXPORT ElError ElBlocksize( ElInt* blocksize );
EL_EXPORT ElError ElSetBlocksize( ElInt blocksize );

EL_EXPORT ElError ElPushBlocksizeStack( ElInt blocksize );
EL_EXPORT ElError ElPopBlocksizeStack();

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
