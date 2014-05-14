/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_ENVIRONMENT_CINT_H
#define EL_ENVIRONMENT_CINT_H

#ifdef __cplusplus
extern "C" {
#endif

void ElPrintVersion( FILE* stream );
void ElPrintConfig( FILE* stream );
void ElPrintCCompilerInfo( FILE* stream );
void ElPrintCxxCompilerInfo( FILE* stream );

void ElInitialize( int* argc, char*** argv );
void ElFinalize();
bool ElInitialized();

bool   ElInput_b( const char* name, const char* desc, bool defaultVal );
int    ElInput_i( const char* name, const char* desc, int defaultVal );
ElInt  ElInput_I( const char* name, const char* desc, ElInt defaultVal );
float  ElInput_s( const char* name, const char* desc, float defaultVal );
double ElInput_d( const char* name, const char* desc, double defaultVal );
const char* ElInput_cstr
( const char* name, const char* desc, const char* defaultVal );

void ElProcessInput();
void ElPrintInputReport();

ElInt ElBlocksize();
void ElSetBlocksize( ElInt blocksize );

void ElPushBlocksizeStack( ElInt blocksize );
void ElPopBlocksizeStack();

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_ENVIRONMENT_CINT_H */
