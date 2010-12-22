/*
   Copyright (c) 2009-2010, Jack Poulson
   All rights reserved.

   This file is part of Elemental.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

    - Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

    - Neither the name of the owner nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/
#ifndef ELEMENTAL_ENVIRONMENT_HPP
#define ELEMENTAL_ENVIRONMENT_HPP 1

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iostream>
#include <memory>
#include <sstream>
#include <stack>
#include <stdexcept>
#include <vector>

#include "elemental/config.hpp"

// If defined, the _OPENMP macro contains the date of the specification
#ifdef _OPENMP
# include <omp.h>
# if _OPENMP >= 200805
#  define COLLAPSE(N) collapse(N)
# else
#  define COLLAPSE(N) 
# endif 
#endif

#ifndef RELEASE
namespace elemental {

void PushCallStack( std::string s );
void PopCallStack();
void DumpCallStack();

}
#endif

#include "elemental/memory.hpp"
#include "elemental/grid.hpp"
#include "elemental/random.hpp"
#include "elemental/timer.hpp"
#include "elemental/types.hpp"
#include "elemental/utilities.hpp"
#include "elemental/wrappers.hpp"

namespace elemental {

void Init( int* argc, char** argv[] );
void Finalize();

// Naive blocksize set and get
int Blocksize();
void SetBlocksize( int blocksize );

void PushBlocksizeStack( int blocksize );
void PopBlocksizeStack();

template<typename R>
R
Abs( R alpha );

#ifndef WITHOUT_COMPLEX
template<typename R>
R
Abs( std::complex<R> alpha );
#endif

template<typename R>
R
FastAbs( R alpha );

#ifndef WITHOUT_COMPLEX
template<typename R>
R
FastAbs( std::complex<R> alpha );
#endif

template<typename R>
R
Conj( R alpha );

#ifndef WITHOUT_COMPLEX
template<typename R>
std::complex<R>
Conj( std::complex<R> alpha );
#endif

} // elemental

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

template<typename R>
inline R
elemental::Abs
( R alpha )
{ return fabs(alpha); }

#ifndef WITHOUT_COMPLEX
template<typename R>
inline R
elemental::Abs
( std::complex<R> alpha )
{ return std::abs( alpha ); }
#endif

template<typename R>
inline R
elemental::FastAbs
( R alpha )
{ return fabs(alpha); }

#ifndef WITHOUT_COMPLEX
template<typename R>
inline R
elemental::FastAbs
( std::complex<R> alpha )
{ return fabs( std::real(alpha) ) + fabs( std::imag(alpha) ); }
#endif

template<typename R>
inline R
elemental::Conj
( R alpha )
{ return alpha; }

#ifndef WITHOUT_COMPLEX
template<typename R>
inline std::complex<R>
elemental::Conj
( std::complex<R> alpha )
{ return std::conj( alpha ); }
#endif

#endif /* ELEMENTAL_ENVIRONMENT_HPP */

