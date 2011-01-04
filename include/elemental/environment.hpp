/*
   Copyright (c) 2009-2011, Jack Poulson
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

// Template conventions:
//   G: general datatype
//
//   T: any ring, e.g., the (Gaussian) integers and the real/complex numbers
//   Z: representation of a real ring, e.g., the integers or real numbers
//   std::complex<Z>: representation of a complex ring, e.g. Gaussian integers
//                    or complex numbers
//
//   F: representation of real or complex number
//   R: representation of real number
//   std::complex<R>: representation of complex number

namespace elemental {

void Init( int* argc, char** argv[] );
void Finalize();

// Naive blocksize set and get
int Blocksize();
void SetBlocksize( int blocksize );

void PushBlocksizeStack( int blocksize );
void PopBlocksizeStack();

template<typename Z>
Z
Abs( Z alpha );

#ifndef WITHOUT_COMPLEX
template<typename Z>
Z
Abs( std::complex<Z> alpha );
#endif

template<typename Z>
Z
FastAbs( Z alpha );

#ifndef WITHOUT_COMPLEX
template<typename Z>
Z
FastAbs( std::complex<Z> alpha );
#endif

template<typename Z>
Z
Conj( Z alpha );

#ifndef WITHOUT_COMPLEX
template<typename Z>
std::complex<Z>
Conj( std::complex<Z> alpha );
#endif

} // elemental

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

template<typename Z>
inline Z 
elemental::Abs
( Z alpha )
{ return std::abs(alpha); }

#ifndef WITHOUT_COMPLEX
template<typename Z>
inline Z
elemental::Abs
( std::complex<Z> alpha )
{ return std::abs( alpha ); }
#endif

template<typename Z>
inline Z
elemental::FastAbs
( Z alpha )
{ return std::abs(alpha); }

#ifndef WITHOUT_COMPLEX
template<typename Z>
inline Z
elemental::FastAbs
( std::complex<Z> alpha )
{ return std::abs( std::real(alpha) ) + std::abs( std::imag(alpha) ); }
#endif

template<typename Z>
inline Z
elemental::Conj
( Z alpha )
{ return alpha; }

#ifndef WITHOUT_COMPLEX
template<typename Z>
inline std::complex<Z>
elemental::Conj
( std::complex<Z> alpha )
{ return std::conj( alpha ); }
#endif

#endif /* ELEMENTAL_ENVIRONMENT_HPP */

