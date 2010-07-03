/*
   This file is part of Elemental, a library for distributed-memory dense 
   linear algebra.

   Copyright (c) 2009-2010 Jack Poulson <jack.poulson@gmail.com>.
   All rights reserved.

   This file is released under the terms of the license contained in the file
   LICENSE-PURE.
*/
#ifndef ELEMENTAL_ENVIRONMENT_HPP
#define ELEMENTAL_ENVIRONMENT_HPP 1

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <memory>
#include <sstream>
#include <stack>
#include <stdexcept>
#include <vector>

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

