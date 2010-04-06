/*
   Copyright 2009-2010 Jack Poulson

   This file is part of Elemental.

   Elemental is free software: you can redistribute it and/or modify it under
   the terms of the GNU Lesser General Public License as published by the
   Free Software Foundation; either version 3 of the License, or 
   (at your option) any later version.

   Elemental is distributed in the hope that it will be useful, but 
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with Elemental. If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef ELEMENTAL_ENVIRONMENT_HPP
#define ELEMENTAL_ENVIRONMENT_HPP 1

#include <cstdlib>
#include <iostream>
#include <memory>
#include <sstream>
#include <stack>
#include <vector>
#include "mpi.h"

#ifndef RELEASE
namespace Elemental
{
    void PushCallStack( std::string s );
    void PopCallStack();
    void DumpCallStack();
}
#endif

#include "Elemental/Memory.hpp"
#include "Elemental/Grid.hpp"
#include "Elemental/Random.hpp"
#include "Elemental/Types.hpp"
#include "Elemental/Utilities.hpp"

namespace Elemental
{
    void Init( int* argc, char** argv[] );
    void Finalize();

    // Naive blocksize set and get
    int Blocksize();
    void SetBlocksize( int blocksize );

    void PushBlocksizeStack( int blocksize );
    void PopBlocksizeStack();

    template<typename R>
    R
    Abs( const R& alpha );

#ifndef WITHOUT_COMPLEX
    template<typename R>
    R
    Abs( const std::complex<R>& alpha );
#endif

    template<typename R>
    R
    FastAbs( const R& alpha );

#ifndef WITHOUT_COMPLEX
    template<typename R>
    R
    FastAbs( const std::complex<R>& alpha );
#endif

    template<typename R>
    R
    Conj( const R& alpha );

#ifndef WITHOUT_COMPLEX
    template<typename R>
    std::complex<R>
    Conj( const std::complex<R>& alpha );
#endif

    template<typename R>
    R
    Imag( const R& alpha );

#ifndef WITHOUT_COMPLEX
    template<typename R>
    R
    Imag( const std::complex<R>& alpha );
#endif

    template<typename R>
    R
    Real( const R& alpha );

#ifndef WITHOUT_COMPLEX
    template<typename R>
    R
    Real( const std::complex<R>& alpha );
#endif
}

/*----------------------------------------------------------------------------*/

template<typename R>
inline R
Elemental::Abs
( const R& alpha )
{ return fabs(alpha); }

#ifndef WITHOUT_COMPLEX
template<typename R>
inline R
Elemental::Abs
( const std::complex<R>& alpha )
{ return std::abs( alpha ); }
#endif

template<typename R>
inline R
Elemental::FastAbs
( const R& alpha )
{ return fabs(alpha); }

#ifndef WITHOUT_COMPLEX
template<typename R>
inline R
Elemental::FastAbs
( const std::complex<R>& alpha )
{ return fabs( std::real(alpha) ) + fabs( std::imag(alpha) ); }
#endif

template<typename R>
inline R
Elemental::Conj
( const R& alpha )
{ return alpha; }

#ifndef WITHOUT_COMPLEX
template<typename R>
inline std::complex<R>
Elemental::Conj
( const std::complex<R>& alpha )
{ return std::conj( alpha ); }
#endif

template<typename R>
inline R
Elemental::Imag
( const R& alpha )
{ return 0; }

#ifndef WITHOUT_COMPLEX
template<typename R>
inline R
Elemental::Imag
( const std::complex<R>& alpha )
{ return std::imag( alpha ); }
#endif

template<typename R>
inline R
Elemental::Real
( const R& alpha )
{ return alpha; }

#ifndef WITHOUT_COMPLEX
template<typename R>
inline R
Elemental::Real
( const std::complex<R>& alpha )
{ return std::real( alpha ); }
#endif

#endif /* ELEMENTAL_ENVIRONMENT_HPP */

