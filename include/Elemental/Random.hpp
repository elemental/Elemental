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
#ifndef ELEMENTAL_RANDOM_HPP
#define ELEMENTAL_RANDOM_HPP 1

#ifndef WITHOUT_COMPLEX
#include <complex>
#endif

namespace Elemental {

// Generate a sample from a uniform PDF over the unit ball about the origin 
// of the vector space implied by the type T
template<typename T> T Random();

}

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

namespace Elemental {

template<>
inline int
Random<int>()
{
    int sample = rand();
    if( sample <= RAND_MAX/3 )
        return -1;
    else if( sample <= (RAND_MAX/3)*2 )
        return 0;
    else
        return +1;
}

template<>
inline float
Random<float>()
{ return ( 2*static_cast<float>(rand())/RAND_MAX-1 ); }

template<>
inline double
Random<double>()
{ return ( 2*static_cast<double>(rand())/RAND_MAX-1 ); }

#ifndef WITHOUT_COMPLEX
template<>
inline std::complex<float>
Random< std::complex<float> >()
{
    float r = Random<float>();
    float angle = 3.141592653589793 * Random<float>();

    std::complex<float> u = std::polar(r,angle);

    return u;
}

template<>
inline std::complex<double>
Random< std::complex<double> >()
{
    double r = Random<double>();
    double angle = 3.141592653589793 * Random<double>();

    std::complex<double> u = std::polar(r,angle);

    return u;
}
#endif

} // Elemental

#endif  /* ELEMENTAL_RANDOM_HPP */

