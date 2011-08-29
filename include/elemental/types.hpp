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
#ifndef ELEMENTAL_TYPES_HPP
#define ELEMENTAL_TYPES_HPP 1

#ifndef WITHOUT_COMPLEX
#include <complex>
#endif
#include <string>

namespace elemental {

typedef unsigned char byte;
 
#ifndef WITHOUT_COMPLEX
typedef std::complex<float>  scomplex; 
typedef std::complex<double> dcomplex;
#endif

enum Diagonal
{
    NON_UNIT,
    UNIT
};

char DiagonalToChar( Diagonal diagonal );
Diagonal CharToDiagonal( char c );

enum Distribution
{
    MC,  // Col of a matrix distribution
    MD,  // Diagonal of a matrix distribution
    MR,  // Row of a matrix distribution
    VC,  // Col-major vector distribution
    VR,  // Row-major vector distribution
    STAR // Do not distribute
};

std::string DistToString( Distribution distribution );
Distribution StringToDist( std::string s );

enum Orientation
{
    NORMAL,
    TRANSPOSE,
    ADJOINT
};

char OrientationToChar( Orientation orientation );
Orientation CharToOrientation( char c );

enum Shape
{
    LOWER,
    UPPER
};

char ShapeToChar( Shape shape );
Shape CharToShape( char c );

enum Side
{
    LEFT,
    RIGHT
};

char SideToChar( Side side );
Side CharToSide( char c );

enum VectorDirection
{
    VERTICAL,
    HORIZONTAL
};

enum ForwardOrBackward
{
    FORWARD,
    BACKWARD
};

enum Conjugation
{
    UNCONJUGATED,
    CONJUGATED
};

enum GridOrder
{
    ROW_MAJOR,
    COLUMN_MAJOR
};

//
// We ensure that all enums are lifted into the elemental namespace so that
// they can be conveniently used as function arguments, e.g., 
//
// using namespace elemental;
// Matrix<double> A;
// ...
// double norm = advanced::Norm( A, ONE_NORM );
//

enum NormType
{
    FROBENIUS_NORM, // Compute the "vector" L_2 norm of the matrix
    INFINITY_NORM,  // Compute the L_oo norm of the matrix
    MAX_NORM,       // Compute the "vector" L_oo norm of the matrix
    ONE_NORM        // Compute the L_1 norm of the matrix
};

enum HermitianGenDefiniteEigType
{
    AXBX=1,
    ABX=2,
    BAX=3
};

enum HermitianTridiagApproach
{
    HERMITIAN_TRIDIAG_NORMAL, // Keep the current grid
    HERMITIAN_TRIDIAG_SQUARE, // Drop to a square process grid
    HERMITIAN_TRIDIAG_DEFAULT // Square grid algorithm only if already square
};

} // elemental

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

inline char elemental::DiagonalToChar( Diagonal diagonal )
{
    char diagonalChar;
    switch( diagonal )
    {
        case NON_UNIT: diagonalChar = 'N'; break;
        default:       diagonalChar = 'U'; break;
    }
    return diagonalChar;
}

inline elemental::Diagonal elemental::CharToDiagonal( char c )
{
    Diagonal diagonal;
    switch( c )
    {
        case 'N': diagonal = NON_UNIT; break;
        case 'U': diagonal = UNIT;     break;
        default:
            throw std::logic_error("CharToDiagonal expects char in {N,U}.");
    }
    return diagonal;
}

inline std::string elemental::DistToString
( Distribution distribution )
{
    std::string distString;
    switch( distribution )
    {
        case MC: distString = "MC"; break;
        case MD: distString = "MD"; break;
        case MR: distString = "MR"; break;
        case VC: distString = "VC"; break;
        case VR: distString = "VR"; break;
        default: distString = "* "; break;
    }
    return distString;
}

inline elemental::Distribution elemental::StringToDist( std::string s )
{
    Distribution distribution;
    if( s == "MC" )
        distribution = MC;
    else if( s == "MD" )
        distribution = MD;
    else if( s == "MR" )
        distribution = MR;
    else if( s == "VC" )
        distribution = VC;
    else if( s == "VR" )
        distribution = VR;
    else if( s == "* " || s == " *" || s == "*" )
        distribution = STAR;
    else
    {
        throw std::logic_error("StringToDist expects string in "
              "{\"MC\",\"MD\",\"MR\",\"VC\",\"VR\",\"* \",\" *\",\"*\"}.");
    }
    return distribution;
}

inline char elemental::OrientationToChar( Orientation orientation )
{
    char orientationChar;
    switch( orientation )
    {
        case NORMAL:    orientationChar = 'N'; break;
        case TRANSPOSE: orientationChar = 'T'; break;
        default:        orientationChar = 'C'; break;
    }
    return orientationChar;
}

inline elemental::Orientation elemental::CharToOrientation( char c )
{
    Orientation orientation;
    switch( c )
    {
        case 'N': orientation = NORMAL;    break;
        case 'T': orientation = TRANSPOSE; break;
        case 'C': orientation = ADJOINT;   break;
        default:
            throw std::logic_error
            ("CharToOrientation expects char in {N,T,C}.");
    }
    return orientation;
}

inline char elemental::ShapeToChar( Shape shape )
{
    char shapeChar;
    switch( shape )
    {
        case LOWER: shapeChar = 'L'; break;
        default:    shapeChar = 'U'; break;
    }
    return shapeChar;
}

inline elemental::Shape elemental::CharToShape( char c )
{
    Shape shape;
    switch( c )
    {
        case 'L': shape = LOWER; break;
        case 'U': shape = UPPER; break;
        default:
            throw std::logic_error("CharToShape expects char in {L,U}.");
    }
    return shape;
}

inline char elemental::SideToChar( Side side )
{
    char sideChar;
    switch( side )
    {
        case LEFT:  sideChar = 'L'; break;
        default:    sideChar = 'R'; break;
    }
    return sideChar;
}
    
inline elemental::Side elemental::CharToSide( char c )
{
    Side side;
    switch( c )
    {
        case 'L': side = LEFT;  break;
        case 'R': side = RIGHT; break;
        default:
            throw std::logic_error("CharToSide expects char in {L,R}.");
    }
    return side;
}

#endif /* ELEMENTAL_TYPES_HPP */

