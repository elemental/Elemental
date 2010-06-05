/*
   This file is part of elemental, a library for distributed-memory dense 
   linear algebra.

   Copyright (C) 2009-2010 Jack Poulson <jack.poulson@gmail.com>

   This program is released under the terms of the license contained in the 
   file LICENSE.
*/
#ifndef ELEMENTAL_TYPES_HPP
#define ELEMENTAL_TYPES_HPP 1

#ifndef WITHOUT_COMPLEX
#include <complex>
#endif
#include <iostream>
#include <limits>
#include <string>

namespace elemental {

#ifndef WITHOUT_COMPLEX
typedef std::complex<float>  scomplex; 
typedef std::complex<double> dcomplex;
#endif

enum Diagonal
{
    NonUnit,
    Unit
};

inline char DiagonalToChar( Diagonal diagonal )
{
    char diagonalChar;
    switch( diagonal )
    {
        case NonUnit: diagonalChar = 'N'; break;
        default:      diagonalChar = 'U'; break;
    }
    return diagonalChar;
}

inline Diagonal CharToDiagonal( char c )
{
    Diagonal diagonal;
    switch( c )
    {
        case 'N': diagonal = NonUnit; break;
        case 'U': diagonal = Unit;    break;
        default:
            throw "CharToDiagonal expects char in {N,U}."; 
    }
    return diagonal;
}

enum Distribution
{
    MC,  // Col of a matrix distribution
    MD,  // Diagonal of a matrix distribution
    MR,  // Row of a matrix distribution
    VC,  // Col-major vector distribution
    VR,  // Row-major vector distribution
    Star // Do not distribute
};

inline std::string DistToString( Distribution distribution )
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

inline Distribution StringToDist( std::string s )
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
        distribution = Star;
    else
    {
        throw "StringToDist expects string in "
              "{\"MC\",\"MD\",\"MR\",\"VC\",\"VR\",\"* \",\" *\",\"*\"}.";
    }
    return distribution;
}

enum Orientation
{
    Normal,
    Transpose,
    ConjugateTranspose
};

inline char OrientationToChar( Orientation orientation )
{
    char orientationChar;
    switch( orientation )
    {
        case Normal:             orientationChar = 'N'; break;
        case Transpose:          orientationChar = 'T'; break;
        default:                 orientationChar = 'C'; break;
    }
    return orientationChar;
}

inline Orientation CharToOrientation( char c )
{
    Orientation orientation;
    switch( c )
    {
        case 'N': orientation = Normal;             break;
        case 'T': orientation = Transpose;          break;
        case 'C': orientation = ConjugateTranspose; break;
        default:
            throw "CharToOrientation expects char in {N,T,C}.";
    }
    return orientation;
}

enum Shape
{
    Lower,
    Upper
};

inline char ShapeToChar( Shape shape )
{
    char shapeChar;
    switch( shape )
    {
        case Lower: shapeChar = 'L'; break;
        default:    shapeChar = 'U'; break;
    }
    return shapeChar;
}

inline Shape CharToShape( char c )
{
    Shape shape;
    switch( c )
    {
        case 'L': shape = Lower; break;
        case 'U': shape = Upper; break;
        default:
            throw "CharToShape expects char in {L,U}.";
    }
    return shape;
}

enum Side
{
    Left,
    Right
};

inline char SideToChar( Side side )
{
    char sideChar;
    switch( side )
    {
        case Left:  sideChar = 'L'; break;
        default:    sideChar = 'R'; break;
    }
    return sideChar;
}
    
inline Side CharToSide( char c )
{
    Side side;
    switch( c )
    {
        case 'L': side = Left;  break;
        case 'R': side = Right; break;
        default:
            throw "CharToSide expects char in {L,R}.";
    }
    return side;
}

} // elemental

#endif /* ELEMENTAL_TYPES_HPP */

