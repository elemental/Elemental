/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename F>
SafeProduct<F>::SafeProduct( Int numEntries )
: rho(1), kappa(0), n(numEntries)
{ }

namespace DistNS {

std::string DistToString( Dist dist )
{
    std::string distString;
    switch( dist )
    {
        case CIRC: distString = "o "; break;
        case MC:   distString = "MC";   break;
        case MD:   distString = "MD";   break;
        case MR:   distString = "MR";   break;
        case VC:   distString = "VC";   break;
        case VR:   distString = "VR";   break;
        default:   distString = "* ";   break;
    }
    return distString;
}

Dist StringToDist( std::string s )
{
    // Most compilers' logic for detecting potentially uninitialized variables
    // is horrendously bad.
    Dist dist=MC;
    if( s == "MC" )
        dist = MC;
    else if( s == "MD" )
        dist = MD;
    else if( s == "MR" )
        dist = MR;
    else if( s == "VC" )
        dist = VC;
    else if( s == "VR" )
        dist = VR;
    else if( s == "* " || s == " *" || s == "*" )
        dist = STAR;
    else if( s == "o " || s == " o" || s == "o" )
        dist = CIRC;
    else
        LogicError
        ("StringToDist expects string in "
         "{\"MC\",\"MD\",\"MR\",\"VC\",\"VR\",\"* \",\" *\",\"*\"}");
    return dist;
}

} // namespace DistNS

namespace LeftOrRightNS {

char LeftOrRightToChar( LeftOrRight side )
{
    char sideChar;
    switch( side )
    {
        case LEFT:  sideChar = 'L'; break;
        default:    sideChar = 'R'; break;
    }
    return sideChar;
}
    
LeftOrRight CharToLeftOrRight( char c )
{
    // Most compilers' logic for detecting potentially uninitialized variables
    // is horrendously bad.
    LeftOrRight side=LEFT;
    switch( c )
    {
        case 'L': side = LEFT;  break;
        case 'R': side = RIGHT; break;
        default:
            LogicError("CharToLeftOrRight expects char in {L,R}");
    }
    return side;
}

} // namespace LeftOrRightNS

namespace OrientationNS {

char OrientationToChar( Orientation orientation )
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

Orientation CharToOrientation( char c )
{
    // Most compilers' logic for detecting potentially uninitialized variables
    // is horrendously bad.
    Orientation orientation=NORMAL;
    switch( c )
    {
        case 'N': orientation = NORMAL;    break;
        case 'T': orientation = TRANSPOSE; break;
        case 'C': orientation = ADJOINT;   break;
        default:
            LogicError
            ("CharToOrientation expects char in {N,T,C}");
    }
    return orientation;
}

} // namespace OrientationNS

namespace UnitOrNonUnitNS {

char UnitOrNonUnitToChar( UnitOrNonUnit diag )
{
    char diagChar;
    switch( diag )
    {
        case NON_UNIT: diagChar = 'N'; break;
        default:       diagChar = 'U'; break;
    }
    return diagChar;
}

UnitOrNonUnit CharToUnitOrNonUnit( char c )
{
    // Most compilers' logic for detecting potentially uninitialized variables
    // is horrendously bad.
    UnitOrNonUnit diag=NON_UNIT;
    switch( c )
    {
        case 'N': diag = NON_UNIT; break;
        case 'U': diag = UNIT;     break;
        default:
            LogicError("CharToUnitOrNonUnit expects char in {N,U}");
    }
    return diag;
}

} // namespace UnitOrNonUnitNS

namespace UpperOrLowerNS {

char UpperOrLowerToChar( UpperOrLower uplo )
{
    char uploChar;
    switch( uplo )
    {
        case LOWER: uploChar = 'L'; break;
        default:    uploChar = 'U'; break;
    }
    return uploChar;
}

UpperOrLower CharToUpperOrLower( char c )
{
    // Most compilers' logic for detecting potentially uninitialized variables
    // is horrendously bad.
    UpperOrLower uplo=LOWER;
    switch( c )
    {
        case 'L': uplo = LOWER; break;
        case 'U': uplo = UPPER; break;
        default:
            LogicError("CharToUpperOrLower expects char in {L,U}");
    }
    return uplo;
}

} // namespace UpperOrLowerNS

template struct SafeProduct<float>;
template struct SafeProduct<double>;
template struct SafeProduct<Complex<float>>;
template struct SafeProduct<Complex<double>>;
#ifdef EL_HAVE_QD
template struct SafeProduct<DoubleDouble>;
template struct SafeProduct<QuadDouble>;
#endif
#ifdef EL_HAVE_QUAD
template struct SafeProduct<Quad>;
template struct SafeProduct<Complex<Quad>>;
#endif
#ifdef EL_HAVE_MPC
template struct SafeProduct<BigFloat>;
#endif

} // namespace El
