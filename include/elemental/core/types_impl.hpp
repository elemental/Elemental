/*
   Copyright (c) 2009-2012, Jack Poulson
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

namespace elem {

template<typename F,typename Int>
inline
SafeProduct<F,Int>::SafeProduct( Int numEntries )
: rho(1), kappa(0), n(numEntries)
{ }

namespace distribution_wrapper {

inline std::string 
DistToString( Distribution distribution )
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

inline Distribution 
StringToDist( std::string s )
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
        throw std::logic_error
        ("StringToDist expects string in "
         "{\"MC\",\"MD\",\"MR\",\"VC\",\"VR\",\"* \",\" *\",\"*\"}");
    }
    return distribution;
}

} // namespace distribution_wrapper

namespace left_or_right_wrapper {

inline char 
LeftOrRightToChar( LeftOrRight side )
{
    char sideChar;
    switch( side )
    {
        case LEFT:  sideChar = 'L'; break;
        default:    sideChar = 'R'; break;
    }
    return sideChar;
}
    
inline LeftOrRight 
CharToLeftOrRight( char c )
{
    LeftOrRight side;
    switch( c )
    {
        case 'L': side = LEFT;  break;
        case 'R': side = RIGHT; break;
        default:
            throw std::logic_error("CharToLeftOrRight expects char in {L,R}");
    }
    return side;
}

} // namespace left_or_right_wrapper

namespace orientation_wrapper {

inline char 
OrientationToChar( Orientation orientation )
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

inline Orientation 
CharToOrientation( char c )
{
    Orientation orientation;
    switch( c )
    {
        case 'N': orientation = NORMAL;    break;
        case 'T': orientation = TRANSPOSE; break;
        case 'C': orientation = ADJOINT;   break;
        default:
            throw std::logic_error
            ("CharToOrientation expects char in {N,T,C}");
    }
    return orientation;
}

} // namespace orientation_wrapper

namespace unit_or_non_unit_wrapper {

inline char 
UnitOrNonUnitToChar( UnitOrNonUnit diag )
{
    char diagChar;
    switch( diag )
    {
        case NON_UNIT: diagChar = 'N'; break;
        default:       diagChar = 'U'; break;
    }
    return diagChar;
}

inline UnitOrNonUnit 
CharToUnitOrNonUnit( char c )
{
    UnitOrNonUnit diag;
    switch( c )
    {
        case 'N': diag = NON_UNIT; break;
        case 'U': diag = UNIT;     break;
        default:
            throw std::logic_error("CharToUnitOrNonUnit expects char in {N,U}");
    }
    return diag;
}

} // namespace unit_or_non_unit_wrapper

namespace upper_or_lower_wrapper {

inline char 
UpperOrLowerToChar( UpperOrLower uplo )
{
    char uploChar;
    switch( uplo )
    {
        case LOWER: uploChar = 'L'; break;
        default:    uploChar = 'U'; break;
    }
    return uploChar;
}

inline UpperOrLower 
CharToUpperOrLower( char c )
{
    UpperOrLower uplo;
    switch( c )
    {
        case 'L': uplo = LOWER; break;
        case 'U': uplo = UPPER; break;
        default:
            throw std::logic_error("CharToUpperOrLower expects char in {L,U}");
    }
    return uplo;
}

} // namespace upper_or_lower_wrapper

} // namespace elem
