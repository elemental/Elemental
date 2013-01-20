/*
   A simple parallel Linear Congruential Generator (LCG)

   Copyright (C) 2011 Jack Poulson, Lexing Ying, and
   The University of Texas at Austin

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef CORE_PLCG_HPP
#define CORE_PLCG_HPP

namespace plcg {

// C++98 does not have 64-bit integer support. We instead store each 16-bits
// as a 32-bit unsigned integer in order to allow for multiplication and 
// addition without overflow.
typedef unsigned UInt32;
struct UInt64 
{ 
    UInt32 d[2]; 
    UInt32& operator[]( int i ) { return d[i]; }
    const UInt32& operator[]( int i ) const { return d[i]; }
};
struct ExpandedUInt64 
{ 
    UInt32 d[4]; 
    UInt32& operator[]( int i ) { return d[i]; }
    const UInt32& operator[]( int i ) const { return d[i]; }
};

UInt32 Lower16Bits( UInt32 a );
UInt32 Upper16Bits( UInt32 a );

ExpandedUInt64 Expand( UInt32 a );
ExpandedUInt64 Expand( UInt64 a );
UInt64 Deflate( ExpandedUInt64 a );

// Carry the upper 16 bits in each of the four 32-bit UInt32's used in the
// expanded representation of a UInt64. The topmost 16 bits are simply zeroed.
void CarryUpper16Bits( ExpandedUInt64& a );

// Add/multiply two expanded UInt64's mod 2^64
ExpandedUInt64 AddWith64BitMod( ExpandedUInt64 a, ExpandedUInt64 b );
ExpandedUInt64 MultiplyWith64BitMod( ExpandedUInt64 a, ExpandedUInt64 b );
ExpandedUInt64 IntegerPowerWith64BitMod( ExpandedUInt64 x, ExpandedUInt64 n );
void Halve( ExpandedUInt64& a );

void SeedSerialLcg( UInt64 globalSeed );
void SeedParallelLcg( UInt32 rank, UInt32 commSize, UInt64 globalSeed );

UInt64 SerialLcg();
UInt64 ParallelLcg();
void ManualLcg( ExpandedUInt64 a, ExpandedUInt64 c, ExpandedUInt64& X );

// For grabbing uniform samples from (0,1]
template<typename R> R SerialUniform();
template<> float SerialUniform<float>();
template<> double SerialUniform<double>();
template<typename R> R ParallelUniform();
template<> float ParallelUniform<float>();
template<> double ParallelUniform<double>();

// For generating Gaussian random variables/vectors
template<typename R>
void SerialBoxMuller( R& X, R& Y );
template<typename R>
void ParallelBoxMuller( R& X, R& Y );
template<typename R>
void SerialGaussianRandomVariable( R& X );
template<typename R>
void ParallelGaussianRandomVariable( R& X );

//----------------------------------------------------------------------------//
// Header implementations                                                     //
//----------------------------------------------------------------------------//

inline UInt32
Lower16Bits( UInt32 a )
{
    return a & 0xFFFF;
}

inline UInt32
Upper16Bits( UInt32 a )
{
    return (a >> 16) & 0xFFFF;
}

inline ExpandedUInt64
Expand( UInt32 a )
{
    ExpandedUInt64 b;
    b[0] = Lower16Bits( a );
    b[1] = Upper16Bits( a );
    b[2] = 0U;
    b[3] = 0U;

    return b;
}

inline ExpandedUInt64
Expand( UInt64 a )
{
    ExpandedUInt64 b;
    b[0] = Lower16Bits( a[0] );
    b[1] = Upper16Bits( a[0] );
    b[2] = Lower16Bits( a[1] );
    b[3] = Upper16Bits( a[1] );

    return b;
}

inline UInt64
Deflate( ExpandedUInt64 a )
{
    UInt64 b;
    b[0] = a[0] + ( a[1] << 16 );
    b[1] = a[2] + ( a[3] << 16 );

    return b;
}

inline void
CarryUpper16Bits( ExpandedUInt64& c )
{
    c[1] += Upper16Bits(c[0]);
    c[0] = Lower16Bits(c[0]);
    c[2] += Upper16Bits(c[1]);
    c[1] = Lower16Bits(c[1]);
    c[3] += Upper16Bits(c[2]);
    c[2] = Lower16Bits(c[2]);
    c[3] = Lower16Bits(c[3]);
}

// Multiply two expanded UInt64's (mod 2^64)
//
// We do so by breaking the 64-bit integers into 16-bit pieces so that
// products can be safely be computed with 32-bit unsigned integers.
//
// a = 2^48 a3 + 2^32 a2 + 2^16 a1 + 2^0 a0,
// b = 2^48 b3 + 2^32 b2 + 2^16 b1 + 2^0 b0,
// where a_j, b_j < 2^16.
//
// Then,
// a b = 
//   2^96 ( a3 b3 ) + 
//   2^80 ( a3 b2 + b3 a2 ) + 
//   2^64 ( a3 b1 + a2 b2 + a1 b3 ) +
//   2^48 ( a3 b0 + a2 b1 + a1 b2 + a0 b3 ) +
//   2^32 ( a2 b0 + a1 b1 + a0 b2 ) + 
//   2^16 ( a1 b0 + a0 b1 ) +
//   2^0  ( a0 b0 )
//
// Since c := a b (mod 2^64), only the last four terms must be computed.
inline ExpandedUInt64 
MultiplyWith64BitMod( ExpandedUInt64 a, ExpandedUInt64 b )
{
    UInt32 temp;
    ExpandedUInt64 c;

    // c = 2^0 a0 b0
    temp = a[0]*b[0];
    c[0] = Lower16Bits( temp );
    c[1] = Upper16Bits( temp );

    // c += 2^16 ( a1 b0 + a0 b1 )
    temp = a[1]*b[0];
    c[1] += Lower16Bits( temp );
    c[2] =  Upper16Bits( temp );
    temp = a[0]*b[1];
    c[1] += Lower16Bits( temp );
    c[2] += Upper16Bits( temp );

    // c += 2^32 ( a2 b0 + a1 b1 + a0 b2 )
    temp = a[2]*b[0];
    c[2] += Lower16Bits( temp );
    c[3] =  Upper16Bits( temp );
    temp = a[1]*b[1];
    c[2] += Lower16Bits( temp );
    c[3] += Upper16Bits( temp );
    temp = a[0]*b[2];
    c[2] += Lower16Bits( temp );
    c[3] += Upper16Bits( temp );

    // c += 2^48 ( a3 b0 + a2 b1 + a1 b2 + a0 b3 )
    temp = a[3]*b[0];
    c[3] += Lower16Bits( temp );
    temp = a[2]*b[1];
    c[3] += Lower16Bits( temp );
    temp = a[1]*b[2];
    c[3] += Lower16Bits( temp );
    temp = a[0]*b[3];
    c[3] += Lower16Bits( temp );

    CarryUpper16Bits( c );
    return c;
}

inline ExpandedUInt64 
AddWith64BitMod( ExpandedUInt64 a, ExpandedUInt64 b )
{
    ExpandedUInt64 c;
    c[0] = a[0] + b[0];
    c[1] = a[1] + b[1];
    c[2] = a[2] + b[2];
    c[3] = a[3] + b[3];

    CarryUpper16Bits( c );
    return c;
}

inline void
ManualLcg( ExpandedUInt64 a, ExpandedUInt64 c, ExpandedUInt64& X )
{
    X = MultiplyWith64BitMod( a, X );
    X = AddWith64BitMod( c, X );
}

// Provide a uniform sample from (0,1]
template<>
inline float SerialUniform<float>()
{
    const UInt64 state = SerialLcg();
    // Use the upper 32-bits of the LCG since they are the most random.
    return (float(state[1])+1.f) / 4294967296.f;
}

// Provide a uniform sample from (0,1]
template<>
inline double SerialUniform<double>()
{
    const UInt64 state = SerialLcg();
    // Use the upper 32-bits of the LCG since they are the most random
    // and we cannot rely on the existence of 64-bit integers in C++.
    return (double(state[1])+1.) / 4294967296.;
}

// Provide a uniform sample from (0,1]
template<>
inline float ParallelUniform<float>()
{
    const UInt64 state = ParallelLcg();
    // Use the upper 32-bits of the LCG since they are the most random.
    return (float(state[1])+1.f) / 4294967296.f;
}

// Provide a uniform sample from (0,1]
template<>
inline double ParallelUniform<double>()
{
    const UInt64 state = ParallelLcg();
    // Use the upper 32-bits of the LCG since they are the most random
    // and we cannot rely on the existence of 64-bit integers in C++.
    return (double(state[1])+1.) / 4294967296.;
}

//
// For generating Gaussian random variables/vectors
//

template<typename R>
inline void
SerialBoxMuller( R& X, R& Y )
{
    const R U = SerialUniform<R>();
    const R V = SerialUniform<R>();
    const R A = sqrt(-2*log(U));
    const R c = cos(2*M_PI*V);
    const R s = sin(2*M_PI*V);
    X = A*c;
    Y = A*s;
}

template<typename R>
inline void
ParallelBoxMuller( R& X, R& Y )
{
    const R U = ParallelUniform<R>();
    const R V = ParallelUniform<R>();
    const R A = sqrt(-2*log(U));
    const R c = cos(2*M_PI*V);
    const R s = sin(2*M_PI*V);
    X = A*c;
    Y = A*s;
}

template<typename R>
inline void
SerialGaussianRandomVariable( R& X )
{
    // Use half of Box-Muller
    const R U = SerialUniform<R>();
    const R V = SerialUniform<R>();
    X = sqrt(-2*log(U)) * cos(2*M_PI*V);
}

template<typename R>
inline void
ParallelGaussianRandomVariable( R& X )
{
    // Use half of Box-Muller
    const R U = ParallelUniform<R>();
    const R V = ParallelUniform<R>();
    X = sqrt(-2*log(U)) * cos(2*M_PI*V);
}

} // namespace plcg

#endif // ifndef CORE_PLCG_HPP
