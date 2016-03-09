/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_IMPORTS_QD_HPP
#define EL_IMPORTS_QD_HPP

#ifdef EL_HAVE_QD
#include "qd/qd_real.h"

namespace El {

#ifdef EL_HAVE_MPC
class BigFloat;
#endif

// The dd_real and qd_real classes unfortunately do not provide a means of
// assignment directly from an integer, which would break the large amount of
// Elemental (user and library) code which makes use of assignments of the
// form "ABuf[i+j*ALDim] = 0;".

// TODO: Move constructors and assignments?

struct DoubleDouble : public dd_real
{
    DoubleDouble() { }
    DoubleDouble( int a ) : dd_real(a) { }
    // TODO: Perform a more accurate conversion
    DoubleDouble( long long int a ): dd_real(double(a)) { }
    DoubleDouble( double a ) : dd_real(a) { }
    DoubleDouble( const dd_real& a ) : dd_real(a) { }
#ifdef EL_HAVE_QUAD
    DoubleDouble( const Quad& a );
#endif
    DoubleDouble( const char* s ) : dd_real(s) { }

    DoubleDouble& operator=( const dd_real& a )
    { dd_real::operator=(a); return *this; }
    DoubleDouble& operator=( double a )
    { dd_real::operator=(a); return *this; }
    DoubleDouble& operator=( int a )
    { dd_real::operator=(a); return *this; }
    // TODO: Perform a more accurate conversion
    DoubleDouble& operator=( long long int a )
    { dd_real::operator=(double(a)); return *this; }
    DoubleDouble& operator=( const char* s )
    { dd_real::operator=(s); return *this; }

    DoubleDouble& operator+=( double a )
    { dd_real::operator+=(a); return *this; }
    DoubleDouble& operator+=( const dd_real& a )
    { dd_real::operator+=(a); return *this; }

    DoubleDouble& operator-=( double a )
    { dd_real::operator-=(a); return *this; }
    DoubleDouble& operator-=( const dd_real& a )
    { dd_real::operator-=(a); return *this; }

    DoubleDouble& operator*=( double a )
    { dd_real::operator*=(a); return *this; }
    DoubleDouble& operator*=( const dd_real& a )
    { dd_real::operator*=(a); return *this; }

    DoubleDouble& operator/=( double a )
    { dd_real::operator/=(a); return *this; }
    DoubleDouble& operator/=( const dd_real& a )
    { dd_real::operator/=(a); return *this; }

    DoubleDouble operator+() const { return *this; }
    DoubleDouble operator-() const { return dd_real::operator-(); }
    // NOTE: It appears to be a bug in QD that dd_real::operator^ is not const
    DoubleDouble operator^( int n ) { return dd_real::operator^(n); }

    // Casting
    explicit operator unsigned() const { return to_int(*this); }
    explicit operator int() const { return to_int(*this); }
    explicit operator long long int() const;
    explicit operator float() const { return to_double(*this); }
    explicit operator double() const { return to_double(*this); }
#ifdef EL_HAVE_QUAD
    explicit operator Quad() const;
#endif
#ifdef EL_HAVE_MPC
    explicit operator BigFloat() const;
#endif
};

inline DoubleDouble operator+( const DoubleDouble& a, const DoubleDouble& b )
{ return static_cast<const dd_real&>(a) + static_cast<const dd_real&>(b); }
inline DoubleDouble operator+( const DoubleDouble& a, double b )
{ return static_cast<const dd_real&>(a) + b; }
inline DoubleDouble operator+( double a, const DoubleDouble& b )
{ return a + static_cast<const dd_real&>(b); }

inline DoubleDouble operator-( const DoubleDouble& a, const DoubleDouble& b )
{ return static_cast<const dd_real&>(a) - static_cast<const dd_real&>(b); }
inline DoubleDouble operator-( const DoubleDouble& a, double b )
{ return static_cast<const dd_real&>(a) - b; }
inline DoubleDouble operator-( double a, const DoubleDouble& b )
{ return a - static_cast<const dd_real&>(b); }

inline DoubleDouble operator*( const DoubleDouble& a, const DoubleDouble& b )
{ return static_cast<const dd_real&>(a) * static_cast<const dd_real&>(b); }
inline DoubleDouble operator*( const DoubleDouble& a, double b )
{ return static_cast<const dd_real&>(a) * b; }
inline DoubleDouble operator*( double a, const DoubleDouble& b )
{ return a * static_cast<const dd_real&>(b); }

inline DoubleDouble operator/( const DoubleDouble& a, const DoubleDouble& b )
{ return static_cast<const dd_real&>(a) / static_cast<const dd_real&>(b); }
inline DoubleDouble operator/( const DoubleDouble& a, double b )
{ return static_cast<const dd_real&>(a) / b; }
inline DoubleDouble operator/( double a, const DoubleDouble& b )
{ return a / static_cast<const dd_real&>(b); }

struct QuadDouble : public qd_real
{
    QuadDouble() { }
    QuadDouble( int a ) : qd_real(a) { }
    // TODO: Perform a more accurate conversion
    QuadDouble( long long int a ) : qd_real(double(a)) { }
    QuadDouble( double a ) : qd_real(a) { }
    QuadDouble( const dd_real& a ) : qd_real(a) { } 
    QuadDouble( const qd_real& a ) : qd_real(a) { } 
#ifdef EL_HAVE_QUAD
    QuadDouble( const Quad& a );
#endif
    QuadDouble( const char* s ) : qd_real(s) { }

    QuadDouble& operator=( const dd_real& a )
    { qd_real::operator=(a); return *this; }
    QuadDouble& operator=( const qd_real& a )
    { qd_real::operator=(a); return *this; }
    QuadDouble& operator=( double a )
    { qd_real::operator=(a); return *this; }
    QuadDouble& operator=( int a )
    { qd_real::operator=(a); return *this; }
    // TODO: Perform a more accurate conversion
    QuadDouble& operator=( long long int a )
    { qd_real::operator=(double(a)); return *this; }
    QuadDouble& operator=( const char* s )
    { qd_real::operator=(s); return *this; }

    QuadDouble& operator+=( double a )
    { qd_real::operator+=(a); return *this; }
    QuadDouble& operator+=( const dd_real& a )
    { qd_real::operator+=(a); return *this; }
    QuadDouble& operator+=( const qd_real& a )
    { qd_real::operator+=(a); return *this; }

    QuadDouble& operator-=( double a )
    { qd_real::operator-=(a); return *this; }
    QuadDouble& operator-=( const dd_real& a )
    { qd_real::operator-=(a); return *this; }
    QuadDouble& operator-=( const qd_real& a )
    { qd_real::operator-=(a); return *this; }

    QuadDouble& operator*=( double a )
    { qd_real::operator*=(a); return *this; }
    QuadDouble& operator*=( const dd_real& a )
    { qd_real::operator*=(a); return *this; }
    QuadDouble& operator*=( const qd_real& a )
    { qd_real::operator*=(a); return *this; }

    QuadDouble& operator/=( double a )
    { qd_real::operator/=(a); return *this; }
    QuadDouble& operator/=( const dd_real& a )
    { qd_real::operator/=(a); return *this; }
    QuadDouble& operator/=( const qd_real& a )
    { qd_real::operator/=(a); return *this; }

    QuadDouble operator+() const { return *this; }
    QuadDouble operator-() const { return qd_real::operator-(); }
    QuadDouble operator^( int n ) const { return qd_real::operator^(n); }

    // Casting
    explicit operator unsigned() const { return to_int(*this); }
    explicit operator int() const { return to_int(*this); }
    explicit operator float() const { return to_double(*this); }
    explicit operator double() const { return to_double(*this); }
    explicit operator long long int() const;
    explicit operator DoubleDouble() const { return to_dd_real(*this); }
#ifdef EL_HAVE_QUAD
    explicit operator Quad() const;
#endif
#ifdef EL_HAVE_MPC
    explicit operator BigFloat() const;
#endif
};

inline QuadDouble operator+( const QuadDouble& a, const QuadDouble& b )
{ return static_cast<const qd_real&>(a) + static_cast<const qd_real&>(b); }
inline QuadDouble operator+( const QuadDouble& a, const DoubleDouble& b )
{ return static_cast<const qd_real&>(a) + static_cast<const dd_real&>(b); }
inline QuadDouble operator+( const DoubleDouble& a, const QuadDouble& b )
{ return static_cast<const dd_real&>(a) + static_cast<const qd_real&>(b); }
inline QuadDouble operator+( const QuadDouble& a, double b )
{ return static_cast<const qd_real&>(a) + b; }
inline QuadDouble operator+( double a, const QuadDouble& b )
{ return a + static_cast<const qd_real&>(b); }

inline QuadDouble operator-( const QuadDouble& a, const QuadDouble& b )
{ return static_cast<const qd_real&>(a) - static_cast<const qd_real&>(b); }
inline QuadDouble operator-( const QuadDouble& a, const DoubleDouble& b )
{ return static_cast<const qd_real&>(a) - static_cast<const dd_real&>(b); }
inline QuadDouble operator-( const DoubleDouble& a, const QuadDouble& b )
{ return static_cast<const dd_real&>(a) - static_cast<const qd_real&>(b); }
inline QuadDouble operator-( const QuadDouble& a, double b )
{ return static_cast<const qd_real&>(a) - b; }
inline QuadDouble operator-( double a, const QuadDouble& b )
{ return a - static_cast<const qd_real&>(b); }

inline QuadDouble operator*( const QuadDouble& a, const QuadDouble& b )
{ return static_cast<const qd_real&>(a) * static_cast<const qd_real&>(b); }
inline QuadDouble operator*( const QuadDouble& a, const DoubleDouble& b )
{ return static_cast<const qd_real&>(a) * static_cast<const dd_real&>(b); }
inline QuadDouble operator*( const DoubleDouble& a, const QuadDouble& b )
{ return static_cast<const dd_real&>(a) * static_cast<const qd_real&>(b); }
inline QuadDouble operator*( const QuadDouble& a, double b )
{ return static_cast<const qd_real&>(a) * b; }
inline QuadDouble operator*( double a, const QuadDouble& b )
{ return a * static_cast<const qd_real&>(b); }

inline QuadDouble operator/( const QuadDouble& a, const QuadDouble& b )
{ return static_cast<const qd_real&>(a) / static_cast<const qd_real&>(b); }
inline QuadDouble operator/( const QuadDouble& a, const DoubleDouble& b )
{ return static_cast<const qd_real&>(a) / static_cast<const dd_real&>(b); }
inline QuadDouble operator/( const DoubleDouble& a, const QuadDouble& b )
{ return static_cast<const dd_real&>(a) / static_cast<const qd_real&>(b); }
inline QuadDouble operator/( const QuadDouble& a, double b )
{ return static_cast<const qd_real&>(a) / b; }
inline QuadDouble operator/( double a, const QuadDouble& b )
{ return a / static_cast<const qd_real&>(b); }

} // namespace El
#endif // ifdef EL_HAVE_QD

#endif // ifndef EL_IMPORTS_QD_HPP
