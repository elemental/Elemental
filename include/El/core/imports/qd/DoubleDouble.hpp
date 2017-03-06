/*
   Copyright (c) 2009-2017, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_IMPORTS_QD_DOUBLEDOUBLE_HPP
#define EL_IMPORTS_QD_DOUBLEDOUBLE_HPP

#ifdef EL_HAVE_QD
#include <qd/qd_real.h>

namespace El {

// The dd_real and qd_real classes unfortunately do not provide a means of
// assignment directly from an integer, which would break the large amount of
// Elemental (user and library) code which makes use of assignments of the
// form "ABuf[i+j*ALDim] = 0;".

// TODO(poulson): Move constructors and assignments?

struct DoubleDouble : public dd_real
{
    DoubleDouble() { }

    template<typename T,
             typename=EnableIf<IsScalar<T>>>
    DoubleDouble( const T& a );

    DoubleDouble( const unsigned& a ): dd_real(double(a)) { }
    DoubleDouble( const int& a ) : dd_real(a) { }
    DoubleDouble( const float& a ) : dd_real(double(a)) { }
    DoubleDouble( const double& a ) : dd_real(a) { }
    DoubleDouble( const dd_real& a ) : dd_real(a) { }
    DoubleDouble( const char* s ) : dd_real(s) { }

    template<typename T>
    DoubleDouble& operator=( const T& a )
    { dd_real::operator=(DoubleDouble(a)); return *this; }

    DoubleDouble& operator=( const dd_real& a )
    { dd_real::operator=(a); return *this; }
    DoubleDouble& operator=( const float& a )
    { dd_real::operator=(double(a)); return *this; }
    DoubleDouble& operator=( const double& a )
    { dd_real::operator=(a); return *this; }
    DoubleDouble& operator=( const unsigned& a )
    { dd_real::operator=(double(a)); return *this; }
    DoubleDouble& operator=( const int& a )
    { dd_real::operator=(double(a)); return *this; }
    DoubleDouble& operator=( const char* s )
    { dd_real::operator=(s); return *this; }

    template<typename T>
    DoubleDouble& operator+=( const T& a )
    { dd_real::operator+=(DoubleDouble(a)); return *this; }

    DoubleDouble& operator+=( const unsigned& a )
    { dd_real::operator+=(double(a)); return *this; }
    DoubleDouble& operator+=( const int& a )
    { dd_real::operator+=(double(a)); return *this; }
    DoubleDouble& operator+=( const float& a )
    { dd_real::operator+=(a); return *this; }
    DoubleDouble& operator+=( const double& a )
    { dd_real::operator+=(a); return *this; }
    DoubleDouble& operator+=( const dd_real& a )
    { dd_real::operator+=(a); return *this; }

    template<typename T>
    DoubleDouble& operator-=( const T& a )
    { dd_real::operator-=(DoubleDouble(a)); return *this; }

    DoubleDouble& operator-=( const unsigned& a )
    { dd_real::operator-=(double(a)); return *this; }
    DoubleDouble& operator-=( const int& a )
    { dd_real::operator-=(double(a)); return *this; }
    DoubleDouble& operator-=( const float& a )
    { dd_real::operator-=(double(a)); return *this; }
    DoubleDouble& operator-=( const double& a )
    { dd_real::operator-=(a); return *this; }
    DoubleDouble& operator-=( const dd_real& a )
    { dd_real::operator-=(a); return *this; }

    template<typename T>
    DoubleDouble& operator*=( const T& a )
    { dd_real::operator*=(DoubleDouble(a)); return *this; }

    DoubleDouble& operator*=( const unsigned& a )
    { dd_real::operator*=(double(a)); return *this; }
    DoubleDouble& operator*=( const int& a )
    { dd_real::operator*=(double(a)); return *this; }
    DoubleDouble& operator*=( const float& a )
    { dd_real::operator*=(double(a)); return *this; }
    DoubleDouble& operator*=( const double& a )
    { dd_real::operator*=(a); return *this; }
    DoubleDouble& operator*=( const dd_real& a )
    { dd_real::operator*=(a); return *this; }

    template<typename T>
    DoubleDouble& operator/=( const T& a )
    { dd_real::operator/=(DoubleDouble(a)); return *this; }

    DoubleDouble& operator/=( const unsigned& a )
    { dd_real::operator/=(double(a)); return *this; }
    DoubleDouble& operator/=( const int& a )
    { dd_real::operator/=(double(a)); return *this; }
    DoubleDouble& operator/=( const float& a )
    { dd_real::operator/=(double(a)); return *this; }
    DoubleDouble& operator/=( const double& a )
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
    explicit operator long double() const;
#ifdef EL_HAVE_QUAD
    explicit operator Quad() const;
#endif
#ifdef EL_HAVE_MPC
    explicit operator BigFloat() const;
#endif
};

template<typename T,
         typename/*=EnableIf<IsScalar<T>>*/>
DoubleDouble::DoubleDouble( const T& a )
{
    T b = a;
    x[0] = double(b);
    b -= x[0];
    x[1] = double(b);
}

inline bool operator==( const DoubleDouble& a, const DoubleDouble& b )
{ return static_cast<const dd_real&>(a) == static_cast<const dd_real&>(b); }
inline bool operator!=( const DoubleDouble& a, const DoubleDouble& b )
{ return !(a==b); }

inline DoubleDouble operator+( const DoubleDouble& a, const DoubleDouble& b )
{ return static_cast<const dd_real&>(a) + static_cast<const dd_real&>(b); }
inline DoubleDouble operator-( const DoubleDouble& a, const DoubleDouble& b )
{ return static_cast<const dd_real&>(a) - static_cast<const dd_real&>(b); }
inline DoubleDouble operator*( const DoubleDouble& a, const DoubleDouble& b )
{ return static_cast<const dd_real&>(a) * static_cast<const dd_real&>(b); }
inline DoubleDouble operator/( const DoubleDouble& a, const DoubleDouble& b )
{ return static_cast<const dd_real&>(a) / static_cast<const dd_real&>(b); }

inline DoubleDouble operator+
( const DoubleDouble& a, const unsigned& b )
{ return static_cast<const dd_real&>(a) + double(b); }
inline DoubleDouble operator-
( const DoubleDouble& a, const unsigned& b )
{ return static_cast<const dd_real&>(a) - double(b); }
inline DoubleDouble operator*
( const DoubleDouble& a, const unsigned& b )
{ return static_cast<const dd_real&>(a) * double(b); }
inline DoubleDouble operator/
( const DoubleDouble& a, const unsigned& b )
{ return static_cast<const dd_real&>(a) / double(b); }
inline DoubleDouble operator+
( const DoubleDouble& a, const unsigned long& b )
{ return a + DoubleDouble(b); }
inline DoubleDouble operator-
( const DoubleDouble& a, const unsigned long& b )
{ return a - DoubleDouble(b); }
inline DoubleDouble operator*
( const DoubleDouble& a, const unsigned long& b )
{ return a * DoubleDouble(b); }
inline DoubleDouble operator/
( const DoubleDouble& a, const unsigned long& b )
{ return a / DoubleDouble(b); }
inline DoubleDouble operator+
( const DoubleDouble& a, const unsigned long long& b )
{ return a + DoubleDouble(b); }
inline DoubleDouble operator-
( const DoubleDouble& a, const unsigned long long& b )
{ return a - DoubleDouble(b); }
inline DoubleDouble operator*
( const DoubleDouble& a, const unsigned long long& b )
{ return a * DoubleDouble(b); }
inline DoubleDouble operator/
( const DoubleDouble& a, const unsigned long long& b )
{ return a / DoubleDouble(b); }
inline DoubleDouble operator+
( const DoubleDouble& a, const int& b )
{ return static_cast<const dd_real&>(a) + double(b); }
inline DoubleDouble operator-
( const DoubleDouble& a, const int& b )
{ return static_cast<const dd_real&>(a) - double(b); }
inline DoubleDouble operator*
( const DoubleDouble& a, const int& b )
{ return static_cast<const dd_real&>(a) * double(b); }
inline DoubleDouble operator/
( const DoubleDouble& a, const int& b )
{ return static_cast<const dd_real&>(a) / double(b); }
inline DoubleDouble operator+
( const DoubleDouble& a, const long int& b )
{ return a + DoubleDouble(b); }
inline DoubleDouble operator-
( const DoubleDouble& a, const long int& b )
{ return a - DoubleDouble(b); }
inline DoubleDouble operator*
( const DoubleDouble& a, const long int& b )
{ return a * DoubleDouble(b); }
inline DoubleDouble operator/
( const DoubleDouble& a, const long int& b )
{ return a / DoubleDouble(b); }
inline DoubleDouble operator+
( const DoubleDouble& a, const long long int& b )
{ return a + DoubleDouble(b); }
inline DoubleDouble operator-
( const DoubleDouble& a, const long long int& b )
{ return a - DoubleDouble(b); }
inline DoubleDouble operator*
( const DoubleDouble& a, const long long int& b )
{ return a * DoubleDouble(b); }
inline DoubleDouble operator/
( const DoubleDouble& a, const long long int& b )
{ return a / DoubleDouble(b); }
inline DoubleDouble operator+
( const DoubleDouble& a, const float& b )
{ return static_cast<const dd_real&>(a) + double(b); }
inline DoubleDouble operator-
( const DoubleDouble& a, const float& b )
{ return static_cast<const dd_real&>(a) - double(b); }
inline DoubleDouble operator*
( const DoubleDouble& a, const float& b )
{ return static_cast<const dd_real&>(a) * double(b); }
inline DoubleDouble operator/
( const DoubleDouble& a, const float& b )
{ return static_cast<const dd_real&>(a) / double(b); }
inline DoubleDouble operator+
( const DoubleDouble& a, const double& b )
{ return static_cast<const dd_real&>(a) + b; }
inline DoubleDouble operator-
( const DoubleDouble& a, const double& b )
{ return static_cast<const dd_real&>(a) - b; }
inline DoubleDouble operator*
( const DoubleDouble& a, const double& b )
{ return static_cast<const dd_real&>(a) * b; }
inline DoubleDouble operator/
( const DoubleDouble& a, const double& b )
{ return static_cast<const dd_real&>(a) / b; }
#ifdef EL_HAVE_QUAD
inline DoubleDouble operator+
( const DoubleDouble& a, const Quad& b )
{ return a + DoubleDouble(b); }
inline DoubleDouble operator-
( const DoubleDouble& a, const Quad& b )
{ return a - DoubleDouble(b); }
inline DoubleDouble operator*
( const DoubleDouble& a, const Quad& b )
{ return a * DoubleDouble(b); }
inline DoubleDouble operator/
( const DoubleDouble& a, const Quad& b )
{ return a / DoubleDouble(b); }
#endif

inline DoubleDouble operator+
( const unsigned& a, const DoubleDouble& b )
{ return double(a) + static_cast<const dd_real&>(b); }
inline DoubleDouble operator-
( const unsigned& a, const DoubleDouble& b )
{ return double(a) - static_cast<const dd_real&>(b); }
inline DoubleDouble operator*
( const unsigned& a, const DoubleDouble& b )
{ return double(a) * static_cast<const dd_real&>(b); }
inline DoubleDouble operator/
( const unsigned& a, const DoubleDouble& b )
{ return double(a) / static_cast<const dd_real&>(b); }
inline DoubleDouble operator+
( const unsigned long& a, const DoubleDouble& b )
{ return DoubleDouble(a) += b; }
inline DoubleDouble operator-
( const unsigned long& a, const DoubleDouble& b )
{ return DoubleDouble(a) -= b; }
inline DoubleDouble operator*
( const unsigned long& a, const DoubleDouble& b )
{ return DoubleDouble(a) *= b; }
inline DoubleDouble operator/
( const unsigned long& a, const DoubleDouble& b )
{ return DoubleDouble(a) /= b; }
inline DoubleDouble operator+
( const unsigned long long& a, const DoubleDouble& b )
{ return DoubleDouble(a) += b; }
inline DoubleDouble operator-
( const unsigned long long& a, const DoubleDouble& b )
{ return DoubleDouble(a) -= b; }
inline DoubleDouble operator*
( const unsigned long long& a, const DoubleDouble& b )
{ return DoubleDouble(a) *= b; }
inline DoubleDouble operator/
( const unsigned long long& a, const DoubleDouble& b )
{ return DoubleDouble(a) /= b; }
inline DoubleDouble operator+
( const int& a, const DoubleDouble& b )
{ return double(a) + static_cast<const dd_real&>(b); }
inline DoubleDouble operator-
( const int& a, const DoubleDouble& b )
{ return double(a) - static_cast<const dd_real&>(b); }
inline DoubleDouble operator*
( const int& a, const DoubleDouble& b )
{ return double(a) * static_cast<const dd_real&>(b); }
inline DoubleDouble operator/
( const int& a, const DoubleDouble& b )
{ return double(a) / static_cast<const dd_real&>(b); }
inline DoubleDouble operator+
( const long int& a, const DoubleDouble& b )
{ return DoubleDouble(a) += b; }
inline DoubleDouble operator-
( const long int& a, const DoubleDouble& b )
{ return DoubleDouble(a) -= b; }
inline DoubleDouble operator*
( const long int& a, const DoubleDouble& b )
{ return DoubleDouble(a) *= b; }
inline DoubleDouble operator/
( const long int& a, const DoubleDouble& b )
{ return DoubleDouble(a) /= b; }
inline DoubleDouble operator+
( const long long int& a, const DoubleDouble& b )
{ return DoubleDouble(a) += b; }
inline DoubleDouble operator-
( const long long int& a, const DoubleDouble& b )
{ return DoubleDouble(a) -= b; }
inline DoubleDouble operator*
( const long long int& a, const DoubleDouble& b )
{ return DoubleDouble(a) *= b; }
inline DoubleDouble operator/
( const long long int& a, const DoubleDouble& b )
{ return DoubleDouble(a) /= b; }
inline DoubleDouble operator+
( const float& a, const DoubleDouble& b )
{ return double(a) + static_cast<const dd_real&>(b); }
inline DoubleDouble operator-
( const float& a, const DoubleDouble& b )
{ return double(a) - static_cast<const dd_real&>(b); }
inline DoubleDouble operator*
( const float& a, const DoubleDouble& b )
{ return double(a) * static_cast<const dd_real&>(b); }
inline DoubleDouble operator/
( const float& a, const DoubleDouble& b )
{ return double(a) / static_cast<const dd_real&>(b); }
inline DoubleDouble operator+
( const double& a, const DoubleDouble& b )
{ return a + static_cast<const dd_real&>(b); }
inline DoubleDouble operator-
( const double& a, const DoubleDouble& b )
{ return a - static_cast<const dd_real&>(b); }
inline DoubleDouble operator*
( const double& a, const DoubleDouble& b )
{ return a * static_cast<const dd_real&>(b); }
inline DoubleDouble operator/
( const double& a, const DoubleDouble& b )
{ return a / static_cast<const dd_real&>(b); }
#ifdef EL_HAVE_QUAD
inline DoubleDouble operator+
( const Quad& a, const DoubleDouble& b )
{ return DoubleDouble(a) += b; }
inline DoubleDouble operator-
( const Quad& a, const DoubleDouble& b )
{ return DoubleDouble(a) -= b; }
inline DoubleDouble operator*
( const Quad& a, const DoubleDouble& b )
{ return DoubleDouble(a) *= b; }
inline DoubleDouble operator/
( const Quad& a, const DoubleDouble& b )
{ return DoubleDouble(a) /= b; }
#endif

} // namespace El
#endif // ifdef EL_HAVE_QD

#endif // ifndef EL_IMPORTS_QD_DOUBLEDOUBLE_HPP
