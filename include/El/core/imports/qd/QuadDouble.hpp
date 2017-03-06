/*
   Copyright (c) 2009-2017, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_IMPORTS_QD_QUADDOUBLE_HPP
#define EL_IMPORTS_QD_QUADDOUBLE_HPP

#ifdef EL_HAVE_QD
#include <qd/qd_real.h>
#include <El/core/imports/qd/DoubleDouble.hpp>

namespace El {

// TODO(poulson): Move constructors and assignments?

struct QuadDouble : public qd_real
{
    QuadDouble() { }

    template<typename T,
             typename=EnableIf<IsScalar<T>>>
    QuadDouble( const T& a );

    QuadDouble( const unsigned& a ) : qd_real(double(a)) { }
    QuadDouble( const int& a ) : qd_real(a) { }
    QuadDouble( const float& a ) : qd_real(double(a)) { }
    QuadDouble( const double& a ) : qd_real(a) { }
    QuadDouble( const dd_real& a ) : qd_real(a) { }
    QuadDouble( const qd_real& a ) : qd_real(a) { }
    QuadDouble( const char* s ) : qd_real(s) { }

    template<typename T>
    QuadDouble& operator=( const T& a )
    { qd_real::operator=(QuadDouble(a)); return *this; }

    QuadDouble& operator=( const dd_real& a )
    { qd_real::operator=(a); return *this; }
    QuadDouble& operator=( const qd_real& a )
    { qd_real::operator=(a); return *this; }
    QuadDouble& operator=( const double& a )
    { qd_real::operator=(a); return *this; }
    QuadDouble& operator=( const float& a )
    { qd_real::operator=(double(a)); return *this; }
    QuadDouble& operator=( const unsigned& a )
    { qd_real::operator=(double(a)); return *this; }
    QuadDouble& operator=( const int& a )
    { qd_real::operator=(a); return *this; }
    QuadDouble& operator=( const char* s )
    { qd_real::operator=(s); return *this; }

    template<typename T>
    QuadDouble& operator+=( const T& a )
    { qd_real::operator+=(QuadDouble(a)); return *this; }

    QuadDouble& operator+=( const unsigned& a )
    { qd_real::operator+=(double(a)); return *this; }
    QuadDouble& operator+=( const int& a )
    { qd_real::operator+=(double(a)); return *this; }
    QuadDouble& operator+=( const float& a )
    { qd_real::operator+=(double(a)); return *this; }
    QuadDouble& operator+=( const double& a )
    { qd_real::operator+=(a); return *this; }
    QuadDouble& operator+=( const dd_real& a )
    { qd_real::operator+=(a); return *this; }
    QuadDouble& operator+=( const qd_real& a )
    { qd_real::operator+=(a); return *this; }

    template<typename T>
    QuadDouble& operator-=( const T& a )
    { qd_real::operator-=(QuadDouble(a)); return *this; }

    QuadDouble& operator-=( const unsigned& a )
    { qd_real::operator-=(double(a)); return *this; }
    QuadDouble& operator-=( const int& a )
    { qd_real::operator-=(double(a)); return *this; }
    QuadDouble& operator-=( const float& a )
    { qd_real::operator-=(double(a)); return *this; }
    QuadDouble& operator-=( const double& a )
    { qd_real::operator-=(a); return *this; }
    QuadDouble& operator-=( const dd_real& a )
    { qd_real::operator-=(a); return *this; }
    QuadDouble& operator-=( const qd_real& a )
    { qd_real::operator-=(a); return *this; }

    template<typename T>
    QuadDouble& operator*=( const T& a )
    { qd_real::operator*=(QuadDouble(a)); return *this; }

    QuadDouble& operator*=( const unsigned& a )
    { qd_real::operator*=(double(a)); return *this; }
    QuadDouble& operator*=( const int& a )
    { qd_real::operator*=(double(a)); return *this; }
    QuadDouble& operator*=( const float& a )
    { qd_real::operator*=(double(a)); return *this; }
    QuadDouble& operator*=( const double& a )
    { qd_real::operator*=(a); return *this; }
    QuadDouble& operator*=( const dd_real& a )
    { qd_real::operator*=(a); return *this; }
    QuadDouble& operator*=( const qd_real& a )
    { qd_real::operator*=(a); return *this; }

    template<typename T>
    QuadDouble& operator/=( const T& a )
    { qd_real::operator/=(QuadDouble(a)); return *this; }

    QuadDouble& operator/=( const unsigned& a )
    { qd_real::operator/=(double(a)); return *this; }
    QuadDouble& operator/=( const int& a )
    { qd_real::operator/=(double(a)); return *this; }
    QuadDouble& operator/=( const float& a )
    { qd_real::operator/=(double(a)); return *this; }
    QuadDouble& operator/=( const double& a )
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
    explicit operator long double() const;
    explicit operator long long int() const;
    explicit operator DoubleDouble() const { return to_dd_real(*this); }
#ifdef EL_HAVE_QUAD
    explicit operator Quad() const;
#endif
#ifdef EL_HAVE_MPC
    explicit operator BigFloat() const;
#endif
};

template<typename T,
         typename/*=EnableIf<IsScalar<T>>*/>
QuadDouble::QuadDouble( const T& a )
{
    T b = a;
    x[0] = double(b);
    for( Int j=1; j<4; ++j )
    {
        b -= x[j-1];
        x[j] = double(b);
    }
}

inline bool operator==( const QuadDouble& a, const QuadDouble& b )
{ return static_cast<const qd_real&>(a) == static_cast<const qd_real&>(b); }
inline bool operator!=( const QuadDouble& a, const QuadDouble& b )
{ return !(a==b); }

inline QuadDouble operator+
( const QuadDouble& a, const QuadDouble& b )
{ return static_cast<const qd_real&>(a) + static_cast<const qd_real&>(b); }
inline QuadDouble operator-
( const QuadDouble& a, const QuadDouble& b )
{ return static_cast<const qd_real&>(a) - static_cast<const qd_real&>(b); }
inline QuadDouble operator*
( const QuadDouble& a, const QuadDouble& b )
{ return static_cast<const qd_real&>(a) * static_cast<const qd_real&>(b); }
inline QuadDouble operator/
( const QuadDouble& a, const QuadDouble& b )
{ return static_cast<const qd_real&>(a) / static_cast<const qd_real&>(b); }

inline QuadDouble operator+
( const QuadDouble& a, const DoubleDouble& b )
{ return static_cast<const qd_real&>(a) + static_cast<const dd_real&>(b); }
inline QuadDouble operator-
( const QuadDouble& a, const DoubleDouble& b )
{ return static_cast<const qd_real&>(a) - static_cast<const dd_real&>(b); }
inline QuadDouble operator*
( const QuadDouble& a, const DoubleDouble& b )
{ return static_cast<const qd_real&>(a) * static_cast<const dd_real&>(b); }
inline QuadDouble operator/
( const QuadDouble& a, const DoubleDouble& b )
{ return static_cast<const qd_real&>(a) / static_cast<const dd_real&>(b); }
#ifdef EL_HAVE_QUAD
inline QuadDouble operator+
( const QuadDouble& a, const Quad& b )
{ return a + QuadDouble(b); }
inline QuadDouble operator-
( const QuadDouble& a, const Quad& b )
{ return a - QuadDouble(b); }
inline QuadDouble operator*
( const QuadDouble& a, const Quad& b )
{ return a * QuadDouble(b); }
inline QuadDouble operator/
( const QuadDouble& a, const Quad& b )
{ return a / QuadDouble(b); }
#endif
inline QuadDouble operator+
( const QuadDouble& a, const double& b )
{ return static_cast<const qd_real&>(a) + b; }
inline QuadDouble operator-
( const QuadDouble& a, const double& b )
{ return static_cast<const qd_real&>(a) - b; }
inline QuadDouble operator*
( const QuadDouble& a, const double& b )
{ return static_cast<const qd_real&>(a) * b; }
inline QuadDouble operator/
( const QuadDouble& a, const double& b )
{ return static_cast<const qd_real&>(a) / b; }
inline QuadDouble operator+
( const QuadDouble& a, const float& b )
{ return a + double(b); }
inline QuadDouble operator-
( const QuadDouble& a, const float& b )
{ return a - double(b); }
inline QuadDouble operator*
( const QuadDouble& a, const float& b )
{ return a * double(b); }
inline QuadDouble operator/
( const QuadDouble& a, const float& b )
{ return a / double(b); }
inline QuadDouble operator+
( const QuadDouble& a, const int& b )
{ return a + double(b); }
inline QuadDouble operator-
( const QuadDouble& a, const int& b )
{ return a - double(b); }
inline QuadDouble operator*
( const QuadDouble& a, const int& b )
{ return a * double(b); }
inline QuadDouble operator/
( const QuadDouble& a, const int& b )
{ return a / double(b); }
inline QuadDouble operator+
( const QuadDouble& a, const long int& b )
{ return a + DoubleDouble(b); }
inline QuadDouble operator-
( const QuadDouble& a, const long int& b )
{ return a - DoubleDouble(b); }
inline QuadDouble operator*
( const QuadDouble& a, const long int& b )
{ return a * DoubleDouble(b); }
inline QuadDouble operator/
( const QuadDouble& a, const long int& b )
{ return a / DoubleDouble(b); }
inline QuadDouble operator+
( const QuadDouble& a, const long long int& b )
{ return a + DoubleDouble(b); }
inline QuadDouble operator-
( const QuadDouble& a, const long long int& b )
{ return a - DoubleDouble(b); }
inline QuadDouble operator*
( const QuadDouble& a, const long long int& b )
{ return a * DoubleDouble(b); }
inline QuadDouble operator/
( const QuadDouble& a, const long long int& b )
{ return a / DoubleDouble(b); }
inline QuadDouble operator+
( const QuadDouble& a, const unsigned& b )
{ return a + double(b); }
inline QuadDouble operator-
( const QuadDouble& a, const unsigned& b )
{ return a - double(b); }
inline QuadDouble operator*
( const QuadDouble& a, const unsigned& b )
{ return a * double(b); }
inline QuadDouble operator/
( const QuadDouble& a, const unsigned& b )
{ return a / double(b); }
inline QuadDouble operator+
( const QuadDouble& a, const unsigned long& b )
{ return a + DoubleDouble(b); }
inline QuadDouble operator-
( const QuadDouble& a, const unsigned long& b )
{ return a - DoubleDouble(b); }
inline QuadDouble operator*
( const QuadDouble& a, const unsigned long& b )
{ return a * DoubleDouble(b); }
inline QuadDouble operator/
( const QuadDouble& a, const unsigned long& b )
{ return a / DoubleDouble(b); }
inline QuadDouble operator+
( const QuadDouble& a, const unsigned long long& b )
{ return a + DoubleDouble(b); }
inline QuadDouble operator-
( const QuadDouble& a, const unsigned long long& b )
{ return a - DoubleDouble(b); }
inline QuadDouble operator*
( const QuadDouble& a, const unsigned long long& b )
{ return a * DoubleDouble(b); }
inline QuadDouble operator/
( const QuadDouble& a, const unsigned long long& b )
{ return a / DoubleDouble(b); }

inline QuadDouble operator+
( const DoubleDouble& a, const QuadDouble& b )
{ return static_cast<const dd_real&>(a) + static_cast<const qd_real&>(b); }
inline QuadDouble operator-
( const DoubleDouble& a, const QuadDouble& b )
{ return static_cast<const dd_real&>(a) - static_cast<const qd_real&>(b); }
inline QuadDouble operator*
( const DoubleDouble& a, const QuadDouble& b )
{ return static_cast<const dd_real&>(a) * static_cast<const qd_real&>(b); }
inline QuadDouble operator/
( const DoubleDouble& a, const QuadDouble& b )
{ return static_cast<const dd_real&>(a) / static_cast<const qd_real&>(b); }
#ifdef EL_HAVE_QUAD
inline QuadDouble operator+
( const Quad& a, const QuadDouble& b )
{ return QuadDouble(a) += b; }
inline QuadDouble operator-
( const Quad& a, const QuadDouble& b )
{ return QuadDouble(a) -= b; }
inline QuadDouble operator*
( const Quad& a, const QuadDouble& b )
{ return QuadDouble(a) *= b; }
inline QuadDouble operator/
( const Quad& a, const QuadDouble& b )
{ return QuadDouble(a) /= b; }
#endif
inline QuadDouble operator+
( const double& a, const QuadDouble& b )
{ return a + static_cast<const qd_real&>(b); }
inline QuadDouble operator-
( const double& a, const QuadDouble& b )
{ return a - static_cast<const qd_real&>(b); }
inline QuadDouble operator*
( const double& a, const QuadDouble& b )
{ return a * static_cast<const qd_real&>(b); }
inline QuadDouble operator/
( const double& a, const QuadDouble& b )
{ return a / static_cast<const qd_real&>(b); }
inline QuadDouble operator+
( const float& a, const QuadDouble& b )
{ return double(a) + static_cast<const qd_real&>(b); }
inline QuadDouble operator-
( const float& a, const QuadDouble& b )
{ return double(a) - static_cast<const qd_real&>(b); }
inline QuadDouble operator*
( const float& a, const QuadDouble& b )
{ return double(a) * static_cast<const qd_real&>(b); }
inline QuadDouble operator/
( const float& a, const QuadDouble& b )
{ return double(a) / static_cast<const qd_real&>(b); }
inline QuadDouble operator+
( const long long int& a, const QuadDouble& b )
{ return QuadDouble(a) += b; }
inline QuadDouble operator-
( const long long int& a, const QuadDouble& b )
{ return QuadDouble(a) -= b; }
inline QuadDouble operator*
( const long long int& a, const QuadDouble& b )
{ return QuadDouble(a) *= b; }
inline QuadDouble operator/
( const long long int& a, const QuadDouble& b )
{ return QuadDouble(a) /= b; }
inline QuadDouble operator+
( const long int& a, const QuadDouble& b )
{ return QuadDouble(a) += b; }
inline QuadDouble operator-
( const long int& a, const QuadDouble& b )
{ return QuadDouble(a) -= b; }
inline QuadDouble operator*
( const long int& a, const QuadDouble& b )
{ return QuadDouble(a) *= b; }
inline QuadDouble operator/
( const long int& a, const QuadDouble& b )
{ return QuadDouble(a) /= b; }
inline QuadDouble operator+
( const int& a, const QuadDouble& b )
{ return double(a) + static_cast<const qd_real&>(b); }
inline QuadDouble operator-
( const int& a, const QuadDouble& b )
{ return double(a) - static_cast<const qd_real&>(b); }
inline QuadDouble operator*
( const int& a, const QuadDouble& b )
{ return double(a) * static_cast<const qd_real&>(b); }
inline QuadDouble operator/
( const int& a, const QuadDouble& b )
{ return double(a) / static_cast<const qd_real&>(b); }
inline QuadDouble operator+
( const unsigned long long& a, const QuadDouble& b )
{ return QuadDouble(a) += b; }
inline QuadDouble operator-
( const unsigned long long& a, const QuadDouble& b )
{ return QuadDouble(a) -= b; }
inline QuadDouble operator*
( const unsigned long long& a, const QuadDouble& b )
{ return QuadDouble(a) *= b; }
inline QuadDouble operator/
( const unsigned long long& a, const QuadDouble& b )
{ return QuadDouble(a) /= b; }
inline QuadDouble operator+
( const unsigned long& a, const QuadDouble& b )
{ return QuadDouble(a) += b; }
inline QuadDouble operator-
( const unsigned long& a, const QuadDouble& b )
{ return QuadDouble(a) -= b; }
inline QuadDouble operator*
( const unsigned long& a, const QuadDouble& b )
{ return QuadDouble(a) *= b; }
inline QuadDouble operator/
( const unsigned long& a, const QuadDouble& b )
{ return QuadDouble(a) /= b; }
inline QuadDouble operator+
( const unsigned& a, const QuadDouble& b )
{ return double(a) + static_cast<const qd_real&>(b); }
inline QuadDouble operator-
( const unsigned& a, const QuadDouble& b )
{ return double(a) - static_cast<const qd_real&>(b); }
inline QuadDouble operator*
( const unsigned& a, const QuadDouble& b )
{ return double(a) * static_cast<const qd_real&>(b); }
inline QuadDouble operator/
( const unsigned& a, const QuadDouble& b )
{ return double(a) / static_cast<const qd_real&>(b); }

} // namespace El
#endif // ifdef EL_HAVE_QD

#endif // ifndef EL_IMPORTS_QD_QUADDOUBLE_HPP
