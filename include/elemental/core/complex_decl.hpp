/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace elem {

// Forward declare this since it is used in Complex
//
// TODO: Figure out how to avoid inlining the Complex friend functions so that
//       this forward declaration is no longer required.
template<typename R>
R Abs( const R& alpha );

// TODO: Think about extending to rings instead of just fields.
template<typename R>
struct Complex 
{
    typedef R BaseType;
    R real, imag;

    Complex();
    Complex( R a );
    Complex( R a, R b );
    Complex( const std::complex<R>& alpha );

    Complex<R>& operator=( const R& alpha );
    Complex<R>& operator+=( const R& alpha );
    Complex<R>& operator-=( const R& alpha );
    Complex<R>& operator*=( const R& alpha );
    Complex<R>& operator/=( const R& alpha );
    Complex<R>& operator=( const Complex<R>& alpha );
    Complex<R>& operator+=( const Complex<R>& alpha );
    Complex<R>& operator-=( const Complex<R>& alpha );
    Complex<R>& operator*=( const Complex<R>& alpha );
    Complex<R>& operator/=( const Complex<R>& alpha );

    // Implement these inline so that we do not have to template them
    // (which forfits implicit conversions, so that we would no longer be 
    // able to type 4*alpha when alpha is a Complex instance)
    //
    // TODO: Figure out how to avoid this...

    friend Complex<R> operator+
    ( const Complex<R>& alpha, const Complex<R>& beta )
    { return Complex<R>(alpha.real+beta.real,alpha.imag+beta.imag); }

    friend Complex<R> operator+
    ( const Complex<R>& alpha, const R& beta )
    { return Complex<R>(alpha.real+beta,alpha.imag); }

    friend Complex<R> operator+
    ( const R& alpha, const Complex<R>& beta )
    { return Complex<R>(alpha+beta.real,beta.imag); }

    friend Complex<R> operator-
    ( const Complex<R>& alpha, const Complex<R>& beta )
    { return Complex<R>(alpha.real-beta.real,alpha.imag-beta.imag); }

    friend Complex<R> operator-
    ( const Complex<R>& alpha, const R& beta )
    { return Complex<R>(alpha.real-beta,alpha.imag); }

    friend Complex<R> operator-
    ( const R& alpha, const Complex<R>& beta )
    { return Complex<R>(alpha-beta.real,-beta.imag); }

    friend Complex<R> operator*
    ( const Complex<R>& alpha, const Complex<R>& beta )
    {
        const R a=alpha.real, b=alpha.imag, c=beta.real, d=beta.imag;
        return Complex<R>(a*c-b*d,a*d+b*c);
    }

    friend Complex<R> operator*
    ( const Complex<R>& alpha, const R& beta )
    { return Complex<R>(alpha.real*beta,alpha.imag*beta); }

    friend Complex<R> operator*
    ( const R& alpha, const Complex<R>& beta )
    { return Complex<R>(alpha*beta.real,alpha*beta.imag); }

    friend Complex<R> operator/
    ( const Complex<R>& alpha, const Complex<R>& beta )
    {
        const R a=alpha.real, b=alpha.imag, c=beta.real, d=beta.imag;
        if( Abs(c) >= Abs(d) )
        {
            const R ratio = d/c;
            const R denom = c + d*ratio;
            const R u = (a+b*ratio)/denom;
            const R v = (b-a*ratio)/denom;
            return Complex<R>(u,v);
        }
        else
        {
            const R ratio = c/d;
            const R denom = c*ratio + d;
            const R u = (a*ratio+b)/denom;
            const R v = (b*ratio-a)/denom;
            return Complex<R>(u,v);
        }
    }

    friend Complex<R> operator/
    ( const Complex<R>& alpha, const R& beta )
    { return Complex<R>(alpha.real/beta,alpha.imag/beta); }

    friend Complex<R> operator/
    ( const R& alpha, const Complex<R>& beta )
    {
        const R c=beta.real, d=beta.imag;
        if( Abs(c) >= Abs(d) )
        {
            const R ratio = d/c;
            const R denom = c + d*ratio;
            const R u = alpha/denom;
            const R v = -alpha*ratio/denom;
            return Complex<R>(u,v);
        }
        else
        {
            const R ratio = c/d;
            const R denom = c*ratio + d;
            const R u = alpha*ratio/denom;
            const R v = -alpha/denom;
            return Complex<R>(u,v);
        }
    }

    friend Complex<R> operator+( const Complex<R>& alpha )
    { return alpha; }

    friend Complex<R> operator-( const Complex<R>& alpha )
    { return Complex<R>(-alpha.real,-alpha.imag); }

    friend bool operator==( const Complex<R>& alpha, const Complex<R>& beta )
    { return alpha.real==beta.real && alpha.imag==beta.imag; }

    friend bool operator==( const Complex<R>& alpha, const R& beta )
    { return alpha.real==beta && alpha.imag==0; }

    friend bool operator==( const R& alpha, const Complex<R>& beta )
    { return alpha==beta.real && 0==beta.imag; }

    friend bool operator!=( const Complex<R>& alpha, const Complex<R>& beta )
    { return alpha.real!=beta.real || alpha.imag!=beta.imag; }

    friend bool operator!=( const Complex<R>& alpha, const R& beta )
    { return alpha.real!=beta || alpha.imag!=0; }

    friend bool operator!=( const R& alpha, const Complex<R>& beta )
    { return alpha!=beta.real || 0!=beta.imag; }

    friend std::ostream& operator<<
    ( std::ostream& os, Complex<R> alpha )
    {
        os << alpha.real << "+" << alpha.imag << "i";
        return os;
    }
};

// For extracting the underlying real datatype, 
// e.g., typename Base<Scalar>::type a = 3.0;
template<typename R>
struct Base { typedef R type; };
template<typename R>
struct Base<Complex<R> > { typedef R type; };

// For querying whether or not a scalar is complex,
// e.g., IsComplex<Scalar>::val
template<typename R>
struct IsComplex { enum { val=0 }; };
template<typename R>
struct IsComplex<Complex<R> > { enum { val=1 }; };

} // namespace elem
