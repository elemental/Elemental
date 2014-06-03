/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

namespace El {

template<typename T,typename S>
void Scale( S alphaS, Matrix<T>& A )
{
    DEBUG_ONLY(CallStackEntry cse("Scale"))
    const T alpha = T(alphaS);
    if( alpha != T(1) )
    {
        if( alpha == T(0) )
            for( Int j=0; j<A.Width(); ++j )
                for( Int i=0; i<A.Height(); ++i )
                    A.Set(i,j,0);
        else
            for( Int j=0; j<A.Width(); ++j )
                blas::Scal( A.Height(), alpha, A.Buffer(0,j), 1 );
    }
}

template<typename Real,typename S>
void Scale( S alphaS, Matrix<Real>& AReal, Matrix<Real>& AImag )
{
    DEBUG_ONLY(CallStackEntry cse("Scale"))
    typedef Complex<Real> C;
    const Int m = AReal.Height();
    const Int n = AReal.Width();
    const C alpha = C(alphaS);
    if( alpha != C(1) )
    {
        if( alpha == C(0) )
        {
            Scale( Real(0), AReal );
            Scale( Real(0), AImag );
        }
        else
        {
            Matrix<Real> aReal, aImag, aRealCopy, aImagCopy;
            for( Int j=0; j<n; ++j )
            {
                aReal = View( aReal, 0, j, m, 1 );
                aImag = View( aImag, 0, j, m, 1 );
                aRealCopy = aReal;
                aImagCopy = aImag;
                Scale( alpha.real(), aReal     );
                Axpy( -alpha.imag(), aImagCopy, aReal );
                Scale( alpha.real(), aImag     );
                Axpy(  alpha.imag(), aRealCopy, aImag );
            }
        }
    }
}

template<typename T,typename S>
void Scale( S alpha, AbstractDistMatrix<T>& A )
{ Scale( alpha, A.Matrix() ); }

template<typename Real,typename S>
void Scale( S alpha, AbstractDistMatrix<Real>& AReal, 
                     AbstractDistMatrix<Real>& AImag )
{ Scale( alpha, AReal.Matrix(), AImag.Matrix() ); }

template<typename T,typename S>
void Scale( S alpha, AbstractBlockDistMatrix<T>& A )
{ Scale( alpha, A.Matrix() ); }

template<typename Real,typename S>
void Scale( S alpha, AbstractBlockDistMatrix<Real>& AReal, 
                     AbstractBlockDistMatrix<Real>& AImag )
{ Scale( alpha, AReal.Matrix(), AImag.Matrix() ); }

#define PROTO_TYPES(T,S) \
  template void Scale( S alpha, Matrix<T>& A ); \
  template void Scale( S alpha, AbstractDistMatrix<T>& A ); \
  template void Scale( S alpha, AbstractBlockDistMatrix<T>& A );

#define PROTO_INT(T) PROTO_TYPES(T,T)

#define PROTO_REAL(T) \
  PROTO_TYPES(T,Int) \
  PROTO_TYPES(T,T) 

#define PROTO_CPX(T) \
  PROTO_TYPES(T,Int) \
  PROTO_TYPES(T,Base<T>) \
  PROTO_TYPES(T,T) \
  template void Scale \
  ( Int alpha, Matrix<Base<T>>& AReal, Matrix<Base<T>>& AImag ); \
  template void Scale \
  ( Base<T> alpha, Matrix<Base<T>>& AReal, Matrix<Base<T>>& AImag ); \
  template void Scale \
  ( T alpha, Matrix<Base<T>>& AReal, Matrix<Base<T>>& AImag ); \
  template void Scale \
  ( Int alpha, AbstractDistMatrix<Base<T>>& AReal, \
               AbstractDistMatrix<Base<T>>& AImag ); \
  template void Scale \
  ( Base<T> alpha, AbstractDistMatrix<Base<T>>& AReal, \
                   AbstractDistMatrix<Base<T>>& AImag ); \
  template void Scale \
  ( T alpha, AbstractDistMatrix<Base<T>>& AReal, \
             AbstractDistMatrix<Base<T>>& AImag ); \
  template void Scale \
  ( Int alpha, AbstractBlockDistMatrix<Base<T>>& AReal, \
               AbstractBlockDistMatrix<Base<T>>& AImag ); \
  template void Scale \
  ( Base<T> alpha, AbstractBlockDistMatrix<Base<T>>& AReal, \
                   AbstractBlockDistMatrix<Base<T>>& AImag ); \
  template void Scale \
  ( T alpha, AbstractBlockDistMatrix<Base<T>>& AReal, \
             AbstractBlockDistMatrix<Base<T>>& AImag );

PROTO_INT(Int);
PROTO_REAL(float);
PROTO_REAL(double);
PROTO_CPX(Complex<float>);
PROTO_CPX(Complex<double>);

} // namespace El
