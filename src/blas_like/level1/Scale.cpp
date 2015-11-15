/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename T,typename S>
void Scale( S alphaS, Matrix<T>& A )
{
    DEBUG_ONLY(CSE cse("Scale"))
    const T alpha = T(alphaS);

    T* ABuf = A.Buffer();
    const Int ALDim = A.LDim();
    const Int height = A.Height();
    const Int width = A.Width();

    // TODO: Use imatcopy if MKL or OpenBLAS is detected

    if( alpha == T(0) )
    {
        for( Int j=0; j<width; ++j )
            for( Int i=0; i<height; ++i )
                ABuf[i+j*ALDim] = 0;
    }
    else if( alpha != T(1) )
    {
        if( height >= width )
        {
            for( Int j=0; j<width; ++j )
                blas::Scal( height, alpha, &ABuf[j*ALDim], 1 );
        }
        else
        {
            for( Int i=0; i<height; ++i )
                blas::Scal( width, alpha, &ABuf[i], ALDim );
        }
    }
}

template<typename Real,typename S,typename>
void Scale( S alphaS, Matrix<Real>& AReal, Matrix<Real>& AImag )
{
    DEBUG_ONLY(CSE cse("Scale"))
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
            const Real alphaReal=alpha.real(), alphaImag=alpha.imag();
            vector<Real> aRealCopy(m);
            Real* ARealBuf = AReal.Buffer();
            Real* AImagBuf = AImag.Buffer();
            const Int ARealLDim = AReal.LDim();
            const Int AImagLDim = AImag.LDim();
            for( Int j=0; j<n; ++j )
            {
                for( Int i=0; i<m; ++i )
                    aRealCopy[i] = ARealBuf[i+j*ARealLDim];

                blas::Scal( m, alphaReal, &ARealBuf[j*ARealLDim], 1 );
                blas::Axpy
                ( m, -alphaImag, &AImagBuf[j*AImagLDim], 1, 
                                 &ARealBuf[j*ARealLDim], 1 );

                blas::Scal( m, alphaReal, &AImagBuf[j*AImagLDim], 1 );
                blas::Axpy
                ( m,  alphaImag, aRealCopy.data(),       1, 
                                 &AImagBuf[j*AImagLDim], 1 );
            }
        }
    }
}

template<typename T,typename S>
void Scale( S alpha, AbstractDistMatrix<T>& A )
{ Scale( alpha, A.Matrix() ); }

template<typename Real,typename S,typename>
void Scale( S alpha, AbstractDistMatrix<Real>& AReal, 
                     AbstractDistMatrix<Real>& AImag )
{ Scale( alpha, AReal.Matrix(), AImag.Matrix() ); }

template<typename T,typename S>
void Scale( S alpha, SparseMatrix<T>& A )
{
    if( alpha == S(0) )
    {
        const Int m = A.Height();
        const Int n = A.Width();
        A.Empty();
        A.Resize( m, n );
    }
    else if( alpha != S(1) )
    {
        T alphaT = alpha; 
        T* valueBuf = A.ValueBuffer();
        const Int numEntries = A.NumEntries();
        for( Int k=0; k<numEntries; ++k )
            valueBuf[k] *= alphaT;
    }
}

template<typename T,typename S>
void Scale( S alpha, DistSparseMatrix<T>& A )
{
    if( alpha == S(0) )
    {
        const Int m = A.Height();
        const Int n = A.Width();
        A.Empty();
        A.Resize( m, n );
    }
    else if( alpha != S(1) )
    {
        T alphaT = alpha;
        T* valueBuf = A.ValueBuffer();
        const Int numLocalEntries = A.NumLocalEntries();
        for( Int k=0; k<numLocalEntries; ++k )
            valueBuf[k] *= alphaT;
    }
}

template<typename T,typename S>
void Scale( S alpha, DistMultiVec<T>& A )
{ Scale( alpha, A.Matrix() ); }

#define PROTO_TYPES(T,S) \
  template void Scale( S alpha, Matrix<T>& A ); \
  template void Scale( S alpha, AbstractDistMatrix<T>& A ); \
  template void Scale( S alpha, SparseMatrix<T>& A ); \
  template void Scale( S alpha, DistSparseMatrix<T>& A ); \
  template void Scale( S alpha, DistMultiVec<T>& A );

#define PROTO_INT(T) PROTO_TYPES(T,T)

#define PROTO_REAL(T) \
  PROTO_TYPES(T,Int) \
  PROTO_TYPES(T,T) 

#define PROTO_COMPLEX(T) \
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
             AbstractDistMatrix<Base<T>>& AImag );

#define EL_ENABLE_QUAD
#include "El/macros/Instantiate.h"

} // namespace El
