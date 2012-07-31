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

#include "./blas-like/internal.hpp"
#include "./blas-like/level1.hpp"
#include "./blas-like/level2.hpp"
#include "./blas-like/level3.hpp"

namespace elem {

//----------------------------------------------------------------------------//
// Local BLAS-like routines: Level 1                                          //
//----------------------------------------------------------------------------//

template<typename T>
inline void
Axpy( T alpha, const Matrix<T>& X, Matrix<T>& Y )
{
#ifndef RELEASE
    PushCallStack("Axpy");
#endif
    // If X and Y are vectors, we can allow one to be a column and the other
    // to be a row. Otherwise we force X and Y to be the same dimension.
    if( (X.Height()==1 || X.Width()==1) && (Y.Height()==1 || Y.Width()==1) )
    {
        const unsigned XLength = ( X.Width()==1 ? X.Height() : X.Width() );
#ifndef RELEASE
        const unsigned YLength = ( Y.Width()==1 ? Y.Height() : Y.Width() );
        if( XLength != YLength )
            throw std::logic_error("Nonconformal Axpy");
#endif
        if( X.Width()==1 && Y.Width()==1 )
        {
            blas::Axpy
            ( XLength, alpha, X.LockedBuffer(0,0), 1, Y.Buffer(0,0), 1 );
        }
        else if( X.Width()==1 )
        {
            blas::Axpy
            ( XLength, alpha, X.LockedBuffer(0,0), 1, Y.Buffer(0,0), Y.LDim() );
        }
        else if( Y.Width()==1 )
        {
            blas::Axpy
            ( XLength, alpha, X.LockedBuffer(0,0), X.LDim(), Y.Buffer(0,0), 1 );
        }
        else
        {
            blas::Axpy
            ( XLength, alpha, 
              X.LockedBuffer(0,0), X.LDim(), Y.Buffer(0,0), Y.LDim() );
        }
    }
    else 
    {
#ifndef RELEASE
        if( X.Height() != Y.Height() || X.Width() != Y.Width() )
            throw std::logic_error("Nonconformal Axpy");
#endif
        if( X.Width() <= X.Height() )
        {
            for( int j=0; j<X.Width(); ++j )
            {
                blas::Axpy
                ( X.Height(), alpha, X.LockedBuffer(0,j), 1, Y.Buffer(0,j), 1 );
            }
        }
        else
        {
            for( int i=0; i<X.Height(); ++i )
            {
                blas::Axpy
                ( X.Width(), alpha, X.LockedBuffer(i,0), X.LDim(),
                                    Y.Buffer(i,0),       Y.LDim() );
            }
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Axpy( typename Base<T>::type alpha, const Matrix<T>& X, Matrix<T>& Y )
{ Axpy( (T)alpha, X, Y ); }

template<typename T>
inline void
Copy( const Matrix<T>& A, Matrix<T>& B )
{
#ifndef RELEASE
    PushCallStack("Copy");
#endif
    B = A;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void 
DiagonalScale
( LeftOrRight side, Orientation orientation, 
  const Matrix<T>& d, Matrix<T>& X )
{
#ifndef RELEASE
    PushCallStack("DiagonalScale");
#endif
    const int m = X.Height();
    const int n = X.Width();
    const int ldim = X.LDim();
    if( side == LEFT )
    {
        for( int i=0; i<m; ++i )
        {
            const T delta = d.Get(i,0);
            T* XBuffer = X.Buffer(i,0);
            if( orientation == ADJOINT )
                for( int j=0; j<n; ++j )
                    XBuffer[j*ldim] *= Conj(delta);
            else
                for( int j=0; j<n; ++j )
                    XBuffer[j*ldim] *= delta;
        }
    }
    else
    {
        for( int j=0; j<n; ++j )
        {
            const T delta = d.Get(j,0);
            T* XBuffer = X.Buffer(0,j);
            if( orientation == ADJOINT )
                for( int i=0; i<m; ++i )
                    XBuffer[i] *= Conj(delta);
            else
                for( int i=0; i<m; ++i )
                    XBuffer[i] *= delta;
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void 
DiagonalScale
( LeftOrRight side, Orientation orientation, 
  const Matrix<typename Base<T>::type>& d, Matrix<T>& X )
{
#ifndef RELEASE
    PushCallStack("DiagonalScale");
#endif
    typedef typename Base<T>::type R;

    const int m = X.Height();
    const int n = X.Width();
    const int ldim = X.LDim();
    if( side == LEFT )
    {
        for( int i=0; i<m; ++i )
        {
            const R delta = d.Get(i,0);
            T* XBuffer = X.Buffer(i,0);
            for( int j=0; j<n; ++j )
                XBuffer[j*ldim] *= delta;
        }
    }
    else
    {
        for( int j=0; j<n; ++j )
        {
            const R delta = d.Get(j,0);
            T* XBuffer = X.Buffer(0,j);
            for( int i=0; i<m; ++i )
                XBuffer[i] *= delta;
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void 
DiagonalSolve
( LeftOrRight side, Orientation orientation, 
  const Matrix<F>& d, Matrix<F>& X, bool checkIfSingular )
{
#ifndef RELEASE
    PushCallStack("DiagonalSolve");
#endif
    const int m = X.Height();
    const int n = X.Width();
    const int ldim = X.LDim();
    if( side == LEFT )
    {
        for( int i=0; i<m; ++i )
        {
            const F delta = d.Get(i,0);
            if( checkIfSingular && delta == (F)0 )
                throw SingularMatrixException();
            const F deltaInv = static_cast<F>(1)/delta;
            F* XBuffer = X.Buffer(i,0);
            if( orientation == ADJOINT )
                for( int j=0; j<n; ++j )
                    XBuffer[j*ldim] *= Conj(deltaInv);
            else
                for( int j=0; j<n; ++j )
                    XBuffer[j*ldim] *= deltaInv;
        }
    }
    else
    {
        for( int j=0; j<n; ++j )
        {
            const F delta = d.Get(j,0);
            if( checkIfSingular && delta == (F)0 )
                throw SingularMatrixException();
            const F deltaInv = static_cast<F>(1)/delta;
            F* XBuffer = X.Buffer(0,j);
            if( orientation == ADJOINT )
                for( int i=0; i<m; ++i )
                    XBuffer[i] *= Conj(deltaInv);
            else
                for( int i=0; i<m; ++i )
                    XBuffer[i] *= deltaInv;
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void 
DiagonalSolve
( LeftOrRight side, Orientation orientation, 
  const Matrix<typename Base<F>::type>& d, Matrix<F>& X, 
  bool checkIfSingular )
{
#ifndef RELEASE
    PushCallStack("DiagonalSolve");
#endif
    typedef typename Base<F>::type R;

    const int m = X.Height();
    const int n = X.Width();
    const int ldim = X.LDim();
    if( side == LEFT )
    {
        for( int i=0; i<m; ++i )
        {
            const R delta = d.Get(i,0);
            if( checkIfSingular && delta == (R)0 )
                throw SingularMatrixException();
            const R deltaInv = static_cast<R>(1)/delta;
            F* XBuffer = X.Buffer(i,0);
            for( int j=0; j<n; ++j )
                XBuffer[j*ldim] *= deltaInv;
        }
    }
    else
    {
        for( int j=0; j<n; ++j )
        {
            const R delta = d.Get(j,0);
            if( checkIfSingular && delta == (R)0 )
                throw SingularMatrixException();
            const R deltaInv = static_cast<R>(1)/delta;
            F* XBuffer = X.Buffer(0,j);
            for( int i=0; i<m; ++i )
                XBuffer[i] *= deltaInv;
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline T
Dot( const Matrix<T>& x, const Matrix<T>& y )
{
#ifndef RELEASE
    PushCallStack("Dot");
    if( (x.Height() != 1 && x.Width() != 1) ||
        (y.Height() != 1 && y.Width() != 1) )
        throw std::logic_error("Expected vector inputs");
    int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
    int yLength = ( y.Width() == 1 ? y.Height() : y.Width() );
    if( xLength != yLength )
        throw std::logic_error("x and y must be the same length");
#endif
    T dotProduct;
    if( x.Width() == 1 && y.Width() == 1 )
    {
        dotProduct = blas::Dot
                     ( x.Height(), x.LockedBuffer(), 1,
                                   y.LockedBuffer(), 1 );
    }
    else if( x.Width() == 1 )
    {
        dotProduct = blas::Dot
                     ( x.Height(), x.LockedBuffer(), 1,
                                   y.LockedBuffer(), y.LDim() );
    }
    else if( y.Width() == 1 )
    {
        dotProduct = blas::Dot
                     ( x.Width(), x.LockedBuffer(), x.LDim(),
                                  y.LockedBuffer(), 1        );
    }
    else
    {
        dotProduct = blas::Dot
                     ( x.Width(), x.LockedBuffer(), x.LDim(),
                                  y.LockedBuffer(), y.LDim() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return dotProduct;
}

template<typename T>
inline T
Dotc( const Matrix<T>& x, const Matrix<T>& y )
{
#ifndef RELEASE
    PushCallStack("Dotc");
    if( (x.Height() != 1 && x.Width() != 1) ||
        (y.Height() != 1 && y.Width() != 1) )
        throw std::logic_error("Expected vector inputs");
    int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
    int yLength = ( y.Width() == 1 ? y.Height() : y.Width() );
    if( xLength != yLength )
        throw std::logic_error("x and y must be the same length");
#endif
    T dotProduct;
    if( x.Width() == 1 && y.Width() == 1 )
    {
        dotProduct = blas::Dotc
                     ( x.Height(), x.LockedBuffer(), 1,
                                   y.LockedBuffer(), 1 );
    }
    else if( x.Width() == 1 )
    {
        dotProduct = blas::Dotc
                     ( x.Height(), x.LockedBuffer(), 1,
                                   y.LockedBuffer(), y.LDim() );
    }
    else if( y.Width() == 1 )
    {
        dotProduct = blas::Dotc
                     ( x.Width(), x.LockedBuffer(), x.LDim(),
                                  y.LockedBuffer(), 1        );
    }
    else
    {
        dotProduct = blas::Dotc
                     ( x.Width(), x.LockedBuffer(), x.LDim(),
                                  y.LockedBuffer(), y.LDim() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return dotProduct;
}

template<typename T>
inline T
Dotu( const Matrix<T>& x, const Matrix<T>& y )
{
#ifndef RELEASE
    PushCallStack("Dotu");
    if( (x.Height() != 1 && x.Width() != 1) ||
        (y.Height() != 1 && y.Width() != 1) )
        throw std::logic_error("Expected vector inputs");
    int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
    int yLength = ( y.Width() == 1 ? y.Height() : y.Width() );
    if( xLength != yLength )
        throw std::logic_error("x and y must be the same length");
#endif
    T dotProduct;
    if( x.Width() == 1 && y.Width() == 1 )
    {
        dotProduct = blas::Dotu
                     ( x.Height(), x.LockedBuffer(), 1,
                                   y.LockedBuffer(), 1 );
    }
    else if( x.Width() == 1 )
    {
        dotProduct = blas::Dotu
                     ( x.Height(), x.LockedBuffer(), 1,
                                   y.LockedBuffer(), y.LDim() );
    }
    else if( y.Width() == 1 )
    {
        dotProduct = blas::Dotu
                     ( x.Width(), x.LockedBuffer(), x.LDim(),
                                  y.LockedBuffer(), 1        );
    }
    else
    {
        dotProduct = blas::Dotu
                     ( x.Width(), x.LockedBuffer(), x.LDim(),
                                  y.LockedBuffer(), y.LDim() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return dotProduct;
}

template<typename F>
inline typename Base<F>::type
Nrm2( const Matrix<F>& x )
{
#ifndef RELEASE
    PushCallStack("Nrm2");
    if( x.Height() != 1 && x.Width() != 1 )
        throw std::logic_error("Expected vector input");
#endif
    typedef typename Base<F>::type R;

    R norm;
    if( x.Width() == 1 )
        norm = blas::Nrm2( x.Height(), x.LockedBuffer(), 1 );
    else
        norm = blas::Nrm2( x.Width(), x.LockedBuffer(), x.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

template<typename T>
inline void
Scale( T alpha, Matrix<T>& X )
{
#ifndef RELEASE
    PushCallStack("Scale");
#endif
    if( alpha != (T)1 )
    {
        if( alpha == (T)0 )
            for( int j=0; j<X.Width(); ++j )
                for( int i=0; i<X.Height(); ++i )
                    X.Set(i,j,0);
        else
            for( int j=0; j<X.Width(); ++j )
                blas::Scal( X.Height(), alpha, X.Buffer(0,j), 1 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Scale( typename Base<T>::type alpha, Matrix<T>& X )
{ Scale( (T)alpha, X ); }

template<typename T>
inline void
Scal( T alpha, Matrix<T>& X )
{ Scale( alpha, X ); }

template<typename T>
inline void
Scal( typename Base<T>::type alpha, Matrix<T>& X )
{ Scale( (T)alpha, X ); }

//----------------------------------------------------------------------------//
// Local BLAS-like routines: Level 1 (extensions)                             //
//----------------------------------------------------------------------------//

template<typename T>
inline void
Adjoint( const Matrix<T>& A, Matrix<T>& B )
{
#ifndef RELEASE
    PushCallStack("Adjoint");
#endif
    const int m = A.Height();
    const int n = A.Width();
    if( B.Viewing() )
    {
        if( B.Height() != n || B.Width() != m )
        {
            std::ostringstream msg;
            msg << "If Adjoint'ing into a view, it must be the right size:\n"
                << "  A ~ " << A.Height() << " x " << A.Width() << "\n"
                << "  B ~ " << B.Height() << " x " << B.Width();
            throw std::logic_error( msg.str().c_str() );
        }
    }
    else
        B.ResizeTo( n, m );

    for( int j=0; j<n; ++j )
        for( int i=0; i<m; ++i )
            B.Set(j,i,Conj(A.Get(i,j)));
#ifndef RELEASE
    PopCallStack();
#endif
}

// Default case is for real datatypes
template<typename Z>
inline void
Conjugate( Matrix<Z>& A )
{ }

// Specialization is to complex datatypes
template<typename Z>
inline void
Conjugate( Matrix<Complex<Z> >& A )
{
#ifndef RELEASE
    PushCallStack("Conjugate (in-place)");
#endif
    const int m = A.Height();
    const int n = A.Width();
    for( int j=0; j<n; ++j )
        for( int i=0; i<m; ++i )
            A.Set(i,j,Conj(A.Get(i,j)));
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Conjugate( const Matrix<T>& A, Matrix<T>& B )
{
#ifndef RELEASE
    PushCallStack("Conjugate");
#endif
    const int m = A.Height();
    const int n = A.Width();
    B.ResizeTo( m, n );
    for( int j=0; j<n; ++j )
        for( int i=0; i<m; ++i )
            B.Set(i,j,Conj(A.Get(i,j)));
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
MakeReal( Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("MakeReal");
#endif
    T* ABuffer = A.Buffer();
    const int height = A.Height();
    const int width = A.Width();
    const int ldim = A.LDim();
    for( int j=0; j<width; ++j )
        for( int i=0; i<height; ++i )
            ABuffer[i+j*ldim] = RealPart(ABuffer[i+j*ldim]);
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V>
inline void
MakeReal( DistMatrix<T,U,V>& A )
{
#ifndef RELEASE
    PushCallStack("MakeReal");
#endif
    MakeReal( A.LocalMatrix() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
MakeTrapezoidal
( LeftOrRight side, UpperOrLower uplo, int offset, Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("MakeTrapezoidal");
#endif
    const int height = A.Height();
    const int width = A.Width();
    const int ldim = A.LDim();
    T* buffer = A.Buffer();

    if( uplo == LOWER )
    {
        if( side == LEFT )
        {
#ifdef HAVE_OPENMP
            #pragma omp parallel for
#endif
            for( int j=std::max(0,offset+1); j<width; ++j )
            {
                const int lastZeroRow = j-offset-1;
                const int numZeroRows = std::min( lastZeroRow+1, height );
                MemZero( &buffer[j*ldim], numZeroRows );
            }
        }
        else
        {
#ifdef HAVE_OPENMP
            #pragma omp parallel for
#endif
            for( int j=std::max(0,offset-height+width+1); j<width; ++j )
            {
                const int lastZeroRow = j-offset+height-width-1;
                const int numZeroRows = std::min( lastZeroRow+1, height );
                MemZero( &buffer[j*ldim], numZeroRows );
            }
        }
    }
    else
    {
        if( side == LEFT )
        {
#ifdef HAVE_OPENMP
            #pragma omp parallel for
#endif
            for( int j=0; j<width; ++j )
            {
                const int firstZeroRow = std::max(j-offset+1,0);
                if( firstZeroRow < height )
                    MemZero
                    ( &buffer[firstZeroRow+j*ldim], height-firstZeroRow );
            }
        }
        else
        {
#ifdef HAVE_OPENMP
            #pragma omp parallel for
#endif
            for( int j=0; j<width; ++j )
            {
                const int firstZeroRow = std::max(j-offset+height-width+1,0);
                if( firstZeroRow < height )
                    MemZero
                    ( &buffer[firstZeroRow+j*ldim], height-firstZeroRow );
            }
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
MakeHermitian( UpperOrLower uplo, Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("MakeHermitian");
#endif
    if( A.Height() != A.Width() )
        throw std::logic_error("Cannot make non-square matrix Hermitian");

    Matrix<T> d;
    A.GetDiagonal( d );
    MakeReal( d );

    if( uplo == LOWER )
        MakeTrapezoidal( LEFT, LOWER, -1, A );
    else
        MakeTrapezoidal( LEFT, UPPER, +1, A );
    Matrix<T> AAdj;
    Adjoint( A, AAdj );
    Axpy( (T)1, AAdj, A );

    A.SetDiagonal( d );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
MakeHermitian( UpperOrLower uplo, DistMatrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("MakeHermitian");
#endif
    const Grid& g = A.Grid();
    if( A.Height() != A.Width() )
        throw std::logic_error("Cannot make non-square matrix Hermitian");

    DistMatrix<T,MD,STAR> d(g);
    A.GetDiagonal( d );
    MakeReal( d );

    if( uplo == LOWER )
        MakeTrapezoidal( LEFT, LOWER, -1, A );
    else
        MakeTrapezoidal( LEFT, UPPER, +1, A );
    DistMatrix<T> AAdj(g);
    Adjoint( A, AAdj );
    Axpy( (T)1, AAdj, A );

    A.SetDiagonal( d );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
MakeSymmetric( UpperOrLower uplo, Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("MakeSymmetric");
#endif
    if( A.Height() != A.Width() )
        throw std::logic_error("Cannot make non-square matrix symmetric");

    Matrix<T> d;
    A.GetDiagonal( d );

    if( uplo == LOWER )
        MakeTrapezoidal( LEFT, LOWER, -1, A );
    else
        MakeTrapezoidal( LEFT, UPPER, +1, A );
    Matrix<T> ATrans;
    Transpose( A, ATrans );
    Axpy( (T)1, ATrans, A );

    A.SetDiagonal( d );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
MakeSymmetric( UpperOrLower uplo, DistMatrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("MakeSymmetric");
#endif
    if( A.Height() != A.Width() )
        throw std::logic_error("Cannot make non-square matrix symmetric");

    const Grid& g = A.Grid();
    DistMatrix<T> d(g);
    A.GetDiagonal( d );

    if( uplo == LOWER )
        MakeTrapezoidal( LEFT, LOWER, -1, A );
    else
        MakeTrapezoidal( LEFT, UPPER, +1, A );
    DistMatrix<T> ATrans(g);
    Transpose( A, ATrans );
    Axpy( (T)1, ATrans, A );

    A.SetDiagonal( d );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
ScaleTrapezoid
( T alpha, LeftOrRight side, UpperOrLower uplo, int offset, Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("ScaleTrapezoid");
#endif
    const int height = A.Height();
    const int width = A.Width();
    const int ldim = A.LDim();
    T* buffer = A.Buffer();

    if( uplo == UPPER )
    {
        if( side == LEFT )
        {
#ifdef HAVE_OPENMP
            #pragma omp parallel for
#endif
            for( int j=std::max(0,offset-1); j<width; ++j )
            {
                const int numRows = j-offset+1;
                for( int i=0; i<numRows; ++i )
                    buffer[i+j*ldim] *= alpha;
            }
        }
        else
        {
#ifdef HAVE_OPENMP
            #pragma omp parallel for
#endif
            for( int j=std::max(0,offset-height+width-1); j<width; ++j )
            {
                const int numRows = j-offset+height-width+1;
                for( int i=0; i<numRows; ++i )
                    buffer[i+j*ldim] *= alpha;
            }
        }
    }
    else
    {
        if( side == LEFT )
        {
#ifdef HAVE_OPENMP
            #pragma omp parallel for
#endif
            for( int j=0; j<width; ++j )
            {
                const int numZeroRows = std::max(j-offset,0);
                for( int i=numZeroRows; i<height; ++i )
                    buffer[i+j*ldim] *= alpha;
            }
        }
        else
        {
#ifdef HAVE_OPENMP
            #pragma omp parallel for
#endif
            for( int j=0; j<width; ++j )
            {
                const int numZeroRows = std::max(j-offset+height-width,0);
                for( int i=numZeroRows; i<height; ++i )
                    buffer[i+j*ldim] *= alpha;
            }
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Transpose( const Matrix<T>& A, Matrix<T>& B )
{
#ifndef RELEASE
    PushCallStack("Transpose");
#endif
    const int m = A.Height();
    const int n = A.Width();
    if( B.Viewing() )
    {
        if( B.Height() != n || B.Width() != m )
        {
            std::ostringstream msg;
            msg << "If Transpose'ing into a view, it must be the right size:\n"
                << "  A ~ " << A.Height() << " x " << A.Width() << "\n"
                << "  B ~ " << B.Height() << " x " << B.Width();
            throw std::logic_error( msg.str().c_str() );
        }
    }
    else
        B.ResizeTo( n, m );

    for( int j=0; j<n; ++j )
        for( int i=0; i<m; ++i )
            B.Set(j,i,A.Get(i,j));
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Zero( Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("Zero");
#endif
    const int height = A.Height();
    const int width = A.Width();
#ifdef HAVE_OPENMP
    #pragma omp parallel for
#endif
    for( int j=0; j<width; ++j )
        MemZero( A.Buffer(0,j), height );
#ifndef RELEASE
    PopCallStack();
#endif
}

//----------------------------------------------------------------------------//
// Local BLAS-like routines: Level 2                                          //
//----------------------------------------------------------------------------//

template<typename T>
inline void
Gemv
( Orientation orientation,
  T alpha, const Matrix<T>& A, const Matrix<T>& x, T beta, Matrix<T>& y )
{
#ifndef RELEASE
    PushCallStack("Gemv");
    if( ( x.Height() != 1 && x.Width() != 1 ) ||
        ( y.Height() != 1 && y.Width() != 1 ) )
    {
        std::ostringstream msg;
        msg << "x and y must be vectors:\n"
            << "  x ~ " << x.Height() << " x " << x.Width() << "\n"
            << "  y ~ " << y.Height() << " x " << y.Width();
        throw std::logic_error( msg.str().c_str() );
    }
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    const int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
    if( orientation == NORMAL )
    {
        if( A.Height() != yLength || A.Width() != xLength )
        {
            std::ostringstream msg;
            msg << "A must conform with x and y:\n"
                << "  A ~ " << A.Height() << " x " << A.Width() << "\n"
                << "  x ~ " << x.Height() << " x " << x.Width() << "\n"
                << "  y ~ " << y.Height() << " x " << y.Width();
            throw std::logic_error( msg.str().c_str() );
        }
    }
    else
    {
        if( A.Width() != yLength || A.Height() != xLength )
        {
            std::ostringstream msg;
            msg << "A must conform with x and y:\n"
                << "  A ~ " << A.Height() << " x " << A.Width() << "\n"
                << "  x ~ " << x.Height() << " x " << x.Width() << "\n"
                << "  y ~ " << y.Height() << " x " << y.Width();
            throw std::logic_error( msg.str().c_str() );
        }
    }
#endif
    const char transChar = OrientationToChar( orientation );
    const int m = A.Height();
    const int n = A.Width();
    const int k = ( transChar == 'N' ? n : m );
    const int incx = ( x.Width()==1 ? 1 : x.LDim() );
    const int incy = ( y.Width()==1 ? 1 : y.LDim() );
    if( k != 0 )
    {
        blas::Gemv
        ( transChar, m, n, 
          alpha, A.LockedBuffer(), A.LDim(), x.LockedBuffer(), incx, 
          beta,  y.Buffer(), incy );
    }
    else
    {
        Scale( beta, y );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Ger( T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("Ger");
    if( ( x.Height() != 1 && x.Width() != 1 ) ||
        ( y.Height() != 1 && y.Width() != 1 ) )
        throw std::logic_error("x and y must be vectors");
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    const int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
    if( xLength != A.Height() || yLength != A.Width() )
    {
        std::ostringstream msg;
        msg << "Nonconformal Ger:\n"
            << "  x ~ " << x.Height() << " x " << x.Width() << "\n"
            << "  y ~ " << y.Height() << " x " << y.Width() << "\n"
            << "  A ~ " << A.Height() << " x " << A.Width();
        throw std::logic_error( msg.str().c_str() );
    }
#endif
    const int m = A.Height(); 
    const int n = A.Width();
    const int incx = ( x.Width()==1 ? 1 : x.LDim() );
    const int incy = ( y.Width()==1 ? 1 : y.LDim() );
    blas::Ger
    ( m, n, alpha, x.LockedBuffer(), incx, y.LockedBuffer(), incy, 
                   A.Buffer(), A.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Gerc( T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("Gerc");
    if( ( x.Height() != 1 && x.Width() != 1 ) ||
        ( y.Height() != 1 && y.Width() != 1 ) )
        throw std::logic_error("x and y must be vectors");
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    const int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
    if( xLength != A.Height() || yLength != A.Width() )
        throw std::logic_error("Nonconformal Gerc");
#endif
    const int m = A.Height(); 
    const int n = A.Width();
    const int incx = ( x.Width()==1 ? 1 : x.LDim() );
    const int incy = ( y.Width()==1 ? 1 : y.LDim() );
    blas::Gerc
    ( m, n, alpha, x.LockedBuffer(), incx, y.LockedBuffer(), incy, 
                   A.Buffer(), A.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Geru( T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("Geru");
    if( ( x.Height() != 1 && x.Width() != 1 ) ||
        ( y.Height() != 1 && y.Width() != 1 ) )
        throw std::logic_error("x and y must be vectors");
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    const int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
    if( xLength != A.Height() || yLength != A.Width() )
        throw std::logic_error("Nonconformal Geru");
#endif
    const int m = A.Height(); 
    const int n = A.Width();
    const int incx = ( x.Width()==1 ? 1 : x.LDim() );
    const int incy = ( y.Width()==1 ? 1 : y.LDim() );
    blas::Geru
    ( m, n, alpha, x.LockedBuffer(), incx, y.LockedBuffer(), incy, 
                   A.Buffer(), A.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Hemv
( UpperOrLower uplo,
  T alpha, const Matrix<T>& A, const Matrix<T>& x, T beta, Matrix<T>& y )
{
#ifndef RELEASE
    PushCallStack("Hemv");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    if( ( x.Height() != 1 && x.Width() != 1 ) ||
        ( y.Height() != 1 && y.Width() != 1 ) )
        throw std::logic_error("x and y must be vectors");
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    const int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
    if( A.Height() != xLength || A.Height() != yLength )
        throw std::logic_error("A must conform with x and y");
#endif
    const char uploChar = UpperOrLowerToChar( uplo );
    const int m = A.Height();
    const int incx = ( x.Width()==1 ? 1 : x.LDim() );
    const int incy = ( y.Width()==1 ? 1 : y.LDim() );
    blas::Hemv
    ( uploChar, m,
      alpha, A.LockedBuffer(), A.LDim(), x.LockedBuffer(), incx,
      beta,  y.Buffer(), incy );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Her( UpperOrLower uplo, T alpha, const Matrix<T>& x, Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("Her");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    if( x.Width() != 1 && x.Height() != 1 )
        throw std::logic_error("x must be a vector");
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    if( xLength != A.Height() )
        throw std::logic_error("x must conform with A");
#endif
    const char uploChar = UpperOrLowerToChar( uplo );
    const int m = A.Height();
    const int incx = ( x.Width()==1 ? 1 : x.LDim() );
    blas::Her
    ( uploChar, m, alpha, x.LockedBuffer(), incx, A.Buffer(), A.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Her2
( UpperOrLower uplo,
  T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("Her2");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    if( (x.Width() != 1 && x.Height() != 1) || 
        (y.Width() != 1 && y.Height() != 1) )
        throw std::logic_error("x and y must be vectors");
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    const int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
    if( xLength != A.Height() || yLength != A.Height() )
        throw std::logic_error("x and y must conform with A");
#endif
    const char uploChar = UpperOrLowerToChar( uplo );
    const int m = A.Height();
    const int incx = ( x.Width()==1 ? 1 : x.LDim() );
    const int incy = ( y.Width()==1 ? 1 : y.LDim() );
    blas::Her2
    ( uploChar, m,
      alpha, x.LockedBuffer(), incx, y.LockedBuffer(), incy,
             A.Buffer(), A.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Symv
( UpperOrLower uplo,
  T alpha, const Matrix<T>& A, const Matrix<T>& x, T beta, Matrix<T>& y )
{
#ifndef RELEASE
    PushCallStack("Symv");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    if( ( x.Height() != 1 && x.Width() != 1 ) ||
        ( y.Height() != 1 && y.Width() != 1 ) )
        throw std::logic_error("x and y must be vectors");
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    const int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
    if( A.Height() != xLength || A.Height() != yLength )
        throw std::logic_error("A must conform with x and y");
#endif
    const char uploChar = UpperOrLowerToChar( uplo );
    const int m = A.Height();
    const int incx = ( x.Width()==1 ? 1 : x.LDim() );
    const int incy = ( y.Width()==1 ? 1 : y.LDim() );
    blas::Symv
    ( uploChar, m, 
      alpha, A.LockedBuffer(), A.LDim(), x.LockedBuffer(), incx, 
      beta,  y.Buffer(), incy );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Syr( UpperOrLower uplo, T alpha, const Matrix<T>& x, Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("Syr");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    if( x.Width() != 1 && x.Height() != 1 )
        throw std::logic_error("x must be a vector");
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    if( xLength != A.Height() )
        throw std::logic_error("x must conform with A");
#endif
    const char uploChar = UpperOrLowerToChar( uplo );
    const int m = A.Height();
    const int incx = ( x.Width()==1 ? 1 : x.LDim() );
    blas::Syr
    ( uploChar, m, alpha, x.LockedBuffer(), incx, A.Buffer(), A.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Syr2
( UpperOrLower uplo, 
  T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("Syr2");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    if( (x.Width() != 1 && x.Height() != 1) || 
        (y.Width() != 1 && y.Height() != 1) )
        throw std::logic_error("x and y must be vectors");
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    const int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
    if( xLength != A.Height() || yLength != A.Height() )
        throw std::logic_error("x and y must conform with A");
#endif
    const char uploChar = UpperOrLowerToChar( uplo );
    const int m = A.Height();
    const int incx = ( x.Width()==1 ? 1 : x.LDim() );
    const int incy = ( y.Width()==1 ? 1 : y.LDim() );
    blas::Syr2
    ( uploChar, m,
      alpha, x.LockedBuffer(), incx, y.LockedBuffer(), incy,
             A.Buffer(), A.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Trmv
( UpperOrLower uplo, Orientation orientation, UnitOrNonUnit diag,
  const Matrix<T>& A, Matrix<T>& x )
{
#ifndef RELEASE
    PushCallStack("Trmv");
    if( x.Height() != 1 && x.Width() != 1 )
        throw std::logic_error("x must be a vector");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    if( xLength != A.Height() )
        throw std::logic_error("x must conform with A");
#endif
    const char uploChar = UpperOrLowerToChar( uplo );
    const char transChar = OrientationToChar( orientation );
    const char diagChar = UnitOrNonUnitToChar( diag );
    const int m = A.Height();
    const int incx = ( x.Width()==1 ? 1 : x.LDim() );
    blas::Trmv
    ( uploChar, transChar, diagChar, m,
      A.LockedBuffer(), A.LDim(), x.Buffer(), incx );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
Trsv
( UpperOrLower uplo, Orientation orientation, UnitOrNonUnit diag,
  const Matrix<F>& A, Matrix<F>& x )
{
#ifndef RELEASE
    PushCallStack("Trsv");
    if( x.Height() != 1 && x.Width() != 1 )
        throw std::logic_error("x must be a vector");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    if( xLength != A.Height() )
        throw std::logic_error("x must conform with A");
#endif
    const char uploChar = UpperOrLowerToChar( uplo );
    const char transChar = OrientationToChar( orientation );
    const char diagChar = UnitOrNonUnitToChar( diag );
    const int m = A.Height();
    const int incx = ( x.Width()==1 ? 1 : x.LDim() );
    blas::Trsv
    ( uploChar, transChar, diagChar, m,
      A.LockedBuffer(), A.LDim(), x.Buffer(), incx );
#ifndef RELEASE
    PopCallStack();
#endif
}

//----------------------------------------------------------------------------//
// Local BLAS-like routines: Level 3                                          //
//----------------------------------------------------------------------------//

template<typename T>
inline void
Gemm
( Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const Matrix<T>& A, const Matrix<T>& B, T beta, Matrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("Gemm");
    if( orientationOfA == NORMAL && orientationOfB == NORMAL )
    {
        if( A.Height() != C.Height() ||
            B.Width()  != C.Width()  ||
            A.Width()  != B.Height() )
            throw std::logic_error("Nonconformal GemmNN");
    }
    else if( orientationOfA == NORMAL )
    {
        if( A.Height() != C.Height() ||
            B.Height() != C.Width()  ||
            A.Width()  != B.Width() )
            throw std::logic_error("Nonconformal GemmN(T/C)");
    }
    else if( orientationOfB == NORMAL )
    {
        if( A.Width()  != C.Height() ||
            B.Width()  != C.Width()  ||
            A.Height() != B.Height() )
            throw std::logic_error("Nonconformal Gemm(T/C)N");
    }
    else
    {
        if( A.Width()  != C.Height() ||
            B.Height() != C.Width()  ||
            A.Height() != B.Width() )
            throw std::logic_error("Nonconformal Gemm(T/C)(T/C)");
    }
#endif
    const char transA = OrientationToChar( orientationOfA );
    const char transB = OrientationToChar( orientationOfB );
    const int m = C.Height();
    const int n = C.Width();
    const int k = ( orientationOfA == NORMAL ? A.Width() : A.Height() );
    if( k != 0 )
    {
        blas::Gemm
        ( transA, transB, m, n, k, 
          alpha, A.LockedBuffer(), A.LDim(), B.LockedBuffer(), B.LDim(),
          beta,  C.Buffer(),       C.LDim() );
    }
    else
    {
        Scale( beta, C );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Hemm
( LeftOrRight side, UpperOrLower uplo,
  T alpha, const Matrix<T>& A, const Matrix<T>& B, T beta, Matrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("Hemm");
#endif
    const char sideChar = LeftOrRightToChar( side );
    const char uploChar = UpperOrLowerToChar( uplo );
    blas::Hemm
    ( sideChar, uploChar, C.Height(), C.Width(),
      alpha, A.LockedBuffer(), A.LDim(), 
             B.LockedBuffer(), B.LDim(),
      beta,  C.Buffer(),       C.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Her2k
( UpperOrLower uplo, Orientation orientation,
  T alpha, const Matrix<T>& A, const Matrix<T>& B, T beta, Matrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("Her2k");
    if( orientation == NORMAL )
    {
        if( A.Height() != C.Height() || A.Height() != C.Width() ||
            B.Height() != C.Height() ||B.Height() != C.Width() )
            throw std::logic_error("Nonconformal Her2k");
    }
    else if( orientation == ADJOINT )
    {
        if( A.Width() != C.Height() || A.Width() != C.Width() ||
            B.Width() != C.Height() || B.Width() != C.Width() )
            throw std::logic_error("Nonconformal Her2k");
    }
    else
        throw std::logic_error
        ("Her2k only accepts NORMAL and ADJOINT options");
#endif
    const char uploChar = UpperOrLowerToChar( uplo );
    const char transChar = OrientationToChar( orientation );
    const int k = ( orientation == NORMAL ? A.Width() : A.Height() );
    blas::Her2k
    ( uploChar, transChar, C.Height(), k, 
      alpha, A.LockedBuffer(), A.LDim(), 
             B.LockedBuffer(), B.LDim(),
      beta,  C.Buffer(),       C.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Herk
( UpperOrLower uplo, Orientation orientation,
  T alpha, const Matrix<T>& A, T beta, Matrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("Herk");
    if( orientation == NORMAL )
    {
        if( A.Height() != C.Height() || A.Height() != C.Width() )
            throw std::logic_error("Nonconformal Herk");
    }
    else if( orientation == ADJOINT )
    {
        if( A.Width() != C.Height() || A.Width() != C.Width() )
            throw std::logic_error("Nonconformal Herk");
    }
    else
        throw std::logic_error("Herk only accepts NORMAL and ADJOINT options.");
#endif
    const char uploChar = UpperOrLowerToChar( uplo );
    const char transChar = OrientationToChar( orientation );
    const int k = ( orientation == NORMAL ? A.Width() : A.Height() );
    blas::Herk
    ( uploChar, transChar, C.Height(), k, 
      alpha, A.LockedBuffer(), A.LDim(), 
      beta,  C.Buffer(),       C.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Symm
( LeftOrRight side, UpperOrLower uplo,
  T alpha, const Matrix<T>& A, const Matrix<T>& B, T beta, Matrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("Symm");
#endif
    const char sideChar = LeftOrRightToChar( side );
    const char uploChar = UpperOrLowerToChar( uplo );
    blas::Symm
    ( sideChar, uploChar, C.Height(), C.Width(),
      alpha, A.LockedBuffer(), A.LDim(), 
             B.LockedBuffer(), B.LDim(),
      beta,  C.Buffer(),       C.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Syr2k
( UpperOrLower uplo, Orientation orientation,
  T alpha, const Matrix<T>& A, const Matrix<T>& B, T beta, Matrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("Syr2k");
    if( orientation == NORMAL )
    {
        if( A.Height() != C.Height() || A.Height() != C.Width() ||
            B.Height() != C.Height() ||B.Height() != C.Width()    )
            throw std::logic_error("Nonconformal Syr2k");
    }
    else if( orientation == TRANSPOSE )
    {
        if( A.Width() != C.Height() || A.Width() != C.Width() ||
            B.Width() != C.Height() || B.Width() != C.Width()   )
            throw std::logic_error("Nonconformal Syr2k");
    }
    else
        throw std::logic_error
        ("Syr2k only accepts NORMAL and TRANSPOSE options");
#endif
    const char uploChar = UpperOrLowerToChar( uplo );
    const char transChar = OrientationToChar( orientation );
    const int k = ( orientation == NORMAL ? A.Width() : A.Height() );
    blas::Syr2k
    ( uploChar, transChar, C.Height(), k, 
      alpha, A.LockedBuffer(), A.LDim(), 
             B.LockedBuffer(), B.LDim(),
      beta,  C.Buffer(),       C.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Syrk
( UpperOrLower uplo, Orientation orientation,
  T alpha, const Matrix<T>& A, T beta, Matrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("Syrk");
    if( orientation == NORMAL )
    {
        if( A.Height() != C.Height() || A.Height() != C.Width() )
            throw std::logic_error("Nonconformal Syrk");
    }
    else if( orientation == TRANSPOSE )
    {
        if( A.Width() != C.Height() || A.Width() != C.Width() )
            throw std::logic_error("Nonconformal Syrk");
    }
    else
        throw std::logic_error
        ("Syrk only accepts NORMAL and TRANSPOSE options");
#endif
    const char uploChar = UpperOrLowerToChar( uplo );
    const char transChar = OrientationToChar( orientation );
    const int k = ( orientation == NORMAL ? A.Width() : A.Height() );
    blas::Syrk
    ( uploChar, transChar, C.Height(), k, 
      alpha, A.LockedBuffer(), A.LDim(), 
      beta,  C.Buffer(),       C.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Trmm
( LeftOrRight side, UpperOrLower uplo, 
  Orientation orientation, UnitOrNonUnit diag,
  T alpha, const Matrix<T>& A, Matrix<T>& B )
{
#ifndef RELEASE
    PushCallStack("Trmm");
    if( A.Height() != A.Width() )
        throw std::logic_error("Triangular matrix must be square");
    if( side == LEFT )
    {
        if( A.Height() != B.Height() )
            throw std::logic_error("Nonconformal Trmm");
    }
    else
    {
        if( A.Height() != B.Width() )
            throw std::logic_error("Nonconformal Trmm");
    }
#endif
    const char sideChar = LeftOrRightToChar( side );
    const char uploChar = UpperOrLowerToChar( uplo );
    const char transChar = OrientationToChar( orientation );
    const char diagChar = UnitOrNonUnitToChar( diag );
    blas::Trmm
    ( sideChar, uploChar, transChar, diagChar, B.Height(), B.Width(),
      alpha, A.LockedBuffer(), A.LDim(), B.Buffer(), B.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
Trsm
( LeftOrRight side, UpperOrLower uplo,
  Orientation orientation, UnitOrNonUnit diag,
  F alpha, const Matrix<F>& A, Matrix<F>& B,
  bool checkIfSingular )
{
#ifndef RELEASE
    PushCallStack("Trsm");
    if( A.Height() != A.Width() )
        throw std::logic_error("Triangular matrix must be square");
    if( side == LEFT )
    {
        if( A.Height() != B.Height() )
            throw std::logic_error("Nonconformal Trsm");
    }
    else
    {
        if( A.Height() != B.Width() )
            throw std::logic_error("Nonconformal Trsm");
    }
#endif
    const char sideChar = LeftOrRightToChar( side );
    const char uploChar = UpperOrLowerToChar( uplo );
    const char transChar = OrientationToChar( orientation );
    const char diagChar = UnitOrNonUnitToChar( diag );
    if( checkIfSingular && diag != UNIT )
    {
        const int n = A.Height();
        for( int j=0; j<n; ++j )
            if( A.Get(j,j) == (F)0 )
                throw SingularMatrixException();
    }
    blas::Trsm
    ( sideChar, uploChar, transChar, diagChar, B.Height(), B.Width(),
      alpha, A.LockedBuffer(), A.LDim(), B.Buffer(), B.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

//----------------------------------------------------------------------------//
// Distributed BLAS-like routines: Level 1                                    //
//----------------------------------------------------------------------------//

template<typename T,Distribution U,Distribution V>
inline void
Axpy( T alpha, const DistMatrix<T,U,V>& X, DistMatrix<T,U,V>& Y )
{
#ifndef RELEASE
    PushCallStack("Axpy");
    if( X.Grid() != Y.Grid() )
        throw std::logic_error
        ("X and Y must be distributed over the same grid");
#endif
    if( X.ColAlignment() == Y.ColAlignment() &&
        X.RowAlignment() == Y.RowAlignment() )
    {
        Axpy( alpha, X.LockedLocalMatrix(), Y.LocalMatrix() );
    }
    else
    {
        DistMatrix<T,U,V> XCopy( X.Grid() );
        XCopy.AlignWith( Y );
        XCopy = X;
        Axpy( alpha, XCopy.LockedLocalMatrix(), Y.LocalMatrix() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V>
inline void
Axpy
( typename Base<T>::type alpha, 
  const DistMatrix<T,U,V>& X, DistMatrix<T,U,V>& Y )
{ Axpy( (T)alpha, X, Y ); }

template<typename T,Distribution U,Distribution V,
                    Distribution W,Distribution Z>
inline void
Copy( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B )
{
#ifndef RELEASE
    PushCallStack("Copy");
#endif
    B = A;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V,
                    Distribution W,Distribution Z>
inline void
DiagonalScale
( LeftOrRight side, Orientation orientation, 
  const DistMatrix<T,U,V>& d, DistMatrix<T,W,Z>& X )
{
#ifndef RELEASE
    PushCallStack("DiagonalScale");
#endif
    if( side == LEFT )
    {
        if( U == W && V == STAR && d.ColAlignment() == X.ColAlignment() )
        {
            DiagonalScale
            ( LEFT, orientation, d.LockedLocalMatrix(), X.LocalMatrix() );
        }
        else
        {
            DistMatrix<T,W,STAR> d_W_STAR( X.Grid() );
            d_W_STAR = d;
            DiagonalScale
            ( LEFT, orientation, 
              d_W_STAR.LockedLocalMatrix(), X.LocalMatrix() );
        }
    }
    else
    {
        if( U == Z && V == STAR && d.ColAlignment() == X.RowAlignment() )
        {
            DiagonalScale
            ( RIGHT, orientation, d.LockedLocalMatrix(), X.LocalMatrix() );
        }
        else
        {
            DistMatrix<T,Z,STAR> d_Z_STAR( X.Grid() );
            d_Z_STAR = d;
            DiagonalScale
            ( RIGHT, orientation, 
              d_Z_STAR.LockedLocalMatrix(), X.LocalMatrix() );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V,
                    Distribution W,Distribution Z>
inline void
DiagonalScale
( LeftOrRight side, Orientation orientation, 
  const DistMatrix<typename Base<T>::type,U,V>& d, DistMatrix<T,W,Z>& X )
{
#ifndef RELEASE
    PushCallStack("DiagonalScale");
#endif
    typedef typename Base<T>::type R;

    if( side == LEFT )
    {
        if( U == W && V == STAR && d.ColAlignment() == X.ColAlignment() )
        {
            DiagonalScale
            ( LEFT, orientation, d.LockedLocalMatrix(), X.LocalMatrix() );
        }
        else
        {
            DistMatrix<R,W,STAR> d_W_STAR( X.Grid() );
            d_W_STAR = d;
            DiagonalScale
            ( LEFT, orientation, 
              d_W_STAR.LockedLocalMatrix(), X.LocalMatrix() );
        }
    }
    else
    {
        if( U == Z && V == STAR && d.ColAlignment() == X.RowAlignment() )
        {
            DiagonalScale
            ( RIGHT, orientation, d.LockedLocalMatrix(), X.LocalMatrix() );
        }
        else
        {
            DistMatrix<R,Z,STAR> d_Z_STAR( X.Grid() );
            d_Z_STAR = d;
            DiagonalScale
            ( RIGHT, orientation, 
              d_Z_STAR.LockedLocalMatrix(), X.LocalMatrix() );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F,Distribution U,Distribution V,
                    Distribution W,Distribution Z>
inline void
DiagonalSolve
( LeftOrRight side, Orientation orientation, 
  const DistMatrix<F,U,V>& d, DistMatrix<F,W,Z>& X,
  bool checkIfSingular )
{
#ifndef RELEASE
    PushCallStack("DiagonalSolve");
#endif
    if( side == LEFT )
    {
        if( U == W && V == STAR && d.ColAlignment() == X.ColAlignment() )
        {
            DiagonalSolve
            ( LEFT, orientation, d.LockedLocalMatrix(), X.LocalMatrix(),
              checkIfSingular );
        }
        else
        {
            DistMatrix<F,W,STAR> d_W_STAR( X.Grid() );
            d_W_STAR = d;
            DiagonalSolve
            ( LEFT, orientation, 
              d_W_STAR.LockedLocalMatrix(), X.LocalMatrix(), checkIfSingular );
        }
    }
    else
    {
        if( U == Z && V == STAR && d.ColAlignment() == X.RowAlignment() )
        {
            DiagonalSolve
            ( RIGHT, orientation, d.LockedLocalMatrix(), X.LocalMatrix(),
              checkIfSingular );
        }
        else
        {
            DistMatrix<F,Z,STAR> d_Z_STAR( X.Grid() );
            d_Z_STAR = d;
            DiagonalSolve
            ( RIGHT, orientation, 
              d_Z_STAR.LockedLocalMatrix(), X.LocalMatrix(), checkIfSingular );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F,Distribution U,Distribution V,
                    Distribution W,Distribution Z>
inline void
DiagonalSolve
( LeftOrRight side, Orientation orientation, 
  const DistMatrix<typename Base<F>::type,U,V>& d, DistMatrix<F,W,Z>& X,
  bool checkIfSingular )
{
#ifndef RELEASE
    PushCallStack("DiagonalSolve");
#endif
    typedef typename Base<F>::type R;

    if( side == LEFT )
    {
        if( U == W && V == STAR && d.ColAlignment() == X.ColAlignment() )
        {
            DiagonalSolve
            ( LEFT, orientation, d.LockedLocalMatrix(), X.LocalMatrix(),
              checkIfSingular );
        }
        else
        {
            DistMatrix<R,W,STAR> d_W_STAR( X.Grid() );
            d_W_STAR = d;
            DiagonalSolve
            ( LEFT, orientation, 
              d_W_STAR.LockedLocalMatrix(), X.LocalMatrix(), checkIfSingular );
        }
    }
    else
    {
        if( U == Z && V == STAR && d.ColAlignment() == X.RowAlignment() )
        {
            DiagonalSolve
            ( RIGHT, orientation, d.LockedLocalMatrix(), X.LocalMatrix(),
              checkIfSingular );
        }
        else
        {
            DistMatrix<R,Z,STAR> d_Z_STAR( X.Grid() );
            d_Z_STAR = d;
            DiagonalSolve
            ( RIGHT, orientation, 
              d_Z_STAR.LockedLocalMatrix(), X.LocalMatrix(), checkIfSingular );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

// Our extended Dotc is equivalent to our extended Dot, 
// but we are burdened with consistency
template<typename T,Distribution U,Distribution V,
                    Distribution W,Distribution Z>
inline T
Dotc( const DistMatrix<T,U,V>& x, const DistMatrix<T,W,Z>& y )
{
#ifndef RELEASE
    PushCallStack("Dotc");
#endif
    T globalDot = Dot( x, y );
#ifndef RELEASE
    PopCallStack();
#endif
    return globalDot;
}

template<typename T,Distribution U,Distribution V>
inline void
Scale( T alpha, DistMatrix<T,U,V>& A )
{ Scale( alpha, A.LocalMatrix() ); }

template<typename T,Distribution U,Distribution V>
inline void
Scale( typename Base<T>::type alpha, DistMatrix<T,U,V>& A )
{ Scale( (T)alpha, A.LocalMatrix() ); }

template<typename T,Distribution U,Distribution V>
inline void
Scal( T alpha, DistMatrix<T,U,V>& A )
{ Scale( alpha, A ); }

template<typename T,Distribution U,Distribution V>
inline void
Scal( typename Base<T>::type alpha, DistMatrix<T,U,V>& A )
{ Scale( (T)alpha, A ); }

//----------------------------------------------------------------------------//
// Distributed BLAS-like routines: Level 1 (extensions)                       //
//----------------------------------------------------------------------------//

template<typename T,Distribution U,Distribution V,
                    Distribution W,Distribution Z>
inline void
Adjoint( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B )
{
#ifndef RELEASE
    PushCallStack("Adjoint");
#endif
    if( B.Viewing() )
    {
        if( A.Height() != B.Width() || A.Width() != B.Height() )
        {
            std::ostringstream msg;
            msg << "If Adjoint'ing into a view, it must be the right size:\n"
                << "  A ~ " << A.Height() << " x " << A.Width() << "\n"
                << "  B ~ " << B.Height() << " x " << B.Width();
            throw std::logic_error( msg.str().c_str() );
        }
    }

    if( U == Z && V == W && 
        A.ColAlignment() == B.RowAlignment() && 
        A.RowAlignment() == B.ColAlignment() )
    {
        Adjoint( A.LockedLocalMatrix(), B.LocalMatrix() );
    }
    else
    {
        DistMatrix<T,Z,W> C( B.Grid() );
        if( B.Viewing() || B.ConstrainedColAlignment() )
            C.AlignRowsWith( B );
        if( B.Viewing() || B.ConstrainedRowAlignment() )
            C.AlignColsWith( B );
        C = A;

        if( !B.Viewing() )
        {
            if( !B.ConstrainedColAlignment() )
                B.AlignColsWith( C );
            if( !B.ConstrainedRowAlignment() )
                B.AlignRowsWith( C );
            B.ResizeTo( A.Width(), A.Height() );
        }
        Adjoint( C.LockedLocalMatrix(), B.LocalMatrix() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V>
inline void
Conjugate( DistMatrix<T,U,V>& A )
{
#ifndef RELEASE
    PushCallStack("Conjugate (in-place)");
#endif
    Conjugate( A.LocalMatrix() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V,
                    Distribution W,Distribution Z>
inline void
Conjugate( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B )
{
#ifndef RELEASE
    PushCallStack("Conjugate");
#endif
    B = A;
    Conjugate( B ); 
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V>
inline void
MakeTrapezoidal
( LeftOrRight side, UpperOrLower uplo, int offset, 
  DistMatrix<T,U,V>& A )
{
#ifndef RELEASE
    PushCallStack("MakeTrapezoidal");
#endif
    const int height = A.Height();
    const int width = A.Width();
    const int localHeight = A.LocalHeight();
    const int localWidth = A.LocalWidth();
    const int colShift = A.ColShift();
    const int rowShift = A.RowShift();
    const int colStride = A.ColStride();
    const int rowStride = A.RowStride();

    T* localBuffer = A.LocalBuffer();
    const int ldim = A.LocalLDim();

    if( uplo == LOWER )
    {

#ifdef HAVE_OPENMP
        #pragma omp parallel for
#endif
        for( int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            const int j = rowShift + jLocal*rowStride;
            const int lastZeroRow =
                ( side==LEFT ? j-offset-1
                             : j-offset+height-width-1 );
            if( lastZeroRow >= 0 )
            {
                const int boundary = std::min( lastZeroRow+1, height );
                const int numZeroRows =
                    RawLocalLength( boundary, colShift, colStride );
                MemZero( &localBuffer[jLocal*ldim], numZeroRows );
            }
        }
    }
    else
    {
#ifdef HAVE_OPENMP
        #pragma omp parallel for
#endif
        for( int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            const int j = rowShift + jLocal*rowStride;
            const int firstZeroRow =
                ( side==LEFT ? std::max(j-offset+1,0)
                             : std::max(j-offset+height-width+1,0) );
            const int numNonzeroRows =
                RawLocalLength(firstZeroRow,colShift,colStride);
            if( numNonzeroRows < localHeight )
            {
                T* col = &localBuffer[numNonzeroRows+jLocal*ldim];
                MemZero( col, localHeight-numNonzeroRows );
            }
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V>
inline void
ScaleTrapezoid
( T alpha, LeftOrRight side, UpperOrLower uplo, int offset, 
  DistMatrix<T,U,V>& A )
{
#ifndef RELEASE
    PushCallStack("ScaleTrapezoid");
#endif
    const int height = A.Height();
    const int width = A.Width();
    const int localHeight = A.LocalHeight();
    const int localWidth = A.LocalWidth();
    const int colShift = A.ColShift();
    const int rowShift = A.RowShift();
    const int colStride = A.ColStride();
    const int rowStride = A.RowStride();

    if( uplo == UPPER )
    {
        T* localBuffer = A.LocalBuffer();
        const int ldim = A.LocalLDim();
#ifdef HAVE_OPENMP
        #pragma omp parallel for
#endif
        for( int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            int j = rowShift + jLocal*rowStride;
            int lastRow = ( side==LEFT ? j-offset : j-offset+height-width );
            int boundary = std::min( lastRow+1, height );
            int numRows = RawLocalLength( boundary, colShift, colStride );
            T* col = &localBuffer[jLocal*ldim];
            for( int iLocal=0; iLocal<numRows; ++iLocal )
                col[iLocal] *= alpha;
        }
    }
    else
    {
        T* localBuffer = A.LocalBuffer();
        const int ldim = A.LocalLDim();
#ifdef HAVE_OPENMP
        #pragma omp parallel for
#endif
        for( int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            int j = rowShift + jLocal*rowStride;
            int firstRow =
                ( side==LEFT ? std::max(j-offset,0)
                             : std::max(j-offset+height-width,0) );
            int numZeroRows = RawLocalLength( firstRow, colShift, colStride );
            T* col = &localBuffer[numZeroRows+jLocal*ldim];
            for( int iLocal=0; iLocal<(localHeight-numZeroRows); ++iLocal )
                col[iLocal] *= alpha;
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V,
                    Distribution W,Distribution Z>
inline void
Transpose( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B )
{
#ifndef RELEASE
    PushCallStack("Transpose");
#endif
    if( B.Viewing() )
    {
        if( A.Height() != B.Width() || A.Width() != B.Height() )
        {
            std::ostringstream msg;
            msg << "If Transpose'ing into a view, it must be the right size:\n"
                << "  A ~ " << A.Height() << " x " << A.Width() << "\n"
                << "  B ~ " << B.Height() << " x " << B.Width();
            throw std::logic_error( msg.str().c_str() );
        }
    }

    if( U == Z && V == W && 
        A.ColAlignment() == B.RowAlignment() && 
        A.RowAlignment() == B.ColAlignment() )
    {
        Transpose( A.LockedLocalMatrix(), B.LocalMatrix() );
    }
    else
    {
        DistMatrix<T,Z,W> C( B.Grid() );
        if( B.Viewing() || B.ConstrainedColAlignment() )
            C.AlignRowsWith( B );
        if( B.Viewing() || B.ConstrainedRowAlignment() )
            C.AlignColsWith( B );
        C = A;

        if( !B.Viewing() )
        {
            if( !B.ConstrainedColAlignment() )
                B.AlignColsWith( C );
            if( !B.ConstrainedRowAlignment() )
                B.AlignRowsWith( C );
            B.ResizeTo( A.Width(), A.Height() );
        }
        Transpose( C.LockedLocalMatrix(), B.LocalMatrix() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V>
inline void
Zero( DistMatrix<T,U,V>& A )
{
#ifndef RELEASE
    PushCallStack("Zero");
#endif
    Zero( A.LocalMatrix() );
#ifndef RELEASE
    PopCallStack();
#endif
}

//----------------------------------------------------------------------------//
// Distributed BLAS-like routines: Level 2                                    //
//----------------------------------------------------------------------------//

// Our extended Ger and Gerc are equivalent
template<typename T>
inline void
Gerc
( T alpha, const DistMatrix<T>& x,
           const DistMatrix<T>& y,
                 DistMatrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("Gerc");
#endif
    Ger( alpha, x, y, A );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
