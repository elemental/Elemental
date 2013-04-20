/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_SKEWHERMITIANEIG_HPP
#define LAPACK_SKEWHERMITIANEIG_HPP

#include "elemental/blas-like/level1/ScaleTrapezoid.hpp"

namespace elem {

//----------------------------------------------------------------------------//
// Grab the full set of eigenvalues                                           //
//----------------------------------------------------------------------------//

template<typename R>
inline void
SkewHermitianEig( UpperOrLower uplo, Matrix<R>& G, Matrix<R>& wImag )
{
#ifndef RELEASE
    PushCallStack("SkewHermitianEig");
#endif
    if( G.Height() != G.Width() )
        throw std::logic_error("SkewHermitian matrices must be square");
    int n = G.Height();
    Matrix<Complex<R> > A( n, n );
    const int ALDim = A.LDim();
    const int GLDim = G.LDim();    
    const Complex<R> negativeImagOne(0,-1);
    const R* GBuffer = G.Buffer();
    Complex<R>* ABuffer = A.Buffer();
    if( uplo == LOWER )
        for( int j=0; j<n; ++j )
            for( int i=j; i<n; ++i )
                ABuffer[i+j*ALDim] = negativeImagOne*GBuffer[i+j*GLDim];
    else
        for( int j=0; j<n; ++j )
            for( int i=0; i<=j; ++i )
                ABuffer[i+j*ALDim] = negativeImagOne*GBuffer[i+j*GLDim];

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, A, wImag );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline void
SkewHermitianEig( UpperOrLower uplo, Matrix<Complex<R> >& G, Matrix<R>& wImag )
{
#ifndef RELEASE
    PushCallStack("SkewHermitianEig");
#endif
    if( G.Height() != G.Width() )
        throw std::logic_error("Skew-Hermitian matrices must be square");
    
    // Make G Hermitian by scaling by -i
    const Complex<R> negativeImagOne(0,-1);
    ScaleTrapezoid( negativeImagOne, LEFT, uplo, 0, G );

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, G, wImag );
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifdef HAVE_PMRRR
template<typename R>
inline void
SkewHermitianEig
( UpperOrLower uplo, 
  DistMatrix<R>& G,
  DistMatrix<R,VR,STAR>& wImag )
{
#ifndef RELEASE
    PushCallStack("SkewHermitianEig");
#endif
    if( G.Height() != G.Width() )
        throw std::logic_error("SkewHermitian matrices must be square");

    int n = G.Height();
    const Grid& grid = G.Grid();

    DistMatrix<Complex<R> > A(grid);
    A.AlignWith( G.DistData() );
    A.ResizeTo( n, n );

    const int localHeight = A.LocalHeight();
    const int localWidth = A.LocalWidth();
    const int ALDim = A.LDim();
    const int GLDim = G.LDim();    
    const Complex<R> negativeImagOne(0,-1);
    const R* GBuffer = G.Buffer();
    Complex<R>* ABuffer = A.Buffer();
    // Just copy the entire local matrix instead of worrying about symmetry
    for( int j=0; j<localWidth; ++j )
        for( int i=0; i<localHeight; ++i )
            ABuffer[i+j*ALDim] = negativeImagOne*GBuffer[i+j*GLDim];

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, A, wImag );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline void
SkewHermitianEig
( UpperOrLower uplo, 
  DistMatrix<Complex<R> >& G,
  DistMatrix<R,VR,STAR>& wImag )
{
#ifndef RELEASE
    PushCallStack("SkewHermitianEig");
#endif
    if( G.Height() != G.Width() )
        throw std::logic_error("Skew-Hermitian matrices must be square");
    
    // Make G Hermitian by scaling by -i
    const Complex<R> negativeImagOne(0,-1);
    ScaleTrapezoid( negativeImagOne, LEFT, uplo, 0, G );

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, G, wImag );
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // ifdef HAVE_PMRRR

//----------------------------------------------------------------------------//
// Grab the full set of eigenpairs                                            //
//----------------------------------------------------------------------------//

template<typename R>
inline void
SkewHermitianEig
( UpperOrLower uplo, Matrix<R>& G, Matrix<R>& wImag, Matrix<Complex<R> >& Z )
{
#ifndef RELEASE
    PushCallStack("SkewHermitianEig");
#endif
    if( G.Height() != G.Width() )
        throw std::logic_error("SkewHermitian matrices must be square");
    int n = G.Height();
    Matrix<Complex<R> > A( n, n );

    const int ALDim = A.LDim();
    const int GLDim = G.LDim();
    const Complex<R> negativeImagOne(0,-1);
    const R* GBuffer = G.Buffer();
    Complex<R>* ABuffer = A.Buffer();
    if( uplo == LOWER )
        for( int j=0; j<n; ++j )
            for( int i=j; i<n; ++i )
                ABuffer[i+j*ALDim] = negativeImagOne*GBuffer[i+j*GLDim];
    else
        for( int j=0; j<n; ++j )
            for( int i=0; i<=j; ++i )
                ABuffer[i+j*ALDim] = negativeImagOne*GBuffer[i+j*GLDim];

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, A, wImag, Z );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline void
SkewHermitianEig
( UpperOrLower uplo, 
  Matrix<Complex<R> >& G, Matrix<R>& wImag, Matrix<Complex<R> >& Z )
{
#ifndef RELEASE
    PushCallStack("SkewHermitianEig");
#endif
    if( G.Height() != G.Width() )
        throw std::logic_error("Skew-Hermitian matrices must be square");

    // Make G Hermitian by scaling by -i
    const Complex<R> negativeImagOne(0,-1);
    ScaleTrapezoid( negativeImagOne, LEFT, uplo, 0, G );

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, G, wImag, Z );
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifdef HAVE_PMRRR
template<typename R>
inline void
SkewHermitianEig
( UpperOrLower uplo, 
  DistMatrix<R>& G,
  DistMatrix<R,VR,STAR>& wImag,
  DistMatrix<Complex<R> >& Z )
{
#ifndef RELEASE
    PushCallStack("SkewHermitianEig");
#endif
    if( G.Height() != G.Width() )
        throw std::logic_error("SkewHermitian matrices must be square");

    int n = G.Height();
    const Grid& grid = G.Grid();

    DistMatrix<Complex<R> > A(grid);
    A.AlignWith( G.DistData() );
    A.ResizeTo( n, n );

    const int localHeight = A.LocalHeight();
    const int localWidth = A.LocalWidth();
    const int ALDim = A.LDim();
    const int GLDim = G.LDim();
    const Complex<R> negativeImagOne(0,-1);
    const R* GBuffer = G.Buffer();
    Complex<R>* ABuffer = A.Buffer();
    // Just copy the entire local matrix instead of worrying about symmetry
    for( int j=0; j<localWidth; ++j )
        for( int i=0; i<localHeight; ++i )
            ABuffer[i+j*ALDim] = negativeImagOne*GBuffer[i+j*GLDim];

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, A, wImag, Z );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline void
SkewHermitianEig
( UpperOrLower uplo, 
  DistMatrix<Complex<R> >& G,
  DistMatrix<R,VR,STAR>& wImag,
  DistMatrix<Complex<R> >& Z )
{
#ifndef RELEASE
    PushCallStack("SkewHermitianEig");
#endif
    if( G.Height() != G.Width() )
        throw std::logic_error("Skew-Hermitian matrices must be square");

    // Make G Hermitian by scaling by -i
    const Complex<R> negativeImagOne(0,-1);
    ScaleTrapezoid( negativeImagOne, LEFT, uplo, 0, G );

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, G, wImag, Z );
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // ifdef HAVE_PMRRR

//----------------------------------------------------------------------------//
// Grab the eigenvalues with indices in a specified range                     //
//----------------------------------------------------------------------------//

template<typename R>
inline void
SkewHermitianEig
( UpperOrLower uplo, Matrix<R>& G, Matrix<R>& wImag, int a, int b )
{
#ifndef RELEASE
    PushCallStack("SkewHermitianEig");
#endif
    if( G.Height() != G.Width() )
        throw std::logic_error("Skew-Hermitian matrices must be square");
    int n = G.Height();

    DistMatrix<Complex<R> > A( n, n );
    const int ALDim = A.LDim();
    const int GLDim = G.LDim();    
    const Complex<R> negativeImagOne(0,-1);
    const R* GBuffer = G.Buffer();
    Complex<R>* ABuffer = A.Buffer();
    if( uplo == LOWER )
        for( int j=0; j<n; ++j )
            for( int i=j; i<n; ++i )
                ABuffer[i+j*ALDim] = negativeImagOne*GBuffer[i+j*GLDim];
    else
        for( int j=0; j<n; ++j )
            for( int i=0; i<=j; ++i )
                ABuffer[i+j*ALDim] = negativeImagOne*GBuffer[i+j*GLDim];

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, A, wImag, a, b );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline void
SkewHermitianEig
( UpperOrLower uplo, Matrix<Complex<R> >& G, Matrix<R>& wImag, int a, int b )
{
#ifndef RELEASE
    PushCallStack("SkewHermitianEig");
#endif
    if( G.Height() != G.Width() )
        throw std::logic_error("Skew-Hermitian matrices must be square");
    
    // Make G Hermitian by scaling by -i
    const Complex<R> negativeImagOne(0,-1);
    ScaleTrapezoid( negativeImagOne, LEFT, uplo, 0, G );

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, G, wImag, a, b );
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifdef HAVE_PMRRR
template<typename R>
inline void
SkewHermitianEig
( UpperOrLower uplo, 
  DistMatrix<R>& G,
  DistMatrix<R,VR,STAR>& wImag,
  int a, int b )
{
#ifndef RELEASE
    PushCallStack("SkewHermitianEig");
#endif
    if( G.Height() != G.Width() )
        throw std::logic_error("Skew-Hermitian matrices must be square");

    int n = G.Height();
    const Grid& grid = G.Grid();

    DistMatrix<Complex<R> > A(grid);
    A.AlignWith( G.DistData() );
    A.ResizeTo( n, n );

    const int localHeight = A.LocalHeight();
    const int localWidth = A.LocalWidth();
    const int ALDim = A.LDim();
    const int GLDim = G.LDim();    
    const Complex<R> negativeImagOne(0,-1);
    const R* GBuffer = G.Buffer();
    Complex<R>* ABuffer = A.Buffer();
    // Just copy the entire local matrix instead of worrying about symmetry
    for( int j=0; j<localWidth; ++j )
        for( int i=0; i<localHeight; ++i )
            ABuffer[i+j*ALDim] = negativeImagOne*GBuffer[i+j*GLDim];

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, A, wImag, a, b );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline void
SkewHermitianEig
( UpperOrLower uplo, 
  DistMatrix<Complex<R> >& G,
  DistMatrix<R,VR,STAR>& wImag,
  int a, int b )
{
#ifndef RELEASE
    PushCallStack("SkewHermitianEig");
#endif
    if( G.Height() != G.Width() )
        throw std::logic_error("Skew-Hermitian matrices must be square");
    
    // Make G Hermitian by scaling by -i
    const Complex<R> negativeImagOne(0,-1);
    ScaleTrapezoid( negativeImagOne, LEFT, uplo, 0, G );

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, G, wImag, a, b );
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // ifdef HAVE_PMRRR

//----------------------------------------------------------------------------//
// Grab the eigenpairs with indices in a specified range                      //
//----------------------------------------------------------------------------//

template<typename R>
inline void
SkewHermitianEig
( UpperOrLower uplo, Matrix<R>& G, Matrix<R>& wImag, Matrix<Complex<R> >& Z,
  int a, int b )
{
#ifndef RELEASE
    PushCallStack("SkewHermitianEig");
#endif
    if( G.Height() != G.Width() )
        throw std::logic_error("Skew-Hermitian matrices must be square");
    int n = G.Height();

    Matrix<Complex<R> > A( n, n );
    const int ALDim = A.LDim();
    const int GLDim = G.LDim();    
    const Complex<R> negativeImagOne(0,-1);
    const R* GBuffer = G.Buffer();
    Complex<R>* ABuffer = A.Buffer();
    if( uplo == LOWER )
        for( int j=0; j<n; ++j )
            for( int i=j; i<n; ++i )
                ABuffer[i+j*ALDim] = negativeImagOne*GBuffer[i+j*GLDim];
    else
        for( int j=0; j<n; ++j )
            for( int i=0; i<=j; ++i )
                ABuffer[i+j*ALDim] = negativeImagOne*GBuffer[i+j*GLDim];

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, A, wImag, Z, a, b );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline void
SkewHermitianEig
( UpperOrLower uplo, 
  Matrix<Complex<R> >& G, Matrix<R>& wImag, Matrix<Complex<R> >& Z,
  int a, int b )
{
#ifndef RELEASE
    PushCallStack("SkewHermitianEig");
#endif
    if( G.Height() != G.Width() )
        throw std::logic_error("Skew-Hermitian matrices must be square");
    
    // Make G Hermitian by scaling by -i
    const Complex<R> negativeImagOne(0,-1);
    ScaleTrapezoid( negativeImagOne, LEFT, uplo, 0, G );

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, G, wImag, Z, a, b );
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifdef HAVE_PMRRR
template<typename R>
inline void
SkewHermitianEig
( UpperOrLower uplo, 
  DistMatrix<R>& G,
  DistMatrix<R,VR,STAR>& wImag,
  DistMatrix<Complex<R> >& Z,
  int a, int b )
{
#ifndef RELEASE
    PushCallStack("SkewHermitianEig");
#endif
    if( G.Height() != G.Width() )
        throw std::logic_error("Skew-Hermitian matrices must be square");

    int n = G.Height();
    const Grid& grid = G.Grid();

    DistMatrix<Complex<R> > A(grid);
    A.AlignWith( G.DistData() );
    A.ResizeTo( n, n );

    const int localHeight = A.LocalHeight();
    const int localWidth = A.LocalWidth();
    const int ALDim = A.LDim();
    const int GLDim = G.LDim();    
    const Complex<R> negativeImagOne(0,-1);
    const R* GBuffer = G.Buffer();
    Complex<R>* ABuffer = A.Buffer();
    // Just copy the entire local matrix instead of worrying about symmetry
    for( int j=0; j<localWidth; ++j )
        for( int i=0; i<localHeight; ++i )
            ABuffer[i+j*ALDim] = negativeImagOne*GBuffer[i+j*GLDim];

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, A, wImag, Z, a, b );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline void
SkewHermitianEig
( UpperOrLower uplo, 
  DistMatrix<Complex<R> >& G,
  DistMatrix<R,VR,STAR>& wImag,
  DistMatrix<Complex<R> >& Z,
  int a, int b )
{
#ifndef RELEASE
    PushCallStack("SkewHermitianEig");
#endif
    if( G.Height() != G.Width() )
        throw std::logic_error("Skew-Hermitian matrices must be square");
    
    // Make G Hermitian by scaling by -i
    const Complex<R> negativeImagOne(0,-1);
    ScaleTrapezoid( negativeImagOne, LEFT, uplo, 0, G );

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, G, wImag, Z, a, b );
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // ifdef HAVE_PMRRR

//----------------------------------------------------------------------------//
// Grab the eigenvalues in the interval i(a,b]                                //
//----------------------------------------------------------------------------//

template<typename R>
inline void
SkewHermitianEig( UpperOrLower uplo, Matrix<R>& G, Matrix<R>& wImag, R a, R b )
{
#ifndef RELEASE
    PushCallStack("SkewHermitianEig");
#endif
    if( G.Height() != G.Width() )
        throw std::logic_error("Skew-Hermitian matrices must be square");
    int n = G.Height();

    Matrix<Complex<R> > A( n, n );
    const int ALDim = A.LDim();
    const int GLDim = G.LDim();    
    const Complex<R> negativeImagOne(0,-1);
    const R* GBuffer = G.Buffer();
    Complex<R>* ABuffer = A.Buffer();
    if( uplo == LOWER )
        for( int j=0; j<n; ++j )
            for( int i=j; i<n; ++i )
                ABuffer[i+j*ALDim] = negativeImagOne*GBuffer[i+j*GLDim];
    else
        for( int j=0; j<n; ++j )
            for( int i=0; i<=j; ++i )
                ABuffer[i+j*ALDim] = negativeImagOne*GBuffer[i+j*GLDim];

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, A, wImag, a, b );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline void
SkewHermitianEig
( UpperOrLower uplo, Matrix<Complex<R> >& G, Matrix<R>& wImag, R a, R b )
{
#ifndef RELEASE
    PushCallStack("SkewHermitianEig");
#endif
    if( G.Height() != G.Width() )
        throw std::logic_error("Skew-Hermitian matrices must be square");
    
    // Make G Hermitian by scaling by -i
    const Complex<R> negativeImagOne(0,-1);
    ScaleTrapezoid( negativeImagOne, LEFT, uplo, 0, G );

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, G, wImag, a, b );
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifdef HAVE_PMRRR
template<typename R>
inline void
SkewHermitianEig
( UpperOrLower uplo, 
  DistMatrix<R>& G,
  DistMatrix<R,VR,STAR>& wImag,
  R a, R b )
{
#ifndef RELEASE
    PushCallStack("SkewHermitianEig");
#endif
    if( G.Height() != G.Width() )
        throw std::logic_error("Skew-Hermitian matrices must be square");

    int n = G.Height();
    const Grid& grid = G.Grid();

    DistMatrix<Complex<R> > A(grid);
    A.AlignWith( G.DistData() );
    A.ResizeTo( n, n );

    const int localHeight = A.LocalHeight();
    const int localWidth = A.LocalWidth();
    const int ALDim = A.LDim();
    const int GLDim = G.LDim();    
    const Complex<R> negativeImagOne(0,-1);
    const R* GBuffer = G.Buffer();
    Complex<R>* ABuffer = A.Buffer();
    // Just copy the entire local matrix instead of worrying about symmetry
    for( int j=0; j<localWidth; ++j )
        for( int i=0; i<localHeight; ++i )
            ABuffer[i+j*ALDim] = negativeImagOne*GBuffer[i+j*GLDim];

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, A, wImag, a, b );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline void
SkewHermitianEig
( UpperOrLower uplo, 
  DistMatrix<Complex<R> >& G,
  DistMatrix<R,VR,STAR>& wImag,
  R a, R b )
{
#ifndef RELEASE
    PushCallStack("SkewHermitianEig");
#endif
    if( G.Height() != G.Width() )
        throw std::logic_error("Skew-Hermitian matrices must be square");
    
    // Make G Hermitian by scaling by -i
    const Complex<R> negativeImagOne(0,-1);
    ScaleTrapezoid( negativeImagOne, LEFT, uplo, 0, G );

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, G, wImag, a, b );
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // ifdef HAVE_PMRRR

//----------------------------------------------------------------------------//
// Grab the eigenpairs with eigenvalues in the interval i(a,b]                //
//----------------------------------------------------------------------------//

template<typename R>
inline void
SkewHermitianEig
( UpperOrLower uplo, Matrix<R>& G, Matrix<R>& wImag, Matrix<Complex<R> >& Z,
  R a, R b )
{
#ifndef RELEASE
    PushCallStack("SkewHermitianEig");
#endif
    if( G.Height() != G.Width() )
        throw std::logic_error("SkewHermitian matrices must be square");
    int n = G.Height();

    Matrix<Complex<R> > A( n, n );
    const int ALDim = A.LDim();
    const int GLDim = G.LDim();    
    const Complex<R> negativeImagOne(0,-1);
    const R* GBuffer = G.Buffer();
    Complex<R>* ABuffer = A.Buffer();
    if( uplo == LOWER )
        for( int j=0; j<n; ++j )
            for( int i=j; i<n; ++i )
                ABuffer[i+j*ALDim] = negativeImagOne*GBuffer[i+j*GLDim];
    else
        for( int j=0; j<n; ++j )
            for( int i=0; i<=j; ++i )
                ABuffer[i+j*ALDim] = negativeImagOne*GBuffer[i+j*GLDim];

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, A, wImag, Z, a, b );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline void
SkewHermitianEig
( UpperOrLower uplo, 
  Matrix<Complex<R> >& G, Matrix<R>& wImag, Matrix<Complex<R> >& Z,
  R a, R b )
{
#ifndef RELEASE
    PushCallStack("SkewHermitianEig");
#endif
    if( G.Height() != G.Width() )
        throw std::logic_error("Skew-Hermitian matrices must be square");
    
    // Make G Hermitian by scaling by -i
    const Complex<R> negativeImagOne(0,-1);
    ScaleTrapezoid( negativeImagOne, LEFT, uplo, 0, G );

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, G, wImag, Z, a, b );
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifdef HAVE_PMRRR
template<typename R>
inline void
SkewHermitianEig
( UpperOrLower uplo, 
  DistMatrix<R>& G,
  DistMatrix<R,VR,STAR>& wImag,
  DistMatrix<Complex<R> >& Z,
  R a, R b )
{
#ifndef RELEASE
    PushCallStack("SkewHermitianEig");
#endif
    if( G.Height() != G.Width() )
        throw std::logic_error("SkewHermitian matrices must be square");

    int n = G.Height();
    const Grid& grid = G.Grid();

    DistMatrix<Complex<R> > A(grid);
    A.AlignWith( G.DistData() );
    A.ResizeTo( n, n );

    const int localHeight = A.LocalHeight();
    const int localWidth = A.LocalWidth();
    const int ALDim = A.LDim();
    const int GLDim = G.LDim();    
    const Complex<R> negativeImagOne(0,-1);
    const R* GBuffer = G.Buffer();
    Complex<R>* ABuffer = A.Buffer();
    // Just copy the entire local matrix instead of worrying about symmetry
    for( int j=0; j<localWidth; ++j )
        for( int i=0; i<localHeight; ++i )
            ABuffer[i+j*ALDim] = negativeImagOne*GBuffer[i+j*GLDim];

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, A, wImag, Z, a, b );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline void
SkewHermitianEig
( UpperOrLower uplo, 
  DistMatrix<Complex<R> >& G,
  DistMatrix<R,VR,STAR>& wImag,
  DistMatrix<Complex<R> >& Z,
  R a, R b )
{
#ifndef RELEASE
    PushCallStack("SkewHermitianEig");
#endif
    if( G.Height() != G.Width() )
        throw std::logic_error("Skew-Hermitian matrices must be square");
    
    // Make G Hermitian by scaling by -i
    const Complex<R> negativeImagOne(0,-1);
    ScaleTrapezoid( negativeImagOne, LEFT, uplo, 0, G );

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, G, wImag, Z, a, b );
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // ifdef HAVE_PMRRR

} // namespace elem

#endif // ifndef LAPACK_SKEWHERMITIANEIG_HPP
