/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_SKEWHERMITIANEIG_HPP
#define ELEM_LAPACK_SKEWHERMITIANEIG_HPP

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
    CallStackEntry entry("SkewHermitianEig");
#endif
    if( G.Height() != G.Width() )
        LogicError("SkewHermitian matrices must be square");
    Int n = G.Height();
    Matrix<Complex<R> > A( n, n );
    const Int ALDim = A.LDim();
    const Int GLDim = G.LDim();    
    const Complex<R> negativeImagOne(0,-1);
    const R* GBuffer = G.Buffer();
    Complex<R>* ABuffer = A.Buffer();
    if( uplo == LOWER )
        for( Int j=0; j<n; ++j )
            for( Int i=j; i<n; ++i )
                ABuffer[i+j*ALDim] = negativeImagOne*GBuffer[i+j*GLDim];
    else
        for( Int j=0; j<n; ++j )
            for( Int i=0; i<=j; ++i )
                ABuffer[i+j*ALDim] = negativeImagOne*GBuffer[i+j*GLDim];

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, A, wImag );
}

template<typename R>
inline void
SkewHermitianEig( UpperOrLower uplo, Matrix<Complex<R> >& G, Matrix<R>& wImag )
{
#ifndef RELEASE
    CallStackEntry entry("SkewHermitianEig");
#endif
    if( G.Height() != G.Width() )
        LogicError("Skew-Hermitian matrices must be square");
    
    // Make G Hermitian by scaling by -i
    const Complex<R> negativeImagOne(0,-1);
    ScaleTrapezoid( negativeImagOne, uplo, G );

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, G, wImag );
}

template<typename R>
inline void
SkewHermitianEig
( UpperOrLower uplo, 
  DistMatrix<R>& G,
  DistMatrix<R,VR,STAR>& wImag )
{
#ifndef RELEASE
    CallStackEntry entry("SkewHermitianEig");
#endif
    EnsurePMRRR();
    if( G.Height() != G.Width() )
        LogicError("SkewHermitian matrices must be square");

    Int n = G.Height();
    const Grid& grid = G.Grid();

    DistMatrix<Complex<R> > A(grid);
    A.AlignWith( G.DistData() );
    A.ResizeTo( n, n );

    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();
    const Int ALDim = A.LDim();
    const Int GLDim = G.LDim();    
    const Complex<R> negativeImagOne(0,-1);
    const R* GBuffer = G.Buffer();
    Complex<R>* ABuffer = A.Buffer();
    // Just copy the entire local matrix instead of worrying about symmetry
    for( Int j=0; j<localWidth; ++j )
        for( Int i=0; i<localHeight; ++i )
            ABuffer[i+j*ALDim] = negativeImagOne*GBuffer[i+j*GLDim];

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, A, wImag );
}

template<typename R>
inline void
SkewHermitianEig
( UpperOrLower uplo, 
  DistMatrix<Complex<R> >& G,
  DistMatrix<R,VR,STAR>& wImag )
{
#ifndef RELEASE
    CallStackEntry entry("SkewHermitianEig");
#endif
    if( G.Height() != G.Width() )
        LogicError("Skew-Hermitian matrices must be square");
    
    // Make G Hermitian by scaling by -i
    const Complex<R> negativeImagOne(0,-1);
    ScaleTrapezoid( negativeImagOne, uplo, G );

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, G, wImag );
}

//----------------------------------------------------------------------------//
// Grab the full set of eigenpairs                                            //
//----------------------------------------------------------------------------//

template<typename R>
inline void
SkewHermitianEig
( UpperOrLower uplo, Matrix<R>& G, Matrix<R>& wImag, Matrix<Complex<R> >& Z )
{
#ifndef RELEASE
    CallStackEntry entry("SkewHermitianEig");
#endif
    if( G.Height() != G.Width() )
        LogicError("SkewHermitian matrices must be square");
    Int n = G.Height();
    Matrix<Complex<R> > A( n, n );

    const Int ALDim = A.LDim();
    const Int GLDim = G.LDim();
    const Complex<R> negativeImagOne(0,-1);
    const R* GBuffer = G.Buffer();
    Complex<R>* ABuffer = A.Buffer();
    if( uplo == LOWER )
        for( Int j=0; j<n; ++j )
            for( Int i=j; i<n; ++i )
                ABuffer[i+j*ALDim] = negativeImagOne*GBuffer[i+j*GLDim];
    else
        for( Int j=0; j<n; ++j )
            for( Int i=0; i<=j; ++i )
                ABuffer[i+j*ALDim] = negativeImagOne*GBuffer[i+j*GLDim];

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, A, wImag, Z );
}

template<typename R>
inline void
SkewHermitianEig
( UpperOrLower uplo, 
  Matrix<Complex<R> >& G, Matrix<R>& wImag, Matrix<Complex<R> >& Z )
{
#ifndef RELEASE
    CallStackEntry entry("SkewHermitianEig");
#endif
    if( G.Height() != G.Width() )
        LogicError("Skew-Hermitian matrices must be square");

    // Make G Hermitian by scaling by -i
    const Complex<R> negativeImagOne(0,-1);
    ScaleTrapezoid( negativeImagOne, uplo, G );

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, G, wImag, Z );
}

template<typename R>
inline void
SkewHermitianEig
( UpperOrLower uplo, 
  DistMatrix<R>& G,
  DistMatrix<R,VR,STAR>& wImag,
  DistMatrix<Complex<R> >& Z )
{
#ifndef RELEASE
    CallStackEntry entry("SkewHermitianEig");
#endif
    EnsurePMRRR();
    if( G.Height() != G.Width() )
        LogicError("SkewHermitian matrices must be square");

    Int n = G.Height();
    const Grid& grid = G.Grid();

    DistMatrix<Complex<R> > A(grid);
    A.AlignWith( G.DistData() );
    A.ResizeTo( n, n );

    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();
    const Int ALDim = A.LDim();
    const Int GLDim = G.LDim();
    const Complex<R> negativeImagOne(0,-1);
    const R* GBuffer = G.Buffer();
    Complex<R>* ABuffer = A.Buffer();
    // Just copy the entire local matrix instead of worrying about symmetry
    for( Int j=0; j<localWidth; ++j )
        for( Int i=0; i<localHeight; ++i )
            ABuffer[i+j*ALDim] = negativeImagOne*GBuffer[i+j*GLDim];

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, A, wImag, Z );
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
    CallStackEntry entry("SkewHermitianEig");
#endif
    EnsurePMRRR();
    if( G.Height() != G.Width() )
        LogicError("Skew-Hermitian matrices must be square");

    // Make G Hermitian by scaling by -i
    const Complex<R> negativeImagOne(0,-1);
    ScaleTrapezoid( negativeImagOne, uplo, G );

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, G, wImag, Z );
}

//----------------------------------------------------------------------------//
// Grab the eigenvalues with indices in a specified range                     //
//----------------------------------------------------------------------------//

template<typename R>
inline void
SkewHermitianEig
( UpperOrLower uplo, Matrix<R>& G, Matrix<R>& wImag, Int a, Int b )
{
#ifndef RELEASE
    CallStackEntry entry("SkewHermitianEig");
#endif
    if( G.Height() != G.Width() )
        LogicError("Skew-Hermitian matrices must be square");
    Int n = G.Height();

    Matrix<Complex<R> > A( n, n );
    const Int ALDim = A.LDim();
    const Int GLDim = G.LDim();    
    const Complex<R> negativeImagOne(0,-1);
    const R* GBuffer = G.Buffer();
    Complex<R>* ABuffer = A.Buffer();
    if( uplo == LOWER )
        for( Int j=0; j<n; ++j )
            for( Int i=j; i<n; ++i )
                ABuffer[i+j*ALDim] = negativeImagOne*GBuffer[i+j*GLDim];
    else
        for( Int j=0; j<n; ++j )
            for( Int i=0; i<=j; ++i )
                ABuffer[i+j*ALDim] = negativeImagOne*GBuffer[i+j*GLDim];

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, A, wImag, a, b );
}

template<typename R>
inline void
SkewHermitianEig
( UpperOrLower uplo, Matrix<Complex<R> >& G, Matrix<R>& wImag, Int a, Int b )
{
#ifndef RELEASE
    CallStackEntry entry("SkewHermitianEig");
#endif
    if( G.Height() != G.Width() )
        LogicError("Skew-Hermitian matrices must be square");
    
    // Make G Hermitian by scaling by -i
    const Complex<R> negativeImagOne(0,-1);
    ScaleTrapezoid( negativeImagOne, uplo, G );

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, G, wImag, a, b );
}

template<typename R>
inline void
SkewHermitianEig
( UpperOrLower uplo, 
  DistMatrix<R>& G,
  DistMatrix<R,VR,STAR>& wImag,
  Int a, Int b )
{
#ifndef RELEASE
    CallStackEntry entry("SkewHermitianEig");
#endif
    EnsurePMRRR();
    if( G.Height() != G.Width() )
        LogicError("Skew-Hermitian matrices must be square");

    Int n = G.Height();
    const Grid& grid = G.Grid();

    DistMatrix<Complex<R> > A(grid);
    A.AlignWith( G.DistData() );
    A.ResizeTo( n, n );

    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();
    const Int ALDim = A.LDim();
    const Int GLDim = G.LDim();    
    const Complex<R> negativeImagOne(0,-1);
    const R* GBuffer = G.Buffer();
    Complex<R>* ABuffer = A.Buffer();
    // Just copy the entire local matrix instead of worrying about symmetry
    for( Int j=0; j<localWidth; ++j )
        for( Int i=0; i<localHeight; ++i )
            ABuffer[i+j*ALDim] = negativeImagOne*GBuffer[i+j*GLDim];

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, A, wImag, a, b );
}

template<typename R>
inline void
SkewHermitianEig
( UpperOrLower uplo, 
  DistMatrix<Complex<R> >& G,
  DistMatrix<R,VR,STAR>& wImag,
  Int a, Int b )
{
#ifndef RELEASE
    CallStackEntry entry("SkewHermitianEig");
#endif
    EnsurePMRRR();
    if( G.Height() != G.Width() )
        LogicError("Skew-Hermitian matrices must be square");
    
    // Make G Hermitian by scaling by -i
    const Complex<R> negativeImagOne(0,-1);
    ScaleTrapezoid( negativeImagOne, uplo, G );

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, G, wImag, a, b );
}

//----------------------------------------------------------------------------//
// Grab the eigenpairs with indices in a specified range                      //
//----------------------------------------------------------------------------//

template<typename R>
inline void
SkewHermitianEig
( UpperOrLower uplo, Matrix<R>& G, Matrix<R>& wImag, Matrix<Complex<R> >& Z,
  Int a, Int b )
{
#ifndef RELEASE
    CallStackEntry entry("SkewHermitianEig");
#endif
    if( G.Height() != G.Width() )
        LogicError("Skew-Hermitian matrices must be square");
    Int n = G.Height();

    Matrix<Complex<R> > A( n, n );
    const Int ALDim = A.LDim();
    const Int GLDim = G.LDim();    
    const Complex<R> negativeImagOne(0,-1);
    const R* GBuffer = G.Buffer();
    Complex<R>* ABuffer = A.Buffer();
    if( uplo == LOWER )
        for( Int j=0; j<n; ++j )
            for( Int i=j; i<n; ++i )
                ABuffer[i+j*ALDim] = negativeImagOne*GBuffer[i+j*GLDim];
    else
        for( Int j=0; j<n; ++j )
            for( Int i=0; i<=j; ++i )
                ABuffer[i+j*ALDim] = negativeImagOne*GBuffer[i+j*GLDim];

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, A, wImag, Z, a, b );
}

template<typename R>
inline void
SkewHermitianEig
( UpperOrLower uplo, 
  Matrix<Complex<R> >& G, Matrix<R>& wImag, Matrix<Complex<R> >& Z,
  Int a, Int b )
{
#ifndef RELEASE
    CallStackEntry entry("SkewHermitianEig");
#endif
    if( G.Height() != G.Width() )
        LogicError("Skew-Hermitian matrices must be square");
    
    // Make G Hermitian by scaling by -i
    const Complex<R> negativeImagOne(0,-1);
    ScaleTrapezoid( negativeImagOne, uplo, G );

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, G, wImag, Z, a, b );
}

template<typename R>
inline void
SkewHermitianEig
( UpperOrLower uplo, 
  DistMatrix<R>& G,
  DistMatrix<R,VR,STAR>& wImag,
  DistMatrix<Complex<R> >& Z,
  Int a, Int b )
{
#ifndef RELEASE
    CallStackEntry entry("SkewHermitianEig");
#endif
    EnsurePMRRR();
    if( G.Height() != G.Width() )
        LogicError("Skew-Hermitian matrices must be square");

    Int n = G.Height();
    const Grid& grid = G.Grid();

    DistMatrix<Complex<R> > A(grid);
    A.AlignWith( G.DistData() );
    A.ResizeTo( n, n );

    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();
    const Int ALDim = A.LDim();
    const Int GLDim = G.LDim();    
    const Complex<R> negativeImagOne(0,-1);
    const R* GBuffer = G.Buffer();
    Complex<R>* ABuffer = A.Buffer();
    // Just copy the entire local matrix instead of worrying about symmetry
    for( Int j=0; j<localWidth; ++j )
        for( Int i=0; i<localHeight; ++i )
            ABuffer[i+j*ALDim] = negativeImagOne*GBuffer[i+j*GLDim];

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, A, wImag, Z, a, b );
}

template<typename R>
inline void
SkewHermitianEig
( UpperOrLower uplo, 
  DistMatrix<Complex<R> >& G,
  DistMatrix<R,VR,STAR>& wImag,
  DistMatrix<Complex<R> >& Z,
  Int a, Int b )
{
#ifndef RELEASE
    CallStackEntry entry("SkewHermitianEig");
#endif
    EnsurePMRRR();
    if( G.Height() != G.Width() )
        LogicError("Skew-Hermitian matrices must be square");
    
    // Make G Hermitian by scaling by -i
    const Complex<R> negativeImagOne(0,-1);
    ScaleTrapezoid( negativeImagOne, uplo, G );

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, G, wImag, Z, a, b );
}

//----------------------------------------------------------------------------//
// Grab the eigenvalues in the interval i(a,b]                                //
//----------------------------------------------------------------------------//

template<typename R>
inline void
SkewHermitianEig( UpperOrLower uplo, Matrix<R>& G, Matrix<R>& wImag, R a, R b )
{
#ifndef RELEASE
    CallStackEntry entry("SkewHermitianEig");
#endif
    if( G.Height() != G.Width() )
        LogicError("Skew-Hermitian matrices must be square");
    Int n = G.Height();

    Matrix<Complex<R> > A( n, n );
    const Int ALDim = A.LDim();
    const Int GLDim = G.LDim();    
    const Complex<R> negativeImagOne(0,-1);
    const R* GBuffer = G.Buffer();
    Complex<R>* ABuffer = A.Buffer();
    if( uplo == LOWER )
        for( Int j=0; j<n; ++j )
            for( Int i=j; i<n; ++i )
                ABuffer[i+j*ALDim] = negativeImagOne*GBuffer[i+j*GLDim];
    else
        for( Int j=0; j<n; ++j )
            for( Int i=0; i<=j; ++i )
                ABuffer[i+j*ALDim] = negativeImagOne*GBuffer[i+j*GLDim];

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, A, wImag, a, b );
}

template<typename R>
inline void
SkewHermitianEig
( UpperOrLower uplo, Matrix<Complex<R> >& G, Matrix<R>& wImag, R a, R b )
{
#ifndef RELEASE
    CallStackEntry entry("SkewHermitianEig");
#endif
    if( G.Height() != G.Width() )
        LogicError("Skew-Hermitian matrices must be square");
    
    // Make G Hermitian by scaling by -i
    const Complex<R> negativeImagOne(0,-1);
    ScaleTrapezoid( negativeImagOne, uplo, G );

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, G, wImag, a, b );
}

template<typename R>
inline void
SkewHermitianEig
( UpperOrLower uplo, 
  DistMatrix<R>& G,
  DistMatrix<R,VR,STAR>& wImag,
  R a, R b )
{
#ifndef RELEASE
    CallStackEntry entry("SkewHermitianEig");
#endif
    EnsurePMRRR();
    if( G.Height() != G.Width() )
        LogicError("Skew-Hermitian matrices must be square");

    Int n = G.Height();
    const Grid& grid = G.Grid();

    DistMatrix<Complex<R> > A(grid);
    A.AlignWith( G.DistData() );
    A.ResizeTo( n, n );

    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();
    const Int ALDim = A.LDim();
    const Int GLDim = G.LDim();    
    const Complex<R> negativeImagOne(0,-1);
    const R* GBuffer = G.Buffer();
    Complex<R>* ABuffer = A.Buffer();
    // Just copy the entire local matrix instead of worrying about symmetry
    for( Int j=0; j<localWidth; ++j )
        for( Int i=0; i<localHeight; ++i )
            ABuffer[i+j*ALDim] = negativeImagOne*GBuffer[i+j*GLDim];

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, A, wImag, a, b );
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
    CallStackEntry entry("SkewHermitianEig");
#endif
    EnsurePMRRR();
    if( G.Height() != G.Width() )
        LogicError("Skew-Hermitian matrices must be square");
    
    // Make G Hermitian by scaling by -i
    const Complex<R> negativeImagOne(0,-1);
    ScaleTrapezoid( negativeImagOne, uplo, G );

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, G, wImag, a, b );
}

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
    CallStackEntry entry("SkewHermitianEig");
#endif
    if( G.Height() != G.Width() )
        LogicError("SkewHermitian matrices must be square");
    Int n = G.Height();

    Matrix<Complex<R> > A( n, n );
    const Int ALDim = A.LDim();
    const Int GLDim = G.LDim();    
    const Complex<R> negativeImagOne(0,-1);
    const R* GBuffer = G.Buffer();
    Complex<R>* ABuffer = A.Buffer();
    if( uplo == LOWER )
        for( Int j=0; j<n; ++j )
            for( Int i=j; i<n; ++i )
                ABuffer[i+j*ALDim] = negativeImagOne*GBuffer[i+j*GLDim];
    else
        for( Int j=0; j<n; ++j )
            for( Int i=0; i<=j; ++i )
                ABuffer[i+j*ALDim] = negativeImagOne*GBuffer[i+j*GLDim];

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, A, wImag, Z, a, b );
}

template<typename R>
inline void
SkewHermitianEig
( UpperOrLower uplo, 
  Matrix<Complex<R> >& G, Matrix<R>& wImag, Matrix<Complex<R> >& Z,
  R a, R b )
{
#ifndef RELEASE
    CallStackEntry entry("SkewHermitianEig");
#endif
    if( G.Height() != G.Width() )
        LogicError("Skew-Hermitian matrices must be square");
    
    // Make G Hermitian by scaling by -i
    const Complex<R> negativeImagOne(0,-1);
    ScaleTrapezoid( negativeImagOne, uplo, G );

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, G, wImag, Z, a, b );
}

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
    CallStackEntry entry("SkewHermitianEig");
#endif
    EnsurePMRRR();
    if( G.Height() != G.Width() )
        LogicError("SkewHermitian matrices must be square");

    Int n = G.Height();
    const Grid& grid = G.Grid();

    DistMatrix<Complex<R> > A(grid);
    A.AlignWith( G.DistData() );
    A.ResizeTo( n, n );

    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();
    const Int ALDim = A.LDim();
    const Int GLDim = G.LDim();    
    const Complex<R> negativeImagOne(0,-1);
    const R* GBuffer = G.Buffer();
    Complex<R>* ABuffer = A.Buffer();
    // Just copy the entire local matrix instead of worrying about symmetry
    for( Int j=0; j<localWidth; ++j )
        for( Int i=0; i<localHeight; ++i )
            ABuffer[i+j*ALDim] = negativeImagOne*GBuffer[i+j*GLDim];

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, A, wImag, Z, a, b );
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
    CallStackEntry entry("SkewHermitianEig");
#endif
    if( G.Height() != G.Width() )
        LogicError("Skew-Hermitian matrices must be square");
    
    // Make G Hermitian by scaling by -i
    const Complex<R> negativeImagOne(0,-1);
    ScaleTrapezoid( negativeImagOne, uplo, G );

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, G, wImag, Z, a, b );
}

} // namespace elem

#endif // ifndef ELEM_LAPACK_SKEWHERMITIANEIG_HPP
