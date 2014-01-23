/*
   Copyright (c) 2009-2014, Jack Poulson
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

template<typename Real>
inline void
SkewHermitianEig
( UpperOrLower uplo, Matrix<Real>& G, Matrix<Real>& wImag, 
  SortType sort=UNSORTED )
{
    DEBUG_ONLY(CallStackEntry cse("SkewHermitianEig"))
    if( G.Height() != G.Width() )
        LogicError("SkewHermitian matrices must be square");
    Int n = G.Height();
    Matrix<Complex<Real>> A( n, n );
    const Int ALDim = A.LDim();
    const Int GLDim = G.LDim();    
    const Complex<Real> negativeImagOne(0,-1);
    const Real* GBuffer = G.Buffer();
    Complex<Real>* ABuffer = A.Buffer();
    if( uplo == LOWER )
        for( Int j=0; j<n; ++j )
            for( Int i=j; i<n; ++i )
                ABuffer[i+j*ALDim] = negativeImagOne*GBuffer[i+j*GLDim];
    else
        for( Int j=0; j<n; ++j )
            for( Int i=0; i<=j; ++i )
                ABuffer[i+j*ALDim] = negativeImagOne*GBuffer[i+j*GLDim];

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, A, wImag, sort );
}

template<typename Real>
inline void
SkewHermitianEig
( UpperOrLower uplo, Matrix<Complex<Real> >& G, 
  Matrix<Real>& wImag, SortType sort=UNSORTED )
{
    DEBUG_ONLY(CallStackEntry cse("SkewHermitianEig"))
    if( G.Height() != G.Width() )
        LogicError("Skew-Hermitian matrices must be square");
    
    // Make G Hermitian by scaling by -i
    const Complex<Real> negativeImagOne(0,-1);
    ScaleTrapezoid( negativeImagOne, uplo, G );

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, G, wImag, sort );
}

template<typename Real>
inline void
SkewHermitianEig
( UpperOrLower uplo, DistMatrix<Real>& G,
  DistMatrix<Real,VR,STAR>& wImag, SortType sort=UNSORTED )
{
    DEBUG_ONLY(CallStackEntry cse("SkewHermitianEig"))
    EnsurePMRRR();
    if( G.Height() != G.Width() )
        LogicError("SkewHermitian matrices must be square");

    Int n = G.Height();
    const Grid& grid = G.Grid();

    DistMatrix<Complex<Real>> A(grid);
    A.AlignWith( G.DistData() );
    A.Resize( n, n );

    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();
    const Int ALDim = A.LDim();
    const Int GLDim = G.LDim();    
    const Complex<Real> negativeImagOne(0,-1);
    const Real* GBuffer = G.Buffer();
    Complex<Real>* ABuffer = A.Buffer();
    // Just copy the entire local matrix instead of worrying about symmetry
    for( Int j=0; j<localWidth; ++j )
        for( Int i=0; i<localHeight; ++i )
            ABuffer[i+j*ALDim] = negativeImagOne*GBuffer[i+j*GLDim];

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, A, wImag, sort );
}

template<typename Real>
inline void
SkewHermitianEig
( UpperOrLower uplo, DistMatrix<Complex<Real> >& G,
  DistMatrix<Real,VR,STAR>& wImag, SortType sort=UNSORTED )
{
    DEBUG_ONLY(CallStackEntry cse("SkewHermitianEig"))
    if( G.Height() != G.Width() )
        LogicError("Skew-Hermitian matrices must be square");
    
    // Make G Hermitian by scaling by -i
    const Complex<Real> negativeImagOne(0,-1);
    ScaleTrapezoid( negativeImagOne, uplo, G );

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, G, wImag, sort );
}

//----------------------------------------------------------------------------//
// Grab the full set of eigenpairs                                            //
//----------------------------------------------------------------------------//

template<typename Real>
inline void
SkewHermitianEig
( UpperOrLower uplo, Matrix<Real>& G, 
  Matrix<Real>& wImag, Matrix<Complex<Real> >& Z,
  SortType sort=UNSORTED )
{
    DEBUG_ONLY(CallStackEntry cse("SkewHermitianEig"))
    if( G.Height() != G.Width() )
        LogicError("SkewHermitian matrices must be square");
    Int n = G.Height();
    Matrix<Complex<Real>> A( n, n );

    const Int ALDim = A.LDim();
    const Int GLDim = G.LDim();
    const Complex<Real> negativeImagOne(0,-1);
    const Real* GBuffer = G.Buffer();
    Complex<Real>* ABuffer = A.Buffer();
    if( uplo == LOWER )
        for( Int j=0; j<n; ++j )
            for( Int i=j; i<n; ++i )
                ABuffer[i+j*ALDim] = negativeImagOne*GBuffer[i+j*GLDim];
    else
        for( Int j=0; j<n; ++j )
            for( Int i=0; i<=j; ++i )
                ABuffer[i+j*ALDim] = negativeImagOne*GBuffer[i+j*GLDim];

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, A, wImag, Z, sort );
}

template<typename Real>
inline void
SkewHermitianEig
( UpperOrLower uplo, Matrix<Complex<Real> >& G, 
  Matrix<Real>& wImag, Matrix<Complex<Real> >& Z,
  SortType sort=UNSORTED )
{
    DEBUG_ONLY(CallStackEntry cse("SkewHermitianEig"))
    if( G.Height() != G.Width() )
        LogicError("Skew-Hermitian matrices must be square");

    // Make G Hermitian by scaling by -i
    const Complex<Real> negativeImagOne(0,-1);
    ScaleTrapezoid( negativeImagOne, uplo, G );

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, G, wImag, Z, sort );
}

template<typename Real>
inline void
SkewHermitianEig
( UpperOrLower uplo, DistMatrix<Real>& G,
  DistMatrix<Real,VR,STAR>& wImag, DistMatrix<Complex<Real> >& Z,
  SortType sort=UNSORTED )
{
    DEBUG_ONLY(CallStackEntry cse("SkewHermitianEig"))
    EnsurePMRRR();
    if( G.Height() != G.Width() )
        LogicError("SkewHermitian matrices must be square");

    Int n = G.Height();
    const Grid& grid = G.Grid();

    DistMatrix<Complex<Real>> A(grid);
    A.AlignWith( G.DistData() );
    A.Resize( n, n );

    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();
    const Int ALDim = A.LDim();
    const Int GLDim = G.LDim();
    const Complex<Real> negativeImagOne(0,-1);
    const Real* GBuffer = G.Buffer();
    Complex<Real>* ABuffer = A.Buffer();
    // Just copy the entire local matrix instead of worrying about symmetry
    for( Int j=0; j<localWidth; ++j )
        for( Int i=0; i<localHeight; ++i )
            ABuffer[i+j*ALDim] = negativeImagOne*GBuffer[i+j*GLDim];

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, A, wImag, Z, sort );
}

template<typename Real>
inline void
SkewHermitianEig
( UpperOrLower uplo, DistMatrix<Complex<Real> >& G,
  DistMatrix<Real,VR,STAR>& wImag, DistMatrix<Complex<Real> >& Z,
  SortType sort=UNSORTED )
{
    DEBUG_ONLY(CallStackEntry cse("SkewHermitianEig"))
    EnsurePMRRR();
    if( G.Height() != G.Width() )
        LogicError("Skew-Hermitian matrices must be square");

    // Make G Hermitian by scaling by -i
    const Complex<Real> negativeImagOne(0,-1);
    ScaleTrapezoid( negativeImagOne, uplo, G );

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, G, wImag, Z, sort );
}

//----------------------------------------------------------------------------//
// Grab the eigenvalues with indices in a specified range                     //
//----------------------------------------------------------------------------//

template<typename Real>
inline void
SkewHermitianEig
( UpperOrLower uplo, Matrix<Real>& G, 
  Matrix<Real>& wImag, Int a, Int b, SortType sort=UNSORTED )
{
    DEBUG_ONLY(CallStackEntry cse("SkewHermitianEig"))
    if( G.Height() != G.Width() )
        LogicError("Skew-Hermitian matrices must be square");
    Int n = G.Height();

    Matrix<Complex<Real> > A( n, n );
    const Int ALDim = A.LDim();
    const Int GLDim = G.LDim();    
    const Complex<Real> negativeImagOne(0,-1);
    const Real* GBuffer = G.Buffer();
    Complex<Real>* ABuffer = A.Buffer();
    if( uplo == LOWER )
        for( Int j=0; j<n; ++j )
            for( Int i=j; i<n; ++i )
                ABuffer[i+j*ALDim] = negativeImagOne*GBuffer[i+j*GLDim];
    else
        for( Int j=0; j<n; ++j )
            for( Int i=0; i<=j; ++i )
                ABuffer[i+j*ALDim] = negativeImagOne*GBuffer[i+j*GLDim];

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, A, wImag, a, b, sort );
}

template<typename Real>
inline void
SkewHermitianEig
( UpperOrLower uplo, Matrix<Complex<Real> >& G, 
  Matrix<Real>& wImag, Int a, Int b, SortType sort=UNSORTED )
{
    DEBUG_ONLY(CallStackEntry cse("SkewHermitianEig"))
    if( G.Height() != G.Width() )
        LogicError("Skew-Hermitian matrices must be square");
    
    // Make G Hermitian by scaling by -i
    const Complex<Real> negativeImagOne(0,-1);
    ScaleTrapezoid( negativeImagOne, uplo, G );

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, G, wImag, a, b, sort );
}

template<typename Real>
inline void
SkewHermitianEig
( UpperOrLower uplo, DistMatrix<Real>& G,
  DistMatrix<Real,VR,STAR>& wImag, Int a, Int b, SortType sort=UNSORTED )
{
    DEBUG_ONLY(CallStackEntry cse("SkewHermitianEig"))
    EnsurePMRRR();
    if( G.Height() != G.Width() )
        LogicError("Skew-Hermitian matrices must be square");

    Int n = G.Height();
    const Grid& grid = G.Grid();

    DistMatrix<Complex<Real> > A(grid);
    A.AlignWith( G.DistData() );
    A.Resize( n, n );

    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();
    const Int ALDim = A.LDim();
    const Int GLDim = G.LDim();    
    const Complex<Real> negativeImagOne(0,-1);
    const Real* GBuffer = G.Buffer();
    Complex<Real>* ABuffer = A.Buffer();
    // Just copy the entire local matrix instead of worrying about symmetry
    for( Int j=0; j<localWidth; ++j )
        for( Int i=0; i<localHeight; ++i )
            ABuffer[i+j*ALDim] = negativeImagOne*GBuffer[i+j*GLDim];

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, A, wImag, a, b, sort );
}

template<typename Real>
inline void
SkewHermitianEig
( UpperOrLower uplo, DistMatrix<Complex<Real> >& G,
  DistMatrix<Real,VR,STAR>& wImag, Int a, Int b, SortType sort=UNSORTED )
{
    DEBUG_ONLY(CallStackEntry cse("SkewHermitianEig"))
    EnsurePMRRR();
    if( G.Height() != G.Width() )
        LogicError("Skew-Hermitian matrices must be square");
    
    // Make G Hermitian by scaling by -i
    const Complex<Real> negativeImagOne(0,-1);
    ScaleTrapezoid( negativeImagOne, uplo, G );

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, G, wImag, a, b, sort );
}

//----------------------------------------------------------------------------//
// Grab the eigenpairs with indices in a specified range                      //
//----------------------------------------------------------------------------//

template<typename Real>
inline void
SkewHermitianEig
( UpperOrLower uplo, Matrix<Real>& G, 
  Matrix<Real>& wImag, Matrix<Complex<Real> >& Z,
  Int a, Int b, SortType sort=UNSORTED )
{
    DEBUG_ONLY(CallStackEntry cse("SkewHermitianEig"))
    if( G.Height() != G.Width() )
        LogicError("Skew-Hermitian matrices must be square");
    Int n = G.Height();

    Matrix<Complex<Real> > A( n, n );
    const Int ALDim = A.LDim();
    const Int GLDim = G.LDim();    
    const Complex<Real> negativeImagOne(0,-1);
    const Real* GBuffer = G.Buffer();
    Complex<Real>* ABuffer = A.Buffer();
    if( uplo == LOWER )
        for( Int j=0; j<n; ++j )
            for( Int i=j; i<n; ++i )
                ABuffer[i+j*ALDim] = negativeImagOne*GBuffer[i+j*GLDim];
    else
        for( Int j=0; j<n; ++j )
            for( Int i=0; i<=j; ++i )
                ABuffer[i+j*ALDim] = negativeImagOne*GBuffer[i+j*GLDim];

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, A, wImag, Z, a, b, sort );
}

template<typename Real>
inline void
SkewHermitianEig
( UpperOrLower uplo, Matrix<Complex<Real> >& G, 
  Matrix<Real>& wImag, Matrix<Complex<Real> >& Z,
  Int a, Int b, SortType sort=UNSORTED )
{
    DEBUG_ONLY(CallStackEntry cse("SkewHermitianEig"))
    if( G.Height() != G.Width() )
        LogicError("Skew-Hermitian matrices must be square");
    
    // Make G Hermitian by scaling by -i
    const Complex<Real> negativeImagOne(0,-1);
    ScaleTrapezoid( negativeImagOne, uplo, G );

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, G, wImag, Z, a, b, sort );
}

template<typename Real>
inline void
SkewHermitianEig
( UpperOrLower uplo, DistMatrix<Real>& G,
  DistMatrix<Real,VR,STAR>& wImag, DistMatrix<Complex<Real> >& Z,
  Int a, Int b, SortType sort=UNSORTED )
{
    DEBUG_ONLY(CallStackEntry cse("SkewHermitianEig"))
    EnsurePMRRR();
    if( G.Height() != G.Width() )
        LogicError("Skew-Hermitian matrices must be square");

    Int n = G.Height();
    const Grid& grid = G.Grid();

    DistMatrix<Complex<Real> > A(grid);
    A.AlignWith( G.DistData() );
    A.Resize( n, n );

    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();
    const Int ALDim = A.LDim();
    const Int GLDim = G.LDim();    
    const Complex<Real> negativeImagOne(0,-1);
    const Real* GBuffer = G.Buffer();
    Complex<Real>* ABuffer = A.Buffer();
    // Just copy the entire local matrix instead of worrying about symmetry
    for( Int j=0; j<localWidth; ++j )
        for( Int i=0; i<localHeight; ++i )
            ABuffer[i+j*ALDim] = negativeImagOne*GBuffer[i+j*GLDim];

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, A, wImag, Z, a, b, sort );
}

template<typename Real>
inline void
SkewHermitianEig
( UpperOrLower uplo, DistMatrix<Complex<Real> >& G,
  DistMatrix<Real,VR,STAR>& wImag, DistMatrix<Complex<Real> >& Z,
  Int a, Int b, SortType sort=UNSORTED )
{
    DEBUG_ONLY(CallStackEntry cse("SkewHermitianEig"))
    EnsurePMRRR();
    if( G.Height() != G.Width() )
        LogicError("Skew-Hermitian matrices must be square");
    
    // Make G Hermitian by scaling by -i
    const Complex<Real> negativeImagOne(0,-1);
    ScaleTrapezoid( negativeImagOne, uplo, G );

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, G, wImag, Z, a, b, sort );
}

//----------------------------------------------------------------------------//
// Grab the eigenvalues in the interval i(a,b]                                //
//----------------------------------------------------------------------------//

template<typename Real>
inline void
SkewHermitianEig
( UpperOrLower uplo, Matrix<Real>& G, 
  Matrix<Real>& wImag, Real a, Real b, SortType sort=UNSORTED )
{
    DEBUG_ONLY(CallStackEntry cse("SkewHermitianEig"))
    if( G.Height() != G.Width() )
        LogicError("Skew-Hermitian matrices must be square");
    Int n = G.Height();

    Matrix<Complex<Real> > A( n, n );
    const Int ALDim = A.LDim();
    const Int GLDim = G.LDim();    
    const Complex<Real> negativeImagOne(0,-1);
    const Real* GBuffer = G.Buffer();
    Complex<Real>* ABuffer = A.Buffer();
    if( uplo == LOWER )
        for( Int j=0; j<n; ++j )
            for( Int i=j; i<n; ++i )
                ABuffer[i+j*ALDim] = negativeImagOne*GBuffer[i+j*GLDim];
    else
        for( Int j=0; j<n; ++j )
            for( Int i=0; i<=j; ++i )
                ABuffer[i+j*ALDim] = negativeImagOne*GBuffer[i+j*GLDim];

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, A, wImag, a, b, sort );
}

template<typename Real>
inline void
SkewHermitianEig
( UpperOrLower uplo, Matrix<Complex<Real> >& G, 
  Matrix<Real>& wImag, Real a, Real b, SortType sort=UNSORTED )
{
    DEBUG_ONLY(CallStackEntry cse("SkewHermitianEig"))
    if( G.Height() != G.Width() )
        LogicError("Skew-Hermitian matrices must be square");
    
    // Make G Hermitian by scaling by -i
    const Complex<Real> negativeImagOne(0,-1);
    ScaleTrapezoid( negativeImagOne, uplo, G );

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, G, wImag, a, b, sort );
}

template<typename Real>
inline void
SkewHermitianEig
( UpperOrLower uplo, DistMatrix<Real>& G,
  DistMatrix<Real,VR,STAR>& wImag, Real a, Real b, SortType sort=UNSORTED )
{
    DEBUG_ONLY(CallStackEntry cse("SkewHermitianEig"))
    EnsurePMRRR();
    if( G.Height() != G.Width() )
        LogicError("Skew-Hermitian matrices must be square");

    Int n = G.Height();
    const Grid& grid = G.Grid();

    DistMatrix<Complex<Real> > A(grid);
    A.AlignWith( G.DistData() );
    A.Resize( n, n );

    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();
    const Int ALDim = A.LDim();
    const Int GLDim = G.LDim();    
    const Complex<Real> negativeImagOne(0,-1);
    const Real* GBuffer = G.Buffer();
    Complex<Real>* ABuffer = A.Buffer();
    // Just copy the entire local matrix instead of worrying about symmetry
    for( Int j=0; j<localWidth; ++j )
        for( Int i=0; i<localHeight; ++i )
            ABuffer[i+j*ALDim] = negativeImagOne*GBuffer[i+j*GLDim];

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, A, wImag, a, b, sort );
}

template<typename Real>
inline void
SkewHermitianEig
( UpperOrLower uplo, DistMatrix<Complex<Real> >& G,
  DistMatrix<Real,VR,STAR>& wImag,
  Real a, Real b, SortType sort=UNSORTED )
{
    DEBUG_ONLY(CallStackEntry cse("SkewHermitianEig"))
    EnsurePMRRR();
    if( G.Height() != G.Width() )
        LogicError("Skew-Hermitian matrices must be square");
    
    // Make G Hermitian by scaling by -i
    const Complex<Real> negativeImagOne(0,-1);
    ScaleTrapezoid( negativeImagOne, uplo, G );

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, G, wImag, a, b, sort );
}

//----------------------------------------------------------------------------//
// Grab the eigenpairs with eigenvalues in the interval i(a,b]                //
//----------------------------------------------------------------------------//

template<typename Real>
inline void
SkewHermitianEig
( UpperOrLower uplo, Matrix<Real>& G, 
  Matrix<Real>& wImag, Matrix<Complex<Real> >& Z,
  Real a, Real b, SortType sort=UNSORTED )
{
    DEBUG_ONLY(CallStackEntry cse("SkewHermitianEig"))
    if( G.Height() != G.Width() )
        LogicError("SkewHermitian matrices must be square");
    Int n = G.Height();

    Matrix<Complex<Real> > A( n, n );
    const Int ALDim = A.LDim();
    const Int GLDim = G.LDim();    
    const Complex<Real> negativeImagOne(0,-1);
    const Real* GBuffer = G.Buffer();
    Complex<Real>* ABuffer = A.Buffer();
    if( uplo == LOWER )
        for( Int j=0; j<n; ++j )
            for( Int i=j; i<n; ++i )
                ABuffer[i+j*ALDim] = negativeImagOne*GBuffer[i+j*GLDim];
    else
        for( Int j=0; j<n; ++j )
            for( Int i=0; i<=j; ++i )
                ABuffer[i+j*ALDim] = negativeImagOne*GBuffer[i+j*GLDim];

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, A, wImag, Z, a, b, sort );
}

template<typename Real>
inline void
SkewHermitianEig
( UpperOrLower uplo, Matrix<Complex<Real> >& G, 
  Matrix<Real>& wImag, Matrix<Complex<Real> >& Z,
  Real a, Real b, SortType sort=UNSORTED )
{
    DEBUG_ONLY(CallStackEntry cse("SkewHermitianEig"))
    if( G.Height() != G.Width() )
        LogicError("Skew-Hermitian matrices must be square");
    
    // Make G Hermitian by scaling by -i
    const Complex<Real> negativeImagOne(0,-1);
    ScaleTrapezoid( negativeImagOne, uplo, G );

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, G, wImag, Z, a, b, sort );
}

template<typename Real>
inline void
SkewHermitianEig
( UpperOrLower uplo, DistMatrix<Real>& G,
  DistMatrix<Real,VR,STAR>& wImag, DistMatrix<Complex<Real> >& Z,
  Real a, Real b, SortType sort=UNSORTED )
{
    DEBUG_ONLY(CallStackEntry cse("SkewHermitianEig"))
    EnsurePMRRR();
    if( G.Height() != G.Width() )
        LogicError("SkewHermitian matrices must be square");

    Int n = G.Height();
    const Grid& grid = G.Grid();

    DistMatrix<Complex<Real> > A(grid);
    A.AlignWith( G.DistData() );
    A.Resize( n, n );

    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();
    const Int ALDim = A.LDim();
    const Int GLDim = G.LDim();    
    const Complex<Real> negativeImagOne(0,-1);
    const Real* GBuffer = G.Buffer();
    Complex<Real>* ABuffer = A.Buffer();
    // Just copy the entire local matrix instead of worrying about symmetry
    for( Int j=0; j<localWidth; ++j )
        for( Int i=0; i<localHeight; ++i )
            ABuffer[i+j*ALDim] = negativeImagOne*GBuffer[i+j*GLDim];

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, A, wImag, Z, a, b, sort );
}

template<typename Real>
inline void
SkewHermitianEig
( UpperOrLower uplo, DistMatrix<Complex<Real> >& G,
  DistMatrix<Real,VR,STAR>& wImag, DistMatrix<Complex<Real> >& Z,
  Real a, Real b, SortType sort=UNSORTED )
{
    DEBUG_ONLY(CallStackEntry cse("SkewHermitianEig"))
    if( G.Height() != G.Width() )
        LogicError("Skew-Hermitian matrices must be square");
    
    // Make G Hermitian by scaling by -i
    const Complex<Real> negativeImagOne(0,-1);
    ScaleTrapezoid( negativeImagOne, uplo, G );

    // Perform the Hermitian eigensolve
    HermitianEig( uplo, G, wImag, Z, a, b, sort );
}

} // namespace elem

#endif // ifndef ELEM_LAPACK_SKEWHERMITIANEIG_HPP
