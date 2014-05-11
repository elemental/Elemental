/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_SKEWHERMITIANEIG_HPP
#define ELEM_SKEWHERMITIANEIG_HPP

#include ELEM_SCALETRAPEZOID_INC

namespace elem {

// Return the full set of eigenvalues
// ==================================

template<typename Real>
inline void
SkewHermitianEig
( UpperOrLower uplo, Matrix<Real>& G, Matrix<Real>& wImag, 
  SortType sort=UNSORTED, 
  const HermitianEigCtrl<Real>& ctrl=HermitianEigCtrl<Real>() )
{
    DEBUG_ONLY(CallStackEntry cse("SkewHermitianEig"))
    Matrix<Complex<Real>> A;
    Copy( G, A );
    ScaleTrapezoid( Complex<Real>(0,-1), uplo, A );
    HermitianEig( uplo, A, wImag, sort, ctrl );
}

template<typename Real>
inline void
SkewHermitianEig
( UpperOrLower uplo, Matrix<Complex<Real>>& G, 
  Matrix<Real>& wImag, SortType sort=UNSORTED,
  const HermitianEigCtrl<Real>& ctrl=HermitianEigCtrl<Real>() )
{
    DEBUG_ONLY(CallStackEntry cse("SkewHermitianEig"))
    ScaleTrapezoid( Complex<Real>(0,-1), uplo, G );
    HermitianEig( uplo, G, wImag, sort, ctrl );
}

template<typename Real>
inline void
SkewHermitianEig
( UpperOrLower uplo, DistMatrix<Real>& G,
  DistMatrix<Real,VR,STAR>& wImag, SortType sort=UNSORTED,
  const HermitianEigCtrl<Real>& ctrl=HermitianEigCtrl<Real>() )
{
    DEBUG_ONLY(CallStackEntry cse("SkewHermitianEig"))
    DistMatrix<Complex<Real>> A(G.Grid());
    Copy( G, A );
    ScaleTrapezoid( Complex<Real>(0,-1), uplo, A );
    HermitianEig( uplo, A, wImag, sort, ctrl );
}

template<typename Real>
inline void
SkewHermitianEig
( UpperOrLower uplo, DistMatrix<Complex<Real>>& G,
  DistMatrix<Real,VR,STAR>& wImag, SortType sort=UNSORTED,
  const HermitianEigCtrl<Real>& ctrl=HermitianEigCtrl<Real>() )
{
    DEBUG_ONLY(CallStackEntry cse("SkewHermitianEig"))
    ScaleTrapezoid( Complex<Real>(0,-1), uplo, G );
    HermitianEig( uplo, G, wImag, sort, ctrl );
}

// Return the full set of eigenpairs
// =================================

template<typename Real>
inline void
SkewHermitianEig
( UpperOrLower uplo, Matrix<Real>& G, 
  Matrix<Real>& wImag, Matrix<Complex<Real>>& Z,
  SortType sort=UNSORTED,
  const HermitianEigCtrl<Real>& ctrl=HermitianEigCtrl<Real>() )
{
    DEBUG_ONLY(CallStackEntry cse("SkewHermitianEig"))
    Matrix<Complex<Real>> A;
    Copy( G, A );
    ScaleTrapezoid( Complex<Real>(0,-1), uplo, A );
    HermitianEig( uplo, A, wImag, Z, sort, ctrl );
}

template<typename Real>
inline void
SkewHermitianEig
( UpperOrLower uplo, Matrix<Complex<Real>>& G, 
  Matrix<Real>& wImag, Matrix<Complex<Real>>& Z,
  SortType sort=UNSORTED,
  const HermitianEigCtrl<Real>& ctrl=HermitianEigCtrl<Real>() )
{
    DEBUG_ONLY(CallStackEntry cse("SkewHermitianEig"))
    ScaleTrapezoid( Complex<Real>(0,-1), uplo, G );
    HermitianEig( uplo, G, wImag, Z, sort, ctrl );
}

template<typename Real>
inline void
SkewHermitianEig
( UpperOrLower uplo, DistMatrix<Real>& G,
  DistMatrix<Real,VR,STAR>& wImag, DistMatrix<Complex<Real>>& Z,
  SortType sort=UNSORTED,
  const HermitianEigCtrl<Real>& ctrl=HermitianEigCtrl<Real>() )
{
    DEBUG_ONLY(CallStackEntry cse("SkewHermitianEig"))
    DistMatrix<Complex<Real>> A(G.Grid());
    Copy( G, A ); 
    ScaleTrapezoid( Complex<Real>(0,-1), uplo, A );
    HermitianEig( uplo, A, wImag, Z, sort, ctrl );
}

template<typename Real>
inline void
SkewHermitianEig
( UpperOrLower uplo, DistMatrix<Complex<Real>>& G,
  DistMatrix<Real,VR,STAR>& wImag, DistMatrix<Complex<Real>>& Z,
  SortType sort=UNSORTED,
  const HermitianEigCtrl<Real>& ctrl=HermitianEigCtrl<Real>() )
{
    DEBUG_ONLY(CallStackEntry cse("SkewHermitianEig"))
    ScaleTrapezoid( Complex<Real>(0,-1), uplo, G );
    HermitianEig( uplo, G, wImag, Z, sort, ctrl );
}

// Return the eigenvalues with indices in a specified range
// ========================================================

template<typename Real>
inline void
SkewHermitianEig
( UpperOrLower uplo, Matrix<Real>& G, 
  Matrix<Real>& wImag, Int a, Int b, SortType sort=UNSORTED,
  const HermitianEigCtrl<Real>& ctrl=HermitianEigCtrl<Real>() )
{
    DEBUG_ONLY(CallStackEntry cse("SkewHermitianEig"))
    Matrix<Complex<Real>> A;
    Copy( G, A );
    ScaleTrapezoid( Complex<Real>(0,-1), uplo, A );
    HermitianEig( uplo, A, wImag, a, b, sort, ctrl );
}

template<typename Real>
inline void
SkewHermitianEig
( UpperOrLower uplo, Matrix<Complex<Real>>& G, 
  Matrix<Real>& wImag, Int a, Int b, SortType sort=UNSORTED,
  const HermitianEigCtrl<Real>& ctrl=HermitianEigCtrl<Real>() )
{
    DEBUG_ONLY(CallStackEntry cse("SkewHermitianEig"))
    ScaleTrapezoid( Complex<Real>(0,-1), uplo, G );
    HermitianEig( uplo, G, wImag, a, b, sort, ctrl );
}

template<typename Real>
inline void
SkewHermitianEig
( UpperOrLower uplo, DistMatrix<Real>& G,
  DistMatrix<Real,VR,STAR>& wImag, Int a, Int b, SortType sort=UNSORTED,
  const HermitianEigCtrl<Real>& ctrl=HermitianEigCtrl<Real>() )
{
    DEBUG_ONLY(CallStackEntry cse("SkewHermitianEig"))
    DistMatrix<Complex<Real>> A(G.Grid());
    Copy( G, A );
    ScaleTrapezoid( Complex<Real>(0,-1), uplo, A );
    HermitianEig( uplo, A, wImag, a, b, sort, ctrl );
}

template<typename Real>
inline void
SkewHermitianEig
( UpperOrLower uplo, DistMatrix<Complex<Real>>& G,
  DistMatrix<Real,VR,STAR>& wImag, Int a, Int b, SortType sort=UNSORTED,
  const HermitianEigCtrl<Real>& ctrl=HermitianEigCtrl<Real>() )
{
    DEBUG_ONLY(CallStackEntry cse("SkewHermitianEig"))
    ScaleTrapezoid( Complex<Real>(0,-1), uplo, G );
    HermitianEig( uplo, G, wImag, a, b, sort, ctrl );
}

// Return the eigenpairs with indices in a specified range
// =======================================================

template<typename Real>
inline void
SkewHermitianEig
( UpperOrLower uplo, Matrix<Real>& G, 
  Matrix<Real>& wImag, Matrix<Complex<Real>>& Z,
  Int a, Int b, SortType sort=UNSORTED,
  const HermitianEigCtrl<Real>& ctrl=HermitianEigCtrl<Real>() )
{
    DEBUG_ONLY(CallStackEntry cse("SkewHermitianEig"))
    Matrix<Complex<Real>> A;
    Copy( G, A );
    ScaleTrapezoid( Complex<Real>(0,-1), uplo, A );
    HermitianEig( uplo, A, wImag, Z, a, b, sort, ctrl );
}

template<typename Real>
inline void
SkewHermitianEig
( UpperOrLower uplo, Matrix<Complex<Real>>& G, 
  Matrix<Real>& wImag, Matrix<Complex<Real>>& Z,
  Int a, Int b, SortType sort=UNSORTED,
  const HermitianEigCtrl<Real>& ctrl=HermitianEigCtrl<Real>() )
{
    DEBUG_ONLY(CallStackEntry cse("SkewHermitianEig"))
    ScaleTrapezoid( Complex<Real>(0,-1), uplo, G );
    HermitianEig( uplo, G, wImag, Z, a, b, sort, ctrl );
}

template<typename Real>
inline void
SkewHermitianEig
( UpperOrLower uplo, DistMatrix<Real>& G,
  DistMatrix<Real,VR,STAR>& wImag, DistMatrix<Complex<Real>>& Z,
  Int a, Int b, SortType sort=UNSORTED,
  const HermitianEigCtrl<Real>& ctrl=HermitianEigCtrl<Real>() )
{
    DEBUG_ONLY(CallStackEntry cse("SkewHermitianEig"))
    DistMatrix<Complex<Real>> A(G.Grid());
    Copy( G, A );
    ScaleTrapezoid( Complex<Real>(0,-1), uplo, A );
    HermitianEig( uplo, A, wImag, Z, a, b, sort, ctrl );
}

template<typename Real>
inline void
SkewHermitianEig
( UpperOrLower uplo, DistMatrix<Complex<Real>>& G,
  DistMatrix<Real,VR,STAR>& wImag, DistMatrix<Complex<Real>>& Z,
  Int a, Int b, SortType sort=UNSORTED,
  const HermitianEigCtrl<Real>& ctrl=HermitianEigCtrl<Real>() )
{
    DEBUG_ONLY(CallStackEntry cse("SkewHermitianEig"))
    ScaleTrapezoid( Complex<Real>(0,-1), uplo, G );
    HermitianEig( uplo, G, wImag, Z, a, b, sort, ctrl );
}

// Return the eigenvalues in the interval i(a,b]
// =============================================

template<typename Real>
inline void
SkewHermitianEig
( UpperOrLower uplo, Matrix<Real>& G, 
  Matrix<Real>& wImag, Real a, Real b, SortType sort=UNSORTED,
  const HermitianEigCtrl<Real>& ctrl=HermitianEigCtrl<Real>() )
{
    DEBUG_ONLY(CallStackEntry cse("SkewHermitianEig"))
    Matrix<Complex<Real>> A;
    Copy( G, A );
    ScaleTrapezoid( Complex<Real>(0,-1), uplo, A );
    HermitianEig( uplo, A, wImag, a, b, sort, ctrl );
}

template<typename Real>
inline void
SkewHermitianEig
( UpperOrLower uplo, Matrix<Complex<Real>>& G, 
  Matrix<Real>& wImag, Real a, Real b, SortType sort=UNSORTED,
  const HermitianEigCtrl<Real>& ctrl=HermitianEigCtrl<Real>() )
{
    DEBUG_ONLY(CallStackEntry cse("SkewHermitianEig"))
    ScaleTrapezoid( Complex<Real>(0,-1), uplo, G );
    HermitianEig( uplo, G, wImag, a, b, sort, ctrl );
}

template<typename Real>
inline void
SkewHermitianEig
( UpperOrLower uplo, DistMatrix<Real>& G,
  DistMatrix<Real,VR,STAR>& wImag, Real a, Real b, SortType sort=UNSORTED,
  const HermitianEigCtrl<Real>& ctrl=HermitianEigCtrl<Real>() )
{
    DEBUG_ONLY(CallStackEntry cse("SkewHermitianEig"))
    DistMatrix<Complex<Real>> A(G.Grid());
    Copy( G, A );
    ScaleTrapezoid( Complex<Real>(0,-1), uplo, A );
    HermitianEig( uplo, A, wImag, a, b, sort, ctrl );
}

template<typename Real>
inline void
SkewHermitianEig
( UpperOrLower uplo, DistMatrix<Complex<Real>>& G,
  DistMatrix<Real,VR,STAR>& wImag,
  Real a, Real b, SortType sort=UNSORTED,
  const HermitianEigCtrl<Real>& ctrl=HermitianEigCtrl<Real>() )
{
    DEBUG_ONLY(CallStackEntry cse("SkewHermitianEig"))
    ScaleTrapezoid( Complex<Real>(0,-1), uplo, G );
    HermitianEig( uplo, G, wImag, a, b, sort, ctrl );
}

// Return the eigenpairs with eigenvalues in the interval i(a,b]
// =============================================================

template<typename Real>
inline void
SkewHermitianEig
( UpperOrLower uplo, Matrix<Real>& G, 
  Matrix<Real>& wImag, Matrix<Complex<Real>>& Z,
  Real a, Real b, SortType sort=UNSORTED,
  const HermitianEigCtrl<Real>& ctrl=HermitianEigCtrl<Real>() )
{
    DEBUG_ONLY(CallStackEntry cse("SkewHermitianEig"))
    Matrix<Complex<Real>> A;
    Copy( G, A );
    ScaleTrapezoid( Complex<Real>(0,-1), uplo, A );
    HermitianEig( uplo, A, wImag, Z, a, b, sort, ctrl );
}

template<typename Real>
inline void
SkewHermitianEig
( UpperOrLower uplo, Matrix<Complex<Real>>& G, 
  Matrix<Real>& wImag, Matrix<Complex<Real>>& Z,
  Real a, Real b, SortType sort=UNSORTED,
  const HermitianEigCtrl<Real>& ctrl=HermitianEigCtrl<Real>() )
{
    DEBUG_ONLY(CallStackEntry cse("SkewHermitianEig"))
    ScaleTrapezoid( Complex<Real>(0,-1), uplo, G );
    HermitianEig( uplo, G, wImag, Z, a, b, sort, ctrl );
}

template<typename Real>
inline void
SkewHermitianEig
( UpperOrLower uplo, DistMatrix<Real>& G,
  DistMatrix<Real,VR,STAR>& wImag, DistMatrix<Complex<Real>>& Z,
  Real a, Real b, SortType sort=UNSORTED,
  const HermitianEigCtrl<Real>& ctrl=HermitianEigCtrl<Real>() )
{
    DEBUG_ONLY(CallStackEntry cse("SkewHermitianEig"))
    DistMatrix<Complex<Real>> A(G.Grid());
    Copy( G, A );
    ScaleTrapezoid( Complex<Real>(0,-1), uplo, A );
    HermitianEig( uplo, A, wImag, Z, a, b, sort, ctrl );
}

template<typename Real>
inline void
SkewHermitianEig
( UpperOrLower uplo, DistMatrix<Complex<Real>>& G,
  DistMatrix<Real,VR,STAR>& wImag, DistMatrix<Complex<Real>>& Z,
  Real a, Real b, SortType sort=UNSORTED,
  const HermitianEigCtrl<Real>& ctrl=HermitianEigCtrl<Real>() )
{
    DEBUG_ONLY(CallStackEntry cse("SkewHermitianEig"))
    ScaleTrapezoid( Complex<Real>(0,-1), uplo, G );
    HermitianEig( uplo, G, wImag, Z, a, b, sort, ctrl );
}

} // namespace elem

#endif // ifndef ELEM_SKEWHERMITIANEIG_HPP
