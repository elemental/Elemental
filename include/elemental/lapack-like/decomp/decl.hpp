/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_DECOMP_DECL_HPP
#define ELEM_DECOMP_DECL_HPP

#include ELEM_CONDENSE_DECL_INC

namespace elem {

namespace hermitian_gen_definite_eig_type_wrapper {
enum HermitianGenDefiniteEigType
{
    AXBX=1,
    ABX=2,
    BAX=3
};
}
using namespace hermitian_gen_definite_eig_type_wrapper;

template<typename Real>
struct HermitianSdcCtrl {
    Int cutoff;
    Int maxInnerIts;
    Int maxOuterIts;
    Real tol;
    Real spreadFactor;
    bool progress;

    HermitianSdcCtrl()
    : cutoff(256), maxInnerIts(2), maxOuterIts(10),
      tol(0), spreadFactor(1e-6),
      progress(false)
    { }
};

template<typename Real>
struct HermitianEigCtrl
{
    HermitianTridiagCtrl tridiagCtrl;
    HermitianSdcCtrl<Real> sdcCtrl;
    bool useSdc;

    HermitianEigCtrl()
    : tridiagCtrl(), sdcCtrl(), useSdc(false)
    { }
};

// Compute the eigenvalues of a Hermitian matrix
template<typename F>
void HermitianEig
( UpperOrLower uplo, Matrix<F>& A, Matrix<BASE(F)>& w, SortType sort=UNSORTED,
  const HermitianEigCtrl<BASE(F)> ctrl=HermitianEigCtrl<BASE(F)>() );
template<typename F>
void HermitianEig
( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A, 
  DistMatrix<BASE(F),STAR,STAR>& w, SortType sort=UNSORTED,
  const HermitianEigCtrl<BASE(F)> ctrl=HermitianEigCtrl<BASE(F)>() );
template<typename F>
void HermitianEig
( UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<BASE(F),VR,STAR>& w, 
  SortType sort=UNSORTED,
  const HermitianEigCtrl<BASE(F)> ctrl=HermitianEigCtrl<BASE(F)>() );

// Compute the full eigenvalue decomposition of a Hermitian matrix
template<typename F>
void HermitianEig
( UpperOrLower uplo, Matrix<F>& A, Matrix<BASE(F)>& w, Matrix<F>& Z, 
  SortType sort=UNSORTED,
  const HermitianEigCtrl<BASE(F)> ctrl=HermitianEigCtrl<BASE(F)>() );
template<typename F>
void HermitianEig
( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A, 
  DistMatrix<BASE(F),STAR,STAR>& w, DistMatrix<F,STAR,STAR>& Z, 
  SortType sort=UNSORTED,
  const HermitianEigCtrl<BASE(F)> ctrl=HermitianEigCtrl<BASE(F)>() );
template<typename F>
void HermitianEig
( UpperOrLower uplo, 
  DistMatrix<F>& A, DistMatrix<BASE(F),VR,STAR>& w, DistMatrix<F>& paddedZ,
  SortType sort=UNSORTED,
  const HermitianEigCtrl<BASE(F)> ctrl=HermitianEigCtrl<BASE(F)>() );

// Compute the eigenvalues of a Hermitian matrix within a selected range
template<typename F>
void HermitianEig
( UpperOrLower uplo, Matrix<F>& A, Matrix<BASE(F)>& w,
  Int lowerBound, Int upperBound, SortType sort=UNSORTED,
  const HermitianEigCtrl<BASE(F)> ctrl=HermitianEigCtrl<BASE(F)>() );
template<typename F>
void HermitianEig
( UpperOrLower uplo, Matrix<F>& A, Matrix<BASE(F)>& w,
  BASE(F) lowerBound, BASE(F) upperBound, SortType sort=UNSORTED,
  const HermitianEigCtrl<BASE(F)> ctrl=HermitianEigCtrl<BASE(F)>() );
template<typename F>
void HermitianEig
( UpperOrLower uplo, 
  DistMatrix<F,STAR,STAR>& A, 
  DistMatrix<BASE(F),STAR,STAR>& w,
  Int lowerBound, Int upperBound, SortType sort=UNSORTED,
  const HermitianEigCtrl<BASE(F)> ctrl=HermitianEigCtrl<BASE(F)>() );
template<typename F>
void HermitianEig
( UpperOrLower uplo, 
  DistMatrix<F,STAR,STAR>& A, 
  DistMatrix<BASE(F),STAR,STAR>& w,
  BASE(F) lowerBound, BASE(F) upperBound, SortType sort=UNSORTED,
  const HermitianEigCtrl<BASE(F)> ctrl=HermitianEigCtrl<BASE(F)>() );
template<typename F>
void HermitianEig
( UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<BASE(F),VR,STAR>& w, 
  Int lowerBound, Int upperBound, SortType sort=UNSORTED,
  const HermitianEigCtrl<BASE(F)> ctrl=HermitianEigCtrl<BASE(F)>() );
template<typename F>
void HermitianEig
( UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<BASE(F),VR,STAR>& w,
  BASE(F) lowerBound, BASE(F) upperBound, SortType sort=UNSORTED,
  const HermitianEigCtrl<BASE(F)> ctrl=HermitianEigCtrl<BASE(F)>() );

// Compute a selected set of eigenpairs of a Hermitian matrix
template<typename F>
void HermitianEig
( UpperOrLower uplo, 
  Matrix<F>& A, Matrix<BASE(F)>& w, Matrix<F>& Z,
  Int lowerBound, Int upperBound, SortType sort=UNSORTED,
  const HermitianEigCtrl<BASE(F)> ctrl=HermitianEigCtrl<BASE(F)>() );
template<typename F>
void HermitianEig
( UpperOrLower uplo, 
  Matrix<F>& A, Matrix<BASE(F)>& w, Matrix<F>& Z,
  BASE(F) lowerBound, BASE(F) upperBound, SortType sort=UNSORTED,
  const HermitianEigCtrl<BASE(F)> ctrl=HermitianEigCtrl<BASE(F)>() );
template<typename F>
void HermitianEig
( UpperOrLower uplo,
  DistMatrix<F,STAR,STAR>& A, 
  DistMatrix<BASE(F),STAR,STAR>& w, 
  DistMatrix<F,STAR,STAR>& Z,
  Int lowerBound, Int upperBound, SortType sort=UNSORTED,
  const HermitianEigCtrl<BASE(F)> ctrl=HermitianEigCtrl<BASE(F)>() );
template<typename F>
void HermitianEig
( UpperOrLower uplo,
  DistMatrix<F,STAR,STAR>& A, 
  DistMatrix<BASE(F),STAR,STAR>& w, 
  DistMatrix<F,STAR,STAR>& Z,
  BASE(F) lowerBound, BASE(F) upperBound, SortType sort=UNSORTED,
  const HermitianEigCtrl<BASE(F)> ctrl=HermitianEigCtrl<BASE(F)>() );
template<typename F>
void HermitianEig
( UpperOrLower uplo, 
  DistMatrix<F>& A, DistMatrix<BASE(F),VR,STAR>& w, DistMatrix<F>& paddedZ,
  Int lowerBound, Int upperBound, SortType sort=UNSORTED,
  const HermitianEigCtrl<BASE(F)> ctrl=HermitianEigCtrl<BASE(F)>() );
template<typename F>
void HermitianEig
( UpperOrLower uplo, 
  DistMatrix<F>& A, DistMatrix<BASE(F),VR,STAR>& w, DistMatrix<F>& paddedZ,
  BASE(F) lowerBound, BASE(F) upperBound, SortType sort=UNSORTED,
  const HermitianEigCtrl<BASE(F)> ctrl=HermitianEigCtrl<BASE(F)>() );

} // namespace elem

#endif // ifndef ELEM_DECOMP_DECL_HPP
