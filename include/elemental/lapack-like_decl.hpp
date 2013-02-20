/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_DECL_HPP
#define LAPACK_DECL_HPP

namespace elem {

template<typename R>
void HermitianTridiag( UpperOrLower uplo, Matrix<R>& A );
template<typename R>
void HermitianTridiag( UpperOrLower uplo, DistMatrix<R>& A );
template<typename R>
void HermitianTridiag
( UpperOrLower uplo, Matrix<Complex<R> >& A, Matrix<Complex<R> >& t );
template<typename R>
void HermitianTridiag
( UpperOrLower uplo, 
  DistMatrix<Complex<R> >& A, DistMatrix<Complex<R>,STAR,STAR>& t );

void HermitianEig
( UpperOrLower uplo, DistMatrix<double>& A, 
  DistMatrix<double,VR,STAR>& w, DistMatrix<double>& paddedZ );
void HermitianEig
( UpperOrLower uplo, DistMatrix<double>& A, 
  DistMatrix<double,VR,STAR>& w, DistMatrix<double>& paddedZ,
  int lowerBound, int upperBound );
void HermitianEig
( UpperOrLower uplo, DistMatrix<double>& A, 
  DistMatrix<double,VR,STAR>& w, DistMatrix<double>& paddedZ,
  double lowerBound, double upperBound );
void HermitianEig
( UpperOrLower uplo, DistMatrix<double>& A,
  DistMatrix<double,VR,STAR>& w );
void HermitianEig
( UpperOrLower uplo, DistMatrix<double>& A,
  DistMatrix<double,VR,STAR>& w,
  int lowerBound, int upperBound );
void HermitianEig
( UpperOrLower uplo, DistMatrix<double>& A,
  DistMatrix<double,VR,STAR>& w,
  double lowerBound, double upperBound );

void HermitianEig
( UpperOrLower uplo, DistMatrix<Complex<double> >& A, 
  DistMatrix<double,VR,STAR>& w, DistMatrix<Complex<double> >& paddedZ );
void HermitianEig
( UpperOrLower uplo, DistMatrix<Complex<double> >& A, 
  DistMatrix<double,VR,STAR>& w, DistMatrix<Complex<double> >& paddedZ,
  int lowerBound, int upperBound );
void HermitianEig
( UpperOrLower uplo, DistMatrix<Complex<double> >& A, 
  DistMatrix<double,VR,STAR>& w, DistMatrix<Complex<double> >& paddedZ,
  double lowerBound, double upperBound );
void HermitianEig
( UpperOrLower uplo, DistMatrix<Complex<double> >& A,
  DistMatrix<double,VR,STAR>& w );
void HermitianEig
( UpperOrLower uplo, DistMatrix<Complex<double> >& A,
  DistMatrix<double,VR,STAR>& w,
  int lowerBound, int upperBound );
void HermitianEig
( UpperOrLower uplo, DistMatrix<Complex<double> >& A,
  DistMatrix<double,VR,STAR>& w,
  double lowerBound, double upperBound );

//
// HermitianGenDefiniteEig (Hermitian Generalized-Definite Eigensolver) 
//
namespace hermitian_gen_definite_eig_type_wrapper {
enum HermitianGenDefiniteEigType
{
    AXBX=1,
    ABX=2,
    BAX=3
};
}
using namespace hermitian_gen_definite_eig_type_wrapper;

//----------------------------------------------------------------------------//
// Tuning parameters                                                          //
//----------------------------------------------------------------------------//

namespace hermitian_tridiag_approach_wrapper {
enum HermitianTridiagApproach
{
    HERMITIAN_TRIDIAG_NORMAL, // Keep the current grid
    HERMITIAN_TRIDIAG_SQUARE, // Drop to a square process grid
    HERMITIAN_TRIDIAG_DEFAULT // Square grid algorithm only if already square
};
}
using namespace hermitian_tridiag_approach_wrapper;

void SetHermitianTridiagApproach( HermitianTridiagApproach approach );
HermitianTridiagApproach GetHermitianTridiagApproach();

// If dropping down to a square grid, the two simplest approaches are to take 
// the first r^2 processes from the original grid (for an r x r grid) and to
// either order them column-major or row-major to form the square grid.
void SetHermitianTridiagGridOrder( GridOrder order );
GridOrder GetHermitianTridiagGridOrder();

} // namespace elem

#endif // ifndef LAPACK_DECL_HPP
