/*
   Copyright (c) 2009-2011, Jack Poulson
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
#ifndef ELEMENTAL_ADVANCED_HPP
#define ELEMENTAL_ADVANCED_HPP 1

#include "elemental/basic.hpp"

namespace elemental {
namespace advanced {

//----------------------------------------------------------------------------//
// ApplyPackedReflectors                                                      //
//                                                                            //
// Applies the accumulated Householder transforms that are stored in the      //
// triangle of H specified by 'shape' to the matrix A.                        //
//                                                                            //
// If 'shape' is set to 'LOWER', then offset determines the diagonal that the //
// transforms are stored above (they are implicitly one on that diagonal).    //
//                                                                            //
// If 'shape' is set to 'UPPER', then offset determines the diagonal that the //
// transforms are stored below (they are implicitly one on that diagonal).    //
//                                                                            //
// 'direction' determines whether the reflectors are stored vertically or     //
// horizontally.                                                              //
//                                                                            //
// 'conjugation' determines whether or not the Householder scalars should be  //
// conjugated.                                                                //
//                                                                            //
// If 'order' is set to forward, then the reflectors are applied              //
// left-to-right, or top-to-bottom, depending on 'direction'. Otherwise, they //
// are applied in the opposite order.
//                                                                            //
// See the above note for QR factorizations regarding the vector 't' and      //
// Householder early-exit conditions.                                         //
//----------------------------------------------------------------------------//

// TODO: Serial versions

template<typename R>
void ApplyPackedReflectors
( Side side, Shape shape, VectorDirection direction, ForwardOrBackward order,
  int offset,
  const DistMatrix<R,MC,MR>& H, 
        DistMatrix<R,MC,MR>& A );

#ifndef WITHOUT_COMPLEX
template<typename R>
void ApplyPackedReflectors
( Side side, Shape shape,
  VectorDirection direction, ForwardOrBackward order, Conjugation conjugation,
  int offset,
  const DistMatrix<std::complex<R>,MC,  MR  >& H,
  const DistMatrix<std::complex<R>,MD,  STAR>& t,
        DistMatrix<std::complex<R>,MC,  MR  >& A );
template<typename R>
void ApplyPackedReflectors
( Side side, Shape shape, 
  VectorDirection direction, ForwardOrBackward order, Conjugation conjugation,
  int offset,
  const DistMatrix<std::complex<R>,MC,  MR  >& H,
  const DistMatrix<std::complex<R>,STAR,STAR>& t,
        DistMatrix<std::complex<R>,MC,  MR  >& A );
#endif

//----------------------------------------------------------------------------//
// ApplyRowPivots                                                             //
//                                                                            //
// Pivot the rows of the matrix A using the pivot vector (or the image and    //
// preimage of the associated permutation).                                   //
//                                                                            //
// SEE: ComposePivots                                                         //
//----------------------------------------------------------------------------//

// TODO: Serial versions

template<typename F>
void ApplyRowPivots
(       DistMatrix<F,  MC,MR  >& A,
  const DistMatrix<int,VC,STAR>& p );

template<typename F>
void ApplyRowPivots
(       DistMatrix<F,  MC,  MR  >& A,
  const DistMatrix<int,STAR,STAR>& p );

template<typename F>
void ApplyRowPivots
(       DistMatrix<F,MC,MR>& A,
  const std::vector<int>& image,
  const std::vector<int>& preimage );

//----------------------------------------------------------------------------//
// Cholesky:                                                                  //
//                                                                            //
// Overwrite a triangle of A with the Cholesky factor of A. 'shape'           //
// determines whether it is the upper or lower triangle.                      //
//----------------------------------------------------------------------------//

template<typename F>
void Cholesky( Shape shape, Matrix<F>& A );

template<typename F>
void Cholesky( Shape shape, DistMatrix<F,MC,MR>& A );

//----------------------------------------------------------------------------//
// CholeskySolve:                                                             //
//                                                                            //
// Overwrites X := inv(A) X, where A is Hermitian positive-definite, using a  //
// Cholesky factorization of A.                                               //
//----------------------------------------------------------------------------//

// TODO: Serial version

template<typename F>
void CholeskySolve
( Shape shape, DistMatrix<F,MC,MR>& A, DistMatrix<F,MC,MR>& X );

//----------------------------------------------------------------------------//
// ComposePivots                                                              //
//                                                                            //
// Explicitly form the image and preimage of the permutation associated with  //
// the given pivot vector.                                                    //
//                                                                            //
// SEE: ApplyRowPivots                                                        //
//----------------------------------------------------------------------------//

void ComposePivots
( const Matrix<int>& p, 
        std::vector<int>& image, 
        std::vector<int>& preimage );

void ComposePivots
( const DistMatrix<int,VC,STAR>& p, 
        std::vector<int>& image, 
        std::vector<int>& preimage );

void ComposePivots
( const DistMatrix<int,STAR,STAR>& p, 
        std::vector<int>& image, 
        std::vector<int>& preimage );

//----------------------------------------------------------------------------//
// Determinant:                                                               //
//                                                                            //
// Return the determinant of the matrix A.                                    //
//                                                                            //
// SafeDeterminant:                                                           //
//                                                                            //
// Return (rho,kappa,n) such that det(A) = rho exp(kappa n).                  //
// This decomposition of the determinant is done in order to reduce the       //
// possibility of (under/over)flow.                                           //
//----------------------------------------------------------------------------//

template<typename F>
F Determinant( Matrix<F>& A );
template<typename F>
SafeProduct<F> SafeDeterminant( Matrix<F>& A );

template<typename F>
F Determinant( DistMatrix<F,MC,MR>& A );
template<typename F>
SafeProduct<F> SafeDeterminant( DistMatrix<F,MC,MR>& A );

// TODO
/*
template<typename F>
F HPDDeterminant( Shape shape, Matrix<F>& A );
template<typename F>
SafeProduct<F> SafeHPDDeterminant( Shape shape, Matrix<F>& A );

template<typename F>
F HPDDeterminant( Shape shape, DistMatrix<F,MC,MR>& A );
template<typename F>
SafeProduct<F> SafeHPDDeterminant( Shape shape, DistMatrix<F,MC,MR>& A );
*/

//----------------------------------------------------------------------------//
// GaussianElimination:                                                       //
//                                                                            //
// Uses an LU factorization with partial pivoting to overwrite B := A^-1 B    //
//----------------------------------------------------------------------------//

// TODO: Serial version

template<typename F>
void GaussianElimination( DistMatrix<F,MC,MR>& A, DistMatrix<F,MC,MR>& B );

//----------------------------------------------------------------------------//
// HermitianGenDefiniteEig (Hermitian Generalized-Definite Eigensolver)       //
//----------------------------------------------------------------------------//

#ifndef WITHOUT_PMRRR
// Grab the full set of eigenpairs of real symmetric A and SPD B
template<typename R>
void HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, Shape shape, 
  DistMatrix<R,MC,  MR>& A, 
  DistMatrix<R,MC,  MR>& B, 
  DistMatrix<R,VR,STAR>& w,
  DistMatrix<R,MC,  MR>& X );
// Grab a partial set of eigenpairs. 
// The partial set is determined by the inclusive zero-indexed range 
//   a,a+1,...,b    ; a >= 0, b < n  
// of the n eigenpairs sorted from smallest to largest eigenvalues.  
template<typename R>
void HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, Shape shape,
  DistMatrix<R,MC,  MR>& A,
  DistMatrix<R,MC,  MR>& B,
  DistMatrix<R,VR,STAR>& w,
  DistMatrix<R,MC,  MR>& X,
  int a, int b );
// Grab a partial set of eigenpairs.
// The partial set is determined by the half-open interval (a,b]
template<typename R>
void HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, Shape shape,
  DistMatrix<R,MC,  MR>& A,
  DistMatrix<R,MC,  MR>& B,
  DistMatrix<R,VR,STAR>& w,
  DistMatrix<R,MC,  MR>& X,
  R a, R b );
// Grab the full set of eigenvalues of real symmetric A and SPD B
template<typename R>
void HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, Shape shape, 
  DistMatrix<R,MC,  MR>& A, 
  DistMatrix<R,MC,  MR>& B, 
  DistMatrix<R,VR,STAR>& w );
// Grab a partial set of eigenvalues. 
// The partial set is determined by the inclusive zero-indexed range 
//   a,a+1,...,b    ; a >= 0, b < n  
// of the n eigenpairs sorted from smallest to largest eigenvalues.  
template<typename R>
void HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, Shape shape,
  DistMatrix<R,MC,  MR>& A,
  DistMatrix<R,MC,  MR>& B,
  DistMatrix<R,VR,STAR>& w,
  int a, int b );
// Grab a partial set of eigenvalues.
// The partial set is determined by the half-open interval (a,b]
template<typename R>
void HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, Shape shape,
  DistMatrix<R,MC,  MR>& A,
  DistMatrix<R,MC,  MR>& B,
  DistMatrix<R,VR,STAR>& w,
  R a, R b );

#ifndef WITHOUT_COMPLEX
// Grab the full set of eigenpairs of complex Hermitian A and HPD B
template<typename R>
void HermitianGenDefiniteEig    
( HermitianGenDefiniteEigType type, Shape shape,
  DistMatrix<std::complex<R>,MC,  MR>& A,
  DistMatrix<std::complex<R>,MC,  MR>& B,
  DistMatrix<             R, VR,STAR>& w,
  DistMatrix<std::complex<R>,MC,  MR>& X );
// Grab a partial set of eigenpairs. 
// The partial set is determined by the inclusive zero-indexed range 
//   a,a+1,...,b    ; a >= 0, b < n  
// of the n eigenpairs sorted from smallest to largest eigenvalues.  
template<typename R>
void HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, Shape shape,
  DistMatrix<std::complex<R>,MC,  MR>& A,
  DistMatrix<std::complex<R>,MC,  MR>& B,
  DistMatrix<             R, VR,STAR>& w,
  DistMatrix<std::complex<R>,MC,  MR>& X,
  int a, int b );
// Grab a partial set of eigenpairs.
// The partial set is determined by the half-open interval (a,b]
template<typename R>
void HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, Shape shape,
  DistMatrix<std::complex<R>,MC,  MR>& A,
  DistMatrix<std::complex<R>,MC,  MR>& B,
  DistMatrix<             R, VR,STAR>& w,
  DistMatrix<std::complex<R>,MC,  MR>& X,
  R a, R b );
// Grab the full set of eigenvalues of complex Hermitian A and HPD B
template<typename R>
void HermitianGenDefiniteEig    
( HermitianGenDefiniteEigType type, Shape shape,
  DistMatrix<std::complex<R>,MC,  MR>& A,
  DistMatrix<std::complex<R>,MC,  MR>& B,
  DistMatrix<             R, VR,STAR>& w );
// Grab a partial set of eigenvalues. 
// The partial set is determined by the inclusive zero-indexed range 
//   a,a+1,...,b    ; a >= 0, b < n  
// of the n eigenpairs sorted from smallest to largest eigenvalues.  
template<typename R>
void HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, Shape shape,
  DistMatrix<std::complex<R>,MC,  MR>& A,
  DistMatrix<std::complex<R>,MC,  MR>& B,
  DistMatrix<             R, VR,STAR>& w,
  int a, int b );
// Grab a partial set of eigenvalues.
// The partial set is determined by the half-open interval (a,b]
template<typename R>
void HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, Shape shape,
  DistMatrix<std::complex<R>,MC,  MR>& A,
  DistMatrix<std::complex<R>,MC,  MR>& B,
  DistMatrix<             R, VR,STAR>& w,
  R a, R b );
#endif // WITHOUT_COMPLEX
#endif // WITHOUT_PMRRR

//----------------------------------------------------------------------------//
// Hegst (HErmitian GEneralized to STandard eigenvalue problem):              //
//                                                                            //
// If side==LEFT,                                                             //
//   reduce the problems                                                      //
//                      A B X = X Lambda to A X = X Lambda                    //
//                      B A X = X Lambda to A X = X Lambda                    //
// If side==RIGHT,                                                            //
//   reduce the problem A X = B X Lambda to A X = X Lambda                    //
//                                                                            //
// D contains the Cholesky factor of B in the triangle corresponding to the   //
// parameter 'shape'.                                                         //
//----------------------------------------------------------------------------//

template<typename F>
void Hegst
( Side side, Shape shape, 
  Matrix<F>& A, const Matrix<F>& B );

template<typename F>
void Hegst
( Side side, Shape shape, 
  DistMatrix<F,MC,MR>& A, const DistMatrix<F,MC,MR>& B );

//----------------------------------------------------------------------------//
// HermitianEig (Hermitian Eigensolver)                                       //
//----------------------------------------------------------------------------//

#ifndef WITHOUT_PMRRR
// Grab the full set of eigenpairs of the real, symmetric matrix A
void HermitianEig
( Shape shape, 
  DistMatrix<double,MC,  MR>& A, 
  DistMatrix<double,VR,STAR>& w,
  DistMatrix<double,MC,  MR>& Z );
// Grab a partial set of eigenpairs of the real, symmetric n x n matrix A. 
// The partial set is determined by the inclusive zero-indexed range 
//   a,a+1,...,b    ; a >= 0, b < n  
// of the n eigenpairs sorted from smallest to largest eigenvalues.  
void HermitianEig
( Shape shape,
  DistMatrix<double,MC,  MR>& A,
  DistMatrix<double,VR,STAR>& w,
  DistMatrix<double,MC,  MR>& Z,
  int a, int b );
// Grab a partial set of eigenpairs of the real, symmetric n x n matrix A. 
// The partial set is determined by the half-open interval (a,b]
void HermitianEig
( Shape shape,
  DistMatrix<double,MC,  MR>& A,
  DistMatrix<double,VR,STAR>& w,
  DistMatrix<double,MC,  MR>& Z,
  double a, double b );
// Grab the full set of eigenvalues of the real, symmetric matrix A
void HermitianEig
( Shape shape,
  DistMatrix<double,MC,  MR>& A,
  DistMatrix<double,VR,STAR>& w );
// Grab a partial set of eigenvalues of the real, symmetric n x n matrix A. 
// The partial set is determined by the inclusive zero-indexed range 
//   a,a+1,...,b    ; a >= 0, b < n  
// of the n eigenpairs sorted from smallest to largest eigenvalues.  
void HermitianEig
( Shape shape,
  DistMatrix<double,MC,  MR>& A,
  DistMatrix<double,VR,STAR>& w,
  int a, int b );
// Grab a partial set of eigenvalues of the real, symmetric n x n matrix A. 
// The partial set is determined by the half-open interval (a,b]
void HermitianEig
( Shape shape,
  DistMatrix<double,MC,  MR>& A,
  DistMatrix<double,VR,STAR>& w,
  double a, double b );
#ifndef WITHOUT_COMPLEX
// Grab the full set of eigenpairs of the complex, Hermitian matrix A
void HermitianEig    
( Shape shape,
  DistMatrix<std::complex<double>,MC,  MR>& A,
  DistMatrix<             double, VR,STAR>& w,
  DistMatrix<std::complex<double>,MC,  MR>& Z );
// Grab a partial set of eigenpairs of the complex, Hermitian n x n matrix A. 
// The partial set is determined by the inclusive zero-indexed range 
//   a,a+1,...,b    ; a >= 0, b < n  
// of the n eigenpairs sorted from smallest to largest eigenvalues.  
void HermitianEig
( Shape shape,
  DistMatrix<std::complex<double>,MC,  MR>& A,
  DistMatrix<             double, VR,STAR>& w,
  DistMatrix<std::complex<double>,MC,  MR>& Z,
  int a, int b );
// Grab a partial set of eigenpairs of the complex, Hermitian n x n matrix A. 
// The partial set is determined by the half-open interval (a,b]
void HermitianEig
( Shape shape,
  DistMatrix<std::complex<double>,MC,  MR>& A,
  DistMatrix<             double, VR,STAR>& w,
  DistMatrix<std::complex<double>,MC,  MR>& Z,
  double a, double b );
// Grab the full set of eigenvalues of the complex, Hermitian matrix A
void HermitianEig
( Shape shape,
  DistMatrix<std::complex<double>,MC,  MR>& A,
  DistMatrix<             double, VR,STAR>& w );
// Grab a partial set of eigenvalues of the complex, Hermitian n x n matrix A. 
// The partial set is determined by the inclusive zero-indexed range 
//   a,a+1,...,b    ; a >= 0, b < n  
// of the n eigenpairs sorted from smallest to largest eigenvalues.  
void HermitianEig
( Shape shape,
  DistMatrix<std::complex<double>,MC,  MR>& A,
  DistMatrix<             double, VR,STAR>& w,
  int a, int b );
// Grab a partial set of eigenvalues of the complex, Hermitian n x n matrix A. 
// The partial set is determined by the half-open interval (a,b]
void HermitianEig
( Shape shape,
  DistMatrix<std::complex<double>,MC,  MR>& A,
  DistMatrix<             double, VR,STAR>& w,
  double a, double b );
#endif // WITHOUT_COMPLEX
#endif // WITHOUT_PMRRR

//----------------------------------------------------------------------------//
// HouseholderSolve:                                                          //
//                                                                            //
// Overwrite X with the solution of inv(A) X or inv(A)^H X, where A need not  //
// be square. NOTE: If the system is underdetermined, then X should be the    //
// size of the solution and the input RHS should be stored in the top m rows  //
// of X, if the underdetermined matrix is m x n.                              //
//----------------------------------------------------------------------------//

// TODO: Serial version

template<typename R>
void HouseholderSolve
( Orientation orientation, DistMatrix<R,MC,MR>& A, DistMatrix<R,MC,MR>& X );
template<typename R>
void HouseholderSolve
( Orientation orientation, 
  DistMatrix<std::complex<R>,MC,MR>& A,
  DistMatrix<std::complex<R>,MC,MR>& X );

//----------------------------------------------------------------------------//
// HPDInverse:                                                              //
//                                                                            //
// Inverts a Hermitian positive-definite matrix.                              //
//----------------------------------------------------------------------------//

// TODO: Serial version

template<typename F>
void HPDInverse( Shape shape, DistMatrix<F,MC,MR>& A );

//----------------------------------------------------------------------------//
// LDLH (LDL^H factorization):                                                //
//                                                                            //
// Overwrite the lower triangle of A with L and d with the diagonal entries   //
// of D, so that A = L D L^H.                                                 //
//                                                                            //
// Partial pivoting is not yet supported.                                     //
//----------------------------------------------------------------------------//

// NOTE: Currently unblocked
template<typename F>
void LDLH( Matrix<F>& A, Matrix<F>& d );

template<typename F>
void LDLH( DistMatrix<F,MC,MR>& A, DistMatrix<F,MC,STAR>& d );

//----------------------------------------------------------------------------//
// LDLT (LDL^T factorization):                                                //
//                                                                            //
// Overwrite the lower triangle of A with L and d with the diagonal entries   //
// of D, so that A = L D L^T.                                                 //
//                                                                            //
// Partial pivoting is not yet supported.                                     //
//----------------------------------------------------------------------------//

// NOTE: Currently unblocked
template<typename F>
void LDLT( Matrix<F>& A, Matrix<F>& d );

template<typename F>
void LDLT( DistMatrix<F,MC,MR>& A, DistMatrix<F,MC,STAR>& d );

//----------------------------------------------------------------------------//
// LU (LU factorization):                                                     //
//                                                                            //
// If a container for a pivot vector is passed in, then A is overwritten with //
// its LU factorization after partial pivoting: P A = L U.                    //
// P is compressed into the vector p by storing the location of the nonzero   //
// element of each row.                                                       //
//                                                                            //
// If pivot vector is given, then A is overwritten with L U. Note that this   //
// version should usually be avoided, as pivoting is usually required for     //
// stability.                                                                 //
//----------------------------------------------------------------------------//

template<typename F>
void LU( Matrix<F>& A );

template<typename F> 
void LU( Matrix<F>& A, Matrix<int>& p );

template<typename F>
void LU( DistMatrix<F,MC,MR>& A );

template<typename F>
void LU( DistMatrix<F,MC,MR>& A, DistMatrix<int,VC,STAR>& p );

//----------------------------------------------------------------------------//
// LQ (LQ factorization):                                                     //
//                                                                            //
// Essentially the adjoint of a QR factorization on the adjoint of the input  //
// matrix.                                                                    //
//----------------------------------------------------------------------------//

// TODO: Serial versions

template<typename R>
void LQ( DistMatrix<R,MC,MR>& A );

#ifndef WITHOUT_COMPLEX
template<typename R>
void LQ
( DistMatrix<std::complex<R>,MC,MR  >& A, 
  DistMatrix<std::complex<R>,MD,STAR>& t );
#endif

//----------------------------------------------------------------------------//
// Norm                                                                       //
//----------------------------------------------------------------------------//

template<typename R>
R Norm( const Matrix<R>& A, NormType type=FROBENIUS_NORM );

template<typename R>
R Norm( const DistMatrix<R,MC,MR>& A, NormType type=FROBENIUS_NORM );

#ifndef WITHOUT_COMPLEX
template<typename R>
R Norm( const Matrix<std::complex<R> >& A, NormType type=FROBENIUS_NORM );

template<typename R>
R Norm
( const DistMatrix<std::complex<R>,MC,MR>& A, NormType type=FROBENIUS_NORM );
#endif

//----------------------------------------------------------------------------//
// HermitianNorm                                                              //
//----------------------------------------------------------------------------//

template<typename R>
R HermitianNorm
( Shape shape, const Matrix<R>& A, 
  NormType type=FROBENIUS_NORM );

template<typename R>
R HermitianNorm
( Shape shape, const DistMatrix<R,MC,MR>& A, 
  NormType type=FROBENIUS_NORM );

#ifndef WITHOUT_COMPLEX
template<typename R>
R HermitianNorm
( Shape shape, const Matrix<std::complex<R> >& A, 
  NormType type=FROBENIUS_NORM );

template<typename R>
R HermitianNorm
( Shape shape, const DistMatrix<std::complex<R>,MC,MR>& A, 
  NormType type=FROBENIUS_NORM );
#endif

//----------------------------------------------------------------------------//
// SymmetricNorm                                                              //
//----------------------------------------------------------------------------//

template<typename R>
R SymmetricNorm
( Shape shape, const Matrix<R>& A, 
  NormType type=FROBENIUS_NORM );

template<typename R>
R SymmetricNorm
( Shape shape, const DistMatrix<R,MC,MR>& A, 
  NormType type=FROBENIUS_NORM );

#ifndef WITHOUT_COMPLEX
template<typename R>
R SymmetricNorm
( Shape shape, const Matrix<std::complex<R> >& A, 
  NormType type=FROBENIUS_NORM );

template<typename R>
R SymmetricNorm
( Shape shape, const DistMatrix<std::complex<R>,MC,MR>& A, 
  NormType type=FROBENIUS_NORM );
#endif

//----------------------------------------------------------------------------//
// QR (QR factorization):                                                     //
//                                                                            //
// Performs a Householder QR factorization that overwrites the upper triangle //
// of A with R and fills the lower triangle with the scaled Householder       //
// transforms used to generate Q (they are implicitly one on the diagonal of  //
// A). The scaling factors for the Householder transforms are stored in t.    //
//                                                                            //
// For the complex case, 't' holds the Householder reflection coefficients    //
// that define the Householder transformation                                 //
//     House(tau,u) = I - tau u u^H                                           //
//                                                                            //
// IMPORTANT NOTE: The LAPACK convention for early-exiting when computing the //
// Householder reflection for a vector a = [ alpha11, a12 ]^T, where          //
// || a12 ||_2 = 0 and Im( alpha11 ) = 0, is to set 'tau' to zero in the      //
// Householder reflector equation:                                            //
//                                                                            //
//   House(tau,u) = I - tau u u^H                                             //
//                                                                            //
// which is not a valid Householder reflection due to the requirement that    //
// u be normalizable. We thus take the approach of setting tau = 2 when       //
// || a12 ||_2 = 0 and Im( alpha11 ) = 0, so that                             //
//                                                                            //
//   House(2,u) a = (I - 2 | 1 | | 1 0 | ) | alpha11 | = | -alpha11 |         //
//                         | 0 |           |    0    |   |     0    |         //
//                                                                            //
// This allows for the computation of the triangular matrix in the Compact WY //
// transform / UT transform to be computed mainly with Level 3 BLAS.          //
//----------------------------------------------------------------------------//

// TODO: Serial versions

template<typename R>
void QR( DistMatrix<R,MC,MR>& A );

#ifndef WITHOUT_COMPLEX
template<typename R>
void QR
( DistMatrix<std::complex<R>,MC,MR  >& A, 
  DistMatrix<std::complex<R>,MD,STAR>& t );
#endif

//----------------------------------------------------------------------------//
// Reflector (Householder reflector):                                         //
//----------------------------------------------------------------------------//

template<typename R>
R Reflector( Matrix<R>& chi, Matrix<R>& x );

#ifndef WITHOUT_COMPLEX
template<typename R>
std::complex<R>
Reflector( Matrix<std::complex<R> >& chi, Matrix<std::complex<R> >& x );
#endif

template<typename F>
F Reflector( DistMatrix<F,MC,MR>& chi, DistMatrix<F,MC,MR>& x );

//----------------------------------------------------------------------------//
// SkewHermitianEig (Skew-Hermitian Eigensolver)                              //
//----------------------------------------------------------------------------//

#ifndef WITHOUT_COMPLEX
#ifndef WITHOUT_PMRRR
// Grab the full set of eigenpairs of the real, skew-symmetric matrix G
void SkewHermitianEig
( Shape shape, 
  DistMatrix<double,              MC,  MR>& G, 
  DistMatrix<double,              VR,STAR>& wImag,
  DistMatrix<std::complex<double>,MC,  MR>& Z );
// Grab a partial set of eigenpairs of the real, skew-symmetric n x n matrix G. 
// The partial set is determined by the inclusive zero-indexed range 
//   a,a+1,...,b    ; a >= 0, b < n  
// of the n eigenpairs sorted from smallest to largest eigenvalues.  
void SkewHermitianEig
( Shape shape,
  DistMatrix<double,              MC,  MR>& G,
  DistMatrix<double,              VR,STAR>& wImag,
  DistMatrix<std::complex<double>,MC,  MR>& Z,
  int a, int b );
// Grab a partial set of eigenpairs of the real, skew-symmetric n x n matrix G. 
// The partial set is determined by the half-open imaginary interval (a,b]
void SkewHermitianEig
( Shape shape,
  DistMatrix<double,              MC,  MR>& G,
  DistMatrix<double,              VR,STAR>& wImag,
  DistMatrix<std::complex<double>,MC,  MR>& Z,
  double a, double b );
// Grab the full set of eigenvalues of the real, skew-symmetric matrix G 
void SkewHermitianEig
( Shape shape,
  DistMatrix<double,MC,  MR>& G,
  DistMatrix<double,VR,STAR>& wImag );
// Grab a partial set of eigenvalues of the real, skew-symmetric n x n matrix G.
// The partial set is determined by the inclusive zero-indexed range 
//   a,a+1,...,b    ; a >= 0, b < n  
// of the n eigenpairs sorted from smallest to largest eigenvalues.  
void SkewHermitianEig
( Shape shape,
  DistMatrix<double,MC,  MR>& G,
  DistMatrix<double,VR,STAR>& wImag,
  int a, int b );
// Grab a partial set of eigenvalues of the real, skew-symmetric n x n matrix G.
// The partial set is determined by the half-open imaginary interval (a,b]
void SkewHermitianEig
( Shape shape,
  DistMatrix<double,MC,  MR>& G,
  DistMatrix<double,VR,STAR>& wImag,
  double a, double b );

// Grab the full set of eigenpairs of the complex, skew-Hermitian matrix G 
void SkewHermitianEig    
( Shape shape,
  DistMatrix<std::complex<double>,MC,  MR>& G,
  DistMatrix<double,              VR,STAR>& wImag,
  DistMatrix<std::complex<double>,MC,  MR>& Z );
// Grab a partial set of eigenpairs of the complex, skew-Hermitian n x n matrix
// G. The partial set is determined by the inclusive zero-indexed range 
//   a,a+1,...,b    ; a >= 0, b < n  
// of the n eigenpairs sorted from smallest to largest eigenvalues.  
void SkewHermitianEig
( Shape shape,
  DistMatrix<std::complex<double>,MC,  MR>& G,
  DistMatrix<double,              VR,STAR>& wImag,
  DistMatrix<std::complex<double>,MC,  MR>& Z,
  int a, int b );
// Grab a partial set of eigenpairs of the complex, skew-Hermitian n x n matrix
// G. The partial set is determined by the half-open interval (a,b]
void SkewHermitianEig
( Shape shape,
  DistMatrix<std::complex<double>,MC,  MR>& G,
  DistMatrix<double,              VR,STAR>& wImag,
  DistMatrix<std::complex<double>,MC,  MR>& Z,
  double a, double b );
// Grab the full set of eigenvalues of the complex, skew-Hermitian matrix G 
void SkewHermitianEig
( Shape shape,
  DistMatrix<std::complex<double>,MC,  MR>& G,
  DistMatrix<double,              VR,STAR>& wImag );
// Grab a partial set of eigenvalues of the complex, skew-Hermitian n x n matrix
// G. The partial set is determined by the inclusive zero-indexed range 
//   a,a+1,...,b    ; a >= 0, b < n  
// of the n eigenpairs sorted from smallest to largest eigenvalues.  
void SkewHermitianEig
( Shape shape,
  DistMatrix<std::complex<double>,MC,  MR>& G,
  DistMatrix<double,              VR,STAR>& wImag,
  int a, int b );
// Grab a partial set of eigenvalues of the complex, skew-Hermitian n x n matrix
// G. The partial set is determined by the half-open imaginary interval (a,b]
void SkewHermitianEig
( Shape shape,
  DistMatrix<std::complex<double>,MC,  MR>& G,
  DistMatrix<double,              VR,STAR>& wImag,
  double a, double b );
#endif // WITHOUT_PMRRR
#endif // WITHOUT_COMPLEX

//----------------------------------------------------------------------------//
// SortEig                                                                    //
//----------------------------------------------------------------------------//

template<typename R>
void SortEig( DistMatrix<R,VR,STAR>& w );

template<typename R>
void SortEig( DistMatrix<R,VR,STAR>& w, DistMatrix<R,MC,MR>& Z );

template<typename R>
void SortEig( DistMatrix<R,VR,STAR>& w, DistMatrix<std::complex<R>,MC,MR>& Z );

//----------------------------------------------------------------------------//
// HermitianTridiag (Reduce Hermitian matrix to tridiagonal form):            //
//                                                                            //
// The diagonal and sub/super-diagonal of A are overwritten with a similar    //
// tridiagonal matrix that is found by successively applying Householder      //
// reflections to zero the matrix outside of the tridiagonal band.            //
//                                                                            //
// 'shape' decided which triangle of A specifies the Hermitian matrix, and on //
// exit the transforms are stored above the super/sub-diagonal and are        //
// implicitly one on the super/sub-diagonal.                                  //
//                                                                            //
// See the above note for QR factorizations detailing 't' and the difference  //
// in Householder transform early-exit approaches for the serial and parallel //
// routines.                                                                  //
//----------------------------------------------------------------------------//

// NOTE: Currently unblocked
template<typename R>
void HermitianTridiag( Shape shape, Matrix<R>& A );

#ifndef WITHOUT_COMPLEX
// NOTE: Currently unblocked
template<typename R>
void HermitianTridiag
( Shape shape, Matrix<std::complex<R> >& A, Matrix<std::complex<R> >& t );
#endif

template<typename R>
void HermitianTridiag( Shape shape, DistMatrix<R,MC,MR>& A );

#ifndef WITHOUT_COMPLEX
template<typename R>
void HermitianTridiag
( Shape shape,
  DistMatrix<std::complex<R>,MC,  MR  >& A,
  DistMatrix<std::complex<R>,STAR,STAR>& t );
#endif

void SetHermitianTridiagApproach( HermitianTridiagApproach approach );
HermitianTridiagApproach GetHermitianTridiagApproach();

// If dropping down to a square grid, the two simplest approaches are to take 
// the first r^2 processes from the original grid (for an r x r grid) and to
// either order them column-major or row-major to form the square grid.
void SetHermitianTridiagGridOrder( GridOrder order );
GridOrder GetHermitianTridiagGridOrder();

//----------------------------------------------------------------------------//
// Trace                                                                      //
//                                                                            //
// Returns the sum of the diagonal entries of a square matrix.                //
//----------------------------------------------------------------------------//

template<typename F>
F Trace( const Matrix<F>& A );

template<typename F>
F Trace( const DistMatrix<F,MC,MR>& A );

//----------------------------------------------------------------------------//
// TriangularInverse                                                          //
//                                                                            //
// Inverts a triangular matrix. 'shape' determines whether A is assumed to be //
// upper or lower triangular, and 'diagonal' determines whether or not A is   //
// to be treated as having a unit diagonal.                                   //
//----------------------------------------------------------------------------//

template<typename F>
void TriangularInverse
( Shape shape, Diagonal diagonal, Matrix<F>& A );

template<typename F>
void TriangularInverse
( Shape shape, Diagonal diagonal, DistMatrix<F,MC,MR>& A  );

} // advanced
} // elemental

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

#include "./advanced/internal.hpp"
#include "./advanced/ApplyPackedReflectors.hpp"
#include "./advanced/Cholesky.hpp"
#include "./advanced/CholeskySolve.hpp"
#include "./advanced/Determinant.hpp"
#include "./advanced/GaussianElimination.hpp"
#include "./advanced/Hegst.hpp"
#include "./advanced/HermitianEig.hpp"
#include "./advanced/HermitianGenDefiniteEig.hpp"
#include "./advanced/HermitianNorm.hpp"
#include "./advanced/HermitianTridiag.hpp"
#include "./advanced/HouseholderSolve.hpp"
#include "./advanced/HPDInverse.hpp"
#include "./advanced/LDL.hpp"
#include "./advanced/LQ.hpp"
#include "./advanced/LU.hpp"
#include "./advanced/Norm.hpp"
#include "./advanced/QR.hpp"
#include "./advanced/Reflector.hpp"
#include "./advanced/SkewHermitianEig.hpp"
#include "./advanced/SortEig.hpp"
#include "./advanced/Trace.hpp"
#include "./advanced/TriangularInverse.hpp"

template<typename F>
inline void
elemental::advanced::Cholesky
( Shape shape, Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("advanced::Cholesky");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
#endif
    const char uplo = ShapeToChar( shape );
    lapack::Cholesky( uplo, A.Height(), A.Buffer(), A.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
elemental::advanced::Hegst
( Side side, Shape shape, Matrix<F>& A, const Matrix<F>& B )
{
#ifndef RELEASE
    PushCallStack("advanced::Hegst");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    if( B.Height() != B.Width() )
        throw std::logic_error("B must be square");
    if( A.Height() != B.Height() )
        throw std::logic_error("A and B must be the same size");
#endif
    const int itype = ( side==LEFT ? 2 : 1 );
    const char uplo = ShapeToChar( shape );
    lapack::Hegst
    ( itype, uplo, A.Height(), 
      A.Buffer(), A.LDim(), B.LockedBuffer(), B.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
elemental::advanced::LU
( Matrix<F>& A, Matrix<int>& p )
{
#ifndef RELEASE
    PushCallStack("advanced::LU");
    if( p.Height() != A.Height() )
        throw std::logic_error("A and p must be the same height");
#endif
    lapack::LU
    ( A.Height(), A.Width(), A.Buffer(), A.LDim(), p.Buffer() );

    // Convert from Fortran to C indexing
    int* pBuffer = p.Buffer();
    const int n = A.Height();
    for( int i=0; i<n; ++i )
        --pBuffer[i];
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline R
elemental::advanced::SymmetricNorm
( Shape shape, const Matrix<R>& A, NormType type )
{ 
#ifndef RELEASE
    PushCallStack("advanced::SymmetricNorm");
#endif
    HermitianNorm( shape, A, type );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline R
elemental::advanced::SymmetricNorm
( Shape shape, const DistMatrix<R,MC,MR>& A, NormType type )
{ 
#ifndef RELEASE
    PushCallStack("advanced::SymmetricNorm");
#endif
    HermitianNorm( shape, A, type );
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
template<typename R>
inline R
elemental::advanced::SymmetricNorm
( Shape shape, const Matrix<std::complex<R> >& A, NormType type )
{ 
#ifndef RELEASE
    PushCallStack("advanced::SymmetricNorm");
#endif
    HermitianNorm( shape, A, type );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline R
elemental::advanced::SymmetricNorm
( Shape shape, const DistMatrix<std::complex<R>,MC,MR>& A, NormType type )
{ 
#ifndef RELEASE
    PushCallStack("advanced::SymmetricNorm");
#endif
    HermitianNorm( shape, A, type );
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX

template<typename F>
inline void
elemental::advanced::TriangularInverse
( Shape shape, Diagonal diagonal, Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("advanced::TriangularInverse");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
#endif
    const char uplo = ShapeToChar( shape );
    const char diag = DiagonalToChar( diagonal );
    lapack::TriangularInverse( uplo, diag, A.Height(), A.Buffer(), A.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

#endif /* ELEMENTAL_ADVANCED_HPP */

