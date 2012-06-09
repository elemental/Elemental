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
#ifndef ELEMENTAL_LAPACKLIKE_HPP
#define ELEMENTAL_LAPACKLIKE_HPP 1

namespace elem {

// TODO: Better organize this file

namespace norm_type_wrapper {
enum NormType
{
    FROBENIUS_NORM, // Compute the "vector" L_2 norm of the matrix
    INFINITY_NORM,  // Compute the L_oo norm of the matrix
    MAX_NORM,       // Compute the "vector" L_oo norm of the matrix
    ONE_NORM        // Compute the L_1 norm of the matrix
};
}
using namespace norm_type_wrapper;

namespace hermitian_gen_definite_eig_type_wrapper {
enum HermitianGenDefiniteEigType
{
    AXBX=1,
    ABX=2,
    BAX=3
};
}
using namespace hermitian_gen_definite_eig_type_wrapper;

namespace hermitian_tridiag_approach_wrapper {
enum HermitianTridiagApproach
{
    HERMITIAN_TRIDIAG_NORMAL, // Keep the current grid
    HERMITIAN_TRIDIAG_SQUARE, // Drop to a square process grid
    HERMITIAN_TRIDIAG_DEFAULT // Square grid algorithm only if already square
};
}
using namespace hermitian_tridiag_approach_wrapper;

//----------------------------------------------------------------------------//
// ApplyPackedReflectors                                                      //
//                                                                            //
// Applies the accumulated Householder transforms that are stored in the      //
// triangle of H specified by 'uplo' to the matrix A.                         //
//                                                                            //
// If 'uplo' is set to 'LOWER', then offset determines the diagonal that the  //
// transforms are stored above (they are implicitly one on that diagonal).    //
//                                                                            //
// If 'uplo' is set to 'UPPER', then offset determines the diagonal that the  //
// transforms are stored below (they are implicitly one on that diagonal).    //
//                                                                            //
// 'dir' determines whether the reflectors are stored vertically or           //
// horizontally.                                                              //
//                                                                            //
// 'conjugation' determines whether or not the Householder scalars should be  //
// conjugated.                                                                //
//                                                                            //
// If 'order' is set to forward, then the reflectors are applied              //
// left-to-right, or top-to-bottom, depending on 'dir'. Otherwise, they       //
// are applied in the opposite order.                                         //
//                                                                            //
// See the above note for QR factorizations regarding the vector 't' and      //
// Householder early-exit conditions.                                         //
//----------------------------------------------------------------------------//

// Serial versions
template<typename R>
void ApplyPackedReflectors
( LeftOrRight side, UpperOrLower uplo, 
  VerticalOrHorizontal dir, ForwardOrBackward order,
  int offset, const Matrix<R>& H, Matrix<R>& A );
template<typename R>
void ApplyPackedReflectors
( LeftOrRight side, UpperOrLower uplo,
  VerticalOrHorizontal dir, ForwardOrBackward order, Conjugation conjugation,
  int offset,
  const Matrix<Complex<R> >& H,
  const Matrix<Complex<R> >& t,
        Matrix<Complex<R> >& A );

// Parallel versions
template<typename R>
void ApplyPackedReflectors
( LeftOrRight side, UpperOrLower uplo, 
  VerticalOrHorizontal dir, ForwardOrBackward order,
  int offset, const DistMatrix<R>& H, DistMatrix<R>& A );
template<typename R>
void ApplyPackedReflectors
( LeftOrRight side, UpperOrLower uplo,
  VerticalOrHorizontal dir, ForwardOrBackward order, Conjugation conjugation,
  int offset,
  const DistMatrix<Complex<R> >& H,
  const DistMatrix<Complex<R>,MD,STAR>& t,
        DistMatrix<Complex<R> >& A );
template<typename R>
void ApplyPackedReflectors
( LeftOrRight side, UpperOrLower uplo, 
  VerticalOrHorizontal dir, ForwardOrBackward order, Conjugation conjugation,
  int offset,
  const DistMatrix<Complex<R> >& H,
  const DistMatrix<Complex<R>,STAR,STAR>& t,
        DistMatrix<Complex<R> >& A );

//----------------------------------------------------------------------------//
// Apply[Inverse]ColumnPivots/Apply[Inverse]RowPivots                         //
//                                                                            //
// Pivot the columns/rows of the matrix A using the pivot vector (or the      //
// image and preimage of the associated permutation).                         //
//                                                                            //
// SEE: ComposePivots                                                         //
//----------------------------------------------------------------------------//

template<typename F>
void ApplyColumnPivots( Matrix<F>& A, const Matrix<int>& p );
template<typename F>
void ApplyInverseColumnPivots( Matrix<F>& A, const Matrix<int>& p );

template<typename F>
void ApplyRowPivots( Matrix<F>& A, const Matrix<int>& p );
template<typename F>
void ApplyInverseRowPivots( Matrix<F>& A, const Matrix<int>& p );

template<typename F>
void ApplyColumnPivots
(       DistMatrix<F>& A,
  const DistMatrix<int,VC,STAR>& p );
template<typename F>
void ApplyInverseColumnPivots
(       DistMatrix<F>& A,
  const DistMatrix<int,VC,STAR>& p );

template<typename F>
void ApplyRowPivots
(       DistMatrix<F>& A,
  const DistMatrix<int,VC,STAR>& p );
template<typename F>
void ApplyInverseRowPivots
(       DistMatrix<F>& A,
  const DistMatrix<int,VC,STAR>& p );

template<typename F>
void ApplyColumnPivots
(       DistMatrix<F>& A,
  const DistMatrix<int,STAR,STAR>& p );
template<typename F>
void ApplyInverseColumnPivots
(       DistMatrix<F>& A,
  const DistMatrix<int,STAR,STAR>& p );

template<typename F>
void ApplyRowPivots
(       DistMatrix<F>& A,
  const DistMatrix<int,STAR,STAR>& p );
template<typename F>
void ApplyInverseRowPivots
(       DistMatrix<F>& A,
  const DistMatrix<int,STAR,STAR>& p );

template<typename F>
void ApplyColumnPivots
(       DistMatrix<F>& A,
  const std::vector<int>& image,
  const std::vector<int>& preimage );
template<typename F>
void ApplyRowPivots
(       DistMatrix<F>& A,
  const std::vector<int>& image,
  const std::vector<int>& preimage );

//----------------------------------------------------------------------------//
// PivotParity                                                                //
//                                                                            //
// Returns true iff the permutation is odd.                                   //
//----------------------------------------------------------------------------//
bool PivotParity( const Matrix<int>& p, int pivotOffset=0 );
bool PivotParity( const DistMatrix<int,VC,STAR>& p, int pivotOffset=0 );

//----------------------------------------------------------------------------//
// Cholesky:                                                                  //
//                                                                            //
// Overwrite a triangle of A with the Cholesky factor of A. 'uplo'            //
// determines whether it is the upper or lower triangle.                      //
//----------------------------------------------------------------------------//

template<typename F>
void Cholesky( UpperOrLower uplo, Matrix<F>& A );

template<typename F>
void Cholesky( UpperOrLower uplo, DistMatrix<F>& A );

//----------------------------------------------------------------------------//
// CholeskySolve:                                                             //
//                                                                            //
// Overwrites B := inv(A) B, where A is Hermitian positive-definite, using a  //
// Cholesky factorization of A.                                               //
//----------------------------------------------------------------------------//

// TODO: Serial version

template<typename F>
void CholeskySolve( UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<F>& B );

//----------------------------------------------------------------------------//
// ComposePivots                                                              //
//                                                                            //
// Explicitly form the image and preimage of the permutation associated with  //
// the given pivot vector.                                                    //
//                                                                            //
// SEE: ApplyColumnPivots/ApplyRowPivots                                      //
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
F Determinant( DistMatrix<F>& A );
template<typename F>
SafeProduct<F> SafeDeterminant( DistMatrix<F>& A );

// TODO
/*
template<typename F>
F HPDDeterminant( UpperOrLower uplo, Matrix<F>& A );
template<typename F>
SafeProduct<F> SafeHPDDeterminant( UpperOrLower uplo, Matrix<F>& A );

template<typename F>
F HPDDeterminant( UpperOrLower uplo, DistMatrix<F>& A );
template<typename F>
SafeProduct<F> SafeHPDDeterminant( UpperOrLower uplo, DistMatrix<F>& A );
*/

//----------------------------------------------------------------------------//
// GaussianElimination:                                                       //
//                                                                            //
// Uses an LU factorization with partial pivoting to overwrite B := A^-1 B    //
//----------------------------------------------------------------------------//

// TODO: Serial version

template<typename F>
void GaussianElimination( DistMatrix<F>& A, DistMatrix<F>& B );

//----------------------------------------------------------------------------//
// HermitianGenDefiniteEig (Hermitian Generalized-Definite Eigensolver)       //
//----------------------------------------------------------------------------//

#ifndef WITHOUT_PMRRR
// Grab the full set of eigenpairs of real symmetric A and SPD B
template<typename R>
void HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, UpperOrLower uplo, 
  DistMatrix<R>& A, 
  DistMatrix<R>& B, 
  DistMatrix<R,VR,STAR>& w,
  DistMatrix<R>& X );
// Grab a partial set of eigenpairs. 
// The partial set is determined by the inclusive zero-indexed range 
//   a,a+1,...,b    ; a >= 0, b < n  
// of the n eigenpairs sorted from smallest to largest eigenvalues.  
template<typename R>
void HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, UpperOrLower uplo,
  DistMatrix<R>& A,
  DistMatrix<R>& B,
  DistMatrix<R,VR,STAR>& w,
  DistMatrix<R>& X,
  int a, int b );
// Grab a partial set of eigenpairs.
// The partial set is determined by the half-open interval (a,b]
template<typename R>
void HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, UpperOrLower uplo,
  DistMatrix<R>& A, DistMatrix<R>& B, DistMatrix<R,VR,STAR>& w,
  DistMatrix<R>& X, R a, R b );
// Grab the full set of eigenvalues of real symmetric A and SPD B
template<typename R>
void HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, UpperOrLower uplo, 
  DistMatrix<R>& A, DistMatrix<R>& B, DistMatrix<R,VR,STAR>& w );
// Grab a partial set of eigenvalues. 
// The partial set is determined by the inclusive zero-indexed range 
//   a,a+1,...,b    ; a >= 0, b < n  
// of the n eigenpairs sorted from smallest to largest eigenvalues.  
template<typename R>
void HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, UpperOrLower uplo,
  DistMatrix<R>& A, DistMatrix<R>& B, DistMatrix<R,VR,STAR>& w,
  int a, int b );
// Grab a partial set of eigenvalues.
// The partial set is determined by the half-open interval (a,b]
template<typename R>
void HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, UpperOrLower uplo,
  DistMatrix<R>& A, DistMatrix<R>& B, DistMatrix<R,VR,STAR>& w,
  R a, R b );
// Grab the full set of eigenpairs of complex Hermitian A and HPD B
template<typename R>
void HermitianGenDefiniteEig    
( HermitianGenDefiniteEigType type, UpperOrLower uplo,
  DistMatrix<Complex<R> >& A,
  DistMatrix<Complex<R> >& B,
  DistMatrix<R,VR,STAR>& w,
  DistMatrix<Complex<R> >& X );
// Grab a partial set of eigenpairs. 
// The partial set is determined by the inclusive zero-indexed range 
//   a,a+1,...,b    ; a >= 0, b < n  
// of the n eigenpairs sorted from smallest to largest eigenvalues.  
template<typename R>
void HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, UpperOrLower uplo,
  DistMatrix<Complex<R> >& A,
  DistMatrix<Complex<R> >& B,
  DistMatrix<R,VR,STAR>& w,
  DistMatrix<Complex<R> >& X,
  int a, int b );
// Grab a partial set of eigenpairs.
// The partial set is determined by the half-open interval (a,b]
template<typename R>
void HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, UpperOrLower uplo,
  DistMatrix<Complex<R> >& A,
  DistMatrix<Complex<R> >& B,
  DistMatrix<R,VR,STAR>& w,
  DistMatrix<Complex<R> >& X,
  R a, R b );
// Grab the full set of eigenvalues of complex Hermitian A and HPD B
template<typename R>
void HermitianGenDefiniteEig    
( HermitianGenDefiniteEigType type, UpperOrLower uplo,
  DistMatrix<Complex<R> >& A,
  DistMatrix<Complex<R> >& B,
  DistMatrix<R,VR,STAR>& w );
// Grab a partial set of eigenvalues. 
// The partial set is determined by the inclusive zero-indexed range 
//   a,a+1,...,b    ; a >= 0, b < n  
// of the n eigenpairs sorted from smallest to largest eigenvalues.  
template<typename R>
void HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, UpperOrLower uplo,
  DistMatrix<Complex<R> >& A,
  DistMatrix<Complex<R> >& B,
  DistMatrix<R,VR,STAR>& w,
  int a, int b );
// Grab a partial set of eigenvalues.
// The partial set is determined by the half-open interval (a,b]
template<typename R>
void HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, UpperOrLower uplo,
  DistMatrix<Complex<R> >& A,
  DistMatrix<Complex<R> >& B,
  DistMatrix<R,VR,STAR>& w,
  R a, R b );
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
// parameter 'uplo'.                                                          //
//----------------------------------------------------------------------------//

template<typename F>
void Hegst
( LeftOrRight side, UpperOrLower uplo, 
  Matrix<F>& A, const Matrix<F>& B );

template<typename F>
void Hegst
( LeftOrRight side, UpperOrLower uplo, 
  DistMatrix<F>& A, const DistMatrix<F>& B );

//----------------------------------------------------------------------------//
// HermitianEig (Hermitian Eigensolver)                                       //
//----------------------------------------------------------------------------//

#ifndef WITHOUT_PMRRR
// Grab the full set of eigenpairs of the real, symmetric matrix A
void HermitianEig
( UpperOrLower uplo, 
  DistMatrix<double>& A, 
  DistMatrix<double,VR,STAR>& w,
  DistMatrix<double>& Z );
// Grab a partial set of eigenpairs of the real, symmetric n x n matrix A. 
// The partial set is determined by the inclusive zero-indexed range 
//   a,a+1,...,b    ; a >= 0, b < n  
// of the n eigenpairs sorted from smallest to largest eigenvalues.  
void HermitianEig
( UpperOrLower uplo,
  DistMatrix<double>& A,
  DistMatrix<double,VR,STAR>& w,
  DistMatrix<double>& Z,
  int a, int b );
// Grab a partial set of eigenpairs of the real, symmetric n x n matrix A. 
// The partial set is determined by the half-open interval (a,b]
void HermitianEig
( UpperOrLower uplo,
  DistMatrix<double>& A,
  DistMatrix<double,VR,STAR>& w,
  DistMatrix<double>& Z,
  double a, double b );
// Grab the full set of eigenvalues of the real, symmetric matrix A
void HermitianEig
( UpperOrLower uplo,
  DistMatrix<double>& A,
  DistMatrix<double,VR,STAR>& w );
// Grab a partial set of eigenvalues of the real, symmetric n x n matrix A. 
// The partial set is determined by the inclusive zero-indexed range 
//   a,a+1,...,b    ; a >= 0, b < n  
// of the n eigenpairs sorted from smallest to largest eigenvalues.  
void HermitianEig
( UpperOrLower uplo,
  DistMatrix<double>& A,
  DistMatrix<double,VR,STAR>& w,
  int a, int b );
// Grab a partial set of eigenvalues of the real, symmetric n x n matrix A. 
// The partial set is determined by the half-open interval (a,b]
void HermitianEig
( UpperOrLower uplo,
  DistMatrix<double>& A,
  DistMatrix<double,VR,STAR>& w,
  double a, double b );
// Grab the full set of eigenpairs of the complex, Hermitian matrix A
void HermitianEig    
( UpperOrLower uplo,
  DistMatrix<Complex<double> >& A,
  DistMatrix<double,VR,STAR>& w,
  DistMatrix<Complex<double> >& Z );
// Grab a partial set of eigenpairs of the complex, Hermitian n x n matrix A. 
// The partial set is determined by the inclusive zero-indexed range 
//   a,a+1,...,b    ; a >= 0, b < n  
// of the n eigenpairs sorted from smallest to largest eigenvalues.  
void HermitianEig
( UpperOrLower uplo,
  DistMatrix<Complex<double> >& A,
  DistMatrix<double,VR,STAR>& w,
  DistMatrix<Complex<double> >& Z,
  int a, int b );
// Grab a partial set of eigenpairs of the complex, Hermitian n x n matrix A. 
// The partial set is determined by the half-open interval (a,b]
void HermitianEig
( UpperOrLower uplo,
  DistMatrix<Complex<double> >& A,
  DistMatrix<double,VR,STAR>& w,
  DistMatrix<Complex<double> >& Z,
  double a, double b );
// Grab the full set of eigenvalues of the complex, Hermitian matrix A
void HermitianEig
( UpperOrLower uplo,
  DistMatrix<Complex<double> >& A,
  DistMatrix<double,VR,STAR>& w );
// Grab a partial set of eigenvalues of the complex, Hermitian n x n matrix A. 
// The partial set is determined by the inclusive zero-indexed range 
//   a,a+1,...,b    ; a >= 0, b < n  
// of the n eigenpairs sorted from smallest to largest eigenvalues.  
void HermitianEig
( UpperOrLower uplo,
  DistMatrix<Complex<double> >& A,
  DistMatrix<double,VR,STAR>& w,
  int a, int b );
// Grab a partial set of eigenvalues of the complex, Hermitian n x n matrix A. 
// The partial set is determined by the half-open interval (a,b]
void HermitianEig
( UpperOrLower uplo,
  DistMatrix<Complex<double> >& A,
  DistMatrix<double,VR,STAR>& w,
  double a, double b );
#endif // WITHOUT_PMRRR

//----------------------------------------------------------------------------//
// [Real/Complex]HermitianFunction                                            //
//                                                                            //
// Uses a Hermitian eigenvalue decomposition and the passed in functor to     //
// modify the eigenvalues of the passed-in matrix.                            //
//                                                                            //
// When the functor's range is real, the resulting matrix will remain         //
// Hermitian, but when the functor's range is complex, the matrix is only     //
// normal (and hence the full matrix must be stored).                         //
//----------------------------------------------------------------------------//

#ifndef WITHOUT_PMRRR
template<typename F,class RealFunctor>
void RealHermitianFunction
( UpperOrLower uplo, DistMatrix<F>& A, const RealFunctor& f );

template<typename R,class ComplexFunctor>
void ComplexHermitianFunction
( UpperOrLower uplo, DistMatrix<Complex<R> >& A, const ComplexFunctor& f );
#endif // WITHOUT_PMRRR

//----------------------------------------------------------------------------//
// HermitianPseudoinverse                                                     //
//                                                                            //
// Specializes RealHermitianFunction routine to compute the pseudoinverse.    //
//----------------------------------------------------------------------------//

#ifndef WITHOUT_PMRRR
template<typename F>
void HermitianPseudoinverse( UpperOrLower uplo, DistMatrix<F>& A );
#endif // WITHOUT_PMRRR

//----------------------------------------------------------------------------//
// HermitianSVD                                                               //
//----------------------------------------------------------------------------//

#ifndef WITHOUT_PMRRR
template<typename F>
void HermitianSVD
( UpperOrLower uplo, 
  DistMatrix<F>& A, DistMatrix<typename Base<F>::type,VR,STAR>& s, 
  DistMatrix<F>& U, DistMatrix<F>& V );

template<typename F>
void HermitianSingularValues
( UpperOrLower uplo, 
  DistMatrix<F>& A, DistMatrix<typename Base<F>::type,VR,STAR>& s );
#endif // WITHOUT_PMRRR

//----------------------------------------------------------------------------//
// HouseholderSolve:                                                          //
//                                                                            //
// Form X as the minimum norm solution of                                     //
//    || op(A) X - B ||_2                                                     //
// using Householder based QR or LQ factorizations, where op(A) can be either //
// A or A^H, depending on the value or 'orientation'.                         //
//----------------------------------------------------------------------------//

// TODO: Serial version

template<typename R>
void HouseholderSolve
( Orientation orientation, 
  DistMatrix<R>& A, const DistMatrix<R>& B, DistMatrix<R>& X );
template<typename R>
void HouseholderSolve
( Orientation orientation, 
        DistMatrix<Complex<R> >& A,
  const DistMatrix<Complex<R> >& B, 
        DistMatrix<Complex<R> >& X );

//----------------------------------------------------------------------------//
// HPDInverse:                                                                //
//                                                                            //
// Inverts a Hermitian positive-definite matrix.                              //
//----------------------------------------------------------------------------//

// TODO: Serial version

template<typename F>
void HPDInverse( UpperOrLower uplo, DistMatrix<F>& A );

//----------------------------------------------------------------------------//
// HPSDCholesky                                                               //
//                                                                            //
// Compute the Cholesky factors of a potentially singular HPSD matrix through //
// the QR or LQ factorization of its matrix square-root:                      //
//    If sqrtA*sqrtA = A, and QR=sqrtA, then sqrtA=sqrtA' implies             //
//               R' Q' Q R = R' R = A,                                        //
// and a similar argument yields L L' = A.                                    //
//----------------------------------------------------------------------------//

#ifndef WITHOUT_PMRRR
template<typename R>
void HPSDCholesky( UpperOrLower uplo, DistMatrix<R>& A );

template<typename R>
void HPSDCholesky( UpperOrLower uplo, DistMatrix<Complex<R> >& A );
#endif // WITHOUT_PMRRR

//----------------------------------------------------------------------------//
// Pseudoinverse                                                              //
//                                                                            //
// Uses a Singular Value Decomposition to form the pseudoinverse of A.        //
//----------------------------------------------------------------------------//

template<typename F>
void Pseudoinverse( DistMatrix<F>& A );

//----------------------------------------------------------------------------//
// SquareRoot                                                                 //
//                                                                            //
// Compute the square-root of a Hermitian positive semi-definite matrix       //
// through its eigenvalue decomposition:                                      //
// A = U Lambda U^H => sqrt(A) = U sqrt(Lambda) U^H                           //
//----------------------------------------------------------------------------//

#ifndef WITHOUT_PMRRR
template<typename F>
void HPSDSquareRoot( UpperOrLower uplo, DistMatrix<F>& A );
#endif // WITHOUT_PMRRR

//----------------------------------------------------------------------------//
// Inverse:                                                                   //
//                                                                            //
// Inverts a general fully-populated matrix.                                  //
//----------------------------------------------------------------------------//

template<typename F>
void Inverse( Matrix<F>& A );

template<typename F>
void Inverse( DistMatrix<F>& A );

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
void LDLH( Matrix<F>& A );
template<typename F>
void LDLH( Matrix<F>& A, Matrix<F>& d );

template<typename F>
void LDLH( DistMatrix<F>& A );
template<typename F>
void LDLH( DistMatrix<F>& A, DistMatrix<F,MC,STAR>& d );

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
void LDLT( Matrix<F>& A );
template<typename F>
void LDLT( Matrix<F>& A, Matrix<F>& d );

template<typename F>
void LDLT( DistMatrix<F>& A );
template<typename F>
void LDLT( DistMatrix<F>& A, DistMatrix<F,MC,STAR>& d );

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

//
// LU without pivoting
//

template<typename F>
void LU( Matrix<F>& A );

template<typename F>
void LU( DistMatrix<F>& A );

//
// LU with partial pivoting
//

template<typename F> 
void LU( Matrix<F>& A, Matrix<int>& p );

template<typename F>
void LU( DistMatrix<F>& A, DistMatrix<int,VC,STAR>& p );

//----------------------------------------------------------------------------//
// LQ (LQ factorization):                                                     //
//                                                                            //
// Essentially the adjoint of a QR factorization on the adjoint of the input  //
// matrix.                                                                    //
//----------------------------------------------------------------------------//

// Implementations which overwrite A with both L and Q (with Q implicitly 
// represented as the product of Householder reflectors).
template<typename Real>
void LQ( Matrix<Real>& A );
template<typename Real>
void LQ( DistMatrix<Real>& A );
template<typename Real>
void LQ( Matrix<Complex<Real> >& A, 
         Matrix<Complex<Real> >& t );
template<typename Real>
void LQ( DistMatrix<Complex<Real> >& A, 
         DistMatrix<Complex<Real>,MD,STAR>& t );

// Implementations which return the explicit Q from the LQ decomposition
template<typename F>
void ExplicitLQ( Matrix<F>& A );
template<typename F>
void ExplicitLQ( DistMatrix<F>& A );

// Implementations which return the explicit LQ decomposition, with 
// Q overwriting A on exit
template<typename F>
void ExplicitLQ( Matrix<F>& L, Matrix<F>& A );
template<typename F>
void ExplicitLQ( DistMatrix<F>& L, DistMatrix<F>& A );

//----------------------------------------------------------------------------//
// Norm                                                                       //
//----------------------------------------------------------------------------//

template<typename F>
typename Base<F>::type 
Norm( const Matrix<F>& A, NormType type=FROBENIUS_NORM );
template<typename F>
typename Base<F>::type 
Norm( const DistMatrix<F>& A, NormType type=FROBENIUS_NORM );

// TODO: provide an option to compute more accurate estimates
template<typename F>
typename Base<F>::type TwoNormLowerBound( const Matrix<F>& A );
template<typename F>
typename Base<F>::type TwoNormUpperBound( const Matrix<F>& A );
template<typename F>
typename Base<F>::type TwoNormLowerBound( const Matrix<F>& A );
template<typename F>
typename Base<F>::type TwoNormUpperBound( const Matrix<F>& A );

//----------------------------------------------------------------------------//
// HermitianNorm                                                              //
//----------------------------------------------------------------------------//

template<typename F>
typename Base<F>::type 
HermitianNorm
( UpperOrLower uplo, const Matrix<F>& A, NormType type=FROBENIUS_NORM );

template<typename F>
typename Base<F>::type
HermitianNorm
( UpperOrLower uplo, const DistMatrix<F>& A, NormType type=FROBENIUS_NORM );

//----------------------------------------------------------------------------//
// SymmetricNorm                                                              //
//----------------------------------------------------------------------------//

template<typename F>
typename Base<F>::type
SymmetricNorm
( UpperOrLower uplo, const Matrix<F>& A, NormType type=FROBENIUS_NORM );

template<typename F>
typename Base<F>::type 
SymmetricNorm
( UpperOrLower uplo, const DistMatrix<F>& A, NormType type=FROBENIUS_NORM );

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

// Implementations which overwrite A with both Q and R (with Q implicitly 
// represented as the product of Householder reflectors).
template<typename Real>
void QR( Matrix<Real>& A );
template<typename Real>
void QR( DistMatrix<Real>& A );
template<typename Real>
void QR( Matrix<Complex<Real> >& A, 
         Matrix<Complex<Real> >& t );
template<typename Real>
void QR( DistMatrix<Complex<Real> >& A, 
         DistMatrix<Complex<Real>,MD,STAR>& t );

// Implementations which return the explicit Q from the QR decomposition
template<typename F>
void ExplicitQR( Matrix<F>& A );
template<typename F>
void ExplicitQR( DistMatrix<F>& A );

// Implementations which return the explicit QR decomposition, with 
// Q overwriting A on exit
template<typename F>
void ExplicitQR( Matrix<F>& A, Matrix<F>& R );
template<typename F>
void ExplicitQR( DistMatrix<F>& A, DistMatrix<F>& R );

//----------------------------------------------------------------------------//
// Reflector (Householder reflector):                                         //
//----------------------------------------------------------------------------//

template<typename R>
R Reflector( Matrix<R>& chi, Matrix<R>& x );

template<typename R>
Complex<R>
Reflector( Matrix<Complex<R> >& chi, Matrix<Complex<R> >& x );

template<typename F>
F Reflector( DistMatrix<F>& chi, DistMatrix<F>& x );

//----------------------------------------------------------------------------//
// SkewHermitianEig (Skew-Hermitian Eigensolver)                              //
//----------------------------------------------------------------------------//

#ifndef WITHOUT_PMRRR
// Grab the full set of eigenpairs of the real, skew-symmetric matrix G
void SkewHermitianEig
( UpperOrLower uplo, 
  DistMatrix<double>& G, 
  DistMatrix<double,VR,STAR>& wImag,
  DistMatrix<Complex<double> >& Z );
// Grab a partial set of eigenpairs of the real, skew-symmetric n x n matrix G. 
// The partial set is determined by the inclusive zero-indexed range 
//   a,a+1,...,b    ; a >= 0, b < n  
// of the n eigenpairs sorted from smallest to largest eigenvalues.  
void SkewHermitianEig
( UpperOrLower uplo,
  DistMatrix<double>& G,
  DistMatrix<double,VR,STAR>& wImag,
  DistMatrix<Complex<double> >& Z,
  int a, int b );
// Grab a partial set of eigenpairs of the real, skew-symmetric n x n matrix G. 
// The partial set is determined by the half-open imaginary interval (a,b]
void SkewHermitianEig
( UpperOrLower uplo,
  DistMatrix<double>& G,
  DistMatrix<double,VR,STAR>& wImag,
  DistMatrix<Complex<double> >& Z,
  double a, double b );
// Grab the full set of eigenvalues of the real, skew-symmetric matrix G 
void SkewHermitianEig
( UpperOrLower uplo,
  DistMatrix<double>& G,
  DistMatrix<double,VR,STAR>& wImag );
// Grab a partial set of eigenvalues of the real, skew-symmetric n x n matrix G.
// The partial set is determined by the inclusive zero-indexed range 
//   a,a+1,...,b    ; a >= 0, b < n  
// of the n eigenpairs sorted from smallest to largest eigenvalues.  
void SkewHermitianEig
( UpperOrLower uplo,
  DistMatrix<double>& G,
  DistMatrix<double,VR,STAR>& wImag,
  int a, int b );
// Grab a partial set of eigenvalues of the real, skew-symmetric n x n matrix G.
// The partial set is determined by the half-open imaginary interval (a,b]
void SkewHermitianEig
( UpperOrLower uplo,
  DistMatrix<double>& G,
  DistMatrix<double,VR,STAR>& wImag,
  double a, double b );

// Grab the full set of eigenpairs of the complex, skew-Hermitian matrix G 
void SkewHermitianEig    
( UpperOrLower uplo,
  DistMatrix<Complex<double> >& G,
  DistMatrix<double,VR,STAR>& wImag,
  DistMatrix<Complex<double> >& Z );
// Grab a partial set of eigenpairs of the complex, skew-Hermitian n x n matrix
// G. The partial set is determined by the inclusive zero-indexed range 
//   a,a+1,...,b    ; a >= 0, b < n  
// of the n eigenpairs sorted from smallest to largest eigenvalues.  
void SkewHermitianEig
( UpperOrLower uplo,
  DistMatrix<Complex<double> >& G,
  DistMatrix<double,VR,STAR>& wImag,
  DistMatrix<Complex<double> >& Z,
  int a, int b );
// Grab a partial set of eigenpairs of the complex, skew-Hermitian n x n matrix
// G. The partial set is determined by the half-open interval (a,b]
void SkewHermitianEig
( UpperOrLower uplo,
  DistMatrix<Complex<double> >& G,
  DistMatrix<double,VR,STAR>& wImag,
  DistMatrix<Complex<double> >& Z,
  double a, double b );
// Grab the full set of eigenvalues of the complex, skew-Hermitian matrix G 
void SkewHermitianEig
( UpperOrLower uplo,
  DistMatrix<Complex<double> >& G,
  DistMatrix<double,VR,STAR>& wImag );
// Grab a partial set of eigenvalues of the complex, skew-Hermitian n x n matrix
// G. The partial set is determined by the inclusive zero-indexed range 
//   a,a+1,...,b    ; a >= 0, b < n  
// of the n eigenpairs sorted from smallest to largest eigenvalues.  
void SkewHermitianEig
( UpperOrLower uplo,
  DistMatrix<Complex<double> >& G,
  DistMatrix<double,VR,STAR>& wImag,
  int a, int b );
// Grab a partial set of eigenvalues of the complex, skew-Hermitian n x n matrix
// G. The partial set is determined by the half-open imaginary interval (a,b]
void SkewHermitianEig
( UpperOrLower uplo,
  DistMatrix<Complex<double> >& G,
  DistMatrix<double,VR,STAR>& wImag,
  double a, double b );
#endif // WITHOUT_PMRRR

//----------------------------------------------------------------------------//
// SortEig                                                                    //
//----------------------------------------------------------------------------//

template<typename R>
void SortEig( DistMatrix<R,VR,STAR>& w );

template<typename R>
void SortEig( DistMatrix<R,VR,STAR>& w, DistMatrix<R>& Z );

template<typename R>
void SortEig( DistMatrix<R,VR,STAR>& w, DistMatrix<Complex<R> >& Z );

//----------------------------------------------------------------------------//
// Bidiag (Reduce general matrix to bidiagonal form):                         //
//                                                                            //
//----------------------------------------------------------------------------//

// NOTE: Currently unblocked
template<typename R>
void Bidiag( Matrix<R>& A );

// NOTE: Currently unblocked
template<typename R>
void Bidiag
( Matrix<Complex<R> >& A, 
  Matrix<Complex<R> >& tP, 
  Matrix<Complex<R> >& tQ );

// WARNING: Not yet finished
template<typename R>
void Bidiag( DistMatrix<R>& A );

// WARNING: Not yet finished
template<typename R>
void Bidiag
( DistMatrix<Complex<R> >& A,
  DistMatrix<Complex<R>,STAR,STAR>& tP,
  DistMatrix<Complex<R>,STAR,STAR>& tQ );

//----------------------------------------------------------------------------//
// HermitianTridiag (Reduce Hermitian matrix to tridiagonal form):            //
//                                                                            //
// The diagonal and sub/super-diagonal of A are overwritten with a similar    //
// tridiagonal matrix that is found by successively applying Householder      //
// reflections to zero the matrix outside of the tridiagonal band.            //
//                                                                            //
// 'uplo' decided which triangle of A specifies the Hermitian matrix, and on  //
// exit the transforms are stored above the super/sub-diagonal and are        //
// implicitly one on the super/sub-diagonal.                                  //
//                                                                            //
// See the above note for QR factorizations detailing 't' and the difference  //
// in Householder transform early-exit approaches for the serial and parallel //
// routines.                                                                  //
//----------------------------------------------------------------------------//

// NOTE: Currently unblocked
template<typename R>
void HermitianTridiag( UpperOrLower uplo, Matrix<R>& A );

// NOTE: Currently unblocked
template<typename R>
void HermitianTridiag
( UpperOrLower uplo, Matrix<Complex<R> >& A, Matrix<Complex<R> >& t );

template<typename R>
void HermitianTridiag( UpperOrLower uplo, DistMatrix<R>& A );

template<typename R>
void HermitianTridiag
( UpperOrLower uplo,
  DistMatrix<Complex<R> >& A,
  DistMatrix<Complex<R>,STAR,STAR>& t );

void SetHermitianTridiagApproach( HermitianTridiagApproach approach );
HermitianTridiagApproach GetHermitianTridiagApproach();

// If dropping down to a square grid, the two simplest approaches are to take 
// the first r^2 processes from the original grid (for an r x r grid) and to
// either order them column-major or row-major to form the square grid.
void SetHermitianTridiagGridOrder( GridOrder order );
GridOrder GetHermitianTridiagGridOrder();

//----------------------------------------------------------------------------//
// Polar decomposition                                                        //
//----------------------------------------------------------------------------//

template<typename F>
void Polar( DistMatrix<F>& A, DistMatrix<F>& P );

//----------------------------------------------------------------------------//
// SVD                                                                        //
//----------------------------------------------------------------------------//

template<typename F>
void SVD
( Matrix<F>& A, Matrix<typename Base<F>::type>& s, Matrix<F>& V, 
  bool useQR=false );
template<typename F>
void SVD
( DistMatrix<F>& A, 
  DistMatrix<typename Base<F>::type,VR,STAR>& s, 
  DistMatrix<F>& V );

template<typename F>
void SingularValues
( Matrix<F>& A, Matrix<typename Base<F>::type>& s );
template<typename F>
void SingularValues
( DistMatrix<F>& A,
  DistMatrix<typename Base<F>::type,VR,STAR>& s );

//----------------------------------------------------------------------------//
// Trace                                                                      //
//                                                                            //
// Returns the sum of the diagonal entries of a square matrix.                //
//----------------------------------------------------------------------------//

template<typename F>
F Trace( const Matrix<F>& A );

template<typename F>
F Trace( const DistMatrix<F>& A );

//----------------------------------------------------------------------------//
// TriangularInverse                                                          //
//                                                                            //
// Inverts a triangular matrix. 'uplo' determines whether A is assumed to be  //
// upper or lower triangular, and 'diag' determines whether or not A is       //
// to be treated as having a unit diagonal.                                   //
//----------------------------------------------------------------------------//

template<typename F>
void TriangularInverse
( UpperOrLower uplo, UnitOrNonUnit diag, Matrix<F>& A );

template<typename F>
void TriangularInverse
( UpperOrLower uplo, UnitOrNonUnit diag, DistMatrix<F>& A  );

//----------------------------------------------------------------------------//
// QR-based Halley iteration for the polar decomposition                      //
//----------------------------------------------------------------------------//

template<typename F>
int QDWH
( Matrix<F>& A, 
  typename Base<F>::type lowerBound, 
  typename Base<F>::type twoNormEstimate );
template<typename F>
int QDWH
( DistMatrix<F>& A, 
  typename Base<F>::type lowerBound, 
  typename Base<F>::type twoNormEstimate );

template<typename F>
int Halley( Matrix<F>& A, typename Base<F>::type twoNormEstimate );
template<typename F>
int Halley( DistMatrix<F>& A, typename Base<F>::type twoNormEstimate );

} // namespace elem

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

#include "./lapack-like/internal.hpp"
#include "./lapack-like/ApplyPackedReflectors.hpp"
#include "./lapack-like/ApplyColumnPivots.hpp"
#include "./lapack-like/ApplyRowPivots.hpp"
#include "./lapack-like/Bidiag.hpp"
#include "./lapack-like/Cholesky.hpp"
#include "./lapack-like/CholeskySolve.hpp"
#include "./lapack-like/ComposePivots.hpp"
#include "./lapack-like/Determinant.hpp"
#include "./lapack-like/ExplicitLQ.hpp"
#include "./lapack-like/ExplicitQR.hpp"
#include "./lapack-like/GaussianElimination.hpp"
#include "./lapack-like/Halley.hpp"
#include "./lapack-like/Hegst.hpp"
#include "./lapack-like/HermitianEig.hpp"
#include "./lapack-like/HermitianFunction.hpp"
#include "./lapack-like/HermitianGenDefiniteEig.hpp"
#include "./lapack-like/HermitianNorm.hpp"
#include "./lapack-like/HermitianPseudoinverse.hpp"
#include "./lapack-like/HermitianSVD.hpp"
#include "./lapack-like/HermitianTridiag.hpp"
#include "./lapack-like/HouseholderSolve.hpp"
#include "./lapack-like/HPDInverse.hpp"
#include "./lapack-like/HPSDCholesky.hpp"
#include "./lapack-like/HPSDSquareRoot.hpp"
#include "./lapack-like/Inverse.hpp"
#include "./lapack-like/LDL.hpp"
#include "./lapack-like/LQ.hpp"
#include "./lapack-like/LU.hpp"
#include "./lapack-like/Norm.hpp"
#include "./lapack-like/PivotParity.hpp"
#include "./lapack-like/Polar.hpp"
#include "./lapack-like/Pseudoinverse.hpp"
#include "./lapack-like/QDWH.hpp"
#include "./lapack-like/QR.hpp"
#include "./lapack-like/Reflector.hpp"
#include "./lapack-like/SkewHermitianEig.hpp"
#include "./lapack-like/SortEig.hpp"
#include "./lapack-like/SVD.hpp"
#include "./lapack-like/Trace.hpp"
#include "./lapack-like/TriangularInverse.hpp"
#include "./lapack-like/TwoNormLowerBound.hpp"
#include "./lapack-like/TwoNormUpperBound.hpp"

namespace elem {

template<typename F>
inline void
Cholesky( UpperOrLower uplo, Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("Cholesky");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
#endif
    const char uploChar = UpperOrLowerToChar( uplo );
    lapack::Cholesky( uploChar, A.Height(), A.Buffer(), A.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
Hegst( LeftOrRight side, UpperOrLower uplo, Matrix<F>& A, const Matrix<F>& B )
{
#ifndef RELEASE
    PushCallStack("Hegst");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    if( B.Height() != B.Width() )
        throw std::logic_error("B must be square");
    if( A.Height() != B.Height() )
        throw std::logic_error("A and B must be the same size");
#endif
    const int itype = ( side==LEFT ? 2 : 1 );
    const char uploChar = UpperOrLowerToChar( uplo );
    lapack::Hegst
    ( itype, uploChar, A.Height(), 
      A.Buffer(), A.LDim(), B.LockedBuffer(), B.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
LU( Matrix<F>& A, Matrix<int>& p )
{
#ifndef RELEASE
    PushCallStack("LU");
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

template<typename F>
inline typename Base<F>::type
SymmetricNorm( UpperOrLower uplo, const Matrix<F>& A, NormType type )
{ 
#ifndef RELEASE
    PushCallStack("SymmetricNorm");
#endif
    typedef typename Base<F>::type R;
    const R norm = HermitianNorm( uplo, A, type );
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

template<typename F>
inline typename Base<F>::type
SymmetricNorm
( UpperOrLower uplo, const DistMatrix<F>& A, NormType type )
{ 
#ifndef RELEASE
    PushCallStack("SymmetricNorm");
#endif
    typedef typename Base<F>::type R;
    const R norm = HermitianNorm( uplo, A, type );
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

template<typename F>
inline void
TriangularInverse
( UpperOrLower uplo, UnitOrNonUnit diag, Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("TriangularInverse");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
#endif
    const char uploChar = UpperOrLowerToChar( uplo );
    const char diagChar = UnitOrNonUnitToChar( diag );
    lapack::TriangularInverse
    ( uploChar, diagChar, A.Height(), A.Buffer(), A.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem

#endif /* ELEMENTAL_LAPACKLIKE_HPP */

