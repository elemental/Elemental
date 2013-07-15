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

// Throws an error if PMRRR was not built along with Elemental
inline void EnsurePMRRR()
{
#ifndef HAVE_PMRRR
    throw std::logic_error("PMRRR is required for this routine");
#endif
}

// Compute the eigenvalues of a Hermitian matrix
template<typename F>
void HermitianEig
( UpperOrLower uplo, Matrix<F>& A, Matrix<BASE(F)>& w );
template<typename F>
void HermitianEig
( UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<BASE(F),VR,STAR>& w );

// Compute the full eigenvalue decomposition of a Hermitian matrix
template<typename F>
void HermitianEig
( UpperOrLower uplo, Matrix<F>& A, Matrix<BASE(F)>& w, Matrix<F>& Z );
template<typename F>
void HermitianEig
( UpperOrLower uplo, 
  DistMatrix<F>& A, DistMatrix<BASE(F),VR,STAR>& w, DistMatrix<F>& paddedZ );

// Compute the eigenvalues of a Hermitian matrix within a selected range
template<typename F>
void HermitianEig
( UpperOrLower uplo, Matrix<F>& A, Matrix<BASE(F)>& w,
  int lowerBound, int upperBound );
template<typename F>
void HermitianEig
( UpperOrLower uplo, Matrix<F>& A, Matrix<BASE(F)>& w,
  BASE(F) lowerBound, BASE(F) upperBound );
template<typename F>
void HermitianEig
( UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<BASE(F),VR,STAR>& w, 
  int lowerBound, int upperBound );
template<typename F>
void HermitianEig
( UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<BASE(F),VR,STAR>& w,
  BASE(F) lowerBound, BASE(F) upperBound );

// Compute a selected set of eigenpairs of a Hermitian matrix
template<typename F>
void HermitianEig
( UpperOrLower uplo, 
  Matrix<F>& A, Matrix<BASE(F)>& w, Matrix<F>& Z,
  int lowerBound, int upperBound );
template<typename F>
void HermitianEig
( UpperOrLower uplo, 
  Matrix<F>& A, Matrix<BASE(F)>& w, Matrix<F>& Z,
  BASE(F) lowerBound, BASE(F) upperBound );
template<typename F>
void HermitianEig
( UpperOrLower uplo, 
  DistMatrix<F>& A, DistMatrix<BASE(F),VR,STAR>& w, DistMatrix<F>& paddedZ,
  int lowerBound, int upperBound );
template<typename F>
void HermitianEig
( UpperOrLower uplo, 
  DistMatrix<F>& A, DistMatrix<BASE(F),VR,STAR>& w, DistMatrix<F>& paddedZ,
  BASE(F) lowerBound, BASE(F) upperBound );

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
// Utilities                                                                  //
//----------------------------------------------------------------------------//

template<typename F>
void PivotFunc
( void* inData, void* outData, int* length, mpi::Datatype* datatype );

template<typename F> mpi::Op PivotOp();
template<> mpi::Op PivotOp<float>();
template<> mpi::Op PivotOp<double>();
template<> mpi::Op PivotOp<scomplex>();
template<> mpi::Op PivotOp<dcomplex>();

template<typename F> void CreatePivotOp();
template<> void CreatePivotOp<float>();
template<> void CreatePivotOp<double>();
template<> void CreatePivotOp<scomplex>();
template<> void CreatePivotOp<dcomplex>();

template<typename T> void DestroyPivotOp();
template<> void DestroyPivotOp<float>();
template<> void DestroyPivotOp<double>();
template<> void DestroyPivotOp<scomplex>();
template<> void DestroyPivotOp<dcomplex>();

} // namespace elem

#endif // ifndef LAPACK_DECL_HPP
