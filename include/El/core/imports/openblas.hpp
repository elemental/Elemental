/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_IMPORTS_OPENBLAS_HPP
#define EL_IMPORTS_OPENBLAS_HPP

#ifdef EL_HAVE_OPENBLAS
namespace El {
namespace openblas {

void omatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  float alpha, const float* A, BlasInt lda, float* B, BlasInt ldb );  
void omatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  double alpha, const double* A, BlasInt lda, double* B, BlasInt ldb );  
void omatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  scomplex alpha, const scomplex* A, BlasInt lda, scomplex* B, BlasInt ldb );  
void omatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  dcomplex alpha, const dcomplex* A, BlasInt lda, dcomplex* B, BlasInt ldb );  

// NOTE: These are filler routines since they are not provided by OpenBLAS
template<typename T>
void omatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  T alpha, const T* A, BlasInt lda, T* B, BlasInt ldb );  

void imatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  float alpha, float* A, BlasInt lda, BlasInt ldb );  
void imatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  double alpha, double* A, BlasInt lda, BlasInt ldb );  
void imatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  scomplex alpha, scomplex* A, BlasInt lda, BlasInt ldb );  
void imatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  dcomplex alpha, dcomplex* A, BlasInt lda, BlasInt ldb );  

// NOTE: These are filler routines since they are not provided by OpenBLAS
template<typename T>
void imatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  T alpha, T* A, BlasInt lda, BlasInt ldb );  

} // namespace openblas
} // namespace El
#endif // ifdef EL_HAVE_OPENBLAS

#endif // ifndef EL_IMPORTS_OPENBLAS_HPP
