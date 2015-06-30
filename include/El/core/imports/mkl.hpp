/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_IMPORTS_MKL_HPP
#define EL_IMPORTS_MKL_HPP

#ifdef EL_HAVE_MKL
namespace El {
namespace mkl {

void csrmv
( Orientation orientation, BlasInt m, BlasInt k,
  float alpha, const char* matDescrA,
  const float* val, const BlasInt* indx, 
  const BlasInt* pntrb, const BlasInt* pntre,
  const float* x, float beta, float* y );
void csrmv
( Orientation orientation, BlasInt m, BlasInt k,
  double alpha, const char* matDescrA,
  const double* val, const BlasInt* indx, 
  const BlasInt* pntrb, const BlasInt* pntre,
  const double* x, double beta, double* y );
void csrmv
( Orientation orientation, BlasInt m, BlasInt k,
  Complex<float> alpha, const char* matDescrA,
  const Complex<float>* val, const BlasInt* indx, 
  const BlasInt* pntrb, const BlasInt* pntre,
  const Complex<float>* x, Complex<float> beta, Complex<float>* y );
void csrmv
( Orientation orientation, BlasInt m, BlasInt k,
  Complex<double> alpha, const char* matDescrA,
  const Complex<double>* val, const BlasInt* indx, 
  const BlasInt* pntrb, const BlasInt* pntre,
  const Complex<double>* x, Complex<double> beta, Complex<double>* y );

} // namespace mkl
} // namespace El
#endif // ifdef EL_HAVE_MKL

#endif // ifndef EL_IMPORTS_MKL_HPP
