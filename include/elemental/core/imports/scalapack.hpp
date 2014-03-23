/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_IMPORTS_SCALAPACK_HPP
#define ELEM_IMPORTS_SCALAPACK_HPP

namespace elem {
namespace scalapack {

void HessenbergSchur
( int n, float* H, const int* desch, scomplex* w );
void HessenbergSchur
( int n, double* H, const int* desch, dcomplex* w );
void HessenbergSchur
( int n, scomplex* H, const int* desch, scomplex* w );
void HessenbergSchur
( int n, dcomplex* H, const int* desch, dcomplex* w );

void HessenbergSchur
( int n, float* H, const int* desch, scomplex* w, 
  float* U, const int* descu );
void HessenbergSchur
( int n, double* H, const int* desch, dcomplex* w, 
  double* U, const int* descu );
void HessenbergSchur
( int n, scomplex* H, const int* desch, scomplex* w, 
  scomplex* U, const int* descu );
void HessenbergSchur
( int n, dcomplex* H, const int* desch, dcomplex* w, 
  dcomplex* U, const int* descu );

} // namespace scalapack
} // namespace elem

#endif // ifndef ELEM_IMPORTS_SCALAPACK_HPP
