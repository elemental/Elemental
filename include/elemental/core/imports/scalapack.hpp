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

#ifdef ELEM_HAVE_SCALAPACK

namespace elem {

namespace blacs {
// Basic Linear Algebra Communication Subprograms
// ==============================================
int Handle( MPI_Comm comm );
int GridInit( int bhandle, bool colMajor, int gridHeight, int gridWidth );
int GridHeight( int context );
int GridWidth( int context );
int GridRow( int context );
int GridCol( int context );
void FreeHandle( int bhandle );
void FreeGrid( int context );
void Exit( bool finished=false );

typedef std::array<int,9> Desc;

} // namespace blacs

namespace scalapack {

// ScaLAPACK proper
// ================

// Hessenberg Schur decomposition via the QR algorithm
// ---------------------------------------------------
// NOTE: In all of these routines, the matrix needs to be explicitly 
//       upper-Hessenberg before the call, otherwise behavior is unpredictable

void HessenbergSchur
( int n, float* H, const int* desch, scomplex* w, bool fullTriangle=false );
void HessenbergSchur
( int n, double* H, const int* desch, dcomplex* w, bool fullTriangle=false );
void HessenbergSchur
( int n, scomplex* H, const int* desch, scomplex* w, bool fullTriangle=false );
void HessenbergSchur
( int n, dcomplex* H, const int* desch, dcomplex* w, bool fullTriangle=false );

void HessenbergSchur
( int n, float* H, const int* desch, scomplex* w, 
  float* Q, const int* descq, bool fullTriangle=true, bool multiplyQ=false );
void HessenbergSchur
( int n, double* H, const int* desch, dcomplex* w, 
  double* Q, const int* descq, bool fullTriangle=true, bool multiplyQ=false );
void HessenbergSchur
( int n, scomplex* H, const int* desch, scomplex* w, 
  scomplex* Q, const int* descq, bool fullTriangle=true, bool multiplyQ=false );
void HessenbergSchur
( int n, dcomplex* H, const int* desch, dcomplex* w, 
  dcomplex* Q, const int* descq, bool fullTriangle=true, bool multiplyQ=false );

// Hessenberg eigenvalues/pairs
// ----------------------------

void HessenbergEig( int n, float* H, const int* desch, scomplex* w );
void HessenbergEig( int n, double* H, const int* desch, dcomplex* w );
void HessenbergEig( int n, scomplex* H, const int* desch, scomplex* w );
void HessenbergEig( int n, dcomplex* H, const int* desch, dcomplex* w );

// TODO: Compute the eigenvectors

} // namespace scalapack
} // namespace elem

#endif // ifdef ELEM_HAVE_SCALAPACK

#endif // ifndef ELEM_IMPORTS_SCALAPACK_HPP
