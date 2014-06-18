/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

// Draw each entry from a normal PDF
template<typename T>
void MakeGaussian( Matrix<T>& A, T mean, Base<T> stddev )
{
    DEBUG_ONLY(CallStackEntry cse("MakeGaussian"))
    EntrywiseFill( A, [=]() { return SampleNormal(mean,stddev); } );
}

template<typename T>
void MakeGaussian( AbstractDistMatrix<T>& A, T mean, Base<T> stddev )
{
    DEBUG_ONLY(CallStackEntry cse("MakeGaussian"))
    if( A.RedundantRank() == 0 )
        MakeGaussian( A.Matrix(), mean, stddev );
    A.BroadcastOver( A.RedundantComm(), 0 );
}

template<typename T>
void MakeGaussian( AbstractBlockDistMatrix<T>& A, T mean, Base<T> stddev )
{
    DEBUG_ONLY(CallStackEntry cse("MakeGaussian"))
    if( A.RedundantRank() == 0 )
        MakeGaussian( A.Matrix(), mean, stddev );
    A.BroadcastOver( A.RedundantComm(), 0 );
}

template<typename T>
void Gaussian( Matrix<T>& A, Int m, Int n, T mean, Base<T> stddev )
{
    DEBUG_ONLY(CallStackEntry cse("Gaussian"))
    A.Resize( m, n );
    MakeGaussian( A, mean, stddev );
}

template<typename T>
void Gaussian
( AbstractDistMatrix<T>& A, Int m, Int n, T mean, Base<T> stddev )
{
    DEBUG_ONLY(CallStackEntry cse("Gaussian"))
    A.Resize( m, n );
    MakeGaussian( A, mean, stddev );
}

template<typename T>
void Gaussian
( AbstractBlockDistMatrix<T>& A, Int m, Int n, T mean, Base<T> stddev )
{
    DEBUG_ONLY(CallStackEntry cse("Gaussian"))
    A.Resize( m, n );
    MakeGaussian( A, mean, stddev );
}

#define PROTO(T) \
  template void MakeGaussian \
  ( Matrix<T>& A, T mean, Base<T> stddev ); \
  template void MakeGaussian \
  ( AbstractDistMatrix<T>& A, T mean, Base<T> stddev ); \
  template void MakeGaussian \
  ( AbstractBlockDistMatrix<T>& A, T mean, Base<T> stddev ); \
  template void Gaussian \
  ( Matrix<T>& A, Int m, Int n, T mean, Base<T> stddev ); \
  template void Gaussian \
  ( AbstractDistMatrix<T>& A, Int m, Int n, T mean, Base<T> stddev ); \
  template void Gaussian \
  ( AbstractBlockDistMatrix<T>& A, Int m, Int n, T mean, Base<T> stddev );

PROTO(float)
PROTO(double)
PROTO(Complex<float>)
PROTO(Complex<double>)

} // namespace El
