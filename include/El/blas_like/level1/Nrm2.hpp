/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_BLAS_NRM2_HPP
#define EL_BLAS_NRM2_HPP

namespace El {

// Forward declarations
template<typename F>
Base<F> FrobeniusNorm( const AbstractDistMatrix<F>& A );
template<typename F>
Base<F> FrobeniusNorm( const DistMultiVec<F>& A );

template<typename F>
Base<F> Nrm2( const Matrix<F>& x )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( x.Height() != 1 && x.Width() != 1 )
          LogicError("Expected vector input");
    )
    Base<F> norm;
    if( x.Width() == 1 )
        norm = blas::Nrm2( x.Height(), x.LockedBuffer(), 1 );
    else
        norm = blas::Nrm2( x.Width(), x.LockedBuffer(), x.LDim() );
    return norm;
}

template<typename F>
Base<F> Nrm2( const AbstractDistMatrix<F>& x )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( x.Height() != 1 && x.Width() != 1 )
          LogicError("x must be a vector");
    )
    return FrobeniusNorm( x );
}

template<typename F>
Base<F> Nrm2( const DistMultiVec<F>& x )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( x.Height() != 1 && x.Width() != 1 )
          LogicError("x must be a vector");
    )
    return FrobeniusNorm( x );
}

#ifdef EL_INSTANTIATE_BLAS_LEVEL1
# define EL_EXTERN
#else
# define EL_EXTERN extern
#endif

#define PROTO(F) \
  EL_EXTERN template Base<F> Nrm2( const Matrix<F>& x ); \
  EL_EXTERN template Base<F> Nrm2( const AbstractDistMatrix<F>& x ); \
  EL_EXTERN template Base<F> Nrm2( const DistMultiVec<F>& x );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

#undef EL_EXTERN

} // namespace El

#endif // ifndef EL_BLAS_NRM2_HPP
