/*
   Copyright (C) 1992-2008 The University of Tennessee
   All rights reserved.

   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is loosely based upon the LAPACK routines dlarfg.f and zlarfg.f.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

#include "./Householder/Col.hpp"
#include "./Householder/Row.hpp"

namespace El {

//
// The LAPACK convention (see {s,d,c,z}larfg [CITATION]) defines tau such that
//
//   H = I - tau [1; v] [1; v]',
//
// but adjoint(H) [chi; x] = [beta; 0]. 
//
// Elemental simply uses H [chi; x] = [beta; 0].
//
// On exit, chi is overwritten with beta, and x is overwritten with v.
//
// Another major difference from LAPACK is in the treatment of the special case 
// of x=0, where LAPACK would put H := I, which is not a valid Householder 
// reflector. We instead use the valid Householder reflector:
//    H [chi; 0] = [-chi; 0],
// which is accomplished by setting tau=2, and v=0.
//

template<typename F>
F LeftReflector( F& chi, Matrix<F>& x )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( x.Height() != 1 && x.Width() != 1 )
          LogicError("x must be a vector");
    )
    if( x.Width() == 1 )
        return lapack::Reflector( x.Height()+1, chi, x.Buffer(), 1 );
    else
        return lapack::Reflector( x.Width()+1, chi, x.Buffer(), x.LDim() );
}

template<typename F>
F LeftReflector( Matrix<F>& chi, Matrix<F>& x )
{
    DEBUG_CSE

    F alpha = chi(0);
    const F tau = LeftReflector( alpha, x );
    chi(0) = alpha;

    return tau;
}

template<typename F>
F LeftReflector( ElementalMatrix<F>& chi, ElementalMatrix<F>& x )
{
    DEBUG_CSE
    DEBUG_ONLY(
      AssertSameGrids( chi, x );
      if( chi.Height() != 1 || chi.Width() != 1 )
          LogicError("chi must be a scalar");
      if( x.Width() != 1 )
          LogicError("x must be a column vector");
      if( chi.Root() != x.Root() )
          LogicError("Roots must be the same");
    )
    F tau;
    if( x.CrossRank() == x.Root() )
    {
        if( x.RowRank() == x.RowAlign() )
            tau = reflector::Col( chi, x );
        mpi::Broadcast( tau, x.RowAlign(), x.RowComm() );
    }
    mpi::Broadcast( tau, x.Root(), x.CrossComm() );
    return tau;
}

template<typename F>
F LeftReflector( F& chi, ElementalMatrix<F>& x )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( x.Width() != 1 )
          LogicError("x must be a column vector");
    )
    F tau;
    if( x.CrossRank() == x.Root() )
    {
        if( x.RowRank() == x.RowAlign() )
            tau = reflector::Col( chi, x );
        mpi::Broadcast( tau, x.RowAlign(), x.RowComm() );
    }
    mpi::Broadcast( tau, x.Root(), x.CrossComm() );
    return tau;
}

//
// Defines tau and v such that
//
//   H = I - tau [1; v] [1, v'],
//
// and [chi x] H = [beta 0]
//

template<typename F>
F RightReflector( Matrix<F>& chi, Matrix<F>& x )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( chi.Height() != 1 || chi.Width() != 1 )
          LogicError("chi must be a scalar");
      if( x.Height() != 1 && x.Width() != 1 )
          LogicError("x must be a vector");
    )
    const F tau = LeftReflector( chi, x );
    // There is no need to conjugate chi, it should be real now
    Conjugate( x );
    return tau;
}

template<typename F>
F RightReflector( F& chi, Matrix<F>& x )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( x.Height() != 1 && x.Width() != 1 )
          LogicError("x must be a vector");
    )
    const F tau = LeftReflector( chi, x );
    // There is no need to conjugate chi, it should be real now
    Conjugate( x );
    return tau;
}

template<typename F>
F RightReflector( ElementalMatrix<F>& chi, ElementalMatrix<F>& x )
{
    DEBUG_CSE
    DEBUG_ONLY(
      AssertSameGrids( chi, x );
      if( chi.Height() != 1 || chi.Width() != 1 )
          LogicError("chi must be a scalar");
      if( x.Height() != 1 )
          LogicError("x must be a row vector");
      if( chi.Root() != x.Root() )
          LogicError("Roots must be the same");
    )
    F tau;
    if( x.CrossRank() == x.Root() )
    {
        if( x.ColRank() == x.ColAlign() )
            tau = reflector::Row( chi, x );
        mpi::Broadcast( tau, x.ColAlign(), x.ColComm() );
    }
    mpi::Broadcast( tau, x.Root(), x.CrossComm() );
    return tau;
}

template<typename F>
F RightReflector( F& chi, ElementalMatrix<F>& x )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( x.Height() != 1 )
          LogicError("x must be a row vector");
    )
    F tau;
    if( x.CrossRank() == x.Root() )
    {
        if( x.ColRank() == x.ColAlign() )
            tau = reflector::Row( chi, x );
        mpi::Broadcast( tau, x.ColAlign(), x.ColComm() );
    }
    mpi::Broadcast( tau, x.Root(), x.CrossComm() );
    return tau;
}

#define PROTO(F) \
  template F LeftReflector( F& chi, Matrix<F>& x ); \
  template F LeftReflector( F& chi, ElementalMatrix<F>& x ); \
  template F LeftReflector( Matrix<F>& chi, Matrix<F>& x ); \
  template F LeftReflector \
  ( ElementalMatrix<F>& chi, ElementalMatrix<F>& x ); \
  template F RightReflector( F& chi, Matrix<F>& x ); \
  template F RightReflector( F& chi, ElementalMatrix<F>& x ); \
  template F RightReflector( Matrix<F>& chi, Matrix<F>& x ); \
  template F RightReflector \
  ( ElementalMatrix<F>& chi, ElementalMatrix<F>& x ); \
  template F reflector::Col( F& chi, ElementalMatrix<F>& x ); \
  template F reflector::Col \
  ( ElementalMatrix<F>& chi, ElementalMatrix<F>& x ); \
  template F reflector::Row( F& chi, ElementalMatrix<F>& x ); \
  template F reflector::Row \
  ( ElementalMatrix<F>& chi, ElementalMatrix<F>& x );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
