/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename F1,typename F2> 
void CauchyLike
( Matrix<F1>& A,
  const vector<F2>& r, const vector<F2>& s,
  const vector<F2>& x, const vector<F2>& y )
{
    DEBUG_ONLY(CSE cse("CauchyLike"))
    const Int m = r.size();
    const Int n = s.size();
    if( x.size() != (Unsigned)m )
        LogicError("x vector was the wrong length");
    if( y.size() != (Unsigned)n )
        LogicError("y vector was the wrong length");
    A.Resize( m, n );

    auto cauchyFill = 
      [&]( Int i, Int j ) -> F1
      {
        DEBUG_ONLY(
          // TODO: Use tolerance instead?
          if( x[i] == y[j] )
              LogicError
              ( "x[", i, "] = y[", j, "] (", x[i],
                ") is not allowed for Cauchy matrices" );
        )
        return F1(r[i]*s[j]/x[i]-y[j]);
      };
    IndexDependentFill( A, function<F1(Int,Int)>(cauchyFill) );
}

template<typename F1,typename F2>
void CauchyLike
( AbstractDistMatrix<F1>& A,
  const vector<F2>& r, const vector<F2>& s, 
  const vector<F2>& x, const vector<F2>& y )
{
    DEBUG_ONLY(CSE cse("CauchyLike"))
    const Int m = r.size();
    const Int n = s.size();
    if( x.size() != (Unsigned)m )
        LogicError("x vector was the wrong length");
    if( y.size() != (Unsigned)n )
        LogicError("y vector was the wrong length");
    A.Resize( m, n );

    auto cauchyFill =
      [&]( Int i, Int j ) -> F1
      {
        DEBUG_ONLY(
          // TODO: Use tolerance instead?
          if( x[i] == y[j] )
              LogicError
              ( "x[", i, "] = y[", j, "] (", x[i],
                ") is not allowed for Cauchy matrices" );
        )
        return F1(r[i]*s[j]/x[i]-y[j]);
      };
    IndexDependentFill( A, function<F1(Int,Int)>(cauchyFill) );
}

#define PROTO_TYPES(F1,F2) \
  template void CauchyLike \
  ( Matrix<F1>& A, \
    const vector<F2>& r, const vector<F2>& s, \
    const vector<F2>& x, const vector<F2>& y ); \
  template void CauchyLike \
  ( AbstractDistMatrix<F1>& A, \
    const vector<F2>& r, const vector<F2>& s, \
    const vector<F2>& x, const vector<F2>& y );

#define PROTO_REAL(F) \
  PROTO_TYPES(F,Int) \
  PROTO_TYPES(F,F)

#define PROTO_COMPLEX(F) \
  PROTO_TYPES(F,Int) \
  PROTO_TYPES(F,Base<F>) \
  PROTO_TYPES(F,F)

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

} // namespace El
