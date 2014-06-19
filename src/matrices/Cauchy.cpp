/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename F1,typename F2> 
void Cauchy( Matrix<F1>& A, const std::vector<F2>& x, const std::vector<F2>& y )
{
    DEBUG_ONLY(CallStackEntry cse("Cauchy"))
    const Int m = x.size();
    const Int n = y.size();
    A.Resize( m, n );
    IndexDependentFill
    ( A, [&]( Int i, Int j )
         {
            DEBUG_ONLY(
                // TODO: Use tolerance instead?
                if( x[i] == y[j] )
                    LogicError
                    ( "x[", i, "] = y[", j, "] (", x[i], 
                      ") is not allowed for Cauchy matrices" );
            ) 
            return F1(1)/F1(x[i]-y[j]);
         } );
}

template<typename F1,typename F2>
void Cauchy
( AbstractDistMatrix<F1>& A, 
  const std::vector<F2>& x, const std::vector<F2>& y )
{
    DEBUG_ONLY(CallStackEntry cse("Cauchy"))
    const Int m = x.size();
    const Int n = y.size();
    A.Resize( m, n );
    IndexDependentFill
    ( A, [&]( Int i, Int j )
         {
            DEBUG_ONLY(
                // TODO: Use tolerance instead?
                if( x[i] == y[j] )
                    LogicError
                    ( "x[", i, "] = y[", j, "] (", x[i], 
                      ") is not allowed for Cauchy matrices" );
            ) 
            return F1(1)/F1(x[i]-y[j]);
         } );
}

template<typename F1,typename F2>
void Cauchy
( AbstractBlockDistMatrix<F1>& A, 
  const std::vector<F2>& x, const std::vector<F2>& y )
{
    DEBUG_ONLY(CallStackEntry cse("Cauchy"))
    const Int m = x.size();
    const Int n = y.size();
    A.Resize( m, n );
    IndexDependentFill
    ( A, [&]( Int i, Int j )
         {
            DEBUG_ONLY(
                // TODO: Use tolerance instead?
                if( x[i] == y[j] )
                    LogicError
                    ( "x[", i, "] = y[", j, "] (", x[i], 
                      ") is not allowed for Cauchy matrices" );
            ) 
            return F1(1)/F1(x[i]-y[j]);
         } );
}

#define PROTO_TYPES(F1,F2) \
  template void Cauchy \
  ( Matrix<F1>& A, const std::vector<F2>& x, const std::vector<F2>& y ); \
  template void Cauchy \
  ( AbstractDistMatrix<F1>& A, \
    const std::vector<F2>& x, const std::vector<F2>& y ); \
  template void Cauchy \
  ( AbstractBlockDistMatrix<F1>& A, \
    const std::vector<F2>& x, const std::vector<F2>& y );

#define PROTO_REAL(F) \
  PROTO_TYPES(F,Int) \
  PROTO_TYPES(F,F)

#define PROTO_CPX(F) \
  PROTO_TYPES(F,Int) \
  PROTO_TYPES(F,Base<F>) \
  PROTO_TYPES(F,F)

PROTO_REAL(float)
PROTO_REAL(double)
PROTO_CPX(Complex<float>)
PROTO_CPX(Complex<double>)

} // namespace El
