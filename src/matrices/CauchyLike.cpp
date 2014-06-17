/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

namespace El {

template<typename F1,typename F2> 
void CauchyLike
( Matrix<F1>& A,
  const std::vector<F2>& r, const std::vector<F2>& s,
  const std::vector<F2>& x, const std::vector<F2>& y )
{
    DEBUG_ONLY(CallStackEntry cse("CauchyLike"))
    const Int m = r.size();
    const Int n = s.size();
    if( x.size() != (Unsigned)m )
        LogicError("x vector was the wrong length");
    if( y.size() != (Unsigned)n )
        LogicError("y vector was the wrong length");
    A.Resize( m, n );

    for( Int j=0; j<n; ++j )
    {
        for( Int i=0; i<m; ++i )
        {
            DEBUG_ONLY(
                // TODO: Use tolerance instead?
                if( x[i] == y[j] )
                    LogicError
                    ( "x[", i, "] = y[", j, "] (", x[i],
                      ") is not allowed for Cauchy-like matrices" );
            )
            A.Set( i, j, r[i]*s[j]/(x[i]-y[j]) );
        }
    }
}

template<typename F1,typename F2>
void CauchyLike
( AbstractDistMatrix<F1>& A,
  const std::vector<F2>& r, const std::vector<F2>& s, 
  const std::vector<F2>& x, const std::vector<F2>& y )
{
    DEBUG_ONLY(CallStackEntry cse("CauchyLike"))
    const Int m = r.size();
    const Int n = s.size();
    if( x.size() != (Unsigned)m )
        LogicError("x vector was the wrong length");
    if( y.size() != (Unsigned)n )
        LogicError("y vector was the wrong length");
    A.Resize( m, n );

    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = A.GlobalCol(jLoc);
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = A.GlobalRow(iLoc);
            DEBUG_ONLY(
                // TODO: Use tolerance instead?
                if( x[i] == y[j] )
                    LogicError
                    ( "x[", i, "] = y[", j, "] (", x[i],
                      ") is not allowed for Cauchy-like matrices" );
            )
            A.SetLocal( iLoc, jLoc, r[i]*s[j]/(x[i]-y[j]) );
        }
    }
}

template<typename F1,typename F2>
void CauchyLike
( AbstractBlockDistMatrix<F1>& A,
  const std::vector<F2>& r, const std::vector<F2>& s, 
  const std::vector<F2>& x, const std::vector<F2>& y )
{
    DEBUG_ONLY(CallStackEntry cse("CauchyLike"))
    const Int m = r.size();
    const Int n = s.size();
    if( x.size() != (Unsigned)m )
        LogicError("x vector was the wrong length");
    if( y.size() != (Unsigned)n )
        LogicError("y vector was the wrong length");
    A.Resize( m, n );

    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = A.GlobalCol(jLoc);
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = A.GlobalRow(iLoc);
            DEBUG_ONLY(
                // TODO: Use tolerance instead?
                if( x[i] == y[j] )
                    LogicError
                    ( "x[", i, "] = y[", j, "] (", x[i],
                      ") is not allowed for Cauchy-like matrices" );
            )
            A.SetLocal( iLoc, jLoc, r[i]*s[j]/(x[i]-y[j]) );
        }
    }
}

#define PROTO_TYPES(F1,F2) \
  template void CauchyLike \
  ( Matrix<F1>& A, \
    const std::vector<F2>& r, const std::vector<F2>& s, \
    const std::vector<F2>& x, const std::vector<F2>& y ); \
  template void CauchyLike \
  ( AbstractDistMatrix<F1>& A, \
    const std::vector<F2>& r, const std::vector<F2>& s, \
    const std::vector<F2>& x, const std::vector<F2>& y ); \
  template void CauchyLike \
  ( AbstractBlockDistMatrix<F1>& A, \
    const std::vector<F2>& r, const std::vector<F2>& s, \
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
