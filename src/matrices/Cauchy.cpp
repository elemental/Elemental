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
void Cauchy( Matrix<F1>& A, const std::vector<F2>& x, const std::vector<F2>& y )
{
    DEBUG_ONLY(CallStackEntry cse("Cauchy"))
    const Int m = x.size();
    const Int n = y.size();
    A.Resize( m, n );

    const F1 one = F1(1);
    for( Int j=0; j<n; ++j )
    {
        for( Int i=0; i<m; ++i )
        {
            DEBUG_ONLY(
                // TODO: Use tolerance instead?
                if( x[i] == y[j] )
                    LogicError
                    ( "x[", i, "] = y[", j, "] (", x[i], 
                      ") is not allowed for Cauchy matrices" );
            ) 
            A.Set( i, j, one/F1(x[i]-y[j]) );
        }
    }
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

    const F1 one = F1(1);
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
                      ") is not allowed for Cauchy matrices" );
            )
            A.SetLocal( iLoc, jLoc, one/F1(x[i]-y[j]) );
        }
    }
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

    const F1 one = F1(1);
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
                      ") is not allowed for Cauchy matrices" );
            )
            A.SetLocal( iLoc, jLoc, one/F1(x[i]-y[j]) );
        }
    }
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
