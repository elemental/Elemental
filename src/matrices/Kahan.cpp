/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

// I haven't decided on the appropriate generalization to complex cosine/sine
// pairs. For now, given phi, we will compute the corresponding partner as the
// real value sqrt(1-|phi|^2)

namespace El {

template<typename F>
void Kahan( Matrix<F>& A, Int n, F phi )
{
    DEBUG_ONLY(CallStackEntry cse("Kahan"))
    A.Resize( n, n );
    const F zeta = Sqrt(F(1)-phi*Conj(phi));
    typedef Base<F> Real;
    auto kahanFill = 
      [=]( Int i, Int j ) -> F
      { if( i == j )      { return      Pow(zeta,Real(i)); }
        else if(  i < j ) { return -phi*Pow(zeta,Real(i)); }
        else              { return F(0);                   } };
    IndexDependentFill( A, function<F(Int,Int)>(kahanFill) );
}

template<typename F>
void Kahan( AbstractDistMatrix<F>& A, Int n, F phi )
{
    DEBUG_ONLY(CallStackEntry cse("Kahan"))
    A.Resize( n, n );
    const F zeta = Sqrt(F(1)-phi*Conj(phi));
    typedef Base<F> Real;
    auto kahanFill = 
      [=]( Int i, Int j ) -> F
      { if( i == j )      { return      Pow(zeta,Real(i)); }
        else if(  i < j ) { return -phi*Pow(zeta,Real(i)); }
        else              { return F(0);                   } };
    IndexDependentFill( A, function<F(Int,Int)>(kahanFill) );
}

template<typename F>
void Kahan( AbstractBlockDistMatrix<F>& A, Int n, F phi )
{
    DEBUG_ONLY(CallStackEntry cse("Kahan"))
    A.Resize( n, n );
    const F zeta = Sqrt(F(1)-phi*Conj(phi));
    typedef Base<F> Real;
    auto kahanFill = 
      [=]( Int i, Int j ) -> F
      { if( i == j )      { return      Pow(zeta,Real(i)); }
        else if(  i < j ) { return -phi*Pow(zeta,Real(i)); }
        else              { return F(0);                   } };
    IndexDependentFill( A, function<F(Int,Int)>(kahanFill) );
}

#define PROTO(F) \
  template void Kahan( Matrix<F>& A, Int n, F phi ); \
  template void Kahan( AbstractDistMatrix<F>& A, Int n, F phi ); \
  template void Kahan( AbstractBlockDistMatrix<F>& A, Int n, F phi ); 

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
