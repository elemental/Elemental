/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_RQ_EXPLICIT_HPP
#define EL_RQ_EXPLICIT_HPP

namespace El {
namespace rq {

template<typename F>
void ExplicitTriang( Matrix<F>& A )
{
    DEBUG_CSE
    Matrix<F> t;
    Matrix<Base<F>> d;
    Householder( A, t, d );
    MakeTrapezoidal( UPPER, A, A.Width()-A.Height() );
}

template<typename F>
void ExplicitTriang( ElementalMatrix<F>& A )
{
    DEBUG_CSE
    DistMatrix<F,MD,STAR> phase(A.Grid());
    DistMatrix<Base<F>,MD,STAR> signature(A.Grid());
    Householder( A, phase, signature );
    MakeTrapezoidal( UPPER, A, A.Width()-A.Height() );
}

// TODO: ExplicitUnitary

// TODO: Explicit

} // namespace rq
} // namespace El

#endif // ifndef EL_RQ_EXPLICIT_HPP
