/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_LQ_EXPLICIT_HPP
#define EL_LQ_EXPLICIT_HPP

namespace El {
namespace lq {

template<typename F>
void ExplicitTriang( Matrix<F>& A )
{
    DEBUG_CSE
    Matrix<F> phase;
    Matrix<Base<F>> signature;
    LQ( A, phase, signature );

    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);
    A.Resize( m, minDim );
    MakeTrapezoidal( LOWER, A );
}

template<typename F>
void ExplicitTriang( ElementalMatrix<F>& A )
{
    DEBUG_CSE
    const Grid& g = A.Grid();
    DistMatrix<F,MD,STAR> phase(g);
    DistMatrix<Base<F>,MD,STAR> signature(g);
    LQ( A, phase, signature );

    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);
    A.Resize( m, minDim );
    MakeTrapezoidal( LOWER, A );
}

template<typename F>
void ExplicitUnitary( Matrix<F>& A )
{
    DEBUG_CSE
    Matrix<F> phase;
    Matrix<Base<F>> signature;
    LQ( A, phase, signature );

    // TODO: Replace this with an in-place expansion of Q
    Matrix<F> Q;
    Identity( Q, A.Height(), A.Width() );
    lq::ApplyQ( RIGHT, NORMAL, A, phase, signature, Q );
    A = Q;
}

template<typename F>
void ExplicitUnitary( ElementalMatrix<F>& APre )
{
    DEBUG_CSE
    const Grid& g = APre.Grid();

    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre );
    auto& A = AProx.Get();

    DistMatrix<F,MD,STAR> phase(g);
    DistMatrix<Base<F>,MD,STAR> signature(g);
    LQ( A, phase, signature );

    // TODO: Replace this with an in-place expansion of Q
    DistMatrix<F> Q(g);
    Q.AlignWith( A );
    Identity( Q, A.Height(), A.Width() );
    lq::ApplyQ( RIGHT, NORMAL, A, phase, signature, Q );
    Copy( Q, APre );
}

template<typename F>
void Explicit( Matrix<F>& L, Matrix<F>& A )
{
    DEBUG_CSE
    Matrix<F> phase;
    Matrix<Base<F>> signature;
    LQ( A, phase, signature );

    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);
    auto AL = A( IR(0,m), IR(0,minDim) );
    L = AL;
    MakeTrapezoidal( LOWER, L );

    // TODO: Replace this with an in-place expansion of Q
    Matrix<F> Q;
    Identity( Q, A.Height(), A.Width() );
    lq::ApplyQ( RIGHT, NORMAL, A, phase, signature, Q );
    A = Q;
}

template<typename F>
void Explicit( ElementalMatrix<F>& L, ElementalMatrix<F>& APre )
{
    DEBUG_CSE
    const Grid& g = APre.Grid();

    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre );
    auto& A = AProx.Get();

    DistMatrix<F,MD,STAR> phase(g);
    DistMatrix<Base<F>,MD,STAR> signature(g);
    LQ( A, phase, signature );

    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);
    auto AL = A( IR(0,m), IR(0,minDim) );
    Copy( AL, L );
    MakeTrapezoidal( LOWER, L );

    // TODO: Replace this with an in-place expansion of Q
    DistMatrix<F> Q(g);
    Identity( Q, A.Height(), A.Width() );
    lq::ApplyQ( RIGHT, NORMAL, A, phase, signature, Q );
    Copy( Q, APre );
}

} // namespace lq
} // namespace El

#endif // ifndef EL_LQ_EXPLICIT_HPP
