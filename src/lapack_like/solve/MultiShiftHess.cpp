/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

// NOTE: These algorithms are adaptations and/or extensions of Alg. 2 from
//       Greg Henry's "The shifted Hessenberg system solve computation".
//       It is important to note that the Givens rotation definition in
//       said paper is the adjoint of the LAPACK definition (as well as 
//       leaving out a conjugation necessary for the complex case).

// TODO: Expose interface which does not solve against Q, since that step
//       is not needed for two-norm pseudospectra

namespace El {
namespace mshs {

template<typename F>
inline void
LN
( F alpha,
  const Matrix<F>& H,
  const Matrix<F>& shifts,
        Matrix<F>& X ) 
{
    DEBUG_ONLY(CSE cse("mshs::LN"))
    X *= alpha;

    const Int m = X.Height();
    const Int n = X.Width();
    if( m == 0 )
        return;

    // Initialize storage for Givens rotations
    typedef Base<F> Real;
    Matrix<Real> C(m,n);
    Matrix<F> S(m,n);

    // Initialize the workspace for shifted columns of H
    Matrix<F> W(m,n);
    for( Int j=0; j<n; ++j )
    {
        MemCopy( W.Buffer(0,j), H.LockedBuffer(), m );
        W.Update( 0, j, -shifts.Get(j,0) );
    }
     
    // Simultaneously find the LQ factorization and solve against L
    for( Int k=0; k<m-1; ++k )
    {
        auto hB = H( IR(k+2,m), IR(k+1) );
        const F etakkp1 = H.Get(k,k+1);
        const F etakp1kp1 = H.Get(k+1,k+1);
        for( Int j=0; j<n; ++j )
        {
            // Find the Givens rotation needed to zero H(k,k+1),
            //   | c        s | | H(k,k)   | = | gamma |
            //   | -conj(s) c | | H(k,k+1) |   | 0     |
            Real c; F s;
            lapack::Givens( W.Get(k,j), etakkp1, &c, &s );
            C.Set( k, j, c );
            S.Set( k, j, s );

            // The new diagonal value of L
            const F lambdakk = c*W.Get(k,j) + s*etakkp1;

            // Divide our current entry of x by the diagonal value of L
            X.Set( k, j, X.Get(k,j)/lambdakk );

            // x(k+1:end) -= x(k) * L(k+1:end,k), where
            // L(k+1:end,k) = c H(k+1:end,k) + s H(k+1:end,k+1). We express this
            // more concisely as xB -= x(k) * ( c wB + s hB ).
            // Note that we carefully handle updating the k+1'th entry since
            // it is shift-dependent.
            const F mu = shifts.Get( j, 0 );
            const F xc = X.Get(k,j)*c;
            const F xs  = X.Get(k,j)*s;
            X.Update( k+1, j, -xc*W.Get(k+1,j)-xs*(etakp1kp1-mu) );
            blas::Axpy
            ( m-(k+2), -xc, W.LockedBuffer(k+2,j), 1, X.Buffer(k+2,j), 1 );
            blas::Axpy
            ( m-(k+2), -xs, hB.LockedBuffer(),     1, X.Buffer(k+2,j), 1 );

            // Change the working vector, wB, from representing a fully-updated
            // portion of the k'th column of H from the end of the last 
            // to a fully-updated portion of the k+1'th column of this iteration
            //
            // w(k+1:end) := -conj(s) H(k+1:end,k) + c H(k+1:end,k+1)
            W.Set( k+1, j, -Conj(s)*W.Get(k+1,j)+c*(etakp1kp1-mu) );
            blas::Scal( m-(k+2), -Conj(s), W.Buffer(k+2,j), 1 );
            blas::Axpy
            ( m-(k+2), F(c), hB.LockedBuffer(), 1, W.Buffer(k+2,j), 1 );
        }
    }
    // Divide x(end) by L(end,end)
    for( Int j=0; j<n; ++j )
        X.Set( m-1, j, X.Get(m-1,j)/W.Get(m-1,j) );

    // Solve against Q
    for( Int j=0; j<n; ++j )        
    {
        F* x = X.Buffer(0,j);
        const Real* c = C.LockedBuffer(0,j);
        const F*    s = S.LockedBuffer(0,j);
        F tau0 = x[m-1];
        for( Int k=m-2; k>=0; --k )
        {
            F tau1 = x[k];
            x[k+1] =       c[k] *tau0 + s[k]*tau1;
            tau0   = -Conj(s[k])*tau0 + c[k]*tau1;
        }
        x[0] = tau0;
    }
}

template<typename F>
inline void
UN
( F alpha,
  const Matrix<F>& H,
  const Matrix<F>& shifts,
        Matrix<F>& X ) 
{
    DEBUG_ONLY(CSE cse("mshs::UN"))
    X *= alpha;

    const Int m = X.Height();
    const Int n = X.Width();
    if( m == 0 )
        return;

    // Initialize storage for Givens rotations
    typedef Base<F> Real;
    Matrix<Real> C(m,n);
    Matrix<F> S(m,n);

    // Initialize the workspace for shifted columns of H
    Matrix<F> W(m,n);
    for( Int j=0; j<n; ++j )
    {
        MemCopy( W.Buffer(0,j), H.LockedBuffer(0,m-1), m );
        W.Update( m-1, j, -shifts.Get(j,0) );
    }
     
    // Simultaneously form the RQ factorization and solve against R
    for( Int k=m-1; k>0; --k )
    {
        auto hT = H( IR(0,k-1), IR(k-1) );
        const F etakkm1 = H.Get(k,k-1);
        const F etakm1km1 = H.Get(k-1,k-1);
        for( Int j=0; j<n; ++j )
        {
            // Find the Givens rotation needed to zero H(k,k-1),
            //   | c        s | | H(k,k)   | = | gamma |
            //   | -conj(s) c | | H(k,k-1) |   | 0     |
            Real c; F s;
            lapack::Givens( W.Get(k,j), etakkm1, &c, &s );
            C.Set( k, j, c );
            S.Set( k, j, s );

            // The new diagonal value of R
            const F rhokk = c*W.Get(k,j) + s*etakkm1;

            // Divide our current entry of x by the diagonal value of R
            X.Set( k, j, X.Get(k,j)/rhokk );

            // x(0:k-1) -= x(k) * R(0:k-1,k), where
            // R(0:k-1,k) = c H(0:k-1,k) + s H(0:k-1,k-1). We express this
            // more concisely as xT -= x(k) * ( c wT + s hT ).
            // Note that we carefully handle updating the k-1'th entry since
            // it is shift-dependent.
            const F mu = shifts.Get( j, 0 );
            const F xc = X.Get(k,j)*c;
            const F xs  = X.Get(k,j)*s;
            blas::Axpy( k-1, -xc, W.LockedBuffer(0,j), 1, X.Buffer(0,j), 1 );
            blas::Axpy( k-1, -xs, hT.LockedBuffer(),   1, X.Buffer(0,j), 1 );
            X.Update( k-1, j, -xc*W.Get(k-1,j)-xs*(etakm1km1-mu) );

            // Change the working vector, wT, from representing a fully-updated
            // portion of the k'th column of H from the end of the last 
            // to a fully-updated portion of the k-1'th column of this iteration
            //
            // w(0:k-1) := -conj(s) H(0:k-1,k) + c H(0:k-1,k-1)
            blas::Scal( k-1, -Conj(s), W.Buffer(0,j), 1 );
            blas::Axpy( k-1, F(c), hT.LockedBuffer(), 1, W.Buffer(0,j), 1 );
            W.Set( k-1, j, -Conj(s)*W.Get(k-1,j)+c*(etakm1km1-mu) );
        }
    }
    // Divide x(0) by R(0,0)
    for( Int j=0; j<n; ++j )
        X.Set( 0, j, X.Get(0,j)/W.Get(0,j) );

    // Solve against Q
    for( Int j=0; j<n; ++j )        
    {
        F* x = X.Buffer(0,j);
        const Real* c = C.LockedBuffer(0,j);
        const F*    s = S.LockedBuffer(0,j);
        F tau0 = x[0];
        for( Int k=1; k<m; ++k )
        {
            F tau1 = x[k];
            x[k-1] =       c[k] *tau0 + s[k]*tau1;
            tau0   = -Conj(s[k])*tau0 + c[k]*tau1;
        }
        x[m-1] = tau0;
    }
}

// NOTE: A [VC,* ] distribution might be most appropriate for the 
//       Hessenberg matrices since whole columns will need to be formed 
//       on every process and this distribution will keep the communication 
//       balanced.

template<typename F>
inline void
LN
( F alpha,
  const ElementalMatrix<F>& H, 
  const ElementalMatrix<F>& shiftsPre,
        ElementalMatrix<F>& XPre ) 
{
    DEBUG_ONLY(CSE cse("mshs::LN"))

    DistMatrixReadWriteProxy<F,F,STAR,VR> XProx( XPre );
    auto& X = XProx.Get();

    ElementalProxyCtrl ctrl; 
    ctrl.colConstrain = true;
    ctrl.colAlign = X.RowAlign();

    DistMatrixReadProxy<F,F,VR,STAR> shiftsProx( shiftsPre, ctrl );
    auto& shifts = shiftsProx.GetLocked();

    X *= alpha;

    const Int m = X.Height();
    const Int nLoc = X.LocalWidth();
    if( m == 0 )
        return;

    // Initialize storage for Givens rotations
    typedef Base<F> Real;
    Matrix<Real> C(m,nLoc);
    Matrix<F> S(m,nLoc);

    // Initialize the workspace for shifted columns of H
    Matrix<F> W(m,nLoc);
    {
        unique_ptr<ElementalMatrix<F>> 
          h0( H.Construct(H.Grid(),H.Root()) );
        LockedView( *h0, H, ALL, IR(0) );
        DistMatrix<F,STAR,STAR> h0_STAR_STAR( *h0 );
        for( Int jLoc=0; jLoc<nLoc; ++jLoc )
        {
            MemCopy( W.Buffer(0,jLoc), h0_STAR_STAR.LockedBuffer(), m );
            W.Update( 0, jLoc, -shifts.GetLocal(jLoc,0) );
        }
    }
     
    // Simultaneously find the LQ factorization and solve against L
    unique_ptr<ElementalMatrix<F>> hB( H.Construct(H.Grid(),H.Root()) );
    DistMatrix<F,STAR,STAR> hB_STAR_STAR( H.Grid() );
    for( Int k=0; k<m-1; ++k )
    {
        LockedView( *hB, H, IR(k+2,m), IR(k+1) );
        hB_STAR_STAR = *hB;
        const F etakkp1 = H.Get(k,k+1);
        const F etakp1kp1 = H.Get(k+1,k+1);
        for( Int jLoc=0; jLoc<nLoc; ++jLoc )
        {
            // Find the Givens rotation needed to zero H(k,k+1),
            //   | c        s | | H(k,k)   | = | gamma |
            //   | -conj(s) c | | H(k,k+1) |   | 0     |
            Real c; F s;
            lapack::Givens( W.Get(k,jLoc), etakkp1, &c, &s );
            C.Set( k, jLoc, c );
            S.Set( k, jLoc, s );

            // The new diagonal value of L
            const F lambdakk = c*W.Get(k,jLoc) + s*etakkp1;

            // Divide our current entry of x by the diagonal value of L
            X.SetLocal( k, jLoc, X.GetLocal(k,jLoc)/lambdakk );

            // x(k+1:end) -= x(k) * L(k+1:end,k), where
            // L(k+1:end,k) = c H(k+1:end,k) + s H(k+1:end,k+1). We express this
            // more concisely as xB -= x(k) * ( c wB + s hB ).
            // Note that we carefully handle updating the k+1'th entry since
            // it is shift-dependent.
            const F mu = shifts.GetLocal( jLoc, 0 );
            const F xc = X.GetLocal(k,jLoc)*c;
            const F xs  = X.GetLocal(k,jLoc)*s;
            X.UpdateLocal( k+1, jLoc, -xc*W.Get(k+1,jLoc)-xs*(etakp1kp1-mu) );
            blas::Axpy
            ( m-(k+2), -xc, W.LockedBuffer(k+2,jLoc), 1, 
                            X.Buffer(k+2,jLoc),       1 );
            blas::Axpy
            ( m-(k+2), -xs, hB_STAR_STAR.LockedBuffer(), 1, 
                            X.Buffer(k+2,jLoc),          1 );

            // Change the working vector, wB, from representing a fully-updated
            // portion of the k'th column of H from the end of the last 
            // to a fully-updated portion of the k+1'th column of this iteration
            //
            // w(k+1:end) := -conj(s) H(k+1:end,k) + c H(k+1:end,k+1)
            W.Set( k+1, jLoc, -Conj(s)*W.Get(k+1,jLoc)+c*(etakp1kp1-mu) );
            blas::Scal( m-(k+2), -Conj(s), W.Buffer(k+2,jLoc), 1 );
            blas::Axpy
            ( m-(k+2), F(c), hB_STAR_STAR.LockedBuffer(), 1, 
                             W.Buffer(k+2,jLoc),          1 );
        }
    }
    // Divide x(end) by L(end,end)
    for( Int jLoc=0; jLoc<nLoc; ++jLoc )
        X.SetLocal( m-1, jLoc, X.GetLocal(m-1,jLoc)/W.Get(m-1,jLoc) );

    // Solve against Q
    for( Int jLoc=0; jLoc<nLoc; ++jLoc ) 
    {
        F* x = X.Buffer(0,jLoc);
        const Real* c = C.LockedBuffer(0,jLoc);
        const F*    s = S.LockedBuffer(0,jLoc);
        F tau0 = x[m-1];
        for( Int k=m-2; k>=0; --k )
        {
            F tau1 = x[k];
            x[k+1] =       c[k] *tau0 + s[k]*tau1;
            tau0   = -Conj(s[k])*tau0 + c[k]*tau1;
        }
        x[0] = tau0;
    }
}

template<typename F>
inline void
UN
( F alpha,
  const ElementalMatrix<F>& H, 
  const ElementalMatrix<F>& shiftsPre,
        ElementalMatrix<F>& XPre ) 
{
    DEBUG_ONLY(CSE cse("mshs::UN"))

    DistMatrixReadWriteProxy<F,F,STAR,VR> XProx( XPre );
    auto& X = XProx.Get();

    ElementalProxyCtrl ctrl;
    ctrl.colConstrain = true;
    ctrl.colAlign = X.RowAlign();

    DistMatrixReadProxy<F,F,VR,STAR> shiftsProx( shiftsPre, ctrl );
    auto& shifts = shiftsProx.GetLocked();

    X *= alpha;

    const Int m = X.Height();
    const Int nLoc = X.LocalWidth();
    if( m == 0 )
        return;

    // Initialize storage for Givens rotations
    typedef Base<F> Real;
    Matrix<Real> C(m,nLoc);
    Matrix<F> S(m,nLoc);

    // Initialize the workspace for shifted columns of H
    Matrix<F> W(m,nLoc);
    {
        unique_ptr<ElementalMatrix<F>> 
          hLast( H.Construct(H.Grid(),H.Root()) );
        LockedView( *hLast, H, ALL, IR(m-1) );
        DistMatrix<F,STAR,STAR> hLast_STAR_STAR( *hLast );
        for( Int jLoc=0; jLoc<nLoc; ++jLoc )
        {
            MemCopy( W.Buffer(0,jLoc), hLast_STAR_STAR.LockedBuffer(), m );
            W.Update( m-1, jLoc, -shifts.GetLocal(jLoc,0) );
        }
    }
     
    // Simultaneously form the RQ factorization and solve against R
    unique_ptr<ElementalMatrix<F>> hT( H.Construct(H.Grid(),H.Root()) );
    DistMatrix<F,STAR,STAR> hT_STAR_STAR( H.Grid() );
    for( Int k=m-1; k>0; --k )
    {
        LockedView( *hT, H, IR(0,k-1), IR(k-1) );
        hT_STAR_STAR = *hT;
        const F etakkm1 = H.Get(k,k-1);
        const F etakm1km1 = H.Get(k-1,k-1);
        for( Int jLoc=0; jLoc<nLoc; ++jLoc )
        {
            // Find the Givens rotation needed to zero H(k,k-1),
            //   | c        s | | H(k,k)   | = | gamma |
            //   | -conj(s) c | | H(k,k-1) |   | 0     |
            Real c; F s;
            lapack::Givens( W.Get(k,jLoc), etakkm1, &c, &s );
            C.Set( k, jLoc, c );
            S.Set( k, jLoc, s );

            // The new diagonal value of R
            const F rhokk = c*W.Get(k,jLoc) + s*etakkm1;

            // Divide our current entry of x by the diagonal value of R
            X.SetLocal( k, jLoc, X.GetLocal(k,jLoc)/rhokk );

            // x(0:k-1) -= x(k) * R(0:k-1,k), where
            // R(0:k-1,k) = c H(0:k-1,k) + s H(0:k-1,k-1). We express this
            // more concisely as xT -= x(k) * ( c wT + s hT ).
            // Note that we carefully handle updating the k-1'th entry since
            // it is shift-dependent.
            const F mu = shifts.GetLocal( jLoc, 0 );
            const F xc = X.GetLocal(k,jLoc)*c;
            const F xs  = X.GetLocal(k,jLoc)*s;
            blas::Axpy
            ( k-1, -xc, W.LockedBuffer(0,jLoc),      1, X.Buffer(0,jLoc), 1 );
            blas::Axpy
            ( k-1, -xs, hT_STAR_STAR.LockedBuffer(), 1, X.Buffer(0,jLoc), 1 );
            X.UpdateLocal( k-1, jLoc, -xc*W.Get(k-1,jLoc)-xs*(etakm1km1-mu) );

            // Change the working vector, wT, from representing a fully-updated
            // portion of the k'th column of H from the end of the last 
            // to a fully-updated portion of the k-1'th column of this iteration
            //
            // w(0:k-1) := -conj(s) H(0:k-1,k) + c H(0:k-1,k-1)
            blas::Scal( k-1, -Conj(s), W.Buffer(0,jLoc), 1 );
            blas::Axpy( k-1, F(c), hT_STAR_STAR.LockedBuffer(), 1, 
                                   W.Buffer(0,jLoc),            1 );
            W.Set( k-1, jLoc, -Conj(s)*W.Get(k-1,jLoc)+c*(etakm1km1-mu) );
        }
    }
    for( Int jLoc=0; jLoc<nLoc; ++jLoc )
        X.SetLocal( 0, jLoc, X.GetLocal(0,jLoc)/W.Get(0,jLoc) );

    // Solve against Q
    for( Int jLoc=0; jLoc<nLoc; ++jLoc )
    {
        F* x = X.Buffer(0,jLoc);
        const Real* c = C.LockedBuffer(0,jLoc);
        const F*    s = S.LockedBuffer(0,jLoc);
        F tau0 = x[0];
        for( Int k=1; k<m; ++k )
        {
            F tau1 = x[k];
            x[k-1] =       c[k] *tau0 + s[k]*tau1;
            tau0   = -Conj(s[k])*tau0 + c[k]*tau1;
        }
        x[m-1] = tau0;
    }
}

// TODO: UT and LT

} // namespace mshs

template<typename F>
void MultiShiftHessSolve
( UpperOrLower uplo,
  Orientation orientation,
  F alpha,
  const Matrix<F>& H,
  const Matrix<F>& shifts,
        Matrix<F>& X )
{
    DEBUG_ONLY(CSE cse("MultiShiftHessSolve"))
    if( uplo == LOWER )
    {
        if( orientation == NORMAL )
            mshs::LN( alpha, H, shifts, X );
        else
            LogicError("This option is not yet supported");
    }
    else
    {
        if( orientation == NORMAL )
            mshs::UN( alpha, H, shifts, X );
        else
            LogicError("This option is not yet supported");
    }
}

template<typename F>
void MultiShiftHessSolve
( UpperOrLower uplo,
  Orientation orientation,
  F alpha,
  const ElementalMatrix<F>& H,
  const ElementalMatrix<F>& shifts, 
        ElementalMatrix<F>& X )
{
    DEBUG_ONLY(CSE cse("MultiShiftHessSolve"))
    if( uplo == LOWER )
    {
        if( orientation == NORMAL )
            mshs::LN( alpha, H, shifts, X );
        else
            LogicError("This option is not yet supported");
    }
    else
    {
        if( orientation == NORMAL )
            mshs::UN( alpha, H, shifts, X );
        else
            LogicError("This option is not yet supported");
    }
}

#define PROTO(F) \
  template void MultiShiftHessSolve \
  ( UpperOrLower uplo, \
    Orientation orientation, \
    F alpha, \
    const Matrix<F>& H, \
    const Matrix<F>& shifts, \
          Matrix<F>& X ); \
  template void MultiShiftHessSolve \
  ( UpperOrLower uplo, \
    Orientation orientation, \
    F alpha, \
    const ElementalMatrix<F>& H, \
    const ElementalMatrix<F>& shifts, \
          ElementalMatrix<F>& X );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

} // namespace El
