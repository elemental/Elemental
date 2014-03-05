/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_MULTISHIFTHESSSOLVE_HPP
#define ELEM_MULTISHIFTHESSSOLVE_HPP

// NOTE: These algorithms are adaptations and/or extensions of Alg. 2 from
//       Greg Henry's "The shifted Hessenberg system solve computation".
//       It is important to note that the Given's rotation definition in
//       said paper is the adjoint of the LAPACK definition (as well as 
//       leaving out a conjugation necessary for the complex case).

namespace elem {
namespace mshs {

template<typename F>
inline void
UN( F alpha, const Matrix<F>& H, const Matrix<F>& shifts, Matrix<F>& X ) 
{
    DEBUG_ONLY(CallStackEntry cse("mshs::UN"))
    Scale( alpha, X );

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
     
    // Simulataneously factor and solve against R
    for( Int k=m-1; k>0; --k )
    {
        auto hT = LockedView( H, 0, k-1, k-1, 1 );
        const F etakkm1 = H.Get(k,k-1);
        const F etakm1km1 = H.Get(k-1,k-1);
        for( Int j=0; j<n; ++j )
        {
            Real c;
            F s, rho;
            lapack::ComputeGivens( W.Get(k,j), etakkm1, &c, &s, &rho );
            C.Set( k, j, c );
            S.Set( k, j, s );

            X.Set( k, j, X.Get(k,j)/(c*W.Get(k,j)+s*etakkm1) );

            auto xT = View( X, 0, j, k-1, 1 );
            auto wT = View( W, 0, j, k-1, 1 );

            const F xc = X.Get(k,j)*c;
            const F xs  = X.Get(k,j)*s;
            Axpy( -xc, wT, xT );
            Axpy( -xs, hT, xT );

            Scale( -Conj(s), wT );
            Axpy( c, hT, wT );

            const F mu = shifts.Get( j, 0 );
            X.Update( k-1, j, -xc*W.Get(k-1,j)-xs*(etakm1km1-mu) );
            W.Set(    k-1, j, c*(etakm1km1-mu)-Conj(s)*W.Get(k-1,j) );
        }
    }
    for( Int j=0; j<n; ++j )
        X.Set( 0, j, X.Get(0,j)/W.Get(0,j) );

    // Solve against Q
    Matrix<F> t1(n,1), t2(n,1);
    for( Int j=0; j<n; ++j )
        t1.Set( j, 0, X.Get(0,j) );
    for( Int k=1; k<m; ++k )
    {
        for( Int j=0; j<n; ++j )        
        {
            t2.Set( j,   0, X.Get(k,j)                                    );
            X.Set(  k-1, j, C.Get(k,j)*t1.Get(j,0)+S.Get(k,j)*t2.Get(j,0) );
            t1.Set( j,   0, C.Get(k,j)*t2.Get(j,0)-
                            Conj(S.Get(k,j))*t1.Get(j,0) );
        }
    }
    for( Int j=0; j<n; ++j )
        X.Set( m-1, j, t1.Get(j,0) );
}

// NOTE: A [VC,* ] distribution might be most appropriate for the 
//       upper-Hessenberg matrix H since whole columns will need to be formed 
//       on every process and this distribution will keep the communication 
//       balanced.
template<typename F,Dist UH,Dist VH,Dist VX>
inline void
UN
( F alpha, const DistMatrix<F,UH,VH>& H, const DistMatrix<F,VX,STAR>& shifts,
  DistMatrix<F,STAR,VX>& X ) 
{
    DEBUG_ONLY(
        CallStackEntry cse("mshs::UN");
        if( shifts.ColAlign() != X.RowAlign() )
            LogicError("shifts and X are not aligned");
    )
    Scale( alpha, X );

    const Int m = X.Height();
    const Int n = X.Width();
    const Int nLoc = X.LocalWidth();
    if( m == 0 )
        return;

    // Initialize storage for Givens rotations
    typedef Base<F> Real;
    const Grid& g = H.Grid();
    DistMatrix<Real,STAR,VX> C(g);
    DistMatrix<F,STAR,VX> S(g);
    C.AlignWith( X );
    S.AlignWith( X );
    C.Resize( m, n );
    S.Resize( m, n );

    // Initialize the workspace for shifted columns of H
    DistMatrix<F,STAR,VX> W(m,n,g);
    {
        auto hLast = LockedView( H, 0, m-1, m, 1 );
        DistMatrix<F,STAR,STAR> hLast_STAR_STAR( hLast );
        for( Int jLoc=0; jLoc<nLoc; ++jLoc )
        {
            MemCopy( W.Buffer(0,jLoc), hLast_STAR_STAR.LockedBuffer(), m );
            W.UpdateLocal( m-1, jLoc, -shifts.GetLocal(jLoc,0) );
        }
    }
     
    // Simulataneously factor and solve against R
    DistMatrix<F,STAR,STAR> hT_STAR_STAR(g);
    for( Int k=m-1; k>0; --k )
    {
        auto hT = LockedView( H, 0, k-1, k-1, 1 );
        hT_STAR_STAR = hT;
        const F etakkm1 = H.Get(k,k-1);
        const F etakm1km1 = H.Get(k-1,k-1);
        for( Int jLoc=0; jLoc<nLoc; ++jLoc )
        {
            Real c;
            F s, rho;
            lapack::ComputeGivens( W.GetLocal(k,jLoc), etakkm1, &c, &s, &rho );
            C.SetLocal( k, jLoc, c );
            S.SetLocal( k, jLoc, s );

            X.SetLocal
            ( k, jLoc, X.GetLocal(k,jLoc)/(c*W.GetLocal(k,jLoc)+s*etakkm1) );

            auto xT = View( X.Matrix(), 0, jLoc, k-1, 1 );
            auto wT = View( W.Matrix(), 0, jLoc, k-1, 1 );

            const F xc = X.GetLocal(k,jLoc)*c;
            const F xs  = X.GetLocal(k,jLoc)*s;
            Axpy( -xc, wT, xT );
            Axpy( -xs, hT_STAR_STAR.LockedMatrix(), xT );

            Scale( -Conj(s), wT );
            Axpy( c, hT_STAR_STAR.LockedMatrix(), wT );

            const F mu = shifts.GetLocal( jLoc, 0 );
            X.UpdateLocal
            ( k-1, jLoc, -xc*W.GetLocal(k-1,jLoc)-xs*(etakm1km1-mu) );
            W.SetLocal
            ( k-1, jLoc, c*(etakm1km1-mu)-Conj(s)*W.GetLocal(k-1,jLoc) );
        }
    }
    for( Int jLoc=0; jLoc<nLoc; ++jLoc )
        X.SetLocal( 0, jLoc, X.GetLocal(0,jLoc)/W.GetLocal(0,jLoc) );

    // Solve against Q
    Matrix<F> t1(nLoc,1), t2(nLoc,1);
    for( Int jLoc=0; jLoc<nLoc; ++jLoc )
        t1.Set( jLoc, 0, X.GetLocal(0,jLoc) );
    for( Int k=1; k<m; ++k )
    {
        for( Int jLoc=0; jLoc<nLoc; ++jLoc )
        {
            t2.Set( jLoc, 0, X.GetLocal(k,jLoc) );
            X.SetLocal( k-1, jLoc, C.GetLocal(k,jLoc)*t1.Get(jLoc,0)+
                                   S.GetLocal(k,jLoc)*t2.Get(jLoc,0) );
            t1.Set( jLoc, 0, C.GetLocal(k,jLoc)*t2.Get(jLoc,0)-
                             Conj(S.GetLocal(k,jLoc))*t1.Get(jLoc,0) );
        }
    }
    for( Int jLoc=0; jLoc<nLoc; ++jLoc )
        X.SetLocal( m-1, jLoc, t1.Get(jLoc,0) );
}

// TODO: UT, LN, and LT

} // namespace mshs

template<typename F>
inline void
MultiShiftHessSolve
( UpperOrLower uplo, Orientation orientation,
  F alpha, const Matrix<F>& H, const Matrix<F>& shifts, Matrix<F>& X )
{
    DEBUG_ONLY(CallStackEntry cse("MultiShiftHessSolve"))
    if( uplo == UPPER )
    {
        if( orientation == NORMAL )
            mshs::UN( alpha, H, shifts, X );
        else
            LogicError("This option is not yet supported");
    }
    else
        LogicError("This option is not yet supported");
}

template<typename F>
inline void
MultiShiftHessSolve
( UpperOrLower uplo, Orientation orientation,
  F alpha, const DistMatrix<F,VC,STAR>& H, const DistMatrix<F,VR,STAR>& shifts, 
  DistMatrix<F>& X )
{
    DEBUG_ONLY(CallStackEntry cse("MultiShiftHessSolve"))
    if( uplo == UPPER )
    {
        if( orientation == NORMAL )
            mshs::UN( alpha, H, shifts, X );
        else
            LogicError("This option is not yet supported");
    }
    else
        LogicError("This option is not yet supported");
}

} // namespace elem

#endif // ifndef ELEM_MULTISHIFTHESSSOLVE_HPP
