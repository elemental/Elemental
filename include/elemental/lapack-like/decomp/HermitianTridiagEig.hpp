/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_HERMITIANTRIDIAGEIG_HPP
#define ELEM_HERMITIANTRIDIAGEIG_HPP

#include ELEM_HERMITIANTRIDIAGEIG_SORT_INC

// NOTE: eReal and ZReal could be packed into their complex counterparts

namespace elem {

// Return the full set of eigenvalues
// ==================================

template<typename Real>
inline void 
HermitianTridiagEig
( Matrix<Real>& d, Matrix<Real>& e, Matrix<Real>& w, SortType sort )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianTridiagEig"))
    const Int n = d.Height();
    w.Resize( n, 1 );
    lapack::SymmetricTridiagEig( n, d.Buffer(), e.Buffer(), w.Buffer() );
    Sort( w, sort );
}

template<typename Real>
inline void 
HermitianTridiagEig
( Matrix<Real>& d, Matrix<Complex<Real>>& e, Matrix<Real>& w, SortType sort )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianTridiagEig"))
    typedef Complex<Real> C;
    const Int n = d.Height();
    Matrix<Real> eReal( n-1, 1 );
    C yLast = 1;
    for( Int j=0; j<n-1; ++j )
    {
        const C psi = e.Get(j,0);
        const Real psiAbs = Abs(psi);
        if( psiAbs == Real(0) )
            yLast = 1;
        else
            yLast = Polar(Real(1),Arg(psi*yLast));
        eReal.Set( j, 0, psiAbs );
    }
    HermitianTridiagEig( d, eReal, w, sort );
}

template<typename Real,Dist U1,Dist V1,Dist U2,Dist V2,Dist U3>
inline void 
HermitianTridiagEig
( const DistMatrix<Real,U1,V1  >& d,
  const DistMatrix<Real,U2,V2  >& e,
        DistMatrix<Real,U3,STAR>& w, SortType sort )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianTridiagEig"))
    const Int n = d.Height();
    w.AlignCols( 0 );
    w.Resize( n, 1 );

    DistMatrix<Real,STAR,STAR> d_STAR_STAR( d );
    DistMatrix<Real,STAR,STAR> e_STAR_STAR( d.Grid() );
    e_STAR_STAR.Resize( n-1, 1, n );
    e_STAR_STAR = e;
    if( w.Participating() )
    {
        std::vector<Real> wVector(n);
        herm_tridiag_eig::Eig
        ( int(n), d_STAR_STAR.Buffer(), e_STAR_STAR.Buffer(), wVector.data(), 
          w.ColComm() );
        for( Int iLoc=0; iLoc<w.LocalHeight(); ++iLoc )
            w.SetLocal( iLoc, 0, wVector[iLoc] );
    }
    Sort( w, sort );
}

template<typename Real,Dist U1,Dist V1,Dist U2,Dist V2,Dist U3>
inline void 
HermitianTridiagEig
( const DistMatrix<Real,         U1,V1  >& d,
  const DistMatrix<Complex<Real>,U2,V2  >& e,
        DistMatrix<Real,         U3,STAR>& w, SortType sort )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianTridiagEig"))
    typedef Complex<Real> C;
    const Int n = d.Height();
    w.AlignCols( 0 );
    w.Resize( n, 1 );

    DistMatrix<Real,STAR,STAR> d_STAR_STAR( d );
    DistMatrix<C,STAR,STAR> e_STAR_STAR( e );

    DistMatrix<Real,STAR,STAR> eReal( n-1, 1, 0, 0, n, d.Grid() );
    if( w.Participating() )
    {
        C yLast = 1;
        for( Int j=0; j<n-1; ++j )
        {
            const C psi = e_STAR_STAR.GetLocal(j,0);
            const Real psiAbs = Abs(psi);
            if( psiAbs == Real(0) )
                yLast = 1;
            else
                yLast = Polar(Real(1),Arg(psi*yLast));
            eReal.SetLocal( j, 0, psiAbs );
        }

        std::vector<Real> wVector(n);
        herm_tridiag_eig::Eig
        ( int(n), d_STAR_STAR.Buffer(), eReal.Buffer(), wVector.data(), 
          w.ColComm() );
        for( Int iLoc=0; iLoc<w.LocalHeight(); ++iLoc )
            w.SetLocal( iLoc, 0, wVector[iLoc] );
    }
    Sort( w, sort );
}

// Return an index range of the eigenvalues
// ========================================

template<typename Real>
inline void 
HermitianTridiagEig
( Matrix<Real>& d, Matrix<Real>& e, Matrix<Real>& w, 
  Int il, Int iu, SortType sort )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianTridiagEig"))
    const Int n = d.Height();
    const Int k = ( n==0 ? 0 : iu-il+1 );
    w.Resize( n, 1 );
    lapack::SymmetricTridiagEig
    ( n, d.Buffer(), e.Buffer(), w.Buffer(), il, iu );
    w.Resize( k, 1 );
    Sort( w, sort );
}

template<typename Real>
inline void 
HermitianTridiagEig
( Matrix<Real>& d, Matrix<Complex<Real>>& e, Matrix<Real>& w, 
  Int il, Int iu, SortType sort )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianTridiagEig"))
    typedef Complex<Real> C;
    const Int n = d.Height();
    C yLast = 1;
    Matrix<Real> eReal( n-1, 1 );
    for( Int j=0; j<n-1; ++j )
    {
        const C psi = e.Get(j,0);
        const Real psiAbs = Abs(psi);
        if( psiAbs == Real(0) )
            yLast = 1;
        else
            yLast = Polar(Real(1),Arg(psi*yLast));
        eReal.Set( j, 0, psiAbs );
    }
    HermitianTridiagEig( d, eReal, w, il, iu, sort );
}

template<typename Real,Dist U1,Dist V1,Dist U2,Dist V2,Dist U3>
inline void 
HermitianTridiagEig
( const DistMatrix<Real,U1,V1  >& d,
  const DistMatrix<Real,U2,V2  >& e,
        DistMatrix<Real,U3,STAR>& w, Int il, Int iu, SortType sort )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianTridiagEig"))
    const Int n = d.Height();
    const Int k = ( n==0 ? 0 : iu-il+1 );
    w.AlignCols( 0 );
    w.Resize( k, 1 );

    DistMatrix<Real,STAR,STAR> d_STAR_STAR( d );
    DistMatrix<Real,STAR,STAR> e_STAR_STAR( d.Grid() );
    e_STAR_STAR.Resize( n-1, 1, n );
    e_STAR_STAR = e;
    if( w.Participating() )
    {
        std::vector<Real> wVector(n);
        herm_tridiag_eig::Eig
        ( int(n), d_STAR_STAR.Buffer(), e_STAR_STAR.Buffer(), wVector.data(), 
          w.ColComm(), int(il), int(iu) );
        for( Int iLoc=0; iLoc<w.LocalHeight(); ++iLoc )
            w.SetLocal( iLoc, 0, wVector[iLoc] );
    }
    Sort( w, sort );
}

template<typename Real,Dist U1,Dist V1,Dist U2,Dist V2,Dist U3>
inline void 
HermitianTridiagEig
( const DistMatrix<Real,         U1,V1  >& d,
  const DistMatrix<Complex<Real>,U2,V2  >& e,
        DistMatrix<Real,         U3,STAR>& w, Int il, Int iu, SortType sort )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianTridiagEig"))
    typedef Complex<Real> C;
    const Int n = d.Height();
    const Int k = ( n==0 ? 0 : iu-il+1 );
    w.AlignCols( 0 );
    w.Resize( k, 1 );

    DistMatrix<Real,STAR,STAR> d_STAR_STAR( d );
    DistMatrix<C,STAR,STAR> e_STAR_STAR( e );

    DistMatrix<Real,STAR,STAR> eReal( n-1, 1, 0, 0, n, d.Grid() );
    if( w.Participating() )
    {
        C yLast = 1;
        for( Int j=0; j<n-1; ++j )
        {
            const C psi = e_STAR_STAR.GetLocal(j,0);
            const Real psiAbs = Abs(psi);
            if( psiAbs == Real(0) )
                yLast = 1;
            else
                yLast = Polar(Real(1),Arg(psi*yLast));
            eReal.SetLocal( j, 0, psiAbs );
        }

        std::vector<Real> wVector(n);
        herm_tridiag_eig::Eig
        ( int(n), d_STAR_STAR.Buffer(), eReal.Buffer(), wVector.data(), 
          w.ColComm(), int(il), int(iu) );
        for( Int iLoc=0; iLoc<w.LocalHeight(); ++iLoc )
            w.SetLocal( iLoc, 0, wVector[iLoc] );
    }
    Sort( w, sort );
}

// Return the eigenvalues in a given interval
// ==========================================

template<typename Real>
inline void 
HermitianTridiagEig
( Matrix<Real>& d, Matrix<Real>& e, Matrix<Real>& w, 
  Real vl, Real vu, SortType sort )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianTridiagEig"))
    const Int n = d.Height();
    w.Resize( n, 1 );
    const Int k = lapack::SymmetricTridiagEig
    ( n, d.Buffer(), e.Buffer(), w.Buffer(), vl, vu );
    w.Resize( k, 1 );
    Sort( w, sort );
}

template<typename Real>
inline void 
HermitianTridiagEig
( Matrix<Real>& d, Matrix<Complex<Real>>& e, Matrix<Real>& w, 
  Real vl, Real vu, SortType sort )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianTridiagEig"))
    typedef Complex<Real> C;
    const Int n = d.Height();
    Matrix<Real> eReal( n-1, 1 );
    C yLast = 1;
    for( Int j=0; j<n-1; ++j )
    {
        const C psi = e.Get(j,0);
        const Real psiAbs = Abs(psi);
        if( psiAbs == Real(0) )
            yLast = 1;
        else
            yLast = Polar(Real(1),Arg(psi*yLast));
        eReal.Set( j, 0, psiAbs );
    }
    HermitianTridiagEig( d, eReal, w, vl, vu, sort );
}

template<typename Real,Dist U1,Dist V1,Dist U2,Dist V2,Dist U3>
inline void 
HermitianTridiagEig
( const DistMatrix<Real,U1,V1  >& d,
  const DistMatrix<Real,U2,V2  >& e,
        DistMatrix<Real,U3,STAR>& w, Real vl, Real vu, SortType sort )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianTridiagEig"))
    const Int n = d.Height();
    w.AlignCols( 0 );

    DistMatrix<Real,STAR,STAR> d_STAR_STAR( d );
    DistMatrix<Real,STAR,STAR> e_STAR_STAR( d.Grid() );
    e_STAR_STAR.Resize( n-1, 1, n );
    e_STAR_STAR = e;
    if( w.Participating() )
    {
        std::vector<Real> wVector(n);
        auto info = herm_tridiag_eig::Eig
        ( int(n), d_STAR_STAR.Buffer(), e_STAR_STAR.Buffer(), wVector.data(), 
          w.ColComm(), vl, vu );
        const Int k = info.numGlobalEigenvalues;
        w.Resize( k, 1 );
        for( Int iLoc=0; iLoc<w.LocalHeight(); ++iLoc )
            w.SetLocal( iLoc, 0, wVector[iLoc] );
    }
    w.MakeConsistent( true );
    Sort( w, sort );
}

template<typename Real,Dist U1,Dist V1,Dist U2,Dist V2,Dist U3>
inline void 
HermitianTridiagEig
( const DistMatrix<Real,         U1,V1  >& d,
  const DistMatrix<Complex<Real>,U2,V2  >& e,
        DistMatrix<Real,         U3,STAR>& w, Real vl, Real vu, SortType sort )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianTridiagEig"))
    typedef Complex<Real> C;
    const Int n = d.Height();
    w.AlignCols( 0 );

    DistMatrix<Real,STAR,STAR> d_STAR_STAR( d );
    DistMatrix<C,STAR,STAR> e_STAR_STAR( e );

    DistMatrix<Real,STAR,STAR> eReal( d.Grid() );
    eReal.Resize( n-1, 1, n );
    if( w.Participating() )
    {
        C yLast = 1;
        for( Int j=0; j<n-1; ++j )
        {
            const C psi = e_STAR_STAR.GetLocal(j,0);
            const Real psiAbs = Abs(psi);
            if( psiAbs == Real(0) )
                yLast = 1;
            else
                yLast = Polar(Real(1),Arg(psi*yLast));
            eReal.SetLocal( j, 0, psiAbs );
        }

        std::vector<Real> wVector(n);
        auto info = herm_tridiag_eig::Eig
        ( int(n), d_STAR_STAR.Buffer(), eReal.Buffer(), wVector.data(), 
          w.ColComm(), vl, vu );
        const Int k = info.numGlobalEigenvalues;
        w.Resize( k, 1 );
        for( Int iLoc=0; iLoc<w.LocalHeight(); ++iLoc )
            w.SetLocal( iLoc, 0, wVector[iLoc] );
    }
    w.MakeConsistent( true );
    Sort( w, sort );
}

// Return the full set of eigenpairs
// =================================

template<typename Real>
inline void 
HermitianTridiagEig
( Matrix<Real>& d, Matrix<Real>& e, Matrix<Real>& w, Matrix<Real>& Z, 
  SortType sort )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianTridiagEig"))
    const Int n = d.Height();
    w.Resize( n, 1 );
    Z.Resize( n, n );
    lapack::SymmetricTridiagEig
    ( n, d.Buffer(), e.Buffer(), w.Buffer(), Z.Buffer(), Z.LDim() );
    herm_eig::Sort( w, Z, sort );
}

template<typename Real>
inline void 
HermitianTridiagEig
( Matrix<Real>& d, Matrix<Complex<Real>>& e, Matrix<Real>& w, 
  Matrix<Complex<Real>>& Z, SortType sort )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianTridiagEig"))
    typedef Complex<Real> C;
    const Int n = d.Height();
    Matrix<Real> eReal( n-1, 1 );
    Matrix<C> y( n, 1 );
    y.Set( 0, 0, 1 );
    for( Int j=0; j<n-1; ++j )
    {
        const C psi = e.Get(j,0);
        const Real psiAbs = Abs(psi);
        if( psiAbs == Real(0) )
            y.Set( j+1, 0, 1 );
        else
            y.Set( j+1, 0, Polar(Real(1),Arg(psi*y.Get(j,0))) );
        eReal.Set( j, 0, psiAbs );
    }
    Matrix<Real> ZReal;
    HermitianTridiagEig( d, eReal, w, ZReal, sort );
    Z.Resize( n, n );
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<n; ++i )
            Z.Set( i, j, y.Get(i,0)*ZReal.Get(i,j) );
}

template<typename Real,Dist U1,Dist V1,Dist U2,Dist V2,Dist U3>
inline void 
HermitianTridiagEig
( const DistMatrix<Real,U1,  V1  >& d,
  const DistMatrix<Real,U2,  V2  >& e,
        DistMatrix<Real,U3,  STAR>& w, 
        DistMatrix<Real,STAR,U3  >& Z, SortType sort )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianTridiagEig"))
    const Int n = d.Height();
    w.AlignCols( 0 );
    w.Resize( n, 1 );
    Z.AlignRows( 0 );
    Z.Resize( n, n );

    DistMatrix<Real,STAR,STAR> d_STAR_STAR( d );
    DistMatrix<Real,STAR,STAR> e_STAR_STAR( d.Grid() );
    e_STAR_STAR.Resize( n-1, 1, n );
    e_STAR_STAR = e;
    if( w.Participating() )
    {
        std::vector<Real> wVector(n);
        herm_tridiag_eig::Eig
        ( int(n), d_STAR_STAR.Buffer(), e_STAR_STAR.Buffer(), wVector.data(), 
          Z.Buffer(), Z.LDim(), w.ColComm() );
        for( Int iLoc=0; iLoc<w.LocalHeight(); ++iLoc )
            w.SetLocal( iLoc, 0, wVector[iLoc] );
    }
    herm_eig::Sort( w, Z, sort );
}

template<typename Real,Dist U1,Dist V1,Dist U2,Dist V2,Dist U3>
inline void 
HermitianTridiagEig
( const DistMatrix<Real,         U1,  V1  >& d,
  const DistMatrix<Complex<Real>,U2,  V2  >& e,
        DistMatrix<Real,         U3,  STAR>& w, 
        DistMatrix<Complex<Real>,STAR,U3  >& Z, SortType sort )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianTridiagEig"))
    typedef Complex<Real> C;
    const Int n = d.Height();
    w.AlignCols( 0 );
    w.Resize( n, 1 );
    Z.AlignRows( 0 );
    Z.Resize( n, n );

    const Grid& g = d.Grid();
    DistMatrix<Real,STAR,STAR> d_STAR_STAR( d );
    DistMatrix<C,STAR,STAR> e_STAR_STAR( e );

    DistMatrix<Real,STAR,STAR> eReal( g );
    eReal.Resize( n-1, 1, n );
    DistMatrix<Real,STAR,U3  > ZReal(g);
    DistMatrix<C,   STAR,STAR> y(n,1,g);
    if( w.Participating() )
    {
        y.SetLocal(0,0,1);
        for( Int j=0; j<n-1; ++j )
        {
            const C psi = e_STAR_STAR.GetLocal(j,0);
            const Real psiAbs = Abs(psi);
            if( psiAbs == Real(0) )
                y.SetLocal( j+1, 0, 1 );
            else
                y.SetLocal( j+1, 0, Polar(Real(1),Arg(psi*y.GetLocal(j,0))) );
            eReal.SetLocal( j, 0, psiAbs );
        }

        std::vector<Real> wVector(n);
        herm_tridiag_eig::Eig
        ( int(n), d_STAR_STAR.Buffer(), eReal.Buffer(), wVector.data(), 
          ZReal.Buffer(), ZReal.LDim(), w.ColComm() );
        for( Int iLoc=0; iLoc<w.LocalHeight(); ++iLoc )
            w.SetLocal( iLoc, 0, wVector[iLoc] );
    }
    herm_eig::Sort( w, ZReal, sort );
    Z.Resize( n, n );
    for( Int jLoc=0; jLoc<Z.LocalWidth(); ++jLoc )
        for( Int i=0; i<n; ++i )
            Z.SetLocal( i, jLoc, y.GetLocal(i,0)*ZReal.GetLocal(i,jLoc) );
}

template<typename F,Dist U1,Dist V1,Dist U2,Dist V2,Dist U3,Dist U4,Dist V4>
inline void 
HermitianTridiagEig
( const DistMatrix<Base<F>,U1,V1  >& d,
  const DistMatrix<F,      U2,V2  >& e,
        DistMatrix<Base<F>,U3,STAR>& w, 
        DistMatrix<F,      U4,V4  >& Z, SortType sort )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianTridiagEig"))
    DistMatrix<F,STAR,U3> Z_STAR_U3( Z );
    HermitianTridiagEig( d, e, w, Z_STAR_U3, sort );
    Z = Z_STAR_U3;
}

// Return an index range of eigenpairs
// ===================================

template<typename Real>
inline void 
HermitianTridiagEig
( Matrix<Real>& d, Matrix<Real>& e, Matrix<Real>& w, 
  Matrix<Real>& Z, Int il, Int iu, SortType sort )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianTridiagEig"))
    const Int n = d.Height();
    const Int k = ( n==0 ? 0 : iu-il+1 );
    w.Resize( n, 1 );
    Z.Resize( n, n );
    lapack::SymmetricTridiagEig
    ( n, d.Buffer(), e.Buffer(), w.Buffer(), Z.Buffer(), Z.LDim(), il, iu );
    w.Resize( k, 1 );
    herm_eig::Sort( w, Z, sort );
}

template<typename Real>
inline void 
HermitianTridiagEig
( Matrix<Real>& d, Matrix<Complex<Real>>& e, Matrix<Real>& w, 
  Matrix<Complex<Real>>& Z, Int il, Int iu, SortType sort )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianTridiagEig"))
    typedef Complex<Real> C;
    const Int n = d.Height();
    const Int k = ( n==0 ? 0 : iu-il+1 );

    Matrix<Real> eReal( n-1, 1 );
    Matrix<C> y( n, 1 );
    y.Set( 0, 0, 1 );
    for( Int j=0; j<n-1; ++j )
    {
        const C psi = e.Get(j,0);
        const Real psiAbs = Abs(psi);
        if( psiAbs == Real(0) )
            y.Set( j+1, 0, 1 );
        else
            y.Set( j+1, 0, Polar(Real(1),Arg(psi*y.Get(j,0))) );
        eReal.Set( j, 0, psiAbs );
    }

    Matrix<Real> ZReal;
    HermitianTridiagEig( d, eReal, w, ZReal, il, iu, sort );

    Z.Resize( n, k );
    for( Int j=0; j<k; ++j )
        for( Int i=0; i<n; ++i )
            Z.Set( i, j, y.Get(i,0)*ZReal.Get(i,j) );
}

template<typename Real,Dist U1,Dist V1,Dist U2,Dist V2,Dist U3>
inline void 
HermitianTridiagEig
( const DistMatrix<Real,U1,  V1  >& d,
  const DistMatrix<Real,U2,  V2  >& e,
        DistMatrix<Real,U3,  STAR>& w, 
        DistMatrix<Real,STAR,U3  >& Z, Int il, Int iu, SortType sort )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianTridiagEig"))
    const Int n = d.Height();
    const Int k = ( n==0 ? 0 : iu-il+1 );
    w.AlignCols( 0 );
    w.Resize( k, 1 );
    Z.AlignRows( 0 );
    Z.Resize( n, k );

    DistMatrix<Real,STAR,STAR> d_STAR_STAR( d );
    DistMatrix<Real,STAR,STAR> e_STAR_STAR( d.Grid() );
    e_STAR_STAR.Resize( n-1, 1, n );
    e_STAR_STAR = e;
    if( w.Participating() )
    {
        std::vector<Real> wVector(n);
        herm_tridiag_eig::Eig
        ( int(n), d_STAR_STAR.Buffer(), e_STAR_STAR.Buffer(), wVector.data(), 
          Z.Buffer(), Z.LDim(), w.ColComm(), int(il), int(iu) );
        for( Int iLoc=0; iLoc<w.LocalHeight(); ++iLoc )
            w.SetLocal( iLoc, 0, wVector[iLoc] );
    }
    herm_eig::Sort( w, Z, sort );
}

template<typename Real,Dist U1,Dist V1,Dist U2,Dist V2,Dist U3>
inline void 
HermitianTridiagEig
( const DistMatrix<Real,         U1,  V1  >& d,
  const DistMatrix<Complex<Real>,U2,  V2  >& e,
        DistMatrix<Real,         U3,  STAR>& w, 
        DistMatrix<Complex<Real>,STAR,U3  >& Z, Int il, Int iu, SortType sort )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianTridiagEig"))
    typedef Complex<Real> C;
    const Int n = d.Height();
    const Int k = ( n==0 ? 0 : iu-il+1 );
    w.AlignCols( 0 );
    w.Resize( k, 1 );
    Z.AlignRows( 0 );
    Z.Resize( n, k );

    const Grid& g = d.Grid();
    DistMatrix<Real,STAR,STAR> d_STAR_STAR( d );
    DistMatrix<C,STAR,STAR> e_STAR_STAR( e );

    DistMatrix<Real,STAR,STAR> eReal( g );
    eReal.Resize( n-1, 1, n );
    DistMatrix<Real,STAR,U3  > ZReal(g);
    DistMatrix<C,   STAR,STAR> y(n,1,g);
    if( w.Participating() )
    {
        y.SetLocal(0,0,1);
        for( Int j=0; j<n-1; ++j )
        {
            const C psi = e_STAR_STAR.GetLocal(j,0);
            const Real psiAbs = Abs(psi);
            if( psiAbs == Real(0) )
                y.SetLocal( j+1, 0, 1 );
            else
                y.SetLocal( j+1, 0, Polar(Real(1),Arg(psi*y.GetLocal(j,0))) );
            eReal.SetLocal( j, 0, psiAbs );
        }

        std::vector<Real> wVector(n);
        herm_tridiag_eig::Eig
        ( int(n), d_STAR_STAR.Buffer(), eReal.Buffer(), wVector.data(), 
          ZReal.Buffer(), ZReal.LDim(), w.ColComm(), int(il), int(iu) );
        for( Int iLoc=0; iLoc<w.LocalHeight(); ++iLoc )
            w.SetLocal( iLoc, 0, wVector[iLoc] );
    }
    herm_eig::Sort( w, ZReal, sort );
    for( Int jLoc=0; jLoc<ZReal.Width(); ++jLoc )
        for( Int i=0; i<n; ++i )
            Z.SetLocal( i, jLoc, y.GetLocal(i,0)*ZReal.GetLocal(i,jLoc) );
}

template<typename F,Dist U1,Dist V1,Dist U2,Dist V2,Dist U3,Dist U4,Dist V4>
inline void 
HermitianTridiagEig
( const DistMatrix<Base<F>,U1,V1  >& d,
  const DistMatrix<F,      U2,V2  >& e,
        DistMatrix<Base<F>,U3,STAR>& w, 
        DistMatrix<F,      U4,V4  >& Z, Int il, Int iu, SortType sort )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianTridiagEig"))
    DistMatrix<F,STAR,U3> Z_STAR_U3( Z );
    HermitianTridiagEig( d, e, w, Z_STAR_U3, il, iu, sort );
    Z = Z_STAR_U3;
}

// Return the eigenpairs with eigenvalues in a given interval
// ==========================================================

template<typename Real>
inline void 
HermitianTridiagEig
( Matrix<Real>& d, Matrix<Real>& e, Matrix<Real>& w, Matrix<Real>& Z,
  Real vl, Real vu, SortType sort )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianTridiagEig"))
    const Int n = d.Height();
    w.Resize( n, 1 );
    Z.Resize( n, n ); // This can be an unnecessary O(n^2) memory usage
    const Int k = lapack::SymmetricTridiagEig
    ( n, d.Buffer(), e.Buffer(), w.Buffer(), Z.Buffer(), Z.LDim(), vl, vu );
    w.Resize( k, 1 );
    Z.Resize( n, k );
    herm_eig::Sort( w, Z, sort );
}

template<typename Real>
inline void 
HermitianTridiagEig
( Matrix<Real>& d, Matrix<Complex<Real>>& e, Matrix<Real>& w, 
  Matrix<Complex<Real>>& Z, Real vl, Real vu, SortType sort )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianTridiagEig"))
    typedef Complex<Real> C;
    const Int n = d.Height();
    Matrix<Real> eReal( n-1, 1 );
    Matrix<C> y( n, 1 );
    y.Set( 0, 0, 1 );
    for( Int j=0; j<n-1; ++j )
    {
        const C psi = e.Get(j,0);
        const Real psiAbs = Abs(psi);
        if( psiAbs == Real(0) )
            y.Set( j+1, 0, 1 );
        else
            y.Set( j+1, 0, Polar(Real(1),Arg(psi*y.Get(j,0))) );
        eReal.Set( j, 0, psiAbs );
    }
    Matrix<Real> ZReal;
    HermitianTridiagEig( d, eReal, w, ZReal, vl, vu, sort );
    const Int k = w.Height();
    Z.Resize( n, k );
    for( Int j=0; j<k; ++j )
        for( Int i=0; i<n; ++i )
            Z.Set( i, j, y.Get(i,0)*ZReal.Get(i,j) );
}

template<typename Real,Dist U1,Dist V1,Dist U2,Dist V2,Dist U3>
inline void 
HermitianTridiagEig
( const DistMatrix<Real,U1,V1  >& d,
  const DistMatrix<Real,U2,V2  >& e,
        DistMatrix<Real,U3,STAR>& w, 
        DistMatrix<Real,STAR,U3>& Z, Real vl, Real vu, SortType sort )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianTridiagEig"))
    const Int n = d.Height();
    w.AlignCols( 0 );
    Z.AlignRows( 0 );

    DistMatrix<Real,STAR,STAR> d_STAR_STAR( d );
    DistMatrix<Real,STAR,STAR> e_STAR_STAR( d.Grid() );
    e_STAR_STAR.Resize( n-1, 1, n );
    e_STAR_STAR = e;
    if( w.Participating() )
    {
        std::vector<Real> dVector(n), eVector(n), wVector(n);
        MemCopy( dVector.data(), d_STAR_STAR.Buffer(), n );
        MemCopy( eVector.data(), e_STAR_STAR.Buffer(), n-1 );
        auto estimate = herm_tridiag_eig::EigEstimate
        ( int(n), dVector.data(), eVector.data(), wVector.data(), w.ColComm(),
          vl, vu );
        SwapClear( dVector );
        SwapClear( eVector );
        const Int kEst = estimate.numGlobalEigenvalues;
        Z.Resize( n, kEst );

        auto info = herm_tridiag_eig::Eig
        ( int(n), d_STAR_STAR.Buffer(), e_STAR_STAR.Buffer(), wVector.data(), 
          Z.Buffer(), Z.LDim(), w.ColComm(), vl, vu );
        const Int k = info.numGlobalEigenvalues;

        w.Resize( k, 1 );
        for( Int iLoc=0; iLoc<w.LocalHeight(); ++iLoc )
            w.SetLocal( iLoc, 0, wVector[iLoc] );
        Z.Resize( n, k );
    }
    w.MakeConsistent( true );
    Z.MakeConsistent( true );
    herm_eig::Sort( w, Z, sort );
}

template<typename Real,Dist U1,Dist V1,Dist U2,Dist V2,Dist U3>
inline void 
HermitianTridiagEig
( const DistMatrix<Real,         U1,  V1  >& d,
  const DistMatrix<Complex<Real>,U2,  V2  >& e,
        DistMatrix<Real,         U3,  STAR>& w, 
        DistMatrix<Complex<Real>,STAR,U3  >& Z, 
  Real vl, Real vu, SortType sort )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianTridiagEig"))
    typedef Complex<Real> C;
    const Int n = d.Height();
    w.AlignCols( 0 );
    Z.AlignRows( 0 );

    const Grid& g = d.Grid();
    DistMatrix<Real,STAR,STAR> d_STAR_STAR( d );
    DistMatrix<C,STAR,STAR> e_STAR_STAR( e );

    DistMatrix<Real,STAR,STAR> eReal( g );
    eReal.Resize( n-1, 1, n );
    DistMatrix<Real,STAR,U3  > ZReal(g);
    DistMatrix<C,   STAR,STAR> y(n,1,g);
    if( w.Participating() )
    {
        y.SetLocal(0,0,1);
        for( Int j=0; j<n-1; ++j )
        {
            const C psi = e_STAR_STAR.GetLocal(j,0);
            const Real psiAbs = Abs(psi);
            if( psiAbs == Real(0) )
                y.SetLocal( j+1, 0, 1 );
            else
                y.SetLocal( j+1, 0, Polar(Real(1),Arg(psi*y.GetLocal(j,0))) );
            eReal.SetLocal( j, 0, psiAbs );
        }

        std::vector<Real> dVector(n), eVector(n), wVector(n);
        MemCopy( dVector.data(), d_STAR_STAR.Buffer(), n );
        MemCopy( eVector.data(), eReal.Buffer(), n-1 );
        auto estimate = herm_tridiag_eig::EigEstimate
        ( int(n), dVector.data(), eVector.data(), wVector.data(), w.ColComm(),
          vl, vu );
        SwapClear( dVector );
        SwapClear( eVector );
        const Int kEst = estimate.numGlobalEigenvalues;
        ZReal.Resize( n, kEst );

        auto info = herm_tridiag_eig::Eig
        ( int(n), d_STAR_STAR.Buffer(), eReal.Buffer(), wVector.data(), 
          ZReal.Buffer(), ZReal.LDim(), w.ColComm(), vl, vu );
        const Int k = info.numGlobalEigenvalues;

        w.Resize( k, 1 );
        for( Int iLoc=0; iLoc<w.LocalHeight(); ++iLoc )
            w.SetLocal( iLoc, 0, wVector[iLoc] );
        ZReal.Resize( n, k );
        Z.Resize( n, k );
    }
    w.MakeConsistent( true );
    Z.MakeConsistent( true );
    ZReal.MakeConsistent( true );
    herm_eig::Sort( w, ZReal, sort );
    for( Int jLoc=0; jLoc<Z.LocalWidth(); ++jLoc ) 
        for( Int i=0; i<n; ++i )
            Z.SetLocal( i, jLoc, y.GetLocal(i,0)*ZReal.GetLocal(i,jLoc) );
}

template<typename Real,Dist U1,Dist V1,Dist U2,Dist V2,Dist U3,Dist U4,Dist V4>
inline void 
HermitianTridiagEig
( const DistMatrix<Real,U1,V1  >& d,
  const DistMatrix<Real,U2,V2  >& e,
        DistMatrix<Real,U3,STAR>& w, 
        DistMatrix<Real,U4,V4  >& Z, Real vl, Real vu, SortType sort )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianTridiagEig"))
    DistMatrix<Real,STAR,U3> Z_STAR_U3( Z );
    HermitianTridiagEig( d, e, w, Z_STAR_U3, vl, vu, sort );
    Z = Z_STAR_U3;
}

template<typename F,Dist U1,Dist V1,Dist U2,Dist V2,Dist U3,Dist U4,Dist V4>
inline void 
HermitianTridiagEig
( const DistMatrix<Base<F>,U1,V1  >& d,
  const DistMatrix<F,      U2,V2  >& e,
        DistMatrix<Base<F>,U3,STAR>& w, 
        DistMatrix<F,      U4,V4  >& Z, Base<F> vl, Base<F> vu, SortType sort )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianTridiagEig"))
    DistMatrix<F,STAR,U3> Z_STAR_U3( Z );
    HermitianTridiagEig( d, e, w, Z_STAR_U3, vl, vu, sort );
    Z = Z_STAR_U3;
}

template<typename Real,Dist U1,Dist V1,Dist U2,Dist V2>
inline Int
HermitianTridiagEigEstimate
( const DistMatrix<Real,U1,V1  >& d,
  const DistMatrix<Real,U2,V2  >& e,
        mpi::Comm wColComm, Real vl, Real vu )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianTridiagEigEstimate"))
    const Int n = d.Height();
    DistMatrix<Real,STAR,STAR> d_STAR_STAR( d );
    DistMatrix<Real,STAR,STAR> e_STAR_STAR( d.Grid() );
    e_STAR_STAR.Resize( n-1, 1, n );
    e_STAR_STAR = e;
    std::vector<Real> dVector(n), eVector(n), wVector(n);
    MemCopy( dVector.data(), d_STAR_STAR.Buffer(), n );
    MemCopy( eVector.data(), e_STAR_STAR.Buffer(), n-1 );
    auto estimate = herm_tridiag_eig::EigEstimate
    ( int(n), dVector.data(), eVector.data(), wVector.data(), wColComm,
      vl, vu );
    return estimate.numGlobalEigenvalues;
}

// Z is assumed to be sufficiently large and properly aligned
template<typename Real,Dist U1,Dist V1,Dist U2,Dist V2,Dist U3>
inline void 
HermitianTridiagEigPostEstimate
( const DistMatrix<Real,U1,V1  >& d,
  const DistMatrix<Real,U2,V2  >& e,
        DistMatrix<Real,U3,STAR>& w, 
        DistMatrix<Real,STAR,U3>& Z, Real vl, Real vu, SortType sort )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianTridiagEigPostEstimate"))
    const Int n = d.Height();
    w.AlignCols( 0 );
    if( Z.RowAlign() != 0 )
        LogicError("Z was not properly aligned");

    DistMatrix<Real,STAR,STAR> d_STAR_STAR( d );
    DistMatrix<Real,STAR,STAR> e_STAR_STAR( d.Grid() );
    e_STAR_STAR.Resize( n-1, 1, n );
    e_STAR_STAR = e;
    if( w.Participating() )
    {
        std::vector<Real> wVector(n);
        auto info = herm_tridiag_eig::Eig
        ( int(n), d_STAR_STAR.Buffer(), e_STAR_STAR.Buffer(), wVector.data(), 
          Z.Buffer(), Z.LDim(), w.ColComm(), vl, vu );
        const Int k = info.numGlobalEigenvalues;

        w.Resize( k, 1 );
        for( Int iLoc=0; iLoc<w.LocalHeight(); ++iLoc )
            w.SetLocal( iLoc, 0, wVector[iLoc] );
        Z.Resize( n, k );
    }
    w.MakeConsistent( true );
    Z.MakeConsistent( true );
    herm_eig::Sort( w, Z, sort );
}

} // namespace elem

#endif // ifndef ELEM_HERMITIANTRIDIAGEIG_HPP
