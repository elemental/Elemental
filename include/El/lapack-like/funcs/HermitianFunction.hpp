/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_HERMITIANFUNCTION_HPP
#define EL_HERMITIANFUNCTION_HPP

#include EL_HERMITIANFROMEVD_INC
#include EL_NORMALFROMEVD_INC

namespace El {

// Modify the eigenvalues of A with the real-valued function f, which will 
// therefore result in a Hermitian matrix, which we store in-place.

template<typename F,class RealFunction>
inline void
RealHermitianFunction( UpperOrLower uplo, Matrix<F>& A, RealFunction func )
{
    DEBUG_ONLY(CallStackEntry cse("RealHermitianFunction"))
    if( A.Height() != A.Width() )
        LogicError("Hermitian matrices must be square");
    typedef Base<F> Real;

    // Get the EVD of A
    Matrix<Real> w;
    Matrix<F> Z;
    HermitianEig( uplo, A, w, Z );

    // Replace w with f(w)
    const Int n = w.Height();
    for( Int i=0; i<n; ++i )
    {
        const Real omega = w.Get(i,0);
        w.Set(i,0,func(omega));
    }

    // A := Z f(Omega) Z^H
    HermitianFromEVD( uplo, A, w, Z );
}

template<typename F,class RealFunction>
inline void
RealHermitianFunction( UpperOrLower uplo, DistMatrix<F>& A, RealFunction func )
{
    DEBUG_ONLY(CallStackEntry cse("RealHermitianFunction"))
    if( A.Height() != A.Width() )
        LogicError("Hermitian matrices must be square");
    typedef Base<F> Real;

    // Get the EVD of A
    const Grid& g = A.Grid();
    DistMatrix<Real,VR,STAR> w(g);
    DistMatrix<F> Z(g);
    HermitianEig( uplo, A, w, Z );

    // Replace w with f(w)
    const Int numLocalEigs = w.LocalHeight();
    for( Int iLoc=0; iLoc<numLocalEigs; ++iLoc )
    {
        const Real omega = w.GetLocal(iLoc,0);
        w.SetLocal(iLoc,0,func(omega));
    }

    // A := Z f(Omega) Z^H
    HermitianFromEVD( uplo, A, w, Z ); 
}

// Modify the eigenvalues of A with the complex-valued function f, which will
// therefore result in a normal (in general, non-Hermitian) matrix, which we 
// store in-place. At some point a version will be written which takes a real
// symmetric matrix as input and produces a complex normal matrix.

template<typename Real,class Function>
inline void
ComplexHermitianFunction
( UpperOrLower uplo, Matrix<Complex<Real>>& A, Function func )
{
    DEBUG_ONLY(CallStackEntry cse("ComplexHermitianFunction"))
    if( A.Height() != A.Width() )
        LogicError("Hermitian matrices must be square");
    typedef Complex<Real> C;

    // Get the EVD of A
    Matrix<Real> w;
    Matrix<C> Z;
    HermitianEig( uplo, A, w, Z );

    // Form f(w)
    const Int n = w.Height();
    Matrix<C> fw( n, 1 );
    for( Int i=0; i<n; ++i )
    {
        const Real omega = w.Get(i,0);
        fw.Set(i,0,func(omega));
    }

    // A := Z f(Omega) Z^H
    NormalFromEVD( A, fw, Z );
}

template<typename Real,class Function>
inline void
ComplexHermitianFunction
( UpperOrLower uplo, DistMatrix<Complex<Real>>& A, Function func )
{
    DEBUG_ONLY(CallStackEntry cse("ComplexHermitianFunction"))
    if( A.Height() != A.Width() )
        LogicError("Hermitian matrices must be square");
    typedef Complex<Real> C;

    // Get the EVD of A
    const Grid& g = A.Grid();
    DistMatrix<Real,VR,STAR> w(g);
    DistMatrix<C> Z(g);
    HermitianEig( uplo, A, w, Z );

    // Form f(w)
    DistMatrix<C,VR,STAR> fw(g);
    fw.AlignWith( w.DistData() );
    fw.Resize( w.Height(), 1 );
    const Int numLocalEigs = w.LocalHeight();
    for( Int iLoc=0; iLoc<numLocalEigs; ++iLoc )
    {
        const Real omega = w.GetLocal(iLoc,0);
        fw.SetLocal(iLoc,0,func(omega));
    }

    // A := Z f(Omega) Z^H
    NormalFromEVD( A, fw, Z );
}

} // namespace El

#endif // ifndef EL_HERMITIANFUNCTION_HPP
