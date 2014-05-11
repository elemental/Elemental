/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_SCHUR_HPP
#define ELEM_SCHUR_HPP

#include "./Schur/CheckReal.hpp"
#include "./Schur/RealToComplex.hpp"
#include "./Schur/QuasiTriangEig.hpp"
#include "./Schur/QR.hpp"
#include "./Schur/SDC.hpp"
#include "./Schur/InverseFreeSDC.hpp"

namespace elem {

template<typename F>
inline void
Schur( Matrix<F>& A, Matrix<Complex<Base<F>>>& w )
{
    DEBUG_ONLY(CallStackEntry cse("Schur"))
    schur::QR( A, w );
}

template<typename F>
inline void
Schur( Matrix<F>& A, Matrix<Complex<Base<F>>>& w, Matrix<F>& Q )
{
    DEBUG_ONLY(CallStackEntry cse("Schur"))
    schur::QR( A, w, Q );
}

template<typename F>
inline void
Schur( DistMatrix<F>& A, DistMatrix<Complex<Base<F>>,VR,STAR>& w )
{
    DEBUG_ONLY(CallStackEntry cse("Schur"))
#ifdef ELEM_HAVE_SCALAPACK
    schur::QR( A, w );
#else
    schur::SDC( A, w );
#endif
}

template<typename F>
inline void
Schur
( DistMatrix<F>& A, DistMatrix<Complex<Base<F>>,VR,STAR>& w, DistMatrix<F>& Q )
{
    DEBUG_ONLY(CallStackEntry cse("Schur"))
#ifdef ELEM_HAVE_SCALAPACK
    schur::QR( A, w, Q );
#else
    schur::SDC( A, w, Q );
#endif
}

} // namespace elem

#endif // ifndef ELEM_SCHUR_HPP
