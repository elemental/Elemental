/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_SCHUR_HPP
#define ELEM_LAPACK_SCHUR_HPP

#include "elemental/lapack-like/Schur/QR.hpp"
#include "elemental/lapack-like/Schur/SDC.hpp"
#include "elemental/lapack-like/Schur/InverseFreeSDC.hpp"

namespace elem {

template<typename F>
inline void
Schur( Matrix<F>& A, Matrix<Complex<BASE(F)>>& w )
{
    DEBUG_ONLY(CallStackEntry cse("Schur"))
    schur::QR( A, w );
}

template<typename F>
inline void
Schur( Matrix<F>& A, Matrix<Complex<BASE(F)>>& w, Matrix<F>& Q )
{
    DEBUG_ONLY(CallStackEntry cse("Schur"))
    schur::QR( A, w, Q );
}

template<typename F>
inline void
Schur( DistMatrix<F>& A, DistMatrix<Complex<BASE(F)>,VR,STAR>& w )
{
    DEBUG_ONLY(CallStackEntry cse("Schur"))
    schur::SDC( A, w );
}

template<typename F>
inline void
Schur
( DistMatrix<F>& A, DistMatrix<Complex<BASE(F)>,VR,STAR>& w, DistMatrix<F>& Q )
{
    DEBUG_ONLY(CallStackEntry cse("Schur"))
    schur::SDC( A, w, Q );
}

} // namespace elem

#endif // ifndef ELEM_LAPACK_SCHUR_HPP
