/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_SCHUR_QR_HPP
#define ELEM_LAPACK_SCHUR_QR_HPP

namespace elem {
namespace schur {

template<typename F>
inline void
QR( Matrix<F>& A, Matrix<Complex<BASE(F)>>& w, bool formATR=false )
{
    DEBUG_ONLY(CallStackEntry cse("schur::qr"))
    const Int n = A.Height();
    w.Resize( n, 1 );
    lapack::Eig( n, A.Buffer(), A.LDim(), w.Buffer(), formATR );
}

template<typename F>
inline void
QR
( Matrix<F>& A, Matrix<Complex<BASE(F)>>& w, Matrix<F>& Q, bool formATR=true )
{
    DEBUG_ONLY(CallStackEntry cse("schur::qr"))
    const Int n = A.Height();
    Q.Resize( n, n );
    w.Resize( n, 1 );
    lapack::Schur
    ( n, A.Buffer(), A.LDim(), Q.Buffer(), Q.LDim(), w.Buffer(), formATR );
}

} // namespace schur
} // namespace elem

#endif // ifndef ELEM_LAPACK_SCHUR_QR_HPP
