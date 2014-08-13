/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename F> 
void GRQ( Matrix<F>& A, Matrix<F>& B )
{
    DEBUG_ONLY(CallStackEntry cse("GRQ"))
    Matrix<F> tA;
    Matrix<Base<F>> dA;
    RQ( A, tA, dA );
    rq::ApplyQ( RIGHT, ADJOINT, A, tA, dA, B );
    MakeTrapezoidal( UPPER, A, Min(A.Height(),A.Width()) );
    QR( B );
}

template<typename F> 
void GRQ( AbstractDistMatrix<F>& APre, AbstractDistMatrix<F>& BPre )
{
    DEBUG_ONLY(CallStackEntry cse("GRQ"))

    const Grid& g = APre.Grid();
    DistMatrix<F> A(g), B(g);
    Copy( APre, A, READ_WRITE_PROXY );
    Copy( BPre, B, READ_WRITE_PROXY );

    DistMatrix<F,MD,STAR> tA(g);
    DistMatrix<Base<F>,MD,STAR> dA(g);
    RQ( A, tA, dA );
    rq::ApplyQ( RIGHT, ADJOINT, A, tA, dA, B );
    MakeTrapezoidal( UPPER, A, Min(A.Height(),A.Width()) );
    QR( B );

    Copy( A, APre, RESTORE_READ_WRITE_PROXY );
    Copy( B, BPre, RESTORE_READ_WRITE_PROXY );
}

template<typename F> 
void GRQ
( Matrix<F>& A, Matrix<F>& tA, Matrix<Base<F>>& dA, 
  Matrix<F>& B, Matrix<F>& tB, Matrix<Base<F>>& dB )
{
    DEBUG_ONLY(CallStackEntry cse("GRQ"))
    RQ( A, tA, dA );
    rq::ApplyQ( RIGHT, ADJOINT, A, tA, dA, B );
    QR( B, tB, dB );
}

template<typename F> 
void GRQ
( AbstractDistMatrix<F>& APre, 
  AbstractDistMatrix<F>& tA, AbstractDistMatrix<Base<F>>& dA,
  AbstractDistMatrix<F>& BPre, 
  AbstractDistMatrix<F>& tB, AbstractDistMatrix<Base<F>>& dB )
{
    DEBUG_ONLY(CallStackEntry cse("GRQ"))
    
    const Grid& g = APre.Grid();
    DistMatrix<F> A(g), B(g); 
    Copy( APre, A, READ_WRITE_PROXY );
    Copy( BPre, B, READ_WRITE_PROXY );

    RQ( A, tA, dA );
    rq::ApplyQ( RIGHT, ADJOINT, A, tA, dA, B );
    QR( B, tB, dB );

    Copy( A, APre, RESTORE_READ_WRITE_PROXY );
    Copy( B, BPre, RESTORE_READ_WRITE_PROXY );
}

#define PROTO(F) \
  template void GRQ( Matrix<F>& A, Matrix<F>& B ); \
  template void GRQ( AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& B ); \
  template void GRQ \
  ( Matrix<F>& A, Matrix<F>& tA, Matrix<Base<F>>& dA, \
    Matrix<F>& B, Matrix<F>& tB, Matrix<Base<F>>& dB ); \
  template void GRQ \
  ( AbstractDistMatrix<F>& A, \
    AbstractDistMatrix<F>& tA, AbstractDistMatrix<Base<F>>& dA, \
    AbstractDistMatrix<F>& B, \
    AbstractDistMatrix<F>& tB, AbstractDistMatrix<Base<F>>& dB );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
