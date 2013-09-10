/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_LDL_VAR3_HPP
#define ELEM_LAPACK_LDL_VAR3_HPP

#include "elemental/blas-like/level1/DiagonalSolve.hpp"
#include "elemental/blas-like/level3/Trsm.hpp"

namespace elem {
namespace ldl {

// Unblocked serial LDL _without_ partial pivoting
template<typename F> 
inline void
Var3Unb( Orientation orientation, Matrix<F>& A, Matrix<F>& d )
{
#ifndef RELEASE
    CallStackEntry entry("ldl::Var3Unb");
    if( A.Height() != A.Width() )
        LogicError("A must be square");
    if( orientation == NORMAL )
        LogicError("Can only perform LDL^T or LDL^H");
#endif
    const Int n = A.Height();
    d.ResizeTo( n, 1 );

    F* ABuffer = A.Buffer();
    F* dBuffer = d.Buffer();
    const Int ldim = A.LDim();
    for( Int j=0; j<n; ++j )
    {
        const Int a21Height = n - (j+1);

        // Extract and store the diagonal of D
        const F alpha11 = ABuffer[j+j*ldim];
        if( alpha11 == F(0) )
            throw SingularMatrixException();
        dBuffer[j] = alpha11; 

        F* RESTRICT a21 = &ABuffer[(j+1)+j*ldim];
        if( orientation == ADJOINT )
        {
            // A22 := A22 - a21 (a21 / alpha11)^H
            for( Int k=0; k<a21Height; ++k )
            {
                const F beta = Conj(a21[k]/alpha11);
                F* RESTRICT A22Col = &ABuffer[(j+1)+(j+1+k)*ldim];
                for( Int i=k; i<a21Height; ++i )
                    A22Col[i] -= a21[i]*beta;
            }
        }
        else
        {
            // A22 := A22 - a21 (a21 / alpha11)^T
            for( Int k=0; k<a21Height; ++k )
            {
                const F beta = a21[k]/alpha11;
                F* RESTRICT A22Col = &ABuffer[(j+1)+(j+1+k)*ldim];
                for( Int i=k; i<a21Height; ++i )
                    A22Col[i] -= a21[i]*beta;
            }
        }
        
        // a21 := a21 / alpha11
        for( Int i=0; i<a21Height; ++i )
            a21[i] /= alpha11;
    }
}

// Blocked serial LDL _without_ partial pivoting
template<typename F>
inline void
Var3( Orientation orientation, Matrix<F>& A, Matrix<F>& d )
{
#ifndef RELEASE
    CallStackEntry entry("ldl::Var3");
    if( A.Height() != A.Width() )
        LogicError("A must be square");
    if( orientation == NORMAL )
        LogicError("Can only perform LDL^T or LDL^H");
#endif
    const Int n = A.Height();
    d.ResizeTo( n, 1 );

    Matrix<F> S21;
    const Int bsize = Blocksize();
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = std::min(bsize,n-k);
        auto A11 = View( A, k,    k,    nb,       nb       );
        auto A21 = View( A, k+nb, k,    n-(k+nb), nb       );
        auto A22 = View( A, k+nb, k+nb, n-(k+nb), n-(k+nb) );
        auto d1  = View( d, k,    0,    nb,       1        );

        ldl::Var3Unb( orientation, A11, d1 );
        Trsm( RIGHT, LOWER, orientation, UNIT, F(1), A11, A21 );
        S21 = A21;
        DiagonalSolve( RIGHT, NORMAL, d1, A21 );
        internal::TrrkNT( LOWER, orientation, F(-1), S21, A21, F(1), A22 );
    }
}

template<typename F>
inline void
Var3( Orientation orientation, DistMatrix<F>& A, DistMatrix<F,MC,STAR>& d )
{
#ifndef RELEASE
    CallStackEntry entry("ldl::Var3");
    if( orientation == NORMAL )
        LogicError("Can only perform LDL^T and LDL^H");
    if( A.Height() != A.Width() )
        LogicError("A must be square");
    if( A.Grid() != d.Grid() )
        LogicError("A and d must use the same grid");
    if( d.Viewing() && d.ColAlignment() != A.ColAlignment() )
        LogicError("d must be aligned with A");
#endif
    const Grid& g = A.Grid();
    if( !d.Viewing() )
        d.AlignWith( A );
    const Int n = A.Height();
    d.ResizeTo( n, 1 );

    DistMatrix<F,STAR,STAR> A11_STAR_STAR(g);
    DistMatrix<F,STAR,STAR> d1_STAR_STAR(g);
    DistMatrix<F,VC,  STAR> A21_VC_STAR(g);
    DistMatrix<F,VR,  STAR> A21_VR_STAR(g);
    DistMatrix<F,STAR,MC  > S21Trans_STAR_MC(g);
    DistMatrix<F,STAR,MR  > A21AdjOrTrans_STAR_MR(g);

    const bool conjugate = ( orientation == ADJOINT );
    const Int bsize = Blocksize();
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = std::min(bsize,n-k);
        auto A11 = View( A, k,    k,    nb,       nb       );
        auto A21 = View( A, k+nb, k,    n-(k+nb), nb       );
        auto A22 = View( A, k+nb, k+nb, n-(k+nb), n-(k+nb) );
        auto d1  = View( d, k,    0,    nb,       1        );

        A11_STAR_STAR = A11;
        LocalLDL( orientation, A11_STAR_STAR, d1_STAR_STAR );
        A11 = A11_STAR_STAR;
        d1 = d1_STAR_STAR;

        A21_VC_STAR.AlignWith( A22 );
        A21_VC_STAR = A21;
        LocalTrsm
        ( RIGHT, LOWER, orientation, UNIT,
          F(1), A11_STAR_STAR, A21_VC_STAR );

        S21Trans_STAR_MC.AlignWith( A22 );
        S21Trans_STAR_MC.TransposeFrom( A21_VC_STAR );
        DiagonalSolve( RIGHT, NORMAL, d1_STAR_STAR, A21_VC_STAR );
        A21_VR_STAR.AlignWith( A22 );
        A21_VR_STAR = A21_VC_STAR;
        A21AdjOrTrans_STAR_MR.AlignWith( A22 );
        A21AdjOrTrans_STAR_MR.TransposeFrom( A21_VR_STAR, conjugate );
        LocalTrrk
        ( LOWER, TRANSPOSE,
          F(-1), S21Trans_STAR_MC, A21AdjOrTrans_STAR_MR, F(1), A22 );

        A21 = A21_VC_STAR;
    }
}

} // namespace ldl
} // namespace elem

#endif // ifndef ELEM_LAPACK_LDL_VAR3_HPP
