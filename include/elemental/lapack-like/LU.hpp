/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_LU_HPP
#define ELEM_LAPACK_LU_HPP

#include "elemental/blas-like/level2/ApplyRowPivots.hpp"
#include "elemental/blas-like/level3/Gemm.hpp"
#include "elemental/blas-like/level3/Trsm.hpp"

#include "elemental/lapack-like/LU/Local.hpp"
#include "elemental/lapack-like/LU/Panel.hpp"
#include "elemental/lapack-like/LU/Full.hpp"

#include "elemental/lapack-like/LU/SolveAfter.hpp"

namespace elem {

template<typename F>
inline void
LocalLU( DistMatrix<F,STAR,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("LocalLU"))
    LU( A.Matrix() );
}

// Performs LU factorization without pivoting

template<typename F> 
inline void
LU( Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("LU"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);
    const Int bsize = Blocksize();
    for( Int k=0; k<minDim; k+=bsize )
    {
        const Int nb = Min(bsize,minDim-k);
        auto A11 = ViewRange( A, k,    k,    k+nb, k+nb );   
        auto A12 = ViewRange( A, k,    k+nb, k+nb, n    );
        auto A21 = ViewRange( A, k+nb, k,    m,    k+nb );
        auto A22 = ViewRange( A, k+nb, k+nb, m,    n    );

        lu::Unb( A11 );
        Trsm( RIGHT, UPPER, NORMAL, NON_UNIT, F(1), A11, A21 );
        Trsm( LEFT, LOWER, NORMAL, UNIT, F(1), A11, A12 );
        Gemm( NORMAL, NORMAL, F(-1), A21, A12, F(1), A22 );
    }
}

template<typename F> 
inline void
LU( DistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("LU"))
    const Grid& g = A.Grid();
    DistMatrix<F,STAR,STAR> A11_STAR_STAR(g);
    DistMatrix<F,MC,  STAR> A21_MC_STAR(g);
    DistMatrix<F,STAR,VR  > A12_STAR_VR(g);
    DistMatrix<F,STAR,MR  > A12_STAR_MR(g);

    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);
    const Int bsize = Blocksize();
    for( Int k=0; k<minDim; k+=bsize )
    {
        const Int nb = Min(bsize,minDim-k);
        auto A11 = ViewRange( A, k,    k,    k+nb, k+nb );          
        auto A12 = ViewRange( A, k,    k+nb, k+nb, n    );
        auto A21 = ViewRange( A, k+nb, k,    m,    k+nb );
        auto A22 = ViewRange( A, k+nb, k+nb, m,    n    );

        A11_STAR_STAR = A11;
        LocalLU( A11_STAR_STAR );
        A11 = A11_STAR_STAR;

        A21_MC_STAR.AlignWith( A22 );
        A21_MC_STAR = A21;
        LocalTrsm
        ( RIGHT, UPPER, NORMAL, NON_UNIT, F(1), A11_STAR_STAR, A21_MC_STAR );
        A21 = A21_MC_STAR;

        // Perhaps we should give up perfectly distributing this operation since
        // it's total contribution is only O(n^2)
        A12_STAR_VR.AlignWith( A22 );
        A12_STAR_VR = A12;
        LocalTrsm
        ( LEFT, LOWER, NORMAL, UNIT, F(1), A11_STAR_STAR, A12_STAR_VR );

        A12_STAR_MR.AlignWith( A22 );
        A12_STAR_MR = A12_STAR_VR;
        LocalGemm( NORMAL, NORMAL, F(-1), A21_MC_STAR, A12_STAR_MR, F(1), A22 );
        A12 = A12_STAR_MR;
    }
}

// Performs LU factorization with partial pivoting

template<typename F> 
inline void
LU( Matrix<F>& A, Matrix<Int>& p )
{
    DEBUG_ONLY(CallStackEntry cse("LU"))
    std::vector<Int> image, preimage;

    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);
    p.Resize( minDim, 1 );
    const Int bsize = Blocksize();
    for( Int k=0; k<minDim; k+=bsize )
    {
        const Int nb = Min(bsize,minDim-k);
        auto A11  = ViewRange( A, k,    k,    k+nb, k+nb );
        auto A12  = ViewRange( A, k,    k+nb, k+nb, n    );
        auto A21  = ViewRange( A, k+nb, k,    m,    k+nb );
        auto A22  = ViewRange( A, k+nb, k+nb, m,    n    );
        auto ABL  = ViewRange( A, k,    0,    m,    k    );
        auto ABRL = ViewRange( A, k,    k,    m,    k+nb );
        auto ABRR = ViewRange( A, k,    k+nb, m,    n    );
        auto p1 = View( p, k, 0, nb, 1 );

        lu::Panel( ABRL, p1, k );

        ComposePivots( p1, k, image, preimage );
        ApplyRowPivots( ABL, image, preimage );
        ApplyRowPivots( ABRR, image, preimage );

        Trsm( LEFT, LOWER, NORMAL, UNIT, F(1), A11, A12 );
        Gemm( NORMAL, NORMAL, F(-1), A21, A12, F(1), A22 );
    }
}

template<typename F> 
inline void
LU( Matrix<F>& A, Matrix<Int>& p, Matrix<Int>& q )
{
    DEBUG_ONLY(CallStackEntry cse("LU"))
    p.Resize( Min(A.Height(),A.Width()), 1 );
    q.Resize( Min(A.Height(),A.Width()), 1 );
    lu::Full( A, p, q );
}

template<typename F> 
inline void
LU( DistMatrix<F>& A, DistMatrix<Int,VC,STAR>& p )
{
    DEBUG_ONLY(
        CallStackEntry cse("LU");
        if( A.Grid() != p.Grid() )
            LogicError("{A,p} must be distributed over the same grid");
    )
    std::vector<Int> image, preimage;

    const Grid& g = A.Grid();
    DistMatrix<F,  STAR,STAR> A11_STAR_STAR(g);
    DistMatrix<F,  MC,  STAR> A21_MC_STAR(g);
    DistMatrix<F,  STAR,VR  > A12_STAR_VR(g);
    DistMatrix<F,  STAR,MR  > A12_STAR_MR(g);
    DistMatrix<Int,STAR,STAR> p1_STAR_STAR(g);

    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);
    p.Resize( minDim, 1 );
    const Int bsize = Blocksize();
    for( Int k=0; k<minDim; k+=bsize )
    {
        const Int nb = Min(bsize,minDim-k);
        auto A11  = ViewRange( A, k,    k,    k+nb, k+nb );
        auto A12  = ViewRange( A, k,    k+nb, k+nb, n    );
        auto A21  = ViewRange( A, k+nb, k,    m,    k+nb );
        auto A22  = ViewRange( A, k+nb, k+nb, m,    n    );
        auto AB   = ViewRange( A, k,    0,    m,    n    );
        auto ABL  = ViewRange( A, k,    0,    m,    k    );
        auto ABRL = ViewRange( A, k,    k,    m,    k+nb );
        auto ABRR = ViewRange( A, k,    k+nb, m,    n    );
        auto p1 = View( p, k, 0, nb, 1 );

        A21_MC_STAR.AlignWith( A22 );
        A21_MC_STAR = A21;
        A11_STAR_STAR = A11;
        p1_STAR_STAR.Resize( p1.Height(), 1 );
        lu::Panel( A11_STAR_STAR, A21_MC_STAR, p1_STAR_STAR, k );
        ComposePivots( p1_STAR_STAR, k, image, preimage );
        ApplyRowPivots( AB, image, preimage );

        // Perhaps we should give up perfectly distributing this operation since
        // it's total contribution is only O(n^2)
        A12_STAR_VR.AlignWith( A22 );
        A12_STAR_VR = A12;
        LocalTrsm
        ( LEFT, LOWER, NORMAL, UNIT, F(1), A11_STAR_STAR, A12_STAR_VR );

        A12_STAR_MR.AlignWith( A22 );
        A12_STAR_MR = A12_STAR_VR;
        LocalGemm( NORMAL, NORMAL, F(-1), A21_MC_STAR, A12_STAR_MR, F(1), A22 );

        A11 = A11_STAR_STAR;
        A12 = A12_STAR_MR;
        A21 = A21_MC_STAR;
        p1 = p1_STAR_STAR;
    }
}

template<typename F> 
inline void
LU( DistMatrix<F>& A, DistMatrix<Int,VC,STAR>& p, DistMatrix<Int,VC,STAR>& q )
{
    DEBUG_ONLY(CallStackEntry cse("LU"))
    p.Resize( Min(A.Height(),A.Width()), 1 );
    q.Resize( Min(A.Height(),A.Width()), 1 );
    lu::Full( A, p, q );
}

} // namespace elem

#endif // ifndef ELEM_LAPACK_LU_HPP
