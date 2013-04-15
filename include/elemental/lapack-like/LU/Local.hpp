/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_LU_LOCAL_HPP
#define LAPACK_LU_LOCAL_HPP

#include "elemental/blas-like/level1/Scale.hpp"
#include "elemental/blas-like/level2/Geru.hpp"

namespace elem {
namespace lu {

// Local LU _without_ partial pivoting

template<typename F> 
inline void
UnbFLAME( Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("UnbFLAME");
#endif
    // Matrix views 
    Matrix<F>
        ATL, ATR,  A00, a01,     A02,  alpha21T,
        ABL, ABR,  a10, alpha11, a12,  a21B,
                   A20, a21,     A22;

    PushBlocksizeStack( 1 );
    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    while( ATL.Height() < A.Height() && ATL.Width() < A.Width() )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ a01,     A02,
         /*************/ /**********************/
               /**/       a10, /**/ alpha11, a12,
          ABL, /**/ ABR,  A20, /**/ a21,     A22 );

        //--------------------------------------------------------------------//
        F alpha = alpha11.Get(0,0);
        if( alpha == F(0) )
            throw SingularMatrixException();
        Scale( 1/alpha, a21 );
        Geru( F(-1), a21, a12, A22 );
        //--------------------------------------------------------------------//

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, a01,     /**/ A02,
               /**/       a10, alpha11, /**/ a12,
         /*************/ /**********************/
          ABL, /**/ ABR,  A20, a21,     /**/ A22 );
    }
    PopBlocksizeStack();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
Unb( Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("Unb");
#endif
    const int m = A.Height();
    const int n = A.Width();
    for( int j=0; j<std::min(m,n); ++j )
    {
        const F alpha = A.Get(j,j);
        if( alpha == F(0) )
            throw SingularMatrixException();

        blas::Scal( m-(j+1), 1/alpha, A.Buffer(j+1,j), 1 );
        blas::Geru
        ( m-(j+1), n-(j+1),
          F(-1), A.LockedBuffer(j+1,j), 1, A.LockedBuffer(j,j+1), A.LDim(),
                 A.Buffer(j+1,j+1), A.LDim() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace lu
} // namespace elem

#endif // ifndef LAPACK_LU_LOCAL_HPP
