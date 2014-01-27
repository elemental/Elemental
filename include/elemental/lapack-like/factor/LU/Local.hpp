/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LU_LOCAL_HPP
#define ELEM_LU_LOCAL_HPP

#include ELEM_SCALE_INC
#include ELEM_GERU_INC

namespace elem {
namespace lu {

// Local LU _without_ partial pivoting

template<typename F> 
inline void
UnbFLAME( Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("lu::UnbFLAME"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);
    for( Int k=0; k<minDim; ++k )
    {
        auto alpha11 = ViewRange( A, k,   k,   k+1, k+1 );
        auto a12     = ViewRange( A, k,   k+1, k+1, n   );
        auto a21     = ViewRange( A, k+1, k,   m,   k+1 );
        auto A22     = ViewRange( A, k+1, k+1, m,   n   );

        F alpha = alpha11.Get(0,0);
        if( alpha == F(0) )
            throw SingularMatrixException();
        Scale( 1/alpha, a21 );
        Geru( F(-1), a21, a12, A22 );
    }
}

template<typename F>
inline void
Unb( Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("lu::Unb"))
    const Int m = A.Height();
    const Int n = A.Width();
    for( Int j=0; j<Min(m,n); ++j )
    {
        const F alpha = A.Get(j,j);
        if( alpha == F(0) )
            throw SingularMatrixException();

        blas::Scal( m-(j+1), F(1)/alpha, A.Buffer(j+1,j), 1 );
        blas::Geru
        ( m-(j+1), n-(j+1),
          F(-1), A.LockedBuffer(j+1,j), 1, A.LockedBuffer(j,j+1), A.LDim(),
                 A.Buffer(j+1,j+1), A.LDim() );
    }
}

} // namespace lu
} // namespace elem

#endif // ifndef ELEM_LU_LOCAL_HPP
