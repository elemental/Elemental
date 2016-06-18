/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_BLAS_ROTATE_HPP
#define EL_BLAS_ROTATE_HPP

namespace El {

template<typename F>
void Rotate( Base<F> c, F s, Matrix<F>& a1, Matrix<F>& a2 )
{
    DEBUG_CSE
    const Int m1 = a1.Height();
    const Int n1 = a1.Width();
    const Int n2 = a2.Width();

    const Int n = ( n1==1 ? m1 : n1 );
    const Int inc1 = ( n1==1 ? 1 : a1.LDim() );
    const Int inc2 = ( n2==1 ? 1 : a2.LDim() );

    blas::Rot( n, a1.Buffer(), inc1, a2.Buffer(), inc2, c, s );
}

template<typename F>
void Rotate
( Base<F> c, F s, AbstractDistMatrix<F>& a1, AbstractDistMatrix<F>& a2 )
{
    DEBUG_CSE
    DistMatrix<F,STAR,STAR> G(2,2,a1.Grid());
    G.Set(0,0,c);
    G.Set(0,1,s);
    G.Set(1,0,-Conj(s));
    G.Set(1,1,c);
    Transform2x2( G, a1, a2 );
}

template<typename F>
void RotateRows( Base<F> c, F s, Matrix<F>& A, Int i1, Int i2 )
{
    DEBUG_CSE
    Matrix<F> G(2,2);
    G.Set(0,0,c);
    G.Set(0,1,s);
    G.Set(1,0,-Conj(s));
    G.Set(1,1,c);
    Transform2x2Rows( G, A, i1, i2 );
}

template<typename F>
void RotateRows( Base<F> c, F s, AbstractDistMatrix<F>& A, Int i1, Int i2 )
{
    DEBUG_CSE
    DistMatrix<F,STAR,STAR> G(2,2,A.Grid());
    G.Set(0,0,c);
    G.Set(0,1,s);
    G.Set(1,0,-Conj(s));
    G.Set(1,1,c);
    Transform2x2Rows( G, A, i1, i2 );
}

template<typename F>
void RotateCols( Base<F> c, F s, Matrix<F>& A, Int i1, Int i2 )
{
    DEBUG_CSE
    Matrix<F> G(2,2);
    G.Set(0,0,c);
    G.Set(0,1,s);
    G.Set(1,0,-Conj(s));
    G.Set(1,1,c);
    Transform2x2Cols( G, A, i1, i2 );
}

template<typename F>
void RotateCols( Base<F> c, F s, AbstractDistMatrix<F>& A, Int i1, Int i2 )
{
    DEBUG_CSE
    DistMatrix<F,STAR,STAR> G(2,2,A.Grid());
    G.Set(0,0,c);
    G.Set(0,1,s);
    G.Set(1,0,-Conj(s));
    G.Set(1,1,c);
    Transform2x2Cols( G, A, i1, i2 );
}

#ifdef EL_INSTANTIATE_BLAS_LEVEL1
# define EL_EXTERN
#else
# define EL_EXTERN extern
#endif

#define PROTO(F) \
  EL_EXTERN template void Rotate \
  ( Base<F> c, F s, Matrix<F>& a1, Matrix<F>& a2 ); \
  EL_EXTERN template void Rotate \
  ( Base<F> c, F s, AbstractDistMatrix<F>& a1, AbstractDistMatrix<F>& a2 ); \
  EL_EXTERN template void RotateRows \
  ( Base<F> c, F s, Matrix<F>& A, Int i1, Int i2 ); \
  EL_EXTERN template void RotateRows \
  ( Base<F> c, F s, AbstractDistMatrix<F>& A, Int i1, Int i2 ); \
  EL_EXTERN template void RotateCols \
  ( Base<F> c, F s, Matrix<F>& A, Int j1, Int j2 ); \
  EL_EXTERN template void RotateCols \
  ( Base<F> c, F s, AbstractDistMatrix<F>& A, Int j1, Int j2 );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

#undef EL_EXTERN

} // namespace El

#endif // ifndef EL_BLAS_ROTATE_HPP
