/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_SPECTRAL_PRODUCT_LANCZOS_HPP
#define EL_SPECTRAL_PRODUCT_LANCZOS_HPP

namespace El {

template<typename F,class ApplyAType,class ApplyAAdjType>
inline void ProductLanczos
(       Int m,
        Int n,
  const ApplyAType& applyA,
  const ApplyAAdjType& applyAAdj,
        Matrix<Base<F>>& T,
        Int basisSize )
{
    DEBUG_ONLY(CSE cse("ProductLanczos"))
    Matrix<F> s;
    if( m >= n )
    {
        auto applyAProd =
          [&]( const Matrix<F>& x, Matrix<F>& y )
          {
              applyA( x, s );
              applyAAdj( s, y );
          };
        Lanczos<F>( n, applyAProd, T, basisSize );
    }
    else
    {
        auto applyAProd =
          [&]( const Matrix<F>& x, Matrix<F>& y )
          {
              applyAAdj( x, s );
              applyA( s, y );
          };
        Lanczos<F>( m, applyAProd, T, basisSize );
    }
}

template<typename F,class ApplyAType,class ApplyAAdjType>
inline Base<F> ProductLanczosDecomp
(       Int m,
        Int n,
  const ApplyAType& applyA,
  const ApplyAAdjType& applyAAdj,
        Matrix<F>& V, 
        Matrix<Base<F>>& T,
        Matrix<F>& v,
        Int basisSize )
{
    DEBUG_ONLY(CSE cse("ProductLanczosDecomp"))
    Matrix<F> s;
    if( m >= n )
    {
        auto applyAProd =
          [&]( const Matrix<F>& x, Matrix<F>& y )
          {
              applyA( x, s );
              applyAAdj( s, y );
          };
        return LanczosDecomp( n, applyAProd, V, T, v, basisSize );
    }
    else
    {
        auto applyAProd =
          [&]( const Matrix<F>& x, Matrix<F>& y )
          {
              applyAAdj( x, s );
              applyA( s, y ); 
          };
        return LanczosDecomp( m, applyAProd, V, T, v, basisSize );
    }
}

template<typename F,class ApplyAType,class ApplyAAdjType>
inline void ProductLanczos
(       Int m,
        Int n,
  const ApplyAType& applyA,
  const ApplyAAdjType& applyAAdj,
        ElementalMatrix<Base<F>>& T,
        Int basisSize )
{
    DEBUG_ONLY(CSE cse("ProductLanczos"))
    mpi::Comm comm = T.Grid().Comm();
    DistMultiVec<F> s(comm);
    if( m >= n )
    {
        auto applyAProd =
          [&]( const DistMultiVec<F>& x, DistMultiVec<F>& y )
          {
              applyA( x, s );
              applyAAdj( s, y );
          };
        Lanczos<F>( n, applyAProd, T, basisSize );
    }
    else
    {
        auto applyAProd =
          [&]( const DistMultiVec<F>& x, DistMultiVec<F>& y )
          {
              applyAAdj( x, s );
              applyA( s, y );
          };
        Lanczos<F>( m, applyAProd, T, basisSize );
    }
}

template<typename F,class ApplyAType,class ApplyAAdjType>
inline Base<F> ProductLanczosDecomp
(       Int m,
        Int n,
  const ApplyAType& applyA,
  const ApplyAAdjType& applyAAdj,
        DistMultiVec<F>& V, 
        ElementalMatrix<Base<F>>& T,
        DistMultiVec<F>& v,
        Int basisSize )
{
    DEBUG_ONLY(CSE cse("ProductLanczosDecomp"))
    mpi::Comm comm = T.Grid().Comm();
    DistMultiVec<F> s(comm);
    if( m >= n )
    {
        auto applyAProd =
          [&]( const DistMultiVec<F>& x, DistMultiVec<F>& y )
          {
              applyA( x, s );
              applyAAdj( s, y );
          };
        return LanczosDecomp( n, applyAProd, V, T, v, basisSize );
    }
    else
    {
        auto applyAProd =
          [&]( const DistMultiVec<F>& x, DistMultiVec<F>& y )
          {
              applyAAdj( x, s );
              applyA( s, y );
          };
        return LanczosDecomp( m, applyAProd, V, T, v, basisSize );
    }
}

} // namespace El

#endif // ifndef EL_SPECTRAL_PRODUCT_LANCZOS
