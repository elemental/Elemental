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

template<typename Field,class ApplyAType,class ApplyAAdjType>
void ProductLanczos
(       Int m,
        Int n,
  const ApplyAType& applyA,
  const ApplyAAdjType& applyAAdj,
        Matrix<Base<Field>>& T,
        Int basisSize )
{
    EL_DEBUG_CSE
    Matrix<Field> s;
    if( m >= n )
    {
        auto applyAProd =
          [&]( const Matrix<Field>& x, Matrix<Field>& y )
          {
              applyA( x, s );
              applyAAdj( s, y );
          };
        Lanczos<Field>( n, applyAProd, T, basisSize );
    }
    else
    {
        auto applyAProd =
          [&]( const Matrix<Field>& x, Matrix<Field>& y )
          {
              applyAAdj( x, s );
              applyA( s, y );
          };
        Lanczos<Field>( m, applyAProd, T, basisSize );
    }
}

template<typename Field,class ApplyAType,class ApplyAAdjType>
Base<Field> ProductLanczosDecomp
(       Int m,
        Int n,
  const ApplyAType& applyA,
  const ApplyAAdjType& applyAAdj,
        Matrix<Field>& V,
        Matrix<Base<Field>>& T,
        Matrix<Field>& v,
        Int basisSize )
{
    EL_DEBUG_CSE
    Matrix<Field> s;
    if( m >= n )
    {
        auto applyAProd =
          [&]( const Matrix<Field>& x, Matrix<Field>& y )
          {
              applyA( x, s );
              applyAAdj( s, y );
          };
        return LanczosDecomp( n, applyAProd, V, T, v, basisSize );
    }
    else
    {
        auto applyAProd =
          [&]( const Matrix<Field>& x, Matrix<Field>& y )
          {
              applyAAdj( x, s );
              applyA( s, y );
          };
        return LanczosDecomp( m, applyAProd, V, T, v, basisSize );
    }
}

template<typename Field,class ApplyAType,class ApplyAAdjType>
void ProductLanczos
(       Int m,
        Int n,
  const ApplyAType& applyA,
  const ApplyAAdjType& applyAAdj,
        AbstractDistMatrix<Base<Field>>& T,
        Int basisSize )
{
    EL_DEBUG_CSE
    mpi::Comm comm = T.Grid().Comm();
    DistMultiVec<Field> s(comm);
    if( m >= n )
    {
        auto applyAProd =
          [&]( const DistMultiVec<Field>& x, DistMultiVec<Field>& y )
          {
              applyA( x, s );
              applyAAdj( s, y );
          };
        Lanczos<Field>( n, applyAProd, T, basisSize );
    }
    else
    {
        auto applyAProd =
          [&]( const DistMultiVec<Field>& x, DistMultiVec<Field>& y )
          {
              applyAAdj( x, s );
              applyA( s, y );
          };
        Lanczos<Field>( m, applyAProd, T, basisSize );
    }
}

template<typename Field,class ApplyAType,class ApplyAAdjType>
Base<Field> ProductLanczosDecomp
(       Int m,
        Int n,
  const ApplyAType& applyA,
  const ApplyAAdjType& applyAAdj,
        DistMultiVec<Field>& V,
        AbstractDistMatrix<Base<Field>>& T,
        DistMultiVec<Field>& v,
        Int basisSize )
{
    EL_DEBUG_CSE
    mpi::Comm comm = T.Grid().Comm();
    DistMultiVec<Field> s(comm);
    if( m >= n )
    {
        auto applyAProd =
          [&]( const DistMultiVec<Field>& x, DistMultiVec<Field>& y )
          {
              applyA( x, s );
              applyAAdj( s, y );
          };
        return LanczosDecomp( n, applyAProd, V, T, v, basisSize );
    }
    else
    {
        auto applyAProd =
          [&]( const DistMultiVec<Field>& x, DistMultiVec<Field>& y )
          {
              applyAAdj( x, s );
              applyA( s, y );
          };
        return LanczosDecomp( m, applyAProd, V, T, v, basisSize );
    }
}

} // namespace El

#endif // ifndef EL_SPECTRAL_PRODUCT_LANCZOS
