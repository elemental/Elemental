/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El-lite.hpp>
#include <El/blas_like/level2.hpp>

namespace El {

template<typename T>
void Syr
( UpperOrLower uplo,
  T alpha,
  const Matrix<T>& x,
        Matrix<T>& A, 
  bool conjugate )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( A.Height() != A.Width() )
          LogicError("A must be square");
      if( x.Width() != 1 && x.Height() != 1 )
          LogicError("x must be a vector");
      const Int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
      if( xLength != A.Height() )
          LogicError("x must conform with A");
    )
    const char uploChar = UpperOrLowerToChar( uplo );
    const Int m = A.Height();
    const Int incx = ( x.Width()==1 ? 1 : x.LDim() );
    if( conjugate )
    {
        blas::Her
        ( uploChar, m, 
          RealPart(alpha), x.LockedBuffer(), incx, A.Buffer(), A.LDim() );
    }
    else
    {
        blas::Syr
        ( uploChar, m, alpha, x.LockedBuffer(), incx, A.Buffer(), A.LDim() );
    }
}

template<typename T>
void Syr
( UpperOrLower uplo,
  T alpha,
  const AbstractDistMatrix<T>& x,
        AbstractDistMatrix<T>& APre,
  bool conjugate )
{
    DEBUG_CSE
    DEBUG_ONLY(
      AssertSameGrids( APre, x );
      if( APre.Height() != APre.Width() )
          LogicError("A must be square");
      const Int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
      if( APre.Height() != xLength )
          LogicError
          ("A must conform with x: \n",DimsString(APre,"A"),"\n",
           DimsString(x,"x"));
    )

    DistMatrixReadWriteProxy<T,T,MC,MR> AProx( APre );
    auto& A = AProx.Get();

    const Grid& g = A.Grid();
    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();

    if( x.Width() == 1 )
    {
        DistMatrix<T,MC,STAR> x_MC_STAR(g);
        x_MC_STAR.AlignWith( A );
        x_MC_STAR = x;

        DistMatrix<T,MR,STAR> x_MR_STAR(g);
        x_MR_STAR.AlignWith( A );
        x_MR_STAR = x_MC_STAR;

        const T* xBuffer = x_MC_STAR.LockedBuffer();
        if( uplo == LOWER )
        {
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const Int j = A.GlobalCol(jLoc);
                const Int heightAboveDiag = A.LocalRowOffset(j);

                const T beta = x_MR_STAR.GetLocal(jLoc,0);
                const T gamma = ( conjugate ? alpha*Conj(beta) : alpha*beta );
                T* ACol = A.Buffer(0,jLoc);
                for( Int iLoc=heightAboveDiag; iLoc<localHeight; ++iLoc )
                    ACol[iLoc] += gamma*xBuffer[iLoc];
            }
        }
        else
        {
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const Int j = A.GlobalCol(jLoc);
                const Int heightToDiag = A.LocalRowOffset(j+1);

                const T beta = x_MR_STAR.GetLocal(jLoc,0);
                const T gamma = ( conjugate ? alpha*Conj(beta) : alpha*beta );
                T* ACol = A.Buffer(0,jLoc);
                for( Int iLoc=0; iLoc<heightToDiag; ++iLoc )
                    ACol[iLoc] += gamma*xBuffer[iLoc];
            }
        }
    }
    else
    {
        DistMatrix<T,STAR,MR> x_STAR_MR(g);
        x_STAR_MR.AlignWith( A );
        x_STAR_MR = x;

        DistMatrix<T,STAR,MC> x_STAR_MC(g);
        x_STAR_MC.AlignWith( A );
        x_STAR_MC = x_STAR_MR;

        const T* xBuffer = x_STAR_MC.LockedBuffer();
        const Int incx = x_STAR_MC.LDim();
        if( uplo == LOWER )
        {
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const Int j = A.GlobalCol(jLoc);
                const Int heightAboveDiag = A.LocalRowOffset(j);

                const T beta = x_STAR_MR.GetLocal(0,jLoc);
                const T gamma = ( conjugate ? alpha*Conj(beta) : alpha*beta );
                T* ACol = A.Buffer(0,jLoc);
                for( Int iLoc=heightAboveDiag; iLoc<localHeight; ++iLoc )
                    ACol[iLoc] += gamma*xBuffer[iLoc*incx];
            }
        }
        else
        {
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const Int j = A.GlobalCol(jLoc);
                const Int heightToDiag = A.LocalRowOffset(j+1);

                const T beta = x_STAR_MR.GetLocal(0,jLoc);
                const T gamma = ( conjugate ? alpha*Conj(beta) : alpha*beta );
                T* ACol = A.Buffer(0,jLoc);
                for( Int iLoc=0; iLoc<heightToDiag; ++iLoc )
                    ACol[iLoc] += gamma*xBuffer[iLoc*incx];
            }
        }
    }
}

#define PROTO(T) \
  template void Syr \
  ( UpperOrLower uplo, T alpha, \
    const Matrix<T>& x, Matrix<T>& A, bool conjugate ); \
  template void Syr \
  ( UpperOrLower uplo, T alpha, \
    const AbstractDistMatrix<T>& x, AbstractDistMatrix<T>& A, bool conjugate );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
