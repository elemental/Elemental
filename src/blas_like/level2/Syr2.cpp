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
void Syr2
( UpperOrLower uplo,
  T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A, 
  bool conjugate )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( A.Height() != A.Width() )
          LogicError("A must be square");
      if( (x.Width() != 1 && x.Height() != 1) ||
          (y.Width() != 1 && y.Height() != 1) )
          LogicError("x and y must be vectors");
      const Int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
      const Int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
      if( xLength != A.Height() || yLength != A.Height() )
          LogicError("x and y must conform with A");
    )
    const char uploChar = UpperOrLowerToChar( uplo );
    const Int m = A.Height();
    const Int incx = ( x.Width()==1 ? 1 : x.LDim() );
    const Int incy = ( y.Width()==1 ? 1 : y.LDim() );
    if( conjugate )
    {
        blas::Her2
        ( uploChar, m,
          alpha, x.LockedBuffer(), incx, y.LockedBuffer(), incy,
                 A.Buffer(), A.LDim() );
    }
    else
    {
        blas::Syr2
        ( uploChar, m,
          alpha, x.LockedBuffer(), incx, y.LockedBuffer(), incy,
                 A.Buffer(), A.LDim() );
    }
}

template<typename T>
void Syr2
( UpperOrLower uplo,
  T alpha, const AbstractDistMatrix<T>& x,
           const AbstractDistMatrix<T>& y,
                 AbstractDistMatrix<T>& APre, bool conjugate )
{
    DEBUG_CSE
    DEBUG_ONLY(
      AssertSameGrids( APre, x, y );
      if( APre.Height() != APre.Width() )
          LogicError("A must be square");
      const Int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
      const Int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
      if( APre.Height() != xLength || APre.Height() != yLength )
          LogicError
          ("A must conform with x: \n",DimsString(APre,"A"),"\n",
           DimsString(x,"x"),"\n",DimsString(y,"y"));
    )

    DistMatrixReadWriteProxy<T,T,MC,MR> AProx( APre );
    auto& A = AProx.Get();

    const Grid& g = A.Grid();
    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();

    if( x.Width() == 1 && y.Width() == 1 )
    {
        DistMatrix<T,MC,STAR> x_MC_STAR(g), y_MC_STAR(g);
        x_MC_STAR.AlignWith( A );
        y_MC_STAR.AlignWith( A );
        x_MC_STAR = x;
        y_MC_STAR = y;

        DistMatrix<T,MR,STAR> x_MR_STAR(g), y_MR_STAR(g);
        x_MR_STAR.AlignWith( A );
        y_MR_STAR.AlignWith( A );
        x_MR_STAR = x_MC_STAR;
        y_MR_STAR = y_MC_STAR;

        const T* xBuffer = x_MC_STAR.LockedBuffer();
        const T* yBuffer = y_MC_STAR.LockedBuffer();
        if( uplo == LOWER )
        {
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const Int j = A.GlobalCol(jLoc);
                const Int heightAboveDiag = A.LocalRowOffset(j);

                const T beta = y_MR_STAR.GetLocal(jLoc,0);
                const T kappa = x_MR_STAR.GetLocal(jLoc,0);
                const T gamma = ( conjugate ? alpha*Conj(beta)  : alpha*beta );
                const T delta = ( conjugate ? Conj(alpha*kappa) : alpha*kappa );
                T* ACol = A.Buffer(0,jLoc);
                for( Int iLoc=heightAboveDiag; iLoc<localHeight; ++iLoc )
                    ACol[iLoc] += gamma*xBuffer[iLoc] + delta*yBuffer[iLoc];
            }
        }
        else
        {
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const Int j = A.GlobalCol(jLoc);
                const Int heightToDiag = A.LocalRowOffset(j+1);

                const T beta = y_MR_STAR.GetLocal(jLoc,0);
                const T kappa = x_MR_STAR.GetLocal(jLoc,0);
                const T gamma = ( conjugate ? alpha*Conj(beta)  : alpha*beta );
                const T delta = ( conjugate ? Conj(alpha*kappa) : alpha*kappa );
                T* ACol = A.Buffer(0,jLoc);
                for( Int iLoc=0; iLoc<heightToDiag; ++iLoc )
                    ACol[iLoc] += gamma*xBuffer[iLoc] + delta*yBuffer[iLoc];
            }
        }
    }
    else if( x.Width() == 1 )
    {
        DistMatrix<T,MC,STAR> x_MC_STAR(g);
        x_MC_STAR.AlignWith( A );
        x_MC_STAR = x;

        DistMatrix<T,MR,STAR> x_MR_STAR(g);
        x_MR_STAR.AlignWith( A );
        x_MR_STAR = x_MC_STAR;

        DistMatrix<T,STAR,MR> y_STAR_MR(g);
        y_STAR_MR.AlignWith( A );
        y_STAR_MR = y;

        DistMatrix<T,STAR,MC> y_STAR_MC(g);
        y_STAR_MC.AlignWith( A );
        y_STAR_MC = y_STAR_MR;

        const T* xBuffer = x_MC_STAR.LockedBuffer();
        const T* yBuffer = y_STAR_MC.LockedBuffer();
        const Int incy = y_STAR_MC.LDim();
        if( uplo == LOWER )
        {
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const Int j = A.GlobalCol(jLoc);
                const Int heightAboveDiag = A.LocalRowOffset(j);

                const T beta = y_STAR_MR.GetLocal(0,jLoc);
                const T kappa = x_MR_STAR.GetLocal(jLoc,0);
                const T gamma = ( conjugate ? alpha*Conj(beta) :  alpha*beta );
                const T delta = ( conjugate ? Conj(alpha*kappa) : alpha*kappa );
                T* ACol = A.Buffer(0,jLoc);
                for( Int iLoc=heightAboveDiag; iLoc<localHeight; ++iLoc )
                    ACol[iLoc] += gamma*xBuffer[iLoc] +
                                  delta*yBuffer[iLoc*incy];
            }
        }
        else
        {
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const Int j = A.GlobalCol(jLoc);
                const Int heightToDiag = A.LocalRowOffset(j+1);

                const T beta = y_STAR_MR.GetLocal(0,jLoc);
                const T kappa = x_MR_STAR.GetLocal(jLoc,0);
                const T gamma = ( conjugate ? alpha*Conj(beta)  : alpha*beta );
                const T delta = ( conjugate ? Conj(alpha*kappa) : alpha*kappa );
                T* ACol = A.Buffer(0,jLoc);
                for( Int iLoc=0; iLoc<heightToDiag; ++iLoc )
                    ACol[iLoc] += gamma*xBuffer[iLoc] +
                                  delta*yBuffer[iLoc*incy];
            }
        }
    }
    else if( y.Width() == 1 )
    {
        DistMatrix<T,STAR,MR> x_STAR_MR(g);
        x_STAR_MR.AlignWith( A );
        x_STAR_MR = x;

        DistMatrix<T,STAR,MC> x_STAR_MC(g);
        x_STAR_MC.AlignWith( A );
        x_STAR_MC = x_STAR_MR;

        DistMatrix<T,MC,STAR> y_MC_STAR(g);
        y_MC_STAR.AlignWith( A );
        y_MC_STAR = y;

        DistMatrix<T,MR,STAR> y_MR_STAR(g);
        y_MR_STAR.AlignWith( A );
        y_MR_STAR = y_MC_STAR;

        const T* xBuffer = x_STAR_MC.LockedBuffer();
        const T* yBuffer = y_MC_STAR.LockedBuffer();
        const Int incx = x_STAR_MC.LDim();
        if( uplo == LOWER )
        {
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const Int j = A.GlobalCol(jLoc);
                const Int heightAboveDiag = A.LocalRowOffset(j);

                const T beta = x_STAR_MR.GetLocal(0,jLoc);
                const T kappa = y_MR_STAR.GetLocal(jLoc,0);
                const T gamma = ( conjugate ? alpha*Conj(beta)  : alpha*beta );
                const T delta = ( conjugate ? Conj(alpha*kappa) : alpha*kappa );
                T* ACol = A.Buffer(0,jLoc);
                for( Int iLoc=heightAboveDiag; iLoc<localHeight; ++iLoc )
                    ACol[iLoc] += gamma*xBuffer[iLoc*incx] +
                                  delta*yBuffer[iLoc]; 
            }
        }
        else
        {
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const Int j = A.GlobalCol(jLoc);
                const Int heightToDiag = A.LocalRowOffset(j+1);

                const T beta = x_STAR_MR.GetLocal(0,jLoc);
                const T kappa = y_MR_STAR.GetLocal(jLoc,0);
                const T gamma = ( conjugate ? alpha*Conj(beta)  : alpha*beta );
                const T delta = ( conjugate ? Conj(alpha*kappa) : alpha*kappa );
                T* ACol = A.Buffer(0,jLoc);
                for( Int iLoc=0; iLoc<heightToDiag; ++iLoc )
                    ACol[iLoc] += gamma*xBuffer[iLoc*incx] +
                                  delta*yBuffer[iLoc];
            }
        }
    }
    else
    {
        DistMatrix<T,STAR,MR> x_STAR_MR(g), y_STAR_MR(g);
        x_STAR_MR.AlignWith( A );
        y_STAR_MR.AlignWith( A );
        y_STAR_MR = y;
        x_STAR_MR = x;

        DistMatrix<T,STAR,MC> x_STAR_MC(g), y_STAR_MC(g);
        x_STAR_MC.AlignWith( A );
        y_STAR_MC.AlignWith( A );
        x_STAR_MC = x_STAR_MR;
        y_STAR_MC = y_STAR_MR;

        const T* xBuffer = x_STAR_MC.LockedBuffer();
        const T* yBuffer = y_STAR_MC.LockedBuffer();
        const Int incx = x_STAR_MC.LDim();
        const Int incy = y_STAR_MC.LDim();
        if( uplo == LOWER )
        {
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const Int j = A.GlobalCol(jLoc);
                const Int heightAboveDiag = A.LocalRowOffset(j);

                const T beta = y_STAR_MR.GetLocal(0,jLoc);
                const T kappa = x_STAR_MR.GetLocal(0,jLoc);
                const T gamma = ( conjugate ? alpha*Conj(beta)  : alpha*beta );
                const T delta = ( conjugate ? Conj(alpha*kappa) : alpha*kappa );
                T* ACol = A.Buffer(0,jLoc);
                for( Int iLoc=heightAboveDiag; iLoc<localHeight; ++iLoc )
                    ACol[iLoc] += gamma*xBuffer[iLoc*incx] +
                                  delta*yBuffer[iLoc*incy];
            }
        }
        else
        {
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const Int j = A.GlobalCol(jLoc);
                const Int heightToDiag = A.LocalRowOffset(j+1);

                const T beta = y_STAR_MR.GetLocal(0,jLoc);
                const T kappa = x_STAR_MR.GetLocal(0,jLoc);
                const T gamma = ( conjugate ? alpha*Conj(beta)  : alpha*beta );
                const T delta = ( conjugate ? Conj(alpha*kappa) : alpha*kappa );
                T* ACol = A.Buffer(0,jLoc);
                for( Int iLoc=0; iLoc<heightToDiag; ++iLoc )
                    ACol[iLoc] += gamma*xBuffer[iLoc*incx] +
                                  delta*yBuffer[iLoc*incy];
            }
        }
    }
}

#define PROTO(T) \
  template void Syr2 \
  ( UpperOrLower uplo, T alpha, \
    const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A, bool conjugate ); \
  template void Syr2 \
  ( UpperOrLower uplo, T alpha, \
    const AbstractDistMatrix<T>& x, const AbstractDistMatrix<T>& y, \
    AbstractDistMatrix<T>& A, bool conjugate );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
