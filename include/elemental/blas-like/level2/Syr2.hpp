/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_BLAS_SYR2_HPP
#define ELEM_BLAS_SYR2_HPP

namespace elem {

template<typename T>
inline void
Syr2
( UpperOrLower uplo,
  T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A, 
  bool conjugate=false )
{
#ifndef RELEASE
    CallStackEntry entry("Syr2");
    if( A.Height() != A.Width() )
        LogicError("A must be square");
    if( (x.Width() != 1 && x.Height() != 1) ||
        (y.Width() != 1 && y.Height() != 1) )
        LogicError("x and y must be vectors");
    const Int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    const Int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
    if( xLength != A.Height() || yLength != A.Height() )
        LogicError("x and y must conform with A");
#endif
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
inline void
Syr2
( UpperOrLower uplo,
  T alpha, const DistMatrix<T>& x,
           const DistMatrix<T>& y,
                 DistMatrix<T>& A,
  bool conjugate=false )
{
#ifndef RELEASE
    CallStackEntry entry("Syr2");
    if( A.Grid() != x.Grid() || x.Grid() != y.Grid() )
        LogicError
        ("{A,x,y} must be distributed over the same grid");
    if( A.Height() != A.Width() )
        LogicError("A must be square");
    const Int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    const Int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
    if( A.Height() != xLength || A.Height() != yLength )
    {
        std::ostringstream msg;
        msg << "A must conform with x: \n"
            << "  A ~ " << A.Height() << " x " << A.Width() << "\n"
            << "  x ~ " << x.Height() << " x " << x.Width() << "\n"
            << "  y ~ " << y.Height() << " x " << y.Width() << "\n";
        LogicError( msg.str() );
    }
#endif
    const Grid& g = A.Grid();

    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int colShift = A.ColShift();
    const Int rowShift = A.RowShift();

    if( x.Width() == 1 && y.Width() == 1 )
    {
        DistMatrix<T,MC,STAR> x_MC_STAR(g), y_MC_STAR(g);
        DistMatrix<T,MR,STAR> x_MR_STAR(g), y_MR_STAR(g);

        x_MC_STAR.AlignWith( A );
        x_MR_STAR.AlignWith( A );
        y_MC_STAR.AlignWith( A );
        y_MR_STAR.AlignWith( A );
        //--------------------------------------------------------------------//
        x_MC_STAR = x;
        x_MR_STAR = x_MC_STAR;
        y_MC_STAR = y;
        y_MR_STAR = y_MC_STAR;

        const T* xBuffer = x_MC_STAR.LockedBuffer();
        const T* yBuffer = y_MC_STAR.LockedBuffer();
        if( uplo == LOWER )
        {
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const Int j = rowShift + jLoc*c;
                const Int heightAboveDiag = Length(j,colShift,r);

                const T beta = y_MR_STAR.GetLocal(jLoc,0);
                const T kappa = x_MR_STAR.GetLocal(jLoc,0);
                const T gamma = ( conjugate ? alpha*Conj(beta) : alpha*beta );
                const T delta = ( conjugate ? alpha*Conj(kappa) : alpha*kappa );
                T* ACol = A.Buffer(0,jLoc);
                for( Int iLoc=heightAboveDiag; iLoc<localHeight; ++iLoc )
                    ACol[iLoc] += gamma*xBuffer[iLoc] + delta*yBuffer[iLoc];
            }
        }
        else
        {
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const Int j = rowShift + jLoc*c;
                const Int heightToDiag = Length(j+1,colShift,r);

                const T beta = y_MR_STAR.GetLocal(jLoc,0);
                const T kappa = x_MR_STAR.GetLocal(jLoc,0);
                const T gamma = ( conjugate ? alpha*Conj(beta) : alpha*beta );
                const T delta = ( conjugate ? alpha*Conj(kappa) : alpha*kappa );
                T* ACol = A.Buffer(0,jLoc);
                for( Int iLoc=0; iLoc<heightToDiag; ++iLoc )
                    ACol[iLoc] += gamma*xBuffer[iLoc] + delta*yBuffer[iLoc];
            }
        }
        //--------------------------------------------------------------------//
    }
    else if( x.Width() == 1 )
    {
        DistMatrix<T,MC,STAR> x_MC_STAR(g);
        DistMatrix<T,MR,STAR> x_MR_STAR(g);
        DistMatrix<T,STAR,MC> y_STAR_MC(g);
        DistMatrix<T,STAR,MR> y_STAR_MR(g);

        x_MC_STAR.AlignWith( A );
        x_MR_STAR.AlignWith( A );
        y_STAR_MC.AlignWith( A );
        y_STAR_MR.AlignWith( A );
        //--------------------------------------------------------------------//
        x_MC_STAR = x;
        x_MR_STAR = x_MC_STAR;
        y_STAR_MR = y;
        y_STAR_MC = y_STAR_MR;

        const T* xBuffer = x_MC_STAR.LockedBuffer();
        const T* yBuffer = y_STAR_MC.LockedBuffer();
        const Int incy = y_STAR_MC.LDim();
        if( uplo == LOWER )
        {
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const Int j = rowShift + jLoc*c;
                const Int heightAboveDiag = Length(j,colShift,r);

                const T beta = y_STAR_MR.GetLocal(0,jLoc);
                const T kappa = x_MR_STAR.GetLocal(jLoc,0);
                const T gamma = ( conjugate ? alpha*Conj(beta) : alpha*beta );
                const T delta = ( conjugate ? alpha*Conj(kappa) : alpha*kappa );
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
                const Int j = rowShift + jLoc*c;
                const Int heightToDiag = Length(j+1,colShift,r);

                const T beta = y_STAR_MR.GetLocal(0,jLoc);
                const T kappa = x_MR_STAR.GetLocal(jLoc,0);
                const T gamma = ( conjugate ? alpha*Conj(beta) : alpha*beta );
                const T delta = ( conjugate ? alpha*Conj(kappa) : alpha*kappa );
                T* ACol = A.Buffer(0,jLoc);
                for( Int iLoc=0; iLoc<heightToDiag; ++iLoc )
                    ACol[iLoc] += gamma*xBuffer[iLoc] +
                                  delta*yBuffer[iLoc*incy];
            }
        }
        //--------------------------------------------------------------------//
    }
    else if( y.Width() == 1 )
    {
        DistMatrix<T,STAR,MC> x_STAR_MC(g);
        DistMatrix<T,STAR,MR> x_STAR_MR(g);
        DistMatrix<T,MC,STAR> y_MC_STAR(g);
        DistMatrix<T,MR,STAR> y_MR_STAR(g);

        x_STAR_MC.AlignWith( A );
        x_STAR_MR.AlignWith( A );
        y_MC_STAR.AlignWith( A );
        y_MR_STAR.AlignWith( A );
        //--------------------------------------------------------------------//
        x_STAR_MR = x;
        x_STAR_MC = x_STAR_MR;
        y_MC_STAR = y;
        y_MR_STAR = y_MC_STAR;

        const T* xBuffer = x_STAR_MC.LockedBuffer();
        const T* yBuffer = y_MC_STAR.LockedBuffer();
        const Int incx = x_STAR_MC.LDim();
        if( uplo == LOWER )
        {
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const Int j = rowShift + jLoc*c;
                const Int heightAboveDiag = Length(j,colShift,r);

                const T beta = x_STAR_MR.GetLocal(0,jLoc);
                const T kappa = y_MR_STAR.GetLocal(jLoc,0);
                const T gamma = ( conjugate ? alpha*Conj(beta) : alpha*beta );
                const T delta = ( conjugate ? alpha*Conj(kappa) : alpha*kappa );
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
                const Int j = rowShift + jLoc*c;
                const Int heightToDiag = Length(j+1,colShift,r);

                const T beta = x_STAR_MR.GetLocal(0,jLoc);
                const T kappa = y_MR_STAR.GetLocal(jLoc,0);
                const T gamma = ( conjugate ? alpha*Conj(beta) : alpha*beta );
                const T delta = ( conjugate ? alpha*Conj(kappa) : alpha*kappa );
                T* ACol = A.Buffer(0,jLoc);
                for( Int iLoc=0; iLoc<heightToDiag; ++iLoc )
                    ACol[iLoc] += gamma*xBuffer[iLoc*incx] +
                                  delta*yBuffer[iLoc];
            }
        }
        //--------------------------------------------------------------------//
    }
    else
    {
        DistMatrix<T,STAR,MC> x_STAR_MC(g), y_STAR_MC(g);
        DistMatrix<T,STAR,MR> x_STAR_MR(g), y_STAR_MR(g);

        x_STAR_MC.AlignWith( A );
        x_STAR_MR.AlignWith( A );
        y_STAR_MC.AlignWith( A );
        y_STAR_MR.AlignWith( A );
        //--------------------------------------------------------------------//
        x_STAR_MR = x;
        x_STAR_MC = x_STAR_MR;
        y_STAR_MR = y;
        y_STAR_MC = y_STAR_MR;

        const T* xBuffer = x_STAR_MC.LockedBuffer();
        const T* yBuffer = y_STAR_MC.LockedBuffer();
        const Int incx = x_STAR_MC.LDim();
        const Int incy = y_STAR_MC.LDim();
        if( uplo == LOWER )
        {
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const Int j = rowShift + jLoc*c;
                const Int heightAboveDiag = Length(j,colShift,r);

                const T beta = y_STAR_MR.GetLocal(0,jLoc);
                const T kappa = x_STAR_MR.GetLocal(0,jLoc);
                const T gamma = ( conjugate ? alpha*Conj(beta) : alpha*beta );
                const T delta = ( conjugate ? alpha*Conj(kappa) : alpha*kappa );
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
                const Int j = rowShift + jLoc*c;
                const Int heightToDiag = Length(j+1,colShift,r);

                const T beta = y_STAR_MR.GetLocal(0,jLoc);
                const T kappa = x_STAR_MR.GetLocal(0,jLoc);
                const T gamma = ( conjugate ? alpha*Conj(beta) : alpha*beta );
                const T delta = ( conjugate ? alpha*Conj(kappa) : alpha*kappa );
                T* ACol = A.Buffer(0,jLoc);
                for( Int iLoc=0; iLoc<heightToDiag; ++iLoc )
                    ACol[iLoc] += gamma*xBuffer[iLoc*incx] +
                                  delta*yBuffer[iLoc*incy];
            }
        }
        //--------------------------------------------------------------------//
    }
}

} // namespace elem

#endif // ifndef ELEM_BLAS_SYR2_HPP
