/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace elem {

template<typename T>
inline void
Her( UpperOrLower uplo, T alpha, const Matrix<T>& x, Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("Her");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    if( x.Width() != 1 && x.Height() != 1 )
        throw std::logic_error("x must be a vector");
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    if( xLength != A.Height() )
        throw std::logic_error("x must conform with A");
#endif
    const char uploChar = UpperOrLowerToChar( uplo );
    const int m = A.Height();
    const int incx = ( x.Width()==1 ? 1 : x.LDim() );
    blas::Her
    ( uploChar, m, alpha, x.LockedBuffer(), incx, A.Buffer(), A.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Her
( UpperOrLower uplo,
  T alpha, const DistMatrix<T>& x,
                 DistMatrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("Her");
    if( A.Grid() != x.Grid() )
        throw std::logic_error("{A,x} must be distributed over the same grid");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    if( A.Height() != xLength )
    {
        std::ostringstream msg;
        msg << "A must conform with x: \n"
            << "  A ~ " << A.Height() << " x " << A.Width() << "\n"
            << "  x ~ " << x.Height() << " x " << x.Width() << "\n";
        throw std::logic_error( msg.str() );
    }
#endif
    const Grid& g = A.Grid();

    const int localHeight = A.LocalHeight();
    const int localWidth = A.LocalWidth();
    const int r = g.Height();
    const int c = g.Width();
    const int colShift = A.ColShift();
    const int rowShift = A.RowShift();

    if( x.Width() == 1 )
    {
        DistMatrix<T,MC,STAR> x_MC_STAR(g);
        DistMatrix<T,MR,STAR> x_MR_STAR(g);

        x_MC_STAR.AlignWith( A );
        x_MR_STAR.AlignWith( A );
        //--------------------------------------------------------------------//
        x_MC_STAR = x;
        x_MR_STAR = x_MC_STAR;

        const T* xLocal = x_MC_STAR.LockedLocalBuffer();
        if( uplo == LOWER )
        {
            for( int jLocal=0; jLocal<localWidth; ++jLocal )
            {
                const int j = rowShift + jLocal*c;
                const int heightAboveDiag = LocalLength(j,colShift,r);

                const T gamma = alpha*Conj(x_MR_STAR.GetLocal(jLocal,0));
                T* ALocalCol = A.LocalBuffer(0,jLocal);
                for( int iLocal=heightAboveDiag; iLocal<localHeight; ++iLocal )
                    ALocalCol[iLocal] += gamma*xLocal[iLocal];
            }
        }
        else
        {
            for( int jLocal=0; jLocal<localWidth; ++jLocal )
            {
                const int j = rowShift + jLocal*c;
                const int heightToDiag = LocalLength(j+1,colShift,r);

                const T gamma = alpha*Conj(x_MR_STAR.GetLocal(jLocal,0));
                T* ALocalCol = A.LocalBuffer(0,jLocal);
                for( int iLocal=0; iLocal<heightToDiag; ++iLocal )
                    ALocalCol[iLocal] += gamma*xLocal[iLocal];
            }
        }
        //--------------------------------------------------------------------//
        x_MC_STAR.FreeAlignments();
        x_MR_STAR.FreeAlignments();
    }
    else
    {
        DistMatrix<T,STAR,MC> x_STAR_MC(g);
        DistMatrix<T,STAR,MR> x_STAR_MR(g);

        x_STAR_MC.AlignWith( A );
        x_STAR_MR.AlignWith( A );
        //--------------------------------------------------------------------//
        x_STAR_MR = x;
        x_STAR_MC = x_STAR_MR;

        const T* xLocal = x_STAR_MC.LockedLocalBuffer();
        const int incx = x_STAR_MC.LocalLDim();
        if( uplo == LOWER )
        {
            for( int jLocal=0; jLocal<localWidth; ++jLocal )
            {
                const int j = rowShift + jLocal*c;
                const int heightAboveDiag = LocalLength(j,colShift,r);

                const T gamma = alpha*Conj(x_STAR_MR.GetLocal(0,jLocal));
                T* ALocalCol = A.LocalBuffer(0,jLocal);
                for( int iLocal=heightAboveDiag; iLocal<localHeight; ++iLocal )
                    ALocalCol[iLocal] += gamma*xLocal[iLocal*incx];
            }
        }
        else
        {
            for( int jLocal=0; jLocal<localWidth; ++jLocal )
            {
                const int j = rowShift + jLocal*c;
                const int heightToDiag = LocalLength(j+1,colShift,r);

                const T gamma = alpha*Conj(x_STAR_MR.GetLocal(0,jLocal));
                T* ALocalCol = A.LocalBuffer(0,jLocal);
                for( int iLocal=0; iLocal<heightToDiag; ++iLocal )
                    ALocalCol[iLocal] += gamma*xLocal[iLocal*incx];
            }
        }
        //--------------------------------------------------------------------//
        x_STAR_MC.FreeAlignments();
        x_STAR_MR.FreeAlignments();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
