/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace elem {

template<typename F>
inline void
DiagonalSolve
( LeftOrRight side, Orientation orientation,
  const Matrix<F>& d, Matrix<F>& X, bool checkIfSingular )
{
#ifndef RELEASE
    PushCallStack("DiagonalSolve");
#endif
    const int m = X.Height();
    const int n = X.Width();
    const int ldim = X.LDim();
    if( side == LEFT )
    {
        for( int i=0; i<m; ++i )
        {
            const F delta = d.Get(i,0);
            if( checkIfSingular && delta == F(0) )
                throw SingularMatrixException();
            const F deltaInv = F(1)/delta;
            F* XBuffer = X.Buffer(i,0);
            if( orientation == ADJOINT )
                for( int j=0; j<n; ++j )
                    XBuffer[j*ldim] *= Conj(deltaInv);
            else
                for( int j=0; j<n; ++j )
                    XBuffer[j*ldim] *= deltaInv;
        }
    }
    else
    {
        for( int j=0; j<n; ++j )
        {
            const F delta = d.Get(j,0);
            if( checkIfSingular && delta == F(0) )
                throw SingularMatrixException();
            const F deltaInv = F(1)/delta;
            F* XBuffer = X.Buffer(0,j);
            if( orientation == ADJOINT )
                for( int i=0; i<m; ++i )
                    XBuffer[i] *= Conj(deltaInv);
            else
                for( int i=0; i<m; ++i )
                    XBuffer[i] *= deltaInv;
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
DiagonalSolve
( LeftOrRight side, Orientation orientation,
  const Matrix<typename Base<F>::type>& d, Matrix<F>& X,
  bool checkIfSingular )
{
#ifndef RELEASE
    PushCallStack("DiagonalSolve");
#endif
    typedef typename Base<F>::type R;

    const int m = X.Height();
    const int n = X.Width();
    const int ldim = X.LDim();
    if( side == LEFT )
    {
        for( int i=0; i<m; ++i )
        {
            const R delta = d.Get(i,0);
            if( checkIfSingular && delta == R(0) )
                throw SingularMatrixException();
            const R deltaInv = R(1)/delta;
            F* XBuffer = X.Buffer(i,0);
            for( int j=0; j<n; ++j )
                XBuffer[j*ldim] *= deltaInv;
        }
    }
    else
    {
        for( int j=0; j<n; ++j )
        {
            const R delta = d.Get(j,0);
            if( checkIfSingular && delta == R(0) )
                throw SingularMatrixException();
            const R deltaInv = R(1)/delta;
            F* XBuffer = X.Buffer(0,j);
            for( int i=0; i<m; ++i )
                XBuffer[i] *= deltaInv;
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F,Distribution U,Distribution V,
                    Distribution W,Distribution Z>
inline void
DiagonalSolve
( LeftOrRight side, Orientation orientation,
  const DistMatrix<F,U,V>& d, DistMatrix<F,W,Z>& X,
  bool checkIfSingular )
{
#ifndef RELEASE
    PushCallStack("DiagonalSolve");
#endif
    if( side == LEFT )
    {
        if( U == W && V == STAR && d.ColAlignment() == X.ColAlignment() )
        {
            DiagonalSolve
            ( LEFT, orientation, d.LockedLocalMatrix(), X.LocalMatrix(),
              checkIfSingular );
        }
        else
        {
            DistMatrix<F,W,STAR> d_W_STAR( X.Grid() );
            d_W_STAR = d;
            DiagonalSolve
            ( LEFT, orientation,
              d_W_STAR.LockedLocalMatrix(), X.LocalMatrix(), checkIfSingular );
        }
    }
    else
    {
        if( U == Z && V == STAR && d.ColAlignment() == X.RowAlignment() )
        {
            DiagonalSolve
            ( RIGHT, orientation, d.LockedLocalMatrix(), X.LocalMatrix(),
              checkIfSingular );
        }
        else
        {
            DistMatrix<F,Z,STAR> d_Z_STAR( X.Grid() );
            d_Z_STAR = d;
            DiagonalSolve
            ( RIGHT, orientation,
              d_Z_STAR.LockedLocalMatrix(), X.LocalMatrix(), checkIfSingular );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F,Distribution U,Distribution V,
                    Distribution W,Distribution Z>
inline void
DiagonalSolve
( LeftOrRight side, Orientation orientation,
  const DistMatrix<typename Base<F>::type,U,V>& d, DistMatrix<F,W,Z>& X,
  bool checkIfSingular )
{
#ifndef RELEASE
    PushCallStack("DiagonalSolve");
#endif
    typedef typename Base<F>::type R;

    if( side == LEFT )
    {
        if( U == W && V == STAR && d.ColAlignment() == X.ColAlignment() )
        {
            DiagonalSolve
            ( LEFT, orientation, d.LockedLocalMatrix(), X.LocalMatrix(),
              checkIfSingular );
        }
        else
        {
            DistMatrix<R,W,STAR> d_W_STAR( X.Grid() );
            d_W_STAR = d;
            DiagonalSolve
            ( LEFT, orientation,
              d_W_STAR.LockedLocalMatrix(), X.LocalMatrix(), checkIfSingular );
        }
    }
    else
    {
        if( U == Z && V == STAR && d.ColAlignment() == X.RowAlignment() )
        {
            DiagonalSolve
            ( RIGHT, orientation, d.LockedLocalMatrix(), X.LocalMatrix(),
              checkIfSingular );
        }
        else
        {
            DistMatrix<R,Z,STAR> d_Z_STAR( X.Grid() );
            d_Z_STAR = d;
            DiagonalSolve
            ( RIGHT, orientation,
              d_Z_STAR.LockedLocalMatrix(), X.LocalMatrix(), checkIfSingular );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
