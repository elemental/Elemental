/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace elem {
namespace internal {

template<typename T>
inline void
GemvT
( Orientation orientation,
  T alpha, const DistMatrix<T>& A,
           const DistMatrix<T>& x,
  T beta,        DistMatrix<T>& y )
{
#ifndef RELEASE
    PushCallStack("internal::GemvT");
    if( A.Grid() != x.Grid() || x.Grid() != y.Grid() )
        throw std::logic_error
        ("{A,x,y} must be distributed over the same grid");
    if( ( x.Width() != 1 && x.Height() != 1 ) ||
        ( y.Width() != 1 && y.Height() != 1 )   )
        throw std::logic_error("GemvT expects x and y to be vectors");
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    const int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
    if( A.Height() != xLength || A.Width() != yLength )
    {
        std::ostringstream msg;
        msg << "Nonconformal GemvT: \n"
            << "  A ~ " << A.Height() << " x " << A.Width() << "\n"
            << "  x ~ " << x.Height() << " x " << x.Width() << "\n"
            << "  y ~ " << y.Height() << " x " << y.Width() << "\n";
        throw std::logic_error( msg.str() );
    }
#endif
    const Grid& g = A.Grid();
    if( x.Width() == 1 && y.Width() == 1 )
    {
        DistMatrix<T,MC,STAR> x_MC_STAR(g);
        DistMatrix<T,MR,STAR> z_MR_STAR(g);
        DistMatrix<T,MR,MC  > z_MR_MC(g);
        DistMatrix<T> z(g);

        Scale( beta, y );
        x_MC_STAR.AlignWith( A );
        z_MR_STAR.AlignWith( A );
        Zeros( A.Width(), 1, z_MR_STAR );
        z_MR_MC.AlignWith( y );
        z.AlignWith( y );
        //--------------------------------------------------------------------//
        x_MC_STAR = x;
        Gemv
        ( orientation,
          alpha, A.LockedLocalMatrix(), 
                 x_MC_STAR.LockedLocalMatrix(),
          T(0),  z_MR_STAR.LocalMatrix() );
        z_MR_MC.SumScatterFrom( z_MR_STAR );
        z = z_MR_MC;
        Axpy( T(1), z, y );
        //--------------------------------------------------------------------//
        x_MC_STAR.FreeAlignments();
        z_MR_STAR.FreeAlignments();
        z_MR_MC.FreeAlignments();
        z.FreeAlignments();
    }
    else if( x.Width() == 1 )
    {
        DistMatrix<T,MC,STAR> x_MC_STAR(g);
        DistMatrix<T,MR,STAR> z_MR_STAR(g);
        DistMatrix<T,MR,MC  > z_MR_MC(g);
        DistMatrix<T> zTrans(g);

        Scale( beta, y );
        x_MC_STAR.AlignWith( A );
        z_MR_STAR.AlignWith( A );
        z_MR_STAR.ResizeTo( A.Width(), 1 );
        z_MR_MC.AlignWith( y );
        zTrans.AlignWith( y );
        //--------------------------------------------------------------------//
        x_MC_STAR = x;
        Gemv
        ( orientation,
          alpha, A.LockedLocalMatrix(),
                 x_MC_STAR.LockedLocalMatrix(),
          T(0),  z_MR_STAR.LocalMatrix() );
        z_MR_MC.SumScatterFrom( z_MR_STAR );
        Transpose( z_MR_MC, zTrans );
        Axpy( T(1), zTrans, y );
        //--------------------------------------------------------------------//
        x_MC_STAR.FreeAlignments();
        z_MR_STAR.FreeAlignments();
        z_MR_MC.FreeAlignments();
        zTrans.FreeAlignments();
    }
    else if( y.Width() == 1 )
    {
        DistMatrix<T,STAR,MC  > x_STAR_MC(g);
        DistMatrix<T,MR,  STAR> z_MR_STAR(g);
        DistMatrix<T,MR,  MC  > z_MR_MC(g);
        DistMatrix<T> z(g);

        Scale( beta, y );
        x_STAR_MC.AlignWith( A );
        z_MR_STAR.AlignWith( A );
        z_MR_STAR.ResizeTo( A.Width(), 1 );
        z_MR_MC.AlignWith( y );
        z.AlignWith( y );
        //--------------------------------------------------------------------//
        x_STAR_MC = x;
        Gemv
        ( orientation,
          alpha, A.LockedLocalMatrix(), 
                 x_STAR_MC.LockedLocalMatrix(),
          T(0),  z_MR_STAR.LocalMatrix() );
        z_MR_MC.SumScatterFrom( z_MR_STAR );
        z = z_MR_MC;
        Axpy( T(1), z, y );
        //--------------------------------------------------------------------//
        x_STAR_MC.FreeAlignments();
        z_MR_STAR.FreeAlignments();
        z_MR_MC.FreeAlignments();
        z.FreeAlignments();
    }
    else
    {
        DistMatrix<T,STAR,MC  > x_STAR_MC(g);
        DistMatrix<T,MR,  STAR> z_MR_STAR(g);
        DistMatrix<T,MR,  MC  > z_MR_MC(g);
        DistMatrix<T> zTrans(g);

        Scale( beta, y );
        x_STAR_MC.AlignWith( A );
        z_MR_STAR.AlignWith( A );
        z_MR_STAR.ResizeTo( A.Width(), 1 );
        z_MR_MC.AlignWith( y );
        zTrans.AlignWith( y );
        //--------------------------------------------------------------------//
        x_STAR_MC = x;
        Gemv
        ( orientation,
          alpha, A.LockedLocalMatrix(),
                 x_STAR_MC.LockedLocalMatrix(),
          T(0),  z_MR_STAR.LocalMatrix() );
        z_MR_MC.SumScatterFrom( z_MR_STAR );
        Transpose( z_MR_MC, zTrans );
        Axpy( T(1), zTrans, y );
        //--------------------------------------------------------------------//
        x_STAR_MC.FreeAlignments();
        z_MR_STAR.FreeAlignments();
        z_MR_MC.FreeAlignments();
        zTrans.FreeAlignments();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
GemvT
( Orientation orientation,
  T alpha, const DistMatrix<T>& A,
           const DistMatrix<T,VC,STAR>& x,
  T beta,        DistMatrix<T,VC,STAR>& y )
{
#ifndef RELEASE
    PushCallStack("internal::GemvT");
    if( A.Grid() != x.Grid() || x.Grid() != y.Grid() )
        throw std::logic_error
        ("{A,x,y} must be distributed over the same grid");
    if( x.Width() != 1 || y.Width() != 1 )
        throw std::logic_error("GemvT expects x and y to be column vectors");
    if( A.Height() != x.Height() || A.Width() != y.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal GemvT: \n"
            << "  A ~ " << A.Height() << " x " << A.Width() << "\n"
            << "  x ~ " << x.Height() << " x " << x.Width() << "\n"
            << "  y ~ " << y.Height() << " x " << y.Width() << "\n";
        throw std::logic_error( msg.str() );
    }
#endif
    const Grid& g = A.Grid();

    DistMatrix<T,MC,STAR> x_MC_STAR(g);
    DistMatrix<T,MR,STAR> z_MR_STAR(g);
    DistMatrix<T,VR,STAR> z_VR_STAR(g);
    DistMatrix<T,VC,STAR> z(g);

    Scale( beta, y );
    x_MC_STAR.AlignWith( A );
    z_MR_STAR.AlignWith( A );
    Zeros( A.Width(), 1, z_MR_STAR );
    z_VR_STAR.AlignWith( A );
    z.AlignWith( y );
    //--------------------------------------------------------------------//
    x_MC_STAR = x;
    Gemv
    ( orientation,
      alpha, A.LockedLocalMatrix(), 
             x_MC_STAR.LockedLocalMatrix(),
      T(0),  z_MR_STAR.LocalMatrix() );
    z_VR_STAR.SumScatterFrom( z_MR_STAR );
    z = z_VR_STAR;
    Axpy( T(1), z, y );
    //--------------------------------------------------------------------//
    x_MC_STAR.FreeAlignments();
    z_MR_STAR.FreeAlignments();
    z_VR_STAR.FreeAlignments();
    z.FreeAlignments();
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace internal
} // namespace elem
