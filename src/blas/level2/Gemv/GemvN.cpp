/*
   This file is part of elemental, a library for distributed-memory dense 
   linear algebra.

   Copyright (C) 2009-2010 Jack Poulson <jack.poulson@gmail.com>

   This program is released under the terms of the license contained in the 
   file LICENSE.
*/
#include "elemental/blas_internal.hpp"
using namespace std;
using namespace elemental;

template<typename T>
void
elemental::blas::internal::GemvN
( T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& x,
  T beta,        DistMatrix<T,MC,MR>& y )
{
#ifndef RELEASE
    PushCallStack("blas::internal::GemvN");
    if( A.GetGrid() != x.GetGrid() || x.GetGrid() != y.GetGrid() )
        throw logic_error( "{A,x,y} must be distributed over the same grid." );
    if( ( x.Width() != 1 && x.Height() != 1 ) ||
        ( y.Width() != 1 && y.Height() != 1 )   )
        throw logic_error( "x and y are assumed to be vectors." );
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    const int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
    if( A.Height() != yLength || A.Width() != xLength )
    {
        ostringstream msg;
        msg << "Nonconformal GemvN: " << endl
            << "  A ~ " << A.Height() << " x " << A.Width() << endl
            << "  x ~ " << x.Height() << " x " << x.Width() << endl 
            << "  y ~ " << y.Height() << " x " << y.Width() << endl;
        throw logic_error( msg.str() );
    }
#endif
    const Grid& g = A.GetGrid();
    if( x.Width() == 1 && y.Width() == 1 )
    {
        // Temporary distributions
        DistMatrix<T,MR,Star> x_MR_Star(g);
        DistMatrix<T,MC,Star> z_MC_Star(g);

        // Start the algorithm
        blas::Scal( beta, y );
        x_MR_Star.AlignWith( A );
        z_MC_Star.AlignWith( A );
        z_MC_Star.ResizeTo( A.Height(), 1 );
        //--------------------------------------------------------------------//
        x_MR_Star = x;
        blas::Gemv
        ( Normal,
          alpha, A.LockedLocalMatrix(), 
                 x_MR_Star.LockedLocalMatrix(),
          (T)0,  z_MC_Star.LocalMatrix() );
        y.SumScatterUpdate( (T)1, z_MC_Star );
        //--------------------------------------------------------------------//
        x_MR_Star.FreeAlignments();
        z_MC_Star.FreeAlignments();
    }
    else if( x.Width() == 1 )
    {
        // Temporary distributions
        DistMatrix<T,MR,Star> x_MR_Star(g);
        DistMatrix<T,MC,Star> z_MC_Star(g);
        DistMatrix<T,MC,MR  > z(g);
        DistMatrix<T,MC,MR  > zTrans(g);

        // Start the algorithm
        blas::Scal( beta, y );
        x_MR_Star.AlignWith( A );
        z_MC_Star.AlignWith( A );
        z_MC_Star.ResizeTo( A.Height(), 1 );
        z.AlignWith( y );
        zTrans.AlignWith( y );
        //--------------------------------------------------------------------//
        x_MR_Star = x;
        blas::Gemv
        ( Normal,
          alpha, A.LockedLocalMatrix(),
                 x_MR_Star.LockedLocalMatrix(),
          (T)0,  z_MC_Star.LocalMatrix() );
        z.SumScatterFrom( z_MC_Star );
        blas::Trans( z, zTrans );
        blas::Axpy( (T)1, zTrans, y );
        //--------------------------------------------------------------------//
        x_MR_Star.FreeAlignments();
        z_MC_Star.FreeAlignments();
        z.FreeAlignments();
        zTrans.FreeAlignments();
    }
    else if( y.Width() == 1 )
    {
        // Temporary distributions
        DistMatrix<T,Star,MR  > x_Star_MR(g);
        DistMatrix<T,MC,  Star> z_MC_Star(g);

        // Start the algorithm
        blas::Scal( beta, y );
        x_Star_MR.AlignWith( A );
        z_MC_Star.AlignWith( A );
        z_MC_Star.ResizeTo( A.Height(), 1 );
        //--------------------------------------------------------------------//
        x_Star_MR = x;
        blas::Gemv
        ( Normal,
          alpha, A.LockedLocalMatrix(), 
                 x_Star_MR.LockedLocalMatrix(),
          (T)0,  z_MC_Star.LocalMatrix() );
        y.SumScatterUpdate( (T)1, z_MC_Star );
        //--------------------------------------------------------------------//
        x_Star_MR.FreeAlignments();
        z_MC_Star.FreeAlignments();
    }
    else
    {
        // Temporary distributions
        DistMatrix<T,Star,MR  > x_Star_MR(g);
        DistMatrix<T,MC,  Star> z_MC_Star(g);
        DistMatrix<T,MC,  MR  > z(g);
        DistMatrix<T,MC,  MR  > zTrans(g);

        // Start the algorithm
        blas::Scal( beta, y );
        x_Star_MR.AlignWith( A );
        z_MC_Star.AlignWith( A );
        z_MC_Star.ResizeTo( A.Height(), 1 );
        z.AlignWith( y );
        zTrans.AlignWith( y );
        //--------------------------------------------------------------------//
        x_Star_MR = x;
        blas::Gemv
        ( Normal,
          alpha, A.LockedLocalMatrix(),
                 x_Star_MR.LockedLocalMatrix(),
          (T)0,  z_MC_Star.LocalMatrix() );
        z.SumScatterFrom( z_MC_Star );
        blas::Trans( z, zTrans );
        blas::Axpy( (T)1, zTrans, y );
        //--------------------------------------------------------------------//
        x_Star_MR.FreeAlignments();
        z_MC_Star.FreeAlignments();
        z.FreeAlignments();
        zTrans.FreeAlignments();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void elemental::blas::internal::GemvN
( float alpha, const DistMatrix<float,MC,MR>& A,
               const DistMatrix<float,MC,MR>& x,
  float beta,        DistMatrix<float,MC,MR>& y );

template void elemental::blas::internal::GemvN
( double alpha, const DistMatrix<double,MC,MR>& A,
                const DistMatrix<double,MC,MR>& x,
  double beta,        DistMatrix<double,MC,MR>& y );

#ifndef WITHOUT_COMPLEX
template void elemental::blas::internal::GemvN
( scomplex alpha, const DistMatrix<scomplex,MC,MR>& A,
                  const DistMatrix<scomplex,MC,MR>& x,
  scomplex beta,        DistMatrix<scomplex,MC,MR>& y );

template void elemental::blas::internal::GemvN
( dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& A,
                  const DistMatrix<dcomplex,MC,MR>& x,
  dcomplex beta,        DistMatrix<dcomplex,MC,MR>& y );
#endif

