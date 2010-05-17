/*
   Copyright 2009-2010 Jack Poulson

   This file is part of Elemental.

   Elemental is free software: you can redistribute it and/or modify it under
   the terms of the GNU Lesser General Public License as published by the
   Free Software Foundation; either version 3 of the License, or 
   (at your option) any later version.

   Elemental is distributed in the hope that it will be useful, but 
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with Elemental. If not, see <http://www.gnu.org/licenses/>.
*/
#include "elemental/blas_internal.hpp"
using namespace std;
using namespace elemental;

template<typename T>
void
elemental::blas::internal::GemvT
( Orientation orientation,
  T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& x,
  T beta,        DistMatrix<T,MC,MR>& y )
{
#ifndef RELEASE
    PushCallStack("blas::internal::GemvT");
    if( A.GetGrid() != x.GetGrid() || x.GetGrid() != y.GetGrid() )
        throw "{A,x,y} must be distributed over the same grid.";
    if( ( x.Width() != 1 && x.Height() != 1 ) ||
        ( y.Width() != 1 && y.Height() != 1 )   )
    {
        throw "GemvT expects x and y to be vectors.";
    }
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    const int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
    if( A.Height() != xLength || A.Width() != yLength )
    {
        ostringstream msg;
        msg << "Nonconformal GemvT: " << endl
            << "  A ~ " << A.Height() << " x " << A.Width() << endl
            << "  x ~ " << x.Height() << " x " << x.Width() << endl
            << "  y ~ " << y.Height() << " x " << y.Width() << endl;
        const string& s = msg.str();
        throw s.c_str();
    }
#endif
    const Grid& grid = A.GetGrid();
    if( x.Width() == 1 && y.Width() == 1 )
    {
        // Temporary distributions
        DistMatrix<T,MC,Star> x_MC_Star(grid);
        DistMatrix<T,MR,Star> z_MR_Star(grid);
        DistMatrix<T,MR,MC  > z_MR_MC(grid);
        DistMatrix<T,MC,MR  > z(grid);

        // Start the algorithm
        blas::Scal( beta, y );
        x_MC_Star.ConformWith( A );
        z_MR_Star.AlignWith( A );
        z_MR_Star.ResizeTo( A.Width(), 1 );
        z_MR_MC.AlignWith( y );
        z.AlignWith( y );
        //--------------------------------------------------------------------//
        x_MC_Star = x;
        blas::Gemv( orientation,
                    alpha, A.LockedLocalMatrix(), 
                           x_MC_Star.LockedLocalMatrix(),
                    (T)0,  z_MR_Star.LocalMatrix()       );
        z_MR_MC.ReduceScatterFrom( z_MR_Star );
        z = z_MR_MC;
        blas::Axpy( (T)1, z, y );
        //--------------------------------------------------------------------//
        x_MC_Star.FreeConstraints();
        z_MR_Star.FreeConstraints();
        z_MR_MC.FreeConstraints();
        z.FreeConstraints();
    }
    else if( x.Width() == 1 )
    {
        // Temporary distributions
        DistMatrix<T,MC,Star> x_MC_Star(grid);
        DistMatrix<T,MR,Star> z_MR_Star(grid);
        DistMatrix<T,MR,MC  > z_MR_MC(grid);
        DistMatrix<T,MC,MR  > zTrans(grid);

        // Start the algorithm
        blas::Scal( beta, y );
        x_MC_Star.ConformWith( A );
        z_MR_Star.AlignWith( A );
        z_MR_Star.ResizeTo( A.Width(), 1 );
        z_MR_MC.AlignWith( y );
        zTrans.AlignWith( y );
        //--------------------------------------------------------------------//
        x_MC_Star = x;
        blas::Gemv( orientation,
                    alpha, A.LockedLocalMatrix(),
                           x_MC_Star.LockedLocalMatrix(),
                    (T)0,  z_MR_Star.LocalMatrix()       );
        z_MR_MC.ReduceScatterFrom( z_MR_Star );
        blas::Trans( z_MR_MC, zTrans );
        blas::Axpy( (T)1, zTrans, y );
        //--------------------------------------------------------------------//
        x_MC_Star.FreeConstraints();
        z_MR_Star.FreeConstraints();
        z_MR_MC.FreeConstraints();
        zTrans.FreeConstraints();
    }
    else if( y.Width() == 1 )
    {
        // Temporary distributions
        DistMatrix<T,Star,MC  > x_Star_MC(grid);
        DistMatrix<T,MR,  Star> z_MR_Star(grid);
        DistMatrix<T,MR,  MC  > z_MR_MC(grid);
        DistMatrix<T,MC,  MR  > z(grid);

        // Start the algorithm
        blas::Scal( beta, y );
        x_Star_MC.ConformWith( A );
        z_MR_Star.AlignWith( A );
        z_MR_Star.ResizeTo( A.Width(), 1 );
        z_MR_MC.AlignWith( y );
        z.AlignWith( y );
        //--------------------------------------------------------------------//
        x_Star_MC = x;
        blas::Gemv( orientation,
                    alpha, A.LockedLocalMatrix(), 
                           x_Star_MC.LockedLocalMatrix(),
                    (T)0,  z_MR_Star.LocalMatrix()       );
        z_MR_MC.ReduceScatterFrom( z_MR_Star );
        z = z_MR_MC;
        blas::Axpy( (T)1, z, y );
        //--------------------------------------------------------------------//
        x_Star_MC.FreeConstraints();
        z_MR_Star.FreeConstraints();
        z_MR_MC.FreeConstraints();
        z.FreeConstraints();
    }
    else
    {
        // Temporary distributions
        DistMatrix<T,Star,MC  > x_Star_MC(grid);
        DistMatrix<T,MR,  Star> z_MR_Star(grid);
        DistMatrix<T,MR,  MC  > z_MR_MC(grid);
        DistMatrix<T,MC,  MR  > zTrans(grid);

        // Start the algorithm
        blas::Scal( beta, y );
        x_Star_MC.ConformWith( A );
        z_MR_Star.AlignWith( A );
        z_MR_Star.ResizeTo( A.Width(), 1 );
        z_MR_MC.AlignWith( y );
        zTrans.AlignWith( y );
        //--------------------------------------------------------------------//
        x_Star_MC = x;
        blas::Gemv( orientation,
                    alpha, A.LockedLocalMatrix(),
                           x_Star_MC.LockedLocalMatrix(),
                    (T)0,  z_MR_Star.LocalMatrix()       );
        z_MR_MC.ReduceScatterFrom( z_MR_Star );
        blas::Trans( z_MR_MC, zTrans );
        blas::Axpy( (T)1, zTrans, y );
        //--------------------------------------------------------------------//
        x_Star_MC.FreeConstraints();
        z_MR_Star.FreeConstraints();
        z_MR_MC.FreeConstraints();
        zTrans.FreeConstraints();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void elemental::blas::internal::GemvT
( Orientation orientation,
  float alpha, const DistMatrix<float,MC,MR>& A,
               const DistMatrix<float,MC,MR>& x,
  float beta,        DistMatrix<float,MC,MR>& y );

template void elemental::blas::internal::GemvT
( Orientation orientation,
  double alpha, const DistMatrix<double,MC,MR>& A,
                const DistMatrix<double,MC,MR>& x,
  double beta,        DistMatrix<double,MC,MR>& y );

#ifndef WITHOUT_COMPLEX
template void elemental::blas::internal::GemvT
( Orientation orientation,
  scomplex alpha, const DistMatrix<scomplex,MC,MR>& A,
                  const DistMatrix<scomplex,MC,MR>& x,
  scomplex beta,        DistMatrix<scomplex,MC,MR>& y );

template void elemental::blas::internal::GemvT
( Orientation orientation,
  dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& A,
                  const DistMatrix<dcomplex,MC,MR>& x,
  dcomplex beta,        DistMatrix<dcomplex,MC,MR>& y );
#endif

