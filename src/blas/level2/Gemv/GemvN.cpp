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
elemental::blas::internal::GemvN
( T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& x,
  T beta,        DistMatrix<T,MC,MR>& y )
{
#ifndef RELEASE
    PushCallStack("blas::internal::GemvN");
    if( A.GetGrid() != x.GetGrid() || x.GetGrid() != y.GetGrid() )
        throw "{A,x,y} must be distributed over the same grid.";
    if( ( x.Width() != 1 && x.Height() != 1 ) ||
        ( y.Width() != 1 && y.Height() != 1 )   )
    {
        throw "x and y are assumed to be vectors.";
    }
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    const int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
    if( A.Height() != yLength || A.Width() != xLength )
    {
        ostringstream msg;
        msg << "Nonconformal GemvN: " << endl
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
        DistMatrix<T,MR,Star> x_MR_Star(grid);
        DistMatrix<T,MC,Star> z_MC_Star(grid);

        // Start the algorithm
        blas::Scal( beta, y );
        x_MR_Star.ConformWith( A );
        z_MC_Star.AlignWith( A );
        z_MC_Star.ResizeTo( A.Height(), 1 );
        //--------------------------------------------------------------------//
        x_MR_Star = x;
        blas::Gemv( Normal,
                    alpha, A.LockedLocalMatrix(), 
                           x_MR_Star.LockedLocalMatrix(),
                    (T)0,  z_MC_Star.LocalMatrix()       );
        y.ReduceScatterUpdate( (T)1, z_MC_Star );
        //--------------------------------------------------------------------//
        x_MR_Star.FreeConstraints();
        z_MC_Star.FreeConstraints();
    }
    else if( x.Width() == 1 )
    {
        // Temporary distributions
        DistMatrix<T,MR,Star> x_MR_Star(grid);
        DistMatrix<T,MC,Star> z_MC_Star(grid);
        DistMatrix<T,MC,MR  > z(grid);
        DistMatrix<T,MC,MR  > zTrans(grid);

        // Start the algorithm
        blas::Scal( beta, y );
        x_MR_Star.ConformWith( A );
        z_MC_Star.AlignWith( A );
        z_MC_Star.ResizeTo( A.Height(), 1 );
        z.AlignWith( y );
        zTrans.AlignWith( y );
        //--------------------------------------------------------------------//
        x_MR_Star = x;
        blas::Gemv( Normal,
                    alpha, A.LockedLocalMatrix(),
                           x_MR_Star.LockedLocalMatrix(),
                    (T)0,  z_MC_Star.LocalMatrix()       );
        z.ReduceScatterFrom( z_MC_Star );
        blas::Trans( z, zTrans );
        blas::Axpy( (T)1, zTrans, y );
        //--------------------------------------------------------------------//
        x_MR_Star.FreeConstraints();
        z_MC_Star.FreeConstraints();
        z.FreeConstraints();
        zTrans.FreeConstraints();
    }
    else if( y.Width() == 1 )
    {
        // Temporary distributions
        DistMatrix<T,Star,MR  > x_Star_MR(grid);
        DistMatrix<T,MC,  Star> z_MC_Star(grid);

        // Start the algorithm
        blas::Scal( beta, y );
        x_Star_MR.ConformWith( A );
        z_MC_Star.AlignWith( A );
        z_MC_Star.ResizeTo( A.Height(), 1 );
        //--------------------------------------------------------------------//
        x_Star_MR = x;
        blas::Gemv( Normal,
                    alpha, A.LockedLocalMatrix(), 
                           x_Star_MR.LockedLocalMatrix(),
                    (T)0,  z_MC_Star.LocalMatrix()       );
        y.ReduceScatterUpdate( (T)1, z_MC_Star );
        //--------------------------------------------------------------------//
        x_Star_MR.FreeConstraints();
        z_MC_Star.FreeConstraints();
    }
    else
    {
        // Temporary distributions
        DistMatrix<T,Star,MR  > x_Star_MR(grid);
        DistMatrix<T,MC,  Star> z_MC_Star(grid);
        DistMatrix<T,MC,  MR  > z(grid);
        DistMatrix<T,MC,  MR  > zTrans(grid);

        // Start the algorithm
        blas::Scal( beta, y );
        x_Star_MR.ConformWith( A );
        z_MC_Star.AlignWith( A );
        z_MC_Star.ResizeTo( A.Height(), 1 );
        z.AlignWith( y );
        zTrans.AlignWith( y );
        //--------------------------------------------------------------------//
        x_Star_MR = x;
        blas::Gemv( Normal,
                    alpha, A.LockedLocalMatrix(),
                           x_Star_MR.LockedLocalMatrix(),
                    (T)0,  z_MC_Star.LocalMatrix()       );
        z.ReduceScatterFrom( z_MC_Star );
        blas::Trans( z, zTrans );
        blas::Axpy( (T)1, zTrans, y );
        //--------------------------------------------------------------------//
        x_Star_MR.FreeConstraints();
        z_MC_Star.FreeConstraints();
        z.FreeConstraints();
        zTrans.FreeConstraints();
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

