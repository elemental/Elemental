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
elemental::blas::internal::SymvColAccumulateU
( T alpha, const DistMatrix<T,MC,MR  >& A,
           const DistMatrix<T,MC,Star>& x_MC_Star,
           const DistMatrix<T,MR,Star>& x_MR_Star,
                 DistMatrix<T,MC,Star>& z_MC_Star,
                 DistMatrix<T,MR,Star>& z_MR_Star )
{
#ifndef RELEASE
    PushCallStack("blas::internal::SymvColAccumulateU");
    if( A.GetGrid() != x_MC_Star.GetGrid() ||
        x_MC_Star.GetGrid() != x_MR_Star.GetGrid() ||
        x_MR_Star.GetGrid() != z_MC_Star.GetGrid() ||
        z_MC_Star.GetGrid() != z_MR_Star.GetGrid()   )
    {
        throw "{A,x,z} must be distributed over the same grid.";
    }
    if( x_MC_Star.Width() != 1 || x_MR_Star.Width() != 1 ||
        z_MC_Star.Width() != 1 || z_MR_Star.Width() != 1   )
    {
        throw "Expected x and z to be column vectors.";
    }
    if( A.Height() != A.Width() || 
        A.Height() != x_MC_Star.Height() ||
        A.Height() != x_MR_Star.Height() ||
        A.Height() != z_MC_Star.Height() ||
        A.Height() != z_MR_Star.Height()   )
    {
        ostringstream msg;
        msg << "Nonconformal SymvColAccumulateU: " << endl
            << "  A ~ " << A.Height() << " x " << A.Width() << endl
            << "  x[MC,* ] ~ " << x_MC_Star.Height() << " x " 
                               << x_MC_Star.Width() << endl
            << "  x[MR,* ] ~ " << x_MR_Star.Height() << " x " 
                               << x_MR_Star.Width() << endl
            << "  z[MC,* ] ~ " << z_MC_Star.Height() << " x " 
                               << z_MC_Star.Width() << endl
            << "  z[MR,* ] ~ " << z_MR_Star.Height() << " x " 
                               << z_MR_Star.Width() << endl;
        const string& s = msg.str();
        throw s.c_str();
    }
    if( x_MC_Star.ColAlignment() != A.ColAlignment() ||
        x_MR_Star.ColAlignment() != A.RowAlignment() ||
        z_MC_Star.ColAlignment() != A.ColAlignment() ||
        z_MR_Star.ColAlignment() != A.RowAlignment()   )
    {
        throw "Partial matrix distributions are misaligned.";
    }
#endif
    const Grid& grid = A.GetGrid();

    // Matrix views
    DistMatrix<T,MC,MR> 
        ATL(grid), ATR(grid),  A00(grid), A01(grid), A02(grid),
        ABL(grid), ABR(grid),  A10(grid), A11(grid), A12(grid),
                               A20(grid), A21(grid), A22(grid);

    DistMatrix<T,MC,MR> D11(grid);

    DistMatrix<T,MC,Star> 
        xT_MC_Star(grid),  x0_MC_Star(grid),
        xB_MC_Star(grid),  x1_MC_Star(grid),
                           x2_MC_Star(grid);

    DistMatrix<T,MR,Star> 
        xT_MR_Star(grid),  x0_MR_Star(grid),
        xB_MR_Star(grid),  x1_MR_Star(grid),
                           x2_MR_Star(grid);

    DistMatrix<T,MC,Star> 
        zT_MC_Star(grid),  z0_MC_Star(grid),
        zB_MC_Star(grid),  z1_MC_Star(grid),
                           z2_MC_Star(grid);

    DistMatrix<T,MR,Star> 
        zT_MR_Star(grid),  z0_MR_Star(grid),
        zB_MR_Star(grid),  z1_MR_Star(grid),
                           z2_MR_Star(grid);

    // We want our local gemvs to be of width blocksize, so we will 
    // temporarily change to max(r,c) times the current blocksize
    const int ratio = max( grid.Height(), grid.Width() );
    PushBlocksizeStack( ratio*Blocksize() );
                 
    LockedPartitionDownDiagonal( A, ATL, ATR,
                                    ABL, ABR );
    LockedPartitionDown( x_MC_Star, xT_MC_Star,
                                    xB_MC_Star );
    LockedPartitionDown( x_MR_Star, xT_MR_Star,
                                    xB_MR_Star );
    PartitionDown( z_MC_Star, zT_MC_Star,
                              zB_MC_Star );
    PartitionDown( z_MR_Star, zT_MR_Star,
                              zB_MR_Star );
    while( ATL.Height() < A.Height() )
    {
        LockedRepartitionDownDiagonal( ATL, /**/ ATR,  A00, /**/ A01, A02,
                                      /*************/ /******************/
                                            /**/       A10, /**/ A11, A12,
                                       ABL, /**/ ABR,  A20, /**/ A21, A22 );

        LockedRepartitionDown( xT_MC_Star,  x0_MC_Star,
                              /**********/ /**********/
                                            x1_MC_Star,
                               xB_MC_Star,  x2_MC_Star );

        LockedRepartitionDown( xT_MR_Star,  x0_MR_Star,
                              /**********/ /**********/
                                            x1_MR_Star,
                               xB_MR_Star,  x2_MR_Star );

        RepartitionDown( zT_MC_Star,  z0_MC_Star,
                        /**********/ /**********/
                                      z1_MC_Star,
                         zB_MC_Star,  z2_MC_Star );

        RepartitionDown( zT_MR_Star,  z0_MR_Star,
                        /**********/ /**********/
                                      z1_MR_Star,
                         zB_MR_Star,  z2_MR_Star );

        D11.AlignWith( A11 );
        //--------------------------------------------------------------------//
        D11 = A11;
        D11.MakeTrapezoidal( Left, Upper );
        blas::Gemv( Normal, 
                    alpha, D11.LockedLocalMatrix(), 
                           x1_MR_Star.LockedLocalMatrix(),
                    (T)1,  z1_MC_Star.LocalMatrix()       );
        D11.MakeTrapezoidal( Left, Upper, 1 );
        blas::Gemv( Transpose,
                    alpha, D11.LockedLocalMatrix(),
                           x1_MC_Star.LockedLocalMatrix(),
                    (T)1,  z1_MR_Star.LocalMatrix()       );

        blas::Gemv( Normal,
                    alpha, A12.LockedLocalMatrix(),
                           x2_MR_Star.LockedLocalMatrix(),
                    (T)1,  z1_MC_Star.LocalMatrix()       );
        blas::Gemv( Transpose,
                    alpha, A12.LockedLocalMatrix(),
                           x1_MC_Star.LockedLocalMatrix(),
                    (T)1,  z2_MR_Star.LocalMatrix()       );
        //--------------------------------------------------------------------//
        D11.FreeConstraints();

        SlideLockedPartitionDownDiagonal( ATL, /**/ ATR,  A00, A01, /**/ A02,
                                               /**/       A10, A11, /**/ A12,
                                         /*************/ /******************/
                                          ABL, /**/ ABR,  A20, A21, /**/ A22 );

        SlideLockedPartitionDown( xT_MC_Star,  x0_MC_Star,
                                               x1_MC_Star,
                                 /**********/ /**********/
                                  xB_MC_Star,  x2_MC_Star );

        SlideLockedPartitionDown( xT_MR_Star,  x0_MR_Star,
                                               x1_MR_Star,
                                 /**********/ /**********/
                                  xB_MR_Star,  x2_MR_Star );
        
        SlidePartitionDown( zT_MC_Star,  z0_MC_Star,
                                         z1_MC_Star,
                           /**********/ /**********/
                            zB_MC_Star,  z2_MC_Star );

        SlidePartitionDown( zT_MR_Star,  z0_MR_Star,
                                         z1_MR_Star,
                           /**********/ /**********/
                            zB_MR_Star,  z2_MR_Star );
    }

    PopBlocksizeStack();

#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::blas::internal::SymvRowAccumulateU
( T alpha, const DistMatrix<T,MC,  MR>& A,
           const DistMatrix<T,Star,MC>& x_Star_MC,
           const DistMatrix<T,Star,MR>& x_Star_MR,
                 DistMatrix<T,Star,MC>& z_Star_MC,
                 DistMatrix<T,Star,MR>& z_Star_MR )
{
#ifndef RELEASE
    PushCallStack("blas::internal::SymvRowAccumulateU");
    if( A.GetGrid() != x_Star_MC.GetGrid() ||
        x_Star_MC.GetGrid() != x_Star_MR.GetGrid() ||
        x_Star_MR.GetGrid() != z_Star_MC.GetGrid() ||
        z_Star_MC.GetGrid() != z_Star_MR.GetGrid()   )
    {
        throw "{A,x,z} must be distributed over the same grid.";
    }
    if( x_Star_MC.Height() != 1 || x_Star_MR.Height() != 1 ||
        z_Star_MC.Height() != 1 || z_Star_MR.Height() != 1   )
    {
        throw "Expected x and z to be row vectors.";
    }
    if( A.Height() != A.Width() || 
        A.Height() != x_Star_MC.Width() ||
        A.Height() != x_Star_MR.Width() ||
        A.Height() != z_Star_MC.Width() ||
        A.Height() != z_Star_MR.Width()   )
    {
        ostringstream msg;
        msg << "Nonconformal SymvRowAccumulateU: " << endl
            << "  A ~ " << A.Height() << " x " << A.Width() << endl
            << "  x[* ,MC] ~ " << x_Star_MC.Height() << " x " 
                               << x_Star_MC.Width() << endl
            << "  x[* ,MR] ~ " << x_Star_MR.Height() << " x " 
                               << x_Star_MR.Width() << endl
            << "  z[* ,MC] ~ " << z_Star_MC.Height() << " x " 
                               << z_Star_MC.Width() << endl
            << "  z[* ,MR] ~ " << z_Star_MR.Height() << " x " 
                               << z_Star_MR.Width() << endl;
        const string& s = msg.str();
        throw s.c_str();
    }
    if( x_Star_MC.RowAlignment() != A.ColAlignment() ||
        x_Star_MR.RowAlignment() != A.RowAlignment() ||
        z_Star_MC.RowAlignment() != A.ColAlignment() ||
        z_Star_MR.RowAlignment() != A.RowAlignment()   )
    {
        throw "Partial matrix distributions are misaligned.";
    }
#endif
    const Grid& grid = A.GetGrid();

    // Matrix views
    DistMatrix<T,MC,MR> 
        ATL(grid), ATR(grid),  A00(grid), A01(grid), A02(grid),
        ABL(grid), ABR(grid),  A10(grid), A11(grid), A12(grid),
                               A20(grid), A21(grid), A22(grid);

    DistMatrix<T,MC,MR> D11(grid);

    DistMatrix<T,Star,MC> 
        xL_Star_MC(grid), xR_Star_MC(grid),
        x0_Star_MC(grid), x1_Star_MC(grid), x2_Star_MC(grid);

    DistMatrix<T,Star,MR> 
        xL_Star_MR(grid), xR_Star_MR(grid),
        x0_Star_MR(grid), x1_Star_MR(grid), x2_Star_MR(grid);

    DistMatrix<T,Star,MC> 
        zL_Star_MC(grid), zR_Star_MC(grid),
        z0_Star_MC(grid), z1_Star_MC(grid), z2_Star_MC(grid);

    DistMatrix<T,Star,MR> 
        zL_Star_MR(grid), zR_Star_MR(grid),
        z0_Star_MR(grid), z1_Star_MR(grid), z2_Star_MR(grid);

    // We want our local gemvs to be of width blocksize, so we will 
    // temporarily change to max(r,c) times the current blocksize
    const int ratio = max( grid.Height(), grid.Width() );
    PushBlocksizeStack( ratio*Blocksize() );
                 
    LockedPartitionDownDiagonal( A, ATL, ATR,
                                    ABL, ABR );
    LockedPartitionRight( x_Star_MC,  xL_Star_MC, xR_Star_MC );
    LockedPartitionRight( x_Star_MR,  xL_Star_MR, xR_Star_MR );
    PartitionRight( z_Star_MC,  zL_Star_MC, zR_Star_MC );
    PartitionRight( z_Star_MR,  zL_Star_MR, zR_Star_MR );
    while( ATL.Height() < A.Height() )
    {
        LockedRepartitionDownDiagonal( ATL, /**/ ATR,  A00, /**/ A01, A02,
                                      /*************/ /******************/
                                            /**/       A10, /**/ A11, A12,
                                       ABL, /**/ ABR,  A20, /**/ A21, A22 );

        LockedRepartitionRight( xL_Star_MC, /**/ xR_Star_MC, 
                                x0_Star_MC, /**/ x1_Star_MC, x2_Star_MC );

        LockedRepartitionRight( xL_Star_MR, /**/ xR_Star_MR, 
                                x0_Star_MR, /**/ x1_Star_MR, x2_Star_MR );

        RepartitionRight( zL_Star_MC, /**/ zR_Star_MC,
                          z0_Star_MC, /**/ z1_Star_MC, z2_Star_MC );

        RepartitionRight( zL_Star_MR, /**/ zR_Star_MR,
                          z0_Star_MR, /**/ z1_Star_MR, z2_Star_MR );

        D11.AlignWith( A11 );
        //--------------------------------------------------------------------//
        D11 = A11;
        D11.MakeTrapezoidal( Left, Upper );
        blas::Gemv( Normal, 
                    alpha, D11.LockedLocalMatrix(), 
                           x1_Star_MR.LockedLocalMatrix(),
                    (T)1,  z1_Star_MC.LocalMatrix()       );
        D11.MakeTrapezoidal( Left, Upper, 1 );
        blas::Gemv( Transpose,
                    alpha, D11.LockedLocalMatrix(),
                           x1_Star_MC.LockedLocalMatrix(),
                    (T)1,  z1_Star_MR.LocalMatrix()       );

        blas::Gemv( Normal,
                    alpha, A12.LockedLocalMatrix(),
                           x2_Star_MR.LockedLocalMatrix(),
                    (T)1,  z1_Star_MC.LocalMatrix()       );
        blas::Gemv( Transpose,
                    alpha, A12.LockedLocalMatrix(),
                           x1_Star_MC.LockedLocalMatrix(),
                    (T)1,  z2_Star_MR.LocalMatrix()       );
        //--------------------------------------------------------------------//
        D11.FreeConstraints();

        SlideLockedPartitionDownDiagonal( ATL, /**/ ATR,  A00, A01, /**/ A02,
                                               /**/       A10, A11, /**/ A12,
                                         /*************/ /******************/
                                          ABL, /**/ ABR,  A20, A21, /**/ A22 );

        SlideLockedPartitionRight( xL_Star_MC,             /**/ xR_Star_MC,
                                   x0_Star_MC, x1_Star_MC, /**/ x2_Star_MC );
        
        SlideLockedPartitionRight( xL_Star_MR,             /**/ xR_Star_MR,
                                   x0_Star_MR, x1_Star_MR, /**/ x2_Star_MR );

        SlidePartitionRight( zL_Star_MC,             /**/ zR_Star_MC,
                             z0_Star_MC, z1_Star_MC, /**/ z2_Star_MC );

        SlidePartitionRight( zL_Star_MR,             /**/ zR_Star_MR,
                             z0_Star_MR, z1_Star_MR, /**/ z2_Star_MR );
    }

    PopBlocksizeStack();

#ifndef RELEASE
    PopCallStack();
#endif
}

template void elemental::blas::internal::SymvColAccumulateU
( float alpha, const DistMatrix<float,MC,MR  >& A,
               const DistMatrix<float,MC,Star>& x_MC_Star,    
               const DistMatrix<float,MR,Star>& x_MR_Star,
                     DistMatrix<float,MC,Star>& z_MC_Star,
                     DistMatrix<float,MR,Star>& z_MR_Star );

template void elemental::blas::internal::SymvRowAccumulateU
( float alpha, const DistMatrix<float,MC,  MR>& A,
               const DistMatrix<float,Star,MC>& x_Star_MC,
               const DistMatrix<float,Star,MR>& x_Star_MR,
                     DistMatrix<float,Star,MC>& z_Star_MC,
                     DistMatrix<float,Star,MR>& z_Star_MR );

template void elemental::blas::internal::SymvColAccumulateU
( double alpha, const DistMatrix<double,MC,MR  >& A,
                const DistMatrix<double,MC,Star>& x_MC_Star,
                const DistMatrix<double,MR,Star>& x_MR_Star,
                      DistMatrix<double,MC,Star>& z_MC_Star,
                      DistMatrix<double,MR,Star>& z_MR_Star );

template void elemental::blas::internal::SymvRowAccumulateU
( double alpha, const DistMatrix<double,MC,  MR>& A,
                const DistMatrix<double,Star,MC>& x_Star_MC,
                const DistMatrix<double,Star,MR>& x_Star_MR,
                      DistMatrix<double,Star,MC>& z_Star_MC,
                      DistMatrix<double,Star,MR>& z_Star_MR );

#ifndef WITHOUT_COMPLEX
template void elemental::blas::internal::SymvColAccumulateU
( scomplex alpha,
  const DistMatrix<scomplex,MC,MR  >& A,
  const DistMatrix<scomplex,MC,Star>& x_MC_Star,
  const DistMatrix<scomplex,MR,Star>& x_MR_Star,
        DistMatrix<scomplex,MC,Star>& z_MC_Star,
        DistMatrix<scomplex,MR,Star>& z_MR_Star );

template void elemental::blas::internal::SymvRowAccumulateU
( scomplex alpha,
  const DistMatrix<scomplex,MC,  MR>& A,
  const DistMatrix<scomplex,Star,MC>& x_Star_MC,
  const DistMatrix<scomplex,Star,MR>& x_Star_MR,
        DistMatrix<scomplex,Star,MC>& z_Star_MC,
        DistMatrix<scomplex,Star,MR>& z_Star_MR );

template void elemental::blas::internal::SymvColAccumulateU
( dcomplex alpha,
  const DistMatrix<dcomplex,MC,MR  >& A,
  const DistMatrix<dcomplex,MC,Star>& x_MC_Star,
  const DistMatrix<dcomplex,MR,Star>& x_MR_Star,
        DistMatrix<dcomplex,MC,Star>& z_MC_Star,
        DistMatrix<dcomplex,MR,Star>& z_MR_Star );

template void elemental::blas::internal::SymvRowAccumulateU
( dcomplex alpha,
  const DistMatrix<dcomplex,MC,  MR>& A,
  const DistMatrix<dcomplex,Star,MC>& x_Star_MC,
  const DistMatrix<dcomplex,Star,MR>& x_Star_MR,
        DistMatrix<dcomplex,Star,MC>& z_Star_MC,
        DistMatrix<dcomplex,Star,MR>& z_Star_MR );
#endif

