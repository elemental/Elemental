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
#include "Elemental/BLASInternal.hpp"
using namespace std;
using namespace Elemental;

namespace {

#ifndef RELEASE
template<typename T>
void 
CheckInput
( const DistMatrix<T,MC,  Star>& A, 
  const DistMatrix<T,Star,MR  >& B,
  const DistMatrix<T,MC,  MR  >& C )
{
    if( A.GetGrid() != B.GetGrid() || B.GetGrid() != C.GetGrid() )
        throw "A, B, and C must be distributed over the same grid.";
    if( A.Height() != C.Height() || B.Width() != C.Width() ||
        A.Width() != B.Height()  || A.Height() != B.Width()   )
    {
        ostringstream msg;
        msg << "Nonconformal TriangularRankK: " << endl
            << "  A[MC,* ] ~ " << A.Height() << " x "
                               << A.Width()  << endl
            << "  B[* ,MR] ~ " << B.Height() << " x "
                               << B.Width()  << endl
            << "  C[MC,MR] ~ " << C.Height() << " x " << C.Width() << endl;
        const string s = msg.str();
        throw s.c_str();
    }
    if( A.ColAlignment() != C.ColAlignment() ||
        B.RowAlignment() != C.RowAlignment()    )
    {
        ostringstream msg;
        msg << "Misaligned TriangularRankK: " << endl
            << "  A[MC,* ] ~ " << A.ColAlignment() << endl
            << "  B[* ,MR] ~ " << B.RowAlignment() << endl
            << "  C[MC,MR] ~ " << C.ColAlignment() << " , " <<
                                  C.RowAlignment() << endl;
        const string s = msg.str();
        throw s.c_str();
    }
}

template<typename T>
void
CheckInput
( const DistMatrix<T,Star,MC  >& A,
  const DistMatrix<T,Star,MR  >& B,
  const DistMatrix<T,MC,  MR  >& C      )
{
    if( A.GetGrid() != B.GetGrid() || B.GetGrid() != C.GetGrid() )
        throw "A, B, and C must be distributed over the same grid.";
    if( A.Width() != C.Height() || B.Width() != C.Width() ||
        A.Height() != B.Height() || A.Width() != B.Width()   )
    {
        ostringstream msg;
        msg << "Nonconformal TriangularRankK: " << endl
            << "  A[* ,MC] ~ " << A.Height() << " x "
                               << A.Width()  << endl
            << "  B[* ,MR] ~ " << B.Height() << " x "
                               << B.Width()  << endl
            << "  C[MC,MR] ~ " << C.Height() << " x " << C.Width() << endl;
        const string s = msg.str();
        throw s.c_str();
    }
    if( A.RowAlignment() != C.ColAlignment() ||
        B.RowAlignment() != C.RowAlignment()    )
    {
        ostringstream msg;
        msg << "Misaligned TriangularRankK: " << endl
            << "  A[* ,MC] ~ " << A.RowAlignment() << endl
            << "  B[* ,MR] ~ " << B.RowAlignment() << endl
            << "  C[MC,MR] ~ " << C.ColAlignment() << " , " <<
                                  C.RowAlignment() << endl;
        const string s = msg.str();
        throw s.c_str();
    }
}

template<typename T>
void
CheckInput
( const DistMatrix<T,Star,MC  >& A,
  const DistMatrix<T,MR,  Star>& B,
  const DistMatrix<T,MC,  MR  >& C      )
{
    if( A.GetGrid() != B.GetGrid() || B.GetGrid() != C.GetGrid() )
        throw "A, B, and C must be distributed over the same grid.";
    if( A.Width() != C.Height() || B.Height() != C.Width() ||
        A.Height() != B.Width() || A.Width() != B.Height()   )
    {
        ostringstream msg;
        msg << "Nonconformal TriangularRankK: " << endl
            << "  A[* ,MC] ~ " << A.Height() << " x "
                               << A.Width()  << endl
            << "  B[MR,* ] ~ " << B.Height() << " x "
                               << B.Width()  << endl
            << "  C[MC,MR] ~ " << C.Height() << " x " << C.Width() << endl;
        const string s = msg.str();
        throw s.c_str();
    }
    if( A.RowAlignment() != C.ColAlignment() ||
        B.ColAlignment() != C.RowAlignment()    )
    {
        ostringstream msg;
        msg << "Misaligned TriangularRankK: " << endl
            << "  A[* ,MC] ~ " << A.RowAlignment() << endl
            << "  B[MR,* ] ~ " << B.ColAlignment() << endl
            << "  C[MC,MR] ~ " << C.ColAlignment() << " , " <<
                                  C.RowAlignment() << endl;
        const string s = msg.str();
        throw s.c_str();
    }
}
#endif

template<typename T>
void
TriangularRankKKernel
( const Shape shape,
  const T alpha, const DistMatrix<T,MC,  Star>& A,
                 const DistMatrix<T,Star,MR  >& B,
  const T beta,        DistMatrix<T,MC,  MR  >& C )
{
#ifndef RELEASE
    PushCallStack("TriangularRankKKernel");
    CheckInput( A, B, C );
#endif
    const Grid& grid = C.GetGrid();

    DistMatrix<T,MC,Star> AT(grid), 
                          AB(grid);

    DistMatrix<T,Star,MR> BL(grid), BR(grid);

    DistMatrix<T,MC,MR> CTL(grid), CTR(grid),
                        CBL(grid), CBR(grid);

    DistMatrix<T,MC,MR> DTL(grid), DBR(grid);

    const unsigned half = C.Height()/2;

    BLAS::Scal( beta, C );

    LockedPartitionDown( A, AT,
                            AB, half );

    LockedPartitionRight( B, BL, BR, half );

    PartitionDownDiagonal( C, CTL, CTR,
                              CBL, CBR, half );

    DTL.AlignWith( CTL );
    DBR.AlignWith( CBR );
    DTL.ResizeTo( CTL.Height(), CTL.Width() );
    DBR.ResizeTo( CBR.Height(), CBR.Width() );
    //------------------------------------------------------------------------//
    if( shape == Lower )
    {
        BLAS::Gemm
        ( Normal, Normal,
          alpha, AB.LockedLocalMatrix(),
                 BL.LockedLocalMatrix(),
          (T)1,  CBL.LocalMatrix()      );
    }
    else
    {
        BLAS::Gemm
        ( Normal, Normal,
          alpha, AT.LockedLocalMatrix(),
                 BR.LockedLocalMatrix(),
          (T)1,  CTR.LocalMatrix()      );
    }

    BLAS::Gemm
    ( Normal, Normal,
      alpha, AT.LockedLocalMatrix(),
             BL.LockedLocalMatrix(),
      (T)0,  DTL.LocalMatrix()      );
    DTL.MakeTrapezoidal( Left, shape );
    BLAS::Axpy( (T)1, DTL, CTL );

    BLAS::Gemm
    ( Normal, Normal,
      alpha, AB.LockedLocalMatrix(),
             BR.LockedLocalMatrix(),
      (T)0,  DBR.LocalMatrix()      );
    DBR.MakeTrapezoidal( Left, shape );
    BLAS::Axpy( (T)1, DBR, CBR );
    //------------------------------------------------------------------------//
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
TriangularRankKKernel
( const Shape shape,
  const T alpha, const DistMatrix<T,Star,MC  >& A,
                 const DistMatrix<T,Star,MR  >& B,
  const T beta,        DistMatrix<T,MC,  MR  >& C )
{
#ifndef RELEASE
    PushCallStack("TriangularRankKKernel");
    CheckInput( A, B, C );
#endif
    const Grid& grid = C.GetGrid();

    DistMatrix<T,Star,MC> AL(grid), AR(grid);

    DistMatrix<T,Star,MR> BL(grid), BR(grid);

    DistMatrix<T,MC,MR> CTL(grid), CTR(grid),
                        CBL(grid), CBR(grid);

    DistMatrix<T,MC,MR> DTL(grid), DBR(grid);

    const unsigned half = C.Height()/2;

    BLAS::Scal( beta, C );

    LockedPartitionRight( A, AL, AR, half );

    LockedPartitionRight( B, BL, BR, half );

    PartitionDownDiagonal( C, CTL, CTR,
                              CBL, CBR, half );

    DTL.AlignWith( CTL );
    DBR.AlignWith( CBR );
    DTL.ResizeTo( CTL.Height(), CTL.Width() );
    DBR.ResizeTo( CBR.Height(), CBR.Width() );
    //------------------------------------------------------------------------//
    if( shape == Lower )
    {
        BLAS::Gemm
        ( Transpose, Normal,
          alpha, AR.LockedLocalMatrix(),
                 BL.LockedLocalMatrix(),
          (T)1,  CBL.LocalMatrix()      );
    }
    else
    {
        BLAS::Gemm
        ( Transpose, Normal,
          alpha, AL.LockedLocalMatrix(),
                 BR.LockedLocalMatrix(),
          (T)1,  CTR.LocalMatrix()      );
    }

    BLAS::Gemm
    ( Transpose, Normal,
      alpha, AL.LockedLocalMatrix(),
             BL.LockedLocalMatrix(),
      (T)0,  DTL.LocalMatrix()      );
    DTL.MakeTrapezoidal( Left, shape );
    BLAS::Axpy( (T)1, DTL, CTL );

    BLAS::Gemm
    ( Transpose, Normal,
      alpha, AR.LockedLocalMatrix(),
             BR.LockedLocalMatrix(),
      (T)0,  DBR.LocalMatrix()      );
    DBR.MakeTrapezoidal( Left, shape );
    BLAS::Axpy( (T)1, DBR, CBR );
    //------------------------------------------------------------------------//
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
TriangularRankKKernel
( const Shape shape,
  const T alpha, const DistMatrix<T,Star,MC  >& A,
                 const DistMatrix<T,MR,  Star>& B,
  const T beta,        DistMatrix<T,MC,  MR  >& C )
{
#ifndef RELEASE
    PushCallStack("TriangularRankKKernel");
    CheckInput( A, B, C );
#endif
    const Grid& grid = C.GetGrid();

    DistMatrix<T,Star,MC> AL(grid), AR(grid);

    DistMatrix<T,MR,Star> BT(grid), 
                          BB(grid);

    DistMatrix<T,MC,MR> CTL(grid), CTR(grid),
                        CBL(grid), CBR(grid);

    DistMatrix<T,MC,MR> DTL(grid), DBR(grid);

    const unsigned half = C.Height()/2;

    BLAS::Scal( beta, C );

    LockedPartitionRight( A, AL, AR, half );

    LockedPartitionDown( B, BT, 
                            BB, half );

    PartitionDownDiagonal( C, CTL, CTR,
                              CBL, CBR, half );

    DTL.AlignWith( CTL );
    DBR.AlignWith( CBR );
    DTL.ResizeTo( CTL.Height(), CTL.Width() );
    DBR.ResizeTo( CBR.Height(), CBR.Width() );
    //------------------------------------------------------------------------//
    if( shape == Lower )
    {
        BLAS::Gemm
        ( Transpose, Transpose,
          alpha, AR.LockedLocalMatrix(),
                 BT.LockedLocalMatrix(),
          (T)1,  CBL.LocalMatrix()      );
    }
    else
    {
        BLAS::Gemm
        ( Transpose, Transpose,
          alpha, AL.LockedLocalMatrix(),
                 BB.LockedLocalMatrix(),
          (T)1,  CTR.LocalMatrix()      );
    }

    BLAS::Gemm
    ( Transpose, Transpose,
      alpha, AL.LockedLocalMatrix(),
             BT.LockedLocalMatrix(),
      (T)0,  DTL.LocalMatrix()      );
    DTL.MakeTrapezoidal( Left, shape );
    BLAS::Axpy( (T)1, DTL, CTL );

    BLAS::Gemm
    ( Transpose, Transpose,
      alpha, AR.LockedLocalMatrix(),
             BB.LockedLocalMatrix(),
      (T)0,  DBR.LocalMatrix()      );
    DBR.MakeTrapezoidal( Left, shape );
    BLAS::Axpy( (T)1, DBR, CBR );
    //------------------------------------------------------------------------//
#ifndef RELEASE
    PopCallStack();
#endif
}

} // anonymous namespace

template<typename T>
void
Elemental::BLAS::Internal::TriangularRankK
( const Shape shape,
  const T alpha, const DistMatrix<T,MC,  Star>& A,
                 const DistMatrix<T,Star,MR  >& B,
  const T beta,        DistMatrix<T,MC,  MR  >& C )
{
#ifndef RELEASE
    PushCallStack("BLAS::Internal::TriangularRankK");
    CheckInput( A, B, C );
#endif
    const Grid& grid = C.GetGrid();

    if( C.Height() < 2*grid.Width()*Blocksize() )
    {
        TriangularRankKKernel
        ( shape, alpha, A, B, beta, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.

        DistMatrix<T,MC,Star> AT(grid),
                              AB(grid);
        DistMatrix<T,Star,MR> BL(grid), BR(grid);
        DistMatrix<T,MC,MR> CTL(grid), CTR(grid),
                            CBL(grid), CBR(grid);

        const unsigned half = C.Height() / 2;

        LockedPartitionDown( A, AT,
                                AB, half );

        LockedPartitionRight( B, BL, BR, half );

        PartitionDownDiagonal( C, CTL, CTR,
                                  CBL, CBR, half );

        if( shape == Lower )
        { 
            BLAS::Gemm
            ( Normal, Normal, 
              alpha, AB.LockedLocalMatrix(),
                     BL.LockedLocalMatrix(),
              beta,  CBL.LocalMatrix()      );
        }
        else
        {
            BLAS::Gemm
            ( Normal, Normal,
              alpha, AT.LockedLocalMatrix(),
                     BR.LockedLocalMatrix(),
              beta,  CTR.LocalMatrix()      );
        }

        // Recurse
        BLAS::Internal::TriangularRankK
        ( shape, alpha, AT, BL, beta, CTL );

        BLAS::Internal::TriangularRankK
        ( shape, alpha, AB, BR, beta, CBR );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::BLAS::Internal::TriangularRankK
( const Shape shape,
  const T alpha, const DistMatrix<T,Star,MC>& A,
                 const DistMatrix<T,Star,MR>& B,
  const T beta,        DistMatrix<T,MC,  MR>& C )
{
#ifndef RELEASE
    PushCallStack("BLAS::Internal::TriangularRankK");
    CheckInput( A, B, C );
#endif
    const Grid& grid = C.GetGrid();

    if( C.Height() < 2*grid.Width()*Blocksize() )
    {
        TriangularRankKKernel
        ( shape, alpha, A, B, beta, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.

        DistMatrix<T,Star,MC> AL(grid), AR(grid);
        DistMatrix<T,Star,MR> BL(grid), BR(grid);
        DistMatrix<T,MC,MR> CTL(grid), CTR(grid),
                            CBL(grid), CBR(grid);

        const unsigned half = C.Height() / 2;

        LockedPartitionRight( A, AL, AR, half );

        LockedPartitionRight( B, BL, BR, half );

        PartitionDownDiagonal( C, CTL, CTR,
                                  CBL, CBR, half );

        if( shape == Lower )
        { 
            BLAS::Gemm
            ( Transpose, Normal, 
              alpha, AR.LockedLocalMatrix(),
                     BL.LockedLocalMatrix(),
              beta,  CBL.LocalMatrix()      );
        }
        else
        {
            BLAS::Gemm
            ( Transpose, Normal,
              alpha, AL.LockedLocalMatrix(),
                     BR.LockedLocalMatrix(),
              beta,  CTR.LocalMatrix()      );
        }

        // Recurse
        BLAS::Internal::TriangularRankK
        ( shape, alpha, AL, BL, beta, CTL );

        BLAS::Internal::TriangularRankK
        ( shape, alpha, AR, BR, beta, CBR );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::BLAS::Internal::TriangularRankK
( const Shape shape,
  const T alpha, const DistMatrix<T,Star,MC  >& A,
                 const DistMatrix<T,MR,  Star>& B,
  const T beta,        DistMatrix<T,MC,  MR  >& C )
{
#ifndef RELEASE
    PushCallStack("BLAS::Internal::TriangularRankK");
    CheckInput( A, B, C );
#endif
    const Grid& grid = C.GetGrid();

    if( C.Height() < 2*grid.Width()*Blocksize() )
    {
        TriangularRankKKernel
        ( shape, alpha, A, B, beta, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.

        DistMatrix<T,Star,MC> AL(grid), AR(grid);
        DistMatrix<T,MR,Star> BT(grid), 
                              BB(grid);
        DistMatrix<T,MC,MR> CTL(grid), CTR(grid),
                            CBL(grid), CBR(grid);

        const unsigned half = C.Height() / 2;

        LockedPartitionRight( A, AL, AR, half );

        LockedPartitionDown( B, BT, 
                                BB, half );

        PartitionDownDiagonal( C, CTL, CTR,
                                  CBL, CBR, half );

        if( shape == Lower )
        { 
            BLAS::Gemm
            ( Transpose, Transpose, 
              alpha, AR.LockedLocalMatrix(),
                     BT.LockedLocalMatrix(),
              beta,  CBL.LocalMatrix()      );
        }
        else
        {
            BLAS::Gemm
            ( Transpose, Transpose,
              alpha, AL.LockedLocalMatrix(),
                     BB.LockedLocalMatrix(),
              beta,  CTR.LocalMatrix()      );
        }

        // Recurse
        BLAS::Internal::TriangularRankK
        ( shape, alpha, AL, BT, beta, CTL );

        BLAS::Internal::TriangularRankK
        ( shape, alpha, AR, BB, beta, CBR );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void
Elemental::BLAS::Internal::TriangularRankK
( const Shape shape, 
  const float alpha, const DistMatrix<float,MC,  Star>& A,
                     const DistMatrix<float,Star,MR  >& B,
  const float beta,        DistMatrix<float,MC,  MR  >& C );

template void
Elemental::BLAS::Internal::TriangularRankK
( const Shape shape, 
  const float alpha, const DistMatrix<float,Star,MC>& A,
                     const DistMatrix<float,Star,MR>& B,
  const float beta,        DistMatrix<float,MC,  MR>& C );

template void
Elemental::BLAS::Internal::TriangularRankK
( const Shape shape, 
  const float alpha, const DistMatrix<float,Star,MC  >& A,
                     const DistMatrix<float,MR,  Star>& B,
  const float beta,        DistMatrix<float,MC,  MR  >& C );

template void
Elemental::BLAS::Internal::TriangularRankK
( const Shape shape,
  const double alpha, const DistMatrix<double,MC,  Star>& A,
                      const DistMatrix<double,Star,MR  >& B,
  const double beta,        DistMatrix<double,MC,  MR  >& C );

template void
Elemental::BLAS::Internal::TriangularRankK
( const Shape shape,
  const double alpha, const DistMatrix<double,Star,MC>& A,
                      const DistMatrix<double,Star,MR>& B,
  const double beta,        DistMatrix<double,MC,  MR>& C );

template void
Elemental::BLAS::Internal::TriangularRankK
( const Shape shape,
  const double alpha, const DistMatrix<double,Star,MC  >& A,
                      const DistMatrix<double,MR,  Star>& B,
  const double beta,        DistMatrix<double,MC,  MR  >& C );

#ifndef WITHOUT_COMPLEX
template void
Elemental::BLAS::Internal::TriangularRankK
( const Shape shape,
  const scomplex alpha, const DistMatrix<scomplex,MC,  Star>& A,
                        const DistMatrix<scomplex,Star,MR  >& B,
  const scomplex beta,        DistMatrix<scomplex,MC,  MR  >& C );

template void
Elemental::BLAS::Internal::TriangularRankK
( const Shape shape,
  const scomplex alpha, const DistMatrix<scomplex,Star,MC>& A,
                        const DistMatrix<scomplex,Star,MR>& B,
  const scomplex beta,        DistMatrix<scomplex,MC,  MR>& C );

template void
Elemental::BLAS::Internal::TriangularRankK
( const Shape shape,
  const scomplex alpha, const DistMatrix<scomplex,Star,MC  >& A,
                        const DistMatrix<scomplex,MR,  Star>& B,
  const scomplex beta,        DistMatrix<scomplex,MC,  MR  >& C );

template void
Elemental::BLAS::Internal::TriangularRankK
( const Shape shape,
  const dcomplex alpha, const DistMatrix<dcomplex,MC,  Star>& A,
                        const DistMatrix<dcomplex,Star,MR  >& B,
  const dcomplex beta,        DistMatrix<dcomplex,MC,  MR  >& C );

template void
Elemental::BLAS::Internal::TriangularRankK
( const Shape shape,
  const dcomplex alpha, const DistMatrix<dcomplex,Star,MC>& A,
                        const DistMatrix<dcomplex,Star,MR>& B,
  const dcomplex beta,        DistMatrix<dcomplex,MC,  MR>& C );

template void
Elemental::BLAS::Internal::TriangularRankK
( const Shape shape,
  const dcomplex alpha, const DistMatrix<dcomplex,Star,MC  >& A,
                        const DistMatrix<dcomplex,MR,  Star>& B,
  const dcomplex beta,        DistMatrix<dcomplex,MC,  MR  >& C );
#endif

