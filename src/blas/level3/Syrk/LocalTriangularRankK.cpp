/*
   Copyright (c) 2009-2010, Jack Poulson
   All rights reserved.

   This file is part of Elemental.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

    - Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

    - Neither the name of the owner nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/
#include "elemental/blas_internal.hpp"
using namespace std;
using namespace elemental;

namespace {

#ifndef RELEASE
template<typename T>
void
CheckInput
( Orientation orientationOfB,
  const DistMatrix<T,MC,Star>& A,
  const DistMatrix<T,MR,Star>& B,
  const DistMatrix<T,MC,MR  >& C )
{
    if( orientationOfB == Normal )
        throw logic_error( "B[MR,* ] must be (Conjugate)Transpose'd." );
    if( A.GetGrid() != B.GetGrid() || B.GetGrid() != C.GetGrid() )
        throw logic_error
        ( "A, B, and C must be distributed over the same grid." );
    if( A.Height() != C.Height() || B.Height() != C.Width() ||
        A.Width() != B.Width() || A.Height() != B.Height() )
    {
        ostringstream msg;
        msg << "Nonconformal LocalTriangularRankK: " << endl
            << "  A[MC,* ] ~ " << A.Height() << " x "
                               << A.Width()  << endl
            << "  B[MR,* ] ~ " << B.Height() << " x "
                               << B.Width()  << endl
            << "  C[MC,MR] ~ " << C.Height() << " x " << C.Width() << endl;
        throw logic_error( msg.str() );
    }
    if( A.ColAlignment() != C.ColAlignment() ||
        B.ColAlignment() != C.RowAlignment() )
    {
        ostringstream msg;
        msg << "Misaligned LocalTriangularRankK: " << endl
            << "  A[MC,* ] ~ " << A.ColAlignment() << endl
            << "  B[MR,* ] ~ " << B.ColAlignment() << endl
            << "  C[MC,MR] ~ " << C.ColAlignment() << " , " <<
                                  C.RowAlignment() << endl;
        throw logic_error( msg.str() );
    }
}

template<typename T>
void
CheckInput
( Orientation orientationOfA,
  Orientation orientationOfB,
  const DistMatrix<T,Star,MC  >& A,
  const DistMatrix<T,MR,  Star>& B,
  const DistMatrix<T,MC,  MR  >& C )
{
    if( orientationOfA == Normal )
        throw logic_error( "A[* ,MC] must be (Conjugate)Transpose'd." );
    if( orientationOfB == Normal )
        throw logic_error( "B[MR,* ] must be (Conjugate)Transpose'd." );
    if( A.GetGrid() != B.GetGrid() || B.GetGrid() != C.GetGrid() )
        throw logic_error
        ( "A, B, and C must be distributed over the same grid." );
    if( A.Width() != C.Height() || B.Height() != C.Width() ||
        A.Height() != B.Width() || A.Width() != B.Height() )
    {
        ostringstream msg;
        msg << "Nonconformal LocalTriangularRankK: " << endl
            << "  A[* ,MC] ~ " << A.Height() << " x "
                               << A.Width()  << endl
            << "  B[MR,* ] ~ " << B.Height() << " x "
                               << B.Width()  << endl
            << "  C[MC,MR] ~ " << C.Height() << " x " << C.Width() << endl;
        throw logic_error( msg.str() );
    }
    if( A.RowAlignment() != C.ColAlignment() ||
        B.ColAlignment() != C.RowAlignment() )
    {
        ostringstream msg;
        msg << "Misaligned LocalTriangularRankK: " << endl
            << "  A[* ,MC] ~ " << A.RowAlignment() << endl
            << "  B[MR,* ] ~ " << B.ColAlignment() << endl
            << "  C[MC,MR] ~ " << C.ColAlignment() << " , " <<
                                  C.RowAlignment() << endl;
        throw logic_error( msg.str() );
    }
}

template<typename T>
void 
CheckInput
( const DistMatrix<T,MC,  Star>& A, 
  const DistMatrix<T,Star,MR  >& B,
  const DistMatrix<T,MC,  MR  >& C )
{
    if( A.GetGrid() != B.GetGrid() || B.GetGrid() != C.GetGrid() )
        throw logic_error
        ( "A, B, and C must be distributed over the same grid." );
    if( A.Height() != C.Height() || B.Width() != C.Width() ||
        A.Width() != B.Height()  || A.Height() != B.Width() )
    {
        ostringstream msg;
        msg << "Nonconformal LocalTriangularRankK: " << endl
            << "  A[MC,* ] ~ " << A.Height() << " x "
                               << A.Width()  << endl
            << "  B[* ,MR] ~ " << B.Height() << " x "
                               << B.Width()  << endl
            << "  C[MC,MR] ~ " << C.Height() << " x " << C.Width() << endl;
        throw logic_error( msg.str() );
    }
    if( A.ColAlignment() != C.ColAlignment() ||
        B.RowAlignment() != C.RowAlignment() )
    {
        ostringstream msg;
        msg << "Misaligned LocalTriangularRankK: " << endl
            << "  A[MC,* ] ~ " << A.ColAlignment() << endl
            << "  B[* ,MR] ~ " << B.RowAlignment() << endl
            << "  C[MC,MR] ~ " << C.ColAlignment() << " , " <<
                                  C.RowAlignment() << endl;
        throw logic_error( msg.str() );
    }
}

template<typename T>
void
CheckInput
( Orientation orientationOfA,
  const DistMatrix<T,Star,MC>& A,
  const DistMatrix<T,Star,MR>& B,
  const DistMatrix<T,MC,  MR>& C )
{
    if( orientationOfA == Normal )
        throw logic_error( "A[* ,MC] must be (Conjugate)Transpose'd." );
    if( A.GetGrid() != B.GetGrid() || B.GetGrid() != C.GetGrid() )
        throw logic_error
        ( "A, B, and C must be distributed over the same grid." );
    if( A.Width() != C.Height() || B.Width() != C.Width() ||
        A.Height() != B.Height() || A.Width() != B.Width() )
    {
        ostringstream msg;
        msg << "Nonconformal LocalTriangularRankK: " << endl
            << "  A[* ,MC] ~ " << A.Height() << " x "
                               << A.Width()  << endl
            << "  B[* ,MR] ~ " << B.Height() << " x "
                               << B.Width()  << endl
            << "  C[MC,MR] ~ " << C.Height() << " x " << C.Width() << endl;
        throw logic_error( msg.str() );
    }
    if( A.RowAlignment() != C.ColAlignment() ||
        B.RowAlignment() != C.RowAlignment() )
    {
        ostringstream msg;
        msg << "Misaligned LocalTriangularRankK: " << endl
            << "  A[* ,MC] ~ " << A.RowAlignment() << endl
            << "  B[* ,MR] ~ " << B.RowAlignment() << endl
            << "  C[MC,MR] ~ " << C.ColAlignment() << " , " <<
                                  C.RowAlignment() << endl;
        throw logic_error( msg.str() );
    }
}
#endif

template<typename T>
void
LocalTriangularRankKKernel
( Shape shape,
  Orientation orientationOfB,
  T alpha, const DistMatrix<T,MC,Star>& A,
           const DistMatrix<T,MR,Star>& B,
  T beta,        DistMatrix<T,MC,MR  >& C )
{
#ifndef RELEASE
    PushCallStack("LocalTriangularRankKKernel");
    CheckInput( orientationOfB, A, B, C );
#endif
    const Grid& g = C.GetGrid();

    DistMatrix<T,MC,Star> AT(g),
                          AB(g);

    DistMatrix<T,MR,Star> BT(g), 
                          BB(g);

    DistMatrix<T,MC,MR> CTL(g), CTR(g),
                        CBL(g), CBR(g);

    DistMatrix<T,MC,MR> DTL(g), DBR(g);

    const unsigned half = C.Height()/2;

    blas::Scal( beta, C );

    LockedPartitionDown
    ( A, AT,
         AB, half );

    LockedPartitionDown
    ( B, BT, 
         BB, half );

    PartitionDownDiagonal
    ( C, CTL, CTR,
         CBL, CBR, half );

    DTL.AlignWith( CTL );
    DBR.AlignWith( CBR );
    DTL.ResizeTo( CTL.Height(), CTL.Width() );
    DBR.ResizeTo( CBR.Height(), CBR.Width() );
    //------------------------------------------------------------------------//
    if( shape == Lower )
    {
        blas::internal::LocalGemm
        ( Normal, orientationOfB, alpha, AB, BT, (T)1, CBL );
    }
    else
    {
        blas::internal::LocalGemm
        ( Normal, orientationOfB, alpha, AT, BB, (T)1, CTR );
    }

    blas::internal::LocalGemm
    ( Normal, orientationOfB, alpha, AT, BT, (T)0, DTL );

    DTL.MakeTrapezoidal( Left, shape );
    blas::Axpy( (T)1, DTL, CTL );

    blas::internal::LocalGemm
    ( Normal, orientationOfB, alpha, AB, BB, (T)0, DBR );

    DBR.MakeTrapezoidal( Left, shape );
    blas::Axpy( (T)1, DBR, CBR );
    //------------------------------------------------------------------------//
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
LocalTriangularRankKKernel
( Shape shape,
  Orientation orientationOfA,
  Orientation orientationOfB,
  T alpha, const DistMatrix<T,Star,MC  >& A,
           const DistMatrix<T,MR,  Star>& B,
  T beta,        DistMatrix<T,MC,  MR  >& C )
{
#ifndef RELEASE
    PushCallStack("LocalTriangularRankKKernel");
    CheckInput( orientationOfA, orientationOfB, A, B, C );
#endif
    const Grid& g = C.GetGrid();

    DistMatrix<T,Star,MC> AL(g), AR(g);

    DistMatrix<T,MR,Star> BT(g), 
                          BB(g);

    DistMatrix<T,MC,MR> CTL(g), CTR(g),
                        CBL(g), CBR(g);

    DistMatrix<T,MC,MR> DTL(g), DBR(g);

    const unsigned half = C.Height()/2;

    blas::Scal( beta, C );

    LockedPartitionRight( A, AL, AR, half );

    LockedPartitionDown
    ( B, BT, 
         BB, half );

    PartitionDownDiagonal
    ( C, CTL, CTR,
         CBL, CBR, half );

    DTL.AlignWith( CTL );
    DBR.AlignWith( CBR );
    DTL.ResizeTo( CTL.Height(), CTL.Width() );
    DBR.ResizeTo( CBR.Height(), CBR.Width() );
    //------------------------------------------------------------------------//
    if( shape == Lower )
    {
        blas::internal::LocalGemm
        ( orientationOfA, orientationOfB, alpha, AR, BT, (T)1, CBL );
    }
    else
    {
        blas::internal::LocalGemm
        ( orientationOfA, orientationOfB, alpha, AL, BB, (T)1, CTR );
    }

    blas::internal::LocalGemm
    ( orientationOfA, orientationOfB, alpha, AL, BT, (T)0, DTL );

    DTL.MakeTrapezoidal( Left, shape );
    blas::Axpy( (T)1, DTL, CTL );

    blas::internal::LocalGemm
    ( orientationOfA, orientationOfB, alpha, AR, BB, (T)0, DBR );

    DBR.MakeTrapezoidal( Left, shape );
    blas::Axpy( (T)1, DBR, CBR );
    //------------------------------------------------------------------------//
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
LocalTriangularRankKKernel
( Shape shape, 
  T alpha, const DistMatrix<T,MC,  Star>& A,
           const DistMatrix<T,Star,MR  >& B,
  T beta,        DistMatrix<T,MC,  MR  >& C )
{
#ifndef RELEASE
    PushCallStack("LocalTriangularRankKKernel");
    CheckInput( A, B, C );
#endif
    const Grid& g = C.GetGrid();

    DistMatrix<T,MC,Star> AT(g), 
                          AB(g);

    DistMatrix<T,Star,MR> BL(g), BR(g);

    DistMatrix<T,MC,MR> CTL(g), CTR(g),
                        CBL(g), CBR(g);

    DistMatrix<T,MC,MR> DTL(g), DBR(g);

    const unsigned half = C.Height()/2;

    blas::Scal( beta, C );

    LockedPartitionDown
    ( A, AT,
         AB, half );

    LockedPartitionRight( B, BL, BR, half );

    PartitionDownDiagonal
    ( C, CTL, CTR,
         CBL, CBR, half );

    DTL.AlignWith( CTL );
    DBR.AlignWith( CBR );
    DTL.ResizeTo( CTL.Height(), CTL.Width() );
    DBR.ResizeTo( CBR.Height(), CBR.Width() );
    //------------------------------------------------------------------------//
    if( shape == Lower )
    {
        blas::internal::LocalGemm
        ( Normal, Normal, alpha, AB, BL, (T)1, CBL );
    }
    else
    {
        blas::internal::LocalGemm
        ( Normal, Normal, alpha, AT, BR, (T)1, CTR );
    }

    blas::internal::LocalGemm
    ( Normal, Normal, alpha, AT, BL, (T)0, DTL );

    DTL.MakeTrapezoidal( Left, shape );
    blas::Axpy( (T)1, DTL, CTL );

    blas::internal::LocalGemm
    ( Normal, Normal, alpha, AB, BR, (T)0, DBR );

    DBR.MakeTrapezoidal( Left, shape );
    blas::Axpy( (T)1, DBR, CBR );
    //------------------------------------------------------------------------//
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
LocalTriangularRankKKernel
( Shape shape,
  Orientation orientationOfA,
  T alpha, const DistMatrix<T,Star,MC>& A,
           const DistMatrix<T,Star,MR>& B,
  T beta,        DistMatrix<T,MC,  MR>& C )
{
#ifndef RELEASE
    PushCallStack("LocalTriangularRankKKernel");
    CheckInput( orientationOfA, A, B, C );
#endif
    const Grid& g = C.GetGrid();

    DistMatrix<T,Star,MC> AL(g), AR(g);

    DistMatrix<T,Star,MR> BL(g), BR(g);

    DistMatrix<T,MC,MR> CTL(g), CTR(g),
                        CBL(g), CBR(g);

    DistMatrix<T,MC,MR> DTL(g), DBR(g);

    const unsigned half = C.Height()/2;

    blas::Scal( beta, C );

    LockedPartitionRight( A, AL, AR, half );
    LockedPartitionRight( B, BL, BR, half );

    PartitionDownDiagonal
    ( C, CTL, CTR,
         CBL, CBR, half );

    DTL.AlignWith( CTL );
    DBR.AlignWith( CBR );
    DTL.ResizeTo( CTL.Height(), CTL.Width() );
    DBR.ResizeTo( CBR.Height(), CBR.Width() );
    //------------------------------------------------------------------------//
    if( shape == Lower )
    {
        blas::internal::LocalGemm
        ( orientationOfA, Normal, alpha, AR, BL, (T)1, CBL );
    }
    else
    {
        blas::internal::LocalGemm
        ( orientationOfA, Normal, alpha, AL, BR, (T)1, CTR );
    }

    blas::internal::LocalGemm
    ( orientationOfA, Normal, alpha, AL, BL, (T)0, DTL );

    DTL.MakeTrapezoidal( Left, shape );
    blas::Axpy( (T)1, DTL, CTL );

    blas::internal::LocalGemm
    ( orientationOfA, Normal, alpha, AR, BR, (T)0, DBR );

    DBR.MakeTrapezoidal( Left, shape );
    blas::Axpy( (T)1, DBR, CBR );
    //------------------------------------------------------------------------//
#ifndef RELEASE
    PopCallStack();
#endif
}

} // anonymous namespace

template<typename T>
void
elemental::blas::internal::LocalTriangularRankK
( Shape shape,
  Orientation orientationOfB,
  T alpha, const DistMatrix<T,MC,Star>& A,
           const DistMatrix<T,MR,Star>& B,
  T beta,        DistMatrix<T,MC,MR  >& C )
{
#ifndef RELEASE
    PushCallStack("blas::internal::LocalTriangularRankK");
    CheckInput( orientationOfB, A, B, C );
#endif
    const Grid& g = C.GetGrid();

    if( C.Height() < 2*g.Width()*Blocksize() )
    {
        LocalTriangularRankKKernel
        ( shape, orientationOfB, alpha, A, B, beta, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.

        DistMatrix<T,MC,Star> AT(g),
                              AB(g);

        DistMatrix<T,MR,Star> BT(g), 
                              BB(g);

        DistMatrix<T,MC,MR> CTL(g), CTR(g),
                            CBL(g), CBR(g);

        const unsigned half = C.Height() / 2;

        LockedPartitionDown
        ( A, AT,
             AB, half );

        LockedPartitionDown
        ( B, BT, 
             BB, half );

        PartitionDownDiagonal
        ( C, CTL, CTR,
             CBL, CBR, half );

        if( shape == Lower )
        { 
            blas::internal::LocalGemm
            ( Normal, orientationOfB, alpha, AB, BT, beta, CBL );
        }
        else
        {
            blas::internal::LocalGemm
            ( Normal, orientationOfB, alpha, AT, BB, beta, CTR );
        }

        // Recurse
        blas::internal::LocalTriangularRankK
        ( shape, orientationOfB, alpha, AT, BT, beta, CTL );

        blas::internal::LocalTriangularRankK
        ( shape, orientationOfB, alpha, AB, BB, beta, CBR );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::blas::internal::LocalTriangularRankK
( Shape shape,
  Orientation orientationOfA,
  Orientation orientationOfB,
  T alpha, const DistMatrix<T,Star,MC  >& A,
           const DistMatrix<T,MR,  Star>& B,
  T beta,        DistMatrix<T,MC,  MR  >& C )
{
#ifndef RELEASE
    PushCallStack("blas::internal::LocalTriangularRankK");
    CheckInput( orientationOfA, orientationOfB, A, B, C );
#endif
    const Grid& g = C.GetGrid();

    if( C.Height() < 2*g.Width()*Blocksize() )
    {
        LocalTriangularRankKKernel
        ( shape, orientationOfA, orientationOfB, alpha, A, B, beta, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.

        DistMatrix<T,Star,MC> AL(g), AR(g);

        DistMatrix<T,MR,Star> BT(g), 
                              BB(g);

        DistMatrix<T,MC,MR> CTL(g), CTR(g),
                            CBL(g), CBR(g);

        const unsigned half = C.Height() / 2;

        LockedPartitionRight( A, AL, AR, half );

        LockedPartitionDown
        ( B, BT, 
             BB, half );

        PartitionDownDiagonal
        ( C, CTL, CTR,
             CBL, CBR, half );

        if( shape == Lower )
        { 
            blas::internal::LocalGemm
            ( orientationOfA, orientationOfB, alpha, AR, BT, beta, CBL );
        }
        else
        {
            blas::internal::LocalGemm
            ( orientationOfA, orientationOfB, alpha, AL, BB, beta, CTR );
        }

        // Recurse
        blas::internal::LocalTriangularRankK
        ( shape, orientationOfA, orientationOfB, alpha, AL, BT, beta, CTL );

        blas::internal::LocalTriangularRankK
        ( shape, orientationOfA, orientationOfB, alpha, AR, BB, beta, CBR );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::blas::internal::LocalTriangularRankK
( Shape shape,
  T alpha, const DistMatrix<T,MC,  Star>& A,
           const DistMatrix<T,Star,MR  >& B,
  T beta,        DistMatrix<T,MC,  MR  >& C )
{
#ifndef RELEASE
    PushCallStack("blas::internal::LocalTriangularRankK");
    CheckInput( A, B, C );
#endif
    const Grid& g = C.GetGrid();

    if( C.Height() < 2*g.Width()*Blocksize() )
    {
        LocalTriangularRankKKernel
        ( shape, alpha, A, B, beta, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.

        DistMatrix<T,MC,Star> AT(g),
                              AB(g);

        DistMatrix<T,Star,MR> BL(g), BR(g);

        DistMatrix<T,MC,MR> CTL(g), CTR(g),
                            CBL(g), CBR(g);

        const unsigned half = C.Height() / 2;

        LockedPartitionDown
        ( A, AT,
             AB, half );

        LockedPartitionRight( B, BL, BR, half );

        PartitionDownDiagonal
        ( C, CTL, CTR,
             CBL, CBR, half );

        if( shape == Lower )
        { 
            blas::internal::LocalGemm
            ( Normal, Normal, alpha, AB, BL, beta, CBL );
        }
        else
        {
            blas::internal::LocalGemm
            ( Normal, Normal, alpha, AT, BR, beta, CTR );
        }

        // Recurse
        blas::internal::LocalTriangularRankK
        ( shape, alpha, AT, BL, beta, CTL );

        blas::internal::LocalTriangularRankK
        ( shape, alpha, AB, BR, beta, CBR );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::blas::internal::LocalTriangularRankK
( Shape shape,
  Orientation orientationOfA,
  T alpha, const DistMatrix<T,Star,MC>& A,
           const DistMatrix<T,Star,MR>& B,
  T beta,        DistMatrix<T,MC,  MR>& C )
{
#ifndef RELEASE
    PushCallStack("blas::internal::LocalTriangularRankK");
    CheckInput( orientationOfA, A, B, C );
#endif
    const Grid& g = C.GetGrid();

    if( C.Height() < 2*g.Width()*Blocksize() )
    {
        LocalTriangularRankKKernel
        ( shape, orientationOfA, alpha, A, B, beta, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.

        DistMatrix<T,Star,MC> AL(g), AR(g);

        DistMatrix<T,Star,MR> BL(g), BR(g);

        DistMatrix<T,MC,MR> CTL(g), CTR(g),
                            CBL(g), CBR(g);

        const unsigned half = C.Height() / 2;

        LockedPartitionRight( A, AL, AR, half );

        LockedPartitionRight( B, BL, BR, half );

        PartitionDownDiagonal
        ( C, CTL, CTR,
             CBL, CBR, half );

        if( shape == Lower )
        { 
            blas::internal::LocalGemm
            ( orientationOfA, Normal, alpha, AR, BL, beta, CBL );
        }
        else
        {
            blas::internal::LocalGemm
            ( orientationOfA, Normal, alpha, AL, BR, beta, CTR );
        }

        // Recurse
        blas::internal::LocalTriangularRankK
        ( shape, orientationOfA, alpha, AL, BL, beta, CTL );

        blas::internal::LocalTriangularRankK
        ( shape, orientationOfA, alpha, AR, BR, beta, CBR );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void
elemental::blas::internal::LocalTriangularRankK
( Shape shape,
  Orientation orientationOfB,
  float alpha, const DistMatrix<float,MC,  Star>& A,
               const DistMatrix<float,MR,  Star>& B,
  float beta,        DistMatrix<float,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRankK
( Shape shape, 
  Orientation orientationOfA,
  Orientation orientationOfB,
  float alpha, const DistMatrix<float,Star,MC  >& A,
               const DistMatrix<float,MR,  Star>& B,
  float beta,        DistMatrix<float,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRankK
( Shape shape, 
  float alpha, const DistMatrix<float,MC,  Star>& A,
               const DistMatrix<float,Star,MR  >& B,
  float beta,        DistMatrix<float,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRankK
( Shape shape, 
  Orientation orientationOfA,
  float alpha, const DistMatrix<float,Star,MC>& A,
               const DistMatrix<float,Star,MR>& B,
  float beta,        DistMatrix<float,MC,  MR>& C );

template void
elemental::blas::internal::LocalTriangularRankK
( Shape shape,
  Orientation orientationOfB,
  double alpha, const DistMatrix<double,MC,  Star>& A,
                const DistMatrix<double,MR,  Star>& B,
  double beta,        DistMatrix<double,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRankK
( Shape shape,
  Orientation orientationOfA,
  Orientation orientationOfB,
  double alpha, const DistMatrix<double,Star,MC  >& A,
                const DistMatrix<double,MR,  Star>& B,
  double beta,        DistMatrix<double,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRankK
( Shape shape,
  double alpha, const DistMatrix<double,MC,  Star>& A,
                const DistMatrix<double,Star,MR  >& B,
  double beta,        DistMatrix<double,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRankK
( Shape shape,
  Orientation orientationOfA,
  double alpha, const DistMatrix<double,Star,MC>& A,
                const DistMatrix<double,Star,MR>& B,
  double beta,        DistMatrix<double,MC,  MR>& C );


#ifndef WITHOUT_COMPLEX
template void
elemental::blas::internal::LocalTriangularRankK
( Shape shape,
  Orientation orientationOfB,
  scomplex alpha, const DistMatrix<scomplex,MC,  Star>& A,
                  const DistMatrix<scomplex,MR,  Star>& B,
  scomplex beta,        DistMatrix<scomplex,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRankK
( Shape shape,
  Orientation orientationOfA,
  Orientation orientationOfB,
  scomplex alpha, const DistMatrix<scomplex,Star,MC  >& A,
                  const DistMatrix<scomplex,MR,  Star>& B,
  scomplex beta,        DistMatrix<scomplex,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRankK
( Shape shape,
  scomplex alpha, const DistMatrix<scomplex,MC,  Star>& A,
                  const DistMatrix<scomplex,Star,MR  >& B,
  scomplex beta,        DistMatrix<scomplex,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRankK
( Shape shape,
  Orientation orientationOfA,
  scomplex alpha, const DistMatrix<scomplex,Star,MC>& A,
                  const DistMatrix<scomplex,Star,MR>& B,
  scomplex beta,        DistMatrix<scomplex,MC,  MR>& C );

template void
elemental::blas::internal::LocalTriangularRankK
( Shape shape,
  Orientation orientationOfB,
  dcomplex alpha, const DistMatrix<dcomplex,MC,  Star>& A,
                  const DistMatrix<dcomplex,MR,  Star>& B,
  dcomplex beta,        DistMatrix<dcomplex,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRankK
( Shape shape,
  Orientation orientationOfA,
  Orientation orientationOfB,
  dcomplex alpha, const DistMatrix<dcomplex,Star,MC  >& A,
                  const DistMatrix<dcomplex,MR,  Star>& B,
  dcomplex beta,        DistMatrix<dcomplex,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRankK
( Shape shape,
  dcomplex alpha, const DistMatrix<dcomplex,MC,  Star>& A,
                  const DistMatrix<dcomplex,Star,MR  >& B,
  dcomplex beta,        DistMatrix<dcomplex,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRankK
( Shape shape,
  Orientation orientationOfA,
  dcomplex alpha, const DistMatrix<dcomplex,Star,MC>& A,
                  const DistMatrix<dcomplex,Star,MR>& B,
  dcomplex beta,        DistMatrix<dcomplex,MC,  MR>& C );
#endif

