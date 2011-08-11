/*
   Copyright (c) 2009-2011, Jack Poulson
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
#include "elemental/basic_internal.hpp"
using namespace std;
using namespace elemental;

// Set up interface for managing tuning parameters
namespace {
int localTriangularRankKFloatBlocksize = 64;
int localTriangularRankKDoubleBlocksize = 64;
#ifndef WITHOUT_COMPLEX
int localTriangularRankKComplexFloatBlocksize = 64;
int localTriangularRankKComplexDoubleBlocksize = 64;
#endif // WITHOUT_COMPLEX
}

template<>
void
elemental::basic::SetLocalTriangularRankKBlocksize<float>
( int blocksize )
{ ::localTriangularRankKFloatBlocksize = blocksize; }

template<>
void
elemental::basic::SetLocalTriangularRankKBlocksize<double>
( int blocksize )
{ ::localTriangularRankKDoubleBlocksize = blocksize; }

#ifndef WITHOUT_COMPLEX
template<>
void
elemental::basic::SetLocalTriangularRankKBlocksize< std::complex<float> >
( int blocksize )
{ ::localTriangularRankKComplexFloatBlocksize = blocksize; }

template<>
void
elemental::basic::SetLocalTriangularRankKBlocksize< std::complex<double> >
( int blocksize )
{ ::localTriangularRankKComplexDoubleBlocksize = blocksize; }
#endif // WITHOUT_COMPLEX

template<>
int
elemental::basic::LocalTriangularRankKBlocksize<float>()
{ return ::localTriangularRankKFloatBlocksize; }

template<>
int
elemental::basic::LocalTriangularRankKBlocksize<double>()
{ return ::localTriangularRankKDoubleBlocksize; }

#ifndef WITHOUT_COMPLEX
template<>
int
elemental::basic::LocalTriangularRankKBlocksize<scomplex>()
{ return ::localTriangularRankKComplexFloatBlocksize; }

template<>
int
elemental::basic::LocalTriangularRankKBlocksize<dcomplex>()
{ return ::localTriangularRankKComplexDoubleBlocksize; }
#endif // WITHOUT_COMPLEX

// Template conventions:
//   G: general datatype
//
//   T: any ring, e.g., the (Gaussian) integers and the real/complex numbers
//   Z: representation of a real ring, e.g., the integers or real numbers
//   std::complex<Z>: representation of a complex ring, e.g. Gaussian integers
//                    or complex numbers
//
//   F: representation of real or complex number
//   R: representation of real number
//   std::complex<R>: representation of complex number

namespace {

#ifndef RELEASE
template<typename T>
void
CheckInput
( Orientation orientationOfB,
  const DistMatrix<T,MC,STAR>& A,
  const DistMatrix<T,MR,STAR>& B,
  const DistMatrix<T,MC,MR  >& C )
{
    if( orientationOfB == NORMAL )
        throw logic_error( "B[MR,* ] must be (Conjugate)Transpose'd." );
    if( A.Grid() != B.Grid() || B.Grid() != C.Grid() )
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
  const DistMatrix<T,STAR,MC  >& A,
  const DistMatrix<T,MR,  STAR>& B,
  const DistMatrix<T,MC,  MR  >& C )
{
    if( orientationOfA == NORMAL )
        throw logic_error( "A[* ,MC] must be (Conjugate)Transpose'd." );
    if( orientationOfB == NORMAL )
        throw logic_error( "B[MR,* ] must be (Conjugate)Transpose'd." );
    if( A.Grid() != B.Grid() || B.Grid() != C.Grid() )
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
( const DistMatrix<T,MC,  STAR>& A, 
  const DistMatrix<T,STAR,MR  >& B,
  const DistMatrix<T,MC,  MR  >& C )
{
    if( A.Grid() != B.Grid() || B.Grid() != C.Grid() )
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
  const DistMatrix<T,STAR,MC>& A,
  const DistMatrix<T,STAR,MR>& B,
  const DistMatrix<T,MC,  MR>& C )
{
    if( orientationOfA == NORMAL )
        throw logic_error( "A[* ,MC] must be (Conjugate)Transpose'd." );
    if( A.Grid() != B.Grid() || B.Grid() != C.Grid() )
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
  T alpha, const DistMatrix<T,MC,STAR>& A,
           const DistMatrix<T,MR,STAR>& B,
  T beta,        DistMatrix<T,MC,MR  >& C )
{
#ifndef RELEASE
    PushCallStack("LocalTriangularRankKKernel");
    CheckInput( orientationOfB, A, B, C );
#endif
    const Grid& g = C.Grid();

    DistMatrix<T,MC,STAR> AT(g),
                          AB(g);

    DistMatrix<T,MR,STAR> BT(g), 
                          BB(g);

    DistMatrix<T,MC,MR> CTL(g), CTR(g),
                        CBL(g), CBR(g);

    DistMatrix<T,MC,MR> DTL(g), DBR(g);

    const unsigned half = C.Height()/2;

    basic::Scal( beta, C );

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
    if( shape == LOWER )
    {
        basic::internal::LocalGemm
        ( NORMAL, orientationOfB, alpha, AB, BT, (T)1, CBL );
    }
    else
    {
        basic::internal::LocalGemm
        ( NORMAL, orientationOfB, alpha, AT, BB, (T)1, CTR );
    }

    basic::internal::LocalGemm
    ( NORMAL, orientationOfB, alpha, AT, BT, (T)0, DTL );

    // TODO: AxpyTrapezoidal?
    DTL.MakeTrapezoidal( LEFT, shape );
    basic::Axpy( (T)1, DTL, CTL );

    basic::internal::LocalGemm
    ( NORMAL, orientationOfB, alpha, AB, BB, (T)0, DBR );

    // TODO: AxpyTrapezoidal?
    DBR.MakeTrapezoidal( LEFT, shape );
    basic::Axpy( (T)1, DBR, CBR );
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
  T alpha, const DistMatrix<T,STAR,MC  >& A,
           const DistMatrix<T,MR,  STAR>& B,
  T beta,        DistMatrix<T,MC,  MR  >& C )
{
#ifndef RELEASE
    PushCallStack("LocalTriangularRankKKernel");
    CheckInput( orientationOfA, orientationOfB, A, B, C );
#endif
    const Grid& g = C.Grid();

    DistMatrix<T,STAR,MC> AL(g), AR(g);

    DistMatrix<T,MR,STAR> BT(g), 
                          BB(g);

    DistMatrix<T,MC,MR> CTL(g), CTR(g),
                        CBL(g), CBR(g);

    DistMatrix<T,MC,MR> DTL(g), DBR(g);

    const unsigned half = C.Height()/2;

    basic::Scal( beta, C );

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
    if( shape == LOWER )
    {
        basic::internal::LocalGemm
        ( orientationOfA, orientationOfB, alpha, AR, BT, (T)1, CBL );
    }
    else
    {
        basic::internal::LocalGemm
        ( orientationOfA, orientationOfB, alpha, AL, BB, (T)1, CTR );
    }

    basic::internal::LocalGemm
    ( orientationOfA, orientationOfB, alpha, AL, BT, (T)0, DTL );

    DTL.MakeTrapezoidal( LEFT, shape );
    basic::Axpy( (T)1, DTL, CTL );

    basic::internal::LocalGemm
    ( orientationOfA, orientationOfB, alpha, AR, BB, (T)0, DBR );

    DBR.MakeTrapezoidal( LEFT, shape );
    basic::Axpy( (T)1, DBR, CBR );
    //------------------------------------------------------------------------//
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
LocalTriangularRankKKernel
( Shape shape, 
  T alpha, const DistMatrix<T,MC,  STAR>& A,
           const DistMatrix<T,STAR,MR  >& B,
  T beta,        DistMatrix<T,MC,  MR  >& C )
{
#ifndef RELEASE
    PushCallStack("LocalTriangularRankKKernel");
    CheckInput( A, B, C );
#endif
    const Grid& g = C.Grid();

    DistMatrix<T,MC,STAR> AT(g), 
                          AB(g);

    DistMatrix<T,STAR,MR> BL(g), BR(g);

    DistMatrix<T,MC,MR> CTL(g), CTR(g),
                        CBL(g), CBR(g);

    DistMatrix<T,MC,MR> DTL(g), DBR(g);

    const unsigned half = C.Height()/2;

    basic::Scal( beta, C );

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
    if( shape == LOWER )
    {
        basic::internal::LocalGemm
        ( NORMAL, NORMAL, alpha, AB, BL, (T)1, CBL );
    }
    else
    {
        basic::internal::LocalGemm
        ( NORMAL, NORMAL, alpha, AT, BR, (T)1, CTR );
    }

    basic::internal::LocalGemm
    ( NORMAL, NORMAL, alpha, AT, BL, (T)0, DTL );

    DTL.MakeTrapezoidal( LEFT, shape );
    basic::Axpy( (T)1, DTL, CTL );

    basic::internal::LocalGemm
    ( NORMAL, NORMAL, alpha, AB, BR, (T)0, DBR );

    DBR.MakeTrapezoidal( LEFT, shape );
    basic::Axpy( (T)1, DBR, CBR );
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
  T alpha, const DistMatrix<T,STAR,MC>& A,
           const DistMatrix<T,STAR,MR>& B,
  T beta,        DistMatrix<T,MC,  MR>& C )
{
#ifndef RELEASE
    PushCallStack("LocalTriangularRankKKernel");
    CheckInput( orientationOfA, A, B, C );
#endif
    const Grid& g = C.Grid();

    DistMatrix<T,STAR,MC> AL(g), AR(g);

    DistMatrix<T,STAR,MR> BL(g), BR(g);

    DistMatrix<T,MC,MR> CTL(g), CTR(g),
                        CBL(g), CBR(g);

    DistMatrix<T,MC,MR> DTL(g), DBR(g);

    const unsigned half = C.Height()/2;

    basic::Scal( beta, C );

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
    if( shape == LOWER )
    {
        basic::internal::LocalGemm
        ( orientationOfA, NORMAL, alpha, AR, BL, (T)1, CBL );
    }
    else
    {
        basic::internal::LocalGemm
        ( orientationOfA, NORMAL, alpha, AL, BR, (T)1, CTR );
    }

    basic::internal::LocalGemm
    ( orientationOfA, NORMAL, alpha, AL, BL, (T)0, DTL );

    DTL.MakeTrapezoidal( LEFT, shape );
    basic::Axpy( (T)1, DTL, CTL );

    basic::internal::LocalGemm
    ( orientationOfA, NORMAL, alpha, AR, BR, (T)0, DBR );

    DBR.MakeTrapezoidal( LEFT, shape );
    basic::Axpy( (T)1, DBR, CBR );
    //------------------------------------------------------------------------//
#ifndef RELEASE
    PopCallStack();
#endif
}

} // anonymous namespace

template<typename T>
void
elemental::basic::internal::LocalTriangularRankK
( Shape shape,
  Orientation orientationOfB,
  T alpha, const DistMatrix<T,MC,STAR>& A,
           const DistMatrix<T,MR,STAR>& B,
  T beta,        DistMatrix<T,MC,MR  >& C )
{
#ifndef RELEASE
    PushCallStack("basic::internal::LocalTriangularRankK");
    CheckInput( orientationOfB, A, B, C );
#endif
    const Grid& g = C.Grid();

    if( C.Height() < g.Width()*LocalTriangularRankKBlocksize<T>() )
    {
        LocalTriangularRankKKernel
        ( shape, orientationOfB, alpha, A, B, beta, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.

        DistMatrix<T,MC,STAR> AT(g),
                              AB(g);

        DistMatrix<T,MR,STAR> BT(g), 
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

        if( shape == LOWER )
        { 
            basic::internal::LocalGemm
            ( NORMAL, orientationOfB, alpha, AB, BT, beta, CBL );
        }
        else
        {
            basic::internal::LocalGemm
            ( NORMAL, orientationOfB, alpha, AT, BB, beta, CTR );
        }

        // Recurse
        basic::internal::LocalTriangularRankK
        ( shape, orientationOfB, alpha, AT, BT, beta, CTL );

        basic::internal::LocalTriangularRankK
        ( shape, orientationOfB, alpha, AB, BB, beta, CBR );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::basic::internal::LocalTriangularRankK
( Shape shape,
  Orientation orientationOfA,
  Orientation orientationOfB,
  T alpha, const DistMatrix<T,STAR,MC  >& A,
           const DistMatrix<T,MR,  STAR>& B,
  T beta,        DistMatrix<T,MC,  MR  >& C )
{
#ifndef RELEASE
    PushCallStack("basic::internal::LocalTriangularRankK");
    CheckInput( orientationOfA, orientationOfB, A, B, C );
#endif
    const Grid& g = C.Grid();

    if( C.Height() < g.Width()*LocalTriangularRankKBlocksize<T>() )
    {
        LocalTriangularRankKKernel
        ( shape, orientationOfA, orientationOfB, alpha, A, B, beta, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.

        DistMatrix<T,STAR,MC> AL(g), AR(g);

        DistMatrix<T,MR,STAR> BT(g), 
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

        if( shape == LOWER )
        { 
            basic::internal::LocalGemm
            ( orientationOfA, orientationOfB, alpha, AR, BT, beta, CBL );
        }
        else
        {
            basic::internal::LocalGemm
            ( orientationOfA, orientationOfB, alpha, AL, BB, beta, CTR );
        }

        // Recurse
        basic::internal::LocalTriangularRankK
        ( shape, orientationOfA, orientationOfB, alpha, AL, BT, beta, CTL );

        basic::internal::LocalTriangularRankK
        ( shape, orientationOfA, orientationOfB, alpha, AR, BB, beta, CBR );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::basic::internal::LocalTriangularRankK
( Shape shape,
  T alpha, const DistMatrix<T,MC,  STAR>& A,
           const DistMatrix<T,STAR,MR  >& B,
  T beta,        DistMatrix<T,MC,  MR  >& C )
{
#ifndef RELEASE
    PushCallStack("basic::internal::LocalTriangularRankK");
    CheckInput( A, B, C );
#endif
    const Grid& g = C.Grid();

    if( C.Height() < g.Width()*LocalTriangularRankKBlocksize<T>() )
    {
        LocalTriangularRankKKernel
        ( shape, alpha, A, B, beta, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.

        DistMatrix<T,MC,STAR> AT(g),
                              AB(g);

        DistMatrix<T,STAR,MR> BL(g), BR(g);

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

        if( shape == LOWER )
        { 
            basic::internal::LocalGemm
            ( NORMAL, NORMAL, alpha, AB, BL, beta, CBL );
        }
        else
        {
            basic::internal::LocalGemm
            ( NORMAL, NORMAL, alpha, AT, BR, beta, CTR );
        }

        // Recurse
        basic::internal::LocalTriangularRankK
        ( shape, alpha, AT, BL, beta, CTL );

        basic::internal::LocalTriangularRankK
        ( shape, alpha, AB, BR, beta, CBR );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::basic::internal::LocalTriangularRankK
( Shape shape,
  Orientation orientationOfA,
  T alpha, const DistMatrix<T,STAR,MC>& A,
           const DistMatrix<T,STAR,MR>& B,
  T beta,        DistMatrix<T,MC,  MR>& C )
{
#ifndef RELEASE
    PushCallStack("basic::internal::LocalTriangularRankK");
    CheckInput( orientationOfA, A, B, C );
#endif
    const Grid& g = C.Grid();

    if( C.Height() < g.Width()*LocalTriangularRankKBlocksize<T>() )
    {
        LocalTriangularRankKKernel
        ( shape, orientationOfA, alpha, A, B, beta, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.

        DistMatrix<T,STAR,MC> AL(g), AR(g);

        DistMatrix<T,STAR,MR> BL(g), BR(g);

        DistMatrix<T,MC,MR> CTL(g), CTR(g),
                            CBL(g), CBR(g);

        const unsigned half = C.Height() / 2;

        LockedPartitionRight( A, AL, AR, half );

        LockedPartitionRight( B, BL, BR, half );

        PartitionDownDiagonal
        ( C, CTL, CTR,
             CBL, CBR, half );

        if( shape == LOWER )
        { 
            basic::internal::LocalGemm
            ( orientationOfA, NORMAL, alpha, AR, BL, beta, CBL );
        }
        else
        {
            basic::internal::LocalGemm
            ( orientationOfA, NORMAL, alpha, AL, BR, beta, CTR );
        }

        // Recurse
        basic::internal::LocalTriangularRankK
        ( shape, orientationOfA, alpha, AL, BL, beta, CTL );

        basic::internal::LocalTriangularRankK
        ( shape, orientationOfA, alpha, AR, BR, beta, CBR );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void
elemental::basic::internal::LocalTriangularRankK
( Shape shape,
  Orientation orientationOfB,
  float alpha, const DistMatrix<float,MC,  STAR>& A,
               const DistMatrix<float,MR,  STAR>& B,
  float beta,        DistMatrix<float,MC,  MR  >& C );

template void
elemental::basic::internal::LocalTriangularRankK
( Shape shape, 
  Orientation orientationOfA,
  Orientation orientationOfB,
  float alpha, const DistMatrix<float,STAR,MC  >& A,
               const DistMatrix<float,MR,  STAR>& B,
  float beta,        DistMatrix<float,MC,  MR  >& C );

template void
elemental::basic::internal::LocalTriangularRankK
( Shape shape, 
  float alpha, const DistMatrix<float,MC,  STAR>& A,
               const DistMatrix<float,STAR,MR  >& B,
  float beta,        DistMatrix<float,MC,  MR  >& C );

template void
elemental::basic::internal::LocalTriangularRankK
( Shape shape, 
  Orientation orientationOfA,
  float alpha, const DistMatrix<float,STAR,MC>& A,
               const DistMatrix<float,STAR,MR>& B,
  float beta,        DistMatrix<float,MC,  MR>& C );

template void
elemental::basic::internal::LocalTriangularRankK
( Shape shape,
  Orientation orientationOfB,
  double alpha, const DistMatrix<double,MC,  STAR>& A,
                const DistMatrix<double,MR,  STAR>& B,
  double beta,        DistMatrix<double,MC,  MR  >& C );

template void
elemental::basic::internal::LocalTriangularRankK
( Shape shape,
  Orientation orientationOfA,
  Orientation orientationOfB,
  double alpha, const DistMatrix<double,STAR,MC  >& A,
                const DistMatrix<double,MR,  STAR>& B,
  double beta,        DistMatrix<double,MC,  MR  >& C );

template void
elemental::basic::internal::LocalTriangularRankK
( Shape shape,
  double alpha, const DistMatrix<double,MC,  STAR>& A,
                const DistMatrix<double,STAR,MR  >& B,
  double beta,        DistMatrix<double,MC,  MR  >& C );

template void
elemental::basic::internal::LocalTriangularRankK
( Shape shape,
  Orientation orientationOfA,
  double alpha, const DistMatrix<double,STAR,MC>& A,
                const DistMatrix<double,STAR,MR>& B,
  double beta,        DistMatrix<double,MC,  MR>& C );


#ifndef WITHOUT_COMPLEX
template void
elemental::basic::internal::LocalTriangularRankK
( Shape shape,
  Orientation orientationOfB,
  scomplex alpha, const DistMatrix<scomplex,MC,  STAR>& A,
                  const DistMatrix<scomplex,MR,  STAR>& B,
  scomplex beta,        DistMatrix<scomplex,MC,  MR  >& C );

template void
elemental::basic::internal::LocalTriangularRankK
( Shape shape,
  Orientation orientationOfA,
  Orientation orientationOfB,
  scomplex alpha, const DistMatrix<scomplex,STAR,MC  >& A,
                  const DistMatrix<scomplex,MR,  STAR>& B,
  scomplex beta,        DistMatrix<scomplex,MC,  MR  >& C );

template void
elemental::basic::internal::LocalTriangularRankK
( Shape shape,
  scomplex alpha, const DistMatrix<scomplex,MC,  STAR>& A,
                  const DistMatrix<scomplex,STAR,MR  >& B,
  scomplex beta,        DistMatrix<scomplex,MC,  MR  >& C );

template void
elemental::basic::internal::LocalTriangularRankK
( Shape shape,
  Orientation orientationOfA,
  scomplex alpha, const DistMatrix<scomplex,STAR,MC>& A,
                  const DistMatrix<scomplex,STAR,MR>& B,
  scomplex beta,        DistMatrix<scomplex,MC,  MR>& C );

template void
elemental::basic::internal::LocalTriangularRankK
( Shape shape,
  Orientation orientationOfB,
  dcomplex alpha, const DistMatrix<dcomplex,MC,  STAR>& A,
                  const DistMatrix<dcomplex,MR,  STAR>& B,
  dcomplex beta,        DistMatrix<dcomplex,MC,  MR  >& C );

template void
elemental::basic::internal::LocalTriangularRankK
( Shape shape,
  Orientation orientationOfA,
  Orientation orientationOfB,
  dcomplex alpha, const DistMatrix<dcomplex,STAR,MC  >& A,
                  const DistMatrix<dcomplex,MR,  STAR>& B,
  dcomplex beta,        DistMatrix<dcomplex,MC,  MR  >& C );

template void
elemental::basic::internal::LocalTriangularRankK
( Shape shape,
  dcomplex alpha, const DistMatrix<dcomplex,MC,  STAR>& A,
                  const DistMatrix<dcomplex,STAR,MR  >& B,
  dcomplex beta,        DistMatrix<dcomplex,MC,  MR  >& C );

template void
elemental::basic::internal::LocalTriangularRankK
( Shape shape,
  Orientation orientationOfA,
  dcomplex alpha, const DistMatrix<dcomplex,STAR,MC>& A,
                  const DistMatrix<dcomplex,STAR,MR>& B,
  dcomplex beta,        DistMatrix<dcomplex,MC,  MR>& C );
#endif

