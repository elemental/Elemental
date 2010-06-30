/*
   This file is part of Elemental, a library for distributed-memory dense 
   linear algebra.

   Copyright (c) 2009-2010 Jack Poulson <jack.poulson@gmail.com>.
   All rights reserved.

   This file is released under the terms of the license contained in the file
   LICENSE-PURE.
*/
#include "elemental/blas_internal.hpp"
using namespace std;
using namespace elemental;

namespace {

#ifndef RELEASE
template<typename T>
void 
CheckInput
( Orientation orientationOfB1,
  Orientation orientationOfB2,
  const DistMatrix<T,MC,Star>& A1, 
  const DistMatrix<T,MC,Star>& A2,
  const DistMatrix<T,MR,Star>& B1,
  const DistMatrix<T,MR,Star>& B2,
  const DistMatrix<T,MC,MR  >& C )
{
    if( orientationOfB1 == Normal )
        throw logic_error( "B1[MR,* ] must be (Conjugate)Transpose'd." );
    if( orientationOfB2 == Normal )
        throw logic_error( "B2[MR,* ] must be (Conjugate)Transpose'd." );
    if( A1.GetGrid() != A2.GetGrid() || A2.GetGrid() != B1.GetGrid() ||
        B1.GetGrid() != B2.GetGrid() || B2.GetGrid() != C.GetGrid() )
        throw logic_error
        ( "A, B, and C must be distributed over the same grid." );
    if( A1.Height() != C.Height() || B1.Height() != C.Width() ||
        A1.Height() != A2.Height() || A1.Width() != A2.Width() ||
        B1.Width() != B2.Width() || B1.Height() != B2.Height() ||
        A1.Width() != B1.Width() )
    {
        ostringstream msg;
        msg << "Nonconformal LocalTriangularRank2K: " << endl
            << "  A1[MC,* ] ~ " << A1.Height() << " x "
                                << A1.Width()  << endl
            << "  A2[MC,* ] ~ " << A2.Height() << " x "
                                << A2.Width()  << endl
            << "  B1[MR,* ] ~ " << B1.Height() << " x "
                                << B1.Width()  << endl
            << "  B2[MR,* ] ~ " << B2.Height() << " x "
                                << B2.Width()  << endl
            << "  C[MC,MR] ~ " << C.Height() << " x " << C.Width() << endl;
        throw logic_error( msg.str() );
    }
    if( A1.ColAlignment() != C.ColAlignment() ||
        B1.ColAlignment() != C.RowAlignment() ||
        A1.ColAlignment() != A2.ColAlignment() ||
        B1.ColAlignment() != B2.ColAlignment() )
    {
        ostringstream msg;
        msg << "Misaligned LocalTriangularRank2K: " << endl
            << "  A1[MC,* ] ~ " << A1.ColAlignment() << endl
            << "  A2[MC,* ] ~ " << A2.ColAlignment() << endl
            << "  B1[MR,* ] ~ " << B1.ColAlignment() << endl
            << "  B2[MR,* ] ~ " << B2.ColAlignment() << endl
            << "  C[MC,MR] ~ " << C.ColAlignment() << " , " <<
                                  C.RowAlignment() << endl;
        throw logic_error( msg.str() );
    }
}

template<typename T>
void 
CheckInput
( Orientation orientationOfA1,
  Orientation orientationOfB1,
  Orientation orientationOfB2,
  const DistMatrix<T,Star,MC  >& A1, 
  const DistMatrix<T,MC,  Star>& A2,
  const DistMatrix<T,MR,  Star>& B1,
  const DistMatrix<T,MR,  Star>& B2,
  const DistMatrix<T,MC,  MR  >& C )
{
    if( orientationOfA1 == Normal )
        throw logic_error( "A1[* ,MC] must be (Conjugate)Transpose'd." );
    if( orientationOfB1 == Normal )
        throw logic_error( "B1[MR,* ] must be (Conjugate)Transpose'd." );
    if( orientationOfB2 == Normal )
        throw logic_error( "B2[MR,* ] must be (Conjugate)Transpose'd." );
    if( A1.GetGrid() != A2.GetGrid() || A2.GetGrid() != B1.GetGrid() ||
        B1.GetGrid() != B2.GetGrid() || B2.GetGrid() != C.GetGrid() )
        throw logic_error
        ( "A, B, and C must be distributed over the same grid." );
    if( A1.Width() != C.Height() || B1.Height() != C.Width() ||
        A1.Width() != A2.Height() || A1.Height() != A2.Width() ||
        B1.Width() != B2.Width() || B1.Height() != B2.Height() ||
        A1.Height() != B1.Width() )
    {
        ostringstream msg;
        msg << "Nonconformal LocalTriangularRank2K: " << endl
            << "  A1[* ,MC] ~ " << A1.Height() << " x "
                                << A1.Width()  << endl
            << "  A2[MC,* ] ~ " << A2.Height() << " x "
                                << A2.Width()  << endl
            << "  B1[MR,* ] ~ " << B1.Height() << " x "
                                << B1.Width()  << endl
            << "  B2[MR,* ] ~ " << B2.Height() << " x "
                                << B2.Width()  << endl
            << "  C[MC,MR] ~ " << C.Height() << " x " << C.Width() << endl;
        throw logic_error( msg.str() );
    }
    if( A1.RowAlignment() != C.ColAlignment() ||
        B1.ColAlignment() != C.RowAlignment() ||
        A1.RowAlignment() != A2.ColAlignment() ||
        B1.ColAlignment() != B2.ColAlignment() )
    {
        ostringstream msg;
        msg << "Misaligned LocalTriangularRank2K: " << endl
            << "  A1[* ,MC] ~ " << A1.RowAlignment() << endl
            << "  A2[MC,* ] ~ " << A2.ColAlignment() << endl
            << "  B1[MR,* ] ~ " << B1.ColAlignment() << endl
            << "  B2[MR,* ] ~ " << B2.ColAlignment() << endl
            << "  C[MC,MR] ~ " << C.ColAlignment() << " , " <<
                                  C.RowAlignment() << endl;
        throw logic_error( msg.str() );
    }
}

template<typename T>
void 
CheckInput
( Orientation orientationOfA2,
  Orientation orientationOfB1,
  Orientation orientationOfB2,
  const DistMatrix<T,MC,  Star>& A1, 
  const DistMatrix<T,Star,MC  >& A2,
  const DistMatrix<T,MR,  Star>& B1,
  const DistMatrix<T,MR,  Star>& B2,
  const DistMatrix<T,MC,  MR  >& C )
{
    if( orientationOfA2 == Normal )
        throw logic_error( "A2[* ,MC] must be (Conjugate)Transpose'd." );
    if( orientationOfB1 == Normal )
        throw logic_error( "B1[MR,* ] must be (Conjugate)Transpose'd." );
    if( orientationOfB2 == Normal )
        throw logic_error( "B2[MR,* ] must be (Conjugate)Transpose'd." );
    if( A1.GetGrid() != A2.GetGrid() || A2.GetGrid() != B1.GetGrid() ||
        B1.GetGrid() != B2.GetGrid() || B2.GetGrid() != C.GetGrid() )
        throw logic_error
        ( "A, B, and C must be distributed over the same grid." );
    if( A1.Height() != C.Height() || B1.Height() != C.Width() ||
        A1.Height() != A2.Width() || A1.Width() != A2.Height() ||
        B1.Width() != B2.Width() || B1.Height() != B2.Height() ||
        A1.Width() != B1.Width() )
    {
        ostringstream msg;
        msg << "Nonconformal LocalTriangularRank2K: " << endl
            << "  A1[MC,* ] ~ " << A1.Height() << " x "
                                << A1.Width()  << endl
            << "  A2[* ,MC] ~ " << A2.Height() << " x "
                                << A2.Width()  << endl
            << "  B1[MR,* ] ~ " << B1.Height() << " x "
                                << B1.Width()  << endl
            << "  B2[MR,* ] ~ " << B2.Height() << " x "
                                << B2.Width()  << endl
            << "  C[MC,MR] ~ " << C.Height() << " x " << C.Width() << endl;
        throw logic_error( msg.str() );
    }
    if( A1.ColAlignment() != C.ColAlignment() ||
        B1.ColAlignment() != C.RowAlignment() ||
        A1.ColAlignment() != A2.RowAlignment() ||
        B1.ColAlignment() != B2.ColAlignment() )
    {
        ostringstream msg;
        msg << "Misaligned LocalTriangularRank2K: " << endl
            << "  A1[MC,* ] ~ " << A1.ColAlignment() << endl
            << "  A2[* ,MC] ~ " << A2.RowAlignment() << endl
            << "  B1[MR,* ] ~ " << B1.ColAlignment() << endl
            << "  B2[MR,* ] ~ " << B2.ColAlignment() << endl
            << "  C[MC,MR] ~ " << C.ColAlignment() << " , " <<
                                  C.RowAlignment() << endl;
        throw logic_error( msg.str() );
    }
}

template<typename T>
void 
CheckInput
( Orientation orientationOfB2,
  const DistMatrix<T,MC,  Star>& A1, 
  const DistMatrix<T,MC,  Star>& A2,
  const DistMatrix<T,Star,MR  >& B1,
  const DistMatrix<T,MR,  Star>& B2,
  const DistMatrix<T,MC,  MR  >& C )
{
    if( orientationOfB2 == Normal )
        throw logic_error( "B2[MR,* ] must be (Conjugate)Transpose'd." );
    if( A1.GetGrid() != A2.GetGrid() || A2.GetGrid() != B1.GetGrid() ||
        B1.GetGrid() != B2.GetGrid() || B2.GetGrid() != C.GetGrid() )
        throw logic_error
        ( "A, B, and C must be distributed over the same grid." );
    if( A1.Height() != C.Height() || B1.Height() != C.Width() ||
        A1.Height() != A2.Height() || A1.Width() != A2.Width() ||
        B1.Height() != B2.Width() || B1.Width() != B2.Height() ||
        A1.Width() != B1.Height() )
    {
        ostringstream msg;
        msg << "Nonconformal LocalTriangularRank2K: " << endl
            << "  A1[MC,* ] ~ " << A1.Height() << " x "
                                << A1.Width()  << endl
            << "  A2[MC,* ] ~ " << A2.Height() << " x "
                                << A2.Width()  << endl
            << "  B1[* ,MR] ~ " << B1.Height() << " x "
                                << B1.Width()  << endl
            << "  B2[MR,* ] ~ " << B2.Height() << " x "
                                << B2.Width()  << endl
            << "  C[MC,MR] ~ " << C.Height() << " x " << C.Width() << endl;
        throw logic_error( msg.str() );
    }
    if( A1.ColAlignment() != C.ColAlignment() ||
        B1.RowAlignment() != C.RowAlignment() ||
        A1.ColAlignment() != A2.ColAlignment() ||
        B1.RowAlignment() != B2.ColAlignment() )
    {
        ostringstream msg;
        msg << "Misaligned LocalTriangularRank2K: " << endl
            << "  A1[MC,* ] ~ " << A1.ColAlignment() << endl
            << "  A2[MC,* ] ~ " << A2.ColAlignment() << endl
            << "  B1[* ,MR] ~ " << B1.RowAlignment() << endl
            << "  B2[MR,* ] ~ " << B2.ColAlignment() << endl
            << "  C[MC,MR] ~ " << C.ColAlignment() << " , " <<
                                  C.RowAlignment() << endl;
        throw logic_error( msg.str() );
    }
}

template<typename T>
void 
CheckInput
( Orientation orientationOfB1,
  const DistMatrix<T,MC,  Star>& A1, 
  const DistMatrix<T,MC,  Star>& A2,
  const DistMatrix<T,MR,  Star>& B1,
  const DistMatrix<T,Star,MR  >& B2,
  const DistMatrix<T,MC,  MR  >& C )
{
    if( orientationOfB1 == Normal )
        throw logic_error( "B1[MR,* ] must be (Conjugate)Transpose'd." );
    if( A1.GetGrid() != A2.GetGrid() || A2.GetGrid() != B1.GetGrid() ||
        B1.GetGrid() != B2.GetGrid() || B2.GetGrid() != C.GetGrid() )
        throw logic_error
        ( "A, B, and C must be distributed over the same grid." );
    if( A1.Height() != C.Height() || B1.Height() != C.Width() ||
        A1.Height() != A2.Height() || A1.Width() != A2.Width() ||
        B1.Width() != B2.Height() || B1.Height() != B2.Width() ||
        A1.Width() != B1.Width() )
    {
        ostringstream msg;
        msg << "Nonconformal LocalTriangularRank2K: " << endl
            << "  A1[MC,* ] ~ " << A1.Height() << " x "
                                << A1.Width()  << endl
            << "  A2[MC,* ] ~ " << A2.Height() << " x "
                                << A2.Width()  << endl
            << "  B1[MR,* ] ~ " << B1.Height() << " x "
                                << B1.Width()  << endl
            << "  B2[* ,MR] ~ " << B2.Height() << " x "
                                << B2.Width()  << endl
            << "  C[MC,MR] ~ " << C.Height() << " x " << C.Width() << endl;
        throw logic_error( msg.str() );
    }
    if( A1.ColAlignment() != C.ColAlignment() ||
        B1.ColAlignment() != C.RowAlignment() ||
        A1.ColAlignment() != A2.ColAlignment() ||
        B1.ColAlignment() != B2.RowAlignment() )
    {
        ostringstream msg;
        msg << "Misaligned LocalTriangularRank2K: " << endl
            << "  A1[MC,* ] ~ " << A1.ColAlignment() << endl
            << "  A2[MC,* ] ~ " << A2.ColAlignment() << endl
            << "  B1[MR,* ] ~ " << B1.ColAlignment() << endl
            << "  B2[* ,MR] ~ " << B2.RowAlignment() << endl
            << "  C[MC,MR] ~ " << C.ColAlignment() << " , " <<
                                  C.RowAlignment() << endl;
        throw logic_error( msg.str() );
    }
}

template<typename T>
void 
CheckInput
( Orientation orientationOfA1,
  Orientation orientationOfA2,
  Orientation orientationOfB1,
  Orientation orientationOfB2,
  const DistMatrix<T,Star,MC  >& A1, 
  const DistMatrix<T,Star,MC  >& A2,
  const DistMatrix<T,MR,  Star>& B1,
  const DistMatrix<T,MR,  Star>& B2,
  const DistMatrix<T,MC,  MR  >& C )
{
    if( orientationOfA1 == Normal )
        throw logic_error( "A1[* ,MC] must be (Conjugate)Transpose'd." );
    if( orientationOfA2 == Normal )
        throw logic_error( "A2[* ,MC] must be (Conjugate)Transpose'd." );
    if( orientationOfB1 == Normal )
        throw logic_error( "B1[MR,* ] must be (Conjugate)Transpose'd." );
    if( orientationOfB2 == Normal )
        throw logic_error( "B2[MR,* ] must be (Conjugate)Transpose'd." );
    if( A1.GetGrid() != A2.GetGrid() || A2.GetGrid() != B1.GetGrid() ||
        B1.GetGrid() != B2.GetGrid() || B2.GetGrid() != C.GetGrid() )
        throw logic_error
        ( "A, B, and C must be distributed over the same grid." );
    if( A1.Width() != C.Height() || B1.Height() != C.Width() ||
        A1.Height() != A2.Height() || A1.Width() != A2.Width() ||
        B1.Height() != B2.Height() || B1.Width() != B2.Width() ||
        A1.Height() != B1.Width() )
    {
        ostringstream msg;
        msg << "Nonconformal LocalTriangularRank2K: " << endl
            << "  A1[* ,MC] ~ " << A1.Height() << " x "
                                << A1.Width()  << endl
            << "  A2[* ,MC] ~ " << A2.Height() << " x "
                                << A2.Width()  << endl
            << "  B1[MR,* ] ~ " << B1.Height() << " x "
                                << B1.Width()  << endl
            << "  B2[MR,* ] ~ " << B2.Height() << " x "
                                << B2.Width()  << endl
            << "  C[MC,MR] ~ " << C.Height() << " x " << C.Width() << endl;
        throw logic_error( msg.str() );
    }
    if( A1.RowAlignment() != C.ColAlignment() ||
        B1.ColAlignment() != C.RowAlignment() ||
        A1.RowAlignment() != A2.RowAlignment() ||
        B1.ColAlignment() != B2.ColAlignment() )
    {
        ostringstream msg;
        msg << "Misaligned LocalTriangularRank2K: " << endl
            << "  A1[* ,MC] ~ " << A1.RowAlignment() << endl
            << "  A2[* ,MC] ~ " << A2.RowAlignment() << endl
            << "  B1[MR,* ] ~ " << B1.ColAlignment() << endl
            << "  B2[MR,* ] ~ " << B2.ColAlignment() << endl
            << "  C[MC,MR] ~ " << C.ColAlignment() << " , " <<
                                  C.RowAlignment() << endl;
        throw logic_error( msg.str() );
    }
}

template<typename T>
void 
CheckInput
( Orientation orientationOfA1,
  Orientation orientationOfB2,
  const DistMatrix<T,Star,MC  >& A1, 
  const DistMatrix<T,MC,  Star>& A2,
  const DistMatrix<T,Star,MR  >& B1,
  const DistMatrix<T,MR,  Star>& B2,
  const DistMatrix<T,MC,  MR  >& C )
{
    if( orientationOfA1 == Normal )
        throw logic_error( "A1[* ,MC] must be (Conjugate)Transpose'd." );
    if( orientationOfB2 == Normal )
        throw logic_error( "B2[MR,* ] must be (Conjugate)Transpose'd." );
    if( A1.GetGrid() != A2.GetGrid() || A2.GetGrid() != B1.GetGrid() ||
        B1.GetGrid() != B2.GetGrid() || B2.GetGrid() != C.GetGrid() )
        throw logic_error
        ( "A, B, and C must be distributed over the same grid." );
    if( A1.Width() != C.Height() || B1.Width() != C.Width() ||
        A1.Width() != A2.Height() || A1.Height() != A2.Width() ||
        B1.Height() != B2.Width() || B1.Width() != B2.Height() ||
        A1.Height() != B1.Height() )
    {
        ostringstream msg;
        msg << "Nonconformal LocalTriangularRank2K: " << endl
            << "  A1[* ,MC] ~ " << A1.Height() << " x "
                                << A1.Width()  << endl
            << "  A2[MC,* ] ~ " << A2.Height() << " x "
                                << A2.Width()  << endl
            << "  B1[* ,MR] ~ " << B1.Height() << " x "
                                << B1.Width()  << endl
            << "  B2[MR,* ] ~ " << B2.Height() << " x "
                                << B2.Width()  << endl
            << "  C[MC,MR] ~ " << C.Height() << " x " << C.Width() << endl;
        throw logic_error( msg.str() );
    }
    if( A1.RowAlignment() != C.ColAlignment() ||
        B1.RowAlignment() != C.RowAlignment() ||
        A1.RowAlignment() != A2.ColAlignment() ||
        B1.RowAlignment() != B2.ColAlignment() )
    {
        ostringstream msg;
        msg << "Misaligned LocalTriangularRank2K: " << endl
            << "  A1[* ,MC] ~ " << A1.RowAlignment() << endl
            << "  A2[MC,* ] ~ " << A2.ColAlignment() << endl
            << "  B1[* ,MR] ~ " << B1.RowAlignment() << endl
            << "  B2[MR,* ] ~ " << B2.ColAlignment() << endl
            << "  C[MC,MR] ~ " << C.ColAlignment() << " , " <<
                                  C.RowAlignment() << endl;
        throw logic_error( msg.str() );
    }
}

template<typename T>
void 
CheckInput
( Orientation orientationOfA1,
  Orientation orientationOfB1,
  const DistMatrix<T,Star,MC  >& A1, 
  const DistMatrix<T,MC,  Star>& A2,
  const DistMatrix<T,MR,  Star>& B1,
  const DistMatrix<T,Star,MR  >& B2,
  const DistMatrix<T,MC,  MR  >& C )
{
    if( orientationOfA1 == Normal )
        throw logic_error( "A1[* ,MC] must be (Conjugate)Transpose'd." );
    if( orientationOfB1 == Normal )
        throw logic_error( "B1[MR,* ] must be (Conjugate)Transpose'd." );
    if( A1.GetGrid() != A2.GetGrid() || A2.GetGrid() != B1.GetGrid() ||
        B1.GetGrid() != B2.GetGrid() || B2.GetGrid() != C.GetGrid() )
        throw logic_error
        ( "A, B, and C must be distributed over the same grid." );
    if( A1.Width() != C.Height() || B1.Height() != C.Width() ||
        A1.Width() != A2.Height() || A1.Height() != A2.Width() ||
        B1.Width() != B2.Height() || B1.Height() != B2.Width() ||
        A1.Height() != B1.Width() )
    {
        ostringstream msg;
        msg << "Nonconformal LocalTriangularRank2K: " << endl
            << "  A1[* ,MC] ~ " << A1.Height() << " x "
                                << A1.Width()  << endl
            << "  A2[MC,* ] ~ " << A2.Height() << " x "
                                << A2.Width()  << endl
            << "  B1[MR,* ] ~ " << B1.Height() << " x "
                                << B1.Width()  << endl
            << "  B2[* ,MR] ~ " << B2.Height() << " x "
                                << B2.Width()  << endl
            << "  C[MC,MR] ~ " << C.Height() << " x " << C.Width() << endl;
        throw logic_error( msg.str() );
    }
    if( A1.RowAlignment() != C.ColAlignment() ||
        B1.ColAlignment() != C.RowAlignment() ||
        A1.RowAlignment() != A2.ColAlignment() ||
        B1.ColAlignment() != B2.RowAlignment() )
    {
        ostringstream msg;
        msg << "Misaligned LocalTriangularRank2K: " << endl
            << "  A1[* ,MC] ~ " << A1.RowAlignment() << endl
            << "  A2[MC,* ] ~ " << A2.ColAlignment() << endl
            << "  B1[MR,* ] ~ " << B1.ColAlignment() << endl
            << "  B2[* ,MR] ~ " << B2.RowAlignment() << endl
            << "  C[MC,MR] ~ " << C.ColAlignment() << " , " <<
                                  C.RowAlignment() << endl;
        throw logic_error( msg.str() );
    }
}

template<typename T>
void 
CheckInput
( Orientation orientationOfA2,
  Orientation orientationOfB2,
  const DistMatrix<T,MC,  Star>& A1,
  const DistMatrix<T,Star,MC  >& A2, 
  const DistMatrix<T,Star,MR  >& B1,
  const DistMatrix<T,MR,  Star>& B2,
  const DistMatrix<T,MC,  MR  >& C )
{
    if( orientationOfA2 == Normal )
        throw logic_error( "A2[* ,MC] must be (Conjugate)Transpose'd." );
    if( orientationOfB2 == Normal )
        throw logic_error( "B2[MR,* ] must be (Conjugate)Transpose'd." );
    if( A1.GetGrid() != A2.GetGrid() || A2.GetGrid() != B1.GetGrid() ||
        B1.GetGrid() != B2.GetGrid() || B2.GetGrid() != C.GetGrid() )
        throw logic_error
        ( "A, B, and C must be distributed over the same grid." );
    if( A1.Height() != C.Height() || B1.Width() != C.Width() ||
        A1.Height() != A2.Width() || A1.Width() != A2.Height() ||
        B1.Height() != B2.Width() || B1.Width() != B2.Height() ||
        A1.Width() != B1.Height() )
    {
        ostringstream msg;
        msg << "Nonconformal LocalTriangularRank2K: " << endl
            << "  A1[MC,* ] ~ " << A1.Height() << " x "
                                << A1.Width()  << endl
            << "  A2[* ,MC] ~ " << A2.Height() << " x "
                                << A2.Width()  << endl
            << "  B1[* ,MR] ~ " << B1.Height() << " x "
                                << B1.Width()  << endl
            << "  B2[MR,* ] ~ " << B2.Height() << " x "
                                << B2.Width()  << endl
            << "  C[MC,MR] ~ " << C.Height() << " x " << C.Width() << endl;
        throw logic_error( msg.str() );
    }
    if( A1.ColAlignment() != C.ColAlignment() ||
        B1.RowAlignment() != C.RowAlignment() ||
        A1.ColAlignment() != A2.RowAlignment() ||
        B1.RowAlignment() != B2.ColAlignment() )
    {
        ostringstream msg;
        msg << "Misaligned LocalTriangularRank2K: " << endl
            << "  A1[MC,* ] ~ " << A1.ColAlignment() << endl
            << "  A2[* ,MC] ~ " << A2.RowAlignment() << endl
            << "  B1[* ,MR] ~ " << B1.RowAlignment() << endl
            << "  B2[MR,* ] ~ " << B2.ColAlignment() << endl
            << "  C[MC,MR] ~ " << C.ColAlignment() << " , " <<
                                  C.RowAlignment() << endl;
        throw logic_error( msg.str() );
    }
}

template<typename T>
void 
CheckInput
( Orientation orientationOfA2,
  Orientation orientationOfB1,
  const DistMatrix<T,MC,  Star>& A1,
  const DistMatrix<T,Star,MC  >& A2, 
  const DistMatrix<T,MR,  Star>& B1,
  const DistMatrix<T,Star,MR  >& B2,
  const DistMatrix<T,MC,  MR  >& C )
{
    if( orientationOfA2 == Normal )
        throw logic_error( "A2[* ,MC] must be (Conjugate)Transpose'd." );
    if( orientationOfB1 == Normal )
        throw logic_error( "B1[MR,* ] must be (Conjugate)Transpose'd." );
    if( A1.GetGrid() != A2.GetGrid() || A2.GetGrid() != B1.GetGrid() ||
        B1.GetGrid() != B2.GetGrid() || B2.GetGrid() != C.GetGrid() )
        throw logic_error
        ( "A, B, and C must be distributed over the same grid." );
    if( A1.Height() != C.Height() || B1.Height() != C.Width() ||
        A1.Height() != A2.Width() || A1.Width() != A2.Height() ||
        B1.Width() != B2.Height() || B1.Height() != B2.Width() ||
        A1.Width() != B1.Width() )
    {
        ostringstream msg;
        msg << "Nonconformal LocalTriangularRank2K: " << endl
            << "  A1[MC,* ] ~ " << A1.Height() << " x "
                                << A1.Width()  << endl
            << "  A2[* ,MC] ~ " << A2.Height() << " x "
                                << A2.Width()  << endl
            << "  B1[MR,* ] ~ " << B1.Height() << " x "
                                << B1.Width()  << endl
            << "  B2[* ,MR] ~ " << B2.Height() << " x "
                                << B2.Width()  << endl
            << "  C[MC,MR] ~ " << C.Height() << " x " << C.Width() << endl;
        throw logic_error( msg.str() );
    }
    if( A1.ColAlignment() != C.ColAlignment() ||
        B1.ColAlignment() != C.RowAlignment() ||
        A1.ColAlignment() != A2.RowAlignment() ||
        B1.ColAlignment() != B2.RowAlignment() )
    {
        ostringstream msg;
        msg << "Misaligned LocalTriangularRank2K: " << endl
            << "  A1[MC,* ] ~ " << A1.ColAlignment() << endl
            << "  A2[* ,MC] ~ " << A2.RowAlignment() << endl
            << "  B1[MR,* ] ~ " << B1.ColAlignment() << endl
            << "  B2[* ,MR] ~ " << B2.RowAlignment() << endl
            << "  C[MC,MR] ~ " << C.ColAlignment() << " , " <<
                                  C.RowAlignment() << endl;
        throw logic_error( msg.str() );
    }
}

template<typename T>
void 
CheckInput
( const DistMatrix<T,MC,  Star>& A1, 
  const DistMatrix<T,MC,  Star>& A2,
  const DistMatrix<T,Star,MR  >& B1,
  const DistMatrix<T,Star,MR  >& B2,
  const DistMatrix<T,MC,  MR  >& C  )
{
    if( A1.GetGrid() != A2.GetGrid() || A2.GetGrid() != B1.GetGrid() ||
        B1.GetGrid() != B2.GetGrid() || B2.GetGrid() != C.GetGrid() )
        throw logic_error
        ( "A, B, and C must be distributed over the same grid." );
    if( A1.Height() != C.Height() || B1.Width() != C.Width() ||
        A1.Height() != A2.Height() || A1.Width() != A2.Width() ||
        B1.Height() != B2.Height() || B1.Width() != B2.Width() ||
        A1.Width() != B1.Height() )
    {
        ostringstream msg;
        msg << "Nonconformal LocalTriangularRank2K: " << endl
            << "  A1[MC,* ] ~ " << A1.Height() << " x "
                                << A1.Width()  << endl
            << "  A2[MC,* ] ~ " << A2.Height() << " x "
                                << A2.Width()  << endl
            << "  B1[* ,MR] ~ " << B1.Height() << " x "
                                << B1.Width()  << endl
            << "  B2[* ,MR] ~ " << B2.Height() << " x "
                                << B2.Width()  << endl
            << "  C[MC,MR] ~ " << C.Height() << " x " << C.Width() << endl;
        throw logic_error( msg.str() );
    }
    if( A1.ColAlignment() != C.ColAlignment() ||
        B1.RowAlignment() != C.RowAlignment() ||
        A1.ColAlignment() != A2.ColAlignment() ||
        B1.RowAlignment() != B2.RowAlignment() )
    {
        ostringstream msg;
        msg << "Misaligned LocalTriangularRank2K: " << endl
            << "  A1[MC,* ] ~ " << A1.ColAlignment() << endl
            << "  A2[MC,* ] ~ " << A2.ColAlignment() << endl
            << "  B1[* ,MR] ~ " << B1.RowAlignment() << endl
            << "  B2[* ,MR] ~ " << B2.RowAlignment() << endl
            << "  C[MC,MR] ~ " << C.ColAlignment() << " , " <<
                                  C.RowAlignment() << endl;
        throw logic_error( msg.str() );
    }
}

template<typename T>
void 
CheckInput
( Orientation orientationOfA1,
  Orientation orientationOfA2,
  Orientation orientationOfB2,
  const DistMatrix<T,Star,MC  >& A1, 
  const DistMatrix<T,Star,MC  >& A2,
  const DistMatrix<T,Star,MR  >& B1,
  const DistMatrix<T,MR,  Star>& B2,
  const DistMatrix<T,MC,  MR  >& C )
{
    if( orientationOfA1 == Normal )
        throw logic_error( "A1[* ,MC] must be (Conjugate)Transpose'd." );
    if( orientationOfA2 == Normal )
        throw logic_error( "A2[* ,MC] must be (Conjugate)Transpose'd." );
    if( orientationOfB2 == Normal )
        throw logic_error( "B2[MR,* ] must be (Conjugate)Transpose'd." );
    if( A1.GetGrid() != A2.GetGrid() || A2.GetGrid() != B1.GetGrid() ||
        B1.GetGrid() != B2.GetGrid() || B2.GetGrid() != C.GetGrid() )
        throw logic_error
        ( "A, B, and C must be distributed over the same grid." );
    if( A1.Width() != C.Height() || B1.Width() != C.Width() ||
        A1.Height() != A2.Height() || A1.Width() != A2.Width() ||
        B1.Width() != B2.Height() || B1.Height() != B2.Width() ||
        A1.Height() != B1.Height() )
    {
        ostringstream msg;
        msg << "Nonconformal LocalTriangularRank2K: " << endl
            << "  A1[* ,MC] ~ " << A1.Height() << " x "
                                << A1.Width()  << endl
            << "  A2[* ,MC] ~ " << A2.Height() << " x "
                                << A2.Width()  << endl
            << "  B1[* ,MR] ~ " << B1.Height() << " x "
                                << B1.Width()  << endl
            << "  B2[MR,* ] ~ " << B2.Height() << " x "
                                << B2.Width()  << endl
            << "  C[MC,MR] ~ " << C.Height() << " x " << C.Width() << endl;
        throw logic_error( msg.str() );
    }
    if( A1.RowAlignment() != C.ColAlignment() ||
        B1.RowAlignment() != C.RowAlignment() ||
        A1.RowAlignment() != A2.RowAlignment() ||
        B1.RowAlignment() != B2.ColAlignment() )
    {
        ostringstream msg;
        msg << "Misaligned LocalTriangularRank2K: " << endl
            << "  A1[* ,MC] ~ " << A1.RowAlignment() << endl
            << "  A2[* ,MC] ~ " << A2.RowAlignment() << endl
            << "  B1[* ,MR] ~ " << B1.RowAlignment() << endl
            << "  B2[MR,* ] ~ " << B2.ColAlignment() << endl
            << "  C[MC,MR] ~ " << C.ColAlignment() << " , " <<
                                  C.RowAlignment() << endl;
        throw logic_error( msg.str() );
    }
}

template<typename T>
void 
CheckInput
( Orientation orientationOfA1,
  Orientation orientationOfA2,
  Orientation orientationOfB1,
  const DistMatrix<T,Star,MC  >& A1, 
  const DistMatrix<T,Star,MC  >& A2,
  const DistMatrix<T,MR,  Star>& B1,
  const DistMatrix<T,Star,MR  >& B2,
  const DistMatrix<T,MC,  MR  >& C )
{
    if( orientationOfA1 == Normal )
        throw logic_error( "A1[* ,MC] must be (Conjugate)Transpose'd." );
    if( orientationOfA2 == Normal )
        throw logic_error( "A2[* ,MC] must be (Conjugate)Transpose'd." );
    if( orientationOfB1 == Normal )
        throw logic_error( "B1[MR,* ] must be (Conjugate)Transpose'd." );
    if( A1.GetGrid() != A2.GetGrid() || A2.GetGrid() != B1.GetGrid() ||
        B1.GetGrid() != B2.GetGrid() || B2.GetGrid() != C.GetGrid() )
        throw logic_error
        ( "A, B, and C must be distributed over the same grid." );
    if( A1.Width() != C.Height() || B1.Height() != C.Width() ||
        A1.Height() != A2.Height() || A1.Width() != A2.Width() ||
        B1.Height() != B2.Width() || B1.Width() != B2.Height() ||
        A1.Height() != B1.Width() )
    {
        ostringstream msg;
        msg << "Nonconformal LocalTriangularRank2K: " << endl
            << "  A1[* ,MC] ~ " << A1.Height() << " x "
                                << A1.Width()  << endl
            << "  A2[* ,MC] ~ " << A2.Height() << " x "
                                << A2.Width()  << endl
            << "  B1[MR,* ] ~ " << B1.Height() << " x "
                                << B1.Width()  << endl
            << "  B2[* ,MR] ~ " << B2.Height() << " x "
                                << B2.Width()  << endl
            << "  C[MC,MR] ~ " << C.Height() << " x " << C.Width() << endl;
        throw logic_error( msg.str() );
    }
    if( A1.RowAlignment() != C.ColAlignment() ||
        B1.ColAlignment() != C.RowAlignment() ||
        A1.RowAlignment() != A2.RowAlignment() ||
        B1.ColAlignment() != B2.RowAlignment() )
    {
        ostringstream msg;
        msg << "Misaligned LocalTriangularRank2K: " << endl
            << "  A1[* ,MC] ~ " << A1.RowAlignment() << endl
            << "  A2[* ,MC] ~ " << A2.RowAlignment() << endl
            << "  B1[MR,* ] ~ " << B1.ColAlignment() << endl
            << "  B2[* ,MR] ~ " << B2.RowAlignment() << endl
            << "  C[MC,MR] ~ " << C.ColAlignment() << " , " <<
                                  C.RowAlignment() << endl;
        throw logic_error( msg.str() );
    }
}

template<typename T>
void 
CheckInput
( Orientation orientationOfA1,
  const DistMatrix<T,Star,MC  >& A1, 
  const DistMatrix<T,MC,  Star>& A2,
  const DistMatrix<T,Star,MR  >& B1,
  const DistMatrix<T,Star,MR  >& B2,
  const DistMatrix<T,MC,  MR  >& C )
{
    if( orientationOfA1 == Normal )
        throw logic_error( "A1[* ,MC] must be (Conjugate)Transpose'd." );
    if( A1.GetGrid() != A2.GetGrid() || A2.GetGrid() != B1.GetGrid() ||
        B1.GetGrid() != B2.GetGrid() || B2.GetGrid() != C.GetGrid() )
        throw logic_error
        ( "A, B, and C must be distributed over the same grid." );
    if( A1.Width() != C.Height() || B1.Width() != C.Width() ||
        A1.Width() != A2.Height() || A1.Height() != A2.Width() ||
        B1.Height() != B2.Height() || B1.Width() != B2.Width() ||
        A1.Height() != B1.Height() )
    {
        ostringstream msg;
        msg << "Nonconformal LocalTriangularRank2K: " << endl
            << "  A1[* ,MC] ~ " << A1.Height() << " x "
                                << A1.Width()  << endl
            << "  A2[MC,* ] ~ " << A2.Height() << " x "
                                << A2.Width()  << endl
            << "  B1[* ,MR] ~ " << B1.Height() << " x "
                                << B1.Width()  << endl
            << "  B2[* ,MR] ~ " << B2.Height() << " x "
                                << B2.Width()  << endl
            << "  C[MC,MR] ~ " << C.Height() << " x " << C.Width() << endl;
        throw logic_error( msg.str() );
    }
    if( A1.RowAlignment() != C.ColAlignment() ||
        B1.RowAlignment() != C.RowAlignment() ||
        A1.RowAlignment() != A2.ColAlignment() ||
        B1.RowAlignment() != B2.RowAlignment() )
    {
        ostringstream msg;
        msg << "Misaligned LocalTriangularRank2K: " << endl
            << "  A1[* ,MC] ~ " << A1.RowAlignment() << endl
            << "  A2[MC,* ] ~ " << A2.ColAlignment() << endl
            << "  B1[* ,MR] ~ " << B1.RowAlignment() << endl
            << "  B2[* ,MR] ~ " << B2.RowAlignment() << endl
            << "  C[MC,MR] ~ " << C.ColAlignment() << " , " <<
                                  C.RowAlignment() << endl;
        throw logic_error( msg.str() );
    }
}

template<typename T>
void 
CheckInput
( Orientation orientationOfA2,
  const DistMatrix<T,MC,  Star>& A1,
  const DistMatrix<T,Star,MC  >& A2, 
  const DistMatrix<T,Star,MR  >& B1,
  const DistMatrix<T,Star,MR  >& B2,
  const DistMatrix<T,MC,  MR  >& C )
{
    if( orientationOfA2 == Normal )
        throw logic_error( "A2[* ,MC] must be (Conjugate)Transpose'd." );
    if( A1.GetGrid() != A2.GetGrid() || A2.GetGrid() != B1.GetGrid() ||
        B1.GetGrid() != B2.GetGrid() || B2.GetGrid() != C.GetGrid() )
        throw logic_error
        ( "A, B, and C must be distributed over the same grid." );
    if( A1.Height() != C.Height() || B1.Width() != C.Width() ||
        A1.Height() != A2.Width() || A1.Width() != A2.Height() ||
        B1.Height() != B2.Height() || B1.Width() != B2.Width() ||
        A1.Width() != B1.Height() )
    {
        ostringstream msg;
        msg << "Nonconformal LocalTriangularRank2K: " << endl
            << "  A1[MC,* ] ~ " << A1.Height() << " x "
                                << A1.Width()  << endl
            << "  A2[* ,MC] ~ " << A2.Height() << " x "
                                << A2.Width()  << endl
            << "  B1[* ,MR] ~ " << B1.Height() << " x "
                                << B1.Width()  << endl
            << "  B2[* ,MR] ~ " << B2.Height() << " x "
                                << B2.Width()  << endl
            << "  C[MC,MR] ~ " << C.Height() << " x " << C.Width() << endl;
        throw logic_error( msg.str() );
    }
    if( A1.ColAlignment() != C.ColAlignment() ||
        B1.RowAlignment() != C.RowAlignment() ||
        A1.ColAlignment() != A2.RowAlignment() ||
        B1.RowAlignment() != B2.RowAlignment() )
    {
        ostringstream msg;
        msg << "Misaligned LocalTriangularRank2K: " << endl
            << "  A1[MC,* ] ~ " << A1.ColAlignment() << endl
            << "  A2[* ,MC] ~ " << A2.RowAlignment() << endl
            << "  B1[* ,MR] ~ " << B1.RowAlignment() << endl
            << "  B2[* ,MR] ~ " << B2.RowAlignment() << endl
            << "  C[MC,MR] ~ " << C.ColAlignment() << " , " <<
                                  C.RowAlignment() << endl;
        throw logic_error( msg.str() );
    }
}

template<typename T>
void 
CheckInput
( Orientation orientationOfA1,
  Orientation orientationOfA2,
  const DistMatrix<T,Star,MC>& A1, 
  const DistMatrix<T,Star,MC>& A2,
  const DistMatrix<T,Star,MR>& B1,
  const DistMatrix<T,Star,MR>& B2,
  const DistMatrix<T,MC,  MR>& C )
{
    if( orientationOfA1 == Normal )
        throw logic_error( "A1[* ,MC] must be (Conjugate)Transpose'd." );
    if( orientationOfA2 == Normal )
        throw logic_error( "A2[* ,MC] must be (Conjugate)Transpose'd." );
    if( A1.GetGrid() != A2.GetGrid() || A2.GetGrid() != B1.GetGrid() ||
        B1.GetGrid() != B2.GetGrid() || B2.GetGrid() != C.GetGrid() )
        throw logic_error
        ( "A, B, and C must be distributed over the same grid." );
    if( A1.Width() != C.Height() || B1.Width() != C.Width() ||
        A1.Height() != A2.Height() || A1.Width() != A2.Width() ||
        B1.Height() != B2.Height() || B1.Width() != B2.Width() ||
        A1.Height() != B1.Height() )
    {
        ostringstream msg;
        msg << "Nonconformal LocalTriangularRank2K: " << endl
            << "  A1[* ,MC] ~ " << A1.Height() << " x "
                                << A1.Width()  << endl
            << "  A2[* ,MC] ~ " << A2.Height() << " x "
                                << A2.Width()  << endl
            << "  B1[* ,MR] ~ " << B1.Height() << " x "
                                << B1.Width()  << endl
            << "  B2[* ,MR] ~ " << B2.Height() << " x "
                                << B2.Width()  << endl
            << "  C[MC,MR] ~ " << C.Height() << " x " << C.Width() << endl;
        throw logic_error( msg.str() );
    }
    if( A1.RowAlignment() != C.ColAlignment() ||
        B1.RowAlignment() != C.RowAlignment() ||
        A1.RowAlignment() != A2.RowAlignment() ||
        B1.RowAlignment() != B2.RowAlignment() )
    {
        ostringstream msg;
        msg << "Misaligned LocalTriangularRank2K: " << endl
            << "  A1[* ,MC] ~ " << A1.RowAlignment() << endl
            << "  A2[* ,MC] ~ " << A2.RowAlignment() << endl
            << "  B1[* ,MR] ~ " << B1.RowAlignment() << endl
            << "  B2[* ,MR] ~ " << B2.RowAlignment() << endl
            << "  C[MC,MR] ~ " << C.ColAlignment() << " , " <<
                                  C.RowAlignment() << endl;
        throw logic_error( msg.str() );
    }
}
#endif 

template<typename T>
void
LocalTriangularRank2KKernel
( Shape shape,
  Orientation orientationOfB1,
  Orientation orientationOfB2,
  T alpha, const DistMatrix<T,MC,Star>& A1,
           const DistMatrix<T,MC,Star>& A2,
           const DistMatrix<T,MR,Star>& B1,
           const DistMatrix<T,MR,Star>& B2,
  T beta,        DistMatrix<T,MC,MR  >& C )
{
#ifndef RELEASE
    PushCallStack("LocalTriangularRank2KKernel");
    CheckInput( orientationOfB1, orientationOfB2, A1, A2, B1, B2, C );
#endif
    const Grid& g = C.GetGrid();

    DistMatrix<T,MC,Star> A1T(g),  A2T(g),
                          A1B(g),  A2B(g);

    DistMatrix<T,MR,Star> B1T(g),  B2T(g),
                          B1B(g),  B2B(g);

    DistMatrix<T,MC,MR> CTL(g), CTR(g),
                        CBL(g), CBR(g);

    DistMatrix<T,MC,MR> DTL(g), DBR(g);

    const unsigned half = C.Height()/2;

    blas::Scal( beta, C );

    LockedPartitionDown
    ( A1, A1T,
          A1B, half );
    LockedPartitionDown
    ( A2, A2T,
          A2B, half );

    LockedPartitionDown
    ( B1, B1T,
          B1B, half );
    LockedPartitionDown
    ( B2, B2T,
          B2B, half );

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
        ( Normal, orientationOfB2, alpha, A1B, B2T, (T)1, CBL );

        blas::internal::LocalGemm
        ( Normal, orientationOfB1, alpha, A2B, B1T, (T)1, CBL );
    }
    else
    {
        blas::internal::LocalGemm
        ( Normal, orientationOfB2, alpha, A1T, B2B, (T)1, CTR );

        blas::internal::LocalGemm
        ( Normal, orientationOfB1, alpha, A2T, B1B, (T)1, CTR );
    }

    blas::internal::LocalGemm
    ( Normal, orientationOfB2, alpha, A1T, B2T, (T)0, DTL );

    blas::internal::LocalGemm
    ( Normal, orientationOfB1, alpha, A2T, B1T, (T)1, DTL );

    DTL.MakeTrapezoidal( Left, shape );
    blas::Axpy( (T)1, DTL, CTL );

    blas::internal::LocalGemm
    ( Normal, orientationOfB2, alpha, A1B, B2B, (T)0, DBR );

    blas::internal::LocalGemm
    ( Normal, orientationOfB1, alpha, A2B, B1B, (T)1, DBR );

    DBR.MakeTrapezoidal( Left, shape );
    blas::Axpy( (T)1, DBR, CBR );
    //------------------------------------------------------------------------//
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
LocalTriangularRank2KKernel
( Shape shape,
  Orientation orientationOfA1,
  Orientation orientationOfB1,
  Orientation orientationOfB2,
  T alpha, const DistMatrix<T,Star,MC  >& A1,
           const DistMatrix<T,MC,  Star>& A2,
           const DistMatrix<T,MR,  Star>& B1,
           const DistMatrix<T,MR,  Star>& B2,
  T beta,        DistMatrix<T,MC,  MR  >& C )
{
#ifndef RELEASE
    PushCallStack("LocalTriangularRank2KKernel");
    CheckInput
    ( orientationOfA1, orientationOfB1, orientationOfB2, A1, A2, B1, B2, C );
#endif
    const Grid& g = C.GetGrid();

    DistMatrix<T,Star,MC> A1L(g), A1R(g);

    DistMatrix<T,MC,Star> A2T(g),
                          A2B(g);

    DistMatrix<T,MR,Star> B1T(g), 
                          B1B(g);

    DistMatrix<T,MR,Star> B2T(g), 
                          B2B(g);

    DistMatrix<T,MC,MR> CTL(g), CTR(g),
                        CBL(g), CBR(g);

    DistMatrix<T,MC,MR> DTL(g), DBR(g);

    const unsigned half = C.Height()/2;

    blas::Scal( beta, C );

    LockedPartitionRight
    ( A1, A1L, A1R, half );
    LockedPartitionDown
    ( A2, A2T,
          A2B, half );

    LockedPartitionDown
    ( B1, B1T,
          B1B, half );
    LockedPartitionDown
    ( B2, B2T,
          B2B, half );

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
        ( orientationOfA1, orientationOfB2, alpha, A1R, B2T, (T)1, CBL );

        blas::internal::LocalGemm
        ( Normal, orientationOfB1, alpha, A2B, B1T, (T)1, CBL );
    }
    else
    {
        blas::internal::LocalGemm
        ( orientationOfA1, orientationOfB2, alpha, A1L, B2B, (T)1, CTR );

        blas::internal::LocalGemm
        ( Normal, orientationOfB1, alpha, A2T, B1B, (T)1, CTR );
    }

    blas::internal::LocalGemm
    ( orientationOfA1, orientationOfB2, alpha, A1L, B2T, (T)0, DTL );

    blas::internal::LocalGemm
    ( Normal, orientationOfB1, alpha, A2T, B1T, (T)1, DTL );

    DTL.MakeTrapezoidal( Left, shape );
    blas::Axpy( (T)1, DTL, CTL );

    blas::internal::LocalGemm
    ( orientationOfA1, orientationOfB2, alpha, A1R, B2B, (T)0, DBR );

    blas::internal::LocalGemm
    ( Normal, orientationOfB1, alpha, A2B, B1B, (T)1, DBR );

    DBR.MakeTrapezoidal( Left, shape );
    blas::Axpy( (T)1, DBR, CBR );
    //------------------------------------------------------------------------//
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
LocalTriangularRank2KKernel
( Shape shape,
  Orientation orientationOfA2,
  Orientation orientationOfB1,
  Orientation orientationOfB2,
  T alpha, const DistMatrix<T,MC,  Star>& A1,
           const DistMatrix<T,Star,MC  >& A2,
           const DistMatrix<T,MR,  Star>& B1,
           const DistMatrix<T,MR,  Star>& B2,
  T beta,        DistMatrix<T,MC,  MR  >& C )
{
#ifndef RELEASE
    PushCallStack("LocalTriangularRank2KKernel");
    CheckInput
    ( orientationOfA2, orientationOfB1, orientationOfB2, A1, A2, B1, B2, C );
#endif
    const Grid& g = C.GetGrid();

    DistMatrix<T,MC,Star> A1T(g),
                          A1B(g);

    DistMatrix<T,Star,MC> A2L(g), A2R(g);

    DistMatrix<T,MR,Star> B1T(g), 
                          B1B(g);

    DistMatrix<T,MR,Star> B2T(g), 
                          B2B(g);

    DistMatrix<T,MC,MR> CTL(g), CTR(g),
                        CBL(g), CBR(g);

    DistMatrix<T,MC,MR> DTL(g), DBR(g);

    const unsigned half = C.Height()/2;

    blas::Scal( beta, C );

    LockedPartitionDown
    ( A1, A1T,
          A1B, half );
    LockedPartitionRight
    ( A2, A2L, A2R, half );

    LockedPartitionDown
    ( B1, B1T,
          B1B, half );
    LockedPartitionDown
    ( B2, B2T,
          B2B, half );

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
        ( Normal, orientationOfB2, alpha, A1B, B2T, (T)1, CBL );

        blas::internal::LocalGemm
        ( orientationOfA2, orientationOfB1, alpha, A2R, B1T, (T)1, CBL );
    }
    else
    {
        blas::internal::LocalGemm
        ( Normal, orientationOfB2, alpha, A1T, B2B, (T)1, CTR );

        blas::internal::LocalGemm
        ( orientationOfA2, orientationOfB1, alpha, A2L, B1B, (T)1, CTR );
    }

    blas::internal::LocalGemm
    ( Normal, orientationOfB2, alpha, A1T, B2T, (T)0, DTL );

    blas::internal::LocalGemm
    ( orientationOfA2, orientationOfB1, alpha, A2L, B1T, (T)1, DTL );

    DTL.MakeTrapezoidal( Left, shape );
    blas::Axpy( (T)1, DTL, CTL );

    blas::internal::LocalGemm
    ( Normal, orientationOfB2, alpha, A1B, B2B, (T)0, DBR );

    blas::internal::LocalGemm
    ( orientationOfA2, orientationOfB1, alpha, A2R, B1B, (T)1, DBR );

    DBR.MakeTrapezoidal( Left, shape );
    blas::Axpy( (T)1, DBR, CBR );
    //------------------------------------------------------------------------//
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
LocalTriangularRank2KKernel
( Shape shape,
  Orientation orientationOfB2,
  T alpha, const DistMatrix<T,MC,  Star>& A1,
           const DistMatrix<T,MC,  Star>& A2,
           const DistMatrix<T,Star,MR  >& B1,
           const DistMatrix<T,MR,  Star>& B2,
  T beta,        DistMatrix<T,MC,  MR  >& C )
{
#ifndef RELEASE
    PushCallStack("LocalTriangularRank2KKernel");
    CheckInput( orientationOfB2, A1, A2, B1, B2, C );
#endif
    const Grid& g = C.GetGrid();

    DistMatrix<T,MC,Star> A1T(g),  A2T(g),
                          A1B(g),  A2B(g);

    DistMatrix<T,Star,MR> B1L(g), B1R(g);

    DistMatrix<T,MR,Star> B2T(g), 
                          B2B(g);

    DistMatrix<T,MC,MR> CTL(g), CTR(g),
                        CBL(g), CBR(g);

    DistMatrix<T,MC,MR> DTL(g), DBR(g);

    const unsigned half = C.Height()/2;

    blas::Scal( beta, C );

    LockedPartitionDown
    ( A1, A1T,
          A1B, half );
    LockedPartitionDown
    ( A2, A2T,
          A2B, half );

    LockedPartitionRight( B1, B1L, B1R, half );
    LockedPartitionDown
    ( B2, B2T,
          B2B, half );

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
        ( Normal, orientationOfB2, alpha, A1B, B2T, (T)1, CBL );

        blas::internal::LocalGemm
        ( Normal, Normal, alpha, A2B, B1L, (T)1, CBL );
    }
    else
    {
        blas::internal::LocalGemm
        ( Normal, orientationOfB2, alpha, A1T, B2B, (T)1, CTR );

        blas::internal::LocalGemm
        ( Normal, Normal, alpha, A2T, B1R, (T)1, CTR );
    }

    blas::internal::LocalGemm
    ( Normal, orientationOfB2, alpha, A1T, B2T, (T)0, DTL );

    blas::internal::LocalGemm
    ( Normal, Normal, alpha, A2T, B1L, (T)1, DTL );

    DTL.MakeTrapezoidal( Left, shape );
    blas::Axpy( (T)1, DTL, CTL );

    blas::internal::LocalGemm
    ( Normal, orientationOfB2, alpha, A1B, B2B, (T)0, DBR );

    blas::internal::LocalGemm
    ( Normal, Normal, alpha, A2B, B1R, (T)1, DBR );

    DBR.MakeTrapezoidal( Left, shape );
    blas::Axpy( (T)1, DBR, CBR );
    //------------------------------------------------------------------------//
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
LocalTriangularRank2KKernel
( Shape shape,
  Orientation orientationOfB1,
  T alpha, const DistMatrix<T,MC,  Star>& A1,
           const DistMatrix<T,MC,  Star>& A2,
           const DistMatrix<T,MR,  Star>& B1,
           const DistMatrix<T,Star,MR  >& B2,
  T beta,        DistMatrix<T,MC,  MR  >& C )
{
#ifndef RELEASE
    PushCallStack("LocalTriangularRank2KKernel");
    CheckInput( orientationOfB1, A1, A2, B1, B2, C );
#endif
    const Grid& g = C.GetGrid();

    DistMatrix<T,MC,Star> A1T(g),  A2T(g),
                          A1B(g),  A2B(g);

    DistMatrix<T,MR,Star> B1T(g), 
                          B1B(g);

    DistMatrix<T,Star,MR> B2L(g), B2R(g);

    DistMatrix<T,MC,MR> CTL(g), CTR(g),
                        CBL(g), CBR(g);

    DistMatrix<T,MC,MR> DTL(g), DBR(g);

    const unsigned half = C.Height()/2;

    blas::Scal( beta, C );

    LockedPartitionDown
    ( A1, A1T,
          A1B, half );
    LockedPartitionDown
    ( A2, A2T,
          A2B, half );

    LockedPartitionDown
    ( B1, B1T,
          B1B, half );

    LockedPartitionRight( B2, B2L, B2R, half );

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
        ( Normal, Normal, alpha, A1B, B2L, (T)1, CBL );

        blas::internal::LocalGemm
        ( Normal, orientationOfB1, alpha, A2B, B1T, (T)1, CBL );
    }
    else
    {
        blas::internal::LocalGemm
        ( Normal, Normal, alpha, A1T, B2R, (T)1, CTR );

        blas::internal::LocalGemm
        ( Normal, orientationOfB1, alpha, A2T, B1B, (T)1, CTR );
    }

    blas::internal::LocalGemm
    ( Normal, Normal, alpha, A1T, B2L, (T)0, DTL );

    blas::internal::LocalGemm
    ( Normal, orientationOfB1, alpha, A2T, B1T, (T)1, DTL );

    DTL.MakeTrapezoidal( Left, shape );
    blas::Axpy( (T)1, DTL, CTL );

    blas::internal::LocalGemm
    ( Normal, Normal, alpha, A1B, B2R, (T)0, DBR );

    blas::internal::LocalGemm
    ( Normal, orientationOfB1, alpha, A2B, B1B, (T)1, DBR );

    DBR.MakeTrapezoidal( Left, shape );
    blas::Axpy( (T)1, DBR, CBR );
    //------------------------------------------------------------------------//
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
LocalTriangularRank2KKernel
( Shape shape,
  Orientation orientationOfA1,
  Orientation orientationOfA2,
  Orientation orientationOfB1,
  Orientation orientationOfB2,
  T alpha, const DistMatrix<T,Star,MC  >& A1,
           const DistMatrix<T,Star,MC  >& A2,
           const DistMatrix<T,MR,  Star>& B1,
           const DistMatrix<T,MR,  Star>& B2,
  T beta,        DistMatrix<T,MC,  MR  >& C )
{
#ifndef RELEASE
    PushCallStack("LocalTriangularRank2KKernel");
    CheckInput
    ( orientationOfA1, orientationOfA2, orientationOfB1, orientationOfB2, 
      A1, A2, B1, B2, C );
#endif
    const Grid& g = C.GetGrid();

    DistMatrix<T,Star,MC> A1L(g), A1R(g),
                          A2L(g), A2R(g);

    DistMatrix<T,MR,Star> B1T(g),  B2T(g),
                          B1B(g),  B2B(g);

    DistMatrix<T,MC,MR> CTL(g), CTR(g),
                        CBL(g), CBR(g);

    DistMatrix<T,MC,MR> DTL(g), DBR(g);

    const unsigned half = C.Height()/2;

    blas::Scal( beta, C );

    LockedPartitionRight( A1, A1L, A1R, half );
    LockedPartitionRight( A2, A2L, A2R, half );

    LockedPartitionDown
    ( B1, B1T,
          B1B, half );

    LockedPartitionDown
    ( B2, B2T,
          B2B, half );

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
        ( orientationOfA1, orientationOfB2, alpha, A1R, B2T, (T)1, CBL );

        blas::internal::LocalGemm
        ( orientationOfA2, orientationOfB1, alpha, A2R, B1T, (T)1, CBL );
    }
    else
    {
        blas::internal::LocalGemm
        ( orientationOfA1, orientationOfB2, alpha, A1L, B2B, (T)1, CTR );

        blas::internal::LocalGemm
        ( orientationOfA2, orientationOfB1, alpha, A2L, B1B, (T)1, CTR );
    }

    blas::internal::LocalGemm
    ( orientationOfA1, orientationOfB2, alpha, A1L, B2T, (T)0, DTL );

    blas::internal::LocalGemm
    ( orientationOfA2, orientationOfB1, alpha, A2L, B1T, (T)1, DTL );

    DTL.MakeTrapezoidal( Left, shape );
    blas::Axpy( (T)1, DTL, CTL );

    blas::internal::LocalGemm
    ( orientationOfA1, orientationOfB2, alpha, A1R, B2B, (T)0, DBR );

    blas::internal::LocalGemm
    ( orientationOfA2, orientationOfB1, alpha, A2R, B1B, (T)1, DBR );

    DBR.MakeTrapezoidal( Left, shape );
    blas::Axpy( (T)1, DBR, CBR );
    //------------------------------------------------------------------------//
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
LocalTriangularRank2KKernel
( Shape shape,
  Orientation orientationOfA1,
  Orientation orientationOfB2,
  T alpha, const DistMatrix<T,Star,MC  >& A1,
           const DistMatrix<T,MC,  Star>& A2,
           const DistMatrix<T,Star,MR  >& B1,
           const DistMatrix<T,MR,  Star>& B2,
  T beta,        DistMatrix<T,MC,  MR  >& C )
{
#ifndef RELEASE
    PushCallStack("LocalTriangularRank2KKernel");
    CheckInput( orientationOfA1, orientationOfB2, A1, A2, B1, B2, C );
#endif
    const Grid& g = C.GetGrid();

    DistMatrix<T,Star,MC> A1L(g), A1R(g);

    DistMatrix<T,MC,Star> A2T(g),
                          A2B(g);

    DistMatrix<T,Star,MR> B1L(g), B1R(g);

    DistMatrix<T,MR,Star> B2T(g),
                          B2B(g);

    DistMatrix<T,MC,MR> CTL(g), CTR(g),
                        CBL(g), CBR(g);

    DistMatrix<T,MC,MR> DTL(g), DBR(g);

    const unsigned half = C.Height()/2;

    blas::Scal( beta, C );

    LockedPartitionRight( A1, A1L, A1R, half );
    LockedPartitionDown
    ( A2, A2T,
          A2B, half );

    LockedPartitionRight( B1, B1L, B1R, half );
    LockedPartitionDown
    ( B2, B2T,
          B2B, half );

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
        ( orientationOfA1, orientationOfB2, alpha, A1R, B2T, (T)1, CBL );

        blas::internal::LocalGemm
        ( Normal, Normal, alpha, A2B, B1L, (T)1, CBL );
    }
    else
    {
        blas::internal::LocalGemm
        ( orientationOfA1, orientationOfB2, alpha, A1L, B2B, (T)1, CTR );

        blas::internal::LocalGemm
        ( Normal, Normal, alpha, A2T, B1R, (T)1, CTR );
    }

    blas::internal::LocalGemm
    ( orientationOfA1, orientationOfB2, alpha, A1L, B2T, (T)0, DTL );

    blas::internal::LocalGemm
    ( Normal, Normal, alpha, A2T, B1L, (T)1, DTL );

    DTL.MakeTrapezoidal( Left, shape );
    blas::Axpy( (T)1, DTL, CTL );

    blas::internal::LocalGemm
    ( orientationOfA1, orientationOfB2, alpha, A1R, B2B, (T)0, DBR );

    blas::internal::LocalGemm
    ( Normal, Normal, alpha, A2B, B1R, (T)1, DBR );

    DBR.MakeTrapezoidal( Left, shape );
    blas::Axpy( (T)1, DBR, CBR );
    //------------------------------------------------------------------------//
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
LocalTriangularRank2KKernel
( Shape shape,
  Orientation orientationOfA1,
  Orientation orientationOfB1,
  T alpha, const DistMatrix<T,Star,MC  >& A1,
           const DistMatrix<T,MC,  Star>& A2,
           const DistMatrix<T,MR,  Star>& B1,
           const DistMatrix<T,Star,MR  >& B2,
  T beta,        DistMatrix<T,MC,  MR  >& C )
{
#ifndef RELEASE
    PushCallStack("LocalTriangularRank2KKernel");
    CheckInput( orientationOfA1, orientationOfB1, A1, A2, B1, B2, C );
#endif
    const Grid& g = C.GetGrid();

    DistMatrix<T,Star,MC> A1L(g), A1R(g);

    DistMatrix<T,MC,Star> A2T(g),
                          A2B(g);

    DistMatrix<T,MR,Star> B1T(g),
                          B1B(g);

    DistMatrix<T,Star,MR> B2L(g), B2R(g);

    DistMatrix<T,MC,MR> CTL(g), CTR(g),
                        CBL(g), CBR(g);

    DistMatrix<T,MC,MR> DTL(g), DBR(g);

    const unsigned half = C.Height()/2;

    blas::Scal( beta, C );

    LockedPartitionRight( A1, A1L, A1R, half );
    LockedPartitionDown
    ( A2, A2T,
          A2B, half );

    LockedPartitionDown
    ( B1, B1T,
          B1B, half );
    LockedPartitionRight( B2, B2L, B2R, half );

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
        ( orientationOfA1, Normal, alpha, A1R, B2L, (T)1, CBL );

        blas::internal::LocalGemm
        ( Normal, orientationOfB1, alpha, A2B, B1T, (T)1, CBL );
    }
    else
    {
        blas::internal::LocalGemm
        ( orientationOfA1, Normal, alpha, A1L, B2R, (T)1, CTR );

        blas::internal::LocalGemm
        ( Normal, orientationOfB1, alpha, A2T, B1B, (T)1, CTR );
    }

    blas::internal::LocalGemm
    ( orientationOfA1, Normal, alpha, A1L, B2L, (T)0, DTL );

    blas::internal::LocalGemm
    ( Normal, orientationOfB1, alpha, A2T, B1T, (T)1, DTL );

    DTL.MakeTrapezoidal( Left, shape );
    blas::Axpy( (T)1, DTL, CTL );

    blas::internal::LocalGemm
    ( orientationOfA1, Normal, alpha, A1R, B2R, (T)0, DBR );

    blas::internal::LocalGemm
    ( Normal, orientationOfB1, alpha, A2B, B1B, (T)1, DBR );

    DBR.MakeTrapezoidal( Left, shape );
    blas::Axpy( (T)1, DBR, CBR );
    //------------------------------------------------------------------------//
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
LocalTriangularRank2KKernel
( Shape shape,
  Orientation orientationOfA2,
  Orientation orientationOfB2,
  T alpha, const DistMatrix<T,MC,  Star>& A1,
           const DistMatrix<T,Star,MC  >& A2,
           const DistMatrix<T,Star,MR  >& B1,
           const DistMatrix<T,MR,  Star>& B2,
  T beta,        DistMatrix<T,MC,  MR  >& C )
{
#ifndef RELEASE
    PushCallStack("LocalTriangularRank2KKernel");
    CheckInput( orientationOfA2, orientationOfB2, A1, A2, B1, B2, C );
#endif
    const Grid& g = C.GetGrid();

    DistMatrix<T,MC,Star> A1T(g),
                          A1B(g);
    
    DistMatrix<T,Star,MC> A2L(g), A2R(g);

    DistMatrix<T,Star,MR> B1L(g), B1R(g);

    DistMatrix<T,MR,Star> B2T(g),
                          B2B(g);

    DistMatrix<T,MC,MR> CTL(g), CTR(g),
                        CBL(g), CBR(g);

    DistMatrix<T,MC,MR> DTL(g), DBR(g);

    const unsigned half = C.Height()/2;

    blas::Scal( beta, C );

    LockedPartitionDown
    ( A1, A1T,
          A1B, half );
    LockedPartitionRight( A2, A2L, A2R, half );

    LockedPartitionRight( B1, B1L, B1R, half );
    LockedPartitionDown
    ( B2, B2T,
          B2B, half );

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
        ( Normal, orientationOfB2, alpha, A1B, B2T, (T)1, CBL );

        blas::internal::LocalGemm
        ( orientationOfA2, Normal, alpha, A2R, B1L, (T)1, CBL );
    }
    else
    {
        blas::internal::LocalGemm
        ( Normal, orientationOfB2, alpha, A1T, B2B, (T)1, CTR );

        blas::internal::LocalGemm
        ( orientationOfA2, Normal, alpha, A2L, B1R, (T)1, CTR );
    }

    blas::internal::LocalGemm
    ( Normal, orientationOfB2, alpha, A1T, B2T, (T)0, DTL );

    blas::internal::LocalGemm
    ( orientationOfA2, Normal, alpha, A2L, B1L, (T)1, DTL );

    DTL.MakeTrapezoidal( Left, shape );
    blas::Axpy( (T)1, DTL, CTL );

    blas::internal::LocalGemm
    ( Normal, orientationOfB2, alpha, A1B, B2B, (T)0, DBR );

    blas::internal::LocalGemm
    ( orientationOfA2, Normal, alpha, A2R, B1R, (T)1, DBR );

    DBR.MakeTrapezoidal( Left, shape );
    blas::Axpy( (T)1, DBR, CBR );
    //------------------------------------------------------------------------//
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
LocalTriangularRank2KKernel
( Shape shape,
  Orientation orientationOfA2,
  Orientation orientationOfB1,
  T alpha, const DistMatrix<T,MC,  Star>& A1,
           const DistMatrix<T,Star,MC  >& A2,
           const DistMatrix<T,MR,  Star>& B1,
           const DistMatrix<T,Star,MR  >& B2,
  T beta,        DistMatrix<T,MC,  MR  >& C )
{
#ifndef RELEASE
    PushCallStack("LocalTriangularRank2KKernel");
    CheckInput( orientationOfA2, orientationOfB1, A1, A2, B1, B2, C );
#endif
    const Grid& g = C.GetGrid();

    DistMatrix<T,MC,Star> A1T(g),
                          A1B(g);
    
    DistMatrix<T,Star,MC> A2L(g), A2R(g);

    DistMatrix<T,MR,Star> B1T(g),
                          B1B(g);

    DistMatrix<T,Star,MR> B2L(g), B2R(g);

    DistMatrix<T,MC,MR> CTL(g), CTR(g),
                        CBL(g), CBR(g);

    DistMatrix<T,MC,MR> DTL(g), DBR(g);

    const unsigned half = C.Height()/2;

    blas::Scal( beta, C );

    LockedPartitionDown
    ( A1, A1T,
          A1B, half );
    LockedPartitionRight( A2, A2L, A2R, half );

    LockedPartitionDown
    ( B1, B1T,
          B1B, half );
    LockedPartitionRight( B2, B2L, B2R, half );

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
        ( Normal, Normal, alpha, A1B, B2L, (T)1, CBL );

        blas::internal::LocalGemm
        ( orientationOfA2, orientationOfB1, alpha, A2R, B1T, (T)1, CBL );
    }
    else
    {
        blas::internal::LocalGemm
        ( Normal, Normal, alpha, A1T, B2R, (T)1, CTR );

        blas::internal::LocalGemm
        ( orientationOfA2, orientationOfB1, alpha, A2L, B1B, (T)1, CTR );
    }

    blas::internal::LocalGemm
    ( Normal, Normal, alpha, A1T, B2L, (T)0, DTL );

    blas::internal::LocalGemm
    ( orientationOfA2, orientationOfB1, alpha, A2L, B1T, (T)1, DTL );

    DTL.MakeTrapezoidal( Left, shape );
    blas::Axpy( (T)1, DTL, CTL );

    blas::internal::LocalGemm
    ( Normal, Normal, alpha, A1B, B2R, (T)0, DBR );

    blas::internal::LocalGemm
    ( orientationOfA2, orientationOfB1, alpha, A2R, B1B, (T)1, DBR );

    DBR.MakeTrapezoidal( Left, shape );
    blas::Axpy( (T)1, DBR, CBR );
    //------------------------------------------------------------------------//
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
LocalTriangularRank2KKernel
( Shape shape,
  T alpha, const DistMatrix<T,MC,  Star>& A1,
           const DistMatrix<T,MC,  Star>& A2,
           const DistMatrix<T,Star,MR  >& B1,
           const DistMatrix<T,Star,MR  >& B2,
  T beta,        DistMatrix<T,MC,  MR  >& C )
{
#ifndef RELEASE
    PushCallStack("LocalTriangularRank2KKernel");
    CheckInput( A1, A2, B1, B2, C );
#endif
    const Grid& g = C.GetGrid();

    DistMatrix<T,MC,Star> A1T(g),  A2T(g),
                          A1B(g),  A2B(g);

    DistMatrix<T,Star,MR> B1L(g), B1R(g),
                          B2L(g), B2R(g);

    DistMatrix<T,MC,MR> CTL(g), CTR(g),
                        CBL(g), CBR(g);

    DistMatrix<T,MC,MR> DTL(g), DBR(g);

    const unsigned half = C.Height()/2;

    blas::Scal( beta, C );

    LockedPartitionDown
    ( A1, A1T,
          A1B, half );
    LockedPartitionDown
    ( A2, A2T,
          A2B, half );

    LockedPartitionRight( B1, B1L, B1R, half );
    LockedPartitionRight( B2, B2L, B2R, half );

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
        ( Normal, Normal, alpha, A1B, B2L, (T)1, CBL );

        blas::internal::LocalGemm
        ( Normal, Normal, alpha, A2B, B1L, (T)1, CBL );
    }
    else
    {
        blas::internal::LocalGemm
        ( Normal, Normal, alpha, A1T, B2R, (T)1, CTR );

        blas::internal::LocalGemm
        ( Normal, Normal, alpha, A2T, B1R, (T)1, CTR );
    }

    blas::internal::LocalGemm
    ( Normal, Normal, alpha, A1T, B2L, (T)0, DTL );

    blas::internal::LocalGemm
    ( Normal, Normal, alpha, A2T, B1L, (T)1, DTL );

    DTL.MakeTrapezoidal( Left, shape );
    blas::Axpy( (T)1, DTL, CTL );

    blas::internal::LocalGemm
    ( Normal, Normal, alpha, A1B, B2R, (T)0, DBR );

    blas::internal::LocalGemm
    ( Normal, Normal, alpha, A2B, B1R, (T)1, DBR );

    DBR.MakeTrapezoidal( Left, shape );
    blas::Axpy( (T)1, DBR, CBR );
    //------------------------------------------------------------------------//
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
LocalTriangularRank2KKernel
( Shape shape,
  Orientation orientationOfA1,
  Orientation orientationOfA2,
  Orientation orientationOfB2,
  T alpha, const DistMatrix<T,Star,MC  >& A1,
           const DistMatrix<T,Star,MC  >& A2,
           const DistMatrix<T,Star,MR  >& B1,
           const DistMatrix<T,MR,  Star>& B2,
  T beta,        DistMatrix<T,MC,  MR>& C )
{
#ifndef RELEASE
    PushCallStack("LocalTriangularRank2KKernel");
    CheckInput
    ( orientationOfA1, orientationOfA2, orientationOfB2, A1, A2, B1, B2, C );
#endif
    const Grid& g = C.GetGrid();

    DistMatrix<T,Star,MC> A1L(g), A1R(g),
                          A2L(g), A2R(g);

    DistMatrix<T,Star,MR> B1L(g), B1R(g);

    DistMatrix<T,MR,Star> B2T(g),
                          B2B(g);

    DistMatrix<T,MC,MR> CTL(g), CTR(g),
                        CBL(g), CBR(g);

    DistMatrix<T,MC,MR> DTL(g), DBR(g);

    const unsigned half = C.Height()/2;

    blas::Scal( beta, C );

    LockedPartitionRight( A1, A1L, A1R, half );
    LockedPartitionRight( A2, A2L, A2R, half );

    LockedPartitionRight( B1, B1L, B1R, half );
    LockedPartitionDown
    ( B2, B2T,
          B2B, half );

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
        ( orientationOfA1, orientationOfB2, alpha, A1R, B2T, (T)1, CBL );
        
        blas::internal::LocalGemm
        ( orientationOfA2, Normal, alpha, A2R, B1L, (T)1, CBL );
    }
    else
    {
        blas::internal::LocalGemm
        ( orientationOfA1, orientationOfB2, alpha, A1L, B2B, (T)1, CTR );

        blas::internal::LocalGemm
        ( orientationOfA2, Normal, alpha, A2L, B1R, (T)1, CTR );
    }

    blas::internal::LocalGemm
    ( orientationOfA1, orientationOfB2, alpha, A1L, B2T, (T)0, DTL );

    blas::internal::LocalGemm
    ( orientationOfA2, Normal, alpha, A2L, B1L, (T)1, DTL );

    DTL.MakeTrapezoidal( Left, shape );
    blas::Axpy( (T)1, DTL, CTL );

    blas::internal::LocalGemm
    ( orientationOfA1, orientationOfB2, alpha, A1R, B2B, (T)0, DBR );

    blas::internal::LocalGemm
    ( orientationOfA2, Normal, alpha, A2R, B1R, (T)1, DBR );

    DBR.MakeTrapezoidal( Left, shape );
    blas::Axpy( (T)1, DBR, CBR );
    //------------------------------------------------------------------------//
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
LocalTriangularRank2KKernel
( Shape shape,
  Orientation orientationOfA1,
  Orientation orientationOfA2,
  Orientation orientationOfB1,
  T alpha, const DistMatrix<T,Star,MC  >& A1,
           const DistMatrix<T,Star,MC  >& A2,
           const DistMatrix<T,MR,  Star>& B1,
           const DistMatrix<T,Star,MR  >& B2,
  T beta,        DistMatrix<T,MC,  MR  >& C )
{
#ifndef RELEASE
    PushCallStack("LocalTriangularRank2KKernel");
    CheckInput
    ( orientationOfA1, orientationOfA2, orientationOfB1, A1, A2, B1, B2, C );
#endif
    const Grid& g = C.GetGrid();

    DistMatrix<T,Star,MC> A1L(g), A1R(g),
                          A2L(g), A2R(g);

    DistMatrix<T,MR,Star> B1T(g),
                          B1B(g);

    DistMatrix<T,Star,MR> B2L(g), B2R(g);

    DistMatrix<T,MC,MR> CTL(g), CTR(g),
                        CBL(g), CBR(g);

    DistMatrix<T,MC,MR> DTL(g), DBR(g);

    const unsigned half = C.Height()/2;

    blas::Scal( beta, C );

    LockedPartitionRight( A1, A1L, A1R, half );
    LockedPartitionRight( A2, A2L, A2R, half );

    LockedPartitionDown
    ( B1, B1T,
          B1B, half );
    LockedPartitionRight( B2, B2L, B2R, half );

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
        ( orientationOfA1, Normal, alpha, A1R, B2L, (T)1, CBL );
        
        blas::internal::LocalGemm
        ( orientationOfA2, orientationOfB1, alpha, A2R, B1T, (T)1, CBL );
    }
    else
    {
        blas::internal::LocalGemm
        ( orientationOfA1, Normal, alpha, A1L, B2R, (T)1, CTR );

        blas::internal::LocalGemm
        ( orientationOfA2, orientationOfB1, alpha, A2L, B1B, (T)1, CTR );
    }

    blas::internal::LocalGemm
    ( orientationOfA1, Normal, alpha, A1L, B2L, (T)0, DTL );

    blas::internal::LocalGemm
    ( orientationOfA2, orientationOfB1, alpha, A2L, B1T, (T)1, DTL );

    DTL.MakeTrapezoidal( Left, shape );
    blas::Axpy( (T)1, DTL, CTL );

    blas::internal::LocalGemm
    ( orientationOfA1, Normal, alpha, A1R, B2R, (T)0, DBR );

    blas::internal::LocalGemm
    ( orientationOfA2, orientationOfB1, alpha, A2R, B1B, (T)1, DBR );

    DBR.MakeTrapezoidal( Left, shape );
    blas::Axpy( (T)1, DBR, CBR );
    //------------------------------------------------------------------------//
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
LocalTriangularRank2KKernel
( Shape shape,
  Orientation orientationOfA1,
  T alpha, const DistMatrix<T,Star,MC  >& A1,
           const DistMatrix<T,MC,  Star>& A2,
           const DistMatrix<T,Star,MR  >& B1,
           const DistMatrix<T,Star,MR  >& B2,
  T beta,        DistMatrix<T,MC,  MR  >& C )
{
#ifndef RELEASE
    PushCallStack("LocalTriangularRank2KKernel");
    CheckInput( orientationOfA1, A1, A2, B1, B2, C );
#endif
    const Grid& g = C.GetGrid();

    DistMatrix<T,Star,MC> A1L(g), A1R(g);

    DistMatrix<T,MC,Star> A2T(g),
                          A2B(g);

    DistMatrix<T,Star,MR> B1L(g), B1R(g),
                          B2L(g), B2R(g);

    DistMatrix<T,MC,MR> CTL(g), CTR(g),
                        CBL(g), CBR(g);

    DistMatrix<T,MC,MR> DTL(g), DBR(g);

    const unsigned half = C.Height()/2;

    blas::Scal( beta, C );

    LockedPartitionRight( A1, A1L, A1R, half );
    LockedPartitionDown
    ( A2, A2T,
          A2B, half );

    LockedPartitionRight( B1, B1L, B1R, half );
    LockedPartitionRight( B2, B2L, B2R, half );

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
        ( orientationOfA1, Normal, alpha, A1R, B2L, (T)1, CBL );
        
        blas::internal::LocalGemm
        ( Normal, Normal, alpha, A2B, B1L, (T)1, CBL );
    }
    else
    {
        blas::internal::LocalGemm
        ( orientationOfA1, Normal, alpha, A1L, B2R, (T)1, CTR );

        blas::internal::LocalGemm
        ( Normal, Normal, alpha, A2T, B1R, (T)1, CTR );
    }

    blas::internal::LocalGemm
    ( orientationOfA1, Normal, alpha, A1L, B2L, (T)0, DTL );

    blas::internal::LocalGemm
    ( Normal, Normal, alpha, A2T, B1L, (T)1, DTL );

    DTL.MakeTrapezoidal( Left, shape );
    blas::Axpy( (T)1, DTL, CTL );

    blas::internal::LocalGemm
    ( orientationOfA1, Normal, alpha, A1R, B2R, (T)0, DBR );

    blas::internal::LocalGemm
    ( Normal, Normal, alpha, A2B, B1R, (T)1, DBR );

    DBR.MakeTrapezoidal( Left, shape );
    blas::Axpy( (T)1, DBR, CBR );
    //------------------------------------------------------------------------//
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
LocalTriangularRank2KKernel
( Shape shape,
  Orientation orientationOfA2,
  T alpha, const DistMatrix<T,MC,  Star>& A1,
           const DistMatrix<T,Star,MC  >& A2,
           const DistMatrix<T,Star,MR  >& B1,
           const DistMatrix<T,Star,MR  >& B2,
  T beta,        DistMatrix<T,MC,  MR  >& C )
{
#ifndef RELEASE
    PushCallStack("LocalTriangularRank2KKernel");
    CheckInput( orientationOfA2, A1, A2, B1, B2, C );
#endif
    const Grid& g = C.GetGrid();

    DistMatrix<T,MC,Star> A1T(g),
                          A1B(g);

    DistMatrix<T,Star,MC> A2L(g), A2R(g);

    DistMatrix<T,Star,MR> B1L(g), B1R(g),
                          B2L(g), B2R(g);

    DistMatrix<T,MC,MR> CTL(g), CTR(g),
                        CBL(g), CBR(g);

    DistMatrix<T,MC,MR> DTL(g), DBR(g);

    const unsigned half = C.Height()/2;

    blas::Scal( beta, C );

    LockedPartitionDown
    ( A1, A1T,
          A1B, half );
    LockedPartitionRight( A2, A2L, A2R, half );

    LockedPartitionRight( B1, B1L, B1R, half );
    LockedPartitionRight( B2, B2L, B2R, half );

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
        ( Normal, Normal, alpha, A1B, B2L, (T)1, CBL );
        
        blas::internal::LocalGemm
        ( orientationOfA2, Normal, alpha, A2R, B1L, (T)1, CBL );
    }
    else
    {
        blas::internal::LocalGemm
        ( Normal, Normal, alpha, A1T, B2R, (T)1, CTR );

        blas::internal::LocalGemm
        ( orientationOfA2, Normal, alpha, A2L, B1R, (T)1, CTR );
    }

    blas::internal::LocalGemm
    ( Normal, Normal, alpha, A1T, B2L, (T)0, DTL );

    blas::internal::LocalGemm
    ( orientationOfA2, Normal, alpha, A2L, B1L, (T)1, DTL );

    DTL.MakeTrapezoidal( Left, shape );
    blas::Axpy( (T)1, DTL, CTL );

    blas::internal::LocalGemm
    ( Normal, Normal, alpha, A1B, B2R, (T)0, DBR );

    blas::internal::LocalGemm
    ( orientationOfA2, Normal, alpha, A2R, B1R, (T)1, DBR );

    DBR.MakeTrapezoidal( Left, shape );
    blas::Axpy( (T)1, DBR, CBR );
    //------------------------------------------------------------------------//
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
LocalTriangularRank2KKernel
( Shape shape,
  Orientation orientationOfA1,
  Orientation orientationOfA2,
  T alpha, const DistMatrix<T,Star,MC>& A1,
           const DistMatrix<T,Star,MC>& A2,
           const DistMatrix<T,Star,MR>& B1,
           const DistMatrix<T,Star,MR>& B2,
  T beta,        DistMatrix<T,MC,  MR>& C )
{
#ifndef RELEASE
    PushCallStack("LocalTriangularRank2KKernel");
    CheckInput( orientationOfA1, orientationOfA2, A1, A2, B1, B2, C );
#endif
    const Grid& g = C.GetGrid();

    DistMatrix<T,Star,MC> A1L(g), A1R(g),
                          A2L(g), A2R(g);

    DistMatrix<T,Star,MR> B1L(g), B1R(g),
                          B2L(g), B2R(g);

    DistMatrix<T,MC,MR> CTL(g), CTR(g),
                        CBL(g), CBR(g);

    DistMatrix<T,MC,MR> DTL(g), DBR(g);

    const unsigned half = C.Height()/2;

    blas::Scal( beta, C );

    LockedPartitionRight( A1, A1L, A1R, half );
    LockedPartitionRight( A2, A2L, A2R, half );

    LockedPartitionRight( B1, B1L, B1R, half );
    LockedPartitionRight( B2, B2L, B2R, half );

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
        ( orientationOfA1, Normal, alpha, A1R, B2L, (T)1, CBL );
        
        blas::internal::LocalGemm
        ( orientationOfA2, Normal, alpha, A2R, B1L, (T)1, CBL );
    }
    else
    {
        blas::internal::LocalGemm
        ( orientationOfA1, Normal, alpha, A1L, B2R, (T)1, CTR );

        blas::internal::LocalGemm
        ( orientationOfA2, Normal, alpha, A2L, B1R, (T)1, CTR );
    }

    blas::internal::LocalGemm
    ( orientationOfA1, Normal, alpha, A1L, B2L, (T)0, DTL );

    blas::internal::LocalGemm
    ( orientationOfA2, Normal, alpha, A2L, B1L, (T)1, DTL );

    DTL.MakeTrapezoidal( Left, shape );
    blas::Axpy( (T)1, DTL, CTL );

    blas::internal::LocalGemm
    ( orientationOfA1, Normal, alpha, A1R, B2R, (T)0, DBR );

    blas::internal::LocalGemm
    ( orientationOfA2, Normal, alpha, A2R, B1R, (T)1, DBR );

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
elemental::blas::internal::LocalTriangularRank2K
( Shape shape,
  Orientation orientationOfB1,
  Orientation orientationOfB2,
  T alpha, const DistMatrix<T,MC,Star>& A1,
           const DistMatrix<T,MC,Star>& A2,
           const DistMatrix<T,MR,Star>& B1,
           const DistMatrix<T,MR,Star>& B2,
  T beta,        DistMatrix<T,MC,MR  >& C  )
{
#ifndef RELEASE
    PushCallStack("blas::internal::LocalTriangularRank2K");
    CheckInput( orientationOfB1, orientationOfB2, A1, A2, B1, B2, C );
#endif
    const Grid& g = C.GetGrid();

    if( C.Height() < 2*g.Width()*Blocksize() )
    {
        LocalTriangularRank2KKernel
        ( shape, orientationOfB1, orientationOfB2, 
          alpha, A1, A2, B1, B2, beta, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.

        DistMatrix<T,MC,Star> A1T(g),  A2T(g),
                              A1B(g),  A2B(g);

        DistMatrix<T,MR,Star> B1T(g),  B2T(g),
                              B1B(g),  B2B(g);

        DistMatrix<T,MC,MR> CTL(g), CTR(g),
                            CBL(g), CBR(g);

        const unsigned half = C.Height() / 2;

        LockedPartitionDown
        ( A1, A1T,
              A1B, half );
        LockedPartitionDown
        ( A2, A2T,
              A2B, half );

        LockedPartitionDown
        ( B1, B1T,
              B1B, half );
        LockedPartitionDown
        ( B2, B2T,
              B2B, half );

        PartitionDownDiagonal
        ( C, CTL, CTR,
             CBL, CBR, half );

        if( shape == Lower )
        { 
            blas::internal::LocalGemm
            ( Normal, orientationOfB2, alpha, A1B, B2T, beta, CBL );

            blas::internal::LocalGemm
            ( Normal, orientationOfB1, alpha, A2B, B1T, (T)1, CBL );
        }
        else
        {
            blas::internal::LocalGemm
            ( Normal, orientationOfB2, alpha, A1T, B2B, beta, CTR );

            blas::internal::LocalGemm
            ( Normal, orientationOfB1, alpha, A2T, B1B, (T)1, CTR );
        }

        // Recurse
        blas::internal::LocalTriangularRank2K
        ( shape, orientationOfB1, orientationOfB2,
          alpha, A1T, A2T, B1T, B2T, beta, CTL );

        blas::internal::LocalTriangularRank2K
        ( shape, orientationOfB1, orientationOfB2,
          alpha, A1B, A2B, B1B, B2B, beta, CBR );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape,
  Orientation orientationOfA1,
  Orientation orientationOfB1,
  Orientation orientationOfB2,
  T alpha, const DistMatrix<T,Star,MC  >& A1,
           const DistMatrix<T,MC,  Star>& A2,
           const DistMatrix<T,MR,  Star>& B1,
           const DistMatrix<T,MR,  Star>& B2,
  T beta,        DistMatrix<T,MC,  MR  >& C  )
{
#ifndef RELEASE
    PushCallStack("blas::internal::LocalTriangularRank2K");
    CheckInput
    ( orientationOfA1, orientationOfB1, orientationOfB2, A1, A2, B1, B2, C );
#endif
    const Grid& g = C.GetGrid();

    if( C.Height() < 2*g.Width()*Blocksize() )
    {
        LocalTriangularRank2KKernel
        ( shape, orientationOfA1, orientationOfB1, orientationOfB2, 
          alpha, A1, A2, B1, B2, beta, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.

        DistMatrix<T,Star,MC> A1L(g), A1R(g);

        DistMatrix<T,MC,Star> A2T(g),
                              A2B(g); 

        DistMatrix<T,MR,Star> B1T(g),  B2T(g),
                              B1B(g),  B2B(g);

        DistMatrix<T,MC,MR> CTL(g), CTR(g),
                            CBL(g), CBR(g);

        const unsigned half = C.Height() / 2;

        LockedPartitionRight( A1, A1L, A1R, half );
        LockedPartitionDown
        ( A2, A2T,
              A2B, half );

        LockedPartitionDown
        ( B1, B1T,
              B1B, half );
        LockedPartitionDown
        ( B2, B2T,
              B2B, half );

        PartitionDownDiagonal
        ( C, CTL, CTR,
             CBL, CBR, half );

        if( shape == Lower )
        { 
            blas::internal::LocalGemm
            ( orientationOfA1, orientationOfB2, alpha, A1R, B2T, beta, CBL );

            blas::internal::LocalGemm
            ( Normal, orientationOfB1, alpha, A2B, B1T, (T)1, CBL );
        }
        else
        {
            blas::internal::LocalGemm
            ( orientationOfA1, orientationOfB2, alpha, A1L, B2B, beta, CTR );

            blas::internal::LocalGemm
            ( Normal, orientationOfB1, alpha, A2T, B1B, (T)1, CTR );
        }

        // Recurse
        blas::internal::LocalTriangularRank2K
        ( shape, orientationOfA1, orientationOfB1, orientationOfB2,
          alpha, A1L, A2T, B1T, B2T, beta, CTL );

        blas::internal::LocalTriangularRank2K
        ( shape, orientationOfA1, orientationOfB1, orientationOfB2,
          alpha, A1R, A2B, B1B, B2B, beta, CBR );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape,
  Orientation orientationOfA2,
  Orientation orientationOfB1,
  Orientation orientationOfB2,
  T alpha, const DistMatrix<T,MC,  Star>& A1,
           const DistMatrix<T,Star,MC  >& A2,
           const DistMatrix<T,MR,  Star>& B1,
           const DistMatrix<T,MR,  Star>& B2,
  T beta,        DistMatrix<T,MC,  MR  >& C  )
{
#ifndef RELEASE
    PushCallStack("blas::internal::LocalTriangularRank2K");
    CheckInput
    ( orientationOfA2, orientationOfB1, orientationOfB2, A1, A2, B1, B2, C );
#endif
    const Grid& g = C.GetGrid();

    if( C.Height() < 2*g.Width()*Blocksize() )
    {
        LocalTriangularRank2KKernel
        ( shape, orientationOfA2, orientationOfB1, orientationOfB2, 
          alpha, A1, A2, B1, B2, beta, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.

        DistMatrix<T,MC,Star> A1T(g),
                              A1B(g); 

        DistMatrix<T,Star,MC> A2L(g), A2R(g);

        DistMatrix<T,MR,Star> B1T(g),  B2T(g),
                              B1B(g),  B2B(g);

        DistMatrix<T,MC,MR> CTL(g), CTR(g),
                            CBL(g), CBR(g);

        const unsigned half = C.Height() / 2;

        LockedPartitionDown
        ( A1, A1T,
              A1B, half );
        LockedPartitionRight( A2, A2L, A2R, half );

        LockedPartitionDown
        ( B1, B1T,
              B1B, half );
        LockedPartitionDown
        ( B2, B2T,
              B2B, half );

        PartitionDownDiagonal
        ( C, CTL, CTR,
             CBL, CBR, half );

        if( shape == Lower )
        { 
            blas::internal::LocalGemm
            ( Normal, orientationOfB2, alpha, A1B, B2T, beta, CBL );

            blas::internal::LocalGemm
            ( orientationOfA2, orientationOfB1, alpha, A2R, B1T, (T)1, CBL );
        }
        else
        {
            blas::internal::LocalGemm
            ( Normal, orientationOfB2, alpha, A1T, B2B, beta, CTR );

            blas::internal::LocalGemm
            ( orientationOfA2, orientationOfB1, alpha, A2L, B1B, (T)1, CTR );
        }

        // Recurse
        blas::internal::LocalTriangularRank2K
        ( shape, orientationOfA2, orientationOfB1, orientationOfB2,
          alpha, A1T, A2L, B1T, B2T, beta, CTL );

        blas::internal::LocalTriangularRank2K
        ( shape, orientationOfA2, orientationOfB1, orientationOfB2,
          alpha, A1B, A2R, B1B, B2B, beta, CBR );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape,
  Orientation orientationOfB2,
  T alpha, const DistMatrix<T,MC,  Star>& A1,
           const DistMatrix<T,MC,  Star>& A2,
           const DistMatrix<T,Star,MR  >& B1,
           const DistMatrix<T,MR,  Star>& B2,
  T beta,        DistMatrix<T,MC,  MR  >& C  )
{
#ifndef RELEASE
    PushCallStack("blas::internal::LocalTriangularRank2K");
    CheckInput( orientationOfB2, A1, A2, B1, B2, C );
#endif
    const Grid& g = C.GetGrid();

    if( C.Height() < 2*g.Width()*Blocksize() )
    {
        LocalTriangularRank2KKernel
        ( shape, orientationOfB2, alpha, A1, A2, B1, B2, beta, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.

        DistMatrix<T,MC,Star> A1T(g),  A2T(g),
                              A1B(g),  A2B(g);

        DistMatrix<T,Star,MR> B1L(g), B1R(g);
        
        DistMatrix<T,MR,Star> B2T(g), 
                              B2B(g);

        DistMatrix<T,MC,MR> CTL(g), CTR(g),
                            CBL(g), CBR(g);

        const unsigned half = C.Height() / 2;

        LockedPartitionDown
        ( A1, A1T,
              A1B, half );
        LockedPartitionDown
        ( A2, A2T,
              A2B, half );

        LockedPartitionRight( B1, B1L, B1R, half );
        LockedPartitionDown
        ( B2, B2T,
              B2B, half );

        PartitionDownDiagonal
        ( C, CTL, CTR,
             CBL, CBR, half );

        if( shape == Lower )
        { 
            blas::internal::LocalGemm
            ( Normal, orientationOfB2, alpha, A1B, B2T, beta, CBL );

            blas::internal::LocalGemm
            ( Normal, Normal, alpha, A2B, B1L, (T)1, CBL );
        }
        else
        {
            blas::internal::LocalGemm
            ( Normal, orientationOfB2, alpha, A1T, B2B, beta, CTR );

            blas::internal::LocalGemm
            ( Normal, Normal, alpha, A2T, B1R, (T)1, CTR );
        }

        // Recurse
        blas::internal::LocalTriangularRank2K
        ( shape, orientationOfB2, alpha, A1T, A2T, B1L, B2T, beta, CTL );

        blas::internal::LocalTriangularRank2K
        ( shape, orientationOfB2, alpha, A1B, A2B, B1R, B2B, beta, CBR );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape,
  Orientation orientationOfB1,
  T alpha, const DistMatrix<T,MC,  Star>& A1,
           const DistMatrix<T,MC,  Star>& A2,
           const DistMatrix<T,MR,  Star>& B1,
           const DistMatrix<T,Star,MR  >& B2,
  T beta,        DistMatrix<T,MC,  MR  >& C  )
{
#ifndef RELEASE
    PushCallStack("blas::internal::LocalTriangularRank2K");
    CheckInput( orientationOfB1, A1, A2, B1, B2, C );
#endif
    const Grid& g = C.GetGrid();

    if( C.Height() < 2*g.Width()*Blocksize() )
    {
        LocalTriangularRank2KKernel
        ( shape, orientationOfB1, alpha, A1, A2, B1, B2, beta, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.

        DistMatrix<T,MC,Star> A1T(g),  A2T(g),
                              A1B(g),  A2B(g);

        DistMatrix<T,MR,Star> B1T(g), 
                              B1B(g);

        DistMatrix<T,Star,MR> B2L(g), B2R(g);

        DistMatrix<T,MC,MR> CTL(g), CTR(g),
                            CBL(g), CBR(g);

        const unsigned half = C.Height() / 2;

        LockedPartitionDown
        ( A1, A1T,
              A1B, half );
        LockedPartitionDown
        ( A2, A2T,
              A2B, half );

        LockedPartitionDown
        ( B1, B1T,
              B1B, half );
        LockedPartitionRight( B2, B2L, B2R, half );

        PartitionDownDiagonal
        ( C, CTL, CTR,
             CBL, CBR, half );

        if( shape == Lower )
        { 
            blas::internal::LocalGemm
            ( Normal, Normal, alpha, A1B, B2L, beta, CBL );

            blas::internal::LocalGemm
            ( Normal, orientationOfB1, alpha, A2B, B1T, (T)1, CBL );
        }
        else
        {
            blas::internal::LocalGemm
            ( Normal, Normal, alpha, A1T, B2R, beta, CTR );

            blas::internal::LocalGemm
            ( Normal, orientationOfB1, alpha, A2T, B1B, (T)1, CTR );
        }

        // Recurse
        blas::internal::LocalTriangularRank2K
        ( shape, orientationOfB1, alpha, A1T, A2T, B1T, B2L, beta, CTL );

        blas::internal::LocalTriangularRank2K
        ( shape, orientationOfB1, alpha, A1B, A2B, B1B, B2R, beta, CBR );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape,
  Orientation orientationOfA1,
  Orientation orientationOfA2,
  Orientation orientationOfB1,
  Orientation orientationOfB2,
  T alpha, const DistMatrix<T,Star,MC  >& A1,
           const DistMatrix<T,Star,MC  >& A2,
           const DistMatrix<T,MR,  Star>& B1,
           const DistMatrix<T,MR,  Star>& B2,
  T beta,        DistMatrix<T,MC,  MR  >& C  )
{
#ifndef RELEASE
    PushCallStack("blas::internal::LocalTriangularRank2K");
    CheckInput
    ( orientationOfA1, orientationOfA2, orientationOfB1, orientationOfB2, 
      A1, A2, B1, B2, C );
#endif
    const Grid& g = C.GetGrid();

    if( C.Height() < 2*g.Width()*Blocksize() )
    {
        LocalTriangularRank2KKernel
        ( shape, 
          orientationOfA1, orientationOfA2, 
          orientationOfB1, orientationOfB2, 
          alpha, A1, A2, B1, B2, beta, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.

        DistMatrix<T,Star,MC> A1L(g), A1R(g),
                              A2L(g), A2R(g);

        DistMatrix<T,MR,Star> B1T(g),  B2T(g),
                              B1B(g),  B2B(g);

        DistMatrix<T,MC,MR> CTL(g), CTR(g),
                            CBL(g), CBR(g);

        const unsigned half = C.Height() / 2;

        LockedPartitionRight( A1, A1L, A1R, half );
        LockedPartitionRight( A2, A2L, A2R, half );

        LockedPartitionDown
        ( B1, B1T,
              B1B, half );
        LockedPartitionDown
        ( B2, B2T,
              B2B, half );

        PartitionDownDiagonal
        ( C, CTL, CTR,
             CBL, CBR, half );

        if( shape == Lower )
        { 
            blas::internal::LocalGemm
            ( orientationOfA1, orientationOfB2, alpha, A1R, B2T, beta, CBL );

            blas::internal::LocalGemm
            ( orientationOfA2, orientationOfB1, alpha, A2R, B1T, (T)1, CBL );
        }
        else
        {
            blas::internal::LocalGemm
            ( orientationOfA1, orientationOfB2, alpha, A1L, B2B, beta, CTR );

            blas::internal::LocalGemm
            ( orientationOfA2, orientationOfB1, alpha, A2L, B1B, (T)1, CTR );
        }

        // Recurse
        blas::internal::LocalTriangularRank2K
        ( shape, 
          orientationOfA1, orientationOfA2, 
          orientationOfB1, orientationOfB2, 
          alpha, A1L, A2L, B1T, B2T, beta, CTL );

        blas::internal::LocalTriangularRank2K
        ( shape, 
          orientationOfA1, orientationOfA2,
          orientationOfB1, orientationOfB2,
          alpha, A1R, A2R, B1B, B2B, beta, CBR );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape,
  Orientation orientationOfA1,
  Orientation orientationOfB2,
  T alpha, const DistMatrix<T,Star,MC  >& A1,
           const DistMatrix<T,MC,  Star>& A2,
           const DistMatrix<T,Star,MR  >& B1,
           const DistMatrix<T,MR,  Star>& B2,
  T beta,        DistMatrix<T,MC,  MR  >& C  )
{
#ifndef RELEASE
    PushCallStack("blas::internal::LocalTriangularRank2K");
    CheckInput
    ( orientationOfA1, orientationOfB2, A1, A2, B1, B2, C );
#endif
    const Grid& g = C.GetGrid();

    if( C.Height() < 2*g.Width()*Blocksize() )
    {
        LocalTriangularRank2KKernel
        ( shape, orientationOfA1, orientationOfB2, 
          alpha, A1, A2, B1, B2, beta, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.

        DistMatrix<T,Star,MC> A1L(g), A1R(g);

        DistMatrix<T,MC,Star> A2T(g),
                              A2B(g); 

        DistMatrix<T,Star,MR> B1L(g), B1R(g);

        DistMatrix<T,MR,Star> B2T(g),
                              B2B(g);

        DistMatrix<T,MC,MR> CTL(g), CTR(g),
                            CBL(g), CBR(g);

        const unsigned half = C.Height() / 2;

        LockedPartitionRight( A1, A1L, A1R, half );
        LockedPartitionDown
        ( A2, A2T,
              A2B, half );

        LockedPartitionRight( B1, B1L, B1R, half );
        LockedPartitionDown
        ( B2, B2T,
              B2B, half );

        PartitionDownDiagonal
        ( C, CTL, CTR,
             CBL, CBR, half );

        if( shape == Lower )
        { 
            blas::internal::LocalGemm
            ( orientationOfA1, orientationOfB2, alpha, A1R, B2T, beta, CBL );

            blas::internal::LocalGemm
            ( Normal, Normal, alpha, A2B, B1L, (T)1, CBL );
        }
        else
        {
            blas::internal::LocalGemm
            ( orientationOfA1, orientationOfB2, alpha, A1L, B2B, beta, CTR );

            blas::internal::LocalGemm
            ( Normal, Normal, alpha, A2T, B1R, (T)1, CTR );
        }

        // Recurse
        blas::internal::LocalTriangularRank2K
        ( shape, orientationOfA1, orientationOfB2,
          alpha, A1L, A2T, B1L, B2T, beta, CTL );

        blas::internal::LocalTriangularRank2K
        ( shape, orientationOfA1, orientationOfB2,
          alpha, A1R, A2B, B1R, B2B, beta, CBR );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape,
  Orientation orientationOfA1,
  Orientation orientationOfB1,
  T alpha, const DistMatrix<T,Star,MC  >& A1,
           const DistMatrix<T,MC,  Star>& A2,
           const DistMatrix<T,MR,  Star>& B1,
           const DistMatrix<T,Star,MR  >& B2,
  T beta,        DistMatrix<T,MC,  MR  >& C  )
{
#ifndef RELEASE
    PushCallStack("blas::internal::LocalTriangularRank2K");
    CheckInput
    ( orientationOfA1, orientationOfB1, A1, A2, B1, B2, C );
#endif
    const Grid& g = C.GetGrid();

    if( C.Height() < 2*g.Width()*Blocksize() )
    {
        LocalTriangularRank2KKernel
        ( shape, orientationOfA1, orientationOfB1,  
          alpha, A1, A2, B1, B2, beta, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.

        DistMatrix<T,Star,MC> A1L(g), A1R(g);

        DistMatrix<T,MC,Star> A2T(g),
                              A2B(g); 

        DistMatrix<T,MR,Star> B1T(g),
                              B1B(g);

        DistMatrix<T,Star,MR> B2L(g), B2R(g);

        DistMatrix<T,MC,MR> CTL(g), CTR(g),
                            CBL(g), CBR(g);

        const unsigned half = C.Height() / 2;

        LockedPartitionRight( A1, A1L, A1R, half );
        LockedPartitionDown
        ( A2, A2T,
              A2B, half );

        LockedPartitionDown
        ( B1, B1T,
              B1B, half );
        LockedPartitionRight( B2, B2L, B2R, half );

        PartitionDownDiagonal
        ( C, CTL, CTR,
             CBL, CBR, half );

        if( shape == Lower )
        { 
            blas::internal::LocalGemm
            ( orientationOfA1, Normal, alpha, A1R, B2L, beta, CBL );

            blas::internal::LocalGemm
            ( Normal, orientationOfB1, alpha, A2B, B1T, (T)1, CBL );
        }
        else
        {
            blas::internal::LocalGemm
            ( orientationOfA1, Normal, alpha, A1L, B2R, beta, CTR );

            blas::internal::LocalGemm
            ( Normal, orientationOfB1, alpha, A2T, B1B, (T)1, CTR );
        }

        // Recurse
        blas::internal::LocalTriangularRank2K
        ( shape, orientationOfA1, orientationOfB1,
          alpha, A1L, A2T, B1T, B2L, beta, CTL );

        blas::internal::LocalTriangularRank2K
        ( shape, orientationOfA1, orientationOfB1,
          alpha, A1R, A2B, B1B, B2R, beta, CBR );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape,
  Orientation orientationOfA2,
  Orientation orientationOfB2,
  T alpha, const DistMatrix<T,MC,  Star>& A1,
           const DistMatrix<T,Star,MC  >& A2,
           const DistMatrix<T,Star,MR  >& B1,
           const DistMatrix<T,MR,  Star>& B2,
  T beta,        DistMatrix<T,MC,  MR  >& C  )
{
#ifndef RELEASE
    PushCallStack("blas::internal::LocalTriangularRank2K");
    CheckInput
    ( orientationOfA2, orientationOfB2, A1, A2, B1, B2, C );
#endif
    const Grid& g = C.GetGrid();

    if( C.Height() < 2*g.Width()*Blocksize() )
    {
        LocalTriangularRank2KKernel
        ( shape, orientationOfA2, orientationOfB2, 
          alpha, A1, A2, B1, B2, beta, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.

        DistMatrix<T,MC,Star> A1T(g),
                              A1B(g); 

        DistMatrix<T,Star,MC> A2L(g), A2R(g);

        DistMatrix<T,Star,MR> B1L(g), B1R(g);

        DistMatrix<T,MR,Star> B2T(g),
                              B2B(g);

        DistMatrix<T,MC,MR> CTL(g), CTR(g),
                            CBL(g), CBR(g);

        const unsigned half = C.Height() / 2;

        LockedPartitionDown
        ( A1, A1T,
              A1B, half );
        LockedPartitionRight( A2, A2L, A2R, half );

        LockedPartitionRight( B1, B1L, B1R, half );
        LockedPartitionDown
        ( B2, B2T,
              B2B, half );

        PartitionDownDiagonal
        ( C, CTL, CTR,
             CBL, CBR, half );

        if( shape == Lower )
        { 
            blas::internal::LocalGemm
            ( Normal, orientationOfB2, alpha, A1B, B2T, beta, CBL );

            blas::internal::LocalGemm
            ( orientationOfA2, Normal, alpha, A2R, B1L, (T)1, CBL );
        }
        else
        {
            blas::internal::LocalGemm
            ( Normal, orientationOfB2, alpha, A1T, B2B, beta, CTR );

            blas::internal::LocalGemm
            ( orientationOfA2, Normal, alpha, A2L, B1R, (T)1, CTR );
        }

        // Recurse
        blas::internal::LocalTriangularRank2K
        ( shape, orientationOfA2, orientationOfB2,
          alpha, A1T, A2L, B1L, B2T, beta, CTL );

        blas::internal::LocalTriangularRank2K
        ( shape, orientationOfA2, orientationOfB2,
          alpha, A1B, A2R, B1R, B2B, beta, CBR );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape,
  Orientation orientationOfA2,
  Orientation orientationOfB1,
  T alpha, const DistMatrix<T,MC,  Star>& A1,
           const DistMatrix<T,Star,MC  >& A2,
           const DistMatrix<T,MR,  Star>& B1,
           const DistMatrix<T,Star,MR  >& B2,
  T beta,        DistMatrix<T,MC,  MR  >& C  )
{
#ifndef RELEASE
    PushCallStack("blas::internal::LocalTriangularRank2K");
    CheckInput
    ( orientationOfA2, orientationOfB1, A1, A2, B1, B2, C );
#endif
    const Grid& g = C.GetGrid();

    if( C.Height() < 2*g.Width()*Blocksize() )
    {
        LocalTriangularRank2KKernel
        ( shape, orientationOfA2, orientationOfB1, 
          alpha, A1, A2, B1, B2, beta, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.

        DistMatrix<T,MC,Star> A1T(g),
                              A1B(g); 

        DistMatrix<T,Star,MC> A2L(g), A2R(g);

        DistMatrix<T,MR,Star> B1T(g),
                              B1B(g);

        DistMatrix<T,Star,MR> B2L(g), B2R(g);

        DistMatrix<T,MC,MR> CTL(g), CTR(g),
                            CBL(g), CBR(g);

        const unsigned half = C.Height() / 2;

        LockedPartitionDown
        ( A1, A1T,
              A1B, half );
        LockedPartitionRight( A2, A2L, A2R, half );

        LockedPartitionDown
        ( B1, B1T,
              B1B, half );
        LockedPartitionRight( B2, B2L, B2R, half );

        PartitionDownDiagonal
        ( C, CTL, CTR,
             CBL, CBR, half );

        if( shape == Lower )
        { 
            blas::internal::LocalGemm
            ( Normal, Normal, alpha, A1B, B2L, beta, CBL );

            blas::internal::LocalGemm
            ( orientationOfA2, orientationOfB1, alpha, A2R, B1T, (T)1, CBL );
        }
        else
        {
            blas::internal::LocalGemm
            ( Normal, Normal, alpha, A1T, B2R, beta, CTR );

            blas::internal::LocalGemm
            ( orientationOfA2, orientationOfB1, alpha, A2L, B1B, (T)1, CTR );
        }

        // Recurse
        blas::internal::LocalTriangularRank2K
        ( shape, orientationOfA2, orientationOfB1, 
          alpha, A1T, A2L, B1T, B2L, beta, CTL );

        blas::internal::LocalTriangularRank2K
        ( shape, orientationOfA2, orientationOfB1,
          alpha, A1B, A2R, B1B, B2R, beta, CBR );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape,
  T alpha, const DistMatrix<T,MC,  Star>& A1,
           const DistMatrix<T,MC,  Star>& A2,
           const DistMatrix<T,Star,MR  >& B1,
           const DistMatrix<T,Star,MR  >& B2,
  T beta,        DistMatrix<T,MC,  MR  >& C )
{
#ifndef RELEASE
    PushCallStack("blas::internal::LocalTriangularRank2K");
    CheckInput( A1, A2, B1, B2, C );
#endif
    const Grid& g = C.GetGrid();

    if( C.Height() < 2*g.Width()*Blocksize() )
    {
        LocalTriangularRank2KKernel
        ( shape, alpha, A1, A2, B1, B2, beta, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.

        DistMatrix<T,MC,Star> A1T(g),  A2T(g),
                              A1B(g),  A2B(g);

        DistMatrix<T,Star,MR> B1L(g), B1R(g),
                              B2L(g), B2R(g);

        DistMatrix<T,MC,MR> CTL(g), CTR(g),
                            CBL(g), CBR(g);

        const unsigned half = C.Height() / 2;

        LockedPartitionDown
        ( A1, A1T,
              A1B, half );
        LockedPartitionDown
        ( A2, A2T,
              A2B, half );

        LockedPartitionRight( B1, B1L, B1R, half );
        LockedPartitionRight( B2, B2L, B2R, half );

        PartitionDownDiagonal
        ( C, CTL, CTR,
             CBL, CBR, half );

        if( shape == Lower )
        { 
            blas::internal::LocalGemm
            ( Normal, Normal, alpha, A1B, B2L, beta, CBL );

            blas::internal::LocalGemm
            ( Normal, Normal, alpha, A2B, B1L, (T)1, CBL );
        }
        else
        {
            blas::internal::LocalGemm
            ( Normal, Normal, alpha, A1T, B2R, beta, CTR );
            
            blas::internal::LocalGemm
            ( Normal, Normal, alpha, A2T, B1R, (T)1, CTR );
        }

        // Recurse
        blas::internal::LocalTriangularRank2K
        ( shape, alpha, A1T, A2T, B1L, B2L, beta, CTL );

        blas::internal::LocalTriangularRank2K
        ( shape, alpha, A1B, A2B, B1R, B2R, beta, CBR );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape,
  Orientation orientationOfA1,
  Orientation orientationOfA2,
  Orientation orientationOfB2,
  T alpha, const DistMatrix<T,Star,MC  >& A1,
           const DistMatrix<T,Star,MC  >& A2,
           const DistMatrix<T,Star,MR  >& B1,
           const DistMatrix<T,MR,  Star>& B2,
  T beta,        DistMatrix<T,MC,  MR>& C  )
{
#ifndef RELEASE
    PushCallStack("blas::internal::LocalTriangularRank2K");
    CheckInput
    ( orientationOfA1, orientationOfA2, orientationOfB2, A1, A2, B1, B2, C );
#endif
    const Grid& g = C.GetGrid();

    if( C.Height() < 2*g.Width()*Blocksize() )
    {
        LocalTriangularRank2KKernel
        ( shape, orientationOfA1, orientationOfA2, orientationOfB2,
          alpha, A1, A2, B1, B2, beta, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.

        DistMatrix<T,Star,MC> A1L(g), A1R(g),
                              A2L(g), A2R(g);

        DistMatrix<T,Star,MR> B1L(g), B1R(g);

        DistMatrix<T,MR,Star> B2T(g),
                              B2B(g);

        DistMatrix<T,MC,MR> CTL(g), CTR(g),
                            CBL(g), CBR(g);

        const unsigned half = C.Height() / 2;

        LockedPartitionRight( A1, A1L, A1R, half );
        LockedPartitionRight( A2, A2L, A2R, half );

        LockedPartitionRight( B1, B1L, B1R, half );
        LockedPartitionDown
        ( B2, B2T, 
              B2B, half );

        PartitionDownDiagonal
        ( C, CTL, CTR,
             CBL, CBR, half );

        if( shape == Lower )
        { 
            blas::internal::LocalGemm
            ( orientationOfA1, orientationOfB2, alpha, A1R, B2T, beta, CBL );

            blas::internal::LocalGemm
            ( orientationOfA2, Normal, alpha, A2R, B1L, (T)1, CBL );
        }
        else
        {
            blas::internal::LocalGemm
            ( orientationOfA1, orientationOfB2, alpha, A1L, B2B, beta, CTR );

            blas::internal::LocalGemm
            ( orientationOfA2, Normal, alpha, A2L, B1R, (T)1, CTR );
        }

        // Recurse
        blas::internal::LocalTriangularRank2K
        ( shape, orientationOfA1, orientationOfA2, orientationOfB2,
          alpha, A1L, A2L, B1L, B2T, beta, CTL );

        blas::internal::LocalTriangularRank2K
        ( shape, orientationOfA1, orientationOfA2, orientationOfB2,
          alpha, A1R, A2R, B1R, B2B, beta, CBR );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape,
  Orientation orientationOfA1,
  Orientation orientationOfA2,
  Orientation orientationOfB1,
  T alpha, const DistMatrix<T,Star,MC  >& A1,
           const DistMatrix<T,Star,MC  >& A2,
           const DistMatrix<T,MR,  Star>& B1,
           const DistMatrix<T,Star,MR  >& B2,
  T beta,        DistMatrix<T,MC,  MR  >& C  )
{
#ifndef RELEASE
    PushCallStack("blas::internal::LocalTriangularRank2K");
    CheckInput
    ( orientationOfA1, orientationOfA2, orientationOfB1, A1, A2, B1, B2, C );
#endif
    const Grid& g = C.GetGrid();

    if( C.Height() < 2*g.Width()*Blocksize() )
    {
        LocalTriangularRank2KKernel
        ( shape, orientationOfA1, orientationOfA2, orientationOfB1,
          alpha, A1, A2, B1, B2, beta, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.

        DistMatrix<T,Star,MC> A1L(g), A1R(g),
                              A2L(g), A2R(g);

        DistMatrix<T,MR,Star> B1T(g),
                              B1B(g);

        DistMatrix<T,Star,MR> B2L(g), B2R(g);

        DistMatrix<T,MC,MR> CTL(g), CTR(g),
                            CBL(g), CBR(g);

        const unsigned half = C.Height() / 2;

        LockedPartitionRight( A1, A1L, A1R, half );
        LockedPartitionRight( A2, A2L, A2R, half );

        LockedPartitionDown
        ( B1, B1T,
              B1B, half );
        LockedPartitionRight( B2, B2L, B2R, half );

        PartitionDownDiagonal
        ( C, CTL, CTR,
             CBL, CBR, half );

        if( shape == Lower )
        { 
            blas::internal::LocalGemm
            ( orientationOfA1, Normal, alpha, A1R, B2L, beta, CBL );

            blas::internal::LocalGemm
            ( orientationOfA2, orientationOfB1, alpha, A2R, B1T, (T)1, CBL );
        }
        else
        {
            blas::internal::LocalGemm
            ( orientationOfA1, Normal, alpha, A1L, B2R, beta, CTR );

            blas::internal::LocalGemm
            ( orientationOfA2, orientationOfB1, alpha, A2L, B1B, (T)1, CTR );
        }

        // Recurse
        blas::internal::LocalTriangularRank2K
        ( shape, orientationOfA1, orientationOfA2, orientationOfB1,
          alpha, A1L, A2L, B1T, B2L, beta, CTL );

        blas::internal::LocalTriangularRank2K
        ( shape, orientationOfA1, orientationOfA2, orientationOfB1,
          alpha, A1R, A2R, B1B, B2R, beta, CBR );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape,
  Orientation orientationOfA1,
  T alpha, const DistMatrix<T,Star,MC  >& A1,
           const DistMatrix<T,MC,  Star>& A2,
           const DistMatrix<T,Star,MR  >& B1,
           const DistMatrix<T,Star,MR  >& B2,
  T beta,        DistMatrix<T,MC,  MR  >& C  )
{
#ifndef RELEASE
    PushCallStack("blas::internal::LocalTriangularRank2K");
    CheckInput( orientationOfA1, A1, A2, B1, B2, C );
#endif
    const Grid& g = C.GetGrid();

    if( C.Height() < 2*g.Width()*Blocksize() )
    {
        LocalTriangularRank2KKernel
        ( shape, orientationOfA1, alpha, A1, A2, B1, B2, beta, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.

        DistMatrix<T,Star,MC> A1L(g), A1R(g);

        DistMatrix<T,MC,Star> A2T(g),
                              A2B(g);

        DistMatrix<T,Star,MR> B1L(g), B1R(g),
                              B2L(g), B2R(g);

        DistMatrix<T,MC,MR> CTL(g), CTR(g),
                            CBL(g), CBR(g);

        const unsigned half = C.Height() / 2;

        LockedPartitionRight( A1, A1L, A1R, half );
        LockedPartitionDown
        ( A2, A2T,
              A2B, half );

        LockedPartitionRight( B1, B1L, B1R, half );
        LockedPartitionRight( B2, B2L, B2R, half );

        PartitionDownDiagonal
        ( C, CTL, CTR,
             CBL, CBR, half );

        if( shape == Lower )
        { 
            blas::internal::LocalGemm
            ( orientationOfA1, Normal, alpha, A1R, B2L, beta, CBL );

            blas::internal::LocalGemm
            ( Normal, Normal, alpha, A2B, B1L, (T)1, CBL );
        }
        else
        {
            blas::internal::LocalGemm
            ( orientationOfA1, Normal, alpha, A1L, B2R, beta, CTR );

            blas::internal::LocalGemm
            ( Normal, Normal, alpha, A2T, B1R, (T)1, CTR );
        }

        // Recurse
        blas::internal::LocalTriangularRank2K
        ( shape, orientationOfA1, alpha, A1L, A2T, B1L, B2L, beta, CTL );

        blas::internal::LocalTriangularRank2K
        ( shape, orientationOfA1, alpha, A1R, A2B, B1R, B2R, beta, CBR );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape,
  Orientation orientationOfA2,
  T alpha, const DistMatrix<T,MC,  Star>& A1,
           const DistMatrix<T,Star,MC  >& A2,
           const DistMatrix<T,Star,MR  >& B1,
           const DistMatrix<T,Star,MR  >& B2,
  T beta,        DistMatrix<T,MC,  MR  >& C  )
{
#ifndef RELEASE
    PushCallStack("blas::internal::LocalTriangularRank2K");
    CheckInput( orientationOfA2, A1, A2, B1, B2, C );
#endif
    const Grid& g = C.GetGrid();

    if( C.Height() < 2*g.Width()*Blocksize() )
    {
        LocalTriangularRank2KKernel
        ( shape, orientationOfA2, alpha, A1, A2, B1, B2, beta, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.

        DistMatrix<T,MC,Star> A1T(g),
                              A1B(g);

        DistMatrix<T,Star,MC> A2L(g), A2R(g);

        DistMatrix<T,Star,MR> B1L(g), B1R(g),
                              B2L(g), B2R(g);

        DistMatrix<T,MC,MR> CTL(g), CTR(g),
                            CBL(g), CBR(g);

        const unsigned half = C.Height() / 2;

        LockedPartitionDown
        ( A1, A1T,
              A1B, half );
        LockedPartitionRight( A2, A2L, A2R, half );

        LockedPartitionRight( B1, B1L, B1R, half );
        LockedPartitionRight( B2, B2L, B2R, half );

        PartitionDownDiagonal
        ( C, CTL, CTR,
             CBL, CBR, half );

        if( shape == Lower )
        { 
            blas::internal::LocalGemm
            ( Normal, Normal, alpha, A1B, B2L, beta, CBL );

            blas::internal::LocalGemm
            ( orientationOfA2, Normal, alpha, A2R, B1L, (T)1, CBL );
        }
        else
        {
            blas::internal::LocalGemm
            ( Normal, Normal, alpha, A1T, B2R, beta, CTR );

            blas::internal::LocalGemm
            ( orientationOfA2, Normal, alpha, A2L, B1R, (T)1, CTR );
        }

        // Recurse
        blas::internal::LocalTriangularRank2K
        ( shape, orientationOfA2, alpha, A1T, A2L, B1L, B2L, beta, CTL );

        blas::internal::LocalTriangularRank2K
        ( shape, orientationOfA2, alpha, A1B, A2R, B1R, B2R, beta, CBR );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape,
  Orientation orientationOfA1,
  Orientation orientationOfA2,
  T alpha, const DistMatrix<T,Star,MC>& A1,
           const DistMatrix<T,Star,MC>& A2,
           const DistMatrix<T,Star,MR>& B1,
           const DistMatrix<T,Star,MR>& B2,
  T beta,        DistMatrix<T,MC,  MR>& C  )
{
#ifndef RELEASE
    PushCallStack("blas::internal::LocalTriangularRank2K");
    CheckInput( orientationOfA1, orientationOfA2, A1, A2, B1, B2, C );
#endif
    const Grid& g = C.GetGrid();

    if( C.Height() < 2*g.Width()*Blocksize() )
    {
        LocalTriangularRank2KKernel
        ( shape, orientationOfA1, orientationOfA2, 
          alpha, A1, A2, B1, B2, beta, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.

        DistMatrix<T,Star,MC> A1L(g), A1R(g),
                              A2L(g), A2R(g);

        DistMatrix<T,Star,MR> B1L(g), B1R(g),
                              B2L(g), B2R(g);

        DistMatrix<T,MC,MR> CTL(g), CTR(g),
                            CBL(g), CBR(g);

        const unsigned half = C.Height() / 2;

        LockedPartitionRight( A1, A1L, A1R, half );
        LockedPartitionRight( A2, A2L, A2R, half );

        LockedPartitionRight( B1, B1L, B1R, half );
        LockedPartitionRight( B2, B2L, B2R, half );

        PartitionDownDiagonal
        ( C, CTL, CTR,
             CBL, CBR, half );

        if( shape == Lower )
        { 
            blas::internal::LocalGemm
            ( orientationOfA1, Normal, alpha, A1R, B2L, beta, CBL );

            blas::internal::LocalGemm
            ( orientationOfA2, Normal, alpha, A2R, B1L, (T)1, CBL );
        }
        else
        {
            blas::internal::LocalGemm
            ( orientationOfA1, Normal, alpha, A1L, B2R, beta, CTR );

            blas::internal::LocalGemm
            ( orientationOfA2, Normal, alpha, A2L, B1R, (T)1, CTR );
        }

        // Recurse
        blas::internal::LocalTriangularRank2K
        ( shape, orientationOfA1, orientationOfA2,
          alpha, A1L, A2L, B1L, B2L, beta, CTL );

        blas::internal::LocalTriangularRank2K
        ( shape, orientationOfA1, orientationOfA2, 
          alpha, A1R, A2R, B1R, B2R, beta, CBR );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape, 
  Orientation orientationOfB1,
  Orientation orientationOfB2,
  float alpha, const DistMatrix<float,MC,  Star>& A1,
               const DistMatrix<float,MC,  Star>& A2,
               const DistMatrix<float,MR,  Star>& B1,
               const DistMatrix<float,MR,  Star>& B2,
  float beta,        DistMatrix<float,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape, 
  Orientation orientationOfA1,
  Orientation orientationOfB1,
  Orientation orientationOfB2,
  float alpha, const DistMatrix<float,Star,MC  >& A1,
               const DistMatrix<float,MC,  Star>& A2,
               const DistMatrix<float,MR,  Star>& B1,
               const DistMatrix<float,MR,  Star>& B2,
  float beta,        DistMatrix<float,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape, 
  Orientation orientationOfA2,
  Orientation orientationOfB1,
  Orientation orientationOfB2,
  float alpha, const DistMatrix<float,MC,  Star>& A1,
               const DistMatrix<float,Star,MC  >& A2,
               const DistMatrix<float,MR,  Star>& B1,
               const DistMatrix<float,MR,  Star>& B2,
  float beta,        DistMatrix<float,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape, 
  Orientation orientationOfB2,
  float alpha, const DistMatrix<float,MC,  Star>& A1,
               const DistMatrix<float,MC,  Star>& A2,
               const DistMatrix<float,Star,MR  >& B1,
               const DistMatrix<float,MR,  Star>& B2,
  float beta,        DistMatrix<float,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape, 
  Orientation orientationOfB1,
  float alpha, const DistMatrix<float,MC,  Star>& A1,
               const DistMatrix<float,MC,  Star>& A2,
               const DistMatrix<float,MR,  Star>& B1,
               const DistMatrix<float,Star,MR  >& B2,
  float beta,        DistMatrix<float,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape, 
  Orientation orientationOfA1,
  Orientation orientationOfA2,
  Orientation orientationOfB1,
  Orientation orientationOfB2,
  float alpha, const DistMatrix<float,Star,MC  >& A1,
               const DistMatrix<float,Star,MC  >& A2,
               const DistMatrix<float,MR,  Star>& B1,
               const DistMatrix<float,MR,  Star>& B2,
  float beta,        DistMatrix<float,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape, 
  Orientation orientationOfA1,
  Orientation orientationOfB2,
  float alpha, const DistMatrix<float,Star,MC  >& A1,
               const DistMatrix<float,MC,  Star>& A2,
               const DistMatrix<float,Star,MR  >& B1,
               const DistMatrix<float,MR,  Star>& B2,
  float beta,        DistMatrix<float,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape, 
  Orientation orientationOfA1,
  Orientation orientationOfB1,
  float alpha, const DistMatrix<float,Star,MC  >& A1,
               const DistMatrix<float,MC,  Star>& A2,
               const DistMatrix<float,MR,  Star>& B1,
               const DistMatrix<float,Star,MR  >& B2,
  float beta,        DistMatrix<float,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape, 
  Orientation orientationOfA2,
  Orientation orientationOfB2,
  float alpha, const DistMatrix<float,MC,  Star>& A1,
               const DistMatrix<float,Star,MC  >& A2,
               const DistMatrix<float,Star,MR  >& B1,
               const DistMatrix<float,MR,  Star>& B2,
  float beta,        DistMatrix<float,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape, 
  Orientation orientationOfA2,
  Orientation orientationOfB1,
  float alpha, const DistMatrix<float,MC,  Star>& A1,
               const DistMatrix<float,Star,MC  >& A2,
               const DistMatrix<float,MR,  Star>& B1,
               const DistMatrix<float,Star,MR  >& B2,
  float beta,        DistMatrix<float,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape, 
  float alpha, const DistMatrix<float,MC,  Star>& A1,
               const DistMatrix<float,MC,  Star>& A2,
               const DistMatrix<float,Star,MR  >& B1,
               const DistMatrix<float,Star,MR  >& B2,
  float beta,        DistMatrix<float,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape, 
  Orientation orientationOfA1,
  Orientation orientationOfA2,
  Orientation orientationOfB2,
  float alpha, const DistMatrix<float,Star,MC  >& A1,
               const DistMatrix<float,Star,MC  >& A2,
               const DistMatrix<float,Star,MR  >& B1,
               const DistMatrix<float,MR,  Star>& B2,
  float beta,        DistMatrix<float,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape, 
  Orientation orientationOfA1,
  Orientation orientationOfA2,
  Orientation orientationOfB1,
  float alpha, const DistMatrix<float,Star,MC  >& A1,
               const DistMatrix<float,Star,MC  >& A2,
               const DistMatrix<float,MR,  Star>& B1,
               const DistMatrix<float,Star,MR  >& B2,
  float beta,        DistMatrix<float,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape, 
  Orientation orientationOfA1,
  float alpha, const DistMatrix<float,Star,MC  >& A1,
               const DistMatrix<float,MC,  Star>& A2,
               const DistMatrix<float,Star,MR  >& B1,
               const DistMatrix<float,Star,MR  >& B2,
  float beta,        DistMatrix<float,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape, 
  Orientation orientationOfA2,
  float alpha, const DistMatrix<float,MC,  Star>& A1,
               const DistMatrix<float,Star,MC  >& A2,
               const DistMatrix<float,Star,MR  >& B1,
               const DistMatrix<float,Star,MR  >& B2,
  float beta,        DistMatrix<float,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape, 
  Orientation orientationOfA1,
  Orientation orientationOfA2,
  float alpha, const DistMatrix<float,Star,MC  >& A1,
               const DistMatrix<float,Star,MC  >& A2,
               const DistMatrix<float,Star,MR  >& B1,
               const DistMatrix<float,Star,MR  >& B2,
  float beta,        DistMatrix<float,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape, 
  Orientation orientationOfB1,
  Orientation orientationOfB2,
  double alpha, const DistMatrix<double,MC,  Star>& A1,
                const DistMatrix<double,MC,  Star>& A2,
                const DistMatrix<double,MR,  Star>& B1,
                const DistMatrix<double,MR,  Star>& B2,
  double beta,        DistMatrix<double,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape, 
  Orientation orientationOfA1,
  Orientation orientationOfB1,
  Orientation orientationOfB2,
  double alpha, const DistMatrix<double,Star,MC  >& A1,
                const DistMatrix<double,MC,  Star>& A2,
                const DistMatrix<double,MR,  Star>& B1,
                const DistMatrix<double,MR,  Star>& B2,
  double beta,        DistMatrix<double,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape, 
  Orientation orientationOfA2,
  Orientation orientationOfB1,
  Orientation orientationOfB2,
  double alpha, const DistMatrix<double,MC,  Star>& A1,
                const DistMatrix<double,Star,MC  >& A2,
                const DistMatrix<double,MR,  Star>& B1,
                const DistMatrix<double,MR,  Star>& B2,
  double beta,        DistMatrix<double,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape, 
  Orientation orientationOfB2,
  double alpha, const DistMatrix<double,MC,  Star>& A1,
                const DistMatrix<double,MC,  Star>& A2,
                const DistMatrix<double,Star,MR  >& B1,
                const DistMatrix<double,MR,  Star>& B2,
  double beta,        DistMatrix<double,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape, 
  Orientation orientationOfB1,
  double alpha, const DistMatrix<double,MC,  Star>& A1,
                const DistMatrix<double,MC,  Star>& A2,
                const DistMatrix<double,MR,  Star>& B1,
                const DistMatrix<double,Star,MR  >& B2,
  double beta,        DistMatrix<double,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape, 
  Orientation orientationOfA1,
  Orientation orientationOfA2,
  Orientation orientationOfB1,
  Orientation orientationOfB2,
  double alpha, const DistMatrix<double,Star,MC  >& A1,
                const DistMatrix<double,Star,MC  >& A2,
                const DistMatrix<double,MR,  Star>& B1,
                const DistMatrix<double,MR,  Star>& B2,
  double beta,        DistMatrix<double,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape, 
  Orientation orientationOfA1,
  Orientation orientationOfB2,
  double alpha, const DistMatrix<double,Star,MC  >& A1,
                const DistMatrix<double,MC,  Star>& A2,
                const DistMatrix<double,Star,MR  >& B1,
                const DistMatrix<double,MR,  Star>& B2,
  double beta,        DistMatrix<double,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape, 
  Orientation orientationOfA1,
  Orientation orientationOfB1,
  double alpha, const DistMatrix<double,Star,MC  >& A1,
                const DistMatrix<double,MC,  Star>& A2,
                const DistMatrix<double,MR,  Star>& B1,
                const DistMatrix<double,Star,MR  >& B2,
  double beta,        DistMatrix<double,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape, 
  Orientation orientationOfA2,
  Orientation orientationOfB2,
  double alpha, const DistMatrix<double,MC,  Star>& A1,
                const DistMatrix<double,Star,MC  >& A2,
                const DistMatrix<double,Star,MR  >& B1,
                const DistMatrix<double,MR,  Star>& B2,
  double beta,        DistMatrix<double,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape, 
  Orientation orientationOfA2,
  Orientation orientationOfB1,
  double alpha, const DistMatrix<double,MC,  Star>& A1,
                const DistMatrix<double,Star,MC  >& A2,
                const DistMatrix<double,MR,  Star>& B1,
                const DistMatrix<double,Star,MR  >& B2,
  double beta,        DistMatrix<double,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape, 
  double alpha, const DistMatrix<double,MC,  Star>& A1,
                const DistMatrix<double,MC,  Star>& A2,
                const DistMatrix<double,Star,MR  >& B1,
                const DistMatrix<double,Star,MR  >& B2,
  double beta,        DistMatrix<double,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape, 
  Orientation orientationOfA1,
  Orientation orientationOfA2,
  Orientation orientationOfB2,
  double alpha, const DistMatrix<double,Star,MC  >& A1,
                const DistMatrix<double,Star,MC  >& A2,
                const DistMatrix<double,Star,MR  >& B1,
                const DistMatrix<double,MR,  Star>& B2,
  double beta,        DistMatrix<double,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape, 
  Orientation orientationOfA1,
  Orientation orientationOfA2,
  Orientation orientationOfB1,
  double alpha, const DistMatrix<double,Star,MC  >& A1,
                const DistMatrix<double,Star,MC  >& A2,
                const DistMatrix<double,MR,  Star>& B1,
                const DistMatrix<double,Star,MR  >& B2,
  double beta,        DistMatrix<double,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape, 
  Orientation orientationOfA1,
  double alpha, const DistMatrix<double,Star,MC  >& A1,
                const DistMatrix<double,MC,  Star>& A2,
                const DistMatrix<double,Star,MR  >& B1,
                const DistMatrix<double,Star,MR  >& B2,
  double beta,        DistMatrix<double,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape, 
  Orientation orientationOfA2,
  double alpha, const DistMatrix<double,MC,  Star>& A1,
                const DistMatrix<double,Star,MC  >& A2,
                const DistMatrix<double,Star,MR  >& B1,
                const DistMatrix<double,Star,MR  >& B2,
  double beta,        DistMatrix<double,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape, 
  Orientation orientationOfA1,
  Orientation orientationOfA2,
  double alpha, const DistMatrix<double,Star,MC  >& A1,
                const DistMatrix<double,Star,MC  >& A2,
                const DistMatrix<double,Star,MR  >& B1,
                const DistMatrix<double,Star,MR  >& B2,
  double beta,        DistMatrix<double,MC,  MR  >& C );

#ifndef WITHOUT_COMPLEX
template void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape, 
  Orientation orientationOfB1,
  Orientation orientationOfB2,
  scomplex alpha, const DistMatrix<scomplex,MC,  Star>& A1,
                  const DistMatrix<scomplex,MC,  Star>& A2,
                  const DistMatrix<scomplex,MR,  Star>& B1,
                  const DistMatrix<scomplex,MR,  Star>& B2,
  scomplex beta,        DistMatrix<scomplex,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape, 
  Orientation orientationOfA1,
  Orientation orientationOfB1,
  Orientation orientationOfB2,
  scomplex alpha, const DistMatrix<scomplex,Star,MC  >& A1,
                  const DistMatrix<scomplex,MC,  Star>& A2,
                  const DistMatrix<scomplex,MR,  Star>& B1,
                  const DistMatrix<scomplex,MR,  Star>& B2,
  scomplex beta,        DistMatrix<scomplex,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape, 
  Orientation orientationOfA2,
  Orientation orientationOfB1,
  Orientation orientationOfB2,
  scomplex alpha, const DistMatrix<scomplex,MC,  Star>& A1,
                  const DistMatrix<scomplex,Star,MC  >& A2,
                  const DistMatrix<scomplex,MR,  Star>& B1,
                  const DistMatrix<scomplex,MR,  Star>& B2,
  scomplex beta,        DistMatrix<scomplex,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape, 
  Orientation orientationOfB2,
  scomplex alpha, const DistMatrix<scomplex,MC,  Star>& A1,
                  const DistMatrix<scomplex,MC,  Star>& A2,
                  const DistMatrix<scomplex,Star,MR  >& B1,
                  const DistMatrix<scomplex,MR,  Star>& B2,
  scomplex beta,        DistMatrix<scomplex,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape, 
  Orientation orientationOfB1,
  scomplex alpha, const DistMatrix<scomplex,MC,  Star>& A1,
                  const DistMatrix<scomplex,MC,  Star>& A2,
                  const DistMatrix<scomplex,MR,  Star>& B1,
                  const DistMatrix<scomplex,Star,MR  >& B2,
  scomplex beta,        DistMatrix<scomplex,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape, 
  Orientation orientationOfA1,
  Orientation orientationOfA2,
  Orientation orientationOfB1,
  Orientation orientationOfB2,
  scomplex alpha, const DistMatrix<scomplex,Star,MC  >& A1,
                  const DistMatrix<scomplex,Star,MC  >& A2,
                  const DistMatrix<scomplex,MR,  Star>& B1,
                  const DistMatrix<scomplex,MR,  Star>& B2,
  scomplex beta,        DistMatrix<scomplex,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape, 
  Orientation orientationOfA1,
  Orientation orientationOfB2,
  scomplex alpha, const DistMatrix<scomplex,Star,MC  >& A1,
                  const DistMatrix<scomplex,MC,  Star>& A2,
                  const DistMatrix<scomplex,Star,MR  >& B1,
                  const DistMatrix<scomplex,MR,  Star>& B2,
  scomplex beta,        DistMatrix<scomplex,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape, 
  Orientation orientationOfA1,
  Orientation orientationOfB1,
  scomplex alpha, const DistMatrix<scomplex,Star,MC  >& A1,
                  const DistMatrix<scomplex,MC,  Star>& A2,
                  const DistMatrix<scomplex,MR,  Star>& B1,
                  const DistMatrix<scomplex,Star,MR  >& B2,
  scomplex beta,        DistMatrix<scomplex,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape, 
  Orientation orientationOfA2,
  Orientation orientationOfB2,
  scomplex alpha, const DistMatrix<scomplex,MC,  Star>& A1,
                  const DistMatrix<scomplex,Star,MC  >& A2,
                  const DistMatrix<scomplex,Star,MR  >& B1,
                  const DistMatrix<scomplex,MR,  Star>& B2,
  scomplex beta,        DistMatrix<scomplex,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape, 
  Orientation orientationOfA2,
  Orientation orientationOfB1,
  scomplex alpha, const DistMatrix<scomplex,MC,  Star>& A1,
                  const DistMatrix<scomplex,Star,MC  >& A2,
                  const DistMatrix<scomplex,MR,  Star>& B1,
                  const DistMatrix<scomplex,Star,MR  >& B2,
  scomplex beta,        DistMatrix<scomplex,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape, 
  scomplex alpha, const DistMatrix<scomplex,MC,  Star>& A1,
                  const DistMatrix<scomplex,MC,  Star>& A2,
                  const DistMatrix<scomplex,Star,MR  >& B1,
                  const DistMatrix<scomplex,Star,MR  >& B2,
  scomplex beta,        DistMatrix<scomplex,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape, 
  Orientation orientationOfA1,
  Orientation orientationOfA2,
  Orientation orientationOfB2,
  scomplex alpha, const DistMatrix<scomplex,Star,MC  >& A1,
                  const DistMatrix<scomplex,Star,MC  >& A2,
                  const DistMatrix<scomplex,Star,MR  >& B1,
                  const DistMatrix<scomplex,MR,  Star>& B2,
  scomplex beta,        DistMatrix<scomplex,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape, 
  Orientation orientationOfA1,
  Orientation orientationOfA2,
  Orientation orientationOfB1,
  scomplex alpha, const DistMatrix<scomplex,Star,MC  >& A1,
                  const DistMatrix<scomplex,Star,MC  >& A2,
                  const DistMatrix<scomplex,MR,  Star>& B1,
                  const DistMatrix<scomplex,Star,MR  >& B2,
  scomplex beta,        DistMatrix<scomplex,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape, 
  Orientation orientationOfA1,
  scomplex alpha, const DistMatrix<scomplex,Star,MC  >& A1,
                  const DistMatrix<scomplex,MC,  Star>& A2,
                  const DistMatrix<scomplex,Star,MR  >& B1,
                  const DistMatrix<scomplex,Star,MR  >& B2,
  scomplex beta,        DistMatrix<scomplex,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape, 
  Orientation orientationOfA2,
  scomplex alpha, const DistMatrix<scomplex,MC,  Star>& A1,
                  const DistMatrix<scomplex,Star,MC  >& A2,
                  const DistMatrix<scomplex,Star,MR  >& B1,
                  const DistMatrix<scomplex,Star,MR  >& B2,
  scomplex beta,        DistMatrix<scomplex,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape, 
  Orientation orientationOfA1,
  Orientation orientationOfA2,
  scomplex alpha, const DistMatrix<scomplex,Star,MC  >& A1,
                  const DistMatrix<scomplex,Star,MC  >& A2,
                  const DistMatrix<scomplex,Star,MR  >& B1,
                  const DistMatrix<scomplex,Star,MR  >& B2,
  scomplex beta,        DistMatrix<scomplex,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape, 
  Orientation orientationOfB1,
  Orientation orientationOfB2,
  dcomplex alpha, const DistMatrix<dcomplex,MC,  Star>& A1,
                  const DistMatrix<dcomplex,MC,  Star>& A2,
                  const DistMatrix<dcomplex,MR,  Star>& B1,
                  const DistMatrix<dcomplex,MR,  Star>& B2,
  dcomplex beta,        DistMatrix<dcomplex,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape, 
  Orientation orientationOfA1,
  Orientation orientationOfB1,
  Orientation orientationOfB2,
  dcomplex alpha, const DistMatrix<dcomplex,Star,MC  >& A1,
                  const DistMatrix<dcomplex,MC,  Star>& A2,
                  const DistMatrix<dcomplex,MR,  Star>& B1,
                  const DistMatrix<dcomplex,MR,  Star>& B2,
  dcomplex beta,        DistMatrix<dcomplex,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape, 
  Orientation orientationOfA2,
  Orientation orientationOfB1,
  Orientation orientationOfB2,
  dcomplex alpha, const DistMatrix<dcomplex,MC,  Star>& A1,
                  const DistMatrix<dcomplex,Star,MC  >& A2,
                  const DistMatrix<dcomplex,MR,  Star>& B1,
                  const DistMatrix<dcomplex,MR,  Star>& B2,
  dcomplex beta,        DistMatrix<dcomplex,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape, 
  Orientation orientationOfB2,
  dcomplex alpha, const DistMatrix<dcomplex,MC,  Star>& A1,
                  const DistMatrix<dcomplex,MC,  Star>& A2,
                  const DistMatrix<dcomplex,Star,MR  >& B1,
                  const DistMatrix<dcomplex,MR,  Star>& B2,
  dcomplex beta,        DistMatrix<dcomplex,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape, 
  Orientation orientationOfB1,
  dcomplex alpha, const DistMatrix<dcomplex,MC,  Star>& A1,
                  const DistMatrix<dcomplex,MC,  Star>& A2,
                  const DistMatrix<dcomplex,MR,  Star>& B1,
                  const DistMatrix<dcomplex,Star,MR  >& B2,
  dcomplex beta,        DistMatrix<dcomplex,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape, 
  Orientation orientationOfA1,
  Orientation orientationOfA2,
  Orientation orientationOfB1,
  Orientation orientationOfB2,
  dcomplex alpha, const DistMatrix<dcomplex,Star,MC  >& A1,
                  const DistMatrix<dcomplex,Star,MC  >& A2,
                  const DistMatrix<dcomplex,MR,  Star>& B1,
                  const DistMatrix<dcomplex,MR,  Star>& B2,
  dcomplex beta,        DistMatrix<dcomplex,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape, 
  Orientation orientationOfA1,
  Orientation orientationOfB2,
  dcomplex alpha, const DistMatrix<dcomplex,Star,MC  >& A1,
                  const DistMatrix<dcomplex,MC,  Star>& A2,
                  const DistMatrix<dcomplex,Star,MR  >& B1,
                  const DistMatrix<dcomplex,MR,  Star>& B2,
  dcomplex beta,        DistMatrix<dcomplex,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape, 
  Orientation orientationOfA1,
  Orientation orientationOfB1,
  dcomplex alpha, const DistMatrix<dcomplex,Star,MC  >& A1,
                  const DistMatrix<dcomplex,MC,  Star>& A2,
                  const DistMatrix<dcomplex,MR,  Star>& B1,
                  const DistMatrix<dcomplex,Star,MR  >& B2,
  dcomplex beta,        DistMatrix<dcomplex,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape, 
  Orientation orientationOfA2,
  Orientation orientationOfB2,
  dcomplex alpha, const DistMatrix<dcomplex,MC,  Star>& A1,
                  const DistMatrix<dcomplex,Star,MC  >& A2,
                  const DistMatrix<dcomplex,Star,MR  >& B1,
                  const DistMatrix<dcomplex,MR,  Star>& B2,
  dcomplex beta,        DistMatrix<dcomplex,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape, 
  Orientation orientationOfA2,
  Orientation orientationOfB1,
  dcomplex alpha, const DistMatrix<dcomplex,MC,  Star>& A1,
                  const DistMatrix<dcomplex,Star,MC  >& A2,
                  const DistMatrix<dcomplex,MR,  Star>& B1,
                  const DistMatrix<dcomplex,Star,MR  >& B2,
  dcomplex beta,        DistMatrix<dcomplex,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape, 
  dcomplex alpha, const DistMatrix<dcomplex,MC,  Star>& A1,
                  const DistMatrix<dcomplex,MC,  Star>& A2,
                  const DistMatrix<dcomplex,Star,MR  >& B1,
                  const DistMatrix<dcomplex,Star,MR  >& B2,
  dcomplex beta,        DistMatrix<dcomplex,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape, 
  Orientation orientationOfA1,
  Orientation orientationOfA2,
  Orientation orientationOfB2,
  dcomplex alpha, const DistMatrix<dcomplex,Star,MC  >& A1,
                  const DistMatrix<dcomplex,Star,MC  >& A2,
                  const DistMatrix<dcomplex,Star,MR  >& B1,
                  const DistMatrix<dcomplex,MR,  Star>& B2,
  dcomplex beta,        DistMatrix<dcomplex,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape, 
  Orientation orientationOfA1,
  Orientation orientationOfA2,
  Orientation orientationOfB1,
  dcomplex alpha, const DistMatrix<dcomplex,Star,MC  >& A1,
                  const DistMatrix<dcomplex,Star,MC  >& A2,
                  const DistMatrix<dcomplex,MR,  Star>& B1,
                  const DistMatrix<dcomplex,Star,MR  >& B2,
  dcomplex beta,        DistMatrix<dcomplex,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape, 
  Orientation orientationOfA1,
  dcomplex alpha, const DistMatrix<dcomplex,Star,MC  >& A1,
                  const DistMatrix<dcomplex,MC,  Star>& A2,
                  const DistMatrix<dcomplex,Star,MR  >& B1,
                  const DistMatrix<dcomplex,Star,MR  >& B2,
  dcomplex beta,        DistMatrix<dcomplex,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape, 
  Orientation orientationOfA2,
  dcomplex alpha, const DistMatrix<dcomplex,MC,  Star>& A1,
                  const DistMatrix<dcomplex,Star,MC  >& A2,
                  const DistMatrix<dcomplex,Star,MR  >& B1,
                  const DistMatrix<dcomplex,Star,MR  >& B2,
  dcomplex beta,        DistMatrix<dcomplex,MC,  MR  >& C );

template void
elemental::blas::internal::LocalTriangularRank2K
( Shape shape, 
  Orientation orientationOfA1,
  Orientation orientationOfA2,
  dcomplex alpha, const DistMatrix<dcomplex,Star,MC  >& A1,
                  const DistMatrix<dcomplex,Star,MC  >& A2,
                  const DistMatrix<dcomplex,Star,MR  >& B1,
                  const DistMatrix<dcomplex,Star,MR  >& B2,
  dcomplex beta,        DistMatrix<dcomplex,MC,  MR  >& C );
#endif

