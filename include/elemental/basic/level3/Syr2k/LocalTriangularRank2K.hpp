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

namespace {
using namespace elemental;

#ifndef RELEASE
template<typename T>
inline void 
CheckInput
( Orientation orientationOfB1,
  Orientation orientationOfB2,
  const DistMatrix<T,MC,STAR>& A1, 
  const DistMatrix<T,MC,STAR>& A2,
  const DistMatrix<T,MR,STAR>& B1,
  const DistMatrix<T,MR,STAR>& B2,
  const DistMatrix<T,MC,MR  >& C )
{
    if( orientationOfB1 == NORMAL )
        throw std::logic_error("B1[MR,* ] must be (Conjugate)Transpose'd");
    if( orientationOfB2 == NORMAL )
        throw std::logic_error("B2[MR,* ] must be (Conjugate)Transpose'd");
    if( A1.Grid() != A2.Grid() || A2.Grid() != B1.Grid() ||
        B1.Grid() != B2.Grid() || B2.Grid() != C.Grid() )
        throw std::logic_error
        ("A, B, and C must be distributed over the same grid");
    if( A1.Height() != C.Height() || B1.Height() != C.Width() ||
        A1.Height() != A2.Height() || A1.Width() != A2.Width() ||
        B1.Width() != B2.Width() || B1.Height() != B2.Height() ||
        A1.Width() != B1.Width() )
    {
        std::ostringstream msg;
        msg << "Nonconformal LocalTriangularRank2K: \n"
            << "  A1[MC,* ] ~ " << A1.Height() << " x "
                                << A1.Width()  << "\n"
            << "  A2[MC,* ] ~ " << A2.Height() << " x "
                                << A2.Width()  << "\n"
            << "  B1[MR,* ] ~ " << B1.Height() << " x "
                                << B1.Width()  << "\n"
            << "  B2[MR,* ] ~ " << B2.Height() << " x "
                                << B2.Width()  << "\n"
            << "  C[MC,MR] ~ " << C.Height() << " x " << C.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( A1.ColAlignment() != C.ColAlignment() ||
        B1.ColAlignment() != C.RowAlignment() ||
        A1.ColAlignment() != A2.ColAlignment() ||
        B1.ColAlignment() != B2.ColAlignment() )
    {
        std::ostringstream msg;
        msg << "Misaligned LocalTriangularRank2K: \n"
            << "  A1[MC,* ] ~ " << A1.ColAlignment() << "\n"
            << "  A2[MC,* ] ~ " << A2.ColAlignment() << "\n"
            << "  B1[MR,* ] ~ " << B1.ColAlignment() << "\n"
            << "  B2[MR,* ] ~ " << B2.ColAlignment() << "\n"
            << "  C[MC,MR] ~ " << C.ColAlignment() << " , " <<
                                  C.RowAlignment() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
}

template<typename T>
inline void 
CheckInput
( Orientation orientationOfA1,
  Orientation orientationOfB1,
  Orientation orientationOfB2,
  const DistMatrix<T,STAR,MC  >& A1, 
  const DistMatrix<T,MC,  STAR>& A2,
  const DistMatrix<T,MR,  STAR>& B1,
  const DistMatrix<T,MR,  STAR>& B2,
  const DistMatrix<T,MC,  MR  >& C )
{
    if( orientationOfA1 == NORMAL )
        throw std::logic_error("A1[* ,MC] must be (Conjugate)Transpose'd");
    if( orientationOfB1 == NORMAL )
        throw std::logic_error("B1[MR,* ] must be (Conjugate)Transpose'd");
    if( orientationOfB2 == NORMAL )
        throw std::logic_error("B2[MR,* ] must be (Conjugate)Transpose'd");
    if( A1.Grid() != A2.Grid() || A2.Grid() != B1.Grid() ||
        B1.Grid() != B2.Grid() || B2.Grid() != C.Grid() )
        throw std::logic_error
        ("A, B, and C must be distributed over the same grid");
    if( A1.Width() != C.Height() || B1.Height() != C.Width() ||
        A1.Width() != A2.Height() || A1.Height() != A2.Width() ||
        B1.Width() != B2.Width() || B1.Height() != B2.Height() ||
        A1.Height() != B1.Width() )
    {
        std::ostringstream msg;
        msg << "Nonconformal LocalTriangularRank2K: \n"
            << "  A1[* ,MC] ~ " << A1.Height() << " x "
                                << A1.Width()  << "\n"
            << "  A2[MC,* ] ~ " << A2.Height() << " x "
                                << A2.Width()  << "\n"
            << "  B1[MR,* ] ~ " << B1.Height() << " x "
                                << B1.Width()  << "\n"
            << "  B2[MR,* ] ~ " << B2.Height() << " x "
                                << B2.Width()  << "\n"
            << "  C[MC,MR] ~ " << C.Height() << " x " << C.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( A1.RowAlignment() != C.ColAlignment() ||
        B1.ColAlignment() != C.RowAlignment() ||
        A1.RowAlignment() != A2.ColAlignment() ||
        B1.ColAlignment() != B2.ColAlignment() )
    {
        std::ostringstream msg;
        msg << "Misaligned LocalTriangularRank2K: \n"
            << "  A1[* ,MC] ~ " << A1.RowAlignment() << "\n"
            << "  A2[MC,* ] ~ " << A2.ColAlignment() << "\n"
            << "  B1[MR,* ] ~ " << B1.ColAlignment() << "\n"
            << "  B2[MR,* ] ~ " << B2.ColAlignment() << "\n"
            << "  C[MC,MR] ~ " << C.ColAlignment() << " , " <<
                                  C.RowAlignment() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
}

template<typename T>
inline void 
CheckInput
( Orientation orientationOfA2,
  Orientation orientationOfB1,
  Orientation orientationOfB2,
  const DistMatrix<T,MC,  STAR>& A1, 
  const DistMatrix<T,STAR,MC  >& A2,
  const DistMatrix<T,MR,  STAR>& B1,
  const DistMatrix<T,MR,  STAR>& B2,
  const DistMatrix<T,MC,  MR  >& C )
{
    if( orientationOfA2 == NORMAL )
        throw std::logic_error("A2[* ,MC] must be (Conjugate)Transpose'd");
    if( orientationOfB1 == NORMAL )
        throw std::logic_error("B1[MR,* ] must be (Conjugate)Transpose'd");
    if( orientationOfB2 == NORMAL )
        throw std::logic_error("B2[MR,* ] must be (Conjugate)Transpose'd");
    if( A1.Grid() != A2.Grid() || A2.Grid() != B1.Grid() ||
        B1.Grid() != B2.Grid() || B2.Grid() != C.Grid() )
        throw std::logic_error
        ("A, B, and C must be distributed over the same grid");
    if( A1.Height() != C.Height() || B1.Height() != C.Width() ||
        A1.Height() != A2.Width() || A1.Width() != A2.Height() ||
        B1.Width() != B2.Width() || B1.Height() != B2.Height() ||
        A1.Width() != B1.Width() )
    {
        std::ostringstream msg;
        msg << "Nonconformal LocalTriangularRank2K: \n"
            << "  A1[MC,* ] ~ " << A1.Height() << " x "
                                << A1.Width()  << "\n"
            << "  A2[* ,MC] ~ " << A2.Height() << " x "
                                << A2.Width()  << "\n"
            << "  B1[MR,* ] ~ " << B1.Height() << " x "
                                << B1.Width()  << "\n"
            << "  B2[MR,* ] ~ " << B2.Height() << " x "
                                << B2.Width()  << "\n"
            << "  C[MC,MR] ~ " << C.Height() << " x " << C.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( A1.ColAlignment() != C.ColAlignment() ||
        B1.ColAlignment() != C.RowAlignment() ||
        A1.ColAlignment() != A2.RowAlignment() ||
        B1.ColAlignment() != B2.ColAlignment() )
    {
        std::ostringstream msg;
        msg << "Misaligned LocalTriangularRank2K: \n"
            << "  A1[MC,* ] ~ " << A1.ColAlignment() << "\n"
            << "  A2[* ,MC] ~ " << A2.RowAlignment() << "\n"
            << "  B1[MR,* ] ~ " << B1.ColAlignment() << "\n"
            << "  B2[MR,* ] ~ " << B2.ColAlignment() << "\n"
            << "  C[MC,MR] ~ " << C.ColAlignment() << " , " <<
                                  C.RowAlignment() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
}

template<typename T>
inline void 
CheckInput
( Orientation orientationOfB2,
  const DistMatrix<T,MC,  STAR>& A1, 
  const DistMatrix<T,MC,  STAR>& A2,
  const DistMatrix<T,STAR,MR  >& B1,
  const DistMatrix<T,MR,  STAR>& B2,
  const DistMatrix<T,MC,  MR  >& C )
{
    if( orientationOfB2 == NORMAL )
        throw std::logic_error("B2[MR,* ] must be (Conjugate)Transpose'd");
    if( A1.Grid() != A2.Grid() || A2.Grid() != B1.Grid() ||
        B1.Grid() != B2.Grid() || B2.Grid() != C.Grid() )
        throw std::logic_error
        ("A, B, and C must be distributed over the same grid");
    if( A1.Height() != C.Height() || B1.Height() != C.Width() ||
        A1.Height() != A2.Height() || A1.Width() != A2.Width() ||
        B1.Height() != B2.Width() || B1.Width() != B2.Height() ||
        A1.Width() != B1.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal LocalTriangularRank2K: \n"
            << "  A1[MC,* ] ~ " << A1.Height() << " x "
                                << A1.Width()  << "\n"
            << "  A2[MC,* ] ~ " << A2.Height() << " x "
                                << A2.Width()  << "\n"
            << "  B1[* ,MR] ~ " << B1.Height() << " x "
                                << B1.Width()  << "\n"
            << "  B2[MR,* ] ~ " << B2.Height() << " x "
                                << B2.Width()  << "\n"
            << "  C[MC,MR] ~ " << C.Height() << " x " << C.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( A1.ColAlignment() != C.ColAlignment() ||
        B1.RowAlignment() != C.RowAlignment() ||
        A1.ColAlignment() != A2.ColAlignment() ||
        B1.RowAlignment() != B2.ColAlignment() )
    {
        std::ostringstream msg;
        msg << "Misaligned LocalTriangularRank2K: \n"
            << "  A1[MC,* ] ~ " << A1.ColAlignment() << "\n"
            << "  A2[MC,* ] ~ " << A2.ColAlignment() << "\n"
            << "  B1[* ,MR] ~ " << B1.RowAlignment() << "\n"
            << "  B2[MR,* ] ~ " << B2.ColAlignment() << "\n"
            << "  C[MC,MR] ~ " << C.ColAlignment() << " , " <<
                                  C.RowAlignment() << "\n";
        throw std::logic_error( msg.str() );
    }
}

template<typename T>
inline void 
CheckInput
( Orientation orientationOfB1,
  const DistMatrix<T,MC,  STAR>& A1, 
  const DistMatrix<T,MC,  STAR>& A2,
  const DistMatrix<T,MR,  STAR>& B1,
  const DistMatrix<T,STAR,MR  >& B2,
  const DistMatrix<T,MC,  MR  >& C )
{
    if( orientationOfB1 == NORMAL )
        throw std::logic_error("B1[MR,* ] must be (Conjugate)Transpose'd");
    if( A1.Grid() != A2.Grid() || A2.Grid() != B1.Grid() ||
        B1.Grid() != B2.Grid() || B2.Grid() != C.Grid() )
        throw std::logic_error
        ("A, B, and C must be distributed over the same grid");
    if( A1.Height() != C.Height() || B1.Height() != C.Width() ||
        A1.Height() != A2.Height() || A1.Width() != A2.Width() ||
        B1.Width() != B2.Height() || B1.Height() != B2.Width() ||
        A1.Width() != B1.Width() )
    {
        std::ostringstream msg;
        msg << "Nonconformal LocalTriangularRank2K: \n"
            << "  A1[MC,* ] ~ " << A1.Height() << " x "
                                << A1.Width()  << "\n"
            << "  A2[MC,* ] ~ " << A2.Height() << " x "
                                << A2.Width()  << "\n"
            << "  B1[MR,* ] ~ " << B1.Height() << " x "
                                << B1.Width()  << "\n"
            << "  B2[* ,MR] ~ " << B2.Height() << " x "
                                << B2.Width()  << "\n"
            << "  C[MC,MR] ~ " << C.Height() << " x " << C.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( A1.ColAlignment() != C.ColAlignment() ||
        B1.ColAlignment() != C.RowAlignment() ||
        A1.ColAlignment() != A2.ColAlignment() ||
        B1.ColAlignment() != B2.RowAlignment() )
    {
        std::ostringstream msg;
        msg << "Misaligned LocalTriangularRank2K: \n"
            << "  A1[MC,* ] ~ " << A1.ColAlignment() << "\n"
            << "  A2[MC,* ] ~ " << A2.ColAlignment() << "\n"
            << "  B1[MR,* ] ~ " << B1.ColAlignment() << "\n"
            << "  B2[* ,MR] ~ " << B2.RowAlignment() << "\n"
            << "  C[MC,MR] ~ " << C.ColAlignment() << " , " <<
                                  C.RowAlignment() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
}

template<typename T>
inline void 
CheckInput
( Orientation orientationOfA1,
  Orientation orientationOfA2,
  Orientation orientationOfB1,
  Orientation orientationOfB2,
  const DistMatrix<T,STAR,MC  >& A1, 
  const DistMatrix<T,STAR,MC  >& A2,
  const DistMatrix<T,MR,  STAR>& B1,
  const DistMatrix<T,MR,  STAR>& B2,
  const DistMatrix<T,MC,  MR  >& C )
{
    if( orientationOfA1 == NORMAL )
        throw std::logic_error("A1[* ,MC] must be (Conjugate)Transpose'd");
    if( orientationOfA2 == NORMAL )
        throw std::logic_error("A2[* ,MC] must be (Conjugate)Transpose'd");
    if( orientationOfB1 == NORMAL )
        throw std::logic_error("B1[MR,* ] must be (Conjugate)Transpose'd");
    if( orientationOfB2 == NORMAL )
        throw std::logic_error("B2[MR,* ] must be (Conjugate)Transpose'd");
    if( A1.Grid() != A2.Grid() || A2.Grid() != B1.Grid() ||
        B1.Grid() != B2.Grid() || B2.Grid() != C.Grid() )
        throw std::logic_error
        ("A, B, and C must be distributed over the same grid");
    if( A1.Width() != C.Height() || B1.Height() != C.Width() ||
        A1.Height() != A2.Height() || A1.Width() != A2.Width() ||
        B1.Height() != B2.Height() || B1.Width() != B2.Width() ||
        A1.Height() != B1.Width() )
    {
        std::ostringstream msg;
        msg << "Nonconformal LocalTriangularRank2K: \n"
            << "  A1[* ,MC] ~ " << A1.Height() << " x "
                                << A1.Width()  << "\n"
            << "  A2[* ,MC] ~ " << A2.Height() << " x "
                                << A2.Width()  << "\n"
            << "  B1[MR,* ] ~ " << B1.Height() << " x "
                                << B1.Width()  << "\n"
            << "  B2[MR,* ] ~ " << B2.Height() << " x "
                                << B2.Width()  << "\n"
            << "  C[MC,MR] ~ " << C.Height() << " x " << C.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( A1.RowAlignment() != C.ColAlignment() ||
        B1.ColAlignment() != C.RowAlignment() ||
        A1.RowAlignment() != A2.RowAlignment() ||
        B1.ColAlignment() != B2.ColAlignment() )
    {
        std::ostringstream msg;
        msg << "Misaligned LocalTriangularRank2K: \n"
            << "  A1[* ,MC] ~ " << A1.RowAlignment() << "\n"
            << "  A2[* ,MC] ~ " << A2.RowAlignment() << "\n"
            << "  B1[MR,* ] ~ " << B1.ColAlignment() << "\n"
            << "  B2[MR,* ] ~ " << B2.ColAlignment() << "\n"
            << "  C[MC,MR] ~ " << C.ColAlignment() << " , " <<
                                  C.RowAlignment() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
}

template<typename T>
inline void 
CheckInput
( Orientation orientationOfA1,
  Orientation orientationOfB2,
  const DistMatrix<T,STAR,MC  >& A1, 
  const DistMatrix<T,MC,  STAR>& A2,
  const DistMatrix<T,STAR,MR  >& B1,
  const DistMatrix<T,MR,  STAR>& B2,
  const DistMatrix<T,MC,  MR  >& C )
{
    if( orientationOfA1 == NORMAL )
        throw std::logic_error("A1[* ,MC] must be (Conjugate)Transpose'd");
    if( orientationOfB2 == NORMAL )
        throw std::logic_error("B2[MR,* ] must be (Conjugate)Transpose'd");
    if( A1.Grid() != A2.Grid() || A2.Grid() != B1.Grid() ||
        B1.Grid() != B2.Grid() || B2.Grid() != C.Grid() )
        throw std::logic_error
        ("A, B, and C must be distributed over the same grid");
    if( A1.Width() != C.Height() || B1.Width() != C.Width() ||
        A1.Width() != A2.Height() || A1.Height() != A2.Width() ||
        B1.Height() != B2.Width() || B1.Width() != B2.Height() ||
        A1.Height() != B1.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal LocalTriangularRank2K: \n"
            << "  A1[* ,MC] ~ " << A1.Height() << " x "
                                << A1.Width()  << "\n"
            << "  A2[MC,* ] ~ " << A2.Height() << " x "
                                << A2.Width()  << "\n"
            << "  B1[* ,MR] ~ " << B1.Height() << " x "
                                << B1.Width()  << "\n"
            << "  B2[MR,* ] ~ " << B2.Height() << " x "
                                << B2.Width()  << "\n"
            << "  C[MC,MR] ~ " << C.Height() << " x " << C.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( A1.RowAlignment() != C.ColAlignment() ||
        B1.RowAlignment() != C.RowAlignment() ||
        A1.RowAlignment() != A2.ColAlignment() ||
        B1.RowAlignment() != B2.ColAlignment() )
    {
        std::ostringstream msg;
        msg << "Misaligned LocalTriangularRank2K: \n"
            << "  A1[* ,MC] ~ " << A1.RowAlignment() << "\n"
            << "  A2[MC,* ] ~ " << A2.ColAlignment() << "\n"
            << "  B1[* ,MR] ~ " << B1.RowAlignment() << "\n"
            << "  B2[MR,* ] ~ " << B2.ColAlignment() << "\n"
            << "  C[MC,MR] ~ " << C.ColAlignment() << " , " <<
                                  C.RowAlignment() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
}

template<typename T>
inline void 
CheckInput
( Orientation orientationOfA1,
  Orientation orientationOfB1,
  const DistMatrix<T,STAR,MC  >& A1, 
  const DistMatrix<T,MC,  STAR>& A2,
  const DistMatrix<T,MR,  STAR>& B1,
  const DistMatrix<T,STAR,MR  >& B2,
  const DistMatrix<T,MC,  MR  >& C )
{
    if( orientationOfA1 == NORMAL )
        throw std::logic_error("A1[* ,MC] must be (Conjugate)Transpose'd");
    if( orientationOfB1 == NORMAL )
        throw std::logic_error("B1[MR,* ] must be (Conjugate)Transpose'd");
    if( A1.Grid() != A2.Grid() || A2.Grid() != B1.Grid() ||
        B1.Grid() != B2.Grid() || B2.Grid() != C.Grid() )
        throw std::logic_error
        ("A, B, and C must be distributed over the same grid");
    if( A1.Width() != C.Height() || B1.Height() != C.Width() ||
        A1.Width() != A2.Height() || A1.Height() != A2.Width() ||
        B1.Width() != B2.Height() || B1.Height() != B2.Width() ||
        A1.Height() != B1.Width() )
    {
        std::ostringstream msg;
        msg << "Nonconformal LocalTriangularRank2K: \n"
            << "  A1[* ,MC] ~ " << A1.Height() << " x "
                                << A1.Width()  << "\n"
            << "  A2[MC,* ] ~ " << A2.Height() << " x "
                                << A2.Width()  << "\n"
            << "  B1[MR,* ] ~ " << B1.Height() << " x "
                                << B1.Width()  << "\n"
            << "  B2[* ,MR] ~ " << B2.Height() << " x "
                                << B2.Width()  << "\n"
            << "  C[MC,MR] ~ " << C.Height() << " x " << C.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( A1.RowAlignment() != C.ColAlignment() ||
        B1.ColAlignment() != C.RowAlignment() ||
        A1.RowAlignment() != A2.ColAlignment() ||
        B1.ColAlignment() != B2.RowAlignment() )
    {
        std::ostringstream msg;
        msg << "Misaligned LocalTriangularRank2K: \n"
            << "  A1[* ,MC] ~ " << A1.RowAlignment() << "\n"
            << "  A2[MC,* ] ~ " << A2.ColAlignment() << "\n"
            << "  B1[MR,* ] ~ " << B1.ColAlignment() << "\n"
            << "  B2[* ,MR] ~ " << B2.RowAlignment() << "\n"
            << "  C[MC,MR] ~ " << C.ColAlignment() << " , " <<
                                  C.RowAlignment() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
}

template<typename T>
inline void 
CheckInput
( Orientation orientationOfA2,
  Orientation orientationOfB2,
  const DistMatrix<T,MC,  STAR>& A1,
  const DistMatrix<T,STAR,MC  >& A2, 
  const DistMatrix<T,STAR,MR  >& B1,
  const DistMatrix<T,MR,  STAR>& B2,
  const DistMatrix<T,MC,  MR  >& C )
{
    if( orientationOfA2 == NORMAL )
        throw std::logic_error("A2[* ,MC] must be (Conjugate)Transpose'd");
    if( orientationOfB2 == NORMAL )
        throw std::logic_error("B2[MR,* ] must be (Conjugate)Transpose'd");
    if( A1.Grid() != A2.Grid() || A2.Grid() != B1.Grid() ||
        B1.Grid() != B2.Grid() || B2.Grid() != C.Grid() )
        throw std::logic_error
        ("A, B, and C must be distributed over the same grid");
    if( A1.Height() != C.Height() || B1.Width() != C.Width() ||
        A1.Height() != A2.Width() || A1.Width() != A2.Height() ||
        B1.Height() != B2.Width() || B1.Width() != B2.Height() ||
        A1.Width() != B1.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal LocalTriangularRank2K: \n"
            << "  A1[MC,* ] ~ " << A1.Height() << " x "
                                << A1.Width()  << "\n"
            << "  A2[* ,MC] ~ " << A2.Height() << " x "
                                << A2.Width()  << "\n"
            << "  B1[* ,MR] ~ " << B1.Height() << " x "
                                << B1.Width()  << "\n"
            << "  B2[MR,* ] ~ " << B2.Height() << " x "
                                << B2.Width()  << "\n"
            << "  C[MC,MR] ~ " << C.Height() << " x " << C.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( A1.ColAlignment() != C.ColAlignment() ||
        B1.RowAlignment() != C.RowAlignment() ||
        A1.ColAlignment() != A2.RowAlignment() ||
        B1.RowAlignment() != B2.ColAlignment() )
    {
        std::ostringstream msg;
        msg << "Misaligned LocalTriangularRank2K: \n"
            << "  A1[MC,* ] ~ " << A1.ColAlignment() << "\n"
            << "  A2[* ,MC] ~ " << A2.RowAlignment() << "\n"
            << "  B1[* ,MR] ~ " << B1.RowAlignment() << "\n"
            << "  B2[MR,* ] ~ " << B2.ColAlignment() << "\n"
            << "  C[MC,MR] ~ " << C.ColAlignment() << " , " <<
                                  C.RowAlignment() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
}

template<typename T>
inline void 
CheckInput
( Orientation orientationOfA2,
  Orientation orientationOfB1,
  const DistMatrix<T,MC,  STAR>& A1,
  const DistMatrix<T,STAR,MC  >& A2, 
  const DistMatrix<T,MR,  STAR>& B1,
  const DistMatrix<T,STAR,MR  >& B2,
  const DistMatrix<T,MC,  MR  >& C )
{
    if( orientationOfA2 == NORMAL )
        throw std::logic_error("A2[* ,MC] must be (Conjugate)Transpose'd");
    if( orientationOfB1 == NORMAL )
        throw std::logic_error("B1[MR,* ] must be (Conjugate)Transpose'd");
    if( A1.Grid() != A2.Grid() || A2.Grid() != B1.Grid() ||
        B1.Grid() != B2.Grid() || B2.Grid() != C.Grid() )
        throw std::logic_error
        ("A, B, and C must be distributed over the same grid");
    if( A1.Height() != C.Height() || B1.Height() != C.Width() ||
        A1.Height() != A2.Width() || A1.Width() != A2.Height() ||
        B1.Width() != B2.Height() || B1.Height() != B2.Width() ||
        A1.Width() != B1.Width() )
    {
        std::ostringstream msg;
        msg << "Nonconformal LocalTriangularRank2K: \n"
            << "  A1[MC,* ] ~ " << A1.Height() << " x "
                                << A1.Width()  << "\n"
            << "  A2[* ,MC] ~ " << A2.Height() << " x "
                                << A2.Width()  << "\n"
            << "  B1[MR,* ] ~ " << B1.Height() << " x "
                                << B1.Width()  << "\n"
            << "  B2[* ,MR] ~ " << B2.Height() << " x "
                                << B2.Width()  << "\n"
            << "  C[MC,MR] ~ " << C.Height() << " x " << C.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( A1.ColAlignment() != C.ColAlignment() ||
        B1.ColAlignment() != C.RowAlignment() ||
        A1.ColAlignment() != A2.RowAlignment() ||
        B1.ColAlignment() != B2.RowAlignment() )
    {
        std::ostringstream msg;
        msg << "Misaligned LocalTriangularRank2K: \n"
            << "  A1[MC,* ] ~ " << A1.ColAlignment() << "\n"
            << "  A2[* ,MC] ~ " << A2.RowAlignment() << "\n"
            << "  B1[MR,* ] ~ " << B1.ColAlignment() << "\n"
            << "  B2[* ,MR] ~ " << B2.RowAlignment() << "\n"
            << "  C[MC,MR] ~ " << C.ColAlignment() << " , " <<
                                  C.RowAlignment() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
}

template<typename T>
inline void 
CheckInput
( const DistMatrix<T,MC,  STAR>& A1, 
  const DistMatrix<T,MC,  STAR>& A2,
  const DistMatrix<T,STAR,MR  >& B1,
  const DistMatrix<T,STAR,MR  >& B2,
  const DistMatrix<T,MC,  MR  >& C  )
{
    if( A1.Grid() != A2.Grid() || A2.Grid() != B1.Grid() ||
        B1.Grid() != B2.Grid() || B2.Grid() != C.Grid() )
        throw std::logic_error
        ("A, B, and C must be distributed over the same grid");
    if( A1.Height() != C.Height() || B1.Width() != C.Width() ||
        A1.Height() != A2.Height() || A1.Width() != A2.Width() ||
        B1.Height() != B2.Height() || B1.Width() != B2.Width() ||
        A1.Width() != B1.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal LocalTriangularRank2K: \n"
            << "  A1[MC,* ] ~ " << A1.Height() << " x "
                                << A1.Width()  << "\n"
            << "  A2[MC,* ] ~ " << A2.Height() << " x "
                                << A2.Width()  << "\n"
            << "  B1[* ,MR] ~ " << B1.Height() << " x "
                                << B1.Width()  << "\n"
            << "  B2[* ,MR] ~ " << B2.Height() << " x "
                                << B2.Width()  << "\n"
            << "  C[MC,MR] ~ " << C.Height() << " x " << C.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( A1.ColAlignment() != C.ColAlignment() ||
        B1.RowAlignment() != C.RowAlignment() ||
        A1.ColAlignment() != A2.ColAlignment() ||
        B1.RowAlignment() != B2.RowAlignment() )
    {
        std::ostringstream msg;
        msg << "Misaligned LocalTriangularRank2K: \n"
            << "  A1[MC,* ] ~ " << A1.ColAlignment() << "\n"
            << "  A2[MC,* ] ~ " << A2.ColAlignment() << "\n"
            << "  B1[* ,MR] ~ " << B1.RowAlignment() << "\n"
            << "  B2[* ,MR] ~ " << B2.RowAlignment() << "\n"
            << "  C[MC,MR] ~ " << C.ColAlignment() << " , " <<
                                  C.RowAlignment() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
}

template<typename T>
inline void 
CheckInput
( Orientation orientationOfA1,
  Orientation orientationOfA2,
  Orientation orientationOfB2,
  const DistMatrix<T,STAR,MC  >& A1, 
  const DistMatrix<T,STAR,MC  >& A2,
  const DistMatrix<T,STAR,MR  >& B1,
  const DistMatrix<T,MR,  STAR>& B2,
  const DistMatrix<T,MC,  MR  >& C )
{
    if( orientationOfA1 == NORMAL )
        throw std::logic_error("A1[* ,MC] must be (Conjugate)Transpose'd");
    if( orientationOfA2 == NORMAL )
        throw std::logic_error("A2[* ,MC] must be (Conjugate)Transpose'd");
    if( orientationOfB2 == NORMAL )
        throw std::logic_error("B2[MR,* ] must be (Conjugate)Transpose'd");
    if( A1.Grid() != A2.Grid() || A2.Grid() != B1.Grid() ||
        B1.Grid() != B2.Grid() || B2.Grid() != C.Grid() )
        throw std::logic_error
        ("A, B, and C must be distributed over the same grid");
    if( A1.Width() != C.Height() || B1.Width() != C.Width() ||
        A1.Height() != A2.Height() || A1.Width() != A2.Width() ||
        B1.Width() != B2.Height() || B1.Height() != B2.Width() ||
        A1.Height() != B1.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal LocalTriangularRank2K: \n"
            << "  A1[* ,MC] ~ " << A1.Height() << " x "
                                << A1.Width()  << "\n"
            << "  A2[* ,MC] ~ " << A2.Height() << " x "
                                << A2.Width()  << "\n"
            << "  B1[* ,MR] ~ " << B1.Height() << " x "
                                << B1.Width()  << "\n"
            << "  B2[MR,* ] ~ " << B2.Height() << " x "
                                << B2.Width()  << "\n"
            << "  C[MC,MR] ~ " << C.Height() << " x " << C.Width() << "\n";
        throw std::logic_error( msg.str() );
    }
    if( A1.RowAlignment() != C.ColAlignment() ||
        B1.RowAlignment() != C.RowAlignment() ||
        A1.RowAlignment() != A2.RowAlignment() ||
        B1.RowAlignment() != B2.ColAlignment() )
    {
        std::ostringstream msg;
        msg << "Misaligned LocalTriangularRank2K: \n"
            << "  A1[* ,MC] ~ " << A1.RowAlignment() << "\n"
            << "  A2[* ,MC] ~ " << A2.RowAlignment() << "\n"
            << "  B1[* ,MR] ~ " << B1.RowAlignment() << "\n"
            << "  B2[MR,* ] ~ " << B2.ColAlignment() << "\n"
            << "  C[MC,MR] ~ " << C.ColAlignment() << " , " <<
                                  C.RowAlignment() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
}

template<typename T>
inline void 
CheckInput
( Orientation orientationOfA1,
  Orientation orientationOfA2,
  Orientation orientationOfB1,
  const DistMatrix<T,STAR,MC  >& A1, 
  const DistMatrix<T,STAR,MC  >& A2,
  const DistMatrix<T,MR,  STAR>& B1,
  const DistMatrix<T,STAR,MR  >& B2,
  const DistMatrix<T,MC,  MR  >& C )
{
    if( orientationOfA1 == NORMAL )
        throw std::logic_error("A1[* ,MC] must be (Conjugate)Transpose'd");
    if( orientationOfA2 == NORMAL )
        throw std::logic_error("A2[* ,MC] must be (Conjugate)Transpose'd");
    if( orientationOfB1 == NORMAL )
        throw std::logic_error("B1[MR,* ] must be (Conjugate)Transpose'd");
    if( A1.Grid() != A2.Grid() || A2.Grid() != B1.Grid() ||
        B1.Grid() != B2.Grid() || B2.Grid() != C.Grid() )
        throw std::logic_error
        ("A, B, and C must be distributed over the same grid");
    if( A1.Width() != C.Height() || B1.Height() != C.Width() ||
        A1.Height() != A2.Height() || A1.Width() != A2.Width() ||
        B1.Height() != B2.Width() || B1.Width() != B2.Height() ||
        A1.Height() != B1.Width() )
    {
        std::ostringstream msg;
        msg << "Nonconformal LocalTriangularRank2K: \n"
            << "  A1[* ,MC] ~ " << A1.Height() << " x "
                                << A1.Width()  << "\n"
            << "  A2[* ,MC] ~ " << A2.Height() << " x "
                                << A2.Width()  << "\n"
            << "  B1[MR,* ] ~ " << B1.Height() << " x "
                                << B1.Width()  << "\n"
            << "  B2[* ,MR] ~ " << B2.Height() << " x "
                                << B2.Width()  << "\n"
            << "  C[MC,MR] ~ " << C.Height() << " x " << C.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( A1.RowAlignment() != C.ColAlignment() ||
        B1.ColAlignment() != C.RowAlignment() ||
        A1.RowAlignment() != A2.RowAlignment() ||
        B1.ColAlignment() != B2.RowAlignment() )
    {
        std::ostringstream msg;
        msg << "Misaligned LocalTriangularRank2K: \n"
            << "  A1[* ,MC] ~ " << A1.RowAlignment() << "\n"
            << "  A2[* ,MC] ~ " << A2.RowAlignment() << "\n"
            << "  B1[MR,* ] ~ " << B1.ColAlignment() << "\n"
            << "  B2[* ,MR] ~ " << B2.RowAlignment() << "\n"
            << "  C[MC,MR] ~ " << C.ColAlignment() << " , " <<
                                  C.RowAlignment() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
}

template<typename T>
inline void 
CheckInput
( Orientation orientationOfA1,
  const DistMatrix<T,STAR,MC  >& A1, 
  const DistMatrix<T,MC,  STAR>& A2,
  const DistMatrix<T,STAR,MR  >& B1,
  const DistMatrix<T,STAR,MR  >& B2,
  const DistMatrix<T,MC,  MR  >& C )
{
    if( orientationOfA1 == NORMAL )
        throw std::logic_error("A1[* ,MC] must be (Conjugate)Transpose'd");
    if( A1.Grid() != A2.Grid() || A2.Grid() != B1.Grid() ||
        B1.Grid() != B2.Grid() || B2.Grid() != C.Grid() )
        throw std::logic_error
        ("A, B, and C must be distributed over the same grid");
    if( A1.Width() != C.Height() || B1.Width() != C.Width() ||
        A1.Width() != A2.Height() || A1.Height() != A2.Width() ||
        B1.Height() != B2.Height() || B1.Width() != B2.Width() ||
        A1.Height() != B1.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal LocalTriangularRank2K: \n"
            << "  A1[* ,MC] ~ " << A1.Height() << " x "
                                << A1.Width()  << "\n"
            << "  A2[MC,* ] ~ " << A2.Height() << " x "
                                << A2.Width()  << "\n"
            << "  B1[* ,MR] ~ " << B1.Height() << " x "
                                << B1.Width()  << "\n"
            << "  B2[* ,MR] ~ " << B2.Height() << " x "
                                << B2.Width()  << "\n"
            << "  C[MC,MR] ~ " << C.Height() << " x " << C.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( A1.RowAlignment() != C.ColAlignment() ||
        B1.RowAlignment() != C.RowAlignment() ||
        A1.RowAlignment() != A2.ColAlignment() ||
        B1.RowAlignment() != B2.RowAlignment() )
    {
        std::ostringstream msg;
        msg << "Misaligned LocalTriangularRank2K: \n"
            << "  A1[* ,MC] ~ " << A1.RowAlignment() << "\n"
            << "  A2[MC,* ] ~ " << A2.ColAlignment() << "\n"
            << "  B1[* ,MR] ~ " << B1.RowAlignment() << "\n"
            << "  B2[* ,MR] ~ " << B2.RowAlignment() << "\n"
            << "  C[MC,MR] ~ " << C.ColAlignment() << " , " <<
                                  C.RowAlignment() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
}

template<typename T>
inline void 
CheckInput
( Orientation orientationOfA2,
  const DistMatrix<T,MC,  STAR>& A1,
  const DistMatrix<T,STAR,MC  >& A2, 
  const DistMatrix<T,STAR,MR  >& B1,
  const DistMatrix<T,STAR,MR  >& B2,
  const DistMatrix<T,MC,  MR  >& C )
{
    if( orientationOfA2 == NORMAL )
        throw std::logic_error("A2[* ,MC] must be (Conjugate)Transpose'd");
    if( A1.Grid() != A2.Grid() || A2.Grid() != B1.Grid() ||
        B1.Grid() != B2.Grid() || B2.Grid() != C.Grid() )
        throw std::logic_error
        ("A, B, and C must be distributed over the same grid");
    if( A1.Height() != C.Height() || B1.Width() != C.Width() ||
        A1.Height() != A2.Width() || A1.Width() != A2.Height() ||
        B1.Height() != B2.Height() || B1.Width() != B2.Width() ||
        A1.Width() != B1.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal LocalTriangularRank2K: \n"
            << "  A1[MC,* ] ~ " << A1.Height() << " x "
                                << A1.Width()  << "\n"
            << "  A2[* ,MC] ~ " << A2.Height() << " x "
                                << A2.Width()  << "\n"
            << "  B1[* ,MR] ~ " << B1.Height() << " x "
                                << B1.Width()  << "\n"
            << "  B2[* ,MR] ~ " << B2.Height() << " x "
                                << B2.Width()  << "\n"
            << "  C[MC,MR] ~ " << C.Height() << " x " << C.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( A1.ColAlignment() != C.ColAlignment() ||
        B1.RowAlignment() != C.RowAlignment() ||
        A1.ColAlignment() != A2.RowAlignment() ||
        B1.RowAlignment() != B2.RowAlignment() )
    {
        std::ostringstream msg;
        msg << "Misaligned LocalTriangularRank2K: \n"
            << "  A1[MC,* ] ~ " << A1.ColAlignment() << "\n"
            << "  A2[* ,MC] ~ " << A2.RowAlignment() << "\n"
            << "  B1[* ,MR] ~ " << B1.RowAlignment() << "\n"
            << "  B2[* ,MR] ~ " << B2.RowAlignment() << "\n"
            << "  C[MC,MR] ~ " << C.ColAlignment() << " , " <<
                                  C.RowAlignment() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
}

template<typename T>
inline void 
CheckInput
( Orientation orientationOfA1,
  Orientation orientationOfA2,
  const DistMatrix<T,STAR,MC>& A1, 
  const DistMatrix<T,STAR,MC>& A2,
  const DistMatrix<T,STAR,MR>& B1,
  const DistMatrix<T,STAR,MR>& B2,
  const DistMatrix<T,MC,  MR>& C )
{
    if( orientationOfA1 == NORMAL )
        throw std::logic_error("A1[* ,MC] must be (Conjugate)Transpose'd");
    if( orientationOfA2 == NORMAL )
        throw std::logic_error("A2[* ,MC] must be (Conjugate)Transpose'd");
    if( A1.Grid() != A2.Grid() || A2.Grid() != B1.Grid() ||
        B1.Grid() != B2.Grid() || B2.Grid() != C.Grid() )
        throw std::logic_error
        ("A, B, and C must be distributed over the same grid");
    if( A1.Width() != C.Height() || B1.Width() != C.Width() ||
        A1.Height() != A2.Height() || A1.Width() != A2.Width() ||
        B1.Height() != B2.Height() || B1.Width() != B2.Width() ||
        A1.Height() != B1.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal LocalTriangularRank2K: \n"
            << "  A1[* ,MC] ~ " << A1.Height() << " x "
                                << A1.Width()  << "\n"
            << "  A2[* ,MC] ~ " << A2.Height() << " x "
                                << A2.Width()  << "\n"
            << "  B1[* ,MR] ~ " << B1.Height() << " x "
                                << B1.Width()  << "\n"
            << "  B2[* ,MR] ~ " << B2.Height() << " x "
                                << B2.Width()  << "\n"
            << "  C[MC,MR] ~ " << C.Height() << " x " << C.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( A1.RowAlignment() != C.ColAlignment() ||
        B1.RowAlignment() != C.RowAlignment() ||
        A1.RowAlignment() != A2.RowAlignment() ||
        B1.RowAlignment() != B2.RowAlignment() )
    {
        std::ostringstream msg;
        msg << "Misaligned LocalTriangularRank2K: \n"
            << "  A1[* ,MC] ~ " << A1.RowAlignment() << "\n"
            << "  A2[* ,MC] ~ " << A2.RowAlignment() << "\n"
            << "  B1[* ,MR] ~ " << B1.RowAlignment() << "\n"
            << "  B2[* ,MR] ~ " << B2.RowAlignment() << "\n"
            << "  C[MC,MR] ~ " << C.ColAlignment() << " , " <<
                                  C.RowAlignment() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
}
#endif 

template<typename T>
inline void
LocalTriangularRank2KKernel
( Shape shape,
  Orientation orientationOfB1,
  Orientation orientationOfB2,
  T alpha, const DistMatrix<T,MC,STAR>& A1,
           const DistMatrix<T,MC,STAR>& A2,
           const DistMatrix<T,MR,STAR>& B1,
           const DistMatrix<T,MR,STAR>& B2,
  T beta,        DistMatrix<T,MC,MR  >& C )
{
#ifndef RELEASE
    PushCallStack("LocalTriangularRank2KKernel");
    CheckInput( orientationOfB1, orientationOfB2, A1, A2, B1, B2, C );
#endif
    const Grid& g = C.Grid();

    DistMatrix<T,MC,STAR> A1T(g),  A2T(g),
                          A1B(g),  A2B(g);

    DistMatrix<T,MR,STAR> B1T(g),  B2T(g),
                          B1B(g),  B2B(g);

    DistMatrix<T,MC,MR> CTL(g), CTR(g),
                        CBL(g), CBR(g);

    DistMatrix<T,MC,MR> DTL(g), DBR(g);

    const unsigned half = C.Height()/2;

    basic::Scal( beta, C );

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
    if( shape == LOWER )
    {
        basic::internal::LocalGemm
        ( NORMAL, orientationOfB2, alpha, A1B, B2T, (T)1, CBL );

        basic::internal::LocalGemm
        ( NORMAL, orientationOfB1, alpha, A2B, B1T, (T)1, CBL );
    }
    else
    {
        basic::internal::LocalGemm
        ( NORMAL, orientationOfB2, alpha, A1T, B2B, (T)1, CTR );

        basic::internal::LocalGemm
        ( NORMAL, orientationOfB1, alpha, A2T, B1B, (T)1, CTR );
    }

    basic::internal::LocalGemm
    ( NORMAL, orientationOfB2, alpha, A1T, B2T, (T)0, DTL );

    basic::internal::LocalGemm
    ( NORMAL, orientationOfB1, alpha, A2T, B1T, (T)1, DTL );

    DTL.MakeTrapezoidal( LEFT, shape );
    basic::Axpy( (T)1, DTL, CTL );

    basic::internal::LocalGemm
    ( NORMAL, orientationOfB2, alpha, A1B, B2B, (T)0, DBR );

    basic::internal::LocalGemm
    ( NORMAL, orientationOfB1, alpha, A2B, B1B, (T)1, DBR );

    DBR.MakeTrapezoidal( LEFT, shape );
    basic::Axpy( (T)1, DBR, CBR );
    //------------------------------------------------------------------------//
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
LocalTriangularRank2KKernel
( Shape shape,
  Orientation orientationOfA1,
  Orientation orientationOfB1,
  Orientation orientationOfB2,
  T alpha, const DistMatrix<T,STAR,MC  >& A1,
           const DistMatrix<T,MC,  STAR>& A2,
           const DistMatrix<T,MR,  STAR>& B1,
           const DistMatrix<T,MR,  STAR>& B2,
  T beta,        DistMatrix<T,MC,  MR  >& C )
{
#ifndef RELEASE
    PushCallStack("LocalTriangularRank2KKernel");
    CheckInput
    ( orientationOfA1, orientationOfB1, orientationOfB2, A1, A2, B1, B2, C );
#endif
    const Grid& g = C.Grid();

    DistMatrix<T,STAR,MC> A1L(g), A1R(g);

    DistMatrix<T,MC,STAR> A2T(g),
                          A2B(g);

    DistMatrix<T,MR,STAR> B1T(g), 
                          B1B(g);

    DistMatrix<T,MR,STAR> B2T(g), 
                          B2B(g);

    DistMatrix<T,MC,MR> CTL(g), CTR(g),
                        CBL(g), CBR(g);

    DistMatrix<T,MC,MR> DTL(g), DBR(g);

    const unsigned half = C.Height()/2;

    basic::Scal( beta, C );

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
    if( shape == LOWER )
    {
        basic::internal::LocalGemm
        ( orientationOfA1, orientationOfB2, alpha, A1R, B2T, (T)1, CBL );

        basic::internal::LocalGemm
        ( NORMAL, orientationOfB1, alpha, A2B, B1T, (T)1, CBL );
    }
    else
    {
        basic::internal::LocalGemm
        ( orientationOfA1, orientationOfB2, alpha, A1L, B2B, (T)1, CTR );

        basic::internal::LocalGemm
        ( NORMAL, orientationOfB1, alpha, A2T, B1B, (T)1, CTR );
    }

    basic::internal::LocalGemm
    ( orientationOfA1, orientationOfB2, alpha, A1L, B2T, (T)0, DTL );

    basic::internal::LocalGemm
    ( NORMAL, orientationOfB1, alpha, A2T, B1T, (T)1, DTL );

    DTL.MakeTrapezoidal( LEFT, shape );
    basic::Axpy( (T)1, DTL, CTL );

    basic::internal::LocalGemm
    ( orientationOfA1, orientationOfB2, alpha, A1R, B2B, (T)0, DBR );

    basic::internal::LocalGemm
    ( NORMAL, orientationOfB1, alpha, A2B, B1B, (T)1, DBR );

    DBR.MakeTrapezoidal( LEFT, shape );
    basic::Axpy( (T)1, DBR, CBR );
    //------------------------------------------------------------------------//
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
LocalTriangularRank2KKernel
( Shape shape,
  Orientation orientationOfA2,
  Orientation orientationOfB1,
  Orientation orientationOfB2,
  T alpha, const DistMatrix<T,MC,  STAR>& A1,
           const DistMatrix<T,STAR,MC  >& A2,
           const DistMatrix<T,MR,  STAR>& B1,
           const DistMatrix<T,MR,  STAR>& B2,
  T beta,        DistMatrix<T,MC,  MR  >& C )
{
#ifndef RELEASE
    PushCallStack("LocalTriangularRank2KKernel");
    CheckInput
    ( orientationOfA2, orientationOfB1, orientationOfB2, A1, A2, B1, B2, C );
#endif
    const Grid& g = C.Grid();

    DistMatrix<T,MC,STAR> A1T(g),
                          A1B(g);

    DistMatrix<T,STAR,MC> A2L(g), A2R(g);

    DistMatrix<T,MR,STAR> B1T(g), 
                          B1B(g);

    DistMatrix<T,MR,STAR> B2T(g), 
                          B2B(g);

    DistMatrix<T,MC,MR> CTL(g), CTR(g),
                        CBL(g), CBR(g);

    DistMatrix<T,MC,MR> DTL(g), DBR(g);

    const unsigned half = C.Height()/2;

    basic::Scal( beta, C );

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
    if( shape == LOWER )
    {
        basic::internal::LocalGemm
        ( NORMAL, orientationOfB2, alpha, A1B, B2T, (T)1, CBL );

        basic::internal::LocalGemm
        ( orientationOfA2, orientationOfB1, alpha, A2R, B1T, (T)1, CBL );
    }
    else
    {
        basic::internal::LocalGemm
        ( NORMAL, orientationOfB2, alpha, A1T, B2B, (T)1, CTR );

        basic::internal::LocalGemm
        ( orientationOfA2, orientationOfB1, alpha, A2L, B1B, (T)1, CTR );
    }

    basic::internal::LocalGemm
    ( NORMAL, orientationOfB2, alpha, A1T, B2T, (T)0, DTL );

    basic::internal::LocalGemm
    ( orientationOfA2, orientationOfB1, alpha, A2L, B1T, (T)1, DTL );

    DTL.MakeTrapezoidal( LEFT, shape );
    basic::Axpy( (T)1, DTL, CTL );

    basic::internal::LocalGemm
    ( NORMAL, orientationOfB2, alpha, A1B, B2B, (T)0, DBR );

    basic::internal::LocalGemm
    ( orientationOfA2, orientationOfB1, alpha, A2R, B1B, (T)1, DBR );

    DBR.MakeTrapezoidal( LEFT, shape );
    basic::Axpy( (T)1, DBR, CBR );
    //------------------------------------------------------------------------//
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
LocalTriangularRank2KKernel
( Shape shape,
  Orientation orientationOfB2,
  T alpha, const DistMatrix<T,MC,  STAR>& A1,
           const DistMatrix<T,MC,  STAR>& A2,
           const DistMatrix<T,STAR,MR  >& B1,
           const DistMatrix<T,MR,  STAR>& B2,
  T beta,        DistMatrix<T,MC,  MR  >& C )
{
#ifndef RELEASE
    PushCallStack("LocalTriangularRank2KKernel");
    CheckInput( orientationOfB2, A1, A2, B1, B2, C );
#endif
    const Grid& g = C.Grid();

    DistMatrix<T,MC,STAR> A1T(g),  A2T(g),
                          A1B(g),  A2B(g);

    DistMatrix<T,STAR,MR> B1L(g), B1R(g);

    DistMatrix<T,MR,STAR> B2T(g), 
                          B2B(g);

    DistMatrix<T,MC,MR> CTL(g), CTR(g),
                        CBL(g), CBR(g);

    DistMatrix<T,MC,MR> DTL(g), DBR(g);

    const unsigned half = C.Height()/2;

    basic::Scal( beta, C );

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
    if( shape == LOWER )
    {
        basic::internal::LocalGemm
        ( NORMAL, orientationOfB2, alpha, A1B, B2T, (T)1, CBL );

        basic::internal::LocalGemm
        ( NORMAL, NORMAL, alpha, A2B, B1L, (T)1, CBL );
    }
    else
    {
        basic::internal::LocalGemm
        ( NORMAL, orientationOfB2, alpha, A1T, B2B, (T)1, CTR );

        basic::internal::LocalGemm
        ( NORMAL, NORMAL, alpha, A2T, B1R, (T)1, CTR );
    }

    basic::internal::LocalGemm
    ( NORMAL, orientationOfB2, alpha, A1T, B2T, (T)0, DTL );

    basic::internal::LocalGemm
    ( NORMAL, NORMAL, alpha, A2T, B1L, (T)1, DTL );

    DTL.MakeTrapezoidal( LEFT, shape );
    basic::Axpy( (T)1, DTL, CTL );

    basic::internal::LocalGemm
    ( NORMAL, orientationOfB2, alpha, A1B, B2B, (T)0, DBR );

    basic::internal::LocalGemm
    ( NORMAL, NORMAL, alpha, A2B, B1R, (T)1, DBR );

    DBR.MakeTrapezoidal( LEFT, shape );
    basic::Axpy( (T)1, DBR, CBR );
    //------------------------------------------------------------------------//
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
LocalTriangularRank2KKernel
( Shape shape,
  Orientation orientationOfB1,
  T alpha, const DistMatrix<T,MC,  STAR>& A1,
           const DistMatrix<T,MC,  STAR>& A2,
           const DistMatrix<T,MR,  STAR>& B1,
           const DistMatrix<T,STAR,MR  >& B2,
  T beta,        DistMatrix<T,MC,  MR  >& C )
{
#ifndef RELEASE
    PushCallStack("LocalTriangularRank2KKernel");
    CheckInput( orientationOfB1, A1, A2, B1, B2, C );
#endif
    const Grid& g = C.Grid();

    DistMatrix<T,MC,STAR> A1T(g),  A2T(g),
                          A1B(g),  A2B(g);

    DistMatrix<T,MR,STAR> B1T(g), 
                          B1B(g);

    DistMatrix<T,STAR,MR> B2L(g), B2R(g);

    DistMatrix<T,MC,MR> CTL(g), CTR(g),
                        CBL(g), CBR(g);

    DistMatrix<T,MC,MR> DTL(g), DBR(g);

    const unsigned half = C.Height()/2;

    basic::Scal( beta, C );

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
    if( shape == LOWER )
    {
        basic::internal::LocalGemm
        ( NORMAL, NORMAL, alpha, A1B, B2L, (T)1, CBL );

        basic::internal::LocalGemm
        ( NORMAL, orientationOfB1, alpha, A2B, B1T, (T)1, CBL );
    }
    else
    {
        basic::internal::LocalGemm
        ( NORMAL, NORMAL, alpha, A1T, B2R, (T)1, CTR );

        basic::internal::LocalGemm
        ( NORMAL, orientationOfB1, alpha, A2T, B1B, (T)1, CTR );
    }

    basic::internal::LocalGemm
    ( NORMAL, NORMAL, alpha, A1T, B2L, (T)0, DTL );

    basic::internal::LocalGemm
    ( NORMAL, orientationOfB1, alpha, A2T, B1T, (T)1, DTL );

    DTL.MakeTrapezoidal( LEFT, shape );
    basic::Axpy( (T)1, DTL, CTL );

    basic::internal::LocalGemm
    ( NORMAL, NORMAL, alpha, A1B, B2R, (T)0, DBR );

    basic::internal::LocalGemm
    ( NORMAL, orientationOfB1, alpha, A2B, B1B, (T)1, DBR );

    DBR.MakeTrapezoidal( LEFT, shape );
    basic::Axpy( (T)1, DBR, CBR );
    //------------------------------------------------------------------------//
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
LocalTriangularRank2KKernel
( Shape shape,
  Orientation orientationOfA1,
  Orientation orientationOfA2,
  Orientation orientationOfB1,
  Orientation orientationOfB2,
  T alpha, const DistMatrix<T,STAR,MC  >& A1,
           const DistMatrix<T,STAR,MC  >& A2,
           const DistMatrix<T,MR,  STAR>& B1,
           const DistMatrix<T,MR,  STAR>& B2,
  T beta,        DistMatrix<T,MC,  MR  >& C )
{
#ifndef RELEASE
    PushCallStack("LocalTriangularRank2KKernel");
    CheckInput
    ( orientationOfA1, orientationOfA2, orientationOfB1, orientationOfB2, 
      A1, A2, B1, B2, C );
#endif
    const Grid& g = C.Grid();

    DistMatrix<T,STAR,MC> A1L(g), A1R(g),
                          A2L(g), A2R(g);

    DistMatrix<T,MR,STAR> B1T(g),  B2T(g),
                          B1B(g),  B2B(g);

    DistMatrix<T,MC,MR> CTL(g), CTR(g),
                        CBL(g), CBR(g);

    DistMatrix<T,MC,MR> DTL(g), DBR(g);

    const unsigned half = C.Height()/2;

    basic::Scal( beta, C );

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
    if( shape == LOWER )
    {
        basic::internal::LocalGemm
        ( orientationOfA1, orientationOfB2, alpha, A1R, B2T, (T)1, CBL );

        basic::internal::LocalGemm
        ( orientationOfA2, orientationOfB1, alpha, A2R, B1T, (T)1, CBL );
    }
    else
    {
        basic::internal::LocalGemm
        ( orientationOfA1, orientationOfB2, alpha, A1L, B2B, (T)1, CTR );

        basic::internal::LocalGemm
        ( orientationOfA2, orientationOfB1, alpha, A2L, B1B, (T)1, CTR );
    }

    basic::internal::LocalGemm
    ( orientationOfA1, orientationOfB2, alpha, A1L, B2T, (T)0, DTL );

    basic::internal::LocalGemm
    ( orientationOfA2, orientationOfB1, alpha, A2L, B1T, (T)1, DTL );

    DTL.MakeTrapezoidal( LEFT, shape );
    basic::Axpy( (T)1, DTL, CTL );

    basic::internal::LocalGemm
    ( orientationOfA1, orientationOfB2, alpha, A1R, B2B, (T)0, DBR );

    basic::internal::LocalGemm
    ( orientationOfA2, orientationOfB1, alpha, A2R, B1B, (T)1, DBR );

    DBR.MakeTrapezoidal( LEFT, shape );
    basic::Axpy( (T)1, DBR, CBR );
    //------------------------------------------------------------------------//
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
LocalTriangularRank2KKernel
( Shape shape,
  Orientation orientationOfA1,
  Orientation orientationOfB2,
  T alpha, const DistMatrix<T,STAR,MC  >& A1,
           const DistMatrix<T,MC,  STAR>& A2,
           const DistMatrix<T,STAR,MR  >& B1,
           const DistMatrix<T,MR,  STAR>& B2,
  T beta,        DistMatrix<T,MC,  MR  >& C )
{
#ifndef RELEASE
    PushCallStack("LocalTriangularRank2KKernel");
    CheckInput( orientationOfA1, orientationOfB2, A1, A2, B1, B2, C );
#endif
    const Grid& g = C.Grid();

    DistMatrix<T,STAR,MC> A1L(g), A1R(g);

    DistMatrix<T,MC,STAR> A2T(g),
                          A2B(g);

    DistMatrix<T,STAR,MR> B1L(g), B1R(g);

    DistMatrix<T,MR,STAR> B2T(g),
                          B2B(g);

    DistMatrix<T,MC,MR> CTL(g), CTR(g),
                        CBL(g), CBR(g);

    DistMatrix<T,MC,MR> DTL(g), DBR(g);

    const unsigned half = C.Height()/2;

    basic::Scal( beta, C );

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
    if( shape == LOWER )
    {
        basic::internal::LocalGemm
        ( orientationOfA1, orientationOfB2, alpha, A1R, B2T, (T)1, CBL );

        basic::internal::LocalGemm
        ( NORMAL, NORMAL, alpha, A2B, B1L, (T)1, CBL );
    }
    else
    {
        basic::internal::LocalGemm
        ( orientationOfA1, orientationOfB2, alpha, A1L, B2B, (T)1, CTR );

        basic::internal::LocalGemm
        ( NORMAL, NORMAL, alpha, A2T, B1R, (T)1, CTR );
    }

    basic::internal::LocalGemm
    ( orientationOfA1, orientationOfB2, alpha, A1L, B2T, (T)0, DTL );

    basic::internal::LocalGemm
    ( NORMAL, NORMAL, alpha, A2T, B1L, (T)1, DTL );

    DTL.MakeTrapezoidal( LEFT, shape );
    basic::Axpy( (T)1, DTL, CTL );

    basic::internal::LocalGemm
    ( orientationOfA1, orientationOfB2, alpha, A1R, B2B, (T)0, DBR );

    basic::internal::LocalGemm
    ( NORMAL, NORMAL, alpha, A2B, B1R, (T)1, DBR );

    DBR.MakeTrapezoidal( LEFT, shape );
    basic::Axpy( (T)1, DBR, CBR );
    //------------------------------------------------------------------------//
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
LocalTriangularRank2KKernel
( Shape shape,
  Orientation orientationOfA1,
  Orientation orientationOfB1,
  T alpha, const DistMatrix<T,STAR,MC  >& A1,
           const DistMatrix<T,MC,  STAR>& A2,
           const DistMatrix<T,MR,  STAR>& B1,
           const DistMatrix<T,STAR,MR  >& B2,
  T beta,        DistMatrix<T,MC,  MR  >& C )
{
#ifndef RELEASE
    PushCallStack("LocalTriangularRank2KKernel");
    CheckInput( orientationOfA1, orientationOfB1, A1, A2, B1, B2, C );
#endif
    const Grid& g = C.Grid();

    DistMatrix<T,STAR,MC> A1L(g), A1R(g);

    DistMatrix<T,MC,STAR> A2T(g),
                          A2B(g);

    DistMatrix<T,MR,STAR> B1T(g),
                          B1B(g);

    DistMatrix<T,STAR,MR> B2L(g), B2R(g);

    DistMatrix<T,MC,MR> CTL(g), CTR(g),
                        CBL(g), CBR(g);

    DistMatrix<T,MC,MR> DTL(g), DBR(g);

    const unsigned half = C.Height()/2;

    basic::Scal( beta, C );

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
    if( shape == LOWER )
    {
        basic::internal::LocalGemm
        ( orientationOfA1, NORMAL, alpha, A1R, B2L, (T)1, CBL );

        basic::internal::LocalGemm
        ( NORMAL, orientationOfB1, alpha, A2B, B1T, (T)1, CBL );
    }
    else
    {
        basic::internal::LocalGemm
        ( orientationOfA1, NORMAL, alpha, A1L, B2R, (T)1, CTR );

        basic::internal::LocalGemm
        ( NORMAL, orientationOfB1, alpha, A2T, B1B, (T)1, CTR );
    }

    basic::internal::LocalGemm
    ( orientationOfA1, NORMAL, alpha, A1L, B2L, (T)0, DTL );

    basic::internal::LocalGemm
    ( NORMAL, orientationOfB1, alpha, A2T, B1T, (T)1, DTL );

    DTL.MakeTrapezoidal( LEFT, shape );
    basic::Axpy( (T)1, DTL, CTL );

    basic::internal::LocalGemm
    ( orientationOfA1, NORMAL, alpha, A1R, B2R, (T)0, DBR );

    basic::internal::LocalGemm
    ( NORMAL, orientationOfB1, alpha, A2B, B1B, (T)1, DBR );

    DBR.MakeTrapezoidal( LEFT, shape );
    basic::Axpy( (T)1, DBR, CBR );
    //------------------------------------------------------------------------//
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
LocalTriangularRank2KKernel
( Shape shape,
  Orientation orientationOfA2,
  Orientation orientationOfB2,
  T alpha, const DistMatrix<T,MC,  STAR>& A1,
           const DistMatrix<T,STAR,MC  >& A2,
           const DistMatrix<T,STAR,MR  >& B1,
           const DistMatrix<T,MR,  STAR>& B2,
  T beta,        DistMatrix<T,MC,  MR  >& C )
{
#ifndef RELEASE
    PushCallStack("LocalTriangularRank2KKernel");
    CheckInput( orientationOfA2, orientationOfB2, A1, A2, B1, B2, C );
#endif
    const Grid& g = C.Grid();

    DistMatrix<T,MC,STAR> A1T(g),
                          A1B(g);
    
    DistMatrix<T,STAR,MC> A2L(g), A2R(g);

    DistMatrix<T,STAR,MR> B1L(g), B1R(g);

    DistMatrix<T,MR,STAR> B2T(g),
                          B2B(g);

    DistMatrix<T,MC,MR> CTL(g), CTR(g),
                        CBL(g), CBR(g);

    DistMatrix<T,MC,MR> DTL(g), DBR(g);

    const unsigned half = C.Height()/2;

    basic::Scal( beta, C );

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
    if( shape == LOWER )
    {
        basic::internal::LocalGemm
        ( NORMAL, orientationOfB2, alpha, A1B, B2T, (T)1, CBL );

        basic::internal::LocalGemm
        ( orientationOfA2, NORMAL, alpha, A2R, B1L, (T)1, CBL );
    }
    else
    {
        basic::internal::LocalGemm
        ( NORMAL, orientationOfB2, alpha, A1T, B2B, (T)1, CTR );

        basic::internal::LocalGemm
        ( orientationOfA2, NORMAL, alpha, A2L, B1R, (T)1, CTR );
    }

    basic::internal::LocalGemm
    ( NORMAL, orientationOfB2, alpha, A1T, B2T, (T)0, DTL );

    basic::internal::LocalGemm
    ( orientationOfA2, NORMAL, alpha, A2L, B1L, (T)1, DTL );

    DTL.MakeTrapezoidal( LEFT, shape );
    basic::Axpy( (T)1, DTL, CTL );

    basic::internal::LocalGemm
    ( NORMAL, orientationOfB2, alpha, A1B, B2B, (T)0, DBR );

    basic::internal::LocalGemm
    ( orientationOfA2, NORMAL, alpha, A2R, B1R, (T)1, DBR );

    DBR.MakeTrapezoidal( LEFT, shape );
    basic::Axpy( (T)1, DBR, CBR );
    //------------------------------------------------------------------------//
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
LocalTriangularRank2KKernel
( Shape shape,
  Orientation orientationOfA2,
  Orientation orientationOfB1,
  T alpha, const DistMatrix<T,MC,  STAR>& A1,
           const DistMatrix<T,STAR,MC  >& A2,
           const DistMatrix<T,MR,  STAR>& B1,
           const DistMatrix<T,STAR,MR  >& B2,
  T beta,        DistMatrix<T,MC,  MR  >& C )
{
#ifndef RELEASE
    PushCallStack("LocalTriangularRank2KKernel");
    CheckInput( orientationOfA2, orientationOfB1, A1, A2, B1, B2, C );
#endif
    const Grid& g = C.Grid();

    DistMatrix<T,MC,STAR> A1T(g),
                          A1B(g);
    
    DistMatrix<T,STAR,MC> A2L(g), A2R(g);

    DistMatrix<T,MR,STAR> B1T(g),
                          B1B(g);

    DistMatrix<T,STAR,MR> B2L(g), B2R(g);

    DistMatrix<T,MC,MR> CTL(g), CTR(g),
                        CBL(g), CBR(g);

    DistMatrix<T,MC,MR> DTL(g), DBR(g);

    const unsigned half = C.Height()/2;

    basic::Scal( beta, C );

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
    if( shape == LOWER )
    {
        basic::internal::LocalGemm
        ( NORMAL, NORMAL, alpha, A1B, B2L, (T)1, CBL );

        basic::internal::LocalGemm
        ( orientationOfA2, orientationOfB1, alpha, A2R, B1T, (T)1, CBL );
    }
    else
    {
        basic::internal::LocalGemm
        ( NORMAL, NORMAL, alpha, A1T, B2R, (T)1, CTR );

        basic::internal::LocalGemm
        ( orientationOfA2, orientationOfB1, alpha, A2L, B1B, (T)1, CTR );
    }

    basic::internal::LocalGemm
    ( NORMAL, NORMAL, alpha, A1T, B2L, (T)0, DTL );

    basic::internal::LocalGemm
    ( orientationOfA2, orientationOfB1, alpha, A2L, B1T, (T)1, DTL );

    DTL.MakeTrapezoidal( LEFT, shape );
    basic::Axpy( (T)1, DTL, CTL );

    basic::internal::LocalGemm
    ( NORMAL, NORMAL, alpha, A1B, B2R, (T)0, DBR );

    basic::internal::LocalGemm
    ( orientationOfA2, orientationOfB1, alpha, A2R, B1B, (T)1, DBR );

    DBR.MakeTrapezoidal( LEFT, shape );
    basic::Axpy( (T)1, DBR, CBR );
    //------------------------------------------------------------------------//
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
LocalTriangularRank2KKernel
( Shape shape,
  T alpha, const DistMatrix<T,MC,  STAR>& A1,
           const DistMatrix<T,MC,  STAR>& A2,
           const DistMatrix<T,STAR,MR  >& B1,
           const DistMatrix<T,STAR,MR  >& B2,
  T beta,        DistMatrix<T,MC,  MR  >& C )
{
#ifndef RELEASE
    PushCallStack("LocalTriangularRank2KKernel");
    CheckInput( A1, A2, B1, B2, C );
#endif
    const Grid& g = C.Grid();

    DistMatrix<T,MC,STAR> A1T(g),  A2T(g),
                          A1B(g),  A2B(g);

    DistMatrix<T,STAR,MR> B1L(g), B1R(g),
                          B2L(g), B2R(g);

    DistMatrix<T,MC,MR> CTL(g), CTR(g),
                        CBL(g), CBR(g);

    DistMatrix<T,MC,MR> DTL(g), DBR(g);

    const unsigned half = C.Height()/2;

    basic::Scal( beta, C );

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
    if( shape == LOWER )
    {
        basic::internal::LocalGemm
        ( NORMAL, NORMAL, alpha, A1B, B2L, (T)1, CBL );

        basic::internal::LocalGemm
        ( NORMAL, NORMAL, alpha, A2B, B1L, (T)1, CBL );
    }
    else
    {
        basic::internal::LocalGemm
        ( NORMAL, NORMAL, alpha, A1T, B2R, (T)1, CTR );

        basic::internal::LocalGemm
        ( NORMAL, NORMAL, alpha, A2T, B1R, (T)1, CTR );
    }

    basic::internal::LocalGemm
    ( NORMAL, NORMAL, alpha, A1T, B2L, (T)0, DTL );

    basic::internal::LocalGemm
    ( NORMAL, NORMAL, alpha, A2T, B1L, (T)1, DTL );

    DTL.MakeTrapezoidal( LEFT, shape );
    basic::Axpy( (T)1, DTL, CTL );

    basic::internal::LocalGemm
    ( NORMAL, NORMAL, alpha, A1B, B2R, (T)0, DBR );

    basic::internal::LocalGemm
    ( NORMAL, NORMAL, alpha, A2B, B1R, (T)1, DBR );

    DBR.MakeTrapezoidal( LEFT, shape );
    basic::Axpy( (T)1, DBR, CBR );
    //------------------------------------------------------------------------//
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
LocalTriangularRank2KKernel
( Shape shape,
  Orientation orientationOfA1,
  Orientation orientationOfA2,
  Orientation orientationOfB2,
  T alpha, const DistMatrix<T,STAR,MC  >& A1,
           const DistMatrix<T,STAR,MC  >& A2,
           const DistMatrix<T,STAR,MR  >& B1,
           const DistMatrix<T,MR,  STAR>& B2,
  T beta,        DistMatrix<T,MC,  MR>& C )
{
#ifndef RELEASE
    PushCallStack("LocalTriangularRank2KKernel");
    CheckInput
    ( orientationOfA1, orientationOfA2, orientationOfB2, A1, A2, B1, B2, C );
#endif
    const Grid& g = C.Grid();

    DistMatrix<T,STAR,MC> A1L(g), A1R(g),
                          A2L(g), A2R(g);

    DistMatrix<T,STAR,MR> B1L(g), B1R(g);

    DistMatrix<T,MR,STAR> B2T(g),
                          B2B(g);

    DistMatrix<T,MC,MR> CTL(g), CTR(g),
                        CBL(g), CBR(g);

    DistMatrix<T,MC,MR> DTL(g), DBR(g);

    const unsigned half = C.Height()/2;

    basic::Scal( beta, C );

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
    if( shape == LOWER )
    {
        basic::internal::LocalGemm
        ( orientationOfA1, orientationOfB2, alpha, A1R, B2T, (T)1, CBL );
        
        basic::internal::LocalGemm
        ( orientationOfA2, NORMAL, alpha, A2R, B1L, (T)1, CBL );
    }
    else
    {
        basic::internal::LocalGemm
        ( orientationOfA1, orientationOfB2, alpha, A1L, B2B, (T)1, CTR );

        basic::internal::LocalGemm
        ( orientationOfA2, NORMAL, alpha, A2L, B1R, (T)1, CTR );
    }

    basic::internal::LocalGemm
    ( orientationOfA1, orientationOfB2, alpha, A1L, B2T, (T)0, DTL );

    basic::internal::LocalGemm
    ( orientationOfA2, NORMAL, alpha, A2L, B1L, (T)1, DTL );

    DTL.MakeTrapezoidal( LEFT, shape );
    basic::Axpy( (T)1, DTL, CTL );

    basic::internal::LocalGemm
    ( orientationOfA1, orientationOfB2, alpha, A1R, B2B, (T)0, DBR );

    basic::internal::LocalGemm
    ( orientationOfA2, NORMAL, alpha, A2R, B1R, (T)1, DBR );

    DBR.MakeTrapezoidal( LEFT, shape );
    basic::Axpy( (T)1, DBR, CBR );
    //------------------------------------------------------------------------//
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
LocalTriangularRank2KKernel
( Shape shape,
  Orientation orientationOfA1,
  Orientation orientationOfA2,
  Orientation orientationOfB1,
  T alpha, const DistMatrix<T,STAR,MC  >& A1,
           const DistMatrix<T,STAR,MC  >& A2,
           const DistMatrix<T,MR,  STAR>& B1,
           const DistMatrix<T,STAR,MR  >& B2,
  T beta,        DistMatrix<T,MC,  MR  >& C )
{
#ifndef RELEASE
    PushCallStack("LocalTriangularRank2KKernel");
    CheckInput
    ( orientationOfA1, orientationOfA2, orientationOfB1, A1, A2, B1, B2, C );
#endif
    const Grid& g = C.Grid();

    DistMatrix<T,STAR,MC> A1L(g), A1R(g),
                          A2L(g), A2R(g);

    DistMatrix<T,MR,STAR> B1T(g),
                          B1B(g);

    DistMatrix<T,STAR,MR> B2L(g), B2R(g);

    DistMatrix<T,MC,MR> CTL(g), CTR(g),
                        CBL(g), CBR(g);

    DistMatrix<T,MC,MR> DTL(g), DBR(g);

    const unsigned half = C.Height()/2;

    basic::Scal( beta, C );

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
    if( shape == LOWER )
    {
        basic::internal::LocalGemm
        ( orientationOfA1, NORMAL, alpha, A1R, B2L, (T)1, CBL );
        
        basic::internal::LocalGemm
        ( orientationOfA2, orientationOfB1, alpha, A2R, B1T, (T)1, CBL );
    }
    else
    {
        basic::internal::LocalGemm
        ( orientationOfA1, NORMAL, alpha, A1L, B2R, (T)1, CTR );

        basic::internal::LocalGemm
        ( orientationOfA2, orientationOfB1, alpha, A2L, B1B, (T)1, CTR );
    }

    basic::internal::LocalGemm
    ( orientationOfA1, NORMAL, alpha, A1L, B2L, (T)0, DTL );

    basic::internal::LocalGemm
    ( orientationOfA2, orientationOfB1, alpha, A2L, B1T, (T)1, DTL );

    DTL.MakeTrapezoidal( LEFT, shape );
    basic::Axpy( (T)1, DTL, CTL );

    basic::internal::LocalGemm
    ( orientationOfA1, NORMAL, alpha, A1R, B2R, (T)0, DBR );

    basic::internal::LocalGemm
    ( orientationOfA2, orientationOfB1, alpha, A2R, B1B, (T)1, DBR );

    DBR.MakeTrapezoidal( LEFT, shape );
    basic::Axpy( (T)1, DBR, CBR );
    //------------------------------------------------------------------------//
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
LocalTriangularRank2KKernel
( Shape shape,
  Orientation orientationOfA1,
  T alpha, const DistMatrix<T,STAR,MC  >& A1,
           const DistMatrix<T,MC,  STAR>& A2,
           const DistMatrix<T,STAR,MR  >& B1,
           const DistMatrix<T,STAR,MR  >& B2,
  T beta,        DistMatrix<T,MC,  MR  >& C )
{
#ifndef RELEASE
    PushCallStack("LocalTriangularRank2KKernel");
    CheckInput( orientationOfA1, A1, A2, B1, B2, C );
#endif
    const Grid& g = C.Grid();

    DistMatrix<T,STAR,MC> A1L(g), A1R(g);

    DistMatrix<T,MC,STAR> A2T(g),
                          A2B(g);

    DistMatrix<T,STAR,MR> B1L(g), B1R(g),
                          B2L(g), B2R(g);

    DistMatrix<T,MC,MR> CTL(g), CTR(g),
                        CBL(g), CBR(g);

    DistMatrix<T,MC,MR> DTL(g), DBR(g);

    const unsigned half = C.Height()/2;

    basic::Scal( beta, C );

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
    if( shape == LOWER )
    {
        basic::internal::LocalGemm
        ( orientationOfA1, NORMAL, alpha, A1R, B2L, (T)1, CBL );
        
        basic::internal::LocalGemm
        ( NORMAL, NORMAL, alpha, A2B, B1L, (T)1, CBL );
    }
    else
    {
        basic::internal::LocalGemm
        ( orientationOfA1, NORMAL, alpha, A1L, B2R, (T)1, CTR );

        basic::internal::LocalGemm
        ( NORMAL, NORMAL, alpha, A2T, B1R, (T)1, CTR );
    }

    basic::internal::LocalGemm
    ( orientationOfA1, NORMAL, alpha, A1L, B2L, (T)0, DTL );

    basic::internal::LocalGemm
    ( NORMAL, NORMAL, alpha, A2T, B1L, (T)1, DTL );

    DTL.MakeTrapezoidal( LEFT, shape );
    basic::Axpy( (T)1, DTL, CTL );

    basic::internal::LocalGemm
    ( orientationOfA1, NORMAL, alpha, A1R, B2R, (T)0, DBR );

    basic::internal::LocalGemm
    ( NORMAL, NORMAL, alpha, A2B, B1R, (T)1, DBR );

    DBR.MakeTrapezoidal( LEFT, shape );
    basic::Axpy( (T)1, DBR, CBR );
    //------------------------------------------------------------------------//
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
LocalTriangularRank2KKernel
( Shape shape,
  Orientation orientationOfA2,
  T alpha, const DistMatrix<T,MC,  STAR>& A1,
           const DistMatrix<T,STAR,MC  >& A2,
           const DistMatrix<T,STAR,MR  >& B1,
           const DistMatrix<T,STAR,MR  >& B2,
  T beta,        DistMatrix<T,MC,  MR  >& C )
{
#ifndef RELEASE
    PushCallStack("LocalTriangularRank2KKernel");
    CheckInput( orientationOfA2, A1, A2, B1, B2, C );
#endif
    const Grid& g = C.Grid();

    DistMatrix<T,MC,STAR> A1T(g),
                          A1B(g);

    DistMatrix<T,STAR,MC> A2L(g), A2R(g);

    DistMatrix<T,STAR,MR> B1L(g), B1R(g),
                          B2L(g), B2R(g);

    DistMatrix<T,MC,MR> CTL(g), CTR(g),
                        CBL(g), CBR(g);

    DistMatrix<T,MC,MR> DTL(g), DBR(g);

    const unsigned half = C.Height()/2;

    basic::Scal( beta, C );

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
    if( shape == LOWER )
    {
        basic::internal::LocalGemm
        ( NORMAL, NORMAL, alpha, A1B, B2L, (T)1, CBL );
        
        basic::internal::LocalGemm
        ( orientationOfA2, NORMAL, alpha, A2R, B1L, (T)1, CBL );
    }
    else
    {
        basic::internal::LocalGemm
        ( NORMAL, NORMAL, alpha, A1T, B2R, (T)1, CTR );

        basic::internal::LocalGemm
        ( orientationOfA2, NORMAL, alpha, A2L, B1R, (T)1, CTR );
    }

    basic::internal::LocalGemm
    ( NORMAL, NORMAL, alpha, A1T, B2L, (T)0, DTL );

    basic::internal::LocalGemm
    ( orientationOfA2, NORMAL, alpha, A2L, B1L, (T)1, DTL );

    DTL.MakeTrapezoidal( LEFT, shape );
    basic::Axpy( (T)1, DTL, CTL );

    basic::internal::LocalGemm
    ( NORMAL, NORMAL, alpha, A1B, B2R, (T)0, DBR );

    basic::internal::LocalGemm
    ( orientationOfA2, NORMAL, alpha, A2R, B1R, (T)1, DBR );

    DBR.MakeTrapezoidal( LEFT, shape );
    basic::Axpy( (T)1, DBR, CBR );
    //------------------------------------------------------------------------//
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
LocalTriangularRank2KKernel
( Shape shape,
  Orientation orientationOfA1,
  Orientation orientationOfA2,
  T alpha, const DistMatrix<T,STAR,MC>& A1,
           const DistMatrix<T,STAR,MC>& A2,
           const DistMatrix<T,STAR,MR>& B1,
           const DistMatrix<T,STAR,MR>& B2,
  T beta,        DistMatrix<T,MC,  MR>& C )
{
#ifndef RELEASE
    PushCallStack("LocalTriangularRank2KKernel");
    CheckInput( orientationOfA1, orientationOfA2, A1, A2, B1, B2, C );
#endif
    const Grid& g = C.Grid();

    DistMatrix<T,STAR,MC> A1L(g), A1R(g),
                          A2L(g), A2R(g);

    DistMatrix<T,STAR,MR> B1L(g), B1R(g),
                          B2L(g), B2R(g);

    DistMatrix<T,MC,MR> CTL(g), CTR(g),
                        CBL(g), CBR(g);

    DistMatrix<T,MC,MR> DTL(g), DBR(g);

    const unsigned half = C.Height()/2;

    basic::Scal( beta, C );

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
    if( shape == LOWER )
    {
        basic::internal::LocalGemm
        ( orientationOfA1, NORMAL, alpha, A1R, B2L, (T)1, CBL );
        
        basic::internal::LocalGemm
        ( orientationOfA2, NORMAL, alpha, A2R, B1L, (T)1, CBL );
    }
    else
    {
        basic::internal::LocalGemm
        ( orientationOfA1, NORMAL, alpha, A1L, B2R, (T)1, CTR );

        basic::internal::LocalGemm
        ( orientationOfA2, NORMAL, alpha, A2L, B1R, (T)1, CTR );
    }

    basic::internal::LocalGemm
    ( orientationOfA1, NORMAL, alpha, A1L, B2L, (T)0, DTL );

    basic::internal::LocalGemm
    ( orientationOfA2, NORMAL, alpha, A2L, B1L, (T)1, DTL );

    DTL.MakeTrapezoidal( LEFT, shape );
    basic::Axpy( (T)1, DTL, CTL );

    basic::internal::LocalGemm
    ( orientationOfA1, NORMAL, alpha, A1R, B2R, (T)0, DBR );

    basic::internal::LocalGemm
    ( orientationOfA2, NORMAL, alpha, A2R, B1R, (T)1, DBR );

    DBR.MakeTrapezoidal( LEFT, shape );
    basic::Axpy( (T)1, DBR, CBR );
    //------------------------------------------------------------------------//
#ifndef RELEASE
    PopCallStack();
#endif
}

} // anonymous namespace

template<typename T>
inline void
elemental::basic::internal::LocalTriangularRank2K
( Shape shape,
  Orientation orientationOfB1,
  Orientation orientationOfB2,
  T alpha, const DistMatrix<T,MC,STAR>& A1,
           const DistMatrix<T,MC,STAR>& A2,
           const DistMatrix<T,MR,STAR>& B1,
           const DistMatrix<T,MR,STAR>& B2,
  T beta,        DistMatrix<T,MC,MR  >& C  )
{
#ifndef RELEASE
    PushCallStack("basic::internal::LocalTriangularRank2K");
    CheckInput( orientationOfB1, orientationOfB2, A1, A2, B1, B2, C );
#endif
    const Grid& g = C.Grid();

    if( C.Height() < g.Width()*LocalTriangularRank2KBlocksize<T>() )
    {
        LocalTriangularRank2KKernel
        ( shape, orientationOfB1, orientationOfB2, 
          alpha, A1, A2, B1, B2, beta, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.

        DistMatrix<T,MC,STAR> A1T(g),  A2T(g),
                              A1B(g),  A2B(g);

        DistMatrix<T,MR,STAR> B1T(g),  B2T(g),
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

        if( shape == LOWER )
        { 
            basic::internal::LocalGemm
            ( NORMAL, orientationOfB2, alpha, A1B, B2T, beta, CBL );

            basic::internal::LocalGemm
            ( NORMAL, orientationOfB1, alpha, A2B, B1T, (T)1, CBL );
        }
        else
        {
            basic::internal::LocalGemm
            ( NORMAL, orientationOfB2, alpha, A1T, B2B, beta, CTR );

            basic::internal::LocalGemm
            ( NORMAL, orientationOfB1, alpha, A2T, B1B, (T)1, CTR );
        }

        // Recurse
        basic::internal::LocalTriangularRank2K
        ( shape, orientationOfB1, orientationOfB2,
          alpha, A1T, A2T, B1T, B2T, beta, CTL );

        basic::internal::LocalTriangularRank2K
        ( shape, orientationOfB1, orientationOfB2,
          alpha, A1B, A2B, B1B, B2B, beta, CBR );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::basic::internal::LocalTriangularRank2K
( Shape shape,
  Orientation orientationOfA1,
  Orientation orientationOfB1,
  Orientation orientationOfB2,
  T alpha, const DistMatrix<T,STAR,MC  >& A1,
           const DistMatrix<T,MC,  STAR>& A2,
           const DistMatrix<T,MR,  STAR>& B1,
           const DistMatrix<T,MR,  STAR>& B2,
  T beta,        DistMatrix<T,MC,  MR  >& C  )
{
#ifndef RELEASE
    PushCallStack("basic::internal::LocalTriangularRank2K");
    CheckInput
    ( orientationOfA1, orientationOfB1, orientationOfB2, A1, A2, B1, B2, C );
#endif
    const Grid& g = C.Grid();

    if( C.Height() < g.Width()*LocalTriangularRank2KBlocksize<T>() )
    {
        LocalTriangularRank2KKernel
        ( shape, orientationOfA1, orientationOfB1, orientationOfB2, 
          alpha, A1, A2, B1, B2, beta, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.

        DistMatrix<T,STAR,MC> A1L(g), A1R(g);

        DistMatrix<T,MC,STAR> A2T(g),
                              A2B(g); 

        DistMatrix<T,MR,STAR> B1T(g),  B2T(g),
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

        if( shape == LOWER )
        { 
            basic::internal::LocalGemm
            ( orientationOfA1, orientationOfB2, alpha, A1R, B2T, beta, CBL );

            basic::internal::LocalGemm
            ( NORMAL, orientationOfB1, alpha, A2B, B1T, (T)1, CBL );
        }
        else
        {
            basic::internal::LocalGemm
            ( orientationOfA1, orientationOfB2, alpha, A1L, B2B, beta, CTR );

            basic::internal::LocalGemm
            ( NORMAL, orientationOfB1, alpha, A2T, B1B, (T)1, CTR );
        }

        // Recurse
        basic::internal::LocalTriangularRank2K
        ( shape, orientationOfA1, orientationOfB1, orientationOfB2,
          alpha, A1L, A2T, B1T, B2T, beta, CTL );

        basic::internal::LocalTriangularRank2K
        ( shape, orientationOfA1, orientationOfB1, orientationOfB2,
          alpha, A1R, A2B, B1B, B2B, beta, CBR );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::basic::internal::LocalTriangularRank2K
( Shape shape,
  Orientation orientationOfA2,
  Orientation orientationOfB1,
  Orientation orientationOfB2,
  T alpha, const DistMatrix<T,MC,  STAR>& A1,
           const DistMatrix<T,STAR,MC  >& A2,
           const DistMatrix<T,MR,  STAR>& B1,
           const DistMatrix<T,MR,  STAR>& B2,
  T beta,        DistMatrix<T,MC,  MR  >& C  )
{
#ifndef RELEASE
    PushCallStack("basic::internal::LocalTriangularRank2K");
    CheckInput
    ( orientationOfA2, orientationOfB1, orientationOfB2, A1, A2, B1, B2, C );
#endif
    const Grid& g = C.Grid();

    if( C.Height() < g.Width()*LocalTriangularRank2KBlocksize<T>() )
    {
        LocalTriangularRank2KKernel
        ( shape, orientationOfA2, orientationOfB1, orientationOfB2, 
          alpha, A1, A2, B1, B2, beta, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.

        DistMatrix<T,MC,STAR> A1T(g),
                              A1B(g); 

        DistMatrix<T,STAR,MC> A2L(g), A2R(g);

        DistMatrix<T,MR,STAR> B1T(g),  B2T(g),
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

        if( shape == LOWER )
        { 
            basic::internal::LocalGemm
            ( NORMAL, orientationOfB2, alpha, A1B, B2T, beta, CBL );

            basic::internal::LocalGemm
            ( orientationOfA2, orientationOfB1, alpha, A2R, B1T, (T)1, CBL );
        }
        else
        {
            basic::internal::LocalGemm
            ( NORMAL, orientationOfB2, alpha, A1T, B2B, beta, CTR );

            basic::internal::LocalGemm
            ( orientationOfA2, orientationOfB1, alpha, A2L, B1B, (T)1, CTR );
        }

        // Recurse
        basic::internal::LocalTriangularRank2K
        ( shape, orientationOfA2, orientationOfB1, orientationOfB2,
          alpha, A1T, A2L, B1T, B2T, beta, CTL );

        basic::internal::LocalTriangularRank2K
        ( shape, orientationOfA2, orientationOfB1, orientationOfB2,
          alpha, A1B, A2R, B1B, B2B, beta, CBR );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::basic::internal::LocalTriangularRank2K
( Shape shape,
  Orientation orientationOfB2,
  T alpha, const DistMatrix<T,MC,  STAR>& A1,
           const DistMatrix<T,MC,  STAR>& A2,
           const DistMatrix<T,STAR,MR  >& B1,
           const DistMatrix<T,MR,  STAR>& B2,
  T beta,        DistMatrix<T,MC,  MR  >& C  )
{
#ifndef RELEASE
    PushCallStack("basic::internal::LocalTriangularRank2K");
    CheckInput( orientationOfB2, A1, A2, B1, B2, C );
#endif
    const Grid& g = C.Grid();

    if( C.Height() < g.Width()*LocalTriangularRank2KBlocksize<T>() )
    {
        LocalTriangularRank2KKernel
        ( shape, orientationOfB2, alpha, A1, A2, B1, B2, beta, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.

        DistMatrix<T,MC,STAR> A1T(g),  A2T(g),
                              A1B(g),  A2B(g);

        DistMatrix<T,STAR,MR> B1L(g), B1R(g);
        
        DistMatrix<T,MR,STAR> B2T(g), 
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

        if( shape == LOWER )
        { 
            basic::internal::LocalGemm
            ( NORMAL, orientationOfB2, alpha, A1B, B2T, beta, CBL );

            basic::internal::LocalGemm
            ( NORMAL, NORMAL, alpha, A2B, B1L, (T)1, CBL );
        }
        else
        {
            basic::internal::LocalGemm
            ( NORMAL, orientationOfB2, alpha, A1T, B2B, beta, CTR );

            basic::internal::LocalGemm
            ( NORMAL, NORMAL, alpha, A2T, B1R, (T)1, CTR );
        }

        // Recurse
        basic::internal::LocalTriangularRank2K
        ( shape, orientationOfB2, alpha, A1T, A2T, B1L, B2T, beta, CTL );

        basic::internal::LocalTriangularRank2K
        ( shape, orientationOfB2, alpha, A1B, A2B, B1R, B2B, beta, CBR );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::basic::internal::LocalTriangularRank2K
( Shape shape,
  Orientation orientationOfB1,
  T alpha, const DistMatrix<T,MC,  STAR>& A1,
           const DistMatrix<T,MC,  STAR>& A2,
           const DistMatrix<T,MR,  STAR>& B1,
           const DistMatrix<T,STAR,MR  >& B2,
  T beta,        DistMatrix<T,MC,  MR  >& C  )
{
#ifndef RELEASE
    PushCallStack("basic::internal::LocalTriangularRank2K");
    CheckInput( orientationOfB1, A1, A2, B1, B2, C );
#endif
    const Grid& g = C.Grid();

    if( C.Height() < g.Width()*LocalTriangularRank2KBlocksize<T>() )
    {
        LocalTriangularRank2KKernel
        ( shape, orientationOfB1, alpha, A1, A2, B1, B2, beta, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.

        DistMatrix<T,MC,STAR> A1T(g),  A2T(g),
                              A1B(g),  A2B(g);

        DistMatrix<T,MR,STAR> B1T(g), 
                              B1B(g);

        DistMatrix<T,STAR,MR> B2L(g), B2R(g);

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

        if( shape == LOWER )
        { 
            basic::internal::LocalGemm
            ( NORMAL, NORMAL, alpha, A1B, B2L, beta, CBL );

            basic::internal::LocalGemm
            ( NORMAL, orientationOfB1, alpha, A2B, B1T, (T)1, CBL );
        }
        else
        {
            basic::internal::LocalGemm
            ( NORMAL, NORMAL, alpha, A1T, B2R, beta, CTR );

            basic::internal::LocalGemm
            ( NORMAL, orientationOfB1, alpha, A2T, B1B, (T)1, CTR );
        }

        // Recurse
        basic::internal::LocalTriangularRank2K
        ( shape, orientationOfB1, alpha, A1T, A2T, B1T, B2L, beta, CTL );

        basic::internal::LocalTriangularRank2K
        ( shape, orientationOfB1, alpha, A1B, A2B, B1B, B2R, beta, CBR );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::basic::internal::LocalTriangularRank2K
( Shape shape,
  Orientation orientationOfA1,
  Orientation orientationOfA2,
  Orientation orientationOfB1,
  Orientation orientationOfB2,
  T alpha, const DistMatrix<T,STAR,MC  >& A1,
           const DistMatrix<T,STAR,MC  >& A2,
           const DistMatrix<T,MR,  STAR>& B1,
           const DistMatrix<T,MR,  STAR>& B2,
  T beta,        DistMatrix<T,MC,  MR  >& C  )
{
#ifndef RELEASE
    PushCallStack("basic::internal::LocalTriangularRank2K");
    CheckInput
    ( orientationOfA1, orientationOfA2, orientationOfB1, orientationOfB2, 
      A1, A2, B1, B2, C );
#endif
    const Grid& g = C.Grid();

    if( C.Height() < g.Width()*LocalTriangularRank2KBlocksize<T>() )
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

        DistMatrix<T,STAR,MC> A1L(g), A1R(g),
                              A2L(g), A2R(g);

        DistMatrix<T,MR,STAR> B1T(g),  B2T(g),
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

        if( shape == LOWER )
        { 
            basic::internal::LocalGemm
            ( orientationOfA1, orientationOfB2, alpha, A1R, B2T, beta, CBL );

            basic::internal::LocalGemm
            ( orientationOfA2, orientationOfB1, alpha, A2R, B1T, (T)1, CBL );
        }
        else
        {
            basic::internal::LocalGemm
            ( orientationOfA1, orientationOfB2, alpha, A1L, B2B, beta, CTR );

            basic::internal::LocalGemm
            ( orientationOfA2, orientationOfB1, alpha, A2L, B1B, (T)1, CTR );
        }

        // Recurse
        basic::internal::LocalTriangularRank2K
        ( shape, 
          orientationOfA1, orientationOfA2, 
          orientationOfB1, orientationOfB2, 
          alpha, A1L, A2L, B1T, B2T, beta, CTL );

        basic::internal::LocalTriangularRank2K
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
inline void
elemental::basic::internal::LocalTriangularRank2K
( Shape shape,
  Orientation orientationOfA1,
  Orientation orientationOfB2,
  T alpha, const DistMatrix<T,STAR,MC  >& A1,
           const DistMatrix<T,MC,  STAR>& A2,
           const DistMatrix<T,STAR,MR  >& B1,
           const DistMatrix<T,MR,  STAR>& B2,
  T beta,        DistMatrix<T,MC,  MR  >& C  )
{
#ifndef RELEASE
    PushCallStack("basic::internal::LocalTriangularRank2K");
    CheckInput
    ( orientationOfA1, orientationOfB2, A1, A2, B1, B2, C );
#endif
    const Grid& g = C.Grid();

    if( C.Height() < g.Width()*LocalTriangularRank2KBlocksize<T>() )
    {
        LocalTriangularRank2KKernel
        ( shape, orientationOfA1, orientationOfB2, 
          alpha, A1, A2, B1, B2, beta, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.

        DistMatrix<T,STAR,MC> A1L(g), A1R(g);

        DistMatrix<T,MC,STAR> A2T(g),
                              A2B(g); 

        DistMatrix<T,STAR,MR> B1L(g), B1R(g);

        DistMatrix<T,MR,STAR> B2T(g),
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

        if( shape == LOWER )
        { 
            basic::internal::LocalGemm
            ( orientationOfA1, orientationOfB2, alpha, A1R, B2T, beta, CBL );

            basic::internal::LocalGemm
            ( NORMAL, NORMAL, alpha, A2B, B1L, (T)1, CBL );
        }
        else
        {
            basic::internal::LocalGemm
            ( orientationOfA1, orientationOfB2, alpha, A1L, B2B, beta, CTR );

            basic::internal::LocalGemm
            ( NORMAL, NORMAL, alpha, A2T, B1R, (T)1, CTR );
        }

        // Recurse
        basic::internal::LocalTriangularRank2K
        ( shape, orientationOfA1, orientationOfB2,
          alpha, A1L, A2T, B1L, B2T, beta, CTL );

        basic::internal::LocalTriangularRank2K
        ( shape, orientationOfA1, orientationOfB2,
          alpha, A1R, A2B, B1R, B2B, beta, CBR );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::basic::internal::LocalTriangularRank2K
( Shape shape,
  Orientation orientationOfA1,
  Orientation orientationOfB1,
  T alpha, const DistMatrix<T,STAR,MC  >& A1,
           const DistMatrix<T,MC,  STAR>& A2,
           const DistMatrix<T,MR,  STAR>& B1,
           const DistMatrix<T,STAR,MR  >& B2,
  T beta,        DistMatrix<T,MC,  MR  >& C  )
{
#ifndef RELEASE
    PushCallStack("basic::internal::LocalTriangularRank2K");
    CheckInput
    ( orientationOfA1, orientationOfB1, A1, A2, B1, B2, C );
#endif
    const Grid& g = C.Grid();

    if( C.Height() < g.Width()*LocalTriangularRank2KBlocksize<T>() )
    {
        LocalTriangularRank2KKernel
        ( shape, orientationOfA1, orientationOfB1,  
          alpha, A1, A2, B1, B2, beta, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.

        DistMatrix<T,STAR,MC> A1L(g), A1R(g);

        DistMatrix<T,MC,STAR> A2T(g),
                              A2B(g); 

        DistMatrix<T,MR,STAR> B1T(g),
                              B1B(g);

        DistMatrix<T,STAR,MR> B2L(g), B2R(g);

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

        if( shape == LOWER )
        { 
            basic::internal::LocalGemm
            ( orientationOfA1, NORMAL, alpha, A1R, B2L, beta, CBL );

            basic::internal::LocalGemm
            ( NORMAL, orientationOfB1, alpha, A2B, B1T, (T)1, CBL );
        }
        else
        {
            basic::internal::LocalGemm
            ( orientationOfA1, NORMAL, alpha, A1L, B2R, beta, CTR );

            basic::internal::LocalGemm
            ( NORMAL, orientationOfB1, alpha, A2T, B1B, (T)1, CTR );
        }

        // Recurse
        basic::internal::LocalTriangularRank2K
        ( shape, orientationOfA1, orientationOfB1,
          alpha, A1L, A2T, B1T, B2L, beta, CTL );

        basic::internal::LocalTriangularRank2K
        ( shape, orientationOfA1, orientationOfB1,
          alpha, A1R, A2B, B1B, B2R, beta, CBR );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::basic::internal::LocalTriangularRank2K
( Shape shape,
  Orientation orientationOfA2,
  Orientation orientationOfB2,
  T alpha, const DistMatrix<T,MC,  STAR>& A1,
           const DistMatrix<T,STAR,MC  >& A2,
           const DistMatrix<T,STAR,MR  >& B1,
           const DistMatrix<T,MR,  STAR>& B2,
  T beta,        DistMatrix<T,MC,  MR  >& C  )
{
#ifndef RELEASE
    PushCallStack("basic::internal::LocalTriangularRank2K");
    CheckInput
    ( orientationOfA2, orientationOfB2, A1, A2, B1, B2, C );
#endif
    const Grid& g = C.Grid();

    if( C.Height() < g.Width()*LocalTriangularRank2KBlocksize<T>() )
    {
        LocalTriangularRank2KKernel
        ( shape, orientationOfA2, orientationOfB2, 
          alpha, A1, A2, B1, B2, beta, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.

        DistMatrix<T,MC,STAR> A1T(g),
                              A1B(g); 

        DistMatrix<T,STAR,MC> A2L(g), A2R(g);

        DistMatrix<T,STAR,MR> B1L(g), B1R(g);

        DistMatrix<T,MR,STAR> B2T(g),
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

        if( shape == LOWER )
        { 
            basic::internal::LocalGemm
            ( NORMAL, orientationOfB2, alpha, A1B, B2T, beta, CBL );

            basic::internal::LocalGemm
            ( orientationOfA2, NORMAL, alpha, A2R, B1L, (T)1, CBL );
        }
        else
        {
            basic::internal::LocalGemm
            ( NORMAL, orientationOfB2, alpha, A1T, B2B, beta, CTR );

            basic::internal::LocalGemm
            ( orientationOfA2, NORMAL, alpha, A2L, B1R, (T)1, CTR );
        }

        // Recurse
        basic::internal::LocalTriangularRank2K
        ( shape, orientationOfA2, orientationOfB2,
          alpha, A1T, A2L, B1L, B2T, beta, CTL );

        basic::internal::LocalTriangularRank2K
        ( shape, orientationOfA2, orientationOfB2,
          alpha, A1B, A2R, B1R, B2B, beta, CBR );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::basic::internal::LocalTriangularRank2K
( Shape shape,
  Orientation orientationOfA2,
  Orientation orientationOfB1,
  T alpha, const DistMatrix<T,MC,  STAR>& A1,
           const DistMatrix<T,STAR,MC  >& A2,
           const DistMatrix<T,MR,  STAR>& B1,
           const DistMatrix<T,STAR,MR  >& B2,
  T beta,        DistMatrix<T,MC,  MR  >& C  )
{
#ifndef RELEASE
    PushCallStack("basic::internal::LocalTriangularRank2K");
    CheckInput
    ( orientationOfA2, orientationOfB1, A1, A2, B1, B2, C );
#endif
    const Grid& g = C.Grid();

    if( C.Height() < g.Width()*LocalTriangularRank2KBlocksize<T>() )
    {
        LocalTriangularRank2KKernel
        ( shape, orientationOfA2, orientationOfB1, 
          alpha, A1, A2, B1, B2, beta, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.

        DistMatrix<T,MC,STAR> A1T(g),
                              A1B(g); 

        DistMatrix<T,STAR,MC> A2L(g), A2R(g);

        DistMatrix<T,MR,STAR> B1T(g),
                              B1B(g);

        DistMatrix<T,STAR,MR> B2L(g), B2R(g);

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

        if( shape == LOWER )
        { 
            basic::internal::LocalGemm
            ( NORMAL, NORMAL, alpha, A1B, B2L, beta, CBL );

            basic::internal::LocalGemm
            ( orientationOfA2, orientationOfB1, alpha, A2R, B1T, (T)1, CBL );
        }
        else
        {
            basic::internal::LocalGemm
            ( NORMAL, NORMAL, alpha, A1T, B2R, beta, CTR );

            basic::internal::LocalGemm
            ( orientationOfA2, orientationOfB1, alpha, A2L, B1B, (T)1, CTR );
        }

        // Recurse
        basic::internal::LocalTriangularRank2K
        ( shape, orientationOfA2, orientationOfB1, 
          alpha, A1T, A2L, B1T, B2L, beta, CTL );

        basic::internal::LocalTriangularRank2K
        ( shape, orientationOfA2, orientationOfB1,
          alpha, A1B, A2R, B1B, B2R, beta, CBR );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::basic::internal::LocalTriangularRank2K
( Shape shape,
  T alpha, const DistMatrix<T,MC,  STAR>& A1,
           const DistMatrix<T,MC,  STAR>& A2,
           const DistMatrix<T,STAR,MR  >& B1,
           const DistMatrix<T,STAR,MR  >& B2,
  T beta,        DistMatrix<T,MC,  MR  >& C )
{
#ifndef RELEASE
    PushCallStack("basic::internal::LocalTriangularRank2K");
    CheckInput( A1, A2, B1, B2, C );
#endif
    const Grid& g = C.Grid();

    if( C.Height() < g.Width()*LocalTriangularRank2KBlocksize<T>() )
    {
        LocalTriangularRank2KKernel
        ( shape, alpha, A1, A2, B1, B2, beta, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.

        DistMatrix<T,MC,STAR> A1T(g),  A2T(g),
                              A1B(g),  A2B(g);

        DistMatrix<T,STAR,MR> B1L(g), B1R(g),
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

        if( shape == LOWER )
        { 
            basic::internal::LocalGemm
            ( NORMAL, NORMAL, alpha, A1B, B2L, beta, CBL );

            basic::internal::LocalGemm
            ( NORMAL, NORMAL, alpha, A2B, B1L, (T)1, CBL );
        }
        else
        {
            basic::internal::LocalGemm
            ( NORMAL, NORMAL, alpha, A1T, B2R, beta, CTR );
            
            basic::internal::LocalGemm
            ( NORMAL, NORMAL, alpha, A2T, B1R, (T)1, CTR );
        }

        // Recurse
        basic::internal::LocalTriangularRank2K
        ( shape, alpha, A1T, A2T, B1L, B2L, beta, CTL );

        basic::internal::LocalTriangularRank2K
        ( shape, alpha, A1B, A2B, B1R, B2R, beta, CBR );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::basic::internal::LocalTriangularRank2K
( Shape shape,
  Orientation orientationOfA1,
  Orientation orientationOfA2,
  Orientation orientationOfB2,
  T alpha, const DistMatrix<T,STAR,MC  >& A1,
           const DistMatrix<T,STAR,MC  >& A2,
           const DistMatrix<T,STAR,MR  >& B1,
           const DistMatrix<T,MR,  STAR>& B2,
  T beta,        DistMatrix<T,MC,  MR>& C  )
{
#ifndef RELEASE
    PushCallStack("basic::internal::LocalTriangularRank2K");
    CheckInput
    ( orientationOfA1, orientationOfA2, orientationOfB2, A1, A2, B1, B2, C );
#endif
    const Grid& g = C.Grid();

    if( C.Height() < g.Width()*LocalTriangularRank2KBlocksize<T>() )
    {
        LocalTriangularRank2KKernel
        ( shape, orientationOfA1, orientationOfA2, orientationOfB2,
          alpha, A1, A2, B1, B2, beta, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.

        DistMatrix<T,STAR,MC> A1L(g), A1R(g),
                              A2L(g), A2R(g);

        DistMatrix<T,STAR,MR> B1L(g), B1R(g);

        DistMatrix<T,MR,STAR> B2T(g),
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

        if( shape == LOWER )
        { 
            basic::internal::LocalGemm
            ( orientationOfA1, orientationOfB2, alpha, A1R, B2T, beta, CBL );

            basic::internal::LocalGemm
            ( orientationOfA2, NORMAL, alpha, A2R, B1L, (T)1, CBL );
        }
        else
        {
            basic::internal::LocalGemm
            ( orientationOfA1, orientationOfB2, alpha, A1L, B2B, beta, CTR );

            basic::internal::LocalGemm
            ( orientationOfA2, NORMAL, alpha, A2L, B1R, (T)1, CTR );
        }

        // Recurse
        basic::internal::LocalTriangularRank2K
        ( shape, orientationOfA1, orientationOfA2, orientationOfB2,
          alpha, A1L, A2L, B1L, B2T, beta, CTL );

        basic::internal::LocalTriangularRank2K
        ( shape, orientationOfA1, orientationOfA2, orientationOfB2,
          alpha, A1R, A2R, B1R, B2B, beta, CBR );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::basic::internal::LocalTriangularRank2K
( Shape shape,
  Orientation orientationOfA1,
  Orientation orientationOfA2,
  Orientation orientationOfB1,
  T alpha, const DistMatrix<T,STAR,MC  >& A1,
           const DistMatrix<T,STAR,MC  >& A2,
           const DistMatrix<T,MR,  STAR>& B1,
           const DistMatrix<T,STAR,MR  >& B2,
  T beta,        DistMatrix<T,MC,  MR  >& C  )
{
#ifndef RELEASE
    PushCallStack("basic::internal::LocalTriangularRank2K");
    CheckInput
    ( orientationOfA1, orientationOfA2, orientationOfB1, A1, A2, B1, B2, C );
#endif
    const Grid& g = C.Grid();

    if( C.Height() < g.Width()*LocalTriangularRank2KBlocksize<T>() )
    {
        LocalTriangularRank2KKernel
        ( shape, orientationOfA1, orientationOfA2, orientationOfB1,
          alpha, A1, A2, B1, B2, beta, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.

        DistMatrix<T,STAR,MC> A1L(g), A1R(g),
                              A2L(g), A2R(g);

        DistMatrix<T,MR,STAR> B1T(g),
                              B1B(g);

        DistMatrix<T,STAR,MR> B2L(g), B2R(g);

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

        if( shape == LOWER )
        { 
            basic::internal::LocalGemm
            ( orientationOfA1, NORMAL, alpha, A1R, B2L, beta, CBL );

            basic::internal::LocalGemm
            ( orientationOfA2, orientationOfB1, alpha, A2R, B1T, (T)1, CBL );
        }
        else
        {
            basic::internal::LocalGemm
            ( orientationOfA1, NORMAL, alpha, A1L, B2R, beta, CTR );

            basic::internal::LocalGemm
            ( orientationOfA2, orientationOfB1, alpha, A2L, B1B, (T)1, CTR );
        }

        // Recurse
        basic::internal::LocalTriangularRank2K
        ( shape, orientationOfA1, orientationOfA2, orientationOfB1,
          alpha, A1L, A2L, B1T, B2L, beta, CTL );

        basic::internal::LocalTriangularRank2K
        ( shape, orientationOfA1, orientationOfA2, orientationOfB1,
          alpha, A1R, A2R, B1B, B2R, beta, CBR );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::basic::internal::LocalTriangularRank2K
( Shape shape,
  Orientation orientationOfA1,
  T alpha, const DistMatrix<T,STAR,MC  >& A1,
           const DistMatrix<T,MC,  STAR>& A2,
           const DistMatrix<T,STAR,MR  >& B1,
           const DistMatrix<T,STAR,MR  >& B2,
  T beta,        DistMatrix<T,MC,  MR  >& C  )
{
#ifndef RELEASE
    PushCallStack("basic::internal::LocalTriangularRank2K");
    CheckInput( orientationOfA1, A1, A2, B1, B2, C );
#endif
    const Grid& g = C.Grid();

    if( C.Height() < g.Width()*LocalTriangularRank2KBlocksize<T>() )
    {
        LocalTriangularRank2KKernel
        ( shape, orientationOfA1, alpha, A1, A2, B1, B2, beta, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.

        DistMatrix<T,STAR,MC> A1L(g), A1R(g);

        DistMatrix<T,MC,STAR> A2T(g),
                              A2B(g);

        DistMatrix<T,STAR,MR> B1L(g), B1R(g),
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

        if( shape == LOWER )
        { 
            basic::internal::LocalGemm
            ( orientationOfA1, NORMAL, alpha, A1R, B2L, beta, CBL );

            basic::internal::LocalGemm
            ( NORMAL, NORMAL, alpha, A2B, B1L, (T)1, CBL );
        }
        else
        {
            basic::internal::LocalGemm
            ( orientationOfA1, NORMAL, alpha, A1L, B2R, beta, CTR );

            basic::internal::LocalGemm
            ( NORMAL, NORMAL, alpha, A2T, B1R, (T)1, CTR );
        }

        // Recurse
        basic::internal::LocalTriangularRank2K
        ( shape, orientationOfA1, alpha, A1L, A2T, B1L, B2L, beta, CTL );

        basic::internal::LocalTriangularRank2K
        ( shape, orientationOfA1, alpha, A1R, A2B, B1R, B2R, beta, CBR );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::basic::internal::LocalTriangularRank2K
( Shape shape,
  Orientation orientationOfA2,
  T alpha, const DistMatrix<T,MC,  STAR>& A1,
           const DistMatrix<T,STAR,MC  >& A2,
           const DistMatrix<T,STAR,MR  >& B1,
           const DistMatrix<T,STAR,MR  >& B2,
  T beta,        DistMatrix<T,MC,  MR  >& C  )
{
#ifndef RELEASE
    PushCallStack("basic::internal::LocalTriangularRank2K");
    CheckInput( orientationOfA2, A1, A2, B1, B2, C );
#endif
    const Grid& g = C.Grid();

    if( C.Height() < g.Width()*LocalTriangularRank2KBlocksize<T>() )
    {
        LocalTriangularRank2KKernel
        ( shape, orientationOfA2, alpha, A1, A2, B1, B2, beta, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.

        DistMatrix<T,MC,STAR> A1T(g),
                              A1B(g);

        DistMatrix<T,STAR,MC> A2L(g), A2R(g);

        DistMatrix<T,STAR,MR> B1L(g), B1R(g),
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

        if( shape == LOWER )
        { 
            basic::internal::LocalGemm
            ( NORMAL, NORMAL, alpha, A1B, B2L, beta, CBL );

            basic::internal::LocalGemm
            ( orientationOfA2, NORMAL, alpha, A2R, B1L, (T)1, CBL );
        }
        else
        {
            basic::internal::LocalGemm
            ( NORMAL, NORMAL, alpha, A1T, B2R, beta, CTR );

            basic::internal::LocalGemm
            ( orientationOfA2, NORMAL, alpha, A2L, B1R, (T)1, CTR );
        }

        // Recurse
        basic::internal::LocalTriangularRank2K
        ( shape, orientationOfA2, alpha, A1T, A2L, B1L, B2L, beta, CTL );

        basic::internal::LocalTriangularRank2K
        ( shape, orientationOfA2, alpha, A1B, A2R, B1R, B2R, beta, CBR );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::basic::internal::LocalTriangularRank2K
( Shape shape,
  Orientation orientationOfA1,
  Orientation orientationOfA2,
  T alpha, const DistMatrix<T,STAR,MC>& A1,
           const DistMatrix<T,STAR,MC>& A2,
           const DistMatrix<T,STAR,MR>& B1,
           const DistMatrix<T,STAR,MR>& B2,
  T beta,        DistMatrix<T,MC,  MR>& C  )
{
#ifndef RELEASE
    PushCallStack("basic::internal::LocalTriangularRank2K");
    CheckInput( orientationOfA1, orientationOfA2, A1, A2, B1, B2, C );
#endif
    const Grid& g = C.Grid();

    if( C.Height() < g.Width()*LocalTriangularRank2KBlocksize<T>() )
    {
        LocalTriangularRank2KKernel
        ( shape, orientationOfA1, orientationOfA2, 
          alpha, A1, A2, B1, B2, beta, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.

        DistMatrix<T,STAR,MC> A1L(g), A1R(g),
                              A2L(g), A2R(g);

        DistMatrix<T,STAR,MR> B1L(g), B1R(g),
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

        if( shape == LOWER )
        { 
            basic::internal::LocalGemm
            ( orientationOfA1, NORMAL, alpha, A1R, B2L, beta, CBL );

            basic::internal::LocalGemm
            ( orientationOfA2, NORMAL, alpha, A2R, B1L, (T)1, CBL );
        }
        else
        {
            basic::internal::LocalGemm
            ( orientationOfA1, NORMAL, alpha, A1L, B2R, beta, CTR );

            basic::internal::LocalGemm
            ( orientationOfA2, NORMAL, alpha, A2L, B1R, (T)1, CTR );
        }

        // Recurse
        basic::internal::LocalTriangularRank2K
        ( shape, orientationOfA1, orientationOfA2,
          alpha, A1L, A2L, B1L, B2L, beta, CTL );

        basic::internal::LocalTriangularRank2K
        ( shape, orientationOfA1, orientationOfA2, 
          alpha, A1R, A2R, B1R, B2R, beta, CBR );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

