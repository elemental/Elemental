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

namespace elemental {
namespace basic {
namespace trr2k {

#ifndef RELEASE
// E := alpha (A B + C D) + beta E
template<typename T>
inline void 
CheckInput
( const DistMatrix<T,MC,  STAR>& A, const DistMatrix<T,STAR,MR>& B, 
  const DistMatrix<T,MC,  STAR>& C, const DistMatrix<T,STAR,MR>& D,
  const DistMatrix<T,MC,  MR  >& E  )
{
    if( A.Grid() != B.Grid() || B.Grid() != C.Grid() ||
        C.Grid() != D.Grid() || D.Grid() != E.Grid() )
        throw std::logic_error
        ("A, B, C, D, and E must be distributed over the same grid");
    if( A.Height() != E.Height() || B.Width() != E.Width() ||
        C.Height() != E.Height() || D.Width() != E.Width() ||
        A.Width()  != B.Height() || C.Width() != D.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal LocalTrr2k: \n"
            << "  A[MC,* ] ~ " << A.Height() << " x "
                               << A.Width()  << "\n"
            << "  B[* ,MR] ~ " << B.Height() << " x "
                               << B.Width()  << "\n"
            << "  C[MC,* ] ~ " << C.Height() << " x "
                               << C.Width()  << "\n"
            << "  D[* ,MR] ~ " << D.Height() << " x "
                               << D.Width()  << "\n"
            << "  E[MC,MR] ~ " << E.Height() << " x " << E.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( A.ColAlignment() != E.ColAlignment() ||
        B.RowAlignment() != E.RowAlignment() ||
        C.ColAlignment() != E.ColAlignment() ||
        D.RowAlignment() != E.RowAlignment() )
    {
        std::ostringstream msg;
        msg << "Misaligned LocalTrr2k: \n"
            << "  A[MC,* ] ~ " << A.ColAlignment() << "\n"
            << "  B[* ,MR] ~ " << B.RowAlignment() << "\n"
            << "  C[MC,* ] ~ " << C.ColAlignment() << "\n"
            << "  D[* ,MR] ~ " << D.RowAlignment() << "\n"
            << "  E[MC,MR] ~ " << E.ColAlignment() << " , " <<
                                  E.RowAlignment() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
}

// E := alpha (A B + C D^{T/H}) + beta E
template<typename T>
inline void 
CheckInput
( Orientation orientationOfD,
  const DistMatrix<T,MC,  STAR>& A, const DistMatrix<T,STAR,MR>& B,
  const DistMatrix<T,MC,  STAR>& C, const DistMatrix<T,MR,STAR>& D,
  const DistMatrix<T,MC,  MR  >& E )
{
    if( orientationOfD == NORMAL )
        throw std::logic_error("D[MR,* ] must be (Conjugate)Transpose'd");
    if( A.Grid() != B.Grid() || B.Grid() != C.Grid() ||
        C.Grid() != D.Grid() || D.Grid() != E.Grid() )
        throw std::logic_error
        ("A, B, C, D, and E must be distributed over the same grid");
    if( A.Height() != E.Height() || B.Width()  != E.Width() ||
        C.Height() != E.Height() || D.Height() != E.Width() ||
        A.Width()  != B.Height() || C.Width()  != D.Width() )
    {
        std::ostringstream msg;
        msg << "Nonconformal LocalTrr2k: \n"
            << "  A[MC,* ] ~ " << A.Height() << " x "
                               << A.Width()  << "\n"
            << "  B[* ,MR] ~ " << B.Height() << " x "
                               << B.Width()  << "\n"
            << "  C[MC,* ] ~ " << C.Height() << " x "
                               << C.Width()  << "\n"
            << "  D[MR,* ] ~ " << D.Height() << " x "
                               << D.Width()  << "\n"
            << "  E[MC,MR] ~ " << E.Height() << " x " << E.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( A.ColAlignment() != E.ColAlignment() ||
        B.RowAlignment() != E.RowAlignment() ||
        C.ColAlignment() != E.ColAlignment() ||
        D.ColAlignment() != E.RowAlignment() )
    {
        std::ostringstream msg;
        msg << "Misaligned LocalTrr2k: \n"
            << "  A[MC,* ] ~ " << A.ColAlignment() << "\n"
            << "  B[* ,MR] ~ " << B.RowAlignment() << "\n"
            << "  C[MC,* ] ~ " << C.ColAlignment() << "\n"
            << "  D[MR,* ] ~ " << D.ColAlignment() << "\n"
            << "  E[MC,MR] ~ " << E.ColAlignment() << " , " <<
                                  E.RowAlignment() << "\n";
        throw std::logic_error( msg.str() );
    }
}

// E := alpha (A B + C^{T/H} D) + beta E
template<typename T>
inline void 
CheckInput
( Orientation orientationOfC,
  const DistMatrix<T,MC,  STAR>& A, const DistMatrix<T,STAR,MR>& B,
  const DistMatrix<T,STAR,MC  >& C, const DistMatrix<T,STAR,MR>& D,
  const DistMatrix<T,MC,  MR  >& E )
{
    if( orientationOfC == NORMAL )
        throw std::logic_error("C[* ,MC] must be (Conjugate)Transpose'd");
    if( A.Grid() != B.Grid() || B.Grid() != C.Grid() ||
        C.Grid() != D.Grid() || D.Grid() != E.Grid() )
        throw std::logic_error
        ("A, B, C, D, and E must be distributed over the same grid");
    if( A.Height() != E.Height() || B.Width()  != E.Width() ||
        C.Width()  != E.Height() || D.Width()  != E.Width() ||
        A.Width()  != B.Height() || C.Height() != D.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal LocalTrr2k: \n"
            << "  A[MC,* ] ~ " << A.Height() << " x "
                               << A.Width()  << "\n"
            << "  B[* ,MR] ~ " << B.Height() << " x "
                               << B.Width()  << "\n"
            << "  C[* ,MC] ~ " << C.Height() << " x "
                               << C.Width()  << "\n"
            << "  D[* ,MR] ~ " << D.Height() << " x "
                               << D.Width()  << "\n"
            << "  E[MC,MR] ~ " << E.Height() << " x " << E.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( A.ColAlignment() != E.ColAlignment() ||
        B.RowAlignment() != E.RowAlignment() ||
        C.RowAlignment() != E.ColAlignment() ||
        D.RowAlignment() != E.RowAlignment() )
    {
        std::ostringstream msg;
        msg << "Misaligned LocalTrr2k: \n"
            << "  A[MC,* ] ~ " << A.ColAlignment() << "\n"
            << "  B[* ,MR] ~ " << B.RowAlignment() << "\n"
            << "  C[* ,MC] ~ " << C.RowAlignment() << "\n"
            << "  D[* ,MR] ~ " << D.RowAlignment() << "\n"
            << "  E[MC,MR] ~ " << E.ColAlignment() << " , " <<
                                  E.RowAlignment() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
}

// E := alpha (A B + C^{T/H} D^{T/H}) + beta E
template<typename T>
inline void 
CheckInput
( Orientation orientationOfC, Orientation orientationOfD,
  const DistMatrix<T,MC,  STAR>& A, const DistMatrix<T,STAR,MR>& B,
  const DistMatrix<T,STAR,MC  >& C, const DistMatrix<T,MR,STAR>& D,
  const DistMatrix<T,MC,  MR  >& E )
{
    if( orientationOfC == NORMAL )
        throw std::logic_error("C[* ,MC] must be (Conjugate)Transpose'd");
    if( orientationOfD == NORMAL )
        throw std::logic_error("D[MR,* ] must be (Conjugate)Transpose'd");
    if( A.Grid() != B.Grid() || B.Grid() != C.Grid() ||
        C.Grid() != D.Grid() || D.Grid() != E.Grid() )
        throw std::logic_error
        ("A, B, C, D, and E must be distributed over the same grid");
    if( A.Height() != E.Height() || B.Width()  != E.Width() ||
        C.Width()  != E.Height() || D.Height() != E.Width() ||
        A.Width()  != B.Height() || C.Height() != D.Width() )
    {
        std::ostringstream msg;
        msg << "Nonconformal LocalTrr2k: \n"
            << "  A[MC,* ] ~ " << A.Height() << " x "
                               << A.Width()  << "\n"
            << "  B[* ,MR] ~ " << B.Height() << " x "
                               << B.Width()  << "\n"
            << "  C[* ,MC] ~ " << C.Height() << " x "
                               << C.Width()  << "\n"
            << "  D[MR,* ] ~ " << D.Height() << " x "
                               << D.Width()  << "\n"
            << "  E[MC,MR] ~ " << E.Height() << " x " << E.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( A.ColAlignment() != E.ColAlignment() ||
        B.RowAlignment() != E.RowAlignment() ||
        C.RowAlignment() != E.ColAlignment() ||
        D.ColAlignment() != E.RowAlignment() )
    {
        std::ostringstream msg;
        msg << "Misaligned LocalTrr2k: \n"
            << "  A[MC,* ] ~ " << A.ColAlignment() << "\n"
            << "  B[* ,MR] ~ " << B.RowAlignment() << "\n"
            << "  C[* ,MC] ~ " << C.RowAlignment() << "\n"
            << "  D[MR,* ] ~ " << D.ColAlignment() << "\n"
            << "  E[MC,MR] ~ " << E.ColAlignment() << " , " <<
                                  E.RowAlignment() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
}

// E := alpha (A B^{T/H} + C D) + beta E
template<typename T>
inline void 
CheckInput
( Orientation orientationOfB,
  const DistMatrix<T,MC,  STAR>& A, const DistMatrix<T,MR,STAR>& B,
  const DistMatrix<T,MC,  STAR>& C, const DistMatrix<T,STAR,MR>& D,
  const DistMatrix<T,MC,  MR  >& E )
{
    if( orientationOfB == NORMAL )
        throw std::logic_error("B[MR,* ] must be (Conjugate)Transpose'd");
    if( A.Grid() != B.Grid() || B.Grid() != C.Grid() ||
        C.Grid() != D.Grid() || D.Grid() != E.Grid() )
        throw std::logic_error
        ("A, B, C, D, and E must be distributed over the same grid");
    if( A.Height() != E.Height() || B.Height() != E.Width() ||
        C.Height() != E.Height() || D.Width()  != E.Width() ||
        A.Width()  != B.Width()  || C.Width()  != D.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal LocalTrr2k: \n"
            << "  A[MC,* ] ~ " << A.Height() << " x "
                               << A.Width()  << "\n"
            << "  B[MR,* ] ~ " << B.Height() << " x "
                               << B.Width()  << "\n"
            << "  C[MC,* ] ~ " << C.Height() << " x "
                               << C.Width()  << "\n"
            << "  D[* ,MR] ~ " << D.Height() << " x "
                               << D.Width()  << "\n"
            << "  E[MC,MR] ~ " << E.Height() << " x " << E.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( A.ColAlignment() != E.ColAlignment() ||
        B.ColAlignment() != E.RowAlignment() ||
        C.ColAlignment() != E.ColAlignment() ||
        D.RowAlignment() != E.RowAlignment() )
    {
        std::ostringstream msg;
        msg << "Misaligned LocalTrr2k: \n"
            << "  A[MC,* ] ~ " << A.ColAlignment() << "\n"
            << "  B[MR,* ] ~ " << B.ColAlignment() << "\n"
            << "  C[MC,* ] ~ " << C.ColAlignment() << "\n"
            << "  D[* ,MR] ~ " << D.RowAlignment() << "\n"
            << "  E[MC,MR] ~ " << E.ColAlignment() << " , " <<
                                  E.RowAlignment() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
}

// E := alpha (A B^{T/H} + C D^{T/H}) + beta E
template<typename T>
inline void 
CheckInput
( Orientation orientationOfB,
  Orientation orientationOfD,
  const DistMatrix<T,MC,STAR>& A, const DistMatrix<T,MR,STAR>& B, 
  const DistMatrix<T,MC,STAR>& C, const DistMatrix<T,MR,STAR>& D,
  const DistMatrix<T,MC,MR  >& E )
{
    if( orientationOfB == NORMAL )
        throw std::logic_error("B[MR,* ] must be (Conjugate)Transpose'd");
    if( orientationOfD == NORMAL )
        throw std::logic_error("D[MR,* ] must be (Conjugate)Transpose'd");
    if( A.Grid() != B.Grid() || B.Grid() != C.Grid() ||
        C.Grid() != D.Grid() || D.Grid() != E.Grid() )
        throw std::logic_error
        ("A, B, C, D, and E must be distributed over the same grid");
    if( A.Height() != E.Height() || B.Height() != E.Width()  ||
        A.Height() != C.Height() || A.Width()  != C.Width()  ||
        B.Width()  != D.Width()  || B.Height() != D.Height() ||
        A.Width()  != B.Width()  || C.Width()  != D.Width() )
    {
        std::ostringstream msg;
        msg << "Nonconformal LocalTrr2k: \n"
            << "  A[MC,* ] ~ " << A.Height() << " x "
                               << A.Width()  << "\n"
            << "  B[MR,* ] ~ " << B.Height() << " x "
                               << B.Width()  << "\n"
            << "  C[MC,* ] ~ " << C.Height() << " x "
                               << C.Width()  << "\n"
            << "  D[MR,* ] ~ " << D.Height() << " x "
                               << D.Width()  << "\n"
            << "  E[MC,MR] ~ " << E.Height() << " x " << E.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( A.ColAlignment() != E.ColAlignment() ||
        B.ColAlignment() != E.RowAlignment() ||
        A.ColAlignment() != C.ColAlignment() ||
        B.ColAlignment() != D.ColAlignment() )
    {
        std::ostringstream msg;
        msg << "Misaligned LocalTrr2k: \n"
            << "  A[MC,* ] ~ " << A.ColAlignment() << "\n"
            << "  B[MR,* ] ~ " << B.ColAlignment() << "\n"
            << "  C[MC,* ] ~ " << C.ColAlignment() << "\n"
            << "  D[MR,* ] ~ " << D.ColAlignment() << "\n"
            << "  E[MC,MR] ~ " << E.ColAlignment() << " , " <<
                                  E.RowAlignment() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
}

// E := alpha (A B^{T/H} + C^{T/H} D) + beta E
template<typename T>
inline void 
CheckInput
( Orientation orientationOfB,
  Orientation orientationOfC,
  const DistMatrix<T,MC,  STAR>& A, const DistMatrix<T,MR,STAR>& B,
  const DistMatrix<T,STAR,MC  >& C, const DistMatrix<T,STAR,MR>& D,
  const DistMatrix<T,MC,  MR  >& E )
{
    if( orientationOfB == NORMAL )
        throw std::logic_error("B[MR,* ] must be (Conjugate)Transpose'd");
    if( orientationOfC == NORMAL )
        throw std::logic_error("C[* ,MC] must be (Conjugate)Transpose'd");
    if( A.Grid() != B.Grid() || B.Grid() != C.Grid() ||
        C.Grid() != D.Grid() || D.Grid() != E.Grid() )
        throw std::logic_error
        ("A, B, C, D, and E must be distributed over the same grid");
    if( A.Height() != E.Height() || B.Height() != E.Width() ||
        C.Width()  != E.Height() || D.Width()  != E.Width() ||
        A.Width()  != B.Width()  || C.Height() != D.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal LocalTrr2k: \n"
            << "  A[MC,* ] ~ " << A.Height() << " x "
                               << A.Width()  << "\n"
            << "  B[MR,* ] ~ " << B.Height() << " x "
                               << B.Width()  << "\n"
            << "  C[* ,MC] ~ " << C.Height() << " x "
                               << C.Width()  << "\n"
            << "  D[* ,MR] ~ " << D.Height() << " x "
                                << D.Width()  << "\n"
            << "  E[MC,MR] ~ " << E.Height() << " x " << E.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( A.ColAlignment() != E.ColAlignment() ||
        B.ColAlignment() != E.RowAlignment() ||
        C.RowAlignment() != E.ColAlignment() ||
        D.RowAlignment() != E.RowAlignment() )
    {
        std::ostringstream msg;
        msg << "Misaligned LocalTrr2k: \n"
            << "  A[MC,* ] ~ " << A.ColAlignment() << "\n"
            << "  B[MR,* ] ~ " << B.ColAlignment() << "\n"
            << "  C[* ,MC] ~ " << C.RowAlignment() << "\n"
            << "  D[* ,MR] ~ " << D.RowAlignment() << "\n"
            << "  E[MC,MR] ~ " << E.ColAlignment() << " , " <<
                                  E.RowAlignment() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
}

// E := alpha (A B^{T/H} + C^{T/H} D^{T/H}) + beta E
template<typename T>
inline void 
CheckInput
( Orientation orientationOfB,
  Orientation orientationOfC,
  Orientation orientationOfD,
  const DistMatrix<T,MC,  STAR>& A, const DistMatrix<T,MR,STAR>& B,
  const DistMatrix<T,STAR,MC  >& C, const DistMatrix<T,MR,STAR>& D,
  const DistMatrix<T,MC,  MR  >& E )
{
    if( orientationOfB == NORMAL )
        throw std::logic_error("B[MR,* ] must be (Conjugate)Transpose'd");
    if( orientationOfC == NORMAL )
        throw std::logic_error("C[* ,MC] must be (Conjugate)Transpose'd");
    if( orientationOfD == NORMAL )
        throw std::logic_error("D[MR,* ] must be (Conjugate)Transpose'd");
    if( A.Grid() != B.Grid() || B.Grid() != C.Grid() ||
        C.Grid() != D.Grid() || D.Grid() != E.Grid() )
        throw std::logic_error
        ("A, B, C, D, and E must be distributed over the same grid");
    if( A.Height() != E.Height() || B.Height() != E.Width() ||
        C.Width()  != E.Height() || D.Height() != E.Width() ||
        A.Width()  != B.Width()  || C.Height() != D.Width() )
    {
        std::ostringstream msg;
        msg << "Nonconformal LocalTrr2k: \n"
            << "  A[MC,* ] ~ " << A.Height() << " x "
                               << A.Width()  << "\n"
            << "  B[MR,* ] ~ " << B.Height() << " x "
                               << B.Width()  << "\n"
            << "  C[* ,MC] ~ " << C.Height() << " x "
                               << C.Width()  << "\n"
            << "  D[MR,* ] ~ " << D.Height() << " x "
                               << D.Width()  << "\n"
            << "  E[MC,MR] ~ " << E.Height() << " x " << E.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( A.ColAlignment() != E.ColAlignment() ||
        B.ColAlignment() != E.RowAlignment() ||
        C.RowAlignment() != E.ColAlignment() ||
        D.ColAlignment() != E.RowAlignment() )
    {
        std::ostringstream msg;
        msg << "Misaligned LocalTrr2k: \n"
            << "  A[MC,* ] ~ " << A.ColAlignment() << "\n"
            << "  B[MR,* ] ~ " << B.ColAlignment() << "\n"
            << "  C[* ,MC] ~ " << C.RowAlignment() << "\n"
            << "  D[MR,* ] ~ " << D.ColAlignment() << "\n"
            << "  E[MC,MR] ~ " << E.ColAlignment() << " , " <<
                                  E.RowAlignment() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
}

// E := alpha (A^{T/H} B + C D) + beta E
template<typename T>
inline void 
CheckInput
( Orientation orientationOfA,
  const DistMatrix<T,STAR,MC  >& A, const DistMatrix<T,STAR,MR>& B, 
  const DistMatrix<T,MC,  STAR>& C, const DistMatrix<T,STAR,MR>& D,
  const DistMatrix<T,MC,  MR  >& E )
{
    if( orientationOfA == NORMAL )
        throw std::logic_error("A[* ,MC] must be (Conjugate)Transpose'd");
    if( A.Grid() != B.Grid() || B.Grid() != C.Grid() ||
        C.Grid() != D.Grid() || D.Grid() != E.Grid() )
        throw std::logic_error
        ("A, B, C, D, and E must be distributed over the same grid");
    if( A.Width()  != E.Height() || B.Width() != E.Width() ||
        C.Height() != E.Height() || D.Width() != E.Width() || 
        A.Height() != B.Height() || C.Width() != D.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal LocalTrr2k: \n"
            << "  A[* ,MC] ~ " << A.Height() << " x "
                               << A.Width()  << "\n"
            << "  B[* ,MR] ~ " << B.Height() << " x "
                               << B.Width()  << "\n"
            << "  C[MC,* ] ~ " << C.Height() << " x "
                               << C.Width()  << "\n"
            << "  D[* ,MR] ~ " << D.Height() << " x "
                               << D.Width()  << "\n"
            << "  E[MC,MR] ~ " << E.Height() << " x " << E.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( A.RowAlignment() != E.ColAlignment() ||
        B.RowAlignment() != E.RowAlignment() ||
        C.ColAlignment() != E.ColAlignment() ||
        D.RowAlignment() != E.RowAlignment() )
    {
        std::ostringstream msg;
        msg << "Misaligned LocalTrr2k: \n"
            << "  A[* ,MC] ~ " << A.RowAlignment() << "\n"
            << "  B[* ,MR] ~ " << B.RowAlignment() << "\n"
            << "  C[MC,* ] ~ " << C.ColAlignment() << "\n"
            << "  D[* ,MR] ~ " << D.RowAlignment() << "\n"
            << "  E[MC,MR] ~ " << E.ColAlignment() << " , " <<
                                  E.RowAlignment() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
}

// E := alpha (A^{T/H} B + C D^{T/H}) + beta E
template<typename T>
inline void 
CheckInput
( Orientation orientationOfA,
  Orientation orientationOfD,
  const DistMatrix<T,STAR,MC  >& A, const DistMatrix<T,STAR,MR>& B,
  const DistMatrix<T,MC,  STAR>& C, const DistMatrix<T,MR,STAR>& D,
  const DistMatrix<T,MC,  MR  >& E )
{
    if( orientationOfA == NORMAL )
        throw std::logic_error("A[* ,MC] must be (Conjugate)Transpose'd");
    if( orientationOfD == NORMAL )
        throw std::logic_error("D[MR,* ] must be (Conjugate)Transpose'd");
    if( A.Grid() != B.Grid() || B.Grid() != C.Grid() ||
        C.Grid() != D.Grid() || D.Grid() != E.Grid() )
        throw std::logic_error
        ("A, B, C, D, and E must be distributed over the same grid");
    if( A.Width()  != E.Height() || B.Width()  != E.Width() ||
        C.Height() != E.Height() || D.Height() != E.Width() ||
        A.Height() != B.Height() || C.Width()  != D.Width() )
    {
        std::ostringstream msg;
        msg << "Nonconformal LocalTrr2k: \n"
            << "  A[* ,MC] ~ " << A.Height() << " x "
                               << A.Width()  << "\n"
            << "  C[MC,* ] ~ " << C.Height() << " x "
                               << C.Width()  << "\n"
            << "  D[MR,* ] ~ " << D.Height() << " x "
                               << D.Width()  << "\n"
            << "  B[* ,MR] ~ " << B.Height() << " x "
                               << B.Width()  << "\n"
            << "  E[MC,MR] ~ " << E.Height() << " x " << E.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( A.RowAlignment() != E.ColAlignment() ||
        B.RowAlignment() != E.RowAlignment() ||
        C.ColAlignment() != E.ColAlignment() ||
        D.ColAlignment() != E.RowAlignment() )
    {
        std::ostringstream msg;
        msg << "Misaligned LocalTrr2k: \n"
            << "  A[* ,MC] ~ " << A.RowAlignment() << "\n"
            << "  B[* ,MR] ~ " << B.RowAlignment() << "\n"
            << "  C[MC,* ] ~ " << C.ColAlignment() << "\n"
            << "  D[MR,* ] ~ " << D.ColAlignment() << "\n"
            << "  E[MC,MR] ~ " << E.ColAlignment() << " , " <<
                                  E.RowAlignment() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
}

// E := alpha (A^{T/H} B + C^{T/H} D) + beta E
template<typename T>
inline void 
CheckInput
( Orientation orientationOfA,
  Orientation orientationOfC,
  const DistMatrix<T,STAR,MC>& A, const DistMatrix<T,STAR,MR>& B, 
  const DistMatrix<T,STAR,MC>& C, const DistMatrix<T,STAR,MR>& D,
  const DistMatrix<T,MC,  MR>& E )
{
    if( orientationOfA == NORMAL )
        throw std::logic_error("A[* ,MC] must be (Conjugate)Transpose'd");
    if( orientationOfC == NORMAL )
        throw std::logic_error("C[* ,MC] must be (Conjugate)Transpose'd");
    if( A.Grid() != B.Grid() || B.Grid() != C.Grid() ||
        C.Grid() != D.Grid() || D.Grid() != E.Grid() )
        throw std::logic_error
        ("A, B, C, D, and E must be distributed over the same grid");
    if( A.Width()  != E.Height() || B.Width()  != E.Width() ||
        C.Width()  != E.Height() || D.Width()  != E.Width() ||
        A.Height() != B.Height() || C.Height() != D.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal LocalTrr2k: \n"
            << "  A[* ,MC] ~ " << A.Height() << " x "
                               << A.Width()  << "\n"
            << "  B[* ,MR] ~ " << B.Height() << " x "
                               << B.Width()  << "\n"
            << "  C[* ,MC] ~ " << C.Height() << " x "
                               << C.Width()  << "\n"
            << "  D[* ,MR] ~ " << D.Height() << " x "
                               << D.Width()  << "\n"
            << "  E[MC,MR] ~ " << E.Height() << " x " << E.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( A.RowAlignment() != E.ColAlignment() ||
        B.RowAlignment() != E.RowAlignment() ||
        C.RowAlignment() != E.ColAlignment() ||
        D.RowAlignment() != E.RowAlignment() )
    {
        std::ostringstream msg;
        msg << "Misaligned LocalTrr2k: \n"
            << "  A[* ,MC] ~ " << A.RowAlignment() << "\n"
            << "  B[* ,MR] ~ " << B.RowAlignment() << "\n"
            << "  C[* ,MC] ~ " << C.RowAlignment() << "\n"
            << "  D[* ,MR] ~ " << D.RowAlignment() << "\n"
            << "  E[MC,MR] ~ " << E.ColAlignment() << " , " <<
                                  E.RowAlignment() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
}

// E := alpha (A^{T/H} B + C^{T/H} D^{T/H}) + beta E
template<typename T>
inline void 
CheckInput
( Orientation orientationOfA,
  Orientation orientationOfC,
  Orientation orientationOfD,
  const DistMatrix<T,STAR,MC>& A, const DistMatrix<T,STAR,MR>& B,
  const DistMatrix<T,STAR,MC>& C, const DistMatrix<T,MR,STAR>& D,
  const DistMatrix<T,MC,  MR>& E )
{
    if( orientationOfA == NORMAL )
        throw std::logic_error("A[* ,MC] must be (Conjugate)Transpose'd");
    if( orientationOfC == NORMAL )
        throw std::logic_error("C[* ,MC] must be (Conjugate)Transpose'd");
    if( orientationOfD == NORMAL )
        throw std::logic_error("D[MR,* ] must be (Conjugate)Transpose'd");
    if( A.Grid() != B.Grid() || B.Grid() != C.Grid() ||
        C.Grid() != D.Grid() || D.Grid() != E.Grid() )
        throw std::logic_error
        ("A, B, C, D, and E must be distributed over the same grid");
    if( A.Width()  != E.Height() || B.Width()  != E.Width() ||
        C.Width()  != E.Height() || D.Height() != E.Width() ||
        A.Height() != B.Height() || C.Height() != D.Width() )
    {
        std::ostringstream msg;
        msg << "Nonconformal LocalTrr2k: \n"
            << "  A[* ,MC] ~ " << A.Height() << " x "
                               << A.Width()  << "\n"
            << "  B[* ,MR] ~ " << B.Height() << " x "
                               << B.Width()  << "\n"
            << "  C[* ,MC] ~ " << C.Height() << " x "
                               << C.Width()  << "\n"
            << "  D[MR,* ] ~ " << D.Height() << " x "
                               << D.Width()  << "\n"
            << "  E[MC,MR] ~ " << E.Height() << " x " << E.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( A.RowAlignment() != E.ColAlignment() ||
        B.RowAlignment() != E.RowAlignment() ||
        C.RowAlignment() != E.ColAlignment() ||
        D.ColAlignment() != E.RowAlignment() )
    {
        std::ostringstream msg;
        msg << "Misaligned LocalTrr2k: \n"
            << "  A[* ,MC] ~ " << A.RowAlignment() << "\n"
            << "  B[* ,MR] ~ " << B.RowAlignment() << "\n"
            << "  C[* ,MC] ~ " << C.RowAlignment() << "\n"
            << "  D[MR,* ] ~ " << D.ColAlignment() << "\n"
            << "  E[MC,MR] ~ " << E.ColAlignment() << " , " <<
                                  E.RowAlignment() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
}

// E := alpha (A^{T/H} B^{T/H} + C D) + beta E
template<typename T>
inline void 
CheckInput
( Orientation orientationOfA,
  Orientation orientationOfB,
  const DistMatrix<T,STAR,MC  >& A, const DistMatrix<T,MR,STAR>& B, 
  const DistMatrix<T,MC,  STAR>& C, const DistMatrix<T,STAR,MR>& D,
  const DistMatrix<T,MC,  MR  >& E )
{
    if( orientationOfA == NORMAL )
        throw std::logic_error("A[* ,MC] must be (Conjugate)Transpose'd");
    if( orientationOfB == NORMAL )
        throw std::logic_error("B[MR,* ] must be (Conjugate)Transpose'd");
    if( A.Grid() != B.Grid() || B.Grid() != C.Grid() ||
        C.Grid() != D.Grid() || D.Grid() != E.Grid() )
        throw std::logic_error
        ("A, B, C, D, and E must be distributed over the same grid");
    if( A.Width()  != E.Height() || B.Height() != E.Width() ||
        C.Height() != E.Height() || D.Width()  != E.Width() ||
        A.Height() != B.Width()  || C.Width()  != D.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal LocalTrr2k: \n"
            << "  A[* ,MC] ~ " << A.Height() << " x "
                               << A.Width()  << "\n"
            << "  B[MR,* ] ~ " << B.Height() << " x "
                               << B.Width()  << "\n"
            << "  C[MC,* ] ~ " << C.Height() << " x "
                               << C.Width()  << "\n"
            << "  D[* ,MR] ~ " << D.Height() << " x "
                               << D.Width()  << "\n"
            << "  E[MC,MR] ~ " << E.Height() << " x " << E.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( A.RowAlignment() != E.ColAlignment() ||
        B.ColAlignment() != E.RowAlignment() ||
        C.ColAlignment() != E.ColAlignment() ||
        D.RowAlignment() != E.RowAlignment() )
    {
        std::ostringstream msg;
        msg << "Misaligned LocalTrr2k: \n"
            << "  A[* ,MC] ~ " << A.RowAlignment() << "\n"
            << "  B[MR,* ] ~ " << B.ColAlignment() << "\n"
            << "  C[MC,* ] ~ " << C.ColAlignment() << "\n"
            << "  D[* ,MR] ~ " << D.RowAlignment() << "\n"
            << "  E[MC,MR] ~ " << E.ColAlignment() << " , " <<
                                  E.RowAlignment() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
}

// E := alpha (A^{T/H} B^{T/H} + C D^{T/H}) + beta E
template<typename T>
inline void 
CheckInput
( Orientation orientationOfA,
  Orientation orientationOfB,
  Orientation orientationOfD,
  const DistMatrix<T,STAR,MC  >& A, const DistMatrix<T,MR,STAR>& B,
  const DistMatrix<T,MC,  STAR>& C, const DistMatrix<T,MR,STAR>& D,
  const DistMatrix<T,MC,  MR  >& E )
{
    if( orientationOfA == NORMAL )
        throw std::logic_error("A[* ,MC] must be (Conjugate)Transpose'd");
    if( orientationOfB == NORMAL )
        throw std::logic_error("B[MR,* ] must be (Conjugate)Transpose'd");
    if( orientationOfD == NORMAL )
        throw std::logic_error("D[MR,* ] must be (Conjugate)Transpose'd");
    if( A.Grid() != B.Grid() || B.Grid() != C.Grid() ||
        C.Grid() != D.Grid() || D.Grid() != E.Grid() )
        throw std::logic_error
        ("A, B, C, D, and E must be distributed over the same grid");
    if( A.Width()  != E.Height() || B.Height() != E.Width() ||
        C.Height() != E.Height() || D.Height() != E.Width() ||
        A.Height() != B.Width()  || C.Width()  != D.Width() )
    {
        std::ostringstream msg;
        msg << "Nonconformal LocalTrr2k: \n"
            << "  A[* ,MC] ~ " << A.Height() << " x "
                               << A.Width()  << "\n"
            << "  B[MR,* ] ~ " << B.Height() << " x "
                               << B.Width()  << "\n"
            << "  C[MC,* ] ~ " << C.Height() << " x "
                               << C.Width()  << "\n"
            << "  D[MR,* ] ~ " << D.Height() << " x "
                               << D.Width()  << "\n"
            << "  E[MC,MR] ~ " << E.Height() << " x " << E.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( A.RowAlignment() != E.ColAlignment() ||
        B.ColAlignment() != E.RowAlignment() ||
        C.ColAlignment() != E.ColAlignment() ||
        D.ColAlignment() != E.RowAlignment() )
    {
        std::ostringstream msg;
        msg << "Misaligned LocalTrr2k: \n"
            << "  A[* ,MC] ~ " << A.RowAlignment() << "\n"
            << "  B[MR,* ] ~ " << B.ColAlignment() << "\n"
            << "  C[MC,* ] ~ " << C.ColAlignment() << "\n"
            << "  D[MR,* ] ~ " << D.ColAlignment() << "\n"
            << "  E[MC,MR] ~ " << E.ColAlignment() << " , " <<
                                  E.RowAlignment() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
}

// E := alpha (A^{T/H} B^{T/H} + C^{T/H} D) + beta E
template<typename T>
inline void 
CheckInput
( Orientation orientationOfA,
  Orientation orientationOfB,
  Orientation orientationOfC,
  const DistMatrix<T,STAR,MC>& A, const DistMatrix<T,MR,STAR>& B,
  const DistMatrix<T,STAR,MC>& C, const DistMatrix<T,STAR,MR>& D,
  const DistMatrix<T,MC,  MR>& E )
{
    if( orientationOfA == NORMAL )
        throw std::logic_error("A[* ,MC] must be (Conjugate)Transpose'd");
    if( orientationOfB == NORMAL )
        throw std::logic_error("B[MR,* ] must be (Conjugate)Transpose'd");
    if( orientationOfC == NORMAL )
        throw std::logic_error("C[* ,MC] must be (Conjugate)Transpose'd");
    if( A.Grid() != B.Grid() || B.Grid() != C.Grid() ||
        C.Grid() != D.Grid() || D.Grid() != E.Grid() )
        throw std::logic_error
        ("A, B, C, D, and E must be distributed over the same grid");
    if( A.Width() != E.Height() || B.Height() != E.Width() ||
        C.Width() != E.Height() || D.Width()  != E.Width() ||
        A.Height() != B.Width() || C.Height() != D.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal LocalTrr2k: \n"
            << "  A[* ,MC] ~ " << A.Height() << " x "
                               << A.Width()  << "\n"
            << "  B[MR,* ] ~ " << B.Height() << " x "
                               << B.Width()  << "\n"
            << "  C[* ,MC] ~ " << C.Height() << " x "
                               << C.Width()  << "\n"
            << "  D[* ,MR] ~ " << D.Height() << " x "
                               << D.Width()  << "\n"
            << "  E[MC,MR] ~ " << E.Height() << " x " << E.Width() << "\n";
        throw std::logic_error( msg.str() );
    }
    if( A.RowAlignment() != E.ColAlignment() ||
        B.ColAlignment() != E.RowAlignment() ||
        C.RowAlignment() != E.ColAlignment() ||
        D.RowAlignment() != E.RowAlignment() )
    {
        std::ostringstream msg;
        msg << "Misaligned LocalTrr2k: \n"
            << "  A[* ,MC] ~ " << A.RowAlignment() << "\n"
            << "  B[MR,* ] ~ " << B.ColAlignment() << "\n"
            << "  C[* ,MC] ~ " << C.RowAlignment() << "\n"
            << "  D[* ,MR] ~ " << D.RowAlignment() << "\n"
            << "  E[MC,MR] ~ " << E.ColAlignment() << " , " <<
                                  E.RowAlignment() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
}

// E := alpha (A^{T/H} B^{T/H} + C^{T/H} D^{T/H}) + beta E
template<typename T>
inline void 
CheckInput
( Orientation orientationOfA, Orientation orientationOfB,
  Orientation orientationOfC, Orientation orientationOfD,
  const DistMatrix<T,STAR,MC>& A, const DistMatrix<T,MR,STAR>& B, 
  const DistMatrix<T,STAR,MC>& C, const DistMatrix<T,MR,STAR>& D,
  const DistMatrix<T,MC,  MR>& E )
{
    if( orientationOfA == NORMAL )
        throw std::logic_error("A[* ,MC] must be (Conjugate)Transpose'd");
    if( orientationOfB == NORMAL )
        throw std::logic_error("B[MR,* ] must be (Conjugate)Transpose'd");
    if( orientationOfC == NORMAL )
        throw std::logic_error("C[* ,MC] must be (Conjugate)Transpose'd");
    if( orientationOfD == NORMAL )
        throw std::logic_error("D[MR,* ] must be (Conjugate)Transpose'd");
    if( A.Grid() != B.Grid() || B.Grid() != C.Grid() ||
        C.Grid() != D.Grid() || D.Grid() != E.Grid() )
        throw std::logic_error
        ("A, B, C, D, and E must be distributed over the same grid");
    if( A.Width()  != E.Height() || B.Height() != E.Width() ||
        C.Width()  != E.Height() || D.Height() != E.Width() ||
        A.Height() != B.Width()  || C.Height() != D.Width() )
    {
        std::ostringstream msg;
        msg << "Nonconformal LocalTrr2k: \n"
            << "  A[* ,MC] ~ " << A.Height() << " x "
                               << A.Width()  << "\n"
            << "  B[MR,* ] ~ " << B.Height() << " x "
                               << B.Width()  << "\n"
            << "  C[* ,MC] ~ " << C.Height() << " x "
                               << C.Width()  << "\n"
            << "  D[MR,* ] ~ " << D.Height() << " x "
                               << D.Width()  << "\n"
            << "  E[MC,MR] ~ " << E.Height() << " x " << E.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( A.RowAlignment() != E.ColAlignment() ||
        B.ColAlignment() != E.RowAlignment() ||
        C.RowAlignment() != E.ColAlignment() ||
        D.ColAlignment() != E.RowAlignment() )
    {
        std::ostringstream msg;
        msg << "Misaligned LocalTrr2k: \n"
            << "  A[* ,MC] ~ " << A.RowAlignment() << "\n"
            << "  B[MR,* ] ~ " << B.ColAlignment() << "\n"
            << "  C[* ,MC] ~ " << C.RowAlignment() << "\n"
            << "  D[MR,* ] ~ " << D.ColAlignment() << "\n"
            << "  E[MC,MR] ~ " << E.ColAlignment() << " , " <<
                                  E.RowAlignment() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
}
#endif // ifndef RELEASE

// E := alpha (A B + C D) + beta E
template<typename T>
inline void
LocalTrr2kKernel
( UpperOrLower uplo,
  T alpha, const DistMatrix<T,MC,  STAR>& A, const DistMatrix<T,STAR,MR>& B,
           const DistMatrix<T,MC,  STAR>& C, const DistMatrix<T,STAR,MR>& D,
  T beta,        DistMatrix<T,MC,  MR  >& E )
{
#ifndef RELEASE
    PushCallStack("LocalTrr2kKernel");
    CheckInput( A, B, C, D, E );
#endif
    const Grid& g = E.Grid();

    DistMatrix<T,MC,STAR> AT(g),  CT(g),
                          AB(g),  CB(g);
    DistMatrix<T,STAR,MR> BL(g), BR(g),
                          DL(g), DR(g);
    DistMatrix<T,MC,MR> ETL(g), ETR(g),
                        EBL(g), EBR(g);
    DistMatrix<T,MC,MR> FTL(g), FBR(g);

    const unsigned half = E.Height()/2;
    E.ScaleTrapezoid( beta, LEFT, uplo );
    LockedPartitionDown
    ( A, AT,
         AB, half );
    LockedPartitionRight( B, BL, BR, half );
    LockedPartitionDown
    ( C, CT,
         CB, half );
    LockedPartitionRight( D, DL, DR, half );
    PartitionDownDiagonal
    ( E, ETL, ETR,
         EBL, EBR, half );

    FTL.AlignWith( ETL );
    FBR.AlignWith( EBR );
    FTL.ResizeTo( ETL.Height(), ETL.Width() );
    FBR.ResizeTo( EBR.Height(), EBR.Width() );
    //------------------------------------------------------------------------//
    if( uplo == LOWER )
    {
        basic::internal::LocalGemm
        ( NORMAL, NORMAL, alpha, AB, BL, (T)1, EBL );
        basic::internal::LocalGemm
        ( NORMAL, NORMAL, alpha, CB, DL, (T)1, EBL );
    }
    else
    {
        basic::internal::LocalGemm
        ( NORMAL, NORMAL, alpha, AT, BR, (T)1, ETR );
        basic::internal::LocalGemm
        ( NORMAL, NORMAL, alpha, CT, DR, (T)1, ETR );
    }

    basic::internal::LocalGemm
    ( NORMAL, NORMAL, alpha, AT, BL, (T)0, FTL );
    basic::internal::LocalGemm
    ( NORMAL, NORMAL, alpha, CT, DL, (T)1, FTL );

    FTL.MakeTrapezoidal( LEFT, uplo );
    basic::Axpy( (T)1, FTL, ETL );

    basic::internal::LocalGemm
    ( NORMAL, NORMAL, alpha, AB, BR, (T)0, FBR );
    basic::internal::LocalGemm
    ( NORMAL, NORMAL, alpha, CB, DR, (T)1, FBR );

    FBR.MakeTrapezoidal( LEFT, uplo );
    basic::Axpy( (T)1, FBR, EBR );
    //------------------------------------------------------------------------//
#ifndef RELEASE
    PopCallStack();
#endif
}

// E := alpha (A B + C D^{T/H}) + beta C
template<typename T>
inline void
LocalTrr2kKernel
( UpperOrLower uplo,
  Orientation orientationOfD,
  T alpha, const DistMatrix<T,MC,  STAR>& A, const DistMatrix<T,STAR,MR>& B,
           const DistMatrix<T,MC,  STAR>& C, const DistMatrix<T,MR,STAR>& D,
  T beta,        DistMatrix<T,MC,  MR  >& E )
{
#ifndef RELEASE
    PushCallStack("LocalTrr2kKernel");
    CheckInput( orientationOfD, A, B, C, D, E );
#endif
    const Grid& g = E.Grid();

    DistMatrix<T,MC,STAR> AT(g),  CT(g),
                          AB(g),  CB(g);
    DistMatrix<T,MR,STAR> DT(g), 
                          DB(g);
    DistMatrix<T,STAR,MR> BL(g), BR(g);
    DistMatrix<T,MC,MR> ETL(g), ETR(g),
                        EBL(g), EBR(g);
    DistMatrix<T,MC,MR> FTL(g), FBR(g);

    const unsigned half = E.Height()/2;
    E.ScaleTrapezoid( beta, LEFT, uplo );
    LockedPartitionDown
    ( A, AT,
         AB, half );
    LockedPartitionRight( B, BL, BR, half );
    LockedPartitionDown
    ( C, CT,
         CB, half );
    LockedPartitionDown
    ( D, DT,
         DB, half );
    PartitionDownDiagonal
    ( E, ETL, ETR,
         EBL, EBR, half );

    FTL.AlignWith( ETL );
    FBR.AlignWith( EBR );
    FTL.ResizeTo( ETL.Height(), ETL.Width() );
    FBR.ResizeTo( EBR.Height(), EBR.Width() );
    //------------------------------------------------------------------------//
    if( uplo == LOWER )
    {
        basic::internal::LocalGemm
        ( NORMAL, NORMAL, alpha, AB, BL, (T)1, EBL );
        basic::internal::LocalGemm
        ( NORMAL, orientationOfD, alpha, CB, DT, (T)1, EBL );
    }
    else
    {
        basic::internal::LocalGemm
        ( NORMAL, NORMAL, alpha, AT, BR, (T)1, ETR );
        basic::internal::LocalGemm
        ( NORMAL, orientationOfD, alpha, CT, DB, (T)1, ETR );
    }

    basic::internal::LocalGemm
    ( NORMAL, NORMAL, alpha, AT, BL, (T)0, FTL );
    basic::internal::LocalGemm
    ( NORMAL, orientationOfD, alpha, CT, DT, (T)1, FTL );

    FTL.MakeTrapezoidal( LEFT, uplo );
    basic::Axpy( (T)1, FTL, ETL );

    basic::internal::LocalGemm
    ( NORMAL, NORMAL, alpha, AB, BR, (T)0, FBR );
    basic::internal::LocalGemm
    ( NORMAL, orientationOfD, alpha, CB, DB, (T)1, FBR );

    FBR.MakeTrapezoidal( LEFT, uplo );
    basic::Axpy( (T)1, FBR, EBR );
    //------------------------------------------------------------------------//
#ifndef RELEASE
    PopCallStack();
#endif
}

// E := alpha (A B + C^{T/H} D) + beta E
template<typename T>
inline void
LocalTrr2kKernel
( UpperOrLower uplo,
  Orientation orientationOfC,
  T alpha, const DistMatrix<T,MC,  STAR>& A, const DistMatrix<T,STAR,MR>& B,
           const DistMatrix<T,STAR,MC  >& C, const DistMatrix<T,STAR,MR>& D,
  T beta,        DistMatrix<T,MC,  MR  >& E )
{
#ifndef RELEASE
    PushCallStack("LocalTrr2kKernel");
    CheckInput( orientationOfC, A, B, C, D, E );
#endif
    const Grid& g = E.Grid();

    DistMatrix<T,MC,STAR> AT(g),
                          AB(g);
    DistMatrix<T,STAR,MC> CL(g), CR(g);
    DistMatrix<T,STAR,MR> BL(g), BR(g),
                          DL(g), DR(g);
    DistMatrix<T,MC,MR> ETL(g), ETR(g),
                        EBL(g), EBR(g);
    DistMatrix<T,MC,MR> FTL(g), FBR(g);

    const unsigned half = E.Height()/2;
    E.ScaleTrapezoid( beta, LEFT, uplo );
    LockedPartitionDown
    ( A, AT,
         AB, half );
    LockedPartitionRight( B, BL, BR, half );
    LockedPartitionRight( C, CL, CR, half );
    LockedPartitionRight( D, DL, DR, half );
    PartitionDownDiagonal
    ( E, ETL, ETR,
         EBL, EBR, half );

    FTL.AlignWith( ETL );
    FBR.AlignWith( EBR );
    FTL.ResizeTo( ETL.Height(), ETL.Width() );
    FBR.ResizeTo( EBR.Height(), EBR.Width() );
    //------------------------------------------------------------------------//
    if( uplo == LOWER )
    {
        basic::internal::LocalGemm
        ( NORMAL, NORMAL, alpha, AB, BL, (T)1, EBL );
        basic::internal::LocalGemm
        ( orientationOfC, NORMAL, alpha, CR, DL, (T)1, EBL );
    }
    else
    {
        basic::internal::LocalGemm
        ( NORMAL, NORMAL, alpha, AT, BR, (T)1, ETR );
        basic::internal::LocalGemm
        ( orientationOfC, NORMAL, alpha, CL, DR, (T)1, ETR );
    }

    basic::internal::LocalGemm
    ( NORMAL, NORMAL, alpha, AT, BL, (T)0, FTL );
    basic::internal::LocalGemm
    ( orientationOfC, NORMAL, alpha, CL, DL, (T)1, FTL );

    FTL.MakeTrapezoidal( LEFT, uplo );
    basic::Axpy( (T)1, FTL, ETL );

    basic::internal::LocalGemm
    ( NORMAL, NORMAL, alpha, AB, BR, (T)0, FBR );
    basic::internal::LocalGemm
    ( orientationOfC, NORMAL, alpha, CR, DR, (T)1, FBR );

    FBR.MakeTrapezoidal( LEFT, uplo );
    basic::Axpy( (T)1, FBR, EBR );
    //------------------------------------------------------------------------//
#ifndef RELEASE
    PopCallStack();
#endif
}

// E := alpha (A B + C^{T/H} D^{T/H}) + beta E
template<typename T>
inline void
LocalTrr2kKernel
( UpperOrLower uplo,
  Orientation orientationOfC, Orientation orientationOfD,
  T alpha, const DistMatrix<T,MC,  STAR>& A, const DistMatrix<T,STAR,MR>& B,
           const DistMatrix<T,STAR,MC  >& C, const DistMatrix<T,MR,STAR>& D,
  T beta,        DistMatrix<T,MC,  MR  >& E )
{
#ifndef RELEASE
    PushCallStack("LocalTrr2kKernel");
    CheckInput( orientationOfC, orientationOfD, A, B, C, D, E );
#endif
    const Grid& g = E.Grid();

    DistMatrix<T,MC,STAR> AT(g),
                          AB(g);
    DistMatrix<T,STAR,MR> BL(g), BR(g);
    DistMatrix<T,STAR,MC> CL(g), CR(g);
    DistMatrix<T,MR,STAR> DT(g),
                          DB(g);
    DistMatrix<T,MC,MR> ETL(g), ETR(g),
                        EBL(g), EBR(g);
    DistMatrix<T,MC,MR> FTL(g), FBR(g);

    const unsigned half = E.Height()/2;
    E.ScaleTrapezoid( beta, LEFT, uplo );
    LockedPartitionDown
    ( A, AT,
         AB, half );
    LockedPartitionRight( B, BL, BR, half );
    LockedPartitionRight( C, CL, CR, half );
    LockedPartitionDown
    ( D, DT,
         DB, half );
    PartitionDownDiagonal
    ( E, ETL, ETR,
         EBL, EBR, half );

    FTL.AlignWith( ETL );
    FBR.AlignWith( EBR );
    FTL.ResizeTo( ETL.Height(), ETL.Width() );
    FBR.ResizeTo( EBR.Height(), EBR.Width() );
    //------------------------------------------------------------------------//
    if( uplo == LOWER )
    {
        basic::internal::LocalGemm
        ( NORMAL, NORMAL, alpha, AB, BL, (T)1, EBL );
        basic::internal::LocalGemm
        ( orientationOfC, orientationOfD, alpha, CR, DT, (T)1, EBL );
    }
    else
    {
        basic::internal::LocalGemm
        ( NORMAL, NORMAL, alpha, AT, BR, (T)1, ETR );
        basic::internal::LocalGemm
        ( orientationOfC, orientationOfD, alpha, CL, DB, (T)1, ETR );
    }

    basic::internal::LocalGemm
    ( NORMAL, NORMAL, alpha, AT, BL, (T)0, FTL );
    basic::internal::LocalGemm
    ( orientationOfC, orientationOfD, alpha, CL, DT, (T)1, FTL );

    FTL.MakeTrapezoidal( LEFT, uplo );
    basic::Axpy( (T)1, FTL, ETL );

    basic::internal::LocalGemm
    ( NORMAL, NORMAL, alpha, AB, BR, (T)0, FBR );
    basic::internal::LocalGemm
    ( orientationOfC, orientationOfD, alpha, CR, DB, (T)1, FBR );

    FBR.MakeTrapezoidal( LEFT, uplo );
    basic::Axpy( (T)1, FBR, EBR );
    //------------------------------------------------------------------------//
#ifndef RELEASE
    PopCallStack();
#endif
}

// E := alpha (A B^{T/H} + C D) + beta C
template<typename T>
inline void
LocalTrr2kKernel
( UpperOrLower uplo,
  Orientation orientationOfB,
  T alpha, const DistMatrix<T,MC,  STAR>& A, const DistMatrix<T,STAR,MR>& B,
           const DistMatrix<T,MC,  STAR>& C, const DistMatrix<T,STAR,MR>& D,
  T beta,        DistMatrix<T,MC,  MR  >& E )
{
#ifndef RELEASE
    PushCallStack("LocalTrr2kKernel");
    CheckInput( orientationOfB, A, B, C, D, E );
#endif
    const Grid& g = E.Grid();

    DistMatrix<T,MC,STAR> AT(g),  CT(g),
                          AB(g),  CB(g);
    DistMatrix<T,MR,STAR> BT(g), 
                          BB(g);
    DistMatrix<T,STAR,MR> DL(g), DR(g);
    DistMatrix<T,MC,MR> ETL(g), ETR(g),
                        EBL(g), EBR(g);
    DistMatrix<T,MC,MR> FTL(g), FBR(g);

    const unsigned half = E.Height()/2;
    E.ScaleTrapezoid( beta, LEFT, uplo );
    LockedPartitionDown
    ( A, AT,
         AB, half );
    LockedPartitionDown
    ( B, BT,
         BB, half );
    LockedPartitionDown
    ( C, CT,
         CB, half );
    LockedPartitionRight( D, DL, DR, half );
    PartitionDownDiagonal
    ( E, ETL, ETR,
         EBL, EBR, half );

    FTL.AlignWith( ETL );
    FBR.AlignWith( EBR );
    FTL.ResizeTo( ETL.Height(), ETL.Width() );
    FBR.ResizeTo( EBR.Height(), EBR.Width() );
    //------------------------------------------------------------------------//
    if( uplo == LOWER )
    {
        basic::internal::LocalGemm
        ( NORMAL, orientationOfB, alpha, AB, BT, (T)1, EBL );
        basic::internal::LocalGemm
        ( NORMAL, NORMAL, alpha, CB, DL, (T)1, EBL );
    }
    else
    {
        basic::internal::LocalGemm
        ( NORMAL, orientationOfB, alpha, AT, BB, (T)1, ETR );
        basic::internal::LocalGemm
        ( NORMAL, NORMAL, alpha, CT, DR, (T)1, ETR );
    }

    basic::internal::LocalGemm
    ( NORMAL, orientationOfB, alpha, AT, BT, (T)0, FTL );
    basic::internal::LocalGemm
    ( NORMAL, NORMAL, alpha, CT, DL, (T)1, FTL );

    FTL.MakeTrapezoidal( LEFT, uplo );
    basic::Axpy( (T)1, FTL, ETL );

    basic::internal::LocalGemm
    ( NORMAL, orientationOfB, alpha, AB, BB, (T)0, FBR );
    basic::internal::LocalGemm
    ( NORMAL, NORMAL, alpha, CB, DR, (T)1, FBR );

    FBR.MakeTrapezoidal( LEFT, uplo );
    basic::Axpy( (T)1, FBR, EBR );
    //------------------------------------------------------------------------//
#ifndef RELEASE
    PopCallStack();
#endif
}

// E := alpha (A B^{T/H} + C D^{T/H}) + beta C
template<typename T>
inline void
LocalTrr2kKernel
( UpperOrLower uplo,
  Orientation orientationOfB,
  Orientation orientationOfD,
  T alpha, const DistMatrix<T,MC,STAR>& A, const DistMatrix<T,MR,STAR>& B,
           const DistMatrix<T,MC,STAR>& C, const DistMatrix<T,MR,STAR>& D,
  T beta,        DistMatrix<T,MC,MR  >& E )
{
#ifndef RELEASE
    PushCallStack("LocalTrr2kKernel");
    CheckInput( orientationOfB, orientationOfD, A, B, C, D, E );
#endif
    const Grid& g = E.Grid();

    DistMatrix<T,MC,STAR> AT(g),  CT(g),
                          AB(g),  CB(g);
    DistMatrix<T,MR,STAR> BT(g),  DT(g),
                          BB(g),  DB(g);
    DistMatrix<T,MC,MR> ETL(g), ETR(g),
                        EBL(g), EBR(g);
    DistMatrix<T,MC,MR> FTL(g), FBR(g);

    const unsigned half = E.Height()/2;
    E.ScaleTrapezoid( beta, LEFT, uplo );
    LockedPartitionDown
    ( A, AT,
         AB, half );
    LockedPartitionDown
    ( B, BT,
         BB, half );
    LockedPartitionDown
    ( C, CT,
         CB, half );
    LockedPartitionDown
    ( D, DT,
         DB, half );
    PartitionDownDiagonal
    ( E, ETL, ETR,
         EBL, EBR, half );

    FTL.AlignWith( ETL );
    FBR.AlignWith( EBR );
    FTL.ResizeTo( ETL.Height(), ETL.Width() );
    FBR.ResizeTo( EBR.Height(), EBR.Width() );
    //------------------------------------------------------------------------//
    if( uplo == LOWER )
    {
        basic::internal::LocalGemm
        ( NORMAL, orientationOfB, alpha, AB, BT, (T)1, EBL );
        basic::internal::LocalGemm
        ( NORMAL, orientationOfD, alpha, CB, DT, (T)1, EBL );
    }
    else
    {
        basic::internal::LocalGemm
        ( NORMAL, orientationOfB, alpha, AT, BB, (T)1, ETR );
        basic::internal::LocalGemm
        ( NORMAL, orientationOfD, alpha, CT, DB, (T)1, ETR );
    }

    basic::internal::LocalGemm
    ( NORMAL, orientationOfB, alpha, AT, BT, (T)0, FTL );
    basic::internal::LocalGemm
    ( NORMAL, orientationOfD, alpha, CT, DT, (T)1, FTL );

    FTL.MakeTrapezoidal( LEFT, uplo );
    basic::Axpy( (T)1, FTL, ETL );

    basic::internal::LocalGemm
    ( NORMAL, orientationOfB, alpha, AB, BB, (T)0, FBR );
    basic::internal::LocalGemm
    ( NORMAL, orientationOfD, alpha, CB, DB, (T)1, FBR );

    FBR.MakeTrapezoidal( LEFT, uplo );
    basic::Axpy( (T)1, FBR, EBR );
    //------------------------------------------------------------------------//
#ifndef RELEASE
    PopCallStack();
#endif
}

// E := alpha (A B^{T/H} + C^{T/H} D) + beta E
template<typename T>
inline void
LocalTrr2kKernel
( UpperOrLower uplo,
  Orientation orientationOfB, Orientation orientationOfC,
  T alpha, const DistMatrix<T,MC,  STAR>& A, const DistMatrix<T,MR,STAR>& B,
           const DistMatrix<T,STAR,MC  >& C, const DistMatrix<T,STAR,MR>& D,
  T beta,        DistMatrix<T,MC,  MR  >& E )
{
#ifndef RELEASE
    PushCallStack("LocalTrr2kKernel");
    CheckInput( orientationOfB, orientationOfC, A, B, C, D, E );
#endif
    const Grid& g = E.Grid();

    DistMatrix<T,MC,STAR> AT(g),
                          AB(g);
    DistMatrix<T,MR,STAR> BT(g),
                          BB(g);
    DistMatrix<T,STAR,MC> CL(g), CR(g);
    DistMatrix<T,STAR,MR> DL(g), DR(g);
    DistMatrix<T,MC,MR> ETL(g), ETR(g),
                        EBL(g), EBR(g);
    DistMatrix<T,MC,MR> FTL(g), FBR(g);

    const unsigned half = E.Height()/2;
    E.ScaleTrapezoid( beta, LEFT, uplo );
    LockedPartitionDown
    ( A, AT,
         AB, half );
    LockedPartitionDown
    ( B, BT,
         BB, half );
    LockedPartitionRight( C, CL, CR, half );
    LockedPartitionRight( D, DL, DR, half );
    PartitionDownDiagonal
    ( E, ETL, ETR,
         EBL, EBR, half );

    FTL.AlignWith( ETL );
    FBR.AlignWith( EBR );
    FTL.ResizeTo( ETL.Height(), ETL.Width() );
    FBR.ResizeTo( EBR.Height(), EBR.Width() );
    //------------------------------------------------------------------------//
    if( uplo == LOWER )
    {
        basic::internal::LocalGemm
        ( NORMAL, orientationOfB, alpha, AB, BT, (T)1, EBL );
        basic::internal::LocalGemm
        ( orientationOfC, NORMAL, alpha, CR, DL, (T)1, EBL );
    }
    else
    {
        basic::internal::LocalGemm
        ( NORMAL, orientationOfB, alpha, AT, BB, (T)1, ETR );
        basic::internal::LocalGemm
        ( orientationOfC, NORMAL, alpha, CL, DR, (T)1, ETR );
    }

    basic::internal::LocalGemm
    ( NORMAL, orientationOfB, alpha, AT, BT, (T)0, FTL );
    basic::internal::LocalGemm
    ( orientationOfC, NORMAL, alpha, CL, DL, (T)1, FTL );

    FTL.MakeTrapezoidal( LEFT, uplo );
    basic::Axpy( (T)1, FTL, ETL );

    basic::internal::LocalGemm
    ( NORMAL, orientationOfB, alpha, AB, BB, (T)0, FBR );
    basic::internal::LocalGemm
    ( orientationOfC, NORMAL, alpha, CR, DR, (T)1, FBR );

    FBR.MakeTrapezoidal( LEFT, uplo );
    basic::Axpy( (T)1, FBR, EBR );
    //------------------------------------------------------------------------//
#ifndef RELEASE
    PopCallStack();
#endif
}

// E := alpha (A B^{T/H} + C^{T/H} D^{T/H}) + beta C
template<typename T>
inline void
LocalTrr2kKernel
( UpperOrLower uplo,
  Orientation orientationOfB,
  Orientation orientationOfC,
  Orientation orientationOfD,
  T alpha, const DistMatrix<T,MC,  STAR>& A, const DistMatrix<T,MR,STAR>& B,
           const DistMatrix<T,STAR,MC  >& C, const DistMatrix<T,MR,STAR>& D,
  T beta,        DistMatrix<T,MC,  MR  >& E )
{
#ifndef RELEASE
    PushCallStack("LocalTrr2kKernel");
    CheckInput( orientationOfB, orientationOfC, orientationOfD, A, B, C, D, E );
#endif
    const Grid& g = E.Grid();

    DistMatrix<T,MC,STAR> AT(g),
                          AB(g);
    DistMatrix<T,MR,STAR> BT(g), 
                          BB(g);
    DistMatrix<T,STAR,MC> CL(g), CR(g);
    DistMatrix<T,MR,STAR> DT(g), 
                          DB(g);
    DistMatrix<T,MC,MR> ETL(g), ETR(g),
                        EBL(g), EBR(g);
    DistMatrix<T,MC,MR> FTL(g), FBR(g);

    const unsigned half = E.Height()/2;
    E.ScaleTrapezoid( beta, LEFT, uplo );
    LockedPartitionDown
    ( A, AT,
         AB, half );
    LockedPartitionDown
    ( B, BT,
         BB, half );
    LockedPartitionRight( C, CL, CR, half );
    LockedPartitionDown
    ( D, DT,
         DB, half );
    PartitionDownDiagonal
    ( E, ETL, ETR,
         EBL, EBR, half );

    FTL.AlignWith( ETL );
    FBR.AlignWith( EBR );
    FTL.ResizeTo( ETL.Height(), ETL.Width() );
    FBR.ResizeTo( EBR.Height(), EBR.Width() );
    //------------------------------------------------------------------------//
    if( uplo == LOWER )
    {
        basic::internal::LocalGemm
        ( NORMAL, orientationOfB, alpha, AB, BT, (T)1, EBL );
        basic::internal::LocalGemm
        ( orientationOfC, orientationOfD, alpha, CR, DT, (T)1, EBL );
    }
    else
    {
        basic::internal::LocalGemm
        ( NORMAL, orientationOfB, alpha, AT, BB, (T)1, ETR );
        basic::internal::LocalGemm
        ( orientationOfC, orientationOfD, alpha, CL, DB, (T)1, ETR );
    }

    basic::internal::LocalGemm
    ( NORMAL, orientationOfB, alpha, AT, BT, (T)0, FTL );
    basic::internal::LocalGemm
    ( orientationOfC, orientationOfD, alpha, CL, DT, (T)1, FTL );

    FTL.MakeTrapezoidal( LEFT, uplo );
    basic::Axpy( (T)1, FTL, ETL );

    basic::internal::LocalGemm
    ( NORMAL, orientationOfB, alpha, AB, BB, (T)0, FBR );
    basic::internal::LocalGemm
    ( orientationOfC, orientationOfD, alpha, CR, DB, (T)1, FBR );

    FBR.MakeTrapezoidal( LEFT, uplo );
    basic::Axpy( (T)1, FBR, EBR );
    //------------------------------------------------------------------------//
#ifndef RELEASE
    PopCallStack();
#endif
}

// E := alpha (A^{T/H} B + C D) + beta E
template<typename T>
inline void
LocalTrr2kKernel
( UpperOrLower uplo,
  Orientation orientationOfA,
  T alpha, const DistMatrix<T,STAR,MC  >& A, const DistMatrix<T,STAR,MR>& B,
           const DistMatrix<T,MC,  STAR>& C, const DistMatrix<T,STAR,MR>& D,
  T beta,        DistMatrix<T,MC,  MR  >& E )
{
#ifndef RELEASE
    PushCallStack("LocalTrr2kKernel");
    CheckInput( orientationOfA, A, B, C, D, E );
#endif
    const Grid& g = E.Grid();

    DistMatrix<T,STAR,MC> AL(g), AR(g);
    DistMatrix<T,MC,STAR> CT(g),
                          CB(g);
    DistMatrix<T,STAR,MR> BL(g), BR(g),
                          DL(g), DR(g);
    DistMatrix<T,MC,MR> ETL(g), ETR(g),
                        EBL(g), EBR(g);
    DistMatrix<T,MC,MR> FTL(g), FBR(g);

    const unsigned half = E.Height()/2;
    E.ScaleTrapezoid( beta, LEFT, uplo );
    LockedPartitionRight( A, AL, AR, half );
    LockedPartitionRight( B, BL, BR, half );
    LockedPartitionDown
    ( C, CT,
         CB, half );
    LockedPartitionRight( D, DL, DR, half );
    PartitionDownDiagonal
    ( E, ETL, ETR,
         EBL, EBR, half );

    FTL.AlignWith( ETL );
    FBR.AlignWith( EBR );
    FTL.ResizeTo( ETL.Height(), ETL.Width() );
    FBR.ResizeTo( EBR.Height(), EBR.Width() );
    //------------------------------------------------------------------------//
    if( uplo == LOWER )
    {
        basic::internal::LocalGemm
        ( orientationOfA, NORMAL, alpha, AR, BL, (T)1, EBL );
        basic::internal::LocalGemm
        ( NORMAL, NORMAL, alpha, CB, DL, (T)1, EBL );
    }
    else
    {
        basic::internal::LocalGemm
        ( orientationOfA, NORMAL, alpha, AL, BR, (T)1, ETR );
        basic::internal::LocalGemm
        ( NORMAL, NORMAL, alpha, CT, DR, (T)1, ETR );
    }

    basic::internal::LocalGemm
    ( orientationOfA, NORMAL, alpha, AL, BL, (T)0, FTL );
    basic::internal::LocalGemm
    ( NORMAL, NORMAL, alpha, CT, DL, (T)1, FTL );

    FTL.MakeTrapezoidal( LEFT, uplo );
    basic::Axpy( (T)1, FTL, ETL );

    basic::internal::LocalGemm
    ( orientationOfA, NORMAL, alpha, AR, BR, (T)0, FBR );
    basic::internal::LocalGemm
    ( NORMAL, NORMAL, alpha, CB, DR, (T)1, FBR );

    FBR.MakeTrapezoidal( LEFT, uplo );
    basic::Axpy( (T)1, FBR, EBR );
    //------------------------------------------------------------------------//
#ifndef RELEASE
    PopCallStack();
#endif
}

// E := alpha (A^{T/H} B + C D^{T/H}) + beta E
template<typename T>
inline void
LocalTrr2kKernel
( UpperOrLower uplo,
  Orientation orientationOfA, Orientation orientationOfD,
  T alpha, const DistMatrix<T,STAR,MC  >& A, const DistMatrix<T,STAR,MR>& B,
           const DistMatrix<T,MC,  STAR>& C, const DistMatrix<T,MR,STAR>& D,
  T beta,        DistMatrix<T,MC,  MR  >& E )
{
#ifndef RELEASE
    PushCallStack("LocalTrr2kKernel");
    CheckInput( orientationOfA, orientationOfD, A, B, C, D, E );
#endif
    const Grid& g = E.Grid();

    DistMatrix<T,STAR,MC> AL(g), AR(g);
    DistMatrix<T,STAR,MR> BL(g), BR(g);
    DistMatrix<T,MC,STAR> CT(g),
                          CB(g);
    DistMatrix<T,MR,STAR> DT(g),
                          DB(g);
    DistMatrix<T,MC,MR> ETL(g), ETR(g),
                        EBL(g), EBR(g);
    DistMatrix<T,MC,MR> FTL(g), FBR(g);

    const unsigned half = E.Height()/2;
    E.ScaleTrapezoid( beta, LEFT, uplo );
    LockedPartitionRight( A, AL, AR, half );
    LockedPartitionRight( B, BL, BR, half );
    LockedPartitionDown
    ( C, CT,
         CB, half );
    LockedPartitionDown
    ( D, DT,
         DB, half );
    PartitionDownDiagonal
    ( E, ETL, ETR,
         EBL, EBR, half );

    FTL.AlignWith( ETL );
    FBR.AlignWith( EBR );
    FTL.ResizeTo( ETL.Height(), ETL.Width() );
    FBR.ResizeTo( EBR.Height(), EBR.Width() );
    //------------------------------------------------------------------------//
    if( uplo == LOWER )
    {
        basic::internal::LocalGemm
        ( orientationOfA, NORMAL, alpha, AR, BL, (T)1, EBL );
        basic::internal::LocalGemm
        ( NORMAL, orientationOfD, alpha, CB, DT, (T)1, EBL );
    }
    else
    {
        basic::internal::LocalGemm
        ( orientationOfA, NORMAL, alpha, AL, BR, (T)1, ETR );
        basic::internal::LocalGemm
        ( NORMAL, orientationOfD, alpha, CT, DB, (T)1, ETR );
    }

    basic::internal::LocalGemm
    ( orientationOfA, NORMAL, alpha, AL, BL, (T)0, FTL );
    basic::internal::LocalGemm
    ( NORMAL, orientationOfD, alpha, CT, DT, (T)1, FTL );

    FTL.MakeTrapezoidal( LEFT, uplo );
    basic::Axpy( (T)1, FTL, ETL );

    basic::internal::LocalGemm
    ( orientationOfA, NORMAL, alpha, AR, BR, (T)0, FBR );
    basic::internal::LocalGemm
    ( NORMAL, orientationOfD, alpha, CB, DB, (T)1, FBR );

    FBR.MakeTrapezoidal( LEFT, uplo );
    basic::Axpy( (T)1, FBR, EBR );
    //------------------------------------------------------------------------//
#ifndef RELEASE
    PopCallStack();
#endif
}

// E := alpha (A^{T/H} B + C^{T/H} D) + beta E
template<typename T>
inline void
LocalTrr2kKernel
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfC,
  T alpha, const DistMatrix<T,STAR,MC>& A, const DistMatrix<T,STAR,MR>& B,
           const DistMatrix<T,STAR,MC>& C, const DistMatrix<T,STAR,MR>& D,
  T beta,        DistMatrix<T,MC,  MR>& E )
{
#ifndef RELEASE
    PushCallStack("LocalTrr2kKernel");
    CheckInput( orientationOfA, orientationOfC, A, B, C, D, E );
#endif
    const Grid& g = E.Grid();

    DistMatrix<T,STAR,MC> AL(g), AR(g),
                          CL(g), CR(g);
    DistMatrix<T,STAR,MR> BL(g), BR(g),
                          DL(g), DR(g);
    DistMatrix<T,MC,MR> ETL(g), ETR(g),
                        EBL(g), EBR(g);
    DistMatrix<T,MC,MR> FTL(g), FBR(g);

    const unsigned half = E.Height()/2;
    E.ScaleTrapezoid( beta, LEFT, uplo );
    LockedPartitionRight( A, AL, AR, half );
    LockedPartitionRight( B, BL, BR, half );
    LockedPartitionRight( C, CL, CR, half );
    LockedPartitionRight( D, DL, DR, half );
    PartitionDownDiagonal
    ( E, ETL, ETR,
         EBL, EBR, half );

    FTL.AlignWith( ETL );
    FBR.AlignWith( EBR );
    FTL.ResizeTo( ETL.Height(), ETL.Width() );
    FBR.ResizeTo( EBR.Height(), EBR.Width() );
    //------------------------------------------------------------------------//
    if( uplo == LOWER )
    {
        basic::internal::LocalGemm
        ( orientationOfA, NORMAL, alpha, AR, BL, (T)1, EBL );
        basic::internal::LocalGemm
        ( orientationOfC, NORMAL, alpha, CR, DL, (T)1, EBL );
    }
    else
    {
        basic::internal::LocalGemm
        ( orientationOfA, NORMAL, alpha, AL, BR, (T)1, ETR );
        basic::internal::LocalGemm
        ( orientationOfC, NORMAL, alpha, CL, DR, (T)1, ETR );
    }

    basic::internal::LocalGemm
    ( orientationOfA, NORMAL, alpha, AL, BL, (T)0, FTL );
    basic::internal::LocalGemm
    ( orientationOfC, NORMAL, alpha, CL, DL, (T)1, FTL );

    FTL.MakeTrapezoidal( LEFT, uplo );
    basic::Axpy( (T)1, FTL, ETL );

    basic::internal::LocalGemm
    ( orientationOfA, NORMAL, alpha, AR, BR, (T)0, FBR );
    basic::internal::LocalGemm
    ( orientationOfC, NORMAL, alpha, CR, DR, (T)1, FBR );

    FBR.MakeTrapezoidal( LEFT, uplo );
    basic::Axpy( (T)1, FBR, EBR );
    //------------------------------------------------------------------------//
#ifndef RELEASE
    PopCallStack();
#endif
}

// E := alpha (A^{T/H} B + C^{T/H} D^{T/H}) + beta E
template<typename T>
inline void
LocalTrr2kKernel
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfC,
  Orientation orientationOfD,
  T alpha, const DistMatrix<T,STAR,MC>& A, const DistMatrix<T,STAR,MR>& B,
           const DistMatrix<T,STAR,MC>& C, const DistMatrix<T,MR,STAR>& D,
  T beta,        DistMatrix<T,MC,  MR>& E )
{
#ifndef RELEASE
    PushCallStack("LocalTrr2kKernel");
    CheckInput( orientationOfA, orientationOfC, orientationOfD, A, B, C, D, E );
#endif
    const Grid& g = E.Grid();

    DistMatrix<T,STAR,MC> AL(g), AR(g),
                          CL(g), CR(g);
    DistMatrix<T,STAR,MR> BL(g), BR(g);
    DistMatrix<T,MR,STAR> DT(g),
                          DB(g);
    DistMatrix<T,MC,MR> ETL(g), ETR(g),
                        EBL(g), EBR(g);
    DistMatrix<T,MC,MR> FTL(g), FBR(g);

    const unsigned half = E.Height()/2;
    E.ScaleTrapezoid( beta, LEFT, uplo );
    LockedPartitionRight( A, AL, AR, half );
    LockedPartitionRight( B, BL, BR, half );
    LockedPartitionRight( C, CL, CR, half );
    LockedPartitionDown
    ( D, DT,
         DB, half );
    PartitionDownDiagonal
    ( E, ETL, ETR,
         EBL, EBR, half );

    FTL.AlignWith( ETL );
    FBR.AlignWith( EBR );
    FTL.ResizeTo( ETL.Height(), ETL.Width() );
    FBR.ResizeTo( EBR.Height(), EBR.Width() );
    //------------------------------------------------------------------------//
    if( uplo == LOWER )
    {
        basic::internal::LocalGemm
        ( orientationOfA, NORMAL, alpha, AR, BL, (T)1, EBL );
        basic::internal::LocalGemm
        ( orientationOfC, orientationOfD, alpha, CR, DT, (T)1, EBL );
    }
    else
    {
        basic::internal::LocalGemm
        ( orientationOfA, NORMAL, alpha, AL, BR, (T)1, ETR );
        basic::internal::LocalGemm
        ( orientationOfC, orientationOfD, alpha, CL, DB, (T)1, ETR );
    }

    basic::internal::LocalGemm
    ( orientationOfA, NORMAL, alpha, AL, BL, (T)0, FTL );
    basic::internal::LocalGemm
    ( orientationOfC, orientationOfD, alpha, CL, DT, (T)1, FTL );

    FTL.MakeTrapezoidal( LEFT, uplo );
    basic::Axpy( (T)1, FTL, ETL );

    basic::internal::LocalGemm
    ( orientationOfA, NORMAL, alpha, AR, BR, (T)0, FBR );
    basic::internal::LocalGemm
    ( orientationOfC, orientationOfD, alpha, CR, DB, (T)1, FBR );

    FBR.MakeTrapezoidal( LEFT, uplo );
    basic::Axpy( (T)1, FBR, EBR );
    //------------------------------------------------------------------------//
#ifndef RELEASE
    PopCallStack();
#endif
}

// E := alpha (A^{T/H} B^{T/H} + C D) + beta E
template<typename T>
inline void
LocalTrr2kKernel
( UpperOrLower uplo,
  Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const DistMatrix<T,STAR,MC  >& A, const DistMatrix<T,MR,STAR>& B,
           const DistMatrix<T,MC,  STAR>& C, const DistMatrix<T,STAR,MR>& D,
  T beta,        DistMatrix<T,MC,  MR  >& E )
{
#ifndef RELEASE
    PushCallStack("LocalTrr2kKernel");
    CheckInput( orientationOfA, orientationOfB, A, B, C, D, E );
#endif
    const Grid& g = E.Grid();

    DistMatrix<T,STAR,MC> AL(g), AR(g);
    DistMatrix<T,MR,STAR> BT(g),
                          BB(g);
    DistMatrix<T,MC,STAR> CT(g),
                          CB(g);
    DistMatrix<T,STAR,MR> DL(g), DR(g);
    DistMatrix<T,MC,MR> ETL(g), ETR(g),
                        EBL(g), EBR(g);
    DistMatrix<T,MC,MR> FTL(g), FBR(g);

    const unsigned half = E.Height()/2;
    E.ScaleTrapezoid( beta, LEFT, uplo );
    LockedPartitionRight( A, AL, AR, half );
    LockedPartitionDown
    ( B, BT,
         BB, half );
    LockedPartitionDown
    ( C, CT,
         CB, half );
    LockedPartitionRight( D, DL, DR, half );
    PartitionDownDiagonal
    ( E, ETL, ETR,
         EBL, EBR, half );

    FTL.AlignWith( ETL );
    FBR.AlignWith( EBR );
    FTL.ResizeTo( ETL.Height(), ETL.Width() );
    FBR.ResizeTo( EBR.Height(), EBR.Width() );
    //------------------------------------------------------------------------//
    if( uplo == LOWER )
    {
        basic::internal::LocalGemm
        ( orientationOfA, orientationOfB, alpha, AR, BT, (T)1, EBL );
        basic::internal::LocalGemm
        ( NORMAL, NORMAL, alpha, CB, DL, (T)1, EBL );
    }
    else
    {
        basic::internal::LocalGemm
        ( orientationOfA, orientationOfB, alpha, AL, BB, (T)1, ETR );
        basic::internal::LocalGemm
        ( NORMAL, NORMAL, alpha, CT, DR, (T)1, ETR );
    }

    basic::internal::LocalGemm
    ( orientationOfA, orientationOfB, alpha, AL, BT, (T)0, FTL );
    basic::internal::LocalGemm
    ( NORMAL, NORMAL, alpha, CT, DL, (T)1, FTL );

    FTL.MakeTrapezoidal( LEFT, uplo );
    basic::Axpy( (T)1, FTL, ETL );

    basic::internal::LocalGemm
    ( orientationOfA, orientationOfB, alpha, AR, BB, (T)0, FBR );
    basic::internal::LocalGemm
    ( NORMAL, NORMAL, alpha, CB, DR, (T)1, FBR );

    FBR.MakeTrapezoidal( LEFT, uplo );
    basic::Axpy( (T)1, FBR, EBR );
    //------------------------------------------------------------------------//
#ifndef RELEASE
    PopCallStack();
#endif
}

// E := alpha (A^{T/H} B^{T/H} + C D^{T/H}) + beta C
template<typename T>
inline void
LocalTrr2kKernel
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfB,
  Orientation orientationOfD,
  T alpha, const DistMatrix<T,STAR,MC  >& A, const DistMatrix<T,MR,STAR>& B,
           const DistMatrix<T,MC,  STAR>& C, const DistMatrix<T,STAR,MR>& D,
  T beta,        DistMatrix<T,MC,  MR  >& E )
{
#ifndef RELEASE
    PushCallStack("LocalTrr2kKernel");
    CheckInput( orientationOfA, orientationOfB, orientationOfB, A, B, C, D, E );
#endif
    const Grid& g = E.Grid();

    DistMatrix<T,STAR,MC> AL(g), AR(g);
    DistMatrix<T,MR,STAR> BT(g), 
                          BB(g);
    DistMatrix<T,MC,STAR> CT(g),
                          CB(g);
    DistMatrix<T,MR,STAR> DT(g), 
                          DB(g);
    DistMatrix<T,MC,MR> ETL(g), ETR(g),
                        EBL(g), EBR(g);
    DistMatrix<T,MC,MR> FTL(g), FBR(g);

    const unsigned half = E.Height()/2;
    E.ScaleTrapezoid( beta, LEFT, uplo );
    LockedPartitionRight( A, AL, AR, half );
    LockedPartitionDown
    ( B, BT,
         BB, half );
    LockedPartitionDown
    ( C, CT,
         CB, half );
    LockedPartitionDown
    ( D, DT,
         DB, half );
    PartitionDownDiagonal
    ( E, ETL, ETR,
         EBL, EBR, half );

    FTL.AlignWith( ETL );
    FBR.AlignWith( EBR );
    FTL.ResizeTo( ETL.Height(), ETL.Width() );
    FBR.ResizeTo( EBR.Height(), EBR.Width() );
    //------------------------------------------------------------------------//
    if( uplo == LOWER )
    {
        basic::internal::LocalGemm
        ( orientationOfA, orientationOfB, alpha, AR, BT, (T)1, EBL );
        basic::internal::LocalGemm
        ( NORMAL, orientationOfD, alpha, CB, DT, (T)1, EBL );
    }
    else
    {
        basic::internal::LocalGemm
        ( orientationOfA, orientationOfB, alpha, AL, BB, (T)1, ETR );
        basic::internal::LocalGemm
        ( NORMAL, orientationOfD, alpha, CT, DB, (T)1, ETR );
    }

    basic::internal::LocalGemm
    ( orientationOfA, orientationOfB, alpha, AL, BT, (T)0, FTL );
    basic::internal::LocalGemm
    ( NORMAL, orientationOfD, alpha, CT, DT, (T)1, FTL );

    FTL.MakeTrapezoidal( LEFT, uplo );
    basic::Axpy( (T)1, FTL, ETL );

    basic::internal::LocalGemm
    ( orientationOfA, orientationOfB, alpha, AR, BB, (T)0, FBR );
    basic::internal::LocalGemm
    ( NORMAL, orientationOfD, alpha, CB, DB, (T)1, FBR );

    FBR.MakeTrapezoidal( LEFT, uplo );
    basic::Axpy( (T)1, FBR, EBR );
    //------------------------------------------------------------------------//
#ifndef RELEASE
    PopCallStack();
#endif
}

// E := alpha (A^{T/H} B^{T/H} + C^{T/H} D) + beta E
template<typename T>
inline void
LocalTrr2kKernel
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfB,
  Orientation orientationOfC,
  T alpha, const DistMatrix<T,STAR,MC>& A, const DistMatrix<T,MR,STAR>& B,
           const DistMatrix<T,STAR,MC>& C, const DistMatrix<T,STAR,MR>& D,
  T beta,        DistMatrix<T,MC,  MR>& E )
{
#ifndef RELEASE
    PushCallStack("LocalTrr2kKernel");
    CheckInput( orientationOfA, orientationOfB, orientationOfC, A, B, C, D, E );
#endif
    const Grid& g = E.Grid();

    DistMatrix<T,STAR,MC> AL(g), AR(g),
                          CL(g), CR(g);
    DistMatrix<T,MR,STAR> BT(g),
                          BB(g);
    DistMatrix<T,STAR,MR> DL(g), DR(g);
    DistMatrix<T,MC,MR> ETL(g), ETR(g),
                        EBL(g), EBR(g);
    DistMatrix<T,MC,MR> FTL(g), FBR(g);

    const unsigned half = E.Height()/2;
    E.ScaleTrapezoid( beta, LEFT, uplo );
    LockedPartitionRight( A, AL, AR, half );
    LockedPartitionDown
    ( B, BT,
         BB, half );
    LockedPartitionRight( C, CL, CR, half );
    LockedPartitionRight( D, DL, DR, half );
    PartitionDownDiagonal
    ( E, ETL, ETR,
         EBL, EBR, half );

    FTL.AlignWith( ETL );
    FBR.AlignWith( EBR );
    FTL.ResizeTo( ETL.Height(), ETL.Width() );
    FBR.ResizeTo( EBR.Height(), EBR.Width() );
    //------------------------------------------------------------------------//
    if( uplo == LOWER )
    {
        basic::internal::LocalGemm
        ( orientationOfA, orientationOfB, alpha, AR, BT, (T)1, EBL );
        basic::internal::LocalGemm
        ( orientationOfC, NORMAL, alpha, CR, DL, (T)1, EBL );
    }
    else
    {
        basic::internal::LocalGemm
        ( orientationOfA, orientationOfB, alpha, AL, BB, (T)1, ETR );
        basic::internal::LocalGemm
        ( orientationOfC, NORMAL, alpha, CL, DR, (T)1, ETR );
    }

    basic::internal::LocalGemm
    ( orientationOfA, orientationOfB, alpha, AL, BT, (T)0, FTL );
    basic::internal::LocalGemm
    ( orientationOfC, NORMAL, alpha, CL, DL, (T)1, FTL );

    FTL.MakeTrapezoidal( LEFT, uplo );
    basic::Axpy( (T)1, FTL, ETL );

    basic::internal::LocalGemm
    ( orientationOfA, orientationOfB, alpha, AR, BB, (T)0, FBR );
    basic::internal::LocalGemm
    ( orientationOfC, NORMAL, alpha, CR, DR, (T)1, FBR );

    FBR.MakeTrapezoidal( LEFT, uplo );
    basic::Axpy( (T)1, FBR, EBR );
    //------------------------------------------------------------------------//
#ifndef RELEASE
    PopCallStack();
#endif
}

// E := alpha (A^{T/H} B^{T/H} + C^{T/H} D^{T/H}) + beta C
template<typename T>
inline void
LocalTrr2kKernel
( UpperOrLower uplo,
  Orientation orientationOfA, Orientation orientationOfB,
  Orientation orientationOfC, Orientation orientationOfD,
  T alpha, const DistMatrix<T,STAR,MC>& A, const DistMatrix<T,MR,STAR>& B,
           const DistMatrix<T,STAR,MC>& C, const DistMatrix<T,MR,STAR>& D,
  T beta,        DistMatrix<T,MC,  MR>& E )
{
#ifndef RELEASE
    PushCallStack("LocalTrr2kKernel");
    CheckInput
    ( orientationOfA, orientationOfB, orientationOfC, orientationOfD, 
      A, B, C, D, E );
#endif
    const Grid& g = E.Grid();

    DistMatrix<T,STAR,MC> AL(g), AR(g),
                          CL(g), CR(g);
    DistMatrix<T,MR,STAR> BT(g),  DT(g),
                          BB(g),  DB(g);
    DistMatrix<T,MC,MR> ETL(g), ETR(g),
                        EBL(g), EBR(g);
    DistMatrix<T,MC,MR> FTL(g), FBR(g);

    const unsigned half = E.Height()/2;
    E.ScaleTrapezoid( beta, LEFT, uplo );
    LockedPartitionRight( A, AL, AR, half );
    LockedPartitionDown
    ( B, BT,
         BB, half );
    LockedPartitionRight( C, CL, CR, half );
    LockedPartitionDown
    ( D, DT,
         DB, half );
    PartitionDownDiagonal
    ( E, ETL, ETR,
         EBL, EBR, half );

    FTL.AlignWith( ETL );
    FBR.AlignWith( EBR );
    FTL.ResizeTo( ETL.Height(), ETL.Width() );
    FBR.ResizeTo( EBR.Height(), EBR.Width() );
    //------------------------------------------------------------------------//
    if( uplo == LOWER )
    {
        basic::internal::LocalGemm
        ( orientationOfA, orientationOfB, alpha, AR, BT, (T)1, EBL );
        basic::internal::LocalGemm
        ( orientationOfC, orientationOfD, alpha, CR, DT, (T)1, EBL );
    }
    else
    {
        basic::internal::LocalGemm
        ( orientationOfA, orientationOfB, alpha, AL, BB, (T)1, ETR );
        basic::internal::LocalGemm
        ( orientationOfC, orientationOfD, alpha, CL, DB, (T)1, ETR );
    }

    basic::internal::LocalGemm
    ( orientationOfA, orientationOfB, alpha, AL, BT, (T)0, FTL );
    basic::internal::LocalGemm
    ( orientationOfC, orientationOfD, alpha, CL, DT, (T)1, FTL );

    FTL.MakeTrapezoidal( LEFT, uplo );
    basic::Axpy( (T)1, FTL, ETL );

    basic::internal::LocalGemm
    ( orientationOfA, orientationOfB, alpha, AR, BB, (T)0, FBR );
    basic::internal::LocalGemm
    ( orientationOfC, orientationOfD, alpha, CR, DB, (T)1, FBR );

    FBR.MakeTrapezoidal( LEFT, uplo );
    basic::Axpy( (T)1, FBR, EBR );
    //------------------------------------------------------------------------//
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace trr2k
} // namespace basic
} // namespace elemental

// E := alpha (A B + C D) + beta E
template<typename T>
inline void
elemental::basic::internal::LocalTrr2k
( UpperOrLower uplo,
  T alpha, const DistMatrix<T,MC,  STAR>& A, const DistMatrix<T,STAR,MR>& B,
           const DistMatrix<T,MC,  STAR>& C, const DistMatrix<T,STAR,MR>& D,
  T beta,        DistMatrix<T,MC,  MR  >& E )
{
    using namespace trr2k;
#ifndef RELEASE
    PushCallStack("basic::internal::LocalTrr2k");
    CheckInput( A, B, C, D, E );
#endif
    const Grid& g = E.Grid();

    if( E.Height() < g.Width()*LocalTrr2kBlocksize<T>() )
    {
        LocalTrr2kKernel( uplo, alpha, A, B, C, D, beta, E );
    }
    else
    {
        // Split E in four roughly equal pieces, perform a large gemm on corner
        // and recurse on ETL and EBR.
        DistMatrix<T,MC,STAR> AT(g),  CT(g),
                              AB(g),  CB(g);
        DistMatrix<T,STAR,MR> BL(g), BR(g),
                              DL(g), DR(g);
        DistMatrix<T,MC,MR> ETL(g), ETR(g),
                            EBL(g), EBR(g);

        const unsigned half = E.Height() / 2;
        LockedPartitionDown
        ( A, AT,
             AB, half );
        LockedPartitionRight( B, BL, BR, half );
        LockedPartitionDown
        ( C, CT,
             CB, half );
        LockedPartitionRight( D, DL, DR, half );
        PartitionDownDiagonal
        ( E, ETL, ETR,
             EBL, EBR, half );

        if( uplo == LOWER )
        { 
            basic::internal::LocalGemm
            ( NORMAL, NORMAL, alpha, AB, BL, beta, EBL );
            basic::internal::LocalGemm
            ( NORMAL, NORMAL, alpha, CB, DL, (T)1, EBL );
        }
        else
        {
            basic::internal::LocalGemm
            ( NORMAL, NORMAL, alpha, AT, BR, beta, ETR );
            basic::internal::LocalGemm
            ( NORMAL, NORMAL, alpha, CT, DR, (T)1, ETR );
        }

        // Recurse
        basic::internal::LocalTrr2k( uplo, alpha, AT, BL, CT, DL, beta, ETL );
        basic::internal::LocalTrr2k( uplo, alpha, AB, BR, CB, DR, beta, EBR );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

// E := alpha (A B + C D^{T/H}) + beta E
template<typename T>
inline void
elemental::basic::internal::LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfD,
  T alpha, const DistMatrix<T,MC,  STAR>& A, const DistMatrix<T,STAR,MR>& B,
           const DistMatrix<T,MC,  STAR>& C, const DistMatrix<T,MR,STAR>& D,
  T beta,        DistMatrix<T,MC,  MR  >& E  )
{
    using namespace trr2k;
#ifndef RELEASE
    PushCallStack("basic::internal::LocalTrr2k");
    CheckInput( orientationOfD, A, B, C, D, E );
#endif
    const Grid& g = E.Grid();

    if( E.Height() < g.Width()*LocalTrr2kBlocksize<T>() )
    {
        LocalTrr2kKernel( uplo, orientationOfD, alpha, A, B, C, D, beta, E );
    }
    else
    {
        // Split E in four roughly equal pieces, perform a large gemm on corner
        // and recurse on ETL and EBR.
        DistMatrix<T,MC,STAR> AT(g),  CT(g),
                              AB(g),  CB(g);
        DistMatrix<T,STAR,MR> BL(g), BR(g);
        DistMatrix<T,MR,STAR> DT(g), 
                              DB(g);
        DistMatrix<T,MC,MR> ETL(g), ETR(g),
                            EBL(g), EBR(g);

        const unsigned half = E.Height() / 2;
        LockedPartitionDown
        ( A, AT,
             AB, half );
        LockedPartitionRight( B, BL, BR, half );
        LockedPartitionDown
        ( C, CT,
             CB, half );
        LockedPartitionDown
        ( D, DT,
             DB, half );
        PartitionDownDiagonal
        ( E, ETL, ETR,
             EBL, EBR, half );

        if( uplo == LOWER )
        { 
            basic::internal::LocalGemm
            ( NORMAL, NORMAL, alpha, AB, BL, (T)1, EBL );
            basic::internal::LocalGemm
            ( NORMAL, orientationOfD, alpha, CB, DT, beta, EBL );
        }
        else
        {
            basic::internal::LocalGemm
            ( NORMAL, NORMAL, alpha, AT, BR, (T)1, ETR );
            basic::internal::LocalGemm
            ( NORMAL, orientationOfD, alpha, CT, DB, beta, ETR );
        }

        // Recurse
        basic::internal::LocalTrr2k
        ( uplo, orientationOfD, alpha, AT, BL, CT, DT, beta, ETL );
        basic::internal::LocalTrr2k
        ( uplo, orientationOfD, alpha, AB, BR, CB, DB, beta, EBR );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

// E := alpha (A B + C^{T/H} D) + beta E
template<typename T>
inline void
elemental::basic::internal::LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfC,
  T alpha, const DistMatrix<T,MC,  STAR>& A, const DistMatrix<T,STAR,MR>& B,
           const DistMatrix<T,STAR,MC  >& C, const DistMatrix<T,STAR,MR>& D,
  T beta,        DistMatrix<T,MC,  MR  >& E  )
{
    using namespace trr2k;
#ifndef RELEASE
    PushCallStack("basic::internal::LocalTrr2k");
    CheckInput( orientationOfC, A, B, C, D, E );
#endif
    const Grid& g = E.Grid();

    if( E.Height() < g.Width()*LocalTrr2kBlocksize<T>() )
    {
        LocalTrr2kKernel( uplo, orientationOfC, alpha, A, B, C, D, beta, E );
    }
    else
    {
        // Split E in four roughly equal pieces, perform a large gemm on corner
        // and recurse on ETL and EBR.
        DistMatrix<T,MC,STAR> AT(g),
                              AB(g);
        DistMatrix<T,STAR,MR> BL(g), BR(g),
                              DL(g), DR(g);
        DistMatrix<T,STAR,MC> CL(g), CR(g);
        DistMatrix<T,MC,MR> ETL(g), ETR(g),
                            EBL(g), EBR(g);

        const unsigned half = E.Height() / 2;
        LockedPartitionDown
        ( A, AT,
             AB, half );
        LockedPartitionRight( B, BL, BR, half );
        LockedPartitionRight( C, CL, CR, half );
        LockedPartitionRight( D, DL, DR, half );
        PartitionDownDiagonal
        ( E, ETL, ETR,
             EBL, EBR, half );

        if( uplo == LOWER )
        { 
            basic::internal::LocalGemm
            ( NORMAL, NORMAL, alpha, AB, BL, beta, EBL );
            basic::internal::LocalGemm
            ( orientationOfC, NORMAL, alpha, CR, DL, (T)1, EBL );
        }
        else
        {
            basic::internal::LocalGemm
            ( NORMAL, NORMAL, alpha, AT, BR, beta, ETR );
            basic::internal::LocalGemm
            ( orientationOfC, NORMAL, alpha, CL, DR, (T)1, ETR );
        }

        // Recurse
        basic::internal::LocalTrr2k
        ( uplo, orientationOfC, alpha, AT, BL, CL, DL, beta, ETL );
        basic::internal::LocalTrr2k
        ( uplo, orientationOfC, alpha, AB, BR, CR, DR, beta, EBR );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

// E := alpha (A B + C^{T/H} D^{T/H}) + beta E
template<typename T>
inline void
elemental::basic::internal::LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfC,
  Orientation orientationOfD,
  T alpha, const DistMatrix<T,MC,  STAR>& A, const DistMatrix<T,STAR,MR>& B,
           const DistMatrix<T,STAR,MC  >& C, const DistMatrix<T,MR,STAR>& D,
  T beta,        DistMatrix<T,MC,  MR  >& E  )
{
    using namespace trr2k;
#ifndef RELEASE
    PushCallStack("basic::internal::LocalTrr2k");
    CheckInput( orientationOfC, orientationOfD, A, B, C, D, E );
#endif
    const Grid& g = E.Grid();

    if( E.Height() < g.Width()*LocalTrr2kBlocksize<T>() )
    {
        LocalTrr2kKernel
        ( uplo, orientationOfC, orientationOfD, alpha, A, B, C, D, beta, E );
    }
    else
    {
        // Split E in four roughly equal pieces, perform a large gemm on corner
        // and recurse on ETL and EBR.
        DistMatrix<T,MC,STAR> AT(g),
                              AB(g); 
        DistMatrix<T,STAR,MR> BL(g), BR(g);
        DistMatrix<T,STAR,MC> CL(g), CR(g);
        DistMatrix<T,MR,STAR> DT(g),
                              DB(g);
        DistMatrix<T,MC,MR> ETL(g), ETR(g),
                            EBL(g), EBR(g);

        const unsigned half = E.Height() / 2;
        LockedPartitionDown
        ( A, AT,
             AB, half );
        LockedPartitionRight( B, BL, BR, half );
        LockedPartitionRight( C, CL, CR, half );
        LockedPartitionDown
        ( D, DT,
             DB, half );
        PartitionDownDiagonal
        ( E, ETL, ETR,
             EBL, EBR, half );

        if( uplo == LOWER )
        { 
            basic::internal::LocalGemm
            ( NORMAL, NORMAL, alpha, AB, BL, beta, EBL );
            basic::internal::LocalGemm
            ( orientationOfC, orientationOfD, alpha, CR, DT, (T)1, EBL );
        }
        else
        {
            basic::internal::LocalGemm
            ( NORMAL, NORMAL, alpha, AT, BR, beta, ETR );
            basic::internal::LocalGemm
            ( orientationOfC, orientationOfD, alpha, CL, DB, (T)1, ETR );
        }

        // Recurse
        basic::internal::LocalTrr2k
        ( uplo, orientationOfC, orientationOfD, 
          alpha, AT, BL, CL, DT, beta, ETL );
        basic::internal::LocalTrr2k
        ( uplo, orientationOfC, orientationOfD,
          alpha, AB, BR, CR, DB, beta, EBR );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

// E := alpha (A B^{T/H} + C D) + beta E
template<typename T>
inline void
elemental::basic::internal::LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfB,
  T alpha, const DistMatrix<T,MC,  STAR>& A, const DistMatrix<T,MR,STAR>& B,
           const DistMatrix<T,MC,  STAR>& C, const DistMatrix<T,STAR,MR>& D,
  T beta,        DistMatrix<T,MC,  MR  >& E  )
{
    using namespace trr2k;
#ifndef RELEASE
    PushCallStack("basic::internal::LocalTrr2k");
    CheckInput( orientationOfB, A, B, C, D, E );
#endif
    const Grid& g = E.Grid();

    if( E.Height() < g.Width()*LocalTrr2kBlocksize<T>() )
    {
        LocalTrr2kKernel( uplo, orientationOfB, alpha, A, B, C, D, beta, E );
    }
    else
    {
        // Split E in four roughly equal pieces, perform a large gemm on corner
        // and recurse on ETL and EBR.
        DistMatrix<T,MC,STAR> AT(g),  CT(g),
                              AB(g),  CB(g);
        DistMatrix<T,MR,STAR> BT(g), 
                              BB(g);
        DistMatrix<T,STAR,MR> DL(g), DR(g);
        DistMatrix<T,MC,MR> ETL(g), ETR(g),
                            EBL(g), EBR(g);

        const unsigned half = E.Height() / 2;
        LockedPartitionDown
        ( A, AT,
             AB, half );
        LockedPartitionDown
        ( B, BT,
             BB, half );
        LockedPartitionDown
        ( C, CT,
             CB, half );
        LockedPartitionRight( D, DL, DR, half );
        PartitionDownDiagonal
        ( E, ETL, ETR,
             EBL, EBR, half );

        if( uplo == LOWER )
        { 
            basic::internal::LocalGemm
            ( NORMAL, orientationOfB, alpha, AB, BT, (T)1, EBL );
            basic::internal::LocalGemm
            ( NORMAL, NORMAL, alpha, CB, DL, beta, EBL );
        }
        else
        {
            basic::internal::LocalGemm
            ( NORMAL, orientationOfB, alpha, AT, BB, (T)1, ETR );
            basic::internal::LocalGemm
            ( NORMAL, NORMAL, alpha, CT, DR, beta, ETR );
        }

        // Recurse
        basic::internal::LocalTrr2k
        ( uplo, orientationOfB, alpha, AT, BT, CT, DL, beta, ETL );
        basic::internal::LocalTrr2k
        ( uplo, orientationOfB, alpha, AB, BB, CB, DR, beta, EBR );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

// E := alpha (A B^{T/H} + C D^{T/H}) + beta E
template<typename T>
inline void
elemental::basic::internal::LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfB,
  Orientation orientationOfD,
  T alpha, const DistMatrix<T,MC,STAR>& A, const DistMatrix<T,MR,STAR>& B,
           const DistMatrix<T,MC,STAR>& C, const DistMatrix<T,MR,STAR>& D,
  T beta,        DistMatrix<T,MC,MR  >& E  )
{
    using namespace trr2k;
#ifndef RELEASE
    PushCallStack("basic::internal::LocalTrr2k");
    CheckInput( orientationOfB, orientationOfD, A, B, C, D, E );
#endif
    const Grid& g = E.Grid();

    if( E.Height() < g.Width()*LocalTrr2kBlocksize<T>() )
    {
        LocalTrr2kKernel
        ( uplo, orientationOfB, orientationOfD, alpha, A, B, C, D, beta, E );
    }
    else
    {
        // Split E in four roughly equal pieces, perform a large gemm on corner
        // and recurse on ETL and EBR.
        DistMatrix<T,MC,STAR> AT(g),  CT(g),
                              AB(g),  CB(g);
        DistMatrix<T,MR,STAR> BT(g),  DT(g),
                              BB(g),  DB(g);
        DistMatrix<T,MC,MR> ETL(g), ETR(g),
                            EBL(g), EBR(g);

        const unsigned half = E.Height() / 2;
        LockedPartitionDown
        ( A, AT,
             AB, half );
        LockedPartitionDown
        ( B, BT,
             BB, half );
        LockedPartitionDown
        ( C, CT,
             CB, half );
        LockedPartitionDown
        ( D, DT,
             DB, half );
        PartitionDownDiagonal
        ( E, ETL, ETR,
             EBL, EBR, half );

        if( uplo == LOWER )
        { 
            basic::internal::LocalGemm
            ( NORMAL, orientationOfB, alpha, AB, BT, beta, EBL );
            basic::internal::LocalGemm
            ( NORMAL, orientationOfD, alpha, CB, DT, (T)1, EBL );
        }
        else
        {
            basic::internal::LocalGemm
            ( NORMAL, orientationOfB, alpha, AT, BB, beta, ETR );
            basic::internal::LocalGemm
            ( NORMAL, orientationOfD, alpha, CT, DB, (T)1, ETR );
        }

        // Recurse
        basic::internal::LocalTrr2k
        ( uplo, orientationOfB, orientationOfD,
          alpha, AT, BT, CT, DT, beta, ETL );
        basic::internal::LocalTrr2k
        ( uplo, orientationOfB, orientationOfD,
          alpha, AB, BB, CB, DB, beta, EBR );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

// E := alpha (A B^{T/H} + C^{T/H} D) + beta E
template<typename T>
inline void
elemental::basic::internal::LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfB,
  Orientation orientationOfC,
  T alpha, const DistMatrix<T,MC,  STAR>& A, const DistMatrix<T,MR,STAR>& B,
           const DistMatrix<T,STAR,MC  >& C, const DistMatrix<T,STAR,MR>& D,
  T beta,        DistMatrix<T,MC,  MR  >& E  )
{
    using namespace trr2k;
#ifndef RELEASE
    PushCallStack("basic::internal::LocalTrr2k");
    CheckInput( orientationOfB, orientationOfC, A, B, C, D, E );
#endif
    const Grid& g = E.Grid();

    if( E.Height() < g.Width()*LocalTrr2kBlocksize<T>() )
    {
        LocalTrr2kKernel
        ( uplo, orientationOfB, orientationOfC, alpha, A, B, C, D, beta, E );
    }
    else
    {
        // Split E in four roughly equal pieces, perform a large gemm on corner
        // and recurse on ETL and EBR.
        DistMatrix<T,MC,STAR> AT(g),
                              AB(g); 
        DistMatrix<T,MR,STAR> BT(g),
                              BB(g);
        DistMatrix<T,STAR,MC> CL(g), CR(g);
        DistMatrix<T,STAR,MR> DL(g), DR(g);
        DistMatrix<T,MC,MR> ETL(g), ETR(g),
                            EBL(g), EBR(g);

        const unsigned half = E.Height() / 2;
        LockedPartitionDown
        ( A, AT,
             AB, half );
        LockedPartitionDown
        ( B, BT,
             BB, half );
        LockedPartitionRight( C, CL, CR, half );
        LockedPartitionRight( D, DL, DR, half );
        PartitionDownDiagonal
        ( E, ETL, ETR,
             EBL, EBR, half );

        if( uplo == LOWER )
        { 
            basic::internal::LocalGemm
            ( NORMAL, orientationOfB, alpha, AB, BT, beta, EBL );
            basic::internal::LocalGemm
            ( orientationOfC, NORMAL, alpha, CR, DL, (T)1, EBL );
        }
        else
        {
            basic::internal::LocalGemm
            ( NORMAL, orientationOfB, alpha, AT, BB, beta, ETR );
            basic::internal::LocalGemm
            ( orientationOfC, NORMAL, alpha, CL, DR, (T)1, ETR );
        }

        // Recurse
        basic::internal::LocalTrr2k
        ( uplo, orientationOfB, orientationOfC,
          alpha, AT, BT, CL, DL, beta, ETL );
        basic::internal::LocalTrr2k
        ( uplo, orientationOfB, orientationOfC,
          alpha, AB, BB, CR, DR, beta, EBR );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

// E := alpha (A B^{T/H} + C^{T/H} D^{T/H}) + beta E
template<typename T>
inline void
elemental::basic::internal::LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfB,
  Orientation orientationOfC,
  Orientation orientationOfD,
  T alpha, const DistMatrix<T,MC,  STAR>& A, const DistMatrix<T,MR,STAR>& B, 
           const DistMatrix<T,STAR,MC  >& C, const DistMatrix<T,MR,STAR>& D,
  T beta,        DistMatrix<T,MC,  MR  >& E  )
{
    using namespace trr2k;
#ifndef RELEASE
    PushCallStack("basic::internal::LocalTrr2k");
    CheckInput( orientationOfB, orientationOfC, orientationOfD, A, B, C, D, E );
#endif
    const Grid& g = E.Grid();

    if( E.Height() < g.Width()*LocalTrr2kBlocksize<T>() )
    {
        LocalTrr2kKernel
        ( uplo, orientationOfB, orientationOfC, orientationOfD, 
          alpha, A, B, C, D, beta, E );
    }
    else
    {
        // Split E in four roughly equal pieces, perform a large gemm on corner
        // and recurse on ETL and EBR.
        DistMatrix<T,MC,STAR> AT(g),
                              AB(g); 
        DistMatrix<T,MR,STAR> BT(g),  DT(g),
                              BB(g),  DB(g);
        DistMatrix<T,STAR,MC> CL(g), CR(g);
        DistMatrix<T,MC,MR> ETL(g), ETR(g),
                            EBL(g), EBR(g);

        const unsigned half = E.Height() / 2;
        LockedPartitionDown
        ( A, AT,
             AB, half );
        LockedPartitionDown
        ( B, BT,
             BB, half );
        LockedPartitionRight( C, CL, CR, half );
        LockedPartitionDown
        ( D, DT,
             DB, half );
        PartitionDownDiagonal
        ( E, ETL, ETR,
             EBL, EBR, half );

        if( uplo == LOWER )
        { 
            basic::internal::LocalGemm
            ( NORMAL, orientationOfB, alpha, AB, BT, beta, EBL );
            basic::internal::LocalGemm
            ( orientationOfC, orientationOfD, alpha, CR, DT, (T)1, EBL );
        }
        else
        {
            basic::internal::LocalGemm
            ( NORMAL, orientationOfB, alpha, AT, BB, beta, ETR );
            basic::internal::LocalGemm
            ( orientationOfC, orientationOfD, alpha, CL, DB, (T)1, ETR );
        }

        // Recurse
        basic::internal::LocalTrr2k
        ( uplo, orientationOfB, orientationOfC, orientationOfD,
          alpha, AT, BT, CL, DT, beta, ETL );
        basic::internal::LocalTrr2k
        ( uplo, orientationOfB, orientationOfC, orientationOfD,
          alpha, AB, BB, CR, DB, beta, EBR );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

// E := alpha (A^{T/H} B + C D) + beta E
template<typename T>
inline void
elemental::basic::internal::LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfA,
  T alpha, const DistMatrix<T,STAR,MC  >& A, const DistMatrix<T,STAR,MR>& B,
           const DistMatrix<T,MC,  STAR>& C, const DistMatrix<T,STAR,MR>& D,
  T beta,        DistMatrix<T,MC,  MR  >& E  )
{
    using namespace trr2k;
#ifndef RELEASE
    PushCallStack("basic::internal::LocalTrr2k");
    CheckInput( orientationOfA, A, B, C, D, E );
#endif
    const Grid& g = E.Grid();

    if( E.Height() < g.Width()*LocalTrr2kBlocksize<T>() )
    {
        LocalTrr2kKernel( uplo, orientationOfA, alpha, A, B, C, D, beta, E );
    }
    else
    {
        // Split E in four roughly equal pieces, perform a large gemm on corner
        // and recurse on ETL and EBR.
        DistMatrix<T,STAR,MC> AL(g), AR(g);
        DistMatrix<T,STAR,MR> BL(g), BR(g),
                              DL(g), DR(g);
        DistMatrix<T,MC,STAR> CT(g),
                              CB(g);
        DistMatrix<T,MC,MR> ETL(g), ETR(g),
                            EBL(g), EBR(g);

        const unsigned half = E.Height() / 2;
        LockedPartitionRight( A, AL, AR, half );
        LockedPartitionRight( B, BL, BR, half );
        LockedPartitionDown
        ( C, CT,
             CB, half );
        LockedPartitionRight( D, DL, DR, half );
        PartitionDownDiagonal
        ( E, ETL, ETR,
             EBL, EBR, half );

        if( uplo == LOWER )
        { 
            basic::internal::LocalGemm
            ( orientationOfA, NORMAL, alpha, AR, BL, beta, EBL );
            basic::internal::LocalGemm
            ( NORMAL, NORMAL, alpha, CB, DL, (T)1, EBL );
        }
        else
        {
            basic::internal::LocalGemm
            ( orientationOfA, NORMAL, alpha, AL, BR, beta, ETR );
            basic::internal::LocalGemm
            ( NORMAL, NORMAL, alpha, CT, DR, (T)1, ETR );
        }

        // Recurse
        basic::internal::LocalTrr2k
        ( uplo, orientationOfA, alpha, AL, BL, CT, DL, beta, ETL );
        basic::internal::LocalTrr2k
        ( uplo, orientationOfA, alpha, AR, BR, CB, DR, beta, EBR );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

// E := alpha (A^{T/H} B + C D^{T/H}) + beta E
template<typename T>
inline void
elemental::basic::internal::LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfD,
  T alpha, const DistMatrix<T,STAR,MC  >& A, const DistMatrix<T,STAR,MR>& B,
           const DistMatrix<T,MC,  STAR>& C, const DistMatrix<T,MR,STAR>& D,
  T beta,        DistMatrix<T,MC,  MR  >& E  )
{
    using namespace trr2k;
#ifndef RELEASE
    PushCallStack("basic::internal::LocalTrr2k");
    CheckInput( orientationOfA, orientationOfD, A, B, C, D, E );
#endif
    const Grid& g = E.Grid();

    if( E.Height() < g.Width()*LocalTrr2kBlocksize<T>() )
    {
        LocalTrr2kKernel
        ( uplo, orientationOfA, orientationOfD, alpha, A, B, C, D, beta, E );
    }
    else
    {
        // Split E in four roughly equal pieces, perform a large gemm on corner
        // and recurse on ETL and EBR.
        DistMatrix<T,STAR,MC> AL(g), AR(g);
        DistMatrix<T,STAR,MR> BL(g), BR(g);
        DistMatrix<T,MC,STAR> CT(g),
                              CB(g); 
        DistMatrix<T,MR,STAR> DT(g),
                              DB(g);
        DistMatrix<T,MC,MR> ETL(g), ETR(g),
                            EBL(g), EBR(g);

        const unsigned half = E.Height() / 2;
        LockedPartitionRight( A, AL, AR, half );
        LockedPartitionRight( B, BL, BR, half );
        LockedPartitionDown
        ( C, CT,
             CB, half );
        LockedPartitionDown
        ( D, DT,
             DB, half );
        PartitionDownDiagonal
        ( E, ETL, ETR,
             EBL, EBR, half );

        if( uplo == LOWER )
        { 
            basic::internal::LocalGemm
            ( orientationOfA, NORMAL, alpha, AR, BL, beta, EBL );
            basic::internal::LocalGemm
            ( NORMAL, orientationOfD, alpha, CB, DT, (T)1, EBL );
        }
        else
        {
            basic::internal::LocalGemm
            ( orientationOfA, NORMAL, alpha, AL, BR, beta, ETR );
            basic::internal::LocalGemm
            ( NORMAL, orientationOfD, alpha, CT, DB, (T)1, ETR );
        }

        // Recurse
        basic::internal::LocalTrr2k
        ( uplo, orientationOfA, orientationOfD,
          alpha, AL, BL, CT, DT, beta, ETL );
        basic::internal::LocalTrr2k
        ( uplo, orientationOfA, orientationOfD,
          alpha, AR, BR, CB, DB, beta, EBR );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

// E := alpha (A^{T/H} B + C^{T/H} D) + beta E
template<typename T>
inline void
elemental::basic::internal::LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfC,
  T alpha, const DistMatrix<T,STAR,MC>& A, const DistMatrix<T,STAR,MR>& B,
           const DistMatrix<T,STAR,MC>& C, const DistMatrix<T,STAR,MR>& D,
  T beta,        DistMatrix<T,MC,  MR>& E  )
{
    using namespace trr2k;
#ifndef RELEASE
    PushCallStack("basic::internal::LocalTrr2k");
    CheckInput( orientationOfA, orientationOfC, A, B, C, D, E );
#endif
    const Grid& g = E.Grid();

    if( E.Height() < g.Width()*LocalTrr2kBlocksize<T>() )
    {
        LocalTrr2kKernel
        ( uplo, orientationOfA, orientationOfC, alpha, A, B, C, D, beta, E );
    }
    else
    {
        // Split E in four roughly equal pieces, perform a large gemm on corner
        // and recurse on ETL and EBR.
        DistMatrix<T,STAR,MC> AL(g), AR(g),
                              CL(g), CR(g);
        DistMatrix<T,STAR,MR> BL(g), BR(g),
                              DL(g), DR(g);
        DistMatrix<T,MC,MR> ETL(g), ETR(g),
                            EBL(g), EBR(g);

        const unsigned half = E.Height() / 2;
        LockedPartitionRight( A, AL, AR, half );
        LockedPartitionRight( B, BL, BR, half );
        LockedPartitionRight( C, CL, CR, half );
        LockedPartitionRight( D, DL, DR, half );
        PartitionDownDiagonal
        ( E, ETL, ETR,
             EBL, EBR, half );

        if( uplo == LOWER )
        { 
            basic::internal::LocalGemm
            ( orientationOfA, NORMAL, alpha, AR, BL, beta, EBL );
            basic::internal::LocalGemm
            ( orientationOfC, NORMAL, alpha, CR, DL, (T)1, EBL );
        }
        else
        {
            basic::internal::LocalGemm
            ( orientationOfA, NORMAL, alpha, AL, BR, beta, ETR );
            basic::internal::LocalGemm
            ( orientationOfC, NORMAL, alpha, CL, DR, (T)1, ETR );
        }

        // Recurse
        basic::internal::LocalTrr2k
        ( uplo, orientationOfA, orientationOfC,
          alpha, AL, BL, CL, DL, beta, ETL );
        basic::internal::LocalTrr2k
        ( uplo, orientationOfA, orientationOfC, 
          alpha, AR, BR, CR, DR, beta, EBR );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

// E := alpha (A^{T/H} B + C^{T/H} D^{T/H}) + beta E
template<typename T>
inline void
elemental::basic::internal::LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfC,
  Orientation orientationOfD,
  T alpha, const DistMatrix<T,STAR,MC>& A, const DistMatrix<T,STAR,MR>& B,
           const DistMatrix<T,STAR,MC>& C, const DistMatrix<T,MR,STAR>& D,
  T beta,        DistMatrix<T,MC,  MR>& E  )
{
    using namespace trr2k;
#ifndef RELEASE
    PushCallStack("basic::internal::LocalTrr2k");
    CheckInput( orientationOfA, orientationOfC, orientationOfD, A, B, C, D, E );
#endif
    const Grid& g = E.Grid();

    if( E.Height() < g.Width()*LocalTrr2kBlocksize<T>() )
    {
        LocalTrr2kKernel
        ( uplo, orientationOfA, orientationOfC, orientationOfD,
          alpha, A, B, C, D, beta, E );
    }
    else
    {
        // Split E in four roughly equal pieces, perform a large gemm on corner
        // and recurse on ETL and EBR.
        DistMatrix<T,STAR,MC> AL(g), AR(g),
                              CL(g), CR(g);
        DistMatrix<T,STAR,MR> BL(g), BR(g);
        DistMatrix<T,MR,STAR> DT(g),
                              DB(g);
        DistMatrix<T,MC,MR> ETL(g), ETR(g),
                            EBL(g), EBR(g);

        const unsigned half = E.Height() / 2;
        LockedPartitionRight( A, AL, AR, half );
        LockedPartitionRight( B, BL, BR, half );
        LockedPartitionRight( C, CL, CR, half );
        LockedPartitionDown
        ( D, DT,
             DB, half );
        PartitionDownDiagonal
        ( E, ETL, ETR,
             EBL, EBR, half );

        if( uplo == LOWER )
        { 
            basic::internal::LocalGemm
            ( orientationOfA, NORMAL, alpha, AR, BL, beta, EBL );
            basic::internal::LocalGemm
            ( orientationOfC, orientationOfD, alpha, CR, DT, (T)1, EBL );
        }
        else
        {
            basic::internal::LocalGemm
            ( orientationOfA, NORMAL, alpha, AL, BR, beta, ETR );
            basic::internal::LocalGemm
            ( orientationOfC, orientationOfD, alpha, CL, DB, (T)1, ETR );
        }

        // Recurse
        basic::internal::LocalTrr2k
        ( uplo, orientationOfA, orientationOfC, orientationOfD,
          alpha, AL, BL, CL, DT, beta, ETL );
        basic::internal::LocalTrr2k
        ( uplo, orientationOfA, orientationOfC, orientationOfD,
          alpha, AR, BR, CR, DB, beta, EBR );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

// E := alpha (A^{T/H} B^{T/H} + C D) + beta E
template<typename T>
inline void
elemental::basic::internal::LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfB,
  T alpha, const DistMatrix<T,STAR,MC  >& A, const DistMatrix<T,MR,STAR>& B,
           const DistMatrix<T,MC,  STAR>& C, const DistMatrix<T,STAR,MR>& D,
  T beta,        DistMatrix<T,MC,  MR  >& E  )
{
    using namespace trr2k;
#ifndef RELEASE
    PushCallStack("basic::internal::LocalTrr2k");
    CheckInput( orientationOfA, orientationOfB, A, B, C, D, E );
#endif
    const Grid& g = E.Grid();

    if( E.Height() < g.Width()*LocalTrr2kBlocksize<T>() )
    {
        LocalTrr2kKernel
        ( uplo, orientationOfA, orientationOfB, alpha, A, B, C, D, beta, E );
    }
    else
    {
        // Split E in four roughly equal pieces, perform a large gemm on corner
        // and recurse on ETL and EBR.
        DistMatrix<T,STAR,MC> AL(g), AR(g);
        DistMatrix<T,MR,STAR> BT(g),
                              BB(g);
        DistMatrix<T,MC,STAR> CT(g),
                              CB(g); 
        DistMatrix<T,STAR,MR> DL(g), DR(g);
        DistMatrix<T,MC,MR> ETL(g), ETR(g),
                            EBL(g), EBR(g);

        const unsigned half = E.Height() / 2;
        LockedPartitionRight( A, AL, AR, half );
        LockedPartitionDown
        ( B, BT,
             BB, half );
        LockedPartitionDown
        ( C, CT,
             CB, half );
        LockedPartitionRight( D, DL, DR, half );
        PartitionDownDiagonal
        ( E, ETL, ETR,
             EBL, EBR, half );

        if( uplo == LOWER )
        { 
            basic::internal::LocalGemm
            ( orientationOfA, orientationOfB, alpha, AR, BT, beta, EBL );
            basic::internal::LocalGemm
            ( NORMAL, NORMAL, alpha, CB, DL, (T)1, EBL );
        }
        else
        {
            basic::internal::LocalGemm
            ( orientationOfA, orientationOfB, alpha, AL, BB, beta, ETR );
            basic::internal::LocalGemm
            ( NORMAL, NORMAL, alpha, CT, DR, (T)1, ETR );
        }

        // Recurse
        basic::internal::LocalTrr2k
        ( uplo, orientationOfA, orientationOfB,
          alpha, AL, BT, CT, DL, beta, ETL );
        basic::internal::LocalTrr2k
        ( uplo, orientationOfA, orientationOfB,
          alpha, AR, BB, CB, DR, beta, EBR );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

// E := alpha (A^{T/H} B^{T/H} + C D^{T/H}) + beta E
template<typename T>
inline void
elemental::basic::internal::LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfB,
  Orientation orientationOfD,
  T alpha, const DistMatrix<T,STAR,MC  >& A, const DistMatrix<T,MR,STAR>& B,
           const DistMatrix<T,MC,  STAR>& C, const DistMatrix<T,MR,STAR>& D,
  T beta,        DistMatrix<T,MC,  MR  >& E  )
{
    using namespace trr2k;
#ifndef RELEASE
    PushCallStack("basic::internal::LocalTrr2k");
    CheckInput( orientationOfA, orientationOfB, orientationOfD, A, B, C, D, E );
#endif
    const Grid& g = E.Grid();

    if( E.Height() < g.Width()*LocalTrr2kBlocksize<T>() )
    {
        LocalTrr2kKernel
        ( uplo, orientationOfA, orientationOfB, orientationOfD, 
          alpha, A, B, C, D, beta, E );
    }
    else
    {
        // Split E in four roughly equal pieces, perform a large gemm on corner
        // and recurse on ETL and EBR.
        DistMatrix<T,STAR,MC> AL(g), AR(g);
        DistMatrix<T,MR,STAR> BT(g),  DT(g),
                              BB(g),  DB(g);
        DistMatrix<T,MC,STAR> CT(g),
                              CB(g); 
        DistMatrix<T,MC,MR> ETL(g), ETR(g),
                            EBL(g), EBR(g);

        const unsigned half = E.Height() / 2;
        LockedPartitionRight( A, AL, AR, half );
        LockedPartitionDown
        ( B, BT,
             BB, half );
        LockedPartitionDown
        ( C, CT,
             CB, half );
        LockedPartitionDown
        ( D, DT,
             DB, half );
        PartitionDownDiagonal
        ( E, ETL, ETR,
             EBL, EBR, half );

        if( uplo == LOWER )
        { 
            basic::internal::LocalGemm
            ( orientationOfA, orientationOfB, alpha, AR, BT, beta, EBL );
            basic::internal::LocalGemm
            ( NORMAL, orientationOfD, alpha, CB, DT, (T)1, EBL );
        }
        else
        {
            basic::internal::LocalGemm
            ( orientationOfA, orientationOfB, alpha, AL, BB, beta, ETR );
            basic::internal::LocalGemm
            ( NORMAL, orientationOfD, alpha, CT, DB, (T)1, ETR );
        }

        // Recurse
        basic::internal::LocalTrr2k
        ( uplo, orientationOfA, orientationOfB, orientationOfD,
          alpha, AL, BT, CT, DT, beta, ETL );
        basic::internal::LocalTrr2k
        ( uplo, orientationOfA, orientationOfB, orientationOfD,
          alpha, AR, BB, CB, DB, beta, EBR );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

// E := alpha (A^{T/H} B^{T/H} + C^{T/H} D) + beta E
template<typename T>
inline void
elemental::basic::internal::LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfB,
  Orientation orientationOfC,
  T alpha, const DistMatrix<T,STAR,MC>& A, const DistMatrix<T,MR,STAR>& B,
           const DistMatrix<T,STAR,MC>& C, const DistMatrix<T,STAR,MR>& D,
  T beta,        DistMatrix<T,MC,  MR>& E  )
{
    using namespace trr2k;
#ifndef RELEASE
    PushCallStack("basic::internal::LocalTrr2k");
    CheckInput( orientationOfA, orientationOfB, orientationOfC, A, B, C, D, E );
#endif
    const Grid& g = E.Grid();

    if( E.Height() < g.Width()*LocalTrr2kBlocksize<T>() )
    {
        LocalTrr2kKernel
        ( uplo, orientationOfA, orientationOfB, orientationOfC,
          alpha, A, B, C, D, beta, E );
    }
    else
    {
        // Split E in four roughly equal pieces, perform a large gemm on corner
        // and recurse on ETL and EBR.
        DistMatrix<T,STAR,MC> AL(g), AR(g),
                              CL(g), CR(g);
        DistMatrix<T,MR,STAR> BT(g),
                              BB(g);
        DistMatrix<T,STAR,MR> DL(g), DR(g);
        DistMatrix<T,MC,MR> ETL(g), ETR(g),
                            EBL(g), EBR(g);

        const unsigned half = E.Height() / 2;
        LockedPartitionRight( A, AL, AR, half );
        LockedPartitionDown
        ( B, BT, 
             BB, half );
        LockedPartitionRight( C, CL, CR, half );
        LockedPartitionRight( D, DL, DR, half );
        PartitionDownDiagonal
        ( E, ETL, ETR,
             EBL, EBR, half );

        if( uplo == LOWER )
        { 
            basic::internal::LocalGemm
            ( orientationOfA, orientationOfB, alpha, AR, BT, beta, EBL );
            basic::internal::LocalGemm
            ( orientationOfC, NORMAL, alpha, CR, DL, (T)1, EBL );
        }
        else
        {
            basic::internal::LocalGemm
            ( orientationOfA, orientationOfB, alpha, AL, BB, beta, ETR );
            basic::internal::LocalGemm
            ( orientationOfC, NORMAL, alpha, CL, DR, (T)1, ETR );
        }

        // Recurse
        basic::internal::LocalTrr2k
        ( uplo, orientationOfA, orientationOfB, orientationOfC,
          alpha, AL, BT, CL, DL, beta, ETL );
        basic::internal::LocalTrr2k
        ( uplo, orientationOfA, orientationOfB, orientationOfC,
          alpha, AR, BB, CR, DR, beta, EBR );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

// E := alpha (A^{T/H} B^{T/H} + C^{T/H} D^{T/H}) + beta E
template<typename T>
inline void
elemental::basic::internal::LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfB,
  Orientation orientationOfC,
  Orientation orientationOfD,
  T alpha, const DistMatrix<T,STAR,MC>& A, const DistMatrix<T,MR,STAR>& B,
           const DistMatrix<T,STAR,MC>& C, const DistMatrix<T,MR,STAR>& D,
  T beta,        DistMatrix<T,MC,  MR>& E  )
{
    using namespace trr2k;
#ifndef RELEASE
    PushCallStack("basic::internal::LocalTrr2k");
    CheckInput
    ( orientationOfA, orientationOfB, orientationOfC, orientationOfD, 
      A, B, C, D, E );
#endif
    const Grid& g = E.Grid();

    if( E.Height() < g.Width()*LocalTrr2kBlocksize<T>() )
    {
        LocalTrr2kKernel
        ( uplo, orientationOfA, orientationOfB, orientationOfC, orientationOfD,
          alpha, A, B, C, D, beta, E );
    }
    else
    {
        // Split E in four roughly equal pieces, perform a large gemm on corner
        // and recurse on ETL and EBR.
        DistMatrix<T,STAR,MC> AL(g), AR(g),
                              CL(g), CR(g);
        DistMatrix<T,MR,STAR> BT(g),  DT(g),
                              BB(g),  DB(g);
        DistMatrix<T,MC,MR> ETL(g), ETR(g),
                            EBL(g), EBR(g);

        const unsigned half = E.Height() / 2;
        LockedPartitionRight( A, AL, AR, half );
        LockedPartitionDown
        ( B, BT,
             BB, half );
        LockedPartitionRight( C, CL, CR, half );
        LockedPartitionDown
        ( D, DT,
             DB, half );
        PartitionDownDiagonal
        ( E, ETL, ETR,
             EBL, EBR, half );

        if( uplo == LOWER )
        { 
            basic::internal::LocalGemm
            ( orientationOfA, orientationOfB, alpha, AR, BT, beta, EBL );
            basic::internal::LocalGemm
            ( orientationOfC, orientationOfD, alpha, CR, DT, (T)1, EBL );
        }
        else
        {
            basic::internal::LocalGemm
            ( orientationOfA, orientationOfB, alpha, AL, BB, beta, ETR );
            basic::internal::LocalGemm
            ( orientationOfC, orientationOfD, alpha, CL, DB, (T)1, ETR );
        }

        // Recurse
        basic::internal::LocalTrr2k
        ( uplo, 
          orientationOfA, orientationOfB, orientationOfC, orientationOfD, 
          alpha, AL, BT, CL, DT, beta, ETL );

        basic::internal::LocalTrr2k
        ( uplo, 
          orientationOfA, orientationOfB, orientationOfC, orientationOfD,
          alpha, AR, BB, CR, DB, beta, EBR );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}
