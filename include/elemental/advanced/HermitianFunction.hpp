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

#ifndef WITHOUT_PMRRR

namespace elemental {
namespace advanced {
namespace hermitian_function {

// A :=  Z f(Omega) Z^T, where f is real-valued
template<typename R,class RealFunctor>
inline void
ReformHermitianMatrix
( UpperOrLower uplo,
        DistMatrix<R,MC,MR>& A,
  const DistMatrix<R,VR,STAR>& w,
  const DistMatrix<R,MC,MR>& Z,
  const RealFunctor& f )
{
#ifndef RELEASE
    PushCallStack("advanced::hermitian_function::ReformHermitianMatrix");
#endif
    const Grid& g = A.Grid();

    DistMatrix<R,MC,MR> ZL(g), ZR(g),
                        Z0(g), Z1(g), Z2(g);
    DistMatrix<R,VR,STAR> wT(g),  w0(g),
                          wB(g),  w1(g),
                                  w2(g);

    DistMatrix<R,MC,  STAR> Z1_MC_STAR(g);
    DistMatrix<R,VR,  STAR> Z1_VR_STAR(g);
    DistMatrix<R,STAR,MR  > Z1Trans_STAR_MR(g);
    DistMatrix<R,STAR,STAR> w1_STAR_STAR(g);

    if( uplo == LOWER )
        A.MakeTrapezoidal( LEFT, UPPER, 1 );
    else
        A.MakeTrapezoidal( LEFT, LOWER, -1 );
    LockedPartitionRight( Z, ZL, ZR, 0 );
    LockedPartitionDown
    ( w, wT,
         wB, 0 );
    while( ZL.Width() < Z.Width() )
    {
        LockedRepartitionRight
        ( ZL, /**/ ZR,
          Z0, /**/ Z1, Z2 );
        LockedRepartitionDown
        ( wT,  w0,
         /**/ /**/
               w1,
          wB,  w2 );

        Z1_MC_STAR.AlignWith( A );
        Z1_VR_STAR.AlignWith( A );
        Z1Trans_STAR_MR.AlignWith( A );
        //--------------------------------------------------------------------//
        Z1_MC_STAR = Z1;
        Z1_VR_STAR = Z1_MC_STAR;
        w1_STAR_STAR = w1;

        // Scale Z1[VR,* ] with the modified eigenvalues
        const int width = Z1_VR_STAR.Width();
        const int localHeight = Z1_VR_STAR.LocalHeight();
        for( int j=0; j<width; ++j )
        {
            const R omega = f(w1_STAR_STAR.GetLocalEntry(j,0));
            R* buffer = Z1_VR_STAR.LocalBuffer(0,j);
            for( int iLocal=0; iLocal<localHeight; ++iLocal )
                buffer[iLocal] *= omega;
        }

        Z1Trans_STAR_MR.TransposeFrom( Z1_VR_STAR );
        basic::internal::LocalTrrk
        ( uplo, (R)1, Z1_MC_STAR, Z1Trans_STAR_MR, (R)1, A );
        //--------------------------------------------------------------------//
        Z1Trans_STAR_MR.FreeAlignments();
        Z1_VR_STAR.FreeAlignments();
        Z1_MC_STAR.FreeAlignments();

        SlideLockedPartitionDown
        ( wT,  w0,
               w1,
         /**/ /**/
          wB,  w2 );
        SlideLockedPartitionRight
        ( ZL,     /**/ ZR,
          Z0, Z1, /**/ Z2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
// A :=  Z f(Omega) Z^T, where f is real-valued
template<typename R,class RealFunctor>
inline void
ReformHermitianMatrix
( UpperOrLower uplo,
        DistMatrix<std::complex<R>,MC,MR>& A,
  const DistMatrix<R,VR,STAR>& w,
  const DistMatrix<std::complex<R>,MC,MR>& Z,
  const RealFunctor& f )
{
#ifndef RELEASE
    PushCallStack("advanced::hermitian_function::ReformHermitianMatrix");
#endif
    const Grid& g = A.Grid();
    typedef std::complex<R> C;

    DistMatrix<C,MC,MR> ZL(g), ZR(g),
                        Z0(g), Z1(g), Z2(g);
    DistMatrix<R,VR,STAR> wT(g),  w0(g),
                          wB(g),  w1(g),
                                  w2(g);

    DistMatrix<C,MC,  STAR> Z1_MC_STAR(g);
    DistMatrix<C,VR,  STAR> Z1_VR_STAR(g);
    DistMatrix<C,STAR,MR  > Z1Adj_STAR_MR(g);
    DistMatrix<R,STAR,STAR> w1_STAR_STAR(g);

    if( uplo == LOWER )
        A.MakeTrapezoidal( LEFT, UPPER, 1 );
    else
        A.MakeTrapezoidal( LEFT, LOWER, -1 );
    LockedPartitionRight( Z, ZL, ZR, 0 );
    LockedPartitionDown
    ( w, wT,
         wB, 0 );
    while( ZL.Width() < Z.Width() )
    {
        LockedRepartitionRight
        ( ZL, /**/ ZR,
          Z0, /**/ Z1, Z2 );
        LockedRepartitionDown
        ( wT,  w0,
         /**/ /**/
               w1,
          wB,  w2 );

        Z1_MC_STAR.AlignWith( A );
        Z1_VR_STAR.AlignWith( A );
        Z1Adj_STAR_MR.AlignWith( A );
        //--------------------------------------------------------------------//
        Z1_MC_STAR = Z1;
        Z1_VR_STAR = Z1_MC_STAR;
        w1_STAR_STAR = w1;

        // Scale Z1[VR,* ] with the modified eigenvalues
        const int width = Z1_VR_STAR.Width();
        const int localHeight = Z1_VR_STAR.LocalHeight();
        for( int j=0; j<width; ++j )
        {
            const R omega = f(w1_STAR_STAR.GetLocalEntry(j,0));
            C* buffer = Z1_VR_STAR.LocalBuffer(0,j);
            for( int iLocal=0; iLocal<localHeight; ++iLocal )
                buffer[iLocal] *= omega;
        }

        Z1Adj_STAR_MR.AdjointFrom( Z1_VR_STAR );
        basic::internal::LocalTrrk
        ( uplo, (C)1, Z1_MC_STAR, Z1Adj_STAR_MR, (C)1, A );
        //--------------------------------------------------------------------//
        Z1Adj_STAR_MR.FreeAlignments();
        Z1_VR_STAR.FreeAlignments();
        Z1_MC_STAR.FreeAlignments();

        SlideLockedPartitionDown
        ( wT,  w0,
               w1,
         /**/ /**/
          wB,  w2 );
        SlideLockedPartitionRight
        ( ZL,     /**/ ZR,
          Z0, Z1, /**/ Z2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

// A :=  Z f(Omega) Z^T, where f is complex-valued
template<typename R,class ComplexFunctor>
inline void
ReformNormalMatrix
(       DistMatrix<std::complex<R>,MC,MR>& A,
  const DistMatrix<R,VR,STAR>& w,
  const DistMatrix<std::complex<R>,MC,MR>& Z,
  const ComplexFunctor& f )
{
#ifndef RELEASE
    PushCallStack("advanced::hermitian_function::ReformNormalMatrix");
#endif
    const Grid& g = A.Grid();
    typedef std::complex<R> C;

    DistMatrix<C,MC,MR> ZL(g), ZR(g),
                        Z0(g), Z1(g), Z2(g);
    DistMatrix<R,VR,STAR> wT(g),  w0(g),
                          wB(g),  w1(g),
                                  w2(g);

    DistMatrix<C,MC,  STAR> Z1_MC_STAR(g);
    DistMatrix<C,VR,  STAR> Z1_VR_STAR(g);
    DistMatrix<C,STAR,MR  > Z1Adj_STAR_MR(g);
    DistMatrix<R,STAR,STAR> w1_STAR_STAR(g);

    A.SetToZero();
    LockedPartitionRight( Z, ZL, ZR, 0 );
    LockedPartitionDown
    ( w, wT,
         wB, 0 );
    while( ZL.Width() < Z.Width() )
    {
        LockedRepartitionRight
        ( ZL, /**/ ZR,
          Z0, /**/ Z1, Z2 );
        LockedRepartitionDown
        ( wT,  w0,
         /**/ /**/
               w1,
          wB,  w2 );

        Z1_MC_STAR.AlignWith( A );
        Z1_VR_STAR.AlignWith( A );
        Z1Adj_STAR_MR.AlignWith( A );
        //--------------------------------------------------------------------//
        Z1_MC_STAR = Z1;
        Z1_VR_STAR = Z1_MC_STAR;
        w1_STAR_STAR = w1;

        // Scale Z1[VR,* ] with the modified eigenvalues
        const int width = Z1_VR_STAR.Width();
        const int localHeight = Z1_VR_STAR.LocalHeight();
        for( int j=0; j<width; ++j )
        {
            const C conjOmega = Conj(f(w1_STAR_STAR.GetLocalEntry(j,0)));
            C* buffer = Z1_VR_STAR.LocalBuffer(0,j);
            for( int iLocal=0; iLocal<localHeight; ++iLocal )
                buffer[iLocal] *= conjOmega;
        }

        Z1Adj_STAR_MR.AdjointFrom( Z1_VR_STAR );
        basic::internal::LocalGemm
        ( NORMAL, NORMAL, (C)1, Z1_MC_STAR, Z1Adj_STAR_MR, (C)1, A );
        //--------------------------------------------------------------------//
        Z1Adj_STAR_MR.FreeAlignments();
        Z1_VR_STAR.FreeAlignments();
        Z1_MC_STAR.FreeAlignments();

        SlideLockedPartitionDown
        ( wT,  w0,
               w1,
         /**/ /**/
          wB,  w2 );
        SlideLockedPartitionRight
        ( ZL,     /**/ ZR,
          Z0, Z1, /**/ Z2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX

} // namespace hermitian_eig
} // namespace advanced
} // namespace elemental

//
// Modify the eigenvalues of A with the real-valued function f, which will 
// therefore result in a Hermitian matrix, which we store in-place.
//

template<typename R,class RealFunctor>
inline void
elemental::advanced::RealHermitianFunction
( UpperOrLower uplo, DistMatrix<R,MC,MR>& A, const RealFunctor& f )
{
#ifndef RELEASE
    PushCallStack("advanced::RealHermitianFunction");
#endif
    if( A.Height() != A.Width() )
        throw std::logic_error("Hermitian matrices must be square");

    // Get the EVD of A
    const Grid& g = A.Grid();
    DistMatrix<R,VR,STAR> w(g);
    DistMatrix<R,MC,MR> Z(g);
    HermitianEig( uplo, A, w, Z );

    // Form the custom outer product, Z f(Omega) Z^T
    hermitian_function::ReformHermitianMatrix( uplo, A, w, Z, f );
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
template<typename R,class RealFunctor>
inline void
elemental::advanced::RealHermitianFunction
( UpperOrLower uplo, 
  DistMatrix<std::complex<R>,MC,MR>& A, const RealFunctor& f )
{
#ifndef RELEASE
    PushCallStack("advanced::RealHermitianFunction");
#endif
    if( A.Height() != A.Width() )
        throw std::logic_error("Hermitian matrices must be square");

    // Get the EVD of A
    const Grid& g = A.Grid();
    DistMatrix<R,VR,STAR> w(g);
    DistMatrix<std::complex<R>,MC,MR> Z(g);
    HermitianEig( uplo, A, w, Z );

    // Form the custom outer product, Z f(Omega) Z^H
    hermitian_function::ReformHermitianMatrix( uplo, A, w, Z, f );
#ifndef RELEASE
    PopCallStack();
#endif
}

//
// Modify the eigenvalues of A with the complex-valued function f, which will
// therefore result in a normal (in general, non-Hermitian) matrix, which we 
// store in-place. At some point a version will be written which takes a real
// symmetric matrix as input and produces a complex normal matrix.
//

template<typename R,class ComplexFunctor>
inline void
elemental::advanced::ComplexHermitianFunction
( UpperOrLower uplo, DistMatrix<std::complex<R>,MC,MR>& A, 
  const ComplexFunctor& f )
{
#ifndef RELEASE
    PushCallStack("advanced::ComplexHermitianFunction");
#endif
    if( A.Height() != A.Width() )
        throw std::logic_error("Hermitian matrices must be square");

    // Get the EVD of A
    const Grid& g = A.Grid();
    DistMatrix<R,VR,STAR> w(g);
    DistMatrix<std::complex<R>,MC,MR> Z(g);
    HermitianEig( uplo, A, w, Z );

    // Form the custom outer product, Z f(Omega) Z^H
    hermitian_function::ReformNormalMatrix( A, w, Z, f );
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX
#endif // WITHOUT_PMRRR
