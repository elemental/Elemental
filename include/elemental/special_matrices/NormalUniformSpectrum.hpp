/*
   Copyright (c) 2009-2012, Jack Poulson
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

namespace elem {

// Draw the spectrum from the specified half-open interval on the real line,
// then rotate it with a random Householder similarity transformation

template<typename R>
inline void
NormalUniformSpectrum
( int n, Matrix<Complex<R> >& A, Complex<R> center, R radius )
{
#ifndef RELEASE
    PushCallStack("NormalUniformSpectrum");
#endif
    A.ResizeTo( n, n );
    MakeNormalUniformSpectrum( A, center, radius );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R,Distribution U,Distribution V>
inline void
NormalUniformSpectrum
( int n, DistMatrix<Complex<R>,U,V>& A, Complex<R> center, R radius )
{
#ifndef RELEASE
    PushCallStack("NormalUniformSpectrum");
#endif
    A.ResizeTo( n, n );
    MakeNormalUniformSpectrum( A, center, radius );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline void
MakeNormalUniformSpectrum
( Matrix<Complex<R> >& A, Complex<R> center, R radius )
{
#ifndef RELEASE
    PushCallStack("MakeNormalUniformSpectrum");
#endif
    typedef Complex<R> C;
    if( A.Height() != A.Width() )
        throw std::logic_error("Cannot make a non-square matrix normal");

    // Sample the diagonal matrix D from the ball B_radius(center)
    // and then rotate it with a random Householder similarity transformation:
    //
    //  (I-2uu^H) D (I-2uu^H)^H = D - 2(u (conj(D) u)^H + (D u) u^H) + 
    //                                (4 u^H D u) u u^H
    //

    // Form d and D
    const int n = A.Height();
    std::vector<C> d( n );
    for( int j=0; j<n; ++j )
        d[j] = center + radius*SampleUnitBall<C>();
    Diagonal( d, A );

    // Form u 
    Matrix<C> u( n, 1 );
    MakeUniform( u );
    const R origNorm = Nrm2( u );
    Scal( (C)1/origNorm, u );

    // Form v := D u
    Matrix<C> v( n, 1 );
    for( int i=0; i<n; ++i )
        v.Set( i, 0, d[i]*u.Get(i,0) );

    // Form w := conj(D) u
    Matrix<C> w( n, 1 );
    for( int i=0; i<n; ++i )
        w.Set( i, 0, Conj(d[i])*u.Get(i,0) );

    // Update A := A - 2(u w^H + v u^H)
    Ger( (C)-2, u, w, A );
    Ger( (C)-2, v, u, A );

    // Form \gamma := 4 u^H (D u) = 4 (u,Du)
    const C gamma = 4*Dot(u,v);

    // Update A := A + gamma u u^H
    Ger( gamma, u, u, A );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R,Distribution U,Distribution V>
inline void
MakeNormalUniformSpectrum
( DistMatrix<Complex<R>,U,V>& A, Complex<R> center, R radius )
{
#ifndef RELEASE
    PushCallStack("MakeNormalUniformSpectrum");
#endif
    typedef Complex<R> C;
    if( A.Height() != A.Width() )
        throw std::logic_error("Cannot make a non-square matrix normal");
    const Grid& grid = A.Grid();
    const bool standardDist = ( U == MC && V == MR );

    // Sample the diagonal matrix D from the ball B_radius(center)
    // and then rotate it with a random Householder similarity transformation:
    //
    //  (I-2uu^H) D (I-2uu^H)^H = D - 2(u (Conj(D) u)^H + (D u) u^H) + 
    //                                (4 u^H D u) u u^H
    //

    // Form d and D
    const int n = A.Height();
    std::vector<C> d( n );
    for( int j=0; j<n; ++j )
        d[j] = center + radius*SampleUnitBall<C>();
    DistMatrix<C,MC,MR> ABackup( grid );
    if( standardDist )
        Diagonal( d, A );
    else
    {
        ABackup.AlignWith( A );
        Diagonal( d, ABackup );
    }

    // Form u 
    DistMatrix<C,MC,MR> u( grid );
    if( standardDist )
        u.AlignWith( A );
    else
        u.AlignWith( ABackup );
    Uniform( n, 1, u );
    const R origNorm = Nrm2( u );
    Scal( (C)1/origNorm, u );

    // Form v := D u
    DistMatrix<C,MC,MR> v( grid );
    if( standardDist )
        v.AlignWith( A );
    else
        v.AlignWith( ABackup );
    v.ResizeTo( n, 1 );
    if( v.LocalWidth() == 1 )
    {
        const int colShift = v.ColShift();
        const int colStride = v.ColStride();
        const int localHeight = v.LocalHeight();
        for( int iLocal=0; iLocal<localHeight; ++iLocal )
        {
            const int i = colShift + iLocal*colStride;
            v.SetLocalEntry( iLocal, 0, d[i]*u.GetLocalEntry(iLocal,0) );
        }
    }

    // Form w := Conj(D) u
    DistMatrix<C,MC,MR> w( grid );
    if( standardDist )
        w.AlignWith( A );
    else
        w.AlignWith( ABackup );
    w.ResizeTo( n, 1 );
    if( w.LocalWidth() == 1 )
    {
        const int colShift = w.ColShift();
        const int colStride = w.ColStride();
        const int localHeight = w.LocalHeight();
        for( int iLocal=0; iLocal<localHeight; ++iLocal )
        {
            const int i = colShift + iLocal*colStride;
            w.SetLocalEntry( iLocal, 0, Conj(d[i])*u.GetLocalEntry(iLocal,0) );
        }
    }

    // Update A := A - 2(u w^H + v u^H)
    if( standardDist )
    {
        Ger( (C)-2, u, w, A );
        Ger( (C)-2, v, u, A );
    }
    else
    {
        Ger( (C)-2, u, w, ABackup );
        Ger( (C)-2, v, u, ABackup );
    }

    // Form \gamma := 4 u^H (D u) = 4 (u,Du)
    const C gamma = 4*Dot(u,v);

    // Update A := A + gamma u u^H
    if( standardDist )
        Ger( gamma, u, u, A );
    else
        Ger( gamma, u, u, ABackup );

    // Copy the result into the correct distribution
    if( !standardDist )
        A = ABackup;
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
