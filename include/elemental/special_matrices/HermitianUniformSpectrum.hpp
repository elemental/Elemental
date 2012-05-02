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

template<typename F>
inline void
HermitianUniformSpectrum
( int n, Matrix<F>& A, 
  typename Base<F>::type lower, typename Base<F>::type upper )
{
#ifndef RELEASE
    PushCallStack("HermitianUniformSpectrum");
#endif
    A.ResizeTo( n, n );
    MakeHermitianUniformSpectrum( A, lower, upper );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F,Distribution U,Distribution V>
inline void
HermitianUniformSpectrum
( int n, DistMatrix<F,U,V>& A, 
  typename Base<F>::type lower, typename Base<F>::type upper )
{
#ifndef RELEASE
    PushCallStack("HermitianUniformSpectrum");
#endif
    A.ResizeTo( n, n );
    MakeHermitianUniformSpectrum( A, lower, upper );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
MakeHermitianUniformSpectrum
( Matrix<F>& A, typename Base<F>::type lower, typename Base<F>::type upper )
{
#ifndef RELEASE
    PushCallStack("MakeHermitianUniformSpectrum");
#endif
    if( A.Height() != A.Width() )
        throw std::logic_error("Cannot make a non-square matrix Hermitian");
    typedef typename Base<F>::type R;
    const bool isComplex = IsComplex<F>::val;

    // Sample the diagonal matrix D from the half-open interval (lower,upper]
    // and then rotate it with a random Householder similarity transformation:
    //
    //  (I-2uu^H) D (I-2uu^H)^H = D - 2(u (D u)^H + (D u) u^H) + 
    //                                (4 u^H D u) u u^H
    //

    // Form d and D
    const int n = A.Height();
    std::vector<F> d( n );
    for( int j=0; j<n; ++j )
        d[j] = lower + (upper-lower)*plcg::SerialUniform<R>();
    Diagonal( d, A );

    // Form u 
    Matrix<F> u( n, 1 );
    MakeUniform( u );
    const R origNorm = Nrm2( u );
    Scal( (F)1/origNorm, u );

    // Form v := D u
    Matrix<F> v( n, 1 );
    for( int i=0; i<n; ++i )
        v.Set( i, 0, d[i]*u.Get(i,0) );

    // Update A := A - 2(u v^H + v u^H)
    Ger( (F)-2, u, v, A );
    Ger( (F)-2, v, u, A );

    // Form \gamma := 4 u^H (D u) = 4 (u,Du)
    const F gamma = 4*Dot(u,v);

    // Update A := A + gamma u u^H
    Ger( gamma, u, u, A );

    // Force the diagonal to be real
    if( isComplex )
    {
        const int n = A.Height();
        for( int j=0; j<n; ++j )
            A.SetImag( j, j, (R)0 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F,Distribution U,Distribution V>
inline void
MakeHermitianUniformSpectrum
( DistMatrix<F,U,V>& A, 
  typename Base<F>::type lower, typename Base<F>::type upper )
{
#ifndef RELEASE
    PushCallStack("MakeHermitianUniformSpectrum");
#endif
    if( A.Height() != A.Width() )
        throw std::logic_error("Cannot make a non-square matrix Hermitian");
    const Grid& grid = A.Grid();
    typedef typename Base<F>::type R;
    const bool isComplex = IsComplex<F>::val;
    const bool standardDist = ( U == MC && V == MR );

    // Sample the diagonal matrix D from the half-open interval (lower,upper]
    // and then rotate it with a random Householder similarity transformation:
    //
    //  (I-2uu^H) D (I-2uu^H)^H = D - 2(u (D u)^H + (D u) u^H) + 
    //                                (4 u^H D u) u u^H
    //

    // Form d and D
    const int n = A.Height();
    std::vector<F> d( n );
    for( int j=0; j<n; ++j )
        d[j] = lower + (upper-lower)*plcg::SerialUniform<R>();
    DistMatrix<F,MC,MR> ABackup( grid );
    if( standardDist )
        Diagonal( d, A );
    else
    {
        ABackup.AlignWith( A );
        Diagonal( d, ABackup );
    }

    // Form u 
    DistMatrix<F,MC,MR> u( grid );
    if( standardDist )
        u.AlignWith( A );
    else
        u.AlignWith( ABackup );
    Uniform( n, 1, u );
    const R origNorm = Nrm2( u );
    Scal( (F)1/origNorm, u );

    // Form v := D u
    DistMatrix<F,MC,MR> v( grid );
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

    // Update A := A - 2(u v^H + v u^H)
    if( standardDist )
    {
        Ger( (F)-2, u, v, A );
        Ger( (F)-2, v, u, A );
    }
    else
    {
        Ger( (F)-2, u, v, ABackup );
        Ger( (F)-2, v, u, ABackup );
    }

    // Form \gamma := 4 u^H (D u) = 4 (u,Du)
    const F gamma = 4*Dot(u,v);

    // Update A := A + gamma u u^H
    if( standardDist )
        Ger( gamma, u, u, A );
    else
        Ger( gamma, u, u, ABackup );

    // Copy the result into the correct distribution
    if( !standardDist )
        A = ABackup;

    // Force the diagonal to be real-valued
    if( isComplex )
    {
        const int localHeight = A.LocalHeight();
        const int localWidth = A.LocalWidth();
        const int colShift = A.ColShift();
        const int rowShift = A.RowShift();
        const int colStride = A.ColStride();
        const int rowStride = A.RowStride();
        for( int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            const int j = rowShift + jLocal*rowStride;
            for( int iLocal=0; iLocal<localHeight; ++iLocal )
            {
                const int i = colShift + iLocal*colStride;
                if( i == j )
                    A.SetImagLocalEntry( iLocal, jLocal, (R)0 );
            }
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
