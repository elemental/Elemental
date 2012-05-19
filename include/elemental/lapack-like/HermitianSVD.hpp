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

#ifndef WITHOUT_PMRRR

namespace elem {

//
// Compute the SVD of a Hermitian matrix through its EVD.
//

template<typename F>
inline void HermitianSVD
( UpperOrLower uplo, 
  DistMatrix<F,MC,MR>& A, DistMatrix<typename Base<F>::type,VR,STAR>& s, 
  DistMatrix<F,MC,MR>& U, DistMatrix<F,MC,MR>& V )
{
#ifndef RELEASE
    PushCallStack("HermitianSVD");
#endif
    typedef typename Base<F>::type R;

    // Grab an eigenvalue decomposition of A
    HermitianEig( uplo, A, s, V ); 
    
    // Redistribute the singular values into an [MR,* ] distribution
    const Grid& grid = A.Grid();
    DistMatrix<R,MR,STAR> s_MR_STAR( grid );
    s_MR_STAR.AlignWith( V );
    s_MR_STAR = s;

    // Set the singular values to the absolute value of the eigenvalues
    const int numLocalVals = s.LocalHeight();
    for( int iLocal=0; iLocal<numLocalVals; ++iLocal )
    {
        const R sigma = s.GetLocalEntry(iLocal,0);
        s.SetLocalEntry(iLocal,0,Abs(sigma));
    }

    // Copy V into U (flipping the sign as necessary)
    U.AlignWith( V );
    U.ResizeTo( V.Height(), V.Width() );
    const int localHeight = V.LocalHeight();
    const int localWidth = V.LocalWidth();
    for( int jLocal=0; jLocal<localWidth; ++jLocal )
    {
        const R sigma = s_MR_STAR.GetLocalEntry( jLocal, 0 );
        F* UCol = U.LocalBuffer( 0, jLocal );
        const F* VCol = V.LockedLocalBuffer( 0, jLocal );
        if( sigma >= 0 )
            for( int iLocal=0; iLocal<localHeight; ++iLocal )
                UCol[iLocal] = VCol[iLocal];
        else
            for( int iLocal=0; iLocal<localHeight; ++iLocal )
                UCol[iLocal] = -VCol[iLocal];
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void HermitianSingularValues
( UpperOrLower uplo, 
  DistMatrix<F,MC,MR>& A, DistMatrix<typename Base<F>::type,VR,STAR>& s )
{
#ifndef RELEASE
    PushCallStack("HermitianSingularValues");
#endif
    typedef typename Base<F>::type R;

    // Grab an eigenvalue decomposition of A
    HermitianEig( uplo, A, s ); 
    
    // Replace the eigenvalues with their absolute values
    const int numLocalVals = s.LocalHeight();
    for( int iLocal=0; iLocal<numLocalVals; ++iLocal )
    {
        const R sigma = s.GetLocalEntry(iLocal,0);
        s.SetLocalEntry(iLocal,0,Abs(sigma));
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem

#endif // WITHOUT_PMRRR
