/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

#ifndef WITHOUT_PMRRR

namespace elem {

//
// Compute the SVD of a Hermitian matrix through its EVD.
//

template<typename F>
inline void HermitianSVD
( UpperOrLower uplo, 
  DistMatrix<F>& A, DistMatrix<typename Base<F>::type,VR,STAR>& s, 
  DistMatrix<F>& U, DistMatrix<F>& V )
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
        const R sigma = s.GetLocal(iLocal,0);
        s.SetLocal(iLocal,0,Abs(sigma));
    }

    // Copy V into U (flipping the sign as necessary)
    U.AlignWith( V );
    U.ResizeTo( V.Height(), V.Width() );
    const int localHeight = V.LocalHeight();
    const int localWidth = V.LocalWidth();
    for( int jLocal=0; jLocal<localWidth; ++jLocal )
    {
        const R sigma = s_MR_STAR.GetLocal( jLocal, 0 );
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
  DistMatrix<F>& A, DistMatrix<typename Base<F>::type,VR,STAR>& s )
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
        const R sigma = s.GetLocal(iLocal,0);
        s.SetLocal(iLocal,0,Abs(sigma));
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem

#endif // WITHOUT_PMRRR
