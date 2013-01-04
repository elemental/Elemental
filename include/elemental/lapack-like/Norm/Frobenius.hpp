/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace elem {
namespace internal {

template<typename F> 
inline typename Base<F>::type
FrobeniusNorm( const Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("internal::FrobeniusNorm");
#endif
    typedef typename Base<F>::type R;

    R scale = 0;
    R scaledSquare = 1;

    const int width = A.Width();
    const int height = A.Height();
    for( int j=0; j<width; ++j )
    {
        for( int i=0; i<height; ++i )
        {
            const R alphaAbs = Abs(A.Get(i,j));
            if( alphaAbs != 0 )
            {
                if( alphaAbs <= scale )
                {
                    const R relScale = alphaAbs/scale;
                    scaledSquare += relScale*relScale;
                }
                else
                {
                    const R relScale = scale/alphaAbs;
                    scaledSquare = scaledSquare*relScale*relScale + 1;
                    scale = alphaAbs;
                }
            }
        }
    }
    const R norm = scale*Sqrt(scaledSquare);
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

template<typename F,Distribution U,Distribution V> 
inline typename Base<F>::type
FrobeniusNorm( const DistMatrix<F,U,V>& A )
{
#ifndef RELEASE
    PushCallStack("internal::FrobeniusNorm");
#endif
    typedef typename Base<F>::type R;
    mpi::Comm comm = NormComm( A );

    R localScale = 0;
    R localScaledSquare = 1;
    const int localHeight = A.LocalHeight();
    const int localWidth = A.LocalWidth();
    for( int jLocal=0; jLocal<localWidth; ++jLocal )
    {
        for( int iLocal=0; iLocal<localHeight; ++iLocal )
        {
            const R alphaAbs = Abs(A.GetLocal(iLocal,jLocal));
            if( alphaAbs != 0 )
            {
                if( alphaAbs <= localScale )
                {
                    const R relScale = alphaAbs/localScale;
                    localScaledSquare += relScale*relScale;
                }
                else
                {
                    const R relScale = localScale/alphaAbs;
                    localScaledSquare = localScaledSquare*relScale*relScale + 1;
                    localScale = alphaAbs; 
                }
            }
        }
    }

    // Find the maximum relative scale
    R scale;
    mpi::AllReduce( &localScale, &scale, 1, mpi::MAX, comm );

    R norm = 0;
    if( scale != 0 )
    {
        // Equilibrate our local scaled sum to the maximum scale
        R relScale = localScale/scale;
        localScaledSquare *= relScale*relScale;

        // The scaled square is now simply the sum of the local contributions
        R scaledSquare;
        mpi::AllReduce( &localScaledSquare, &scaledSquare, 1, mpi::SUM, comm );

        norm = scale*Sqrt(scaledSquare);
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

} // namespace internal
} // namespace elem
