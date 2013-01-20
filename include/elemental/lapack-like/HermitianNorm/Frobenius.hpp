/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_HERMITIANNORM_FROBENIUS_HPP
#define LAPACK_HERMITIANNORM_FROBENIUS_HPP

namespace elem {
namespace internal {

template<typename F>
inline typename Base<F>::type 
HermitianFrobeniusNorm( UpperOrLower uplo, const Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("internal::HermitianFrobeniusNorm");
#endif
    typedef typename Base<F>::type R;

    if( A.Height() != A.Width() )
        throw std::logic_error("Hermitian matrices must be square.");

    R scale = 0;
    R scaledSquare = 1;
    const int height = A.Height();
    const int width = A.Width();
    if( uplo == UPPER )
    {
        for( int j=0; j<width; ++j )
        {
            for( int i=0; i<j; ++i )
            {
                const R alphaAbs = Abs(A.Get(i,j));
                if( alphaAbs != 0 )
                {
                    if( alphaAbs <= scale )
                    {
                        const R relScale = alphaAbs/scale;
                        scaledSquare += 2*relScale*relScale;
                    }
                    else
                    {
                        const R relScale = scale/alphaAbs;
                        scaledSquare = scaledSquare*relScale*relScale + 2;
                        scale = alphaAbs;
                    }
                }
            }
            const R alphaAbs = Abs(A.Get(j,j));
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
    else
    {
        for( int j=0; j<width; ++j )
        {
            for( int i=j+1; i<height; ++i )
            {
                const R alphaAbs = Abs(A.Get(i,j));
                if( alphaAbs != 0 )
                {
                    if( alphaAbs <= scale )
                    {
                        const R relScale = alphaAbs/scale;
                        scaledSquare += 2*relScale*relScale;
                    }
                    else
                    {
                        const R relScale = scale/alphaAbs;
                        scaledSquare = scaledSquare*relScale*relScale + 2;
                        scale = alphaAbs;
                    }
                }
            }
            const R alphaAbs = Abs(A.Get(j,j));
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

template<typename F>
inline typename Base<F>::type
HermitianFrobeniusNorm
( UpperOrLower uplo, const DistMatrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("internal::HermitianFrobeniusNorm");
#endif
    typedef typename Base<F>::type R;

    if( A.Height() != A.Width() )
        throw std::logic_error("Hermitian matrices must be square.");

    const int r = A.Grid().Height();
    const int c = A.Grid().Width();
    const int colShift = A.ColShift();
    const int rowShift = A.RowShift();

    R localScale = 0;
    R localScaledSquare = 1;
    const int localWidth = A.LocalWidth();
    if( uplo == UPPER )
    {
        for( int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            int j = rowShift + jLocal*c;
            int numUpperRows = LocalLength(j+1,colShift,r);
            for( int iLocal=0; iLocal<numUpperRows; ++iLocal )
            {
                int i = colShift + iLocal*r;
                const R alphaAbs = Abs(A.GetLocal(iLocal,jLocal));
                if( alphaAbs != 0 )
                {
                    if( alphaAbs <= localScale )
                    {
                        const R relScale = alphaAbs/localScale;
                        if( i != j )
                            localScaledSquare += 2*relScale*relScale;
                        else
                            localScaledSquare += relScale*relScale;
                    }
                    else
                    {
                        const R relScale = localScale/alphaAbs;
                        if( i != j )
                            localScaledSquare = 
                                localScaledSquare*relScale*relScale + 2;
                        else
                            localScaledSquare = 
                                localScaledSquare*relScale*relScale + 1;
                        localScale = alphaAbs;
                    }
                }
            }
        }
    }
    else
    {
        for( int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            int j = rowShift + jLocal*c;
            int numStrictlyUpperRows = LocalLength(j,colShift,r);
            for( int iLocal=numStrictlyUpperRows; 
                 iLocal<A.LocalHeight(); ++iLocal )
            {
                int i = colShift + iLocal*r;
                const R alphaAbs = Abs(A.GetLocal(iLocal,jLocal));
                if( alphaAbs != 0 )
                {
                    if( alphaAbs <= localScale )
                    {
                        const R relScale = alphaAbs/localScale;
                        if( i != j )
                            localScaledSquare += 2*relScale*relScale;
                        else
                            localScaledSquare += relScale*relScale;
                    }
                    else
                    {
                        const R relScale = localScale/alphaAbs;
                        if( i != j )
                            localScaledSquare = 
                                localScaledSquare*relScale*relScale + 2;
                        else
                            localScaledSquare =
                                localScaledSquare*relScale*relScale + 1; 
                        localScale = alphaAbs;
                    }
                }
            }
        }
    }

    // Find the maximum relative scale
    R scale;
    mpi::AllReduce( &localScale, &scale, 1, mpi::MAX, A.Grid().VCComm() );

    R norm = 0;
    if( scale != 0 )
    {
        // Equilibrate our local scaled sum to the maximum scale
        R relScale = localScale/scale;
        localScaledSquare *= relScale*relScale;

        // The scaled square is now simply the sum of the local contributions
        R scaledSquare;
        mpi::AllReduce
        ( &localScaledSquare, &scaledSquare, 1, mpi::SUM, A.Grid().VCComm() );

        norm = scale*Sqrt(scaledSquare);
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

} // namespace internal
} // namespace elem

#endif // ifndef LAPACK_HERMITIANNORM_FROBENIUS_HPP
