/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_NORM_FROBENIUS_HPP
#define EL_NORM_FROBENIUS_HPP

namespace El {

template<typename F> 
Base<F> FrobeniusNorm( const Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("FrobeniusNorm"))
    typedef Base<F> R;
    R scale = 0;
    R scaledSquare = 1;
    const Int width = A.Width();
    const Int height = A.Height();
    for( Int j=0; j<width; ++j )
    {
        for( Int i=0; i<height; ++i )
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
    return scale*Sqrt(scaledSquare);
}

template<typename F>
Base<F> HermitianFrobeniusNorm( UpperOrLower uplo, const Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianFrobeniusNorm"))
    if( A.Height() != A.Width() )
        LogicError("Hermitian matrices must be square.");

    typedef Base<F> R;
    R scale = 0;
    R scaledSquare = 1;
    const Int height = A.Height();
    const Int width = A.Width();
    if( uplo == UPPER )
    {
        for( Int j=0; j<width; ++j )
        {
            for( Int i=0; i<j; ++i )
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
        for( Int j=0; j<width; ++j )
        {
            for( Int i=j+1; i<height; ++i )
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
    return scale*Sqrt(scaledSquare);
}

template<typename F>
Base<F> SymmetricFrobeniusNorm( UpperOrLower uplo, const Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("SymmetricFrobeniusNorm"))
    return HermitianFrobeniusNorm( uplo, A );
}

template<typename F> 
Base<F> FrobeniusNorm( const AbstractDistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("FrobeniusNorm"))
    typedef Base<F> Real;
    Real norm;
    if( A.Participating() )
    {
        Real locScale=0, locScaledSquare=1;
        const Int localHeight = A.LocalHeight();
        const Int localWidth = A.LocalWidth();
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
            {
                const Real alphaAbs = Abs(A.GetLocal(iLoc,jLoc));
                if( alphaAbs != 0 )
                {
                    if( alphaAbs <= locScale )
                    {
                        const Real relScale = alphaAbs/locScale;
                        locScaledSquare += relScale*relScale;
                    }
                    else
                    {
                        const Real relScale = locScale/alphaAbs;
                        locScaledSquare = locScaledSquare*relScale*relScale + 1;
                        locScale = alphaAbs; 
                    }
                }
            }
        }

        // Find the maximum relative scale
        mpi::Comm comm = A.DistComm();
        const Real scale = mpi::AllReduce( locScale, mpi::MAX, comm );

        norm = 0;
        if( scale != 0 )
        {
            // Equilibrate our local scaled sum to the maximum scale
            Real relScale = locScale/scale;
            locScaledSquare *= relScale*relScale;

            // The scaled square is now the sum of the local contributions
            const Real scaledSquare = mpi::AllReduce( locScaledSquare, comm );
            norm = scale*Sqrt(scaledSquare);
        }
    }
    mpi::Broadcast( norm, A.Root(), A.CrossComm() );
    return norm;
}

template<typename F>
Base<F> HermitianFrobeniusNorm
( UpperOrLower uplo, const AbstractDistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianFrobeniusNorm"))
    if( A.Height() != A.Width() )
        LogicError("Hermitian matrices must be square.");

    typedef Base<F> Real;
    Real norm;
    if( A.Participating() )
    {
        Real localScale = 0;
        Real localScaledSquare = 1;
        const Int localWidth = A.LocalWidth();
        if( uplo == UPPER )
        {
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const Int j = A.GlobalCol(jLoc);
                const Int numUpperRows = A.LocalRowOffset(j+1);
                for( Int iLoc=0; iLoc<numUpperRows; ++iLoc )
                {
                    const Int i = A.GlobalRow(iLoc);
                    const Real alphaAbs = Abs(A.GetLocal(iLoc,jLoc));
                    if( alphaAbs != 0 )
                    {
                        if( alphaAbs <= localScale )
                        {
                            const Real relScale = alphaAbs/localScale;
                            if( i != j )
                                localScaledSquare += 2*relScale*relScale;
                            else
                                localScaledSquare += relScale*relScale;
                        }
                        else
                        {
                            const Real relScale = localScale/alphaAbs;
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
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const Int j = A.GlobalCol(jLoc);
                const Int numStrictlyUpperRows = A.LocalRowOffset(j);
                for( Int iLoc=numStrictlyUpperRows;
                     iLoc<A.LocalHeight(); ++iLoc )
                {
                    const Int i = A.GlobalRow(iLoc);
                    const Real alphaAbs = Abs(A.GetLocal(iLoc,jLoc));
                    if( alphaAbs != 0 )
                    {
                        if( alphaAbs <= localScale )
                        {
                            const Real relScale = alphaAbs/localScale;
                            if( i != j )
                                localScaledSquare += 2*relScale*relScale;
                            else
                                localScaledSquare += relScale*relScale;
                        }
                        else
                        {
                            const Real relScale = localScale/alphaAbs;
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
        const Real scale = mpi::AllReduce( localScale, mpi::MAX, A.DistComm() );

        norm = 0;
        if( scale != 0 )
        {
            // Equilibrate our local scaled sum to the maximum scale
            Real relScale = localScale/scale;
            localScaledSquare *= relScale*relScale;

            // The scaled square is now the sum of the local contributions
            const Real scaledSquare = 
                mpi::AllReduce( localScaledSquare, A.Grid().VCComm() );
            norm = scale*Sqrt(scaledSquare);
        }
    }
    mpi::Broadcast( norm, A.Root(), A.CrossComm() );
    return norm;
}

template<typename F>
Base<F> SymmetricFrobeniusNorm
( UpperOrLower uplo, const AbstractDistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("SymmetricFrobeniusNorm"))
    return HermitianFrobeniusNorm( uplo, A );
}

} // namespace El

#endif // ifndef EL_NORM_FROBENIUS_HPP
