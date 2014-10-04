/*
   Copyright (c) 2009-2014, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, Stanford University, and the
   Georgia Insitute of Technology.
   All rights reserved.
 
   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_SPARSEDIRECT_NUMERIC_LDL_FRONTBLOCK_HPP
#define EL_SPARSEDIRECT_NUMERIC_LDL_FRONTBLOCK_HPP

namespace El {

template<typename F> 
inline void FrontBlockLDL
( Matrix<F>& AL, Matrix<F>& ABR, bool conjugate, bool intraPiv )
{
    DEBUG_ONLY(CallStackEntry cse("FrontBlockLDL"))
    Matrix<F> ATL, ABL;
    PartitionDown( AL, ATL, ABL, AL.Width() );
    
    // Make a copy of the original contents of ABL
    Matrix<F> BBL( ABL );

    if( intraPiv )
    {
        Matrix<Int> p;
        Matrix<F> dSub;
        // TODO: Expose the pivot type as an option?
        LDL( ATL, dSub, p, conjugate, BUNCH_KAUFMAN_A );

        // Solve against ABL and update ABR
        // NOTE: This does not exploit symmetry
        ldl::SolveAfter( ATL, dSub, p, ABL, conjugate );
        const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );
        Gemm( NORMAL, orientation, F(-1), ABL, BBL, F(1), ABR );

        // Copy the original contents of ABL back
        ABL = BBL;

        // Finish inverting ATL
        TriangularInverse( LOWER, UNIT, ATL );
        Trdtrmm( LOWER, ATL, dSub, conjugate );
        // TODO: SymmetricPermutation
        MakeSymmetric( LOWER, ATL, conjugate );
        PermuteRows( ATL, p );
        PermuteCols( ATL, p );
    }
    else
    {
        // Call the standard routine
        FrontLDL( AL, ABR, conjugate );

        // Copy the original contents of ABL back
        ABL = BBL;

        // Finish inverting ATL
        TriangularInverse( LOWER, UNIT, ATL );
        Trdtrmm( LOWER, ATL, conjugate );
        MakeSymmetric( LOWER, ATL, conjugate );
    }
}

template<typename F>
inline void FrontBlockLDL
( DistMatrix<F>& AL, DistMatrix<F>& ABR, bool conjugate, bool intraPiv )
{
    DEBUG_ONLY(CallStackEntry cse("FrontBlockLDL"))
    const Grid& g = AL.Grid();
    DistMatrix<F> ATL(g), ABL(g);
    PartitionDown( AL, ATL, ABL, AL.Width() );

    // Make a copy of the original contents of ABL
    DistMatrix<F> BBL( ABL );

    if( intraPiv )
    {
        DistMatrix<Int,VC,STAR> p( ATL.Grid() );
        DistMatrix<F,MD,STAR> dSub( ATL.Grid() );
        // TODO: Expose the pivot type as an option?
        LDL( ATL, dSub, p, conjugate, BUNCH_KAUFMAN_A );

        // Solve against ABL and update ABR
        // NOTE: This update does not exploit symmetry
        ldl::SolveAfter( ATL, dSub, p, ABL, conjugate );
        const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );
        Gemm( NORMAL, orientation, F(-1), ABL, BBL, F(1), ABR );

        // Copy the original contents of ABL back
        ABL = BBL;

        // Finish inverting ATL
        TriangularInverse( LOWER, UNIT, ATL );
        Trdtrmm( LOWER, ATL, dSub, conjugate );
        ApplyInverseSymmetricPivots( LOWER, ATL, p, conjugate );
    }
    else
    {
        // Call the standard routine
        FrontLDL( AL, ABR, conjugate );

        // Copy the original contents of ABL back
        ABL = BBL;

        // Finish inverting ATL
        TriangularInverse( LOWER, UNIT, ATL );
        Trdtrmm( LOWER, ATL, conjugate );
    }
    MakeSymmetric( LOWER, ATL, conjugate );
}

} // namespace El

#endif // ifndef EL_SPARSEDIRECT_NUMERIC_LDL_FRONTBLOCK_HPP
