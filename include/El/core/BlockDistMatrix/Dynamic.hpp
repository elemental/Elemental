/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_BLOCKDISTMATRIX_DYNAMIC_DECL_HPP
#define EL_BLOCKDISTMATRIX_DYNAMIC_DECL_HPP

namespace El {

// This is designed to prevent the combinatorial explosion of the number of
// external interface routines (e.g., for Python) needed for routines such as 
// 'Copy', which allow for different distributions for the input and output 
// matrices. The idea is to use the U and V Dist member variables to branch
// to various dynamic_cast's.
template<typename T> 
struct DynamicBlockDistMatrix
{
    Dist U=MC, V=MR;
    AbstractBlockDistMatrix<T>* ABDM=nullptr;

    ~DynamicBlockDistMatrix() { delete ABDM; }

    void SetDistribution( Dist UNew, Dist VNew, const Grid& grid )
    {
        if( ABDM != nullptr )
        {
            if( U == UNew && V == VNew && ABDM->Grid() == grid )
                return;
            else
                delete ABDM;
        }

        U = UNew;
        V = VNew;
        if( UNew == CIRC && VNew == CIRC )
            ABDM = new BlockDistMatrix<T,CIRC,CIRC>(grid);
        else if( UNew == MC && VNew == MR )
            ABDM = new BlockDistMatrix<T,MC,MR>(grid);
        else if( UNew == MC && VNew == STAR )
            ABDM = new BlockDistMatrix<T,MC,STAR>(grid);
        else if( UNew == MD && VNew == STAR )
            ABDM = new BlockDistMatrix<T,MD,STAR>(grid);
        else if( UNew == MR && VNew == MC )
            ABDM = new BlockDistMatrix<T,MR,MC>(grid);
        else if( UNew == MR && VNew == STAR )
            ABDM = new BlockDistMatrix<T,MR,STAR>(grid);
        else if( UNew == STAR && VNew == MC )
            ABDM = new BlockDistMatrix<T,STAR,MC>(grid);
        else if( UNew == STAR && VNew == MD )
            ABDM = new BlockDistMatrix<T,STAR,MD>(grid);
        else if( UNew == STAR && VNew == MR )
            ABDM = new BlockDistMatrix<T,STAR,MR>(grid);
        else if( UNew == STAR && VNew == STAR )
            ABDM = new BlockDistMatrix<T,STAR,STAR>(grid);
        else if( UNew == STAR && VNew == VC )
            ABDM = new BlockDistMatrix<T,STAR,VC>(grid);
        else if( UNew == STAR && VNew == VR )
            ABDM = new BlockDistMatrix<T,STAR,VR>(grid);
        else if( UNew == VC && VNew == STAR )
            ABDM = new BlockDistMatrix<T,VC,STAR>(grid);
        else if( UNew == VR && VNew == STAR )
            ABDM = new BlockDistMatrix<T,VR,STAR>(grid);
    }
};

} // namespace El

#endif // ifndef EL_BLOCKDISTMATRIX_DYNAMIC_DECL_HPP
