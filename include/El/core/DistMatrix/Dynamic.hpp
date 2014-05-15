/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_DISTMATRIX_DYNAMIC_DECL_HPP
#define EL_DISTMATRIX_DYNAMIC_DECL_HPP

namespace El {

// This is designed to prevent the combinatorial explosion of the number of
// external interface routines (e.g., for Python) needed for routines such as 
// 'Copy', which allow for different distributions for the input and output 
// matrices. The idea is to use the U and V Dist member variables to branch
// to various dynamic_cast's.
template<typename T> 
struct DynamicDistMatrix
{
    Dist U=MC, V=MR;
    AbstractDistMatrix<T>* ADM=nullptr;

    ~DynamicDistMatrix() { delete ADM; }

    void SetDistribution( Dist UNew, Dist VNew, const Grid& grid )
    {
        if( ADM != nullptr )
        {
            if( U == UNew && V == VNew && ADM->Grid() == grid )
                return;
            else
                delete ADM;
        }

        U = UNew;
        V = VNew;
        if( UNew == CIRC && VNew == CIRC )
            ADM = new DistMatrix<T,CIRC,CIRC>(grid);
        else if( UNew == MC && VNew == MR )
            ADM = new DistMatrix<T,MC,MR>(grid);
        else if( UNew == MC && VNew == STAR )
            ADM = new DistMatrix<T,MC,STAR>(grid);
        else if( UNew == MD && VNew == STAR )
            ADM = new DistMatrix<T,MD,STAR>(grid);
        else if( UNew == MR && VNew == MC )
            ADM = new DistMatrix<T,MR,MC>(grid);
        else if( UNew == MR && VNew == STAR )
            ADM = new DistMatrix<T,MR,STAR>(grid);
        else if( UNew == STAR && VNew == MC )
            ADM = new DistMatrix<T,STAR,MC>(grid);
        else if( UNew == STAR && VNew == MD )
            ADM = new DistMatrix<T,STAR,MD>(grid);
        else if( UNew == STAR && VNew == MR )
            ADM = new DistMatrix<T,STAR,MR>(grid);
        else if( UNew == STAR && VNew == STAR )
            ADM = new DistMatrix<T,STAR,STAR>(grid);
        else if( UNew == STAR && VNew == VC )
            ADM = new DistMatrix<T,STAR,VC>(grid);
        else if( UNew == STAR && VNew == VR )
            ADM = new DistMatrix<T,STAR,VR>(grid);
        else if( UNew == VC && VNew == STAR )
            ADM = new DistMatrix<T,VC,STAR>(grid);
        else if( UNew == VR && VNew == STAR )
            ADM = new DistMatrix<T,VR,STAR>(grid);
    }
};

} // namespace El

#endif // ifndef EL_DISTMATRIX_DYNAMIC_DECL_HPP
