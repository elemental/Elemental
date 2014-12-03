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
#ifndef EL_OPTIMIZATION_REGLDL_PROCESSLOCALTREE_HPP
#define EL_OPTIMIZATION_REGLDL_PROCESSLOCALTREE_HPP

namespace El {
namespace reg_ldl {

template<typename F> 
inline void ProcessLocalTree
( DistSymmInfo& info, DistSymmFrontTree<F>& L, Base<F> pivTol,
  const DistNodalMultiVec<Base<F>>& regCand, DistNodalMultiVec<Base<F>>& reg )
{
    DEBUG_ONLY(CallStackEntry cse("reg_ldl::ProcessLocalTree"))

    const int numLocalNodes = info.localNodes.size();
    for( int s=0; s<numLocalNodes; ++s )
    {
        SymmNodeInfo& node = info.localNodes[s];
        const int updateSize = node.lowerStruct.size();
        SymmFront<F>& front = L.localFronts[s];
        Matrix<F>& frontL = front.frontL;
        Matrix<F>& frontBR = front.work;
        frontBR.Empty();
        DEBUG_ONLY(
            if( frontL.Height() != node.size+updateSize ||
                frontL.Width() != node.size )
                LogicError("Front was not the proper size");
        )

        // Add updates from children (if they exist)
        // TODO: Generalize to an arbitrary number of children
        Zeros( frontBR, updateSize, updateSize );
        const int numChildren = node.children.size();
        if( numChildren == 2 )
        {
            const int leftInd = node.children[0];
            const int rightInd = node.children[1];
            Matrix<F>& leftUpdate = L.localFronts[leftInd].work;
            Matrix<F>& rightUpdate = L.localFronts[rightInd].work;

            // Add the left child's update matrix
            const int leftUpdateSize = leftUpdate.Height();
            for( int jChild=0; jChild<leftUpdateSize; ++jChild )
            {
                const int jFront = node.leftRelInds[jChild];
                for( int iChild=jChild; iChild<leftUpdateSize; ++iChild )
                {
                    const int iFront = node.leftRelInds[iChild];
                    const F value = leftUpdate.Get(iChild,jChild);
                    DEBUG_ONLY(
                        if( iFront < jFront )
                            LogicError("Tried to update upper triangle");
                    )
                    if( jFront < node.size )
                        frontL.Update( iFront, jFront, value );
                    else if( iFront >= node.size )
                        frontBR.Update
                        ( iFront-node.size, jFront-node.size, value );
                }
            }
            leftUpdate.Empty();

            // Add the right child's update matrix
            const int rightUpdateSize = rightUpdate.Height();
            for( int jChild=0; jChild<rightUpdateSize; ++jChild )
            {
                const int jFront = node.rightRelInds[jChild];
                for( int iChild=jChild; iChild<rightUpdateSize; ++iChild )
                {
                    const int iFront = node.rightRelInds[iChild];
                    const F value = rightUpdate.Get(iChild,jChild);
                    DEBUG_ONLY(
                        if( iFront < jFront )
                            LogicError("Tried to update upper triangle");
                    )
                    if( jFront < node.size )
                        frontL.Update( iFront, jFront, value );
                    else if( iFront >= node.size )
                        frontBR.Update
                        ( iFront-node.size, jFront-node.size, value );
                }
            }
            rightUpdate.Empty();
        }

        // Call the custom partial LDL
        ProcessFront
        ( frontL, frontBR, pivTol, 
          regCand.localNodes[s], reg.localNodes[s] );
        frontL.GetDiagonal( front.diag );
        SetDiagonal( frontL, F(1) );
    }
}

} // namespace reg_ldl
} // namespace El

#endif // ifndef EL_OPTIMIZATION_REGLDL_PROCESSLOCALTREE_HPP
