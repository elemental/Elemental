/*
   Copyright (c) 2009-2015, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, Stanford University, and the
   Georgia Insitute of Technology.
   All rights reserved.
 
   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_SPARSEDIRECT_NUMERIC_LOWERSOLVE_LOCAL_HPP
#define EL_SPARSEDIRECT_NUMERIC_LOWERSOLVE_LOCAL_HPP

namespace El {

template<typename F> 
inline void LocalLowerForwardSolve
( const DistSymmInfo& info, 
  const DistSymmFrontTree<F>& L, DistNodalMultiVec<F>& X )
{
    DEBUG_ONLY(CallStackEntry cse("LocalLowerForwardSolve"))
    const Int numLocalNodes = info.localNodes.size();
    const Int width = X.Width();

    const SymmFrontType frontType = L.frontType;
    if( Unfactored(frontType) )
        LogicError("Nonsensical front type for solve");
    const bool blocked = BlockFactorization( frontType );
    const bool pivoted = PivotedFactorization( frontType );

    for( Int s=0; s<numLocalNodes; ++s )
    {
        const SymmNodeInfo& node = info.localNodes[s];
        const SymmFront<F>& front = L.localFronts[s];
        const Matrix<F>& frontL = front.frontL;
        Matrix<F>& W = front.work;

        // Set up a workspace
        W.Resize( frontL.Height(), width );
        Matrix<F> WT, WB;
        PartitionDown( W, WT, WB, node.size );
        WT = X.localNodes[s];
        Zero( WB );

        // Update using the children (if they exist)
        const Int numChildren = node.children.size();
        if( numChildren == 2 )
        {
            const Int leftInd = node.children[0];
            const Int rightInd = node.children[1];
            Matrix<F>& leftWork = L.localFronts[leftInd].work;
            Matrix<F>& rightWork = L.localFronts[rightInd].work;
            const Int leftNodeSize = info.localNodes[leftInd].size;
            const Int rightNodeSize = info.localNodes[rightInd].size;
            const Int leftUpdateSize = leftWork.Height()-leftNodeSize;
            const Int rightUpdateSize = rightWork.Height()-rightNodeSize;

            // Add the left child's update onto ours
            auto leftUpdate = 
                LockedView( leftWork, leftNodeSize, 0, leftUpdateSize, width );
            for( Int iChild=0; iChild<leftUpdateSize; ++iChild )
            {
                const Int iFront = node.leftRelInds[iChild]; 
                for( Int j=0; j<width; ++j )
                    W.Update( iFront, j, leftUpdate.Get(iChild,j) );
            }
            leftWork.Empty();

            // Add the right child's update onto ours
            auto rightUpdate =
                LockedView
                ( rightWork, rightNodeSize, 0, rightUpdateSize, width );
            for( Int iChild=0; iChild<rightUpdateSize; ++iChild )
            {
                const Int iFront = node.rightRelInds[iChild];
                for( Int j=0; j<width; ++j )
                    W.Update( iFront, j, rightUpdate.Get(iChild,j) );
            }
            rightWork.Empty();
        }
        // else numChildren == 0

        // Solve against this front
        if( blocked )
            FrontBlockLowerForwardSolve( frontL, W );
        else if( pivoted )
            FrontIntraPivLowerForwardSolve( frontL, front.piv, W );
        else
            FrontLowerForwardSolve( frontL, W );

        // Store this node's portion of the result
        X.localNodes[s] = WT;
    }
}

// This is an exact copy of the DistNodalMultiVec version...
template<typename F> 
inline void LocalLowerForwardSolve
( const DistSymmInfo& info, 
  const DistSymmFrontTree<F>& L, DistNodalMatrix<F>& X )
{
    DEBUG_ONLY(CallStackEntry cse("LocalLowerForwardSolve"))
    const Int numLocalNodes = info.localNodes.size();
    const Int width = X.Width();

    const SymmFrontType frontType = L.frontType;
    if( Unfactored(frontType) )
        LogicError("Nonsensical front type for solve");
    const bool blocked = BlockFactorization( frontType );
    const bool pivoted = PivotedFactorization( frontType );

    for( Int s=0; s<numLocalNodes; ++s )
    {
        const SymmNodeInfo& node = info.localNodes[s];
        const SymmFront<F>& front = L.localFronts[s];
        const Matrix<F>& frontL = front.frontL;
        Matrix<F>& W = front.work;

        // Set up a workspace
        W.Resize( frontL.Height(), width );
        Matrix<F> WT, WB;
        PartitionDown( W, WT, WB, node.size );
        WT = X.localNodes[s];
        Zero( WB );

        // Update using the children (if they exist)
        const Int numChildren = node.children.size();
        if( numChildren == 2 )
        {
            const Int leftInd = node.children[0];
            const Int rightInd = node.children[1];
            Matrix<F>& leftWork = L.localFronts[leftInd].work;
            Matrix<F>& rightWork = L.localFronts[rightInd].work;
            const Int leftNodeSize = info.localNodes[leftInd].size;
            const Int rightNodeSize = info.localNodes[rightInd].size;
            const Int leftUpdateSize = leftWork.Height()-leftNodeSize;
            const Int rightUpdateSize = rightWork.Height()-rightNodeSize;

            // Add the left child's update onto ours
            auto leftUpdate =
                LockedView( leftWork, leftNodeSize, 0, leftUpdateSize, width );
            for( Int iChild=0; iChild<leftUpdateSize; ++iChild )
            {
                const Int iFront = node.leftRelInds[iChild]; 
                for( Int j=0; j<width; ++j )
                    W.Update( iFront, j, leftUpdate.Get(iChild,j) );
            }
            leftWork.Empty();

            // Add the right child's update onto ours
            auto rightUpdate =
                LockedView
                ( rightWork, rightNodeSize, 0, rightUpdateSize, width );
            for( Int iChild=0; iChild<rightUpdateSize; ++iChild )
            {
                const Int iFront = node.rightRelInds[iChild];
                for( Int j=0; j<width; ++j )
                    W.Update( iFront, j, rightUpdate.Get(iChild,j) );
            }
            rightWork.Empty();
        }
        // else numChildren == 0

        // Solve against this front
        if( blocked )
            FrontBlockLowerForwardSolve( frontL, W );
        else if( pivoted )
            FrontIntraPivLowerForwardSolve( frontL, front.piv, W );
        else
            FrontLowerForwardSolve( frontL, W );

        // Store this node's portion of the result
        X.localNodes[s] = WT;
    }
}

template<typename F> 
inline void LocalLowerBackwardSolve
( const DistSymmInfo& info, 
  const DistSymmFrontTree<F>& L, DistNodalMultiVec<F>& X,
  bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("LocalLowerBackwardSolve"))
    const Int numLocalNodes = info.localNodes.size();
    const Int width = X.Width();

    const SymmFrontType frontType = L.frontType;
    if( Unfactored(frontType) )
        LogicError("Nonsensical front type for solve");
    const bool blocked = BlockFactorization( frontType );
    const bool pivoted = PivotedFactorization( frontType );

    for( Int s=numLocalNodes-2; s>=0; --s )
    {
        const SymmNodeInfo& node = info.localNodes[s];
        const SymmFront<F>& front = L.localFronts[s];
        const Matrix<F>& frontL = front.frontL;
        Matrix<F>& W = front.work;

        // Set up a workspace
        W.Resize( frontL.Height(), width );
        Matrix<F> WT, WB;
        PartitionDown( W, WT, WB, node.size );
        WT = X.localNodes[s];

        // Update using the parent
        const Int parent = node.parent;
        DEBUG_ONLY(
            if( parent < 0 )
                LogicError("Parent index was negative: ",parent);
            if( parent >= numLocalNodes )  
                LogicError
                ("Parent index was too large: ",parent," >= ",numLocalNodes);
        )
        Matrix<F>& parentWork = L.localFronts[parent].work;
        const SymmNodeInfo& parentNode = info.localNodes[parent];
        const Int currentUpdateSize = WB.Height();
        const std::vector<Int>& parentRelInds = 
          ( node.onLeft ? parentNode.leftRelInds : parentNode.rightRelInds );
        for( Int iCurrent=0; iCurrent<currentUpdateSize; ++iCurrent )
        {
            const Int iParent = parentRelInds[iCurrent];
            for( Int j=0; j<width; ++j )
                WB.Set( iCurrent, j, parentWork.Get(iParent,j) );
        }

        // The left child is numbered lower than the right child, so 
        // we can safely free the parent's work if we are the left child
        if( node.onLeft )
        {
            parentWork.Empty();
            if( parent == numLocalNodes-1 )
                L.distFronts[0].work1d.Empty();
        }

        // Solve against this front
        if( blocked )
            FrontBlockLowerBackwardSolve( frontL, W, conjugate );
        else if( pivoted )
            FrontIntraPivLowerBackwardSolve( frontL, front.piv, W, conjugate );
        else
            FrontLowerBackwardSolve( frontL, W, conjugate );

        // Store this node's portion of the result
        X.localNodes[s] = WT;
    }

    // Ensure that all of the temporary buffers are freed (this is overkill)
    L.distFronts[0].work1d.Empty();
    for( Int s=0; s<numLocalNodes; ++s )
        L.localFronts[s].work.Empty();
}

template<typename F> 
inline void LocalLowerBackwardSolve
( const DistSymmInfo& info, 
  const DistSymmFrontTree<F>& L, DistNodalMatrix<F>& X,
  bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("LocalLowerBackwardSolve"))
    const Int numLocalNodes = info.localNodes.size();
    const Int width = X.Width();

    const SymmFrontType frontType = L.frontType;
    if( Unfactored(frontType) )
        LogicError("Nonsensical front type for solve");
    const bool blocked = BlockFactorization( frontType );
    const bool pivoted = PivotedFactorization( frontType );

    for( Int s=numLocalNodes-2; s>=0; --s )
    {
        const SymmNodeInfo& node = info.localNodes[s];
        const SymmFront<F>& front = L.localFronts[s];
        const Matrix<F>& frontL = front.frontL;
        Matrix<F>& W = front.work;

        // Set up a workspace
        W.Resize( frontL.Height(), width );
        Matrix<F> WT, WB;
        PartitionDown( W, WT, WB, node.size );
        WT = X.localNodes[s];

        // Update using the parent
        const Int parent = node.parent;
        DEBUG_ONLY(
            if( parent < 0 )
                LogicError("Parent index was negative: ",parent);
            if( parent >= numLocalNodes )  
                LogicError
                ("Parent index was too large: ",parent," >= ",numLocalNodes);
        )
        Matrix<F>& parentWork = L.localFronts[parent].work;
        const SymmNodeInfo& parentNode = info.localNodes[parent];
        const Int currentUpdateSize = WB.Height();
        const std::vector<Int>& parentRelInds = 
          ( node.onLeft ? parentNode.leftRelInds : parentNode.rightRelInds );
        for( Int iCurrent=0; iCurrent<currentUpdateSize; ++iCurrent )
        {
            const Int iParent = parentRelInds[iCurrent];
            for( Int j=0; j<width; ++j )
                WB.Set( iCurrent, j, parentWork.Get(iParent,j) );
        }

        // The left child is numbered lower than the right child, so 
        // we can safely free the parent's work if we are the left child
        if( node.onLeft )
        {
            parentWork.Empty();
            if( parent == numLocalNodes-1 )
                L.distFronts[0].work2d.Empty();
        }

        // Solve against this front
        if( blocked )
            FrontBlockLowerBackwardSolve( frontL, W, conjugate );
        else if( pivoted )
            FrontIntraPivLowerBackwardSolve( frontL, front.piv, W, conjugate );
        else
            FrontLowerBackwardSolve( frontL, W, conjugate );

        // Store this node's portion of the result
        X.localNodes[s] = WT;
    }

    // Ensure that all of the temporary buffers are freed (this is overkill)
    L.distFronts[0].work2d.Empty();
    for( Int s=0; s<numLocalNodes; ++s )
        L.localFronts[s].work.Empty();
}

} // namespace El

#endif // ifndef EL_SPARSEDIRECT_NUMERIC_LOWERSOLVE_LOCAL_HPP
