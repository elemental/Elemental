/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_BLOCKDISTMATRIX_HPP
#define EL_BLOCKDISTMATRIX_HPP

namespace El {

struct BlockDistData
{
    Dist colDist, rowDist;
    Int blockHeight, blockWidth;
    int colAlign, rowAlign; 
    Int colCut, rowCut;
    int root;  // relevant for [o ,o ]/[MD,* ]/[* ,MD]
    const Grid* grid;

    BlockDistData() { }

    template<typename T>
    BlockDistData( const AbstractBlockDistMatrix<T>& A )
    : colDist(A.ColDist()), rowDist(A.RowDist()), 
      blockHeight(A.BlockHeight()), blockWidth(A.BlockWidth()),
      colAlign(A.ColAlign()), rowAlign(A.RowAlign()),
      colCut(A.ColCut()), rowCut(A.RowCut()),
      root(A.Root()), grid(&A.Grid())
    { }
};
inline bool operator==( const BlockDistData& A, const BlockDistData& B )
{ return A.colDist     == B.colDist &&
         A.rowDist     == B.rowDist &&
         A.blockHeight == B.blockHeight &&
         A.blockWidth  == B.blockWidth &&
         A.colAlign    == B.colAlign &&
         A.rowAlign    == B.rowAlign &&
         A.root        == B.root &&
         A.grid        == B.grid; }


} // namespace El

#include "./BlockDistMatrix/Abstract.hpp"
#include "./BlockDistMatrix/CIRC_CIRC.hpp"
#include "./BlockDistMatrix/MC_MR.hpp"
#include "./BlockDistMatrix/MC_STAR.hpp"
#include "./BlockDistMatrix/MD_STAR.hpp"
#include "./BlockDistMatrix/MR_MC.hpp"
#include "./BlockDistMatrix/MR_STAR.hpp"
#include "./BlockDistMatrix/STAR_MC.hpp"
#include "./BlockDistMatrix/STAR_MD.hpp"
#include "./BlockDistMatrix/STAR_MR.hpp"
#include "./BlockDistMatrix/STAR_STAR.hpp"
#include "./BlockDistMatrix/STAR_VC.hpp"
#include "./BlockDistMatrix/STAR_VR.hpp"
#include "./BlockDistMatrix/VC_STAR.hpp"
#include "./BlockDistMatrix/VR_STAR.hpp"

namespace El {

#ifdef EL_HAVE_SCALAPACK
template<typename T>
inline typename blacs::Desc
FillDesc( const BlockDistMatrix<T>& A, int context )
{
    if( A.ColCut() != 0 || A.RowCut() != 0 )
        LogicError("Cannot produce a meaningful descriptor if nonzero cut");
    typename blacs::Desc desc = 
        { 1, context, A.Height(), A.Width(), A.BlockHeight(), A.BlockWidth(), 
          A.ColAlign(), A.RowAlign(), A.LDim() };
    return desc;
}
#endif

template<typename T>
inline void AssertSameGrids( const AbstractBlockDistMatrix<T>& A ) { }

template<typename T1,typename T2>
inline void AssertSameGrids
( const AbstractBlockDistMatrix<T1>& A1, const AbstractBlockDistMatrix<T2>& A2 )
{
    if( A1.Grid() != A2.Grid() )
        LogicError("Grids did not match");
}

template<typename T1,typename T2,typename... Args>
inline void AssertSameGrids
( const AbstractBlockDistMatrix<T1>& A1, const AbstractBlockDistMatrix<T2>& A2,
  Args&... args )
{
    if( A1.Grid() != A2.Grid() )
        LogicError("Grids did not match");
    AssertSameGrids( A2, args... );
}

template<typename T>
inline void AssertSameDists( const AbstractBlockDistMatrix<T>& A ) { }

template<typename T>
inline void AssertSameDists
( const AbstractBlockDistMatrix<T>& A1, const AbstractBlockDistMatrix<T>& A2 ) 
{
    if( A1.ColDist() != A2.ColDist() || A1.RowDist() != A2.RowDist() )
        LogicError("Distributions did not match");
}

template<typename T,typename... Args>
inline void AssertSameDists
( const AbstractBlockDistMatrix<T>& A1, const AbstractBlockDistMatrix<T>& A2, 
  Args&... args )
{
    if( A1.ColDist() != A2.ColDist() || A1.RowDist() != A2.RowDist() )
        LogicError("Distributions did not match");
    AssertSameDists( A2, args... );
}

} // namespace El

#endif // ifndef EL_BLOCKDISTMATRIX_HPP
