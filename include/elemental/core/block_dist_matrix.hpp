/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_BLOCKDISTMATRIX_HPP
#define ELEM_BLOCKDISTMATRIX_HPP

namespace elem {

struct BlockDistData
{
    Dist colDist, rowDist;
    Int blockHeight, blockWidth;
    Int colAlign, rowAlign; 
    Int colCut, rowCut;
    Int root;  // relevant for [o ,o ]/[MD,* ]/[* ,MD]
    const Grid* grid;

    template<typename T,Dist U,Dist V>
    BlockDistData( const GeneralBlockDistMatrix<T,U,V>& A )
    : colDist(U), rowDist(V), 
      blockHeight(A.BlockHeight()), blockWidth(A.BlockWidth()),
      colAlign(A.ColAlign()), rowAlign(A.RowAlign()),
      colCut(A.ColCut()), rowCut(A.RowCut()),
      root(A.Root()), grid(&A.Grid())
    { }
};

} // namespace elem

#include "./block_dist_matrix/abstract.hpp"
#include "./block_dist_matrix/general.hpp"
#include "./block_dist_matrix/circ_circ.hpp"
#include "./block_dist_matrix/mc_mr.hpp"
#include "./block_dist_matrix/mc_star.hpp"
#include "./block_dist_matrix/md_star.hpp"
#include "./block_dist_matrix/mr_mc.hpp"
#include "./block_dist_matrix/mr_star.hpp"
#include "./block_dist_matrix/star_mc.hpp"
#include "./block_dist_matrix/star_md.hpp"
#include "./block_dist_matrix/star_mr.hpp"
#include "./block_dist_matrix/star_star.hpp"
#include "./block_dist_matrix/star_vc.hpp"
#include "./block_dist_matrix/star_vr.hpp"
#include "./block_dist_matrix/vc_star.hpp"
#include "./block_dist_matrix/vr_star.hpp"

namespace elem {

#ifdef ELEM_HAVE_SCALAPACK
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

} // namespace elem

#endif // ifndef ELEM_BLOCKDISTMATRIX_HPP
