/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_BLOCKDISTMATRIX_STAR_MR_DECL_HPP
#define ELEM_BLOCKDISTMATRIX_STAR_MR_DECL_HPP

namespace elem {

// Partial specialization to A[* ,MR].
//
// The columns of these distributed matrices will be replicated on all 
// processes (*), and the rows will be distributed like "Matrix Rows" (MR).
// Thus the rows will be distributed among rows of the process grid.
template<typename T>
class BlockDistMatrix<T,STAR,MR> : public GeneralBlockDistMatrix<T,STAR,MR>
{
    // TODO
};

} // namespace elem

#endif // ifndef ELEM_BLOCKDISTMATRIX_STAR_MR_DECL_HPP
