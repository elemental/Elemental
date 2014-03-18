/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_BLOCKDISTMATRIX_STAR_STAR_DECL_HPP
#define ELEM_BLOCKDISTMATRIX_STAR_STAR_DECL_HPP

namespace elem {

// Partial specialization to A[* ,* ].
//
// The entire matrix is replicated across all processes.
template<typename T>
class BlockDistMatrix<T,STAR,STAR> : public GeneralBlockDistMatrix<T,STAR,STAR>
{
    // TODO
};

} // namespace elem

#endif // ifndef ELEM_BLOCKDISTMATRIX_STAR_STAR_DECL_HPP
