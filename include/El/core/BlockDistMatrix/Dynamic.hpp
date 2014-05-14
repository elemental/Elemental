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
};

} // namespace El

#endif // ifndef EL_BLOCKDISTMATRIX_DYNAMIC_DECL_HPP
