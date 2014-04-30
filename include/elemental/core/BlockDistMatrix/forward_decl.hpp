/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_BLOCKDISTMATRIX_FORWARD_DECL_HPP
#define ELEM_BLOCKDISTMATRIX_FORWARD_DECL_HPP

namespace elem {

template<typename T>
class AbstractBlockDistMatrix;

template<typename T,Dist U=MC,Dist V=MR>
class GeneralBlockDistMatrix;

template<typename T,Dist U=MC,Dist V=MR>
class BlockDistMatrix;

} // namespace elem

#endif // ifndef ELEM_BLOCKDISTMATRIX_FORWARD_DECL_HPP
