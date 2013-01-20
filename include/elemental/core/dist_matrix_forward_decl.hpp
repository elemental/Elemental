/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef CORE_DISTMATRIX_FORWARD_DECL_HPP
#define CORE_DISTMATRIX_FORWARD_DECL_HPP

namespace elem {

template<typename T,typename Int=int>
class AbstractDistMatrix;

template<typename T,Distribution ColDist=MC,Distribution RowDist=MR,
         typename Int=int>
class DistMatrix;

} // namespace elem

#endif // ifndef CORE_DISTMATRIX_FORWARD_DECL_HPP
