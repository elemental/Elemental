/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_DISTMATRIX_FORWARD_DECL_HPP
#define EL_DISTMATRIX_FORWARD_DECL_HPP

namespace El {

template<typename T>
class AbstractDistMatrix;

template<typename T,Dist U=MC,Dist V=MR>
class GeneralDistMatrix;

template<typename T,Dist U=MC,Dist V=MR>
class DistMatrix;

template<typename T>
class DynamicDistMatrix;

} // namespace El

#endif // ifndef EL_DISTMATRIX_FORWARD_DECL_HPP
