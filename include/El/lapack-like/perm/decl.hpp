/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_PERM_DECL_HPP
#define EL_PERM_DECL_HPP

namespace El {

bool PermutationParity( const Matrix<Int>& origPerm );
template<Dist UPerm>
bool PermutationParity( const DistMatrix<Int,UPerm,STAR>& origPerm );

bool PivotParity( const Matrix<Int>& p, Int pivotOffset=0 );
// TODO: Generalize implementation?
bool PivotParity( const DistMatrix<Int,VC,STAR>& p, Int pivotOffset=0 );

} // namespace El

#endif // ifndef EL_PERM_DECL_HPP
