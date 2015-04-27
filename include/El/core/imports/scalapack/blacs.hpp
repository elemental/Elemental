/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_IMPORTS_SCALAPACK_BLACS_HPP
#define EL_IMPORTS_SCALAPACK_BLACS_HPP

#ifdef EL_HAVE_SCALAPACK

namespace El {
namespace blacs {

int Handle( MPI_Comm comm );
int GridInit( int bhandle, bool colMajor, int gridHeight, int gridWidth );
int GridHeight( int context );
int GridWidth( int context );
int GridRow( int context );
int GridCol( int context );
void FreeHandle( int bhandle );
void FreeGrid( int context );
void Exit( bool finished=false );

typedef typename std::array<int,9> Desc;

} // namespace blacs
} // namespace El

#endif // ifdef EL_HAVE_SCALAPACK
#endif // ifndef EL_IMPORTS_SCALAPACK_BLACS_HPP
