/*
   Copyright (c) 2009-2013, Jack Poulson
                      2013, Michael C. Grant
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
%module elem_mpi

%include "common.swg"

/*
 * MPI
 */

%import  "std_except.i"
%include "elemental/core/imports/mpi.hpp"

// The communication routines have not yet been exposed to SWIG. We have to decide how
// to do that. For instance, for Python, do we typemap all of the arrays to NumPy objects?
// And what to do about the ValueInt<T> type? NumPy record objects, perhaps?
