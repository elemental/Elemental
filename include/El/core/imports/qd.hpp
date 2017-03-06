/*
   Copyright (c) 2009-2017, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_IMPORTS_QD_HPP
#define EL_IMPORTS_QD_HPP

#ifdef EL_HAVE_QD
#include <qd/qd_real.h>

#include <El/core/imports/qd/DoubleDouble.hpp>
#include <El/core/imports/qd/QuadDouble.hpp>

namespace El {

// To be called internally by Elemental
void InitializeQD();
void FinalizeQD();

} // namespace El
#endif // ifdef EL_HAVE_QD

#endif // ifndef EL_IMPORTS_QD_HPP
