/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_BLAS_LEVEL2_CREFLECT_C_HPP
#define EL_BLAS_LEVEL2_CREFLECT_C_HPP

namespace El {

template<typename T>
ElSymvCtrl CReflect( const SymvCtrl<T>& ctrl )
{
    ElSymvCtrl ctrlC;
    ctrlC.bsize = ctrl.bsize;
    ctrlC.avoidTrmvBasedLocalSymv = ctrl.avoidTrmvBasedLocalSymv;
    return ctrlC;
}

template<typename T>
SymvCtrl<T> CReflect( const ElSymvCtrl& ctrlC )
{
    SymvCtrl<T> ctrl;
    ctrl.bsize = ctrlC.bsize;
    ctrl.avoidTrmvBasedLocalSymv = ctrlC.avoidTrmvBasedLocalSymv;
    return ctrl;
}

} // namespace El

#endif // ifndef EL_BLAS_LEVEL2_CREFLECT_C_HPP
