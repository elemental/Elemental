/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_LAPACK_CONDENSE_CREFLECT_C_HPP
#define EL_LAPACK_CONDENSE_CREFLECT_C_HPP

namespace El {

inline ElHermitianTridiagApproach
CReflect( HermitianTridiagApproach approach )
{ return static_cast<ElHermitianTridiagApproach>( approach ); }

inline HermitianTridiagApproach
CReflect( ElHermitianTridiagApproach approach )
{ return static_cast<HermitianTridiagApproach>( approach ); }

template<typename Field>
ElHermitianTridiagCtrl
CReflect( const HermitianTridiagCtrl<Field>& ctrl )
{
    ElHermitianTridiagCtrl ctrlC;
    ctrlC.approach = CReflect(ctrl.approach);
    ctrlC.order = CReflect(ctrl.order);
    ctrlC.symvCtrl = CReflect(ctrl.symvCtrl);
    return ctrlC;
}

template<typename Field>
HermitianTridiagCtrl<Field>
CReflect( const ElHermitianTridiagCtrl& ctrlC )
{
    HermitianTridiagCtrl<Field> ctrl;
    ctrl.approach = CReflect(ctrlC.approach);
    ctrl.order = CReflect(ctrlC.order);
    ctrl.symvCtrl = CReflect<Field>(ctrlC.symvCtrl);
    return ctrl;
}

} // namespace El

#endif // ifndef EL_LAPACK_CONDENSE_CREFLECT_C_HPP
