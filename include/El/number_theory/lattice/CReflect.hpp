/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_LATTICE_CREFLECT_C_HPP
#define EL_LATTICE_CREFLECT_C_HPP

namespace El {

// Lattice
// -------
inline ElLLLVariant CReflect( LLLVariant var )
{ return static_cast<ElLLLVariant>(var); }
inline LLLVariant CReflect( ElLLLVariant var )
{ return static_cast<LLLVariant>(var); }

inline LLLInfo<float> CReflect( const ElLLLInfo_s& infoC )
{
    LLLInfo<float> info;
    info.delta = infoC.delta;
    info.eta = infoC.eta;
    info.rank = infoC.rank;
    info.nullity = infoC.nullity;
    info.numSwaps = infoC.numSwaps;
    info.firstSwap = infoC.firstSwap;
    info.logVol = infoC.logVol;
    return info;
}

inline LLLInfo<double> CReflect( const ElLLLInfo_d& infoC )
{
    LLLInfo<double> info;
    info.delta = infoC.delta;
    info.eta = infoC.eta;
    info.rank = infoC.rank;
    info.nullity = infoC.nullity;
    info.numSwaps = infoC.numSwaps;
    info.firstSwap = infoC.firstSwap;
    info.logVol = infoC.logVol;
    return info;
}

inline ElLLLInfo_s CReflect( const LLLInfo<float>& info )
{
    ElLLLInfo_s infoC;
    infoC.delta = info.delta;
    infoC.eta = info.eta;
    infoC.rank = info.rank;
    infoC.nullity = info.nullity;
    infoC.numSwaps = info.numSwaps;
    infoC.firstSwap = info.firstSwap;
    infoC.logVol = info.logVol;
    return infoC;
}

inline ElLLLInfo_d CReflect( const LLLInfo<double>& info )
{
    ElLLLInfo_d infoC;
    infoC.delta = info.delta;
    infoC.eta = info.eta;
    infoC.rank = info.rank;
    infoC.nullity = info.nullity;
    infoC.numSwaps = info.numSwaps;
    infoC.firstSwap = info.firstSwap;
    infoC.logVol = info.logVol;
    return infoC;
}

inline LLLCtrl<float> CReflect( const ElLLLCtrl_s& ctrlC )
{
    LLLCtrl<float> ctrl;
    ctrl.delta = ctrlC.delta;
    ctrl.eta = ctrlC.eta;
    ctrl.variant = CReflect(ctrlC.variant);
    ctrl.recursive = ctrlC.recursive;
    ctrl.cutoff = ctrlC.cutoff;
    ctrl.presort = ctrlC.presort;
    ctrl.smallestFirst = ctrlC.smallestFirst;
    ctrl.reorthogTol = ctrlC.reorthogTol;
    ctrl.numOrthog = ctrlC.numOrthog;
    ctrl.zeroTol = ctrlC.zeroTol;
    ctrl.blockingThresh = ctrlC.blockingThresh;
    ctrl.progress = ctrlC.progress;
    ctrl.time = ctrlC.time;
    ctrl.jumpstart = ctrlC.jumpstart;
    ctrl.startCol = ctrlC.startCol;
    return ctrl;
}

inline LLLCtrl<double> CReflect( const ElLLLCtrl_d& ctrlC )
{
    LLLCtrl<double> ctrl;
    ctrl.delta = ctrlC.delta;
    ctrl.eta = ctrlC.eta;
    ctrl.variant = CReflect(ctrlC.variant);
    ctrl.recursive = ctrlC.recursive;
    ctrl.cutoff = ctrlC.cutoff;
    ctrl.presort = ctrlC.presort;
    ctrl.smallestFirst = ctrlC.smallestFirst;
    ctrl.reorthogTol = ctrlC.reorthogTol;
    ctrl.numOrthog = ctrlC.numOrthog;
    ctrl.zeroTol = ctrlC.zeroTol;
    ctrl.blockingThresh = ctrlC.blockingThresh;
    ctrl.progress = ctrlC.progress;
    ctrl.time = ctrlC.time;
    ctrl.jumpstart = ctrlC.jumpstart;
    ctrl.startCol = ctrlC.startCol;
    return ctrl;
}

inline ElLLLCtrl_s CReflect( const LLLCtrl<float>& ctrl )
{
    ElLLLCtrl_s ctrlC;
    ctrlC.delta = ctrl.delta;
    ctrlC.eta = ctrl.eta;
    ctrlC.variant = CReflect(ctrl.variant);
    ctrlC.recursive = ctrl.recursive;
    ctrlC.cutoff = ctrl.cutoff;
    ctrlC.presort = ctrl.presort;
    ctrlC.smallestFirst = ctrl.smallestFirst;
    ctrlC.reorthogTol = ctrl.reorthogTol;
    ctrlC.numOrthog = ctrl.numOrthog;
    ctrlC.zeroTol = ctrl.zeroTol;
    ctrlC.blockingThresh = ctrl.blockingThresh;
    ctrlC.progress = ctrl.progress;
    ctrlC.time = ctrl.time;
    ctrlC.jumpstart = ctrl.jumpstart;
    ctrlC.startCol = ctrl.startCol;
    return ctrlC;
}

inline ElLLLCtrl_d CReflect( const LLLCtrl<double>& ctrl )
{
    ElLLLCtrl_d ctrlC;
    ctrlC.delta = ctrl.delta;
    ctrlC.eta = ctrl.eta;
    ctrlC.variant = CReflect(ctrl.variant);
    ctrlC.recursive = ctrl.recursive;
    ctrlC.cutoff = ctrl.cutoff;
    ctrlC.presort = ctrl.presort;
    ctrlC.smallestFirst = ctrl.smallestFirst;
    ctrlC.reorthogTol = ctrl.reorthogTol;
    ctrlC.numOrthog = ctrl.numOrthog;
    ctrlC.zeroTol = ctrl.zeroTol;
    ctrlC.blockingThresh = ctrl.blockingThresh;
    ctrlC.progress = ctrl.progress;
    ctrlC.time = ctrl.time;
    ctrlC.jumpstart = ctrl.jumpstart;
    ctrlC.startCol = ctrl.startCol;
    return ctrlC;
}

} // namespace El

#endif // ifndef EL_LATTICE_CREFLECT_C_HPP
