/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_LAPACK_FACTOR_CREFLECT_C_HPP
#define EL_LAPACK_FACTOR_CREFLECT_C_HPP

namespace El {

inline ElLDLPivotType CReflect( LDLPivotType pivotType )
{ return static_cast<ElLDLPivotType>( pivotType ); }

inline LDLPivotType CReflect( ElLDLPivotType pivotType )
{ return static_cast<LDLPivotType>( pivotType ); }

inline ElLDLPivotCtrl_s CReflect( const LDLPivotCtrl<float>& ctrl )
{
    ElLDLPivotCtrl_s ctrlC;
    ctrlC.pivotType = CReflect(ctrl.pivotType);
    ctrlC.gamma = ctrl.gamma;
    return ctrlC;
}
inline ElLDLPivotCtrl_d CReflect( const LDLPivotCtrl<double>& ctrl )
{
    ElLDLPivotCtrl_d ctrlC;
    ctrlC.pivotType = CReflect(ctrl.pivotType);
    ctrlC.gamma = ctrl.gamma;
    return ctrlC;
}

inline LDLPivotCtrl<float> CReflect( const ElLDLPivotCtrl_s& ctrlC )
{
    LDLPivotCtrl<float> ctrl;
    ctrl.pivotType = CReflect(ctrlC.pivotType);
    ctrl.gamma = ctrlC.gamma;
    return ctrl;
}
inline LDLPivotCtrl<double> CReflect( const ElLDLPivotCtrl_d& ctrlC )
{
    LDLPivotCtrl<double> ctrl;
    ctrl.pivotType = CReflect(ctrlC.pivotType);
    ctrl.gamma = ctrlC.gamma;
    return ctrl;
}

inline ElLDLPivot CReflect( const LDLPivot& pivot )
{
    ElLDLPivot pivotC;
    pivotC.nb = pivot.nb;
    pivotC.from[0] = pivot.from[0];
    pivotC.from[1] = pivot.from[1];
    return pivotC;
}

inline LDLPivot CReflect( const ElLDLPivot& pivotC )
{
    LDLPivot pivot;
    pivot.nb = pivotC.nb;
    pivot.from[0] = pivotC.from[0];
    pivot.from[1] = pivotC.from[1];
    return pivot;
}

inline ElInertiaType CReflect( const InertiaType& inertia )
{ 
    ElInertiaType inertiaC;
    inertiaC.numPositive = inertia.numPositive;
    inertiaC.numNegative = inertia.numNegative;
    inertiaC.numZero = inertia.numZero;
    return inertiaC;
}

inline InertiaType CReflect( const ElInertiaType& inertiaC )
{ 
    InertiaType inertia;
    inertia.numPositive = inertiaC.numPositive;
    inertia.numNegative = inertiaC.numNegative;
    inertia.numZero = inertiaC.numZero;
    return inertia;
}

inline ElQRCtrl_s CReflect( const QRCtrl<float>& ctrl )
{ 
    ElQRCtrl_s ctrlC;
    ctrlC.colPiv = ctrl.colPiv;
    ctrlC.boundRank = ctrl.boundRank;
    ctrlC.maxRank = ctrl.maxRank;
    ctrlC.adaptive = ctrl.adaptive;
    ctrlC.tol = ctrl.tol;
    ctrlC.alwaysRecomputeNorms = ctrl.alwaysRecomputeNorms;
    ctrlC.smallestFirst = ctrl.smallestFirst;
    return ctrlC;
}
inline ElQRCtrl_d CReflect( const QRCtrl<double>& ctrl )
{ 
    ElQRCtrl_d ctrlC;
    ctrlC.colPiv = ctrl.colPiv;
    ctrlC.boundRank = ctrl.boundRank;
    ctrlC.maxRank = ctrl.maxRank;
    ctrlC.adaptive = ctrl.adaptive;
    ctrlC.tol = ctrl.tol;
    ctrlC.alwaysRecomputeNorms = ctrl.alwaysRecomputeNorms;
    ctrlC.smallestFirst = ctrl.smallestFirst;
    return ctrlC;
}

inline QRCtrl<float> CReflect( const ElQRCtrl_s& ctrlC )
{ 
    QRCtrl<float> ctrl;
    ctrl.colPiv = ctrlC.colPiv;
    ctrl.boundRank = ctrlC.boundRank;
    ctrl.maxRank = ctrlC.maxRank;
    ctrl.adaptive = ctrlC.adaptive;
    ctrl.tol = ctrlC.tol;
    ctrl.alwaysRecomputeNorms = ctrlC.alwaysRecomputeNorms;
    ctrl.smallestFirst = ctrlC.smallestFirst;
    return ctrl;
}
inline QRCtrl<double> CReflect( const ElQRCtrl_d& ctrlC )
{ 
    QRCtrl<double> ctrl;
    ctrl.colPiv = ctrlC.colPiv;
    ctrl.boundRank = ctrlC.boundRank;
    ctrl.maxRank = ctrlC.maxRank;
    ctrl.adaptive = ctrlC.adaptive;
    ctrl.tol = ctrlC.tol;
    ctrl.alwaysRecomputeNorms = ctrlC.alwaysRecomputeNorms;
    ctrl.smallestFirst = ctrlC.smallestFirst;
    return ctrl;
}

inline RegSolveAlg CReflect( ElRegSolveAlg alg ) 
{ return static_cast<RegSolveAlg>(alg); }
inline ElRegSolveAlg CReflect( RegSolveAlg alg )
{ return static_cast<ElRegSolveAlg>(alg); }

inline ElRegSolveCtrl_s CReflect( const RegSolveCtrl<float>& ctrl )
{
    ElRegSolveCtrl_s ctrlC;
    ctrlC.alg          = CReflect(ctrl.alg);
    ctrlC.relTol       = ctrl.relTol;
    ctrlC.relTolRefine = ctrl.relTolRefine;
    ctrlC.maxRefineIts = ctrl.maxRefineIts;
    ctrlC.restart      = ctrl.restart;
    ctrlC.progress     = ctrl.progress;
    ctrlC.time         = ctrl.time;
    return ctrlC;
}

inline ElRegSolveCtrl_d CReflect( const RegSolveCtrl<double>& ctrl )
{
    ElRegSolveCtrl_d ctrlC;
    ctrlC.alg          = CReflect(ctrl.alg);
    ctrlC.relTol       = ctrl.relTol;
    ctrlC.relTolRefine = ctrl.relTolRefine;
    ctrlC.maxRefineIts = ctrl.maxRefineIts;
    ctrlC.restart      = ctrl.restart;
    ctrlC.progress     = ctrl.progress;
    ctrlC.time         = ctrl.time;
    return ctrlC;
}

inline RegSolveCtrl<float> CReflect( const ElRegSolveCtrl_s& ctrlC )
{
    RegSolveCtrl<float> ctrl;
    ctrl.alg          = CReflect(ctrlC.alg);
    ctrl.relTol       = ctrlC.relTol;
    ctrl.relTolRefine = ctrlC.relTolRefine;
    ctrl.maxRefineIts = ctrlC.maxRefineIts;
    ctrl.restart      = ctrlC.restart;
    ctrl.progress     = ctrlC.progress;
    ctrl.time         = ctrlC.time;
    return ctrl;
}

inline RegSolveCtrl<double> CReflect( const ElRegSolveCtrl_d& ctrlC )
{
    RegSolveCtrl<double> ctrl;
    ctrl.alg          = CReflect(ctrlC.alg);
    ctrl.relTol       = ctrlC.relTol;
    ctrl.relTolRefine = ctrlC.relTolRefine;
    ctrl.maxRefineIts = ctrlC.maxRefineIts;
    ctrl.restart      = ctrlC.restart;
    ctrl.progress     = ctrlC.progress;
    ctrl.time         = ctrlC.time;
    return ctrl;
}

} // namespace El

#endif // ifndef EL_LAPACK_FACTOR_CREFLECT_C_HPP
