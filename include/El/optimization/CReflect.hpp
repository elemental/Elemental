/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_OPTIMIZATION_CREFLECT_C_HPP
#define EL_OPTIMIZATION_CREFLECT_C_HPP

namespace El {

// Optimization
// ------------
inline ElRegularization CReflect( Regularization penalty )
{ return static_cast<ElRegularization>(penalty); }
inline Regularization CReflect( ElRegularization penalty )
{ return static_cast<Regularization>(penalty); }

inline ElKKTSystem CReflect( KKTSystem system )
{ return static_cast<ElKKTSystem>(system); }
inline KKTSystem CReflect( ElKKTSystem system )
{ return static_cast<KKTSystem>(system); }

/* Mehrotra's Predictor-Corrector IPM
   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ */
inline ElMehrotraCtrl_s CReflect( const MehrotraCtrl<float>& ctrl )
{
    ElMehrotraCtrl_s ctrlC;
    ctrlC.primalInit    = ctrl.primalInit;
    ctrlC.dualInit      = ctrl.dualInit;
    ctrlC.minTol        = ctrl.minTol;
    ctrlC.targetTol     = ctrl.targetTol;
    ctrlC.maxIts        = ctrl.maxIts;
    ctrlC.maxStepRatio  = ctrl.maxStepRatio;
    ctrlC.system        = CReflect(ctrl.system);
    ctrlC.mehrotra      = ctrl.mehrotra;
    ctrlC.forceSameStep = ctrl.forceSameStep;
    ctrlC.solveCtrl     = CReflect(ctrl.solveCtrl);
    ctrlC.resolveReg    = ctrl.resolveReg;
    ctrlC.outerEquil    = ctrl.outerEquil;
    ctrlC.basisSize     = ctrl.basisSize;
    ctrlC.print         = ctrl.print;
    ctrlC.time          = ctrl.time;
    ctrlC.wSafeMaxNorm  = ctrl.wSafeMaxNorm;
    ctrlC.wMaxLimit     = ctrl.wMaxLimit;
    ctrlC.ruizEquilTol  = ctrl.ruizEquilTol;
    ctrlC.ruizMaxIter   = ctrl.ruizMaxIter;
    ctrlC.diagEquilTol  = ctrl.diagEquilTol;
    ctrlC.checkResiduals = ctrl.checkResiduals;

    ctrlC.reg0Tmp = ctrl.reg0Tmp;
    ctrlC.reg1Tmp = ctrl.reg1Tmp;
    ctrlC.reg2Tmp = ctrl.reg2Tmp;
    ctrlC.reg0Perm = ctrl.reg0Perm;
    ctrlC.reg1Perm = ctrl.reg1Perm;
    ctrlC.reg2Perm = ctrl.reg2Perm;

    return ctrlC;
}
inline ElMehrotraCtrl_d CReflect( const MehrotraCtrl<double>& ctrl )
{
    ElMehrotraCtrl_d ctrlC;
    ctrlC.primalInit    = ctrl.primalInit;
    ctrlC.dualInit      = ctrl.dualInit;
    ctrlC.minTol        = ctrl.minTol;
    ctrlC.targetTol     = ctrl.targetTol;
    ctrlC.maxIts        = ctrl.maxIts;
    ctrlC.maxStepRatio  = ctrl.maxStepRatio;
    ctrlC.system        = CReflect(ctrl.system);
    ctrlC.mehrotra      = ctrl.mehrotra;
    ctrlC.forceSameStep = ctrl.forceSameStep;
    ctrlC.solveCtrl     = CReflect(ctrl.solveCtrl);
    ctrlC.resolveReg    = ctrl.resolveReg;
    ctrlC.outerEquil    = ctrl.outerEquil;
    ctrlC.basisSize     = ctrl.basisSize;
    ctrlC.print         = ctrl.print;
    ctrlC.time          = ctrl.time;
    ctrlC.wSafeMaxNorm  = ctrl.wSafeMaxNorm;
    ctrlC.wMaxLimit     = ctrl.wMaxLimit;
    ctrlC.ruizEquilTol  = ctrl.ruizEquilTol;
    ctrlC.ruizMaxIter   = ctrl.ruizMaxIter;
    ctrlC.diagEquilTol  = ctrl.diagEquilTol;
    ctrlC.checkResiduals = ctrl.checkResiduals;

    ctrlC.reg0Tmp = ctrl.reg0Tmp;
    ctrlC.reg1Tmp = ctrl.reg1Tmp;
    ctrlC.reg2Tmp = ctrl.reg2Tmp;
    ctrlC.reg0Perm = ctrl.reg0Perm;
    ctrlC.reg1Perm = ctrl.reg1Perm;
    ctrlC.reg2Perm = ctrl.reg2Perm;

    return ctrlC;
}
inline MehrotraCtrl<float> CReflect( const ElMehrotraCtrl_s& ctrlC )
{
    MehrotraCtrl<float> ctrl;
    ctrl.primalInit    = ctrlC.primalInit;
    ctrl.dualInit      = ctrlC.dualInit;
    ctrl.minTol        = ctrlC.minTol;
    ctrl.targetTol     = ctrlC.targetTol;
    ctrl.maxIts        = ctrlC.maxIts;
    ctrl.maxStepRatio  = ctrlC.maxStepRatio;
    ctrl.system        = CReflect(ctrlC.system);
    ctrl.mehrotra      = ctrlC.mehrotra;
    ctrl.forceSameStep = ctrlC.forceSameStep;
    ctrl.solveCtrl     = CReflect(ctrlC.solveCtrl);
    ctrl.resolveReg    = ctrlC.resolveReg;
    ctrl.outerEquil    = ctrlC.outerEquil;
    ctrl.basisSize     = ctrlC.basisSize;
    ctrl.print         = ctrlC.print;
    ctrl.time          = ctrlC.time;
    ctrl.wSafeMaxNorm  = ctrlC.wSafeMaxNorm;
    ctrl.wMaxLimit     = ctrlC.wMaxLimit;
    ctrl.ruizEquilTol  = ctrlC.ruizEquilTol;
    ctrl.ruizMaxIter   = ctrlC.ruizMaxIter;
    ctrl.diagEquilTol  = ctrlC.diagEquilTol;
    ctrl.checkResiduals = ctrlC.checkResiduals;

    ctrl.reg0Tmp = ctrlC.reg0Tmp;
    ctrl.reg1Tmp = ctrlC.reg1Tmp;
    ctrl.reg2Tmp = ctrlC.reg2Tmp;
    ctrl.reg0Perm = ctrlC.reg0Perm;
    ctrl.reg1Perm = ctrlC.reg1Perm;
    ctrl.reg2Perm = ctrlC.reg2Perm;

    return ctrl;
}
inline MehrotraCtrl<double> CReflect( const ElMehrotraCtrl_d& ctrlC )
{
    MehrotraCtrl<double> ctrl;
    ctrl.primalInit    = ctrlC.primalInit;
    ctrl.dualInit      = ctrlC.dualInit;
    ctrl.minTol        = ctrlC.minTol;
    ctrl.targetTol     = ctrlC.targetTol;
    ctrl.maxIts        = ctrlC.maxIts;
    ctrl.maxStepRatio  = ctrlC.maxStepRatio;
    ctrl.system        = CReflect(ctrlC.system);
    ctrl.mehrotra      = ctrlC.mehrotra;
    ctrl.forceSameStep = ctrlC.forceSameStep;
    ctrl.solveCtrl     = CReflect(ctrlC.solveCtrl);
    ctrl.resolveReg    = ctrlC.resolveReg;
    ctrl.outerEquil    = ctrlC.outerEquil;
    ctrl.basisSize     = ctrlC.basisSize;
    ctrl.print         = ctrlC.print;
    ctrl.time          = ctrlC.time;
    ctrl.wSafeMaxNorm  = ctrlC.wSafeMaxNorm;
    ctrl.wMaxLimit     = ctrlC.wMaxLimit;
    ctrl.ruizEquilTol  = ctrlC.ruizEquilTol;
    ctrl.ruizMaxIter   = ctrlC.ruizMaxIter;
    ctrl.diagEquilTol  = ctrlC.diagEquilTol;
    ctrl.checkResiduals = ctrlC.checkResiduals;

    ctrl.reg0Tmp = ctrlC.reg0Tmp;
    ctrl.reg1Tmp = ctrlC.reg1Tmp;
    ctrl.reg2Tmp = ctrlC.reg2Tmp;
    ctrl.reg0Perm = ctrlC.reg0Perm;
    ctrl.reg1Perm = ctrlC.reg1Perm;
    ctrl.reg2Perm = ctrlC.reg2Perm;

    return ctrl;
}

/* Alternating Direction Method of Multipliers
   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ */
inline ElADMMCtrl_s CReflect( const ADMMCtrl<float>& ctrl )
{
    ElADMMCtrl_s ctrlC;
    ctrlC.rho     = ctrl.rho;
    ctrlC.alpha   = ctrl.alpha;
    ctrlC.maxIter = ctrl.maxIter;
    ctrlC.absTol  = ctrl.absTol;
    ctrlC.relTol  = ctrl.relTol;
    ctrlC.inv     = ctrl.inv;
    ctrlC.print   = ctrl.print;
    return ctrlC;
}
inline ElADMMCtrl_d CReflect( const ADMMCtrl<double>& ctrl )
{
    ElADMMCtrl_d ctrlC;
    ctrlC.rho     = ctrl.rho;
    ctrlC.alpha   = ctrl.alpha;
    ctrlC.maxIter = ctrl.maxIter;
    ctrlC.absTol  = ctrl.absTol;
    ctrlC.relTol  = ctrl.relTol;
    ctrlC.inv     = ctrl.inv;
    ctrlC.print   = ctrl.print;
    return ctrlC;
}
inline ADMMCtrl<float> CReflect( const ElADMMCtrl_s& ctrlC )
{
    ADMMCtrl<float> ctrl;
    ctrl.rho     = ctrlC.rho;
    ctrl.alpha   = ctrlC.alpha;
    ctrl.maxIter = ctrlC.maxIter;
    ctrl.absTol  = ctrlC.absTol;
    ctrl.relTol  = ctrlC.relTol;
    ctrl.inv     = ctrlC.inv;
    ctrl.print   = ctrlC.print;
    return ctrl;
}
inline ADMMCtrl<double> CReflect( const ElADMMCtrl_d& ctrlC )
{
    ADMMCtrl<double> ctrl;
    ctrl.rho     = ctrlC.rho;
    ctrl.alpha   = ctrlC.alpha;
    ctrl.maxIter = ctrlC.maxIter;
    ctrl.absTol  = ctrlC.absTol;
    ctrl.relTol  = ctrlC.relTol;
    ctrl.inv     = ctrlC.inv;
    ctrl.print   = ctrlC.print;
    return ctrl;
}

/* Linear programs
   ^^^^^^^^^^^^^^^ */
inline ElLPApproach CReflect( LPApproach approach )
{ return static_cast<ElLPApproach>(approach); }
inline LPApproach CReflect( ElLPApproach approach )
{ return static_cast<LPApproach>(approach); }

/* Direct conic form
   """"""""""""""""" */
inline ElLPDirectCtrl_s CReflect( const lp::direct::Ctrl<float>& ctrl )
{
    ElLPDirectCtrl_s ctrlC;
    ctrlC.approach     = CReflect(ctrl.approach);
    ctrlC.admmCtrl     = CReflect(ctrl.admmCtrl);
    ctrlC.mehrotraCtrl = CReflect(ctrl.mehrotraCtrl);
    return ctrlC;
}
inline ElLPDirectCtrl_d CReflect( const lp::direct::Ctrl<double>& ctrl )
{
    ElLPDirectCtrl_d ctrlC;
    ctrlC.approach     = CReflect(ctrl.approach);
    ctrlC.admmCtrl     = CReflect(ctrl.admmCtrl);
    ctrlC.mehrotraCtrl = CReflect(ctrl.mehrotraCtrl);
    return ctrlC;
}
inline lp::direct::Ctrl<float> CReflect( const ElLPDirectCtrl_s& ctrlC )
{
    lp::direct::Ctrl<float> ctrl(false);
    ctrl.approach     = CReflect(ctrlC.approach);
    ctrl.admmCtrl     = CReflect(ctrlC.admmCtrl);
    ctrl.mehrotraCtrl = CReflect(ctrlC.mehrotraCtrl);
    return ctrl;
}
inline lp::direct::Ctrl<double> CReflect( const ElLPDirectCtrl_d& ctrlC )
{
    lp::direct::Ctrl<double> ctrl(false);
    ctrl.approach     = CReflect(ctrlC.approach);
    ctrl.admmCtrl     = CReflect(ctrlC.admmCtrl);
    ctrl.mehrotraCtrl = CReflect(ctrlC.mehrotraCtrl);
    return ctrl;
}

/* Affine conic form
   """"""""""""""""" */
inline ElLPAffineCtrl_s CReflect( const lp::affine::Ctrl<float>& ctrl )
{
    ElLPAffineCtrl_s ctrlC;
    ctrlC.approach     = CReflect(ctrl.approach);
    ctrlC.mehrotraCtrl = CReflect(ctrl.mehrotraCtrl);
    return ctrlC;
}
inline ElLPAffineCtrl_d CReflect( const lp::affine::Ctrl<double>& ctrl )
{
    ElLPAffineCtrl_d ctrlC;
    ctrlC.approach     = CReflect(ctrl.approach);
    ctrlC.mehrotraCtrl = CReflect(ctrl.mehrotraCtrl);
    return ctrlC;
}
inline lp::affine::Ctrl<float> CReflect( const ElLPAffineCtrl_s& ctrlC )
{
    lp::affine::Ctrl<float> ctrl;
    ctrl.approach     = CReflect(ctrlC.approach);
    ctrl.mehrotraCtrl = CReflect(ctrlC.mehrotraCtrl);
    return ctrl;
}
inline lp::affine::Ctrl<double> CReflect( const ElLPAffineCtrl_d& ctrlC )
{
    lp::affine::Ctrl<double> ctrl;
    ctrl.approach     = CReflect(ctrlC.approach);
    ctrl.mehrotraCtrl = CReflect(ctrlC.mehrotraCtrl);
    return ctrl;
}

/* Quadratic programs
   ^^^^^^^^^^^^^^^^^^ */
inline ElQPApproach CReflect( QPApproach approach )
{ return static_cast<ElQPApproach>(approach); }
inline QPApproach CReflect( ElQPApproach approach )
{ return static_cast<QPApproach>(approach); }

/* Direct conic form
   """"""""""""""""" */
inline ElQPDirectCtrl_s CReflect( const qp::direct::Ctrl<float>& ctrl )
{
    ElQPDirectCtrl_s ctrlC;
    ctrlC.approach     = CReflect(ctrl.approach);
    ctrlC.mehrotraCtrl = CReflect(ctrl.mehrotraCtrl);
    return ctrlC;
}
inline ElQPDirectCtrl_d CReflect( const qp::direct::Ctrl<double>& ctrl )
{
    ElQPDirectCtrl_d ctrlC;
    ctrlC.approach     = CReflect(ctrl.approach);
    ctrlC.mehrotraCtrl = CReflect(ctrl.mehrotraCtrl);
    return ctrlC;
}
inline qp::direct::Ctrl<float> CReflect( const ElQPDirectCtrl_s& ctrlC )
{
    qp::direct::Ctrl<float> ctrl;
    ctrl.approach     = CReflect(ctrlC.approach);
    ctrl.mehrotraCtrl = CReflect(ctrlC.mehrotraCtrl);
    return ctrl;
}
inline qp::direct::Ctrl<double> CReflect( const ElQPDirectCtrl_d& ctrlC )
{
    qp::direct::Ctrl<double> ctrl;
    ctrl.approach     = CReflect(ctrlC.approach);
    ctrl.mehrotraCtrl = CReflect(ctrlC.mehrotraCtrl);
    return ctrl;
}

/* Affine conic form
   """"""""""""""""" */
inline ElQPAffineCtrl_s CReflect( const qp::affine::Ctrl<float>& ctrl )
{
    ElQPAffineCtrl_s ctrlC;
    ctrlC.approach     = CReflect(ctrl.approach);
    ctrlC.mehrotraCtrl = CReflect(ctrl.mehrotraCtrl);
    return ctrlC;
}
inline ElQPAffineCtrl_d CReflect( const qp::affine::Ctrl<double>& ctrl )
{
    ElQPAffineCtrl_d ctrlC;
    ctrlC.approach     = CReflect(ctrl.approach);
    ctrlC.mehrotraCtrl = CReflect(ctrl.mehrotraCtrl);
    return ctrlC;
}
inline qp::affine::Ctrl<float> CReflect( const ElQPAffineCtrl_s& ctrlC )
{
    qp::affine::Ctrl<float> ctrl;
    ctrl.approach     = CReflect(ctrlC.approach);
    ctrl.mehrotraCtrl = CReflect(ctrlC.mehrotraCtrl);
    return ctrl;
}
inline qp::affine::Ctrl<double> CReflect( const ElQPAffineCtrl_d& ctrlC )
{
    qp::affine::Ctrl<double> ctrl;
    ctrl.approach     = CReflect(ctrlC.approach);
    ctrl.mehrotraCtrl = CReflect(ctrlC.mehrotraCtrl);
    return ctrl;
}

/* Second-order cone programs
   ^^^^^^^^^^^^^^^^^^^^^^^^^^ */
inline ElSOCPApproach CReflect( SOCPApproach approach )
{ return static_cast<ElSOCPApproach>(approach); }
inline SOCPApproach CReflect( ElSOCPApproach approach )
{ return static_cast<SOCPApproach>(approach); }

/* Direct conic form
   """"""""""""""""" */
inline ElSOCPDirectCtrl_s CReflect( const socp::direct::Ctrl<float>& ctrl )
{
    ElSOCPDirectCtrl_s ctrlC;
    ctrlC.approach     = CReflect(ctrl.approach);
    ctrlC.mehrotraCtrl = CReflect(ctrl.mehrotraCtrl);
    return ctrlC;
}
inline ElSOCPDirectCtrl_d CReflect( const socp::direct::Ctrl<double>& ctrl )
{
    ElSOCPDirectCtrl_d ctrlC;
    ctrlC.approach     = CReflect(ctrl.approach);
    ctrlC.mehrotraCtrl = CReflect(ctrl.mehrotraCtrl);
    return ctrlC;
}
inline socp::direct::Ctrl<float> CReflect( const ElSOCPDirectCtrl_s& ctrlC )
{
    socp::direct::Ctrl<float> ctrl;
    ctrl.approach     = CReflect(ctrlC.approach);
    ctrl.mehrotraCtrl = CReflect(ctrlC.mehrotraCtrl);
    return ctrl;
}
inline socp::direct::Ctrl<double> CReflect( const ElSOCPDirectCtrl_d& ctrlC )
{
    socp::direct::Ctrl<double> ctrl;
    ctrl.approach     = CReflect(ctrlC.approach);
    ctrl.mehrotraCtrl = CReflect(ctrlC.mehrotraCtrl);
    return ctrl;
}

/* Affine conic form
   """"""""""""""""" */
inline ElSOCPAffineCtrl_s CReflect( const socp::affine::Ctrl<float>& ctrl )
{
    ElSOCPAffineCtrl_s ctrlC;
    ctrlC.approach     = CReflect(ctrl.approach);
    ctrlC.mehrotraCtrl = CReflect(ctrl.mehrotraCtrl);
    return ctrlC;
}
inline ElSOCPAffineCtrl_d CReflect( const socp::affine::Ctrl<double>& ctrl )
{
    ElSOCPAffineCtrl_d ctrlC;
    ctrlC.approach     = CReflect(ctrl.approach);
    ctrlC.mehrotraCtrl = CReflect(ctrl.mehrotraCtrl);
    return ctrlC;
}
inline socp::affine::Ctrl<float> CReflect( const ElSOCPAffineCtrl_s& ctrlC )
{
    socp::affine::Ctrl<float> ctrl;
    ctrl.approach     = CReflect(ctrlC.approach);
    ctrl.mehrotraCtrl = CReflect(ctrlC.mehrotraCtrl);
    return ctrl;
}
inline socp::affine::Ctrl<double> CReflect( const ElSOCPAffineCtrl_d& ctrlC )
{
    socp::affine::Ctrl<double> ctrl;
    ctrl.approach     = CReflect(ctrlC.approach);
    ctrl.mehrotraCtrl = CReflect(ctrlC.mehrotraCtrl);
    return ctrl;
}

// Models
// ^^^^^^

// Basis Pursuit
// """""""""""""
inline ElBPADMMCtrl_s CReflect( const bp::ADMMCtrl<float>& ctrl )
{
    ElBPADMMCtrl_s ctrlC;
    ctrlC.rho      = ctrl.rho;
    ctrlC.alpha    = ctrl.alpha;
    ctrlC.maxIter  = ctrl.maxIter;
    ctrlC.absTol   = ctrl.absTol;
    ctrlC.relTol   = ctrl.relTol;
    ctrlC.usePinv  = ctrl.usePinv;
    ctrlC.pinvTol  = ctrl.pinvTol;
    ctrlC.progress = ctrl.progress;
    return ctrlC;
}

inline ElBPADMMCtrl_d CReflect( const bp::ADMMCtrl<double>& ctrl )
{
    ElBPADMMCtrl_d ctrlC;
    ctrlC.rho      = ctrl.rho;
    ctrlC.alpha    = ctrl.alpha;
    ctrlC.maxIter  = ctrl.maxIter;
    ctrlC.absTol   = ctrl.absTol;
    ctrlC.relTol   = ctrl.relTol;
    ctrlC.usePinv  = ctrl.usePinv;
    ctrlC.pinvTol  = ctrl.pinvTol;
    ctrlC.progress = ctrl.progress;
    return ctrlC;
}

inline bp::ADMMCtrl<float> CReflect( const ElBPADMMCtrl_s& ctrlC )
{
    bp::ADMMCtrl<float> ctrl;
    ctrl.rho      = ctrlC.rho;
    ctrl.alpha    = ctrlC.alpha;
    ctrl.maxIter  = ctrlC.maxIter;
    ctrl.absTol   = ctrlC.absTol;
    ctrl.relTol   = ctrlC.relTol;
    ctrl.usePinv  = ctrlC.usePinv;
    ctrl.pinvTol  = ctrlC.pinvTol;
    ctrl.progress = ctrlC.progress;
    return ctrl;
}

inline bp::ADMMCtrl<double> CReflect( const ElBPADMMCtrl_d& ctrlC )
{
    bp::ADMMCtrl<double> ctrl;
    ctrl.rho      = ctrlC.rho;
    ctrl.alpha    = ctrlC.alpha;
    ctrl.maxIter  = ctrlC.maxIter;
    ctrl.absTol   = ctrlC.absTol;
    ctrl.relTol   = ctrlC.relTol;
    ctrl.usePinv  = ctrlC.usePinv;
    ctrl.pinvTol  = ctrlC.pinvTol;
    ctrl.progress = ctrlC.progress;
    return ctrl;
}

inline ElBPCtrl_s CReflect( const BPCtrl<float>& ctrl )
{
    ElBPCtrl_s ctrlC;
    ctrlC.useIPM      = ctrl.useIPM;
    ctrlC.useSOCP     = ctrl.useSOCP;
    ctrlC.admmCtrl    = CReflect(ctrl.admmCtrl);
    ctrlC.lpIPMCtrl   = CReflect(ctrl.lpIPMCtrl);
    ctrlC.socpIPMCtrl = CReflect(ctrl.socpIPMCtrl);
    return ctrlC;
}

inline ElBPCtrl_d CReflect( const BPCtrl<double>& ctrl )
{
    ElBPCtrl_d ctrlC;
    ctrlC.useIPM      = ctrl.useIPM;
    ctrlC.useSOCP     = ctrl.useSOCP;
    ctrlC.admmCtrl    = CReflect(ctrl.admmCtrl);
    ctrlC.lpIPMCtrl   = CReflect(ctrl.lpIPMCtrl);
    ctrlC.socpIPMCtrl = CReflect(ctrl.socpIPMCtrl);
    return ctrlC;
}

inline ElBPCtrl_c CReflect( const BPCtrl<Complex<float>>& ctrl )
{
    ElBPCtrl_c ctrlC;
    ctrlC.ipmCtrl = CReflect(ctrl.ipmCtrl);
    return ctrlC;
}

inline ElBPCtrl_z CReflect( const BPCtrl<Complex<double>>& ctrl )
{
    ElBPCtrl_z ctrlC;
    ctrlC.ipmCtrl = CReflect(ctrl.ipmCtrl);
    return ctrlC;
}

inline BPCtrl<float> CReflect( const ElBPCtrl_s& ctrlC )
{
    BPCtrl<float> ctrl(false);
    ctrl.useIPM      = ctrlC.useIPM;
    ctrl.useSOCP     = ctrlC.useSOCP;
    ctrl.admmCtrl    = CReflect(ctrlC.admmCtrl);
    ctrl.lpIPMCtrl   = CReflect(ctrlC.lpIPMCtrl);
    ctrl.socpIPMCtrl = CReflect(ctrlC.socpIPMCtrl);
    return ctrl;
}

inline BPCtrl<double> CReflect( const ElBPCtrl_d& ctrlC )
{
    BPCtrl<double> ctrl(false);
    ctrl.useIPM      = ctrlC.useIPM;
    ctrl.useSOCP     = ctrlC.useSOCP;
    ctrl.admmCtrl    = CReflect(ctrlC.admmCtrl);
    ctrl.lpIPMCtrl   = CReflect(ctrlC.lpIPMCtrl);
    ctrl.socpIPMCtrl = CReflect(ctrlC.socpIPMCtrl);
    return ctrl;
}

inline BPCtrl<Complex<float>> CReflect( const ElBPCtrl_c& ctrlC )
{
    BPCtrl<Complex<float>> ctrl;
    ctrl.ipmCtrl = CReflect(ctrlC.ipmCtrl);
    return ctrl;
}

inline BPCtrl<Complex<double>> CReflect( const ElBPCtrl_z& ctrlC )
{
    BPCtrl<Complex<double>> ctrl;
    ctrl.ipmCtrl = CReflect(ctrlC.ipmCtrl);
    return ctrl;
}


// BPDN / LASSO
// """"""""""""

inline ElBPDNADMMCtrl_s CReflect( const bpdn::ADMMCtrl<float>& ctrl )
{
    ElBPDNADMMCtrl_s ctrlC;
    ctrlC.rho      = ctrl.rho;
    ctrlC.alpha    = ctrl.alpha;
    ctrlC.maxIter  = ctrl.maxIter;
    ctrlC.absTol   = ctrl.absTol;
    ctrlC.relTol   = ctrl.relTol;
    ctrlC.inv      = ctrl.inv;
    ctrlC.progress = ctrl.progress;
    return ctrlC;
}

inline ElBPDNADMMCtrl_d CReflect( const bpdn::ADMMCtrl<double>& ctrl )
{
    ElBPDNADMMCtrl_d ctrlC;
    ctrlC.rho      = ctrl.rho;
    ctrlC.alpha    = ctrl.alpha;
    ctrlC.maxIter  = ctrl.maxIter;
    ctrlC.absTol   = ctrl.absTol;
    ctrlC.relTol   = ctrl.relTol;
    ctrlC.inv      = ctrl.inv;
    ctrlC.progress = ctrl.progress;
    return ctrlC;
}

inline bpdn::ADMMCtrl<float> CReflect( const ElBPDNADMMCtrl_s& ctrlC )
{
    bpdn::ADMMCtrl<float> ctrl;
    ctrl.rho      = ctrlC.rho;
    ctrl.alpha    = ctrlC.alpha;
    ctrl.maxIter  = ctrlC.maxIter;
    ctrl.absTol   = ctrlC.absTol;
    ctrl.relTol   = ctrlC.relTol;
    ctrl.inv      = ctrlC.inv;
    ctrl.progress = ctrlC.progress;
    return ctrl;
}

inline bpdn::ADMMCtrl<double> CReflect( const ElBPDNADMMCtrl_d& ctrlC )
{
    bpdn::ADMMCtrl<double> ctrl;
    ctrl.rho      = ctrlC.rho;
    ctrl.alpha    = ctrlC.alpha;
    ctrl.maxIter  = ctrlC.maxIter;
    ctrl.absTol   = ctrlC.absTol;
    ctrl.relTol   = ctrlC.relTol;
    ctrl.inv      = ctrlC.inv;
    ctrl.progress = ctrlC.progress;
    return ctrl;
}

inline ElBPDNCtrl_s CReflect( const BPDNCtrl<float>& ctrl )
{
    ElBPDNCtrl_s ctrlC;
    ctrlC.useIPM   = ctrl.useIPM;
    ctrlC.admmCtrl = CReflect(ctrl.admmCtrl);
    ctrlC.ipmCtrl  = CReflect(ctrl.ipmCtrl);
    return ctrlC;
}

inline ElBPDNCtrl_d CReflect( const BPDNCtrl<double>& ctrl )
{
    ElBPDNCtrl_d ctrlC;
    ctrlC.useIPM   = ctrl.useIPM;
    ctrlC.admmCtrl = CReflect(ctrl.admmCtrl);
    ctrlC.ipmCtrl  = CReflect(ctrl.ipmCtrl);
    return ctrlC;
}

inline BPDNCtrl<float> CReflect( const ElBPDNCtrl_s& ctrlC )
{
    BPDNCtrl<float> ctrl;
    ctrl.useIPM   = ctrlC.useIPM;
    ctrl.admmCtrl = CReflect(ctrlC.admmCtrl);
    ctrl.ipmCtrl  = CReflect(ctrlC.ipmCtrl);
    return ctrl;
}

inline BPDNCtrl<double> CReflect( const ElBPDNCtrl_d& ctrlC )
{
    BPDNCtrl<double> ctrl;
    ctrl.useIPM   = ctrlC.useIPM;
    ctrl.admmCtrl = CReflect(ctrlC.admmCtrl);
    ctrl.ipmCtrl  = CReflect(ctrlC.ipmCtrl);
    return ctrl;
}

// Non-negative Least Squares
// """"""""""""""""""""""""""

inline NNLSApproach CReflect( ElNNLSApproach approach ) 
{ return static_cast<NNLSApproach>(approach); }
inline ElNNLSApproach CReflect( NNLSApproach approach )
{ return static_cast<ElNNLSApproach>(approach); }

inline ElNNLSCtrl_s CReflect( const NNLSCtrl<float>& ctrl )
{
    ElNNLSCtrl_s ctrlC;
    ctrlC.approach = CReflect(ctrl.approach);
    ctrlC.admmCtrl = CReflect(ctrl.admmCtrl);
    ctrlC.qpCtrl   = CReflect(ctrl.qpCtrl);
    ctrlC.socpCtrl = CReflect(ctrl.socpCtrl);
    return ctrlC;
}

inline ElNNLSCtrl_d CReflect( const NNLSCtrl<double>& ctrl )
{
    ElNNLSCtrl_d ctrlC;
    ctrlC.approach = CReflect(ctrl.approach);
    ctrlC.admmCtrl = CReflect(ctrl.admmCtrl);
    ctrlC.qpCtrl   = CReflect(ctrl.qpCtrl);
    ctrlC.socpCtrl = CReflect(ctrl.socpCtrl);
    return ctrlC;
}

inline NNLSCtrl<float> CReflect( const ElNNLSCtrl_s& ctrlC )
{
    NNLSCtrl<float> ctrl;
    ctrl.approach = CReflect(ctrlC.approach);
    ctrl.admmCtrl = CReflect(ctrlC.admmCtrl);
    ctrl.qpCtrl   = CReflect(ctrlC.qpCtrl);
    ctrl.socpCtrl = CReflect(ctrlC.socpCtrl);
    return ctrl;
}

inline NNLSCtrl<double> CReflect( const ElNNLSCtrl_d& ctrlC )
{
    NNLSCtrl<double> ctrl;
    ctrl.approach = CReflect(ctrlC.approach);
    ctrl.admmCtrl = CReflect(ctrlC.admmCtrl);
    ctrl.qpCtrl   = CReflect(ctrlC.qpCtrl);
    ctrl.socpCtrl = CReflect(ctrlC.socpCtrl);
    return ctrl;
}

// Non-negative Matrix Factorization
// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
inline ElNMFCtrl_s CReflect( const NMFCtrl<float>& ctrl )
{
    ElNMFCtrl_s ctrlC;
    ctrlC.nnlsCtrl = CReflect(ctrl.nnlsCtrl);
    ctrlC.maxIter = ctrl.maxIter;
    return ctrlC;
}

inline ElNMFCtrl_d CReflect( const NMFCtrl<double>& ctrl )
{
    ElNMFCtrl_d ctrlC;
    ctrlC.nnlsCtrl = CReflect(ctrl.nnlsCtrl);
    ctrlC.maxIter = ctrl.maxIter;
    return ctrlC;
}

inline NMFCtrl<float> CReflect( const ElNMFCtrl_s& ctrlC )
{
    NMFCtrl<float> ctrl;
    ctrl.nnlsCtrl = CReflect(ctrlC.nnlsCtrl);
    ctrl.maxIter = ctrlC.maxIter;
    return ctrl;
}

inline NMFCtrl<double> CReflect( const ElNMFCtrl_d& ctrlC )
{
    NMFCtrl<double> ctrl;
    ctrl.nnlsCtrl = CReflect(ctrlC.nnlsCtrl);
    ctrl.maxIter = ctrlC.maxIter;
    return ctrl;
}

// Robust Principal Component Analysis
// """""""""""""""""""""""""""""""""""

inline ElRPCACtrl_s CReflect( const RPCACtrl<float>& ctrl )
{
    ElRPCACtrl_s ctrlC;
    ctrlC.useALM      = ctrl.useALM;
    ctrlC.usePivQR    = ctrl.usePivQR;
    ctrlC.progress    = ctrl.progress;
    ctrlC.numPivSteps = ctrl.numPivSteps;
    ctrlC.maxIts      = ctrl.maxIts;
    ctrlC.tau         = ctrl.tau;
    ctrlC.beta        = ctrl.beta;
    ctrlC.rho         = ctrl.rho;
    ctrlC.tol         = ctrl.tol;
    return ctrlC;
}

inline ElRPCACtrl_d CReflect( const RPCACtrl<double>& ctrl )
{
    ElRPCACtrl_d ctrlC;
    ctrlC.useALM      = ctrl.useALM;
    ctrlC.usePivQR    = ctrl.usePivQR;
    ctrlC.progress    = ctrl.progress;
    ctrlC.numPivSteps = ctrl.numPivSteps;
    ctrlC.maxIts      = ctrl.maxIts;
    ctrlC.tau         = ctrl.tau;
    ctrlC.beta        = ctrl.beta;
    ctrlC.rho         = ctrl.rho;
    ctrlC.tol         = ctrl.tol;
    return ctrlC;
}

inline RPCACtrl<float> CReflect( const ElRPCACtrl_s& ctrlC )
{
    RPCACtrl<float> ctrl;
    ctrl.useALM      = ctrlC.useALM;
    ctrl.usePivQR    = ctrlC.usePivQR;
    ctrl.progress    = ctrlC.progress;
    ctrl.numPivSteps = ctrlC.numPivSteps;
    ctrl.maxIts      = ctrlC.maxIts;
    ctrl.tau         = ctrlC.tau;
    ctrl.beta        = ctrlC.beta;
    ctrl.rho         = ctrlC.rho;
    ctrl.tol         = ctrlC.tol;
    return ctrl;
}

inline RPCACtrl<double> CReflect( const ElRPCACtrl_d& ctrlC )
{
    RPCACtrl<double> ctrl;
    ctrl.useALM      = ctrlC.useALM;
    ctrl.usePivQR    = ctrlC.usePivQR;
    ctrl.progress    = ctrlC.progress;
    ctrl.numPivSteps = ctrlC.numPivSteps;
    ctrl.maxIts      = ctrlC.maxIts;
    ctrl.tau         = ctrlC.tau;
    ctrl.beta        = ctrlC.beta;
    ctrl.rho         = ctrlC.rho;
    ctrl.tol         = ctrlC.tol;
    return ctrl;
}

// Sparse inverse covariance selection
// """""""""""""""""""""""""""""""""""

inline ElSparseInvCovCtrl_s CReflect( const SparseInvCovCtrl<float>& ctrl )
{
    ElSparseInvCovCtrl_s ctrlC;
    ctrlC.rho      = ctrl.rho;
    ctrlC.alpha    = ctrl.alpha;
    ctrlC.maxIter  = ctrl.maxIter;
    ctrlC.absTol   = ctrl.absTol;
    ctrlC.relTol   = ctrl.relTol;
    ctrlC.progress = ctrl.progress;
    return ctrlC;
}

inline ElSparseInvCovCtrl_d CReflect( const SparseInvCovCtrl<double>& ctrl )
{
    ElSparseInvCovCtrl_d ctrlC;
    ctrlC.rho      = ctrl.rho;
    ctrlC.alpha    = ctrl.alpha;
    ctrlC.maxIter  = ctrl.maxIter;
    ctrlC.absTol   = ctrl.absTol;
    ctrlC.relTol   = ctrl.relTol;
    ctrlC.progress = ctrl.progress;
    return ctrlC;
}

inline SparseInvCovCtrl<float> CReflect( const ElSparseInvCovCtrl_s& ctrlC )
{
    SparseInvCovCtrl<float> ctrl;
    ctrl.rho      = ctrlC.rho;
    ctrl.alpha    = ctrlC.alpha;
    ctrl.maxIter  = ctrlC.maxIter;
    ctrl.absTol   = ctrlC.absTol;
    ctrl.relTol   = ctrlC.relTol;
    ctrl.progress = ctrlC.progress;
    return ctrl;
}

inline SparseInvCovCtrl<double> CReflect( const ElSparseInvCovCtrl_d& ctrlC )
{
    SparseInvCovCtrl<double> ctrl;
    ctrl.rho      = ctrlC.rho;
    ctrl.alpha    = ctrlC.alpha;
    ctrl.maxIter  = ctrlC.maxIter;
    ctrl.absTol   = ctrlC.absTol;
    ctrl.relTol   = ctrlC.relTol;
    ctrl.progress = ctrlC.progress;
    return ctrl;
}

/* Model Fit
   """"""""" */
inline ElModelFitCtrl_s CReflect( const ModelFitCtrl<float>& ctrl )
{
    ElModelFitCtrl_s ctrlC;
    ctrlC.rho      = ctrl.rho;
    ctrlC.maxIter  = ctrl.maxIter;
    ctrlC.inv      = ctrl.inv;
    ctrlC.progress = ctrl.progress;
    return ctrlC;
}

inline ElModelFitCtrl_d CReflect( const ModelFitCtrl<double>& ctrl )
{
    ElModelFitCtrl_d ctrlC;
    ctrlC.rho      = ctrl.rho;
    ctrlC.maxIter  = ctrl.maxIter;
    ctrlC.inv      = ctrl.inv;
    ctrlC.progress = ctrl.progress;
    return ctrlC;
}

inline ModelFitCtrl<float> CReflect( const ElModelFitCtrl_s& ctrlC )
{
    ModelFitCtrl<float> ctrl;
    ctrl.rho      = ctrlC.rho;
    ctrl.maxIter  = ctrlC.maxIter;
    ctrl.inv      = ctrlC.inv;
    ctrl.progress = ctrlC.progress;
    return ctrl;
}

inline ModelFitCtrl<double> CReflect( const ElModelFitCtrl_d& ctrlC )
{
    ModelFitCtrl<double> ctrl;
    ctrl.rho      = ctrlC.rho;
    ctrl.maxIter  = ctrlC.maxIter;
    ctrl.inv      = ctrlC.inv;
    ctrl.progress = ctrlC.progress;
    return ctrl;
}

/* Support Vector Machine
   """""""""""""""""""""" */
inline ElSVMCtrl_s CReflect( const SVMCtrl<float>& ctrl )
{
    ElSVMCtrl_s ctrlC;
    ctrlC.useIPM       = ctrl.useIPM;
    ctrlC.modelFitCtrl = CReflect(ctrl.modelFitCtrl);
    ctrlC.ipmCtrl      = CReflect(ctrl.ipmCtrl);
    return ctrlC;
}

inline ElSVMCtrl_d CReflect( const SVMCtrl<double>& ctrl )
{
    ElSVMCtrl_d ctrlC;
    ctrlC.useIPM       = ctrl.useIPM;
    ctrlC.modelFitCtrl = CReflect(ctrl.modelFitCtrl);
    ctrlC.ipmCtrl      = CReflect(ctrl.ipmCtrl);
    return ctrlC;
}

inline SVMCtrl<float> CReflect( const ElSVMCtrl_s& ctrlC )
{
    SVMCtrl<float> ctrl;
    ctrl.useIPM       = ctrlC.useIPM;
    ctrl.modelFitCtrl = CReflect(ctrlC.modelFitCtrl);
    ctrl.ipmCtrl      = CReflect(ctrlC.ipmCtrl);
    return ctrl;
}

inline SVMCtrl<double> CReflect( const ElSVMCtrl_d& ctrlC )
{
    SVMCtrl<double> ctrl;
    ctrl.useIPM       = ctrlC.useIPM;
    ctrl.modelFitCtrl = CReflect(ctrlC.modelFitCtrl);
    ctrl.ipmCtrl      = CReflect(ctrlC.ipmCtrl);
    return ctrl;
}

} // namespace El

#endif // ifndef EL_OPTIMIZATION_CREFLECT_C_HPP
