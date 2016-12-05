/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_LAPACK_SPECTRAL_CREFLECT_C_HPP
#define EL_LAPACK_SPECTRAL_CREFLECT_C_HPP

namespace El {

/* Cubic secular */
inline ElFlipOrClip CReflect( FlipOrClip negativeFix )
{ return static_cast<ElFlipOrClip>(negativeFix); }

inline FlipOrClip CReflect( ElFlipOrClip negativeFix )
{ return static_cast<FlipOrClip>(negativeFix); }

inline ElCubicSecularCtrl CReflect( const CubicSecularCtrl& ctrl )
{
    ElCubicSecularCtrl ctrlC;
    ctrlC.maxIterations = ctrl.maxIterations;
    ctrlC.negativeFix = CReflect( ctrl.negativeFix );
    return ctrlC;
}

inline CubicSecularCtrl CReflect( const ElCubicSecularCtrl& ctrlC )
{
    CubicSecularCtrl ctrl;
    ctrl.maxIterations = ctrlC.maxIterations;
    ctrl.negativeFix = CReflect( ctrlC.negativeFix );
    return ctrl;
}

/* Secular EVD */
inline ElSecularEVDCtrl_s
CReflect( const SecularEVDCtrl<float>& ctrl )
{
    ElSecularEVDCtrl_s ctrlC;
    ctrlC.maxIterations = ctrl.maxIterations;
    ctrlC.sufficientDecay = ctrl.sufficientDecay;
    ctrlC.negativeFix = CReflect( ctrl.negativeFix );
    ctrlC.penalizeDerivative = ctrl.penalizeDerivative;
    ctrlC.progress = ctrl.progress;
    ctrlC.cubicCtrl = CReflect( ctrl.cubicCtrl );
    return ctrlC;
}

inline ElSecularEVDCtrl_d
CReflect( const SecularEVDCtrl<double>& ctrl )
{
    ElSecularEVDCtrl_d ctrlC;
    ctrlC.maxIterations = ctrl.maxIterations;
    ctrlC.sufficientDecay = ctrl.sufficientDecay;
    ctrlC.negativeFix = CReflect( ctrl.negativeFix );
    ctrlC.penalizeDerivative = ctrl.penalizeDerivative;
    ctrlC.progress = ctrl.progress;
    ctrlC.cubicCtrl = CReflect( ctrl.cubicCtrl );
    return ctrlC;
}

inline SecularEVDCtrl<float>
CReflect( const ElSecularEVDCtrl_s& ctrlC )
{
    SecularEVDCtrl<float> ctrl;
    ctrl.maxIterations = ctrlC.maxIterations;
    ctrl.sufficientDecay = ctrlC.sufficientDecay;
    ctrl.negativeFix = CReflect( ctrlC.negativeFix );
    ctrl.penalizeDerivative = ctrlC.penalizeDerivative;
    ctrl.progress = ctrlC.progress;
    ctrl.cubicCtrl = CReflect( ctrlC.cubicCtrl );
    return ctrl;
}

inline SecularEVDCtrl<double>
CReflect( const ElSecularEVDCtrl_d& ctrlC )
{
    SecularEVDCtrl<double> ctrl;
    ctrl.maxIterations = ctrlC.maxIterations;
    ctrl.sufficientDecay = ctrlC.sufficientDecay;
    ctrl.negativeFix = CReflect( ctrlC.negativeFix );
    ctrl.penalizeDerivative = ctrlC.penalizeDerivative;
    ctrl.progress = ctrlC.progress;
    ctrl.cubicCtrl = CReflect( ctrlC.cubicCtrl );
    return ctrl;
}

/* Secular SVD */
inline ElSecularSVDCtrl_s
CReflect( const SecularSVDCtrl<float>& ctrl )
{
    ElSecularSVDCtrl_s ctrlC;
    ctrlC.maxIterations = ctrl.maxIterations;
    ctrlC.sufficientDecay = ctrl.sufficientDecay;
    ctrlC.negativeFix = CReflect( ctrl.negativeFix );
    ctrlC.penalizeDerivative = ctrl.penalizeDerivative;
    ctrlC.progress = ctrl.progress;
    ctrlC.cubicCtrl = CReflect( ctrl.cubicCtrl );
    return ctrlC;
}

inline ElSecularSVDCtrl_d
CReflect( const SecularSVDCtrl<double>& ctrl )
{
    ElSecularSVDCtrl_d ctrlC;
    ctrlC.maxIterations = ctrl.maxIterations;
    ctrlC.sufficientDecay = ctrl.sufficientDecay;
    ctrlC.negativeFix = CReflect( ctrl.negativeFix );
    ctrlC.penalizeDerivative = ctrl.penalizeDerivative;
    ctrlC.progress = ctrl.progress;
    ctrlC.cubicCtrl = CReflect( ctrl.cubicCtrl );
    return ctrlC;
}

inline SecularSVDCtrl<float>
CReflect( const ElSecularSVDCtrl_s& ctrlC )
{
    SecularSVDCtrl<float> ctrl;
    ctrl.maxIterations = ctrlC.maxIterations;
    ctrl.sufficientDecay = ctrlC.sufficientDecay;
    ctrl.negativeFix = CReflect( ctrlC.negativeFix );
    ctrl.penalizeDerivative = ctrlC.penalizeDerivative;
    ctrl.progress = ctrlC.progress;
    ctrl.cubicCtrl = CReflect( ctrlC.cubicCtrl );
    return ctrl;
}

inline SecularSVDCtrl<double>
CReflect( const ElSecularSVDCtrl_d& ctrlC )
{
    SecularSVDCtrl<double> ctrl;
    ctrl.maxIterations = ctrlC.maxIterations;
    ctrl.sufficientDecay = ctrlC.sufficientDecay;
    ctrl.negativeFix = CReflect( ctrlC.negativeFix );
    ctrl.penalizeDerivative = ctrlC.penalizeDerivative;
    ctrl.progress = ctrlC.progress;
    ctrl.cubicCtrl = CReflect( ctrlC.cubicCtrl );
    return ctrl;
}

/* HermitianEig */

inline ElHermitianTridiagEigAlg CReflect
( const HermitianTridiagEigAlg& alg )
{ return static_cast<ElHermitianTridiagEigAlg>(alg); }
inline HermitianTridiagEigAlg CReflect
( const ElHermitianTridiagEigAlg& alg )
{ return static_cast<HermitianTridiagEigAlg>(alg); }

/* HermitianEigSubset */
inline ElHermitianEigSubset_s CReflect
( const HermitianEigSubset<float>& subset )
{
    ElHermitianEigSubset_s subsetC;
    subsetC.indexSubset = subset.indexSubset;
    subsetC.lowerIndex = subset.lowerIndex;
    subsetC.upperIndex = subset.upperIndex;
    subsetC.rangeSubset = subset.rangeSubset;
    subsetC.lowerBound = subset.lowerBound;
    subsetC.upperBound = subset.upperBound;
    return subsetC;
}
inline ElHermitianEigSubset_d CReflect
( const HermitianEigSubset<double>& subset )
{
    ElHermitianEigSubset_d subsetC;
    subsetC.indexSubset = subset.indexSubset;
    subsetC.lowerIndex = subset.lowerIndex;
    subsetC.upperIndex = subset.upperIndex;
    subsetC.rangeSubset = subset.rangeSubset;
    subsetC.lowerBound = subset.lowerBound;
    subsetC.upperBound = subset.upperBound;
    return subsetC;
}

inline HermitianEigSubset<float> CReflect
( const ElHermitianEigSubset_s& subsetC )
{
    HermitianEigSubset<float> subset;
    subset.indexSubset = subsetC.indexSubset;
    subset.lowerIndex = subsetC.lowerIndex;
    subset.upperIndex = subsetC.upperIndex;
    subset.rangeSubset = subsetC.rangeSubset;
    subset.lowerBound = subsetC.lowerBound;
    subset.upperBound = subsetC.upperBound;
    return subset;
}
inline HermitianEigSubset<double> CReflect
( const ElHermitianEigSubset_d& subsetC )
{
    HermitianEigSubset<double> subset;
    subset.indexSubset = subsetC.indexSubset;
    subset.lowerIndex = subsetC.lowerIndex;
    subset.upperIndex = subsetC.upperIndex;
    subset.rangeSubset = subsetC.rangeSubset;
    subset.lowerBound = subsetC.lowerBound;
    subset.upperBound = subsetC.upperBound;
    return subset;
}

/* herm_tridiag_eig::QRCtrl */
inline ElHermitianTridiagEigQRCtrl CReflect
( const herm_tridiag_eig::QRCtrl& ctrl )
{
    ElHermitianTridiagEigQRCtrl ctrlC;
    ctrlC.maxIterPerEig = ctrl.maxIterPerEig;
    ctrlC.demandConverged = ctrl.demandConverged;
    ctrlC.fullAccuracyTwoByTwo = ctrl.fullAccuracyTwoByTwo;
    ctrlC.broadcast = ctrl.broadcast;
    return ctrlC;
}

inline herm_tridiag_eig::QRCtrl CReflect
( const ElHermitianTridiagEigQRCtrl& ctrlC )
{
    herm_tridiag_eig::QRCtrl ctrl;
    ctrl.maxIterPerEig = ctrlC.maxIterPerEig;
    ctrl.demandConverged = ctrlC.demandConverged;
    ctrl.fullAccuracyTwoByTwo = ctrlC.fullAccuracyTwoByTwo;
    ctrl.broadcast = ctrlC.broadcast;
    return ctrl;
}

/* herm_tridiag_eig::DCCtrl */
inline herm_tridiag_eig::DCCtrl<float> CReflect
( const ElHermitianTridiagEigDCCtrl_s& ctrlC )
{
    herm_tridiag_eig::DCCtrl<float> ctrl;
    ctrl.secularCtrl = CReflect( ctrlC.secularCtrl );
    ctrl.deflationFudge = ctrlC.deflationFudge;
    ctrl.cutoff = ctrlC.cutoff;
    ctrl.exploitStructure = ctrlC.exploitStructure;
    return ctrl;
}

inline herm_tridiag_eig::DCCtrl<double> CReflect
( const ElHermitianTridiagEigDCCtrl_d& ctrlC )
{
    herm_tridiag_eig::DCCtrl<double> ctrl;
    ctrl.secularCtrl = CReflect( ctrlC.secularCtrl );
    ctrl.deflationFudge = ctrlC.deflationFudge;
    ctrl.cutoff = ctrlC.cutoff;
    ctrl.exploitStructure = ctrlC.exploitStructure;
    return ctrl;
}

inline ElHermitianTridiagEigDCCtrl_s CReflect
( const herm_tridiag_eig::DCCtrl<float>& ctrl )
{
    ElHermitianTridiagEigDCCtrl_s ctrlC;
    ctrlC.secularCtrl = CReflect( ctrl.secularCtrl );
    ctrlC.deflationFudge = ctrl.deflationFudge;
    ctrlC.cutoff = ctrl.cutoff;
    ctrlC.exploitStructure = ctrl.exploitStructure;
    return ctrlC;
}
inline ElHermitianTridiagEigDCCtrl_d CReflect
( const herm_tridiag_eig::DCCtrl<double>& ctrl )
{
    ElHermitianTridiagEigDCCtrl_d ctrlC;
    ctrlC.secularCtrl = CReflect( ctrl.secularCtrl );
    ctrlC.deflationFudge = ctrl.deflationFudge;
    ctrlC.cutoff = ctrl.cutoff;
    ctrlC.exploitStructure = ctrl.exploitStructure;
    return ctrlC;
}

/* HermitianTridiagEigCtrl */
inline ElHermitianTridiagEigCtrl_s CReflect
( const HermitianTridiagEigCtrl<float>& ctrl )
{
    ElHermitianTridiagEigCtrl_s ctrlC;
    ctrlC.wantEigVecs = ctrl.wantEigVecs;
    ctrlC.accumulateEigVecs = ctrl.accumulateEigVecs;
    ctrlC.sort = CReflect(ctrl.sort);
    ctrlC.subset = CReflect(ctrl.subset);
    ctrlC.progress = ctrl.progress;
    ctrlC.alg = CReflect(ctrl.alg);
    ctrlC.qrCtrl = CReflect(ctrl.qrCtrl);
    ctrlC.dcCtrl = CReflect(ctrl.dcCtrl);
    return ctrlC;
}

inline ElHermitianTridiagEigCtrl_d CReflect
( const HermitianTridiagEigCtrl<double>& ctrl )
{
    ElHermitianTridiagEigCtrl_d ctrlC;
    ctrlC.wantEigVecs = ctrl.wantEigVecs;
    ctrlC.accumulateEigVecs = ctrl.accumulateEigVecs;
    ctrlC.sort = CReflect(ctrl.sort);
    ctrlC.subset = CReflect(ctrl.subset);
    ctrlC.progress = ctrl.progress;
    ctrlC.alg = CReflect(ctrl.alg);
    ctrlC.qrCtrl = CReflect(ctrl.qrCtrl);
    ctrlC.dcCtrl = CReflect(ctrl.dcCtrl);
    return ctrlC;
}

inline HermitianTridiagEigCtrl<float> CReflect
( const ElHermitianTridiagEigCtrl_s& ctrlC )
{
    HermitianTridiagEigCtrl<float> ctrl;
    ctrl.wantEigVecs = ctrlC.wantEigVecs;
    ctrl.accumulateEigVecs = ctrlC.accumulateEigVecs;
    ctrl.sort = CReflect(ctrlC.sort);
    ctrl.subset = CReflect(ctrlC.subset);
    ctrl.progress = ctrlC.progress;
    ctrl.alg = CReflect(ctrlC.alg);
    ctrl.qrCtrl = CReflect(ctrlC.qrCtrl);
    ctrl.dcCtrl = CReflect(ctrlC.dcCtrl);
    return ctrl;
}

inline HermitianTridiagEigCtrl<double> CReflect
( const ElHermitianTridiagEigCtrl_d& ctrlC )
{
    HermitianTridiagEigCtrl<double> ctrl;
    ctrl.wantEigVecs = ctrlC.wantEigVecs;
    ctrl.accumulateEigVecs = ctrlC.accumulateEigVecs;
    ctrl.sort = CReflect(ctrlC.sort);
    ctrl.subset = CReflect(ctrlC.subset);
    ctrl.progress = ctrlC.progress;
    ctrl.alg = CReflect(ctrlC.alg);
    ctrl.qrCtrl = CReflect(ctrlC.qrCtrl);
    ctrl.dcCtrl = CReflect(ctrlC.dcCtrl);
    return ctrl;
}

/* Pencil */
inline ElPencil CReflect( Pencil pencil )
{ return static_cast<ElPencil>(pencil); }

inline Pencil CReflect( ElPencil pencil )
{ return static_cast<Pencil>(pencil); }

/* HermitianSDCCtrl */
inline ElHermitianSDCCtrl_s CReflect( const HermitianSDCCtrl<float>& ctrl )
{
    ElHermitianSDCCtrl_s ctrlC;
    ctrlC.cutoff = ctrl.cutoff;
    ctrlC.maxInnerIts = ctrl.maxInnerIts;
    ctrlC.maxOuterIts = ctrl.maxOuterIts;
    ctrlC.tol = ctrl.tol;
    ctrlC.spreadFactor = ctrl.spreadFactor;
    ctrlC.progress = ctrl.progress;
    return ctrlC;
}
inline ElHermitianSDCCtrl_d CReflect( const HermitianSDCCtrl<double>& ctrl )
{
    ElHermitianSDCCtrl_d ctrlC;
    ctrlC.cutoff = ctrl.cutoff;
    ctrlC.maxInnerIts = ctrl.maxInnerIts;
    ctrlC.maxOuterIts = ctrl.maxOuterIts;
    ctrlC.tol = ctrl.tol;
    ctrlC.spreadFactor = ctrl.spreadFactor;
    ctrlC.progress = ctrl.progress;
    return ctrlC;
}

inline HermitianSDCCtrl<float> CReflect( const ElHermitianSDCCtrl_s& ctrlC )
{
    HermitianSDCCtrl<float> ctrl;
    ctrl.cutoff = ctrlC.cutoff;
    ctrl.maxInnerIts = ctrlC.maxInnerIts;
    ctrl.maxOuterIts = ctrlC.maxOuterIts;
    ctrl.tol = ctrlC.tol;
    ctrl.spreadFactor = ctrlC.spreadFactor;
    ctrl.progress = ctrlC.progress;
    return ctrl;
}
inline HermitianSDCCtrl<double> CReflect( const ElHermitianSDCCtrl_d& ctrlC )
{
    HermitianSDCCtrl<double> ctrl;
    ctrl.cutoff = ctrlC.cutoff;
    ctrl.maxInnerIts = ctrlC.maxInnerIts;
    ctrl.maxOuterIts = ctrlC.maxOuterIts;
    ctrl.tol = ctrlC.tol;
    ctrl.spreadFactor = ctrlC.spreadFactor;
    ctrl.progress = ctrlC.progress;
    return ctrl;
}

/* HermitianEigCtrl */
inline ElHermitianEigCtrl_s CReflect( const HermitianEigCtrl<float>& ctrl )
{
    ElHermitianEigCtrl_s ctrlC;
    ctrlC.tridiagCtrl = CReflect( ctrl.tridiagCtrl );
    ctrlC.tridiagEigCtrl = CReflect( ctrl.tridiagEigCtrl );
    ctrlC.sdcCtrl = CReflect( ctrl.sdcCtrl );
    ctrlC.useSDC = ctrl.useSDC;
    return ctrlC;
}

inline ElHermitianEigCtrl_d CReflect( const HermitianEigCtrl<double>& ctrl )
{
    ElHermitianEigCtrl_d ctrlC;
    ctrlC.tridiagCtrl = CReflect( ctrl.tridiagCtrl );
    ctrlC.tridiagEigCtrl = CReflect( ctrl.tridiagEigCtrl );
    ctrlC.sdcCtrl = CReflect( ctrl.sdcCtrl );
    ctrlC.useSDC = ctrl.useSDC;
    return ctrlC;
}
inline ElHermitianEigCtrl_c
CReflect( const HermitianEigCtrl<Complex<float>>& ctrl )
{
    ElHermitianEigCtrl_c ctrlC;
    ctrlC.tridiagCtrl = CReflect( ctrl.tridiagCtrl );
    ctrlC.tridiagEigCtrl = CReflect( ctrl.tridiagEigCtrl );
    ctrlC.sdcCtrl = CReflect( ctrl.sdcCtrl );
    ctrlC.useSDC = ctrl.useSDC;
    return ctrlC;
}
inline ElHermitianEigCtrl_z
CReflect( const HermitianEigCtrl<Complex<double>>& ctrl )
{
    ElHermitianEigCtrl_z ctrlC;
    ctrlC.tridiagCtrl = CReflect( ctrl.tridiagCtrl );
    ctrlC.tridiagEigCtrl = CReflect( ctrl.tridiagEigCtrl );
    ctrlC.sdcCtrl = CReflect( ctrl.sdcCtrl );
    ctrlC.useSDC = ctrl.useSDC;
    return ctrlC;
}

inline HermitianEigCtrl<float> CReflect( const ElHermitianEigCtrl_s& ctrlC )
{
    HermitianEigCtrl<float> ctrl;
    ctrl.tridiagCtrl = CReflect<float>( ctrlC.tridiagCtrl );
    ctrl.tridiagEigCtrl = CReflect( ctrlC.tridiagEigCtrl );
    ctrl.sdcCtrl = CReflect( ctrlC.sdcCtrl );
    ctrl.useSDC = ctrlC.useSDC;
    return ctrl;
}
inline HermitianEigCtrl<double> CReflect( const ElHermitianEigCtrl_d& ctrlC )
{
    HermitianEigCtrl<double> ctrl;
    ctrl.tridiagCtrl = CReflect<double>( ctrlC.tridiagCtrl );
    ctrl.tridiagEigCtrl = CReflect( ctrlC.tridiagEigCtrl );
    ctrl.sdcCtrl = CReflect( ctrlC.sdcCtrl );
    ctrl.useSDC = ctrlC.useSDC;
    return ctrl;
}
inline HermitianEigCtrl<Complex<float>>
CReflect( const ElHermitianEigCtrl_c& ctrlC )
{
    HermitianEigCtrl<Complex<float>> ctrl;
    ctrl.tridiagCtrl = CReflect<Complex<float>>( ctrlC.tridiagCtrl );
    ctrl.tridiagEigCtrl = CReflect( ctrlC.tridiagEigCtrl );
    ctrl.sdcCtrl = CReflect( ctrlC.sdcCtrl );
    ctrl.useSDC = ctrlC.useSDC;
    return ctrl;
}
inline HermitianEigCtrl<Complex<double>>
CReflect( const ElHermitianEigCtrl_z& ctrlC )
{
    HermitianEigCtrl<Complex<double>> ctrl;
    ctrl.tridiagCtrl = CReflect<Complex<double>>( ctrlC.tridiagCtrl );
    ctrl.tridiagEigCtrl = CReflect( ctrlC.tridiagEigCtrl );
    ctrl.sdcCtrl = CReflect( ctrlC.sdcCtrl );
    ctrl.useSDC = ctrlC.useSDC;
    return ctrl;
}

/* QDWHCtrl */
inline ElQDWHCtrl CReflect( const QDWHCtrl& ctrl )
{
    ElQDWHCtrl ctrlC;
    ctrlC.colPiv = ctrl.colPiv;
    ctrlC.maxIts = ctrl.maxIts;
    return ctrlC;
}

inline QDWHCtrl CReflect( const ElQDWHCtrl& ctrlC )
{
    QDWHCtrl ctrl;
    ctrl.colPiv = ctrlC.colPiv;
    ctrl.maxIts = ctrlC.maxIts;
    return ctrl;
}

/* PolarCtrl */
inline ElPolarCtrl CReflect( const PolarCtrl& ctrl )
{
    ElPolarCtrl ctrlC;
    ctrlC.qdwh = ctrl.qdwh;
    ctrlC.qdwhCtrl = CReflect(ctrl.qdwhCtrl);
    return ctrlC;
}

inline PolarCtrl CReflect( const ElPolarCtrl& ctrlC )
{
    PolarCtrl ctrl;
    ctrl.qdwh = ctrlC.qdwh;
    ctrl.qdwhCtrl = CReflect(ctrlC.qdwhCtrl);
    return ctrl;
}

/* QDWHInfo */
inline ElQDWHInfo CReflect( const QDWHInfo& info )
{
    ElQDWHInfo infoC;
    infoC.numIts = info.numIts;
    infoC.numQRIts = info.numQRIts;
    infoC.numCholIts = info.numCholIts;
    return infoC;
}

inline QDWHInfo CReflect( const ElQDWHInfo& infoC )
{
    QDWHInfo info;
    info.numIts = infoC.numIts;
    info.numQRIts = infoC.numQRIts;
    info.numCholIts = infoC.numCholIts;
    return info;
}

/* PolarInfo */
inline ElPolarInfo CReflect( const PolarInfo& info )
{
    ElPolarInfo infoC;
    infoC.qdwhInfo = CReflect(info.qdwhInfo);
    return infoC;
}

inline PolarInfo CReflect( const ElPolarInfo& infoC )
{
    PolarInfo info;
    info.qdwhInfo = CReflect(infoC.qdwhInfo);
    return info;
}

/* Bidiag SVD */
inline ElSVDApproach CReflect( SVDApproach approach )
{ return static_cast<ElSVDApproach>( approach ); }

inline SVDApproach CReflect( ElSVDApproach approach )
{ return static_cast<SVDApproach>( approach ); }

inline ElSingularValueToleranceType
CReflect( SingularValueToleranceType tolType )
{ return static_cast<ElSingularValueToleranceType>( tolType ); }

inline SingularValueToleranceType
CReflect( ElSingularValueToleranceType tolType )
{ return static_cast<SingularValueToleranceType>( tolType ); }

/* bidiag_svd::QRCtrl */
inline bidiag_svd::QRCtrl CReflect( const ElBidiagSVDQRCtrl& ctrlC )
{
    bidiag_svd::QRCtrl ctrl;
    ctrl.maxIterPerVal = ctrlC.maxIterPerVal;
    ctrl.demandConverged = ctrlC.demandConverged;
    ctrl.looseMinSingValEst = ctrlC.looseMinSingValEst;
    ctrl.useFLAME = ctrlC.useFLAME;
    ctrl.useLAPACK = ctrlC.useLAPACK;
    ctrl.broadcast = ctrlC.broadcast;
    return ctrl;
}

inline ElBidiagSVDQRCtrl CReflect( const bidiag_svd::QRCtrl& ctrl )
{
    ElBidiagSVDQRCtrl ctrlC;
    ctrlC.maxIterPerVal = ctrl.maxIterPerVal;
    ctrlC.demandConverged = ctrl.demandConverged;
    ctrlC.looseMinSingValEst = ctrl.looseMinSingValEst;
    ctrlC.useFLAME = ctrl.useFLAME;
    ctrlC.useLAPACK = ctrl.useLAPACK;
    ctrlC.broadcast = ctrl.broadcast;
    return ctrlC;
}

/* bidiag_svd::DCCtrl */
inline bidiag_svd::DCCtrl<float> CReflect( const ElBidiagSVDDCCtrl_s& ctrlC )
{
    bidiag_svd::DCCtrl<float> ctrl;
    ctrl.secularCtrl = CReflect( ctrlC.secularCtrl );
    ctrl.deflationFudge = ctrlC.deflationFudge;
    ctrl.cutoff = ctrlC.cutoff;
    ctrl.exploitStructure = ctrlC.exploitStructure;
    return ctrl;
}

inline bidiag_svd::DCCtrl<double> CReflect( const ElBidiagSVDDCCtrl_d& ctrlC )
{
    bidiag_svd::DCCtrl<double> ctrl;
    ctrl.secularCtrl = CReflect( ctrlC.secularCtrl );
    ctrl.deflationFudge = ctrlC.deflationFudge;
    ctrl.cutoff = ctrlC.cutoff;
    ctrl.exploitStructure = ctrlC.exploitStructure;
    return ctrl;
}

inline ElBidiagSVDDCCtrl_s CReflect( const bidiag_svd::DCCtrl<float>& ctrl )
{
    ElBidiagSVDDCCtrl_s ctrlC;
    ctrlC.secularCtrl = CReflect( ctrl.secularCtrl );
    ctrlC.deflationFudge = ctrl.deflationFudge;
    ctrlC.cutoff = ctrl.cutoff;
    ctrlC.exploitStructure = ctrl.exploitStructure;
    return ctrlC;
}
inline ElBidiagSVDDCCtrl_d CReflect( const bidiag_svd::DCCtrl<double>& ctrl )
{
    ElBidiagSVDDCCtrl_d ctrlC;
    ctrlC.secularCtrl = CReflect( ctrl.secularCtrl );
    ctrlC.deflationFudge = ctrl.deflationFudge;
    ctrlC.cutoff = ctrl.cutoff;
    ctrlC.exploitStructure = ctrl.exploitStructure;
    return ctrlC;
}

inline BidiagSVDCtrl<float> CReflect( const ElBidiagSVDCtrl_s& ctrlC )
{
    BidiagSVDCtrl<float> ctrl;
    ctrl.wantU = ctrlC.wantU;
    ctrl.wantV = ctrlC.wantV;
    ctrl.accumulateU = ctrlC.accumulateU;
    ctrl.accumulateV = ctrlC.accumulateV;
    ctrl.approach = CReflect(ctrlC.approach);
    ctrl.tolType = CReflect(ctrlC.tolType);
    ctrl.tol = ctrlC.tol;
    ctrl.progress = ctrlC.progress;
    ctrl.useQR = ctrlC.useQR;
    ctrl.qrCtrl = CReflect(ctrlC.qrCtrl);
    ctrl.dcCtrl = CReflect(ctrlC.dcCtrl);
    return ctrl;
}

inline BidiagSVDCtrl<double> CReflect( const ElBidiagSVDCtrl_d& ctrlC )
{
    BidiagSVDCtrl<double> ctrl;
    ctrl.wantU = ctrlC.wantU;
    ctrl.wantV = ctrlC.wantV;
    ctrl.accumulateU = ctrlC.accumulateU;
    ctrl.accumulateV = ctrlC.accumulateV;
    ctrl.approach = CReflect(ctrlC.approach);
    ctrl.tolType = CReflect(ctrlC.tolType);
    ctrl.tol = ctrlC.tol;
    ctrl.progress = ctrlC.progress;
    ctrl.useQR = ctrlC.useQR;
    ctrl.qrCtrl = CReflect(ctrlC.qrCtrl);
    ctrl.dcCtrl = CReflect(ctrlC.dcCtrl);
    return ctrl;
}

inline ElBidiagSVDCtrl_s CReflect( const BidiagSVDCtrl<float>& ctrl )
{
    ElBidiagSVDCtrl_s ctrlC;
    ctrlC.wantU = ctrl.wantU;
    ctrlC.wantV = ctrl.wantV;
    ctrlC.accumulateU = ctrl.accumulateU;
    ctrlC.accumulateV = ctrl.accumulateV;
    ctrlC.tolType = CReflect(ctrl.tolType);
    ctrlC.tol = ctrl.tol;
    ctrlC.progress = ctrl.progress;
    ctrlC.useQR = ctrl.useQR;
    ctrlC.qrCtrl = CReflect(ctrl.qrCtrl);
    ctrlC.dcCtrl = CReflect(ctrl.dcCtrl);
    return ctrlC;
}

inline ElBidiagSVDCtrl_d CReflect( const BidiagSVDCtrl<double>& ctrl )
{
    ElBidiagSVDCtrl_d ctrlC;
    ctrlC.wantU = ctrl.wantU;
    ctrlC.wantV = ctrl.wantV;
    ctrlC.accumulateU = ctrl.accumulateU;
    ctrlC.accumulateV = ctrl.accumulateV;
    ctrlC.tolType = CReflect(ctrl.tolType);
    ctrlC.tol = ctrl.tol;
    ctrlC.progress = ctrl.progress;
    ctrlC.useQR = ctrl.useQR;
    ctrlC.qrCtrl = CReflect(ctrl.qrCtrl);
    ctrlC.dcCtrl = CReflect(ctrl.dcCtrl);
    return ctrlC;
}

/* SVDCtrl */
inline SVDCtrl<float> CReflect( const ElSVDCtrl_s& ctrlC )
{ SVDCtrl<float> ctrl;
    ctrl.overwrite = ctrlC.overwrite;
    ctrl.time = ctrlC.time;
    ctrl.useLAPACK = ctrlC.useLAPACK;
    ctrl.useScaLAPACK = ctrlC.useScaLAPACK;
    ctrl.valChanRatio = ctrlC.valChanRatio;
    ctrl.fullChanRatio = ctrlC.fullChanRatio;
    ctrl.bidiagSVDCtrl = CReflect(ctrlC.bidiagSVDCtrl);
    return ctrl;
}

inline SVDCtrl<double> CReflect( const ElSVDCtrl_d& ctrlC )
{
    SVDCtrl<double> ctrl;
    ctrl.overwrite = ctrlC.overwrite;
    ctrl.time = ctrlC.time;
    ctrl.useLAPACK = ctrlC.useLAPACK;
    ctrl.useScaLAPACK = ctrlC.useScaLAPACK;
    ctrl.valChanRatio = ctrlC.valChanRatio;
    ctrl.fullChanRatio = ctrlC.fullChanRatio;
    ctrl.bidiagSVDCtrl = CReflect(ctrlC.bidiagSVDCtrl);
    return ctrl;
}

inline ElSVDCtrl_s CReflect( const SVDCtrl<float>& ctrl )
{
    ElSVDCtrl_s ctrlC;
    ctrlC.overwrite = ctrl.overwrite;
    ctrlC.time = ctrl.time;
    ctrlC.useLAPACK = ctrl.useLAPACK;
    ctrlC.useScaLAPACK = ctrl.useScaLAPACK;
    ctrlC.valChanRatio = ctrl.valChanRatio;
    ctrlC.fullChanRatio = ctrl.fullChanRatio;
    ctrlC.bidiagSVDCtrl = CReflect(ctrl.bidiagSVDCtrl);
    return ctrlC;
}

inline ElSVDCtrl_d CReflect( const SVDCtrl<double>& ctrl )
{
    ElSVDCtrl_d ctrlC;
    ctrlC.overwrite = ctrl.overwrite;
    ctrlC.time = ctrl.time;
    ctrlC.useLAPACK = ctrl.useLAPACK;
    ctrlC.useScaLAPACK = ctrl.useScaLAPACK;
    ctrlC.valChanRatio = ctrl.valChanRatio;
    ctrlC.fullChanRatio = ctrl.fullChanRatio;
    ctrlC.bidiagSVDCtrl = CReflect(ctrl.bidiagSVDCtrl);
    return ctrlC;
}

/* HessenbergSchurCtrl */
inline ElHessenbergSchurAlg CReflect( const HessenbergSchurAlg& alg )
{ return static_cast<ElHessenbergSchurAlg>(alg); }
inline HessenbergSchurAlg CReflect( const ElHessenbergSchurAlg& alg )
{ return static_cast<HessenbergSchurAlg>(alg); }

inline ElHessenbergSchurCtrl CReflect( const HessenbergSchurCtrl& ctrl )
{
    ElHessenbergSchurCtrl ctrlC;
    ctrlC.winBeg = ctrl.winBeg;
    ctrlC.winEnd = ctrl.winEnd;
    ctrlC.fullTriangle = ctrl.fullTriangle;
    ctrlC.wantSchurVecs = ctrl.wantSchurVecs;
    ctrlC.accumulateSchurVecs = ctrl.accumulateSchurVecs;
    ctrlC.demandConverged = ctrl.demandConverged;

    ctrlC.alg = CReflect(ctrl.alg);
    ctrlC.recursiveAED = ctrl.recursiveAED;
    ctrlC.accumulateReflections = ctrl.accumulateReflections;
    ctrlC.sortShifts = ctrl.sortShifts;

    ctrlC.progress = ctrl.progress;

    ctrlC.minMultiBulgeSize = ctrl.minMultiBulgeSize;
    ctrlC.minDistMultiBulgeSize = ctrl.minDistMultiBulgeSize;
    auto numShiftsRes = ctrl.numShifts.target<ElInt(*)(ElInt,ElInt)>();
    if( numShiftsRes )
        ctrlC.numShifts = *numShiftsRes;
    else
        RuntimeError("Could not convert numShifts to C function pointer");

    auto deflationSizeRes =
      ctrl.deflationSize.target<ElInt(*)(ElInt,ElInt,ElInt)>();
    if( deflationSizeRes )
        ctrlC.deflationSize = *deflationSizeRes;
    else
        RuntimeError("Could not convert deflationSize to C function pointer");

    auto sufficientDeflationRes =
      ctrl.sufficientDeflation.target<ElInt(*)(ElInt)>();
    if( sufficientDeflationRes )
        ctrlC.sufficientDeflation = *sufficientDeflationRes;
    else
        RuntimeError
        ("Could not convert sufficientDeflation to C function pointer");

    ctrlC.scalapack = ctrl.scalapack;
    ctrlC.blockHeight = ctrl.blockHeight;
    auto numBulgesPerBlockRes =
      ctrl.numBulgesPerBlock.target<ElInt(*)(ElInt)>();
    if( numBulgesPerBlockRes )
        ctrlC.numBulgesPerBlock = *numBulgesPerBlockRes;
    else
        RuntimeError
        ("Could not convert numBulgesPerBlock to C function pointer");

    return ctrlC;
}

inline HessenbergSchurCtrl CReflect( const ElHessenbergSchurCtrl& ctrlC )
{
    HessenbergSchurCtrl ctrl;
    ctrl.winBeg = ctrlC.winBeg;
    ctrl.winEnd = ctrlC.winEnd;
    ctrl.fullTriangle = ctrlC.fullTriangle;
    ctrl.wantSchurVecs = ctrlC.wantSchurVecs;
    ctrl.accumulateSchurVecs = ctrlC.accumulateSchurVecs;
    ctrl.demandConverged = ctrlC.demandConverged;

    ctrl.alg = CReflect(ctrlC.alg);
    ctrl.recursiveAED = ctrlC.recursiveAED;
    ctrl.accumulateReflections = ctrlC.accumulateReflections;
    ctrl.sortShifts = ctrlC.sortShifts;

    ctrl.progress = ctrlC.progress;

    ctrl.minMultiBulgeSize = ctrlC.minMultiBulgeSize;
    ctrl.minDistMultiBulgeSize = ctrlC.minDistMultiBulgeSize;
    ctrl.numShifts = ctrlC.numShifts;
    ctrl.deflationSize = ctrlC.deflationSize;
    ctrl.sufficientDeflation = ctrlC.sufficientDeflation;

    ctrl.scalapack = ctrlC.scalapack;
    ctrl.blockHeight = ctrlC.blockHeight;
    ctrl.numBulgesPerBlock = ctrlC.numBulgesPerBlock;

    return ctrl;
}

/* SDCCtrl */
inline ElSDCCtrl_s CReflect( const SDCCtrl<float>& ctrl )
{
    ElSDCCtrl_s ctrlC;
    ctrlC.cutoff = ctrl.cutoff;
    ctrlC.maxInnerIts = ctrl.maxInnerIts;
    ctrlC.maxOuterIts = ctrl.maxOuterIts;
    ctrlC.tol = ctrl.tol;
    ctrlC.spreadFactor = ctrl.spreadFactor;
    ctrlC.random = ctrl.random;
    ctrlC.progress = ctrl.progress;
    return ctrlC;
}
inline ElSDCCtrl_d CReflect( const SDCCtrl<double>& ctrl )
{
    ElSDCCtrl_d ctrlC;
    ctrlC.cutoff = ctrl.cutoff;
    ctrlC.maxInnerIts = ctrl.maxInnerIts;
    ctrlC.maxOuterIts = ctrl.maxOuterIts;
    ctrlC.tol = ctrl.tol;
    ctrlC.spreadFactor = ctrl.spreadFactor;
    ctrlC.random = ctrl.random;
    ctrlC.progress = ctrl.progress;
    return ctrlC;
}

inline SDCCtrl<float> CReflect( const ElSDCCtrl_s& ctrlC )
{
    SDCCtrl<float> ctrl;
    ctrl.cutoff = ctrlC.cutoff;
    ctrl.maxInnerIts = ctrlC.maxInnerIts;
    ctrl.maxOuterIts = ctrlC.maxOuterIts;
    ctrl.tol = ctrlC.tol;
    ctrl.spreadFactor = ctrlC.spreadFactor;
    ctrl.random = ctrlC.random;
    ctrl.progress = ctrlC.progress;
    return ctrl;
}
inline SDCCtrl<double> CReflect( const ElSDCCtrl_d& ctrlC )
{
    SDCCtrl<double> ctrl;
    ctrl.cutoff = ctrlC.cutoff;
    ctrl.maxInnerIts = ctrlC.maxInnerIts;
    ctrl.maxOuterIts = ctrlC.maxOuterIts;
    ctrl.tol = ctrlC.tol;
    ctrl.spreadFactor = ctrlC.spreadFactor;
    ctrl.random = ctrlC.random;
    ctrl.progress = ctrlC.progress;
    return ctrl;
}

/* SchurCtrl */
inline ElSchurCtrl_s CReflect( const SchurCtrl<float>& ctrl )
{
    ElSchurCtrl_s ctrlC;
    ctrlC.useSDC = ctrl.useSDC;
    ctrlC.hessSchurCtrl = CReflect( ctrl.hessSchurCtrl );
    ctrlC.sdcCtrl = CReflect( ctrl.sdcCtrl );
    ctrlC.time = ctrl.time;
    return ctrlC;
}
inline ElSchurCtrl_d CReflect( const SchurCtrl<double>& ctrl )
{
    ElSchurCtrl_d ctrlC;
    ctrlC.useSDC = ctrl.useSDC;
    ctrlC.hessSchurCtrl = CReflect( ctrl.hessSchurCtrl );
    ctrlC.sdcCtrl = CReflect( ctrl.sdcCtrl );
    ctrlC.time = ctrl.time;
    return ctrlC;
}

inline SchurCtrl<float> CReflect( const ElSchurCtrl_s& ctrlC )
{
    SchurCtrl<float> ctrl;
    ctrl.useSDC = ctrlC.useSDC;
    ctrl.hessSchurCtrl = CReflect( ctrlC.hessSchurCtrl );
    ctrl.sdcCtrl = CReflect( ctrlC.sdcCtrl );
    ctrl.time = ctrlC.time;
    return ctrl;
}
inline SchurCtrl<double> CReflect( const ElSchurCtrl_d& ctrlC )
{
    SchurCtrl<double> ctrl;
    ctrl.useSDC = ctrlC.useSDC;
    ctrl.hessSchurCtrl = CReflect( ctrlC.hessSchurCtrl );
    ctrl.sdcCtrl = CReflect( ctrlC.sdcCtrl );
    ctrl.time = ctrlC.time;
    return ctrl;
}

inline ElPseudospecNorm CReflect( PseudospecNorm psNorm )
{ return static_cast<ElPseudospecNorm>(psNorm); }
inline PseudospecNorm CReflect( ElPseudospecNorm psNorm )
{ return static_cast<PseudospecNorm>(psNorm); }

inline ElSnapshotCtrl CReflect( const SnapshotCtrl& ctrl )
{
    ElSnapshotCtrl ctrlC;
    ctrlC.realSize = ctrl.realSize;
    ctrlC.imagSize = ctrl.imagSize;
    ctrlC.imgSaveFreq = ctrl.imgSaveFreq;
    ctrlC.numSaveFreq = ctrl.numSaveFreq;
    ctrlC.imgDispFreq = ctrl.imgDispFreq;
    ctrlC.imgSaveCount = ctrl.imgSaveCount;
    ctrlC.numSaveCount = ctrl.numSaveCount;
    ctrlC.imgDispCount = ctrl.imgDispCount;
    ctrlC.imgBase = CReflect(ctrl.imgBase);
    ctrlC.numBase = CReflect(ctrl.numBase);
    ctrlC.imgFormat = CReflect(ctrl.imgFormat);
    ctrlC.numFormat = CReflect(ctrl.numFormat);
    ctrlC.itCounts = ctrl.itCounts;
    return ctrlC;
}
inline SnapshotCtrl CReflect( const ElSnapshotCtrl& ctrlC )
{
    SnapshotCtrl ctrl;
    ctrl.realSize = ctrlC.realSize;
    ctrl.imagSize = ctrlC.imagSize;
    ctrl.imgSaveFreq = ctrlC.imgSaveFreq;
    ctrl.numSaveFreq = ctrlC.numSaveFreq;
    ctrl.imgDispFreq = ctrlC.imgDispFreq;
    ctrl.imgSaveCount = ctrlC.imgSaveCount;
    ctrl.numSaveCount = ctrlC.numSaveCount;
    ctrl.imgDispCount = ctrlC.imgDispCount;
    ctrl.imgBase = CReflect(ctrlC.imgBase);
    ctrl.numBase = CReflect(ctrlC.numBase);
    ctrl.imgFormat = CReflect(ctrlC.imgFormat);
    ctrl.numFormat = CReflect(ctrlC.numFormat);
    ctrl.itCounts = ctrlC.itCounts;
    return ctrl;
}

inline ElPseudospecCtrl_s CReflect( const PseudospecCtrl<float>& ctrl )
{
    ElPseudospecCtrl_s ctrlC;
    ctrlC.norm = CReflect(ctrl.norm);
    ctrlC.blockWidth = ctrl.norm;
    ctrlC.schur = ctrl.schur;
    ctrlC.forceComplexSchur = ctrl.forceComplexSchur;
    ctrlC.forceComplexPs = ctrl.forceComplexPs;
    ctrlC.schurCtrl = CReflect(ctrl.schurCtrl);
    ctrlC.maxIts = ctrl.maxIts;
    ctrlC.tol = ctrl.tol;
    ctrlC.deflate = ctrl.deflate;
    ctrlC.arnoldi = ctrl.arnoldi;
    ctrlC.basisSize = ctrl.basisSize;
    ctrlC.reorthog = ctrl.reorthog;
    ctrlC.progress = ctrl.progress;
    ctrlC.snapCtrl = CReflect(ctrl.snapCtrl);
    return ctrlC;
}
inline ElPseudospecCtrl_d CReflect( const PseudospecCtrl<double>& ctrl )
{
    ElPseudospecCtrl_d ctrlC;
    ctrlC.norm = CReflect(ctrl.norm);
    ctrlC.blockWidth = ctrl.norm;
    ctrlC.schur = ctrl.schur;
    ctrlC.forceComplexSchur = ctrl.forceComplexSchur;
    ctrlC.forceComplexPs = ctrl.forceComplexPs;
    ctrlC.schurCtrl = CReflect(ctrl.schurCtrl);
    ctrlC.maxIts = ctrl.maxIts;
    ctrlC.tol = ctrl.tol;
    ctrlC.deflate = ctrl.deflate;
    ctrlC.arnoldi = ctrl.arnoldi;
    ctrlC.basisSize = ctrl.basisSize;
    ctrlC.reorthog = ctrl.reorthog;
    ctrlC.progress = ctrl.progress;
    ctrlC.snapCtrl = CReflect(ctrl.snapCtrl);
    return ctrlC;
}

inline PseudospecCtrl<float> CReflect( const ElPseudospecCtrl_s& ctrlC )
{
    PseudospecCtrl<float> ctrl;
    ctrl.norm = CReflect(ctrlC.norm);
    ctrl.blockWidth = ctrlC.norm;
    ctrl.schur = ctrlC.schur;
    ctrl.forceComplexSchur = ctrlC.forceComplexSchur;
    ctrl.forceComplexPs = ctrlC.forceComplexPs;
    ctrl.schurCtrl = CReflect(ctrlC.schurCtrl);
    ctrl.maxIts = ctrlC.maxIts;
    ctrl.tol = ctrlC.tol;
    ctrl.deflate = ctrlC.deflate;
    ctrl.arnoldi = ctrlC.arnoldi;
    ctrl.basisSize = ctrlC.basisSize;
    ctrl.reorthog = ctrlC.reorthog;
    ctrl.progress = ctrlC.progress;
    ctrl.snapCtrl = CReflect(ctrlC.snapCtrl);
    return ctrl;
}
inline PseudospecCtrl<double> CReflect( const ElPseudospecCtrl_d& ctrlC )
{
    PseudospecCtrl<double> ctrl;
    ctrl.norm = CReflect(ctrlC.norm);
    ctrl.blockWidth = ctrlC.norm;
    ctrl.schur = ctrlC.schur;
    ctrl.forceComplexSchur = ctrlC.forceComplexSchur;
    ctrl.forceComplexPs = ctrlC.forceComplexPs;
    ctrl.schurCtrl = CReflect(ctrlC.schurCtrl);
    ctrl.maxIts = ctrlC.maxIts;
    ctrl.tol = ctrlC.tol;
    ctrl.deflate = ctrlC.deflate;
    ctrl.arnoldi = ctrlC.arnoldi;
    ctrl.basisSize = ctrlC.basisSize;
    ctrl.reorthog = ctrlC.reorthog;
    ctrl.progress = ctrlC.progress;
    ctrl.snapCtrl = CReflect(ctrlC.snapCtrl);
    return ctrl;
}

inline ElSpectralBox_s CReflect( const SpectralBox<float>& box )
{
    ElSpectralBox_s boxC;
    boxC.center = CReflect(box.center);
    boxC.realWidth = box.realWidth;
    boxC.imagWidth = box.imagWidth;
    return boxC;
}

inline ElSpectralBox_d CReflect( const SpectralBox<double>& box )
{
    ElSpectralBox_d boxC;
    boxC.center = CReflect(box.center);
    boxC.realWidth = box.realWidth;
    boxC.imagWidth = box.imagWidth;
    return boxC;
}

inline SpectralBox<float> CReflect( const ElSpectralBox_s& boxC )
{
    SpectralBox<float> box;
    box.center = CReflect(boxC.center);
    box.realWidth = CReflect(boxC.realWidth);
    box.imagWidth = CReflect(boxC.imagWidth);
    return box;
}

inline SpectralBox<double> CReflect( const ElSpectralBox_d& boxC )
{
    SpectralBox<double> box;
    box.center = CReflect(boxC.center);
    box.realWidth = CReflect(boxC.realWidth);
    box.imagWidth = CReflect(boxC.imagWidth);
    return box;
}

} // namespace El

#endif // ifndef EL_LAPACK_SPECTRAL_CREFLECT_C_HPP
