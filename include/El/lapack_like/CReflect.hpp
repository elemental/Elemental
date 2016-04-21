/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_LAPACK_CREFLECT_C_HPP
#define EL_LAPACK_CREFLECT_C_HPP

namespace El {

// LAPACK-like
// -----------

inline ElSortType CReflect( SortType type )
{ return static_cast<ElSortType>(type); }

inline SortType CReflect( ElSortType type )
{ return static_cast<SortType>(type); }

// Permutations
// ^^^^^^^^^^^^

inline ElPermutationMeta CReflect( const PermutationMeta& meta )
{
    ElPermutationMeta metaC;    

    metaC.align = meta.align;
    metaC.comm = meta.comm.comm;

    const Int commSize = mpi::Size( meta.comm );
    metaC.sendCounts = new int[commSize];
    metaC.sendDispls = new int[commSize];
    metaC.recvCounts = new int[commSize];
    metaC.recvDispls = new int[commSize];
    MemCopy( metaC.sendCounts, meta.sendCounts.data(), commSize );
    MemCopy( metaC.sendDispls, meta.sendDispls.data(), commSize ); 
    MemCopy( metaC.recvCounts, meta.recvCounts.data(), commSize );
    MemCopy( metaC.recvDispls, meta.recvDispls.data(), commSize );

    metaC.numSendIdx = meta.sendIdx.size();
    metaC.numRecvIdx = meta.recvIdx.size();
    metaC.sendIdx   = new int[metaC.numSendIdx];
    metaC.sendRanks = new int[metaC.numSendIdx];
    metaC.recvIdx   = new int[metaC.numRecvIdx];
    metaC.recvRanks = new int[metaC.numRecvIdx];
    MemCopy( metaC.sendIdx,   meta.sendIdx.data(),   metaC.numSendIdx );
    MemCopy( metaC.sendRanks, meta.sendRanks.data(), metaC.numSendIdx );
    MemCopy( metaC.recvIdx,   meta.recvIdx.data(),   metaC.numRecvIdx );
    MemCopy( metaC.recvRanks, meta.recvRanks.data(), metaC.numRecvIdx );

    return metaC;
}

inline PermutationMeta CReflect( const ElPermutationMeta& metaC )
{
    PermutationMeta meta;

    meta.align = metaC.align;
    meta.comm = metaC.comm;

    int commSize;
    MPI_Comm_size( metaC.comm, &commSize );
    meta.sendCounts = 
        vector<int>( metaC.sendCounts, metaC.sendCounts+commSize );
    meta.sendDispls = 
        vector<int>( metaC.sendDispls, metaC.sendDispls+commSize );
    meta.recvCounts =
        vector<int>( metaC.recvCounts, metaC.recvCounts+commSize );
    meta.recvDispls =
        vector<int>( metaC.recvDispls, metaC.recvDispls+commSize );

    meta.sendIdx = 
        vector<int>( metaC.sendIdx, metaC.sendIdx+metaC.numSendIdx );
    meta.sendRanks =
        vector<int>( metaC.sendRanks, metaC.sendRanks+metaC.numSendIdx );
    meta.recvIdx =
        vector<int>( metaC.recvIdx, metaC.recvIdx+metaC.numRecvIdx );
    meta.recvRanks =
        vector<int>( metaC.recvRanks, metaC.recvRanks+metaC.numRecvIdx );

    return meta;
}

inline Permutation* CReflect( ElPermutation p )
{ return EL_RC(Permutation*,p); }
inline ElPermutation CReflect( Permutation* p )
{ return (ElPermutation)EL_RC(struct ElPermutationDummy*,p); }

inline DistPermutation* CReflect( ElDistPermutation p )
{ return EL_RC(DistPermutation*,p); }
inline ElDistPermutation CReflect( DistPermutation* p )
{ return (ElDistPermutation)EL_RC(struct ElDistPermutationDummy*,p); }

inline const Permutation* CReflect( ElConstPermutation p )
{ return EL_RC(const Permutation*,p); }
inline ElConstPermutation CReflect( const Permutation* p )
{ return (ElConstPermutation)EL_RC(const struct ElPermutationDummy*,p); }

inline const DistPermutation* CReflect( ElConstDistPermutation p )
{ return EL_RC(const DistPermutation*,p); }
inline ElConstDistPermutation CReflect( const DistPermutation* p )
{ return (ElConstDistPermutation)
         EL_RC(const struct ElDistPermutationDummy*,p); }

// Condensed form
// ^^^^^^^^^^^^^^
inline ElHermitianTridiagApproach 
CReflect( HermitianTridiagApproach approach )
{ return static_cast<ElHermitianTridiagApproach>( approach ); }

inline HermitianTridiagApproach 
CReflect( ElHermitianTridiagApproach approach )
{ return static_cast<HermitianTridiagApproach>( approach ); }

template<typename F>
inline ElHermitianTridiagCtrl
CReflect( const HermitianTridiagCtrl<F>& ctrl )
{ 
    ElHermitianTridiagCtrl ctrlC;
    ctrlC.approach = CReflect(ctrl.approach);
    ctrlC.order = CReflect(ctrl.order);
    ctrlC.symvCtrl = CReflect(ctrl.symvCtrl);
    return ctrlC;
}

template<typename F>
inline HermitianTridiagCtrl<F>
CReflect( const ElHermitianTridiagCtrl& ctrlC )
{ 
    HermitianTridiagCtrl<F> ctrl;
    ctrl.approach = CReflect(ctrlC.approach);
    ctrl.order = CReflect(ctrlC.order);
    ctrl.symvCtrl = CReflect<F>(ctrlC.symvCtrl);
    return ctrl;
}

// Decompositions
// ^^^^^^^^^^^^^^

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

/* HermitianEigCtrl */
inline ElHermitianEigCtrl_s CReflect( const HermitianEigCtrl<float>& ctrl )
{
    ElHermitianEigCtrl_s ctrlC;
    ctrlC.tridiagCtrl = CReflect( ctrl.tridiagCtrl );
    ctrlC.sdcCtrl = CReflect( ctrl.sdcCtrl );
    ctrlC.useSDC = ctrl.useSDC;
    return ctrlC;
}

inline ElHermitianEigCtrl_d CReflect( const HermitianEigCtrl<double>& ctrl )
{
    ElHermitianEigCtrl_d ctrlC;
    ctrlC.tridiagCtrl = CReflect( ctrl.tridiagCtrl );
    ctrlC.sdcCtrl = CReflect( ctrl.sdcCtrl );
    ctrlC.useSDC = ctrl.useSDC;
    return ctrlC;
}
inline ElHermitianEigCtrl_c 
CReflect( const HermitianEigCtrl<Complex<float>>& ctrl )
{
    ElHermitianEigCtrl_c ctrlC;
    ctrlC.tridiagCtrl = CReflect( ctrl.tridiagCtrl );
    ctrlC.sdcCtrl = CReflect( ctrl.sdcCtrl );
    ctrlC.useSDC = ctrl.useSDC;
    return ctrlC;
}
inline ElHermitianEigCtrl_z
CReflect( const HermitianEigCtrl<Complex<double>>& ctrl )
{
    ElHermitianEigCtrl_z ctrlC;
    ctrlC.tridiagCtrl = CReflect( ctrl.tridiagCtrl );
    ctrlC.sdcCtrl = CReflect( ctrl.sdcCtrl );
    ctrlC.useSDC = ctrl.useSDC;
    return ctrlC;
}

inline HermitianEigCtrl<float> CReflect( const ElHermitianEigCtrl_s& ctrlC )
{
    HermitianEigCtrl<float> ctrl;
    ctrl.tridiagCtrl = CReflect<float>( ctrlC.tridiagCtrl );
    ctrl.sdcCtrl = CReflect( ctrlC.sdcCtrl );
    ctrl.useSDC = ctrlC.useSDC;
    return ctrl;
}
inline HermitianEigCtrl<double> CReflect( const ElHermitianEigCtrl_d& ctrlC )
{
    HermitianEigCtrl<double> ctrl;
    ctrl.tridiagCtrl = CReflect<double>( ctrlC.tridiagCtrl );
    ctrl.sdcCtrl = CReflect( ctrlC.sdcCtrl );
    ctrl.useSDC = ctrlC.useSDC;
    return ctrl;
}
inline HermitianEigCtrl<Complex<float>> 
CReflect( const ElHermitianEigCtrl_c& ctrlC )
{
    HermitianEigCtrl<Complex<float>> ctrl;
    ctrl.tridiagCtrl = CReflect<Complex<float>>( ctrlC.tridiagCtrl );
    ctrl.sdcCtrl = CReflect( ctrlC.sdcCtrl );
    ctrl.useSDC = ctrlC.useSDC;
    return ctrl;
}
inline HermitianEigCtrl<Complex<double>> 
CReflect( const ElHermitianEigCtrl_z& ctrlC )
{
    HermitianEigCtrl<Complex<double>> ctrl;
    ctrl.tridiagCtrl = CReflect<Complex<double>>( ctrlC.tridiagCtrl );
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

/* SVDCtrl */
inline ElSVDApproach CReflect( SVDApproach approach )
{ return static_cast<ElSVDApproach>( approach ); }

inline SVDApproach CReflect( ElSVDApproach approach )
{ return static_cast<SVDApproach>( approach ); }

inline SVDCtrl<float> CReflect( const ElSVDCtrl_s& ctrlC )
{
    SVDCtrl<float> ctrl;
    ctrl.approach = CReflect(ctrlC.approach);
    ctrl.overwrite = ctrlC.overwrite;
    ctrl.avoidComputingU = ctrlC.avoidComputingU;
    ctrl.avoidComputingV = ctrlC.avoidComputingV;
    ctrl.time = ctrlC.time;
    ctrl.avoidLibflame = ctrlC.avoidLibflame;

    ctrl.seqQR = ctrlC.seqQR;
    ctrl.valChanRatio = ctrlC.valChanRatio;
    ctrl.fullChanRatio = ctrlC.fullChanRatio;
    ctrl.relative = ctrlC.relative;
    ctrl.tol = ctrlC.tol;

    return ctrl;
}

inline SVDCtrl<double> CReflect( const ElSVDCtrl_d& ctrlC )
{
    SVDCtrl<double> ctrl;
    ctrl.approach = CReflect(ctrlC.approach);
    ctrl.overwrite = ctrlC.overwrite;
    ctrl.avoidComputingU = ctrlC.avoidComputingU;
    ctrl.avoidComputingV = ctrlC.avoidComputingV;
    ctrl.time = ctrlC.time;
    ctrl.avoidLibflame = ctrlC.avoidLibflame;

    ctrl.seqQR = ctrlC.seqQR;
    ctrl.valChanRatio = ctrlC.valChanRatio;
    ctrl.fullChanRatio = ctrlC.fullChanRatio;
    ctrl.relative = ctrlC.relative;
    ctrl.tol = ctrlC.tol;

    return ctrl;
}

inline ElSVDCtrl_s CReflect( const SVDCtrl<float>& ctrl )
{
    ElSVDCtrl_s ctrlC;
    ctrlC.approach = CReflect(ctrl.approach);
    ctrlC.overwrite = ctrl.overwrite;
    ctrlC.avoidComputingU = ctrl.avoidComputingU;
    ctrlC.avoidComputingV = ctrl.avoidComputingV;
    ctrlC.time = ctrl.time;
    ctrlC.avoidLibflame = ctrl.avoidLibflame;

    ctrlC.seqQR = ctrl.seqQR;
    ctrlC.valChanRatio = ctrl.valChanRatio;
    ctrlC.fullChanRatio = ctrl.fullChanRatio;
    ctrlC.relative = ctrl.relative;
    ctrlC.tol = ctrl.tol;

    return ctrlC;
}

inline ElSVDCtrl_d CReflect( const SVDCtrl<double>& ctrl )
{
    ElSVDCtrl_d ctrlC;
    ctrlC.approach = CReflect(ctrl.approach);
    ctrlC.overwrite = ctrl.overwrite;
    ctrlC.avoidComputingU = ctrl.avoidComputingU;
    ctrlC.avoidComputingV = ctrl.avoidComputingV;
    ctrlC.time = ctrl.time;
    ctrlC.avoidLibflame = ctrl.avoidLibflame;

    ctrlC.seqQR = ctrl.seqQR;
    ctrlC.valChanRatio = ctrl.valChanRatio;
    ctrlC.fullChanRatio = ctrl.fullChanRatio;
    ctrlC.relative = ctrl.relative;
    ctrlC.tol = ctrl.tol;

    return ctrlC;
}

/* HessQRCtrl */
inline ElHessQRCtrl CReflect( const HessQRCtrl& ctrl )
{
    ElHessQRCtrl ctrlC;
    ctrlC.distAED = ctrl.distAED;
    ctrlC.blockHeight = ctrl.blockHeight;
    ctrlC.blockWidth = ctrl.blockWidth;
    return ctrlC;
}

inline HessQRCtrl CReflect( const ElHessQRCtrl& ctrlC )
{
    HessQRCtrl ctrl;
    ctrl.distAED = ctrlC.distAED;
    ctrl.blockHeight = ctrlC.blockHeight;
    ctrl.blockWidth = ctrlC.blockWidth;
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
    ctrlC.qrCtrl = CReflect( ctrl.qrCtrl );
    ctrlC.sdcCtrl = CReflect( ctrl.sdcCtrl );
    ctrlC.time = ctrl.time;
    return ctrlC;
}
inline ElSchurCtrl_d CReflect( const SchurCtrl<double>& ctrl )
{
    ElSchurCtrl_d ctrlC;
    ctrlC.useSDC = ctrl.useSDC;
    ctrlC.qrCtrl = CReflect( ctrl.qrCtrl );
    ctrlC.sdcCtrl = CReflect( ctrl.sdcCtrl );
    ctrlC.time = ctrl.time;
    return ctrlC;
}

inline SchurCtrl<float> CReflect( const ElSchurCtrl_s& ctrlC )
{
    SchurCtrl<float> ctrl;
    ctrl.useSDC = ctrlC.useSDC;
    ctrl.qrCtrl = CReflect( ctrlC.qrCtrl );
    ctrl.sdcCtrl = CReflect( ctrlC.sdcCtrl );
    ctrl.time = ctrlC.time;
    return ctrl;
}
inline SchurCtrl<double> CReflect( const ElSchurCtrl_d& ctrlC )
{
    SchurCtrl<double> ctrl;
    ctrl.useSDC = ctrlC.useSDC;
    ctrl.qrCtrl = CReflect( ctrlC.qrCtrl );
    ctrl.sdcCtrl = CReflect( ctrlC.sdcCtrl );
    ctrl.time = ctrlC.time;
    return ctrl;
}

// Factorizations
// ^^^^^^^^^^^^^^
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

inline ElLeastSquaresCtrl_s CReflect( const LeastSquaresCtrl<float>& ctrl )
{
    ElLeastSquaresCtrl_s ctrlC;
    ctrlC.scaleTwoNorm = ctrl.scaleTwoNorm;
    ctrlC.basisSize    = ctrl.basisSize;
    ctrlC.alpha        = ctrl.alpha; 
    ctrlC.solveCtrl    = CReflect(ctrl.solveCtrl);
    ctrlC.equilibrate  = ctrl.equilibrate;
    ctrlC.progress     = ctrl.progress;
    ctrlC.time         = ctrl.time;
    return ctrlC;
}

inline ElLeastSquaresCtrl_d CReflect( const LeastSquaresCtrl<double>& ctrl )
{
    ElLeastSquaresCtrl_d ctrlC;
    ctrlC.scaleTwoNorm = ctrl.scaleTwoNorm;
    ctrlC.basisSize    = ctrl.basisSize;
    ctrlC.alpha        = ctrl.alpha; 
    ctrlC.solveCtrl    = CReflect(ctrl.solveCtrl);
    ctrlC.equilibrate  = ctrl.equilibrate;
    ctrlC.progress     = ctrl.progress;
    ctrlC.time         = ctrl.time;
    return ctrlC;
}

inline LeastSquaresCtrl<float> CReflect( const ElLeastSquaresCtrl_s& ctrlC )
{
    LeastSquaresCtrl<float> ctrl;
    ctrl.scaleTwoNorm = ctrlC.scaleTwoNorm;
    ctrl.basisSize    = ctrlC.basisSize;
    ctrl.alpha        = ctrlC.alpha; 
    ctrl.solveCtrl    = CReflect(ctrlC.solveCtrl);
    ctrl.equilibrate  = ctrlC.equilibrate;
    ctrl.progress     = ctrlC.progress;
    ctrl.time         = ctrlC.time;
    return ctrl;
}

inline LeastSquaresCtrl<double> CReflect( const ElLeastSquaresCtrl_d& ctrlC )
{
    LeastSquaresCtrl<double> ctrl;
    ctrl.scaleTwoNorm = ctrlC.scaleTwoNorm;
    ctrl.basisSize    = ctrlC.basisSize;
    ctrl.alpha        = ctrlC.alpha; 
    ctrl.solveCtrl    = CReflect(ctrlC.solveCtrl);
    ctrl.equilibrate  = ctrlC.equilibrate;
    ctrl.progress     = ctrlC.progress;
    ctrl.time         = ctrlC.time;
    return ctrl;
}

// Properties
// ^^^^^^^^^^
inline ElNormType CReflect( NormType type )
{ return static_cast<ElNormType>(type); }
inline NormType CReflect( ElNormType type )
{ return static_cast<NormType>(type); }

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

// Solvers
// ^^^^^^^
inline ElTikhonovAlg CReflect( TikhonovAlg alg )
{ return static_cast<ElTikhonovAlg>(alg); }
inline TikhonovAlg CReflect( ElTikhonovAlg alg )
{ return static_cast<TikhonovAlg>(alg); }

inline ElRidgeAlg CReflect( RidgeAlg alg )
{ return static_cast<ElRidgeAlg>(alg); }
inline RidgeAlg CReflect( ElRidgeAlg alg )
{ return static_cast<RidgeAlg>(alg); }

} // namespace El

#endif // ifndef EL_LAPACK_CREFLECT_C_HPP
