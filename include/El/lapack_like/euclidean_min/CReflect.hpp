/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_LAPACK_EUCLIDEAN_MIN_CREFLECT_C_HPP
#define EL_LAPACK_EUCLIDEAN_MIN_CREFLECT_C_HPP

namespace El {

inline ElLeastSquaresCtrl_s CReflect( const LeastSquaresCtrl<float>& ctrl )
{
    ElLeastSquaresCtrl_s ctrlC;
    ctrlC.scaleTwoNorm = ctrl.scaleTwoNorm;
    ctrlC.basisSize    = ctrl.basisSize;
    ctrlC.alpha        = ctrl.alpha; 
    ctrlC.reg0Tmp      = ctrl.reg0Tmp;
    ctrlC.reg0Perm     = ctrl.reg0Perm;
    ctrlC.reg1Tmp      = ctrl.reg1Tmp;
    ctrlC.reg1Perm     = ctrl.reg1Perm;
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
    ctrlC.reg0Tmp      = ctrl.reg0Tmp;
    ctrlC.reg0Perm     = ctrl.reg0Perm;
    ctrlC.reg1Tmp      = ctrl.reg1Tmp;
    ctrlC.reg1Perm     = ctrl.reg1Perm;
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
    ctrl.reg0Tmp      = ctrlC.reg0Tmp;
    ctrl.reg0Perm     = ctrlC.reg0Perm;
    ctrl.reg1Tmp      = ctrlC.reg1Tmp;
    ctrl.reg1Perm     = ctrlC.reg1Perm;
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
    ctrl.reg0Tmp      = ctrlC.reg0Tmp;
    ctrl.reg0Perm     = ctrlC.reg0Perm;
    ctrl.reg1Tmp      = ctrlC.reg1Tmp;
    ctrl.reg1Perm     = ctrlC.reg1Perm;
    ctrl.solveCtrl    = CReflect(ctrlC.solveCtrl);
    ctrl.equilibrate  = ctrlC.equilibrate;
    ctrl.progress     = ctrlC.progress;
    ctrl.time         = ctrlC.time;
    return ctrl;
}

inline ElTikhonovAlg CReflect( TikhonovAlg alg )
{ return static_cast<ElTikhonovAlg>(alg); }
inline TikhonovAlg CReflect( ElTikhonovAlg alg )
{ return static_cast<TikhonovAlg>(alg); }

inline ElRidgeAlg CReflect( RidgeAlg alg )
{ return static_cast<ElRidgeAlg>(alg); }
inline RidgeAlg CReflect( ElRidgeAlg alg )
{ return static_cast<RidgeAlg>(alg); }

} // namespace El

#endif // ifndef EL_LAPACK_EUCLIDEAN_MIN_CREFLECT_C_HPP
