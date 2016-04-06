/*
   Copyright (c) 2009-2016, Jack Poulson and Tim Moon
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace El {
namespace safemstrsm {

/* Determine machine dependent parameters to control overflow
 *   Note: LAPACK uses more complicated parameters to handle 
 *   issues that can happen on Cray machines.
 */
template<typename Real>
inline pair<Real,Real>
OverflowParameters()
{
    const Real underflow = lapack::MachineSafeMin<Real>();
    const Real overflow = lapack::MachineOverflowThreshold<Real>();
    const Real ulp  = lapack::MachinePrecision<Real>();
    const Real smallNum = Max( underflow/ulp, 1/(overflow*ulp) );
    const Real bigNum = 1/smallNum;
    return pair<Real,Real>(smallNum,bigNum);
}

}
}
