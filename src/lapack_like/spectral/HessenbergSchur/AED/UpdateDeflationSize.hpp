/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_HESS_SCHUR_AED_UPDATE_DEFLATION_SIZE_HPP
#define EL_HESS_SCHUR_AED_UPDATE_DEFLATION_SIZE_HPP

namespace El {
namespace hess_schur {
namespace aed {

// Intelligently choose a deflation window size
// --------------------------------------------
// Cf. LAPACK's DLAQR0 for the high-level approach

template<typename Field>
void UpdateDeflationSize
( Int& deflationSize,
  Int& decreaseLevel,
  Int deflationSizeRec,
  Int numIterSinceDeflation,
  Int numStaleIterBeforeExceptional,
  const Matrix<Field>& hSubIter )
{
    const Int iterWinSize = hSubIter.Height() + 1;
    if( numIterSinceDeflation < numStaleIterBeforeExceptional )
    {
        // Use the recommendation if possible
        deflationSize = Min( iterWinSize, deflationSizeRec );
    }
    else
    {
        // Double the size if possible
        deflationSize = Min( iterWinSize, 2*deflationSize );
    }
    if( deflationSize > iterWinSize-2 )
    {
        // Go ahead and increase by at most one to use the full window
        deflationSize = iterWinSize;
    }
    else
    {
        const Int deflationBeg = iterWinSize - deflationSize;
        if( Abs(hSubIter(deflationBeg-1)) > Abs(hSubIter(deflationBeg-2)) )
            ++deflationSize;
    }
    if( numIterSinceDeflation < numStaleIterBeforeExceptional )
    {
        decreaseLevel = -1;
    }
    else if( decreaseLevel >= 0 || deflationSize == iterWinSize )
    {
        ++decreaseLevel;
        if( deflationSize-decreaseLevel < 2 )
            decreaseLevel = 0;
        deflationSize -= decreaseLevel;
    }
}

} // namespace aed
} // namespace hess_schur
} // namespace El

#endif // ifndef EL_HESS_SCHUR_AED_UPDATE_DEFLATION_SIZE_HPP
