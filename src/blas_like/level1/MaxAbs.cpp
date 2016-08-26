/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El-lite.hpp>
#include <El/blas_like/level1.hpp>

namespace El {

    template<typename Real,typename>
    Real MaxAbs( const Matrix<Real>& A )
    {
        DEBUG_CSE
        const Int m = A.Height();
        const Int n = A.Width();
        const Real* ABuf = A.LockedBuffer();
        const Int ALDim = A.LDim();

        Real value = 0;
        for( Int j=0; j<n; ++j )
            for( Int i=0; i<m; ++i )
                value = Max(std::abs(value),ABuf[i+j*ALDim]);
        return value;
    }

    template<typename Real,typename>
    Real MaxAbs( const AbstractDistMatrix<Real>& A )
    {
        DEBUG_CSE
        DEBUG_ONLY(
          if( !A.Grid().InGrid() )
              LogicError("Viewing processes are not allowed");
        )
        Real value = 0;
        if( A.Participating() )
        {
            // Store the index/value of the local pivot candidate
            const Int mLocal = A.LocalHeight();
            const Int nLocal = A.LocalWidth();
            const Real* ABuf = A.LockedBuffer();
            const Int ALDim = A.LDim();
            for( Int jLoc=0; jLoc<nLocal; ++jLoc )
                for( Int iLoc=0; iLoc<mLocal; ++iLoc )
                    value = Max(value,ABuf[iLoc+jLoc*ALDim]);

            value = mpi::AllReduce( std::abs(value), mpi::MAX, A.DistComm() );
        }
        mpi::Broadcast( value, A.Root(), A.CrossComm() );
        return value;
    }

}
