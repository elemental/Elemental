/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {

namespace svp {

// See Algorithm 2 from:
//
//   Nicolas Gama, Phong Q. Nguyen, and Oded Regev,
//   "Lattice enumeration using extreme pruning", Eurocrypt 2010.
//
// Note that our algorithm uses 'R' to denote the upper triangular matrix from
// the QR factorization of B and 'upperBounds' to denote the sequence of
// upper-bounds
//
//     upperBounds(0) <= upperBounds(1) <= ... <= upperBounds(n-1)
//
// such that upperBounds(j) bounds the two-norm of the last j+1 entries of
// R v = (diag(d) N) v.
//
// See the documentation in <El/number_theory/lattice/Enumerate.hpp> for more
// details on the complex generalization of GNR enumeration.

namespace gnr_enum {

template<typename F>
Base<F> Helper
( const Matrix<Base<F>>& d,
  const Matrix<F>& N,
  const Matrix<Base<F>>& upperBounds,
        Matrix<F>& v,
  const EnumCtrl<Base<F>>& ctrl )
{
    EL_DEBUG_CSE
    typedef Base<F> Real;
    const Int m = N.Height();
    const Int n = N.Width();
    if( n > m )
        LogicError("Expected height(N) >= width(N)");

    Matrix<F> partialSums;
    Zeros( partialSums, n+1, n );

    // The 'sumIndices' vector is indexed differently from the 'r' from
    // Gama/Nguyen/Regev, and (r_0,r_1,r_2,...,r_n)=(0,1,2,...,n) becomes
    // (-1,0,1,...,n-1).
    Matrix<Int> sumIndices;
    Zeros( sumIndices, n+1, 1 );
    for( Int j=0; j<=n; ++j )
        sumIndices(j) = j-1;

    // Note: We maintain the norms rather than their squares
    Matrix<Real> partialNorms;
    Zeros( partialNorms, n+1, 1 );

    Zeros( v, n, 1 );
    if( n == 0 )
        return Real(0);
    v(0) = F(1);
    Int lastNonzero = 0; // will be set to -1 if all sumIndices are zero

    Matrix<F> centers;
    Zeros( centers, n, 1 );

    vector<SpiralState<F>> spiralStates(n);

    Int k=0;
    F* vBuf = &v(0);
    while( true )
    {
        const F entry = d(k)*(vBuf[k] - centers(k));
        const Real partialNorm = SafeNorm( partialNorms(k+1), entry );
        partialNorms(k) = partialNorm;
        if( partialNorm < upperBounds((n-1)-k) )
        {
            if( k == 0 )
            {
                // Success
                return partialNorm;
            }
            else
            {
                // Move down the tree
                --k;
                sumIndices(k) = Max(sumIndices(k),sumIndices(k+1));

                F* s = &partialSums(0,k);
                for( Int i=sumIndices(k+1); i>=k+1; --i )
                    s[i] = s[i+1] + N(k,i)*vBuf[i];

                centers(k) = -partialSums(k+1,k);
                vBuf[k] = Round(centers(k));
                spiralStates[k].Initialize( centers(k) );
            }
        }
        else
        {
            // Move up the tree
            ++k;
            if( k == n )
            {
                // Return an arbitrary value greater than upperBounds(n-1)
                return 2*upperBounds(n-1)+1;
            }
            sumIndices(k) = k; // indicate that (i,j) are not synchronized
            if( k > lastNonzero )
            {
                // NOTE: We should have k == lastNonzero+1 

                // Seed a constrained spiral out from zero
                spiralStates[k].Initialize( true );
                vBuf[k] = spiralStates[k].Step();
                lastNonzero = k;
                if( ctrl.innerProgress )
                    Output("lastNonzero: ",lastNonzero);
            }
            else
            {
                vBuf[k] = spiralStates[k].Step();
            }
        }
    }
}

template<typename F>
Base<F> TransposedHelper
( const Matrix<Base<F>>& d,
  const Matrix<F>& NTrans,
  const Matrix<Base<F>>& upperBounds,
        Matrix<F>& v,
  const EnumCtrl<Base<F>>& ctrl )
{
    EL_DEBUG_CSE
    typedef Base<F> Real;
    const Int m = NTrans.Width();
    const Int n = NTrans.Height();
    if( n > m )
        LogicError("Expected height(N) >= width(N)");

    Matrix<F> partialSums;
    Zeros( partialSums, n+1, n );

    // The 'sumIndices' vector is indexed differently from the 'r' from
    // Gama/Nguyen/Regev, and (r_0,r_1,r_2,...,r_n)=(0,1,2,...,n) becomes
    // (-1,0,1,...,n-1).
    Matrix<Int> sumIndices;
    Zeros( sumIndices, n+1, 1 );
    for( Int j=0; j<=n; ++j )
        sumIndices(j) = j-1;

    // Note: We maintain the norms rather than their squares
    Matrix<Real> partialNorms;
    Zeros( partialNorms, n+1, 1 );

    Zeros( v, n, 1 );
    if( n == 0 )
        return Real(0);
    v(0) = F(1);
    Int lastNonzero = 0; // -1 if all sumIndices are zero

    Matrix<F> centers;
    Zeros( centers, n, 1 );

    vector<SpiralState<F>> spiralStates(n);

    Int k=0;
    F* vBuf = &v(0);
    while( true )
    {
        const F entry = d(k)*(vBuf[k] - centers(k));
        const Real partialNorm = SafeNorm( partialNorms(k+1), entry );
        partialNorms(k) = partialNorm;
        if( partialNorm < upperBounds((n-1)-k) )
        {
            if( k == 0 )
            {
                // Success
                return partialNorm;
            }
            else
            {
                // Move down the tree
                --k;
                sumIndices(k) = Max(sumIndices(k),sumIndices(k+1));

                      F* s = &partialSums(0,k);
                const F* nBuf = &NTrans(0,k);
                for( Int i=sumIndices(k+1); i>=k+1; --i )
                    s[i] = s[i+1] + nBuf[i]*vBuf[i];

                centers(k) = -partialSums(k+1,k);
                vBuf[k] = Round(centers(k));
                spiralStates[k].Initialize( centers(k) );
            }
        }
        else
        {
            // Move up the tree
            ++k;
            if( k == n )
            {
                // Return an arbitrary value > than upperBounds(n-1)
                return 2*upperBounds(n-1)+1;
            }
            sumIndices(k) = k; // indicate that (i,j) are not synchronized
            if( k > lastNonzero )
            {
                // NOTE: We should have k == lastNonzero+1 

                // Seed a constrained spiral out from zero
                spiralStates[k].Initialize( true );
                vBuf[k] = spiralStates[k].Step();
                lastNonzero = k;
                if( ctrl.innerProgress )
                    Output("lastNonzero: ",lastNonzero);
            }
            else
            {
                vBuf[k] = spiralStates[k].Step();
            }
        }
    }
}

} // namespace gnr_enum

template<typename F>
Base<F> GNREnumeration
( const Matrix<Base<F>>& d,
  const Matrix<F>& N,
  const Matrix<Base<F>>& upperBounds,
        Matrix<F>& v,
  const EnumCtrl<Base<F>>& ctrl )
{
    EL_DEBUG_CSE
    if( ctrl.explicitTranspose )
    {
        Matrix<F> NTrans;
        Transpose( N, NTrans );
        return gnr_enum::TransposedHelper( d, NTrans, upperBounds, v, ctrl );
    }
    else
    {
        return gnr_enum::Helper( d, N, upperBounds, v, ctrl );
    }
}

} // namespace svp

#define PROTO(F) \
  template Base<F> svp::GNREnumeration \
  ( const Matrix<Base<F>>& d, \
    const Matrix<F>& N, \
    const Matrix<Base<F>>& u, \
          Matrix<F>& v, \
    const EnumCtrl<Base<F>>& ctrl );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
