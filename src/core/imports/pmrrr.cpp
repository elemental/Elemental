/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El-lite.hpp>

extern "C" {

int pmrrr
( const char* jobz,  // 'N' ~ only eigenvalues, 'V' ~ also eigenvectors
  const char* range, // 'A'~all eigenpairs, 'V'~interval (vl,vu], 'I'~il-iu
  const int* n,      // size of matrix
        double* d,   // full diagonal of tridiagonal matrix [length n]
        double* e,   // full subdiagonal in first n-1 entries [length n]
  const double* vl,  // if range=='V', compute eigenpairs in (vl,vu]
  const double* vu,
  const int* il, // if range=='I', compute il-iu eigenpairs
  const int* iu,
  int* tryrac, // if nonzero, try for high relative accuracy
  MPI_Comm comm,
  int* nz,        // number of locally computed eigenvectors
  int* offset,    // the first eigenpair computed by our process
  double* w,      // eigenvalues corresponding to local eigenvectors [length nz]
  double* Z,      // local eigenvectors [size ldz x nz]
  const int* ldz, // leading dimension of Z
  int* ZSupp      // support of eigenvectors [length 2n]
);

} // extern "C"

namespace El {
namespace herm_tridiag_eig {

// Return upper bounds on the number of (local) eigenvalues in the given range,
// (lowerBound,upperBound]
Estimate EigEstimate
( int n, double* d, double* e, double* w, mpi::Comm comm, 
  double lowerBound, double upperBound )
{
    EL_DEBUG_CSE
    Estimate estimate;
    char jobz='C';
    char range='V';
    int il, iu;
    int highAccuracy=0;
    int nz, offset;
    int ldz=1;
    vector<int> ZSupport(2*n);
    int retval = pmrrr
    ( &jobz, &range, &n, d, e, &lowerBound, &upperBound, &il, &iu, 
      &highAccuracy, comm.comm, &nz, &offset, w, 0, &ldz, ZSupport.data() );
    if( retval != 0 )
        RuntimeError("pmrrr returned ",retval);

    estimate.numLocalEigenvalues = nz;
    estimate.numGlobalEigenvalues = mpi::AllReduce( nz, comm );
    return estimate;
}

// Compute all of the eigenvalues
Info Eig( int n, double* d, double* e, double* w, mpi::Comm comm )
{
    EL_DEBUG_CSE
    Info info;
    char jobz='N';
    char range='A';
    double vl, vu;
    int il, iu;
    int highAccuracy=0; 
    int nz, offset;
    int ldz=1;
    vector<int> ZSupport(2*n);
    int retval = pmrrr
    ( &jobz, &range, &n, d, e, &vl, &vu, &il, &iu, &highAccuracy, comm.comm,
      &nz, &offset, w, 0, &ldz, ZSupport.data() );
    if( retval != 0 )
        RuntimeError("pmrrr returned ",retval);

    info.numLocalEigenvalues=nz;
    info.firstLocalEigenvalue=offset;
    info.numGlobalEigenvalues=n;
    return info;
}

// Compute all of the eigenpairs
Info Eig
( int n, double* d, double* e, double* w, double* Z, int ldz, mpi::Comm comm )
{
    EL_DEBUG_CSE
    Info info;
    char jobz='V';
    char range='A';
    double vl, vu;
    int il, iu;
    int highAccuracy=0; 
    int nz, offset;
    vector<int> ZSupport(2*n);
    int retval = pmrrr
    ( &jobz, &range, &n, d, e, &vl, &vu, &il, &iu, &highAccuracy, comm.comm,
      &nz, &offset, w, Z, &ldz, ZSupport.data() );
    if( retval != 0 )
        RuntimeError("pmrrr returned ",retval);

    info.numLocalEigenvalues=nz;
    info.firstLocalEigenvalue=offset;
    info.numGlobalEigenvalues=n;
    return info;
}

// Compute all of the eigenvalues in (lowerBound,upperBound]
Info Eig
( int n, double* d, double* e, double* w, mpi::Comm comm, 
  double lowerBound, double upperBound )
{
    EL_DEBUG_CSE
    Info info;
    char jobz='N';
    char range='V';
    int il, iu;
    int highAccuracy=0; 
    int nz, offset;
    int ldz=1;
    vector<int> ZSupport(2*n);
    int retval = pmrrr
    ( &jobz, &range, &n, d, e, &lowerBound, &upperBound, &il, &iu, 
      &highAccuracy, comm.comm, &nz, &offset, w, 0, &ldz, ZSupport.data() );
    if( retval != 0 )
        RuntimeError("pmrrr returned ",retval);

    info.numLocalEigenvalues=nz;
    info.firstLocalEigenvalue=offset;
    mpi::AllReduce( &nz, &info.numGlobalEigenvalues, 1, mpi::SUM, comm );
    return info;
}

// Compute all of the eigenpairs with eigenvalues in (lowerBound,upperBound]
Info Eig
( int n, double* d, double* e, double* w, double* Z, int ldz, mpi::Comm comm, 
  double lowerBound, double upperBound )
{
    EL_DEBUG_CSE
    Info info;
    char jobz='V';
    char range='V';
    int il, iu;
    int highAccuracy=0; 
    int nz, offset;
    vector<int> ZSupport(2*n);
    int retval = pmrrr
    ( &jobz, &range, &n, d, e, &lowerBound, &upperBound, &il, &iu, 
      &highAccuracy, comm.comm, &nz, &offset, w, Z, &ldz, ZSupport.data() );
    if( retval != 0 )
        RuntimeError("pmrrr returned ",retval);

    info.numLocalEigenvalues=nz;
    info.firstLocalEigenvalue=offset;
    mpi::AllReduce( &nz, &info.numGlobalEigenvalues, 1, mpi::SUM, comm );
    return info;
}

// Compute all of the eigenvalues with indices in [lowerBound,upperBound]
Info Eig
( int n, double* d, double* e, double* w, mpi::Comm comm, 
  int lowerBound, int upperBound )
{
    EL_DEBUG_CSE
    Info info;
    ++lowerBound;
    ++upperBound;
    char jobz='N';
    char range='I';
    double vl, vu;
    int highAccuracy=0; 
    int nz, offset;
    int ldz=1;
    vector<int> ZSupport(2*n);
    int retval = pmrrr
    ( &jobz, &range, &n, d, e, &vl, &vu, &lowerBound, &upperBound, 
      &highAccuracy, comm.comm, &nz, &offset, w, 0, &ldz, ZSupport.data() );
    if( retval != 0 )
        RuntimeError("pmrrr returned ",retval);

    info.numLocalEigenvalues=nz;
    info.firstLocalEigenvalue=offset;
    info.numGlobalEigenvalues=(upperBound-lowerBound)+1;
    return info;
}

// Compute all of the eigenpairs with eigenvalues indices in 
// [lowerBound,upperBound]
Info Eig
( int n, double* d, double* e, double* w, double* Z, int ldz, mpi::Comm comm, 
  int lowerBound, int upperBound )
{
    EL_DEBUG_CSE
    Info info;
    ++lowerBound;
    ++upperBound;
    char jobz='V';
    char range='I';
    double vl, vu;
    int highAccuracy=0; 
    int nz, offset;
    vector<int> ZSupport(2*n);
    int retval = pmrrr
    ( &jobz, &range, &n, d, e, &vl, &vu, &lowerBound, &upperBound, 
      &highAccuracy, comm.comm, &nz, &offset, w, Z, &ldz, ZSupport.data() );
    if( retval != 0 )
        RuntimeError("pmrrr returned ",retval);

    info.numLocalEigenvalues=nz;
    info.firstLocalEigenvalue=offset;
    info.numGlobalEigenvalues=(upperBound-lowerBound)+1;
    return info;
}

} // namespace herm_tridiag_eig
} // namespace El
