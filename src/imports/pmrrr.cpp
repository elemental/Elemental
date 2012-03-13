/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

    - Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

    - Neither the name of the owner nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/
#include "elemental/core/environment.hpp"

#ifndef WITHOUT_PMRRR

extern "C" {

int PMRRR
( const char* jobz,  // 'N' ~ only eigenvalues, 'V' ~ also eigenvectors
  const char* range, // 'A'~all eigenpairs, 'V'~interval (vl,vu], 'I'~il-iu
  const int* n,      // size of matrix
  const double* d,   // full diagonal of tridiagonal matrix [length n]
  const double* e,   // full subdiagonal in first n-1 entries [length n]
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

namespace elem {
namespace pmrrr {

// Return upper bounds on the number of (local) eigenvalues in the given range,
// (lowerBound,upperBound]
Estimate EigEstimate
( int n, const double* d, const double* e, double* w, 
  mpi::Comm comm, double lowerBound, double upperBound )
{
#ifndef RELEASE
    PushCallStack("pmrrr::EigEstimate");
#endif
    char jobz='C';
    char range='V';
    int il, iu;
    int highAccuracy=0;
    int nz, offset;
    int ldz=1;
    std::vector<int> ZSupport(2*n);
    int retval = PMRRR
    ( &jobz, &range, &n, d, e, &lowerBound, &upperBound, &il, &iu, 
      &highAccuracy, comm, &nz, &offset, w, 0, &ldz, &ZSupport[0] );
    if( retval != 0 )
    {
        std::ostringstream msg;
        msg << "PMRRR returned " << retval; 
        throw std::runtime_error( msg.str().c_str() );
    }

    Estimate estimate;
    estimate.numLocalEigenvalues = nz;
    mpi::AllReduce( &nz, &estimate.numGlobalEigenvalues, 1, mpi::SUM, comm );

#ifndef RELEASE
    PopCallStack();
#endif
    return estimate;
}

// Compute all of the eigenvalues
Info Eig( int n, const double* d, const double* e, double* w, mpi::Comm comm )
{
#ifndef RELEASE
    PushCallStack("pmrrr::Eig");
#endif
    char jobz='N';
    char range='A';
    double vl, vu;
    int il, iu;
    int highAccuracy=0; 
    int nz, offset;
    int ldz=1;
    std::vector<int> ZSupport(2*n);
    int retval = PMRRR
    ( &jobz, &range, &n, d, e, &vl, &vu, &il, &iu, &highAccuracy, comm,
      &nz, &offset, w, 0, &ldz, &ZSupport[0] );
    if( retval != 0 )
    {
        std::ostringstream msg;        
        msg << "PMRRR returned " << retval;
        throw std::runtime_error( msg.str().c_str() );
    }

    Info info;
    info.numLocalEigenvalues=nz;
    info.firstLocalEigenvalue=offset;
    info.numGlobalEigenvalues=n;

#ifndef RELEASE
    PopCallStack();
#endif
    return info;
}

// Compute all of the eigenpairs
Info Eig
( int n, const double* d, const double* e, double* w, double* Z, int ldz, 
  mpi::Comm comm )
{
#ifndef RELEASE
    PushCallStack("pmrrr::Eig");
#endif
    char jobz='V';
    char range='A';
    double vl, vu;
    int il, iu;
    int highAccuracy=0; 
    int nz, offset;
    std::vector<int> ZSupport(2*n);
    int retval = PMRRR
    ( &jobz, &range, &n, d, e, &vl, &vu, &il, &iu, &highAccuracy, comm,
      &nz, &offset, w, Z, &ldz, &ZSupport[0] );
    if( retval != 0 )
    {
        std::ostringstream msg;        
        msg << "PMRRR returned " << retval;
        throw std::runtime_error( msg.str().c_str() );
    }

    Info info;
    info.numLocalEigenvalues=nz;
    info.firstLocalEigenvalue=offset;
    info.numGlobalEigenvalues=n;

#ifndef RELEASE
    PopCallStack();
#endif
    return info;
}

// Compute all of the eigenvalues in (lowerBound,upperBound]
Info Eig
( int n, const double* d, const double* e, double* w, 
  mpi::Comm comm, double lowerBound, double upperBound )
{
#ifndef RELEASE
    PushCallStack("pmrrr::Eig");
#endif
    char jobz='N';
    char range='V';
    int il, iu;
    int highAccuracy=0; 
    int nz, offset;
    int ldz=1;
    std::vector<int> ZSupport(2*n);
    int retval = PMRRR
    ( &jobz, &range, &n, d, e, &lowerBound, &upperBound, &il, &iu, 
      &highAccuracy, comm, &nz, &offset, w, 0, &ldz, &ZSupport[0] );
    if( retval != 0 )
    {
        std::ostringstream msg;        
        msg << "PMRRR returned " << retval;
        throw std::runtime_error( msg.str().c_str() );
    }

    Info info;
    info.numLocalEigenvalues=nz;
    info.firstLocalEigenvalue=offset;
    mpi::AllReduce( &nz, &info.numGlobalEigenvalues, 1, mpi::SUM, comm );

#ifndef RELEASE
    PopCallStack();
#endif
    return info;
}

// Compute all of the eigenpairs with eigenvalues in (lowerBound,upperBound]
Info Eig
( int n, const double* d, const double* e, double* w, double* Z, int ldz, 
  mpi::Comm comm, double lowerBound, double upperBound )
{
#ifndef RELEASE
    PushCallStack("pmrrr::Eig");
#endif
    char jobz='V';
    char range='V';
    int il, iu;
    int highAccuracy=0; 
    int nz, offset;
    std::vector<int> ZSupport(2*n);
    int retval = PMRRR
    ( &jobz, &range, &n, d, e, &lowerBound, &upperBound, &il, &iu, 
      &highAccuracy, comm, &nz, &offset, w, Z, &ldz, &ZSupport[0] );
    if( retval != 0 )
    {
        std::ostringstream msg;        
        msg << "PMRRR returned " << retval;
        throw std::runtime_error( msg.str().c_str() );
    }

    Info info;
    info.numLocalEigenvalues=nz;
    info.firstLocalEigenvalue=offset;
    mpi::AllReduce( &nz, &info.numGlobalEigenvalues, 1, mpi::SUM, comm );

#ifndef RELEASE
    PopCallStack();
#endif
    return info;
}

// Compute all of the eigenvalues with indices in [lowerBound,upperBound]
Info Eig
( int n, const double* d, const double* e, double* w, 
  mpi::Comm comm, int lowerBound, int upperBound )
{
#ifndef RELEASE
    PushCallStack("pmrrr::Eig");
#endif
    ++lowerBound;
    ++upperBound;
    char jobz='N';
    char range='I';
    double vl, vu;
    int highAccuracy=0; 
    int nz, offset;
    int ldz=1;
    std::vector<int> ZSupport(2*n);
    int retval = PMRRR
    ( &jobz, &range, &n, d, e, &vl, &vu, &lowerBound, &upperBound, 
      &highAccuracy, comm, &nz, &offset, w, 0, &ldz, &ZSupport[0] );
    if( retval != 0 )
    {
        std::ostringstream msg;        
        msg << "PMRRR returned " << retval;
        throw std::runtime_error( msg.str().c_str() );
    }

    Info info;
    info.numLocalEigenvalues=nz;
    info.firstLocalEigenvalue=offset;
    info.numGlobalEigenvalues=(upperBound-lowerBound)+1;

#ifndef RELEASE
    PopCallStack();
#endif
    return info;
}

// Compute all of the eigenpairs with eigenvalues indices in 
// [lowerBound,upperBound]
Info Eig
( int n, const double* d, const double* e, double* w, double* Z, int ldz, 
  mpi::Comm comm, int lowerBound, int upperBound )
{
#ifndef RELEASE
    PushCallStack("pmrrr::Eig");
#endif
    ++lowerBound;
    ++upperBound;
    char jobz='V';
    char range='I';
    double vl, vu;
    int highAccuracy=0; 
    int nz, offset;
    std::vector<int> ZSupport(2*n);
    int retval = PMRRR
    ( &jobz, &range, &n, d, e, &vl, &vu, &lowerBound, &upperBound, 
      &highAccuracy, comm, &nz, &offset, w, Z, &ldz, &ZSupport[0] );
    if( retval != 0 )
    {
        std::ostringstream msg;        
        msg << "PMRRR returned " << retval;
        throw std::runtime_error( msg.str().c_str() );
    }

    Info info;
    info.numLocalEigenvalues=nz;
    info.firstLocalEigenvalue=offset;
    info.numGlobalEigenvalues=(upperBound-lowerBound)+1;

#ifndef RELEASE
    PopCallStack();
#endif
    return info;
}

} // namespace pmrrr
} // namespace elem

#endif // WITHOUT_PMRRR
