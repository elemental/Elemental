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
#ifndef PMRRR_HPP
#define PMRRR_HPP 1

namespace elem {
namespace pmrrr {

struct Estimate {
    int numLocalEigenvalues;
    int numGlobalEigenvalues;
};

// Return an upper bound on the number of (local) eigenvalues in the given range
Estimate EigEstimate
( int n, const double* d, const double* e, double* w, 
  mpi::Comm comm, double lowerBound, double upperBound );


struct Info {
    int numLocalEigenvalues;
    int numGlobalEigenvalues;

    int firstLocalEigenvalue;
};

// Compute all of the eigenvalues
Info Eig
( int n, const double* d, const double* e, double* w, mpi::Comm comm );

// Compute all of the eigenpairs
Info Eig
( int n, const double* d, const double* e, double* w, double* Z, int ldz, 
  mpi::Comm comm );

// Compute all of the eigenvalues in [lowerBound,upperBound)
Info Eig
( int n, const double* d, const double* e, double* w,  
  mpi::Comm comm, double lowerBound, double upperBound );

// Compute all of the eigenpairs with eigenvalues in [lowerBound,upperBound)
Info Eig
( int n, const double* d, const double* e, double* w, double* Z, int ldz, 
  mpi::Comm comm, double lowerBound, double upperBound );

// Compute all of the eigenvalues with indices in [lowerBound,upperBound)
Info Eig
( int n, const double* d, const double* e, double* w,
  mpi::Comm comm, int lowerBound, int upperBound );

// Compute all of the eigenpairs with ordered eigenvalue indices in 
// [lowerBound,upperBound)
Info Eig
( int n, const double* d, const double* e, double* w,
  double* Z, int ldz, mpi::Comm comm, int lowerBound, int upperBound );

} // namespace pmrrr
} // namespace elem

#endif // PMRRR_HPP
