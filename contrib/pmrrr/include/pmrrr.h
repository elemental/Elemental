/* Copyright (c) 2010, RWTH Aachen University
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or 
 * without modification, are permitted provided that the following
 * conditions are met:
 *   * Redistributions of source code must retain the above 
 *     copyright notice, this list of conditions and the following
 *     disclaimer.
 *   * Redistributions in binary form must reproduce the above 
 *     copyright notice, this list of conditions and the following 
 *     disclaimer in the documentation and/or other materials 
 *     provided with the distribution.
 *   * Neither the name of the RWTH Aachen University nor the
 *     names of its contributors may be used to endorse or promote 
 *     products derived from this software without specific prior 
 *     written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL RWTH 
 * AACHEN UNIVERSITY BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF 
 * USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT 
 * OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF 
 * SUCH DAMAGE.
 *
 * Coded by Matthias Petschow (petschow@aices.rwth-aachen.de),
 * August 2010, Version 0.6
 *
 * This code was the result of a collaboration between 
 * Matthias Petschow and Paolo Bientinesi. When you use this 
 * code, kindly reference a paper related to this work.
 *
 */

#ifndef PPMRRR_H
#define PPMRRR_H

#include "mpi.h"
#include "global.h"

/* Parallel computation of all or a subset of eigenvalues and 
 * optionally eigenvectors of a symmetric tridiagonal matrix based on 
 * the algorithm of Multiple Relatively Robust Representations (MRRR). 
 * The routine targets hybrid architectures consisting of multiple SMP 
 * nodes. It also runs in fully distributed mode, with each node 
 * having only one processor, and fully SMP mode, in which case no 
 * message passing is required. The implementation is based on 
 * LAPACK's routine 'dstemr'.
 *
 * Function prototype: */

int pmrrr(char *jobz, char *range, int *n, double  *D,
	  double *E, double *vl, double *vu, int *il, int *iu,
	  int *tryrac, MPI_Comm comm, int *nz, int *offset,
	  double *W, double *Z, int *ldz, int *Zsupp);

/* Arguments:
 * ----------
 *
 * INPUTS: 
 * -------
 * jobz              "N" or "n" - compute only eigenvalues
 *                   "V" or "v" - compute also eigenvectors
 *                   "C" or "c" - count the maximal number of 
 *                                locally computed eigenvectors
 * range             "A" or "a" - all
 *                   "V" or "v" - by interval: (VL,VU]
 *                   "I" or "i" - by index:     IL-IU
 * n                 matrix size
 * ldz               must be set on input to the leading dimension 
 *                   of of eigenvector matrix Z; this is often equal 
 *                   to matrix size n (not changed on output)
 *
 * INPUT + OUTPUT: 
 * ---------------
 * D (double[n])     Diagonal elements of tridiagonal T.
 *                   (On output the array will be overwritten).
 * E (double[n])     Off-diagonal elements of tridiagonal T.
 *                   First n-1 elements contain off-diagonals,
 *                   the last element can have an abitrary value. 
 *                   (On output the array will be overwritten.)
 * vl                If range="V", lower bound of interval
 *                   (vl,vu], on output refined.
 *                   If range="A" or "I" not referenced as input.
 *                   On output the interval (vl,vu] contains ALL
 *                   the computed eigenvalues.
 * vu                If range="V", upper bound of interval
 *                   (vl,vu], on output refined.
 *                   If range="A" or "I" not referenced as input.
 *                   On output the interval (vl,vu] contains ALL
 *                   the computed eigenvalues.
 * il                If range="I", lower index (1-based indexing) of 
 *                   the subset 'il' to 'iu'.
 *                   If range="A" or "V" not referenced as input.
 *                   On output the eigenvalues with index il to iu are 
 *                   computed by ALL processes.
 * iu                If range="I", upper index (1-based indexing) of 
 *                   the subset 'il' to 'iu'.
 *                   If range="A" or "V" not referenced as input.
 *                   On output the eigenvalues with index il to iu are 
 *                   computed by ALL processes.
 * tryrac            0 - do not try to achieve high relative accuracy.
 *                   1 - relative accuracy will be attempted; 
 *                       on output it is set to zero if high relative 
 *                       accuracy is not achieved.
 * comm              MPI communicator; commonly: MPI_COMM_WORLD.
 *
 * OUTPUT: 
 * -------
 * nz                Number of eigenvalues and eigenvectors computed 
 *                   locally.
 *                   If jobz="C", 'nz' will be set to the maximal
 *                   number of locally computed eigenvectors such 
 *                   that double[n*nz] will provide enough memory 
 *                   for the local eigenvectors;  this is only 
 *                   important in case of range="V" since 
 *                   '#eigenpairs' are not known in advance
 * offset            Index, relative to the computed eigenvalues, of 
 *                   the smallest eigenvalue computed locally
 *                   (0-based indexing).
 * W (double[n])     Locally computed eigenvalues;
 *                   The first nz entries contain the eigenvalues 
 *                   computed locally; the first entry contains the 
 *                   'offset + 1'-th eigenvalue computed and the 
 *                   'offset + il'-th eigenvalue of the input matrix 
 *                   (1-based indexing in both cases).
 *                   In some situations it is desirable to have all 
 *                   computed eigenvalues in W, instead of only 
 *                   those computed locally. In this case the 
 *                   routine 'PMR_comm_eigvals' can be called right 
 *                   after 'pmrrr' returns (see below).
 * Z                 Locally computed eigenvectors.
 * (double[n*nz])    Enough space must be provided to store the
 *                   vectors. 'nz' should be bigger or equal 
 *                   to ceil('#eigenpairs'/'#processes'), where 
 *                   '#eigenpairs' is 'n' in case of range="A" and
 *                  'iu-il+1' in case of range="I". Alternatively, 
 *                   and for range="V" 'nz' can be obtained 
 *                   by running the routine with jobz="C". 
 * Zsupp             Support of eigenvectors, which is given by
 * (double[2*n])     Z[2*i] to Z[2*i+1] for the i-th eigenvector
 *                   stored locally (1-based indexing in both 
 *                   cases).
 *
 * RETURN VALUE: 
 * -------------
 *                 0 - success  
 *                 1 - wrong input parameter
 *                 2 - misc errors  
 *
 * The Fortran interface takes an additinal integer argument INFO
 * to retrieve the return value. 
 * An example call in Fortran looks therefore like
 *
 * CALL PMRRR('V', 'A', N, D, E, VL, VU, IL, IU, TRYRAC, 
 *            MPI_COMM_WORLD, NZ, MYFIRST, W, Z, LDZ, ZSUPP, INFO)
 *
 * Some additinal comments:
 * ------------------------
 * + In the case '#processors'='#nodes'*'#cores' << 16 the routine
 *   might be optimal in terms of performance.
 *
 */


/* Set the number of threads in case PMR_NUM_THREADS is not 
 * specified */
#define DEFAULT_NUM_THREADS 4

/* Call LAPACK's dstemr in every process to compute all desiered 
 * eigenpairs redundantly (and discard the once that would usually 
 * not be computed by the process) if n < DSTEMR_IF_SMALLER; 
 * default: 4 */ 
#define DSTEMR_IF_SMALLER   4

/* Make sure that eigenpairs are sorted globally; if set to false
 * they are in most cases sorted, but it is not double checked and 
 * can therefore not be guaranteed; default: true */
#define ASSERT_SORTED_EIGENPAIRS true

/* Set flag if Rayleigh Quotient Correction should be used, 
 * which is usually faster; default: true */
#define TRY_RQC          true

/* Maximum numver of iterations of inverse iteration;
 * default: 10 */
#define MAXITER            10

/* Set the min. relative gap for an eigenvalue to be considered 
 * well separated, that is a singleton; this is a very important 
 * parameter of the computation; default: 10e-3 */
#define MIN_RELGAP       1e-3

/* Set the maximal allowed element growth for being accepted as 
 * an RRR, that is if max. pivot < MAX_GROWTH * 'spectral diameter'
 * the RRR is accepted; default: 64.0 */
#define MAX_GROWTH         64.0

/* Set how many iterations should be executed to find the root 
 * representation; default: 6 */
#define MAX_TRY_RRRR       10


/*
 * Routine to communicate eigenvalues such that every process has
 * all computed eigenvalues (iu-il+1) in W; this routine is designed 
 * to be called right after 'pmrrr'.
 */
int PMR_comm_eigvals(MPI_Comm comm, int *nz, int *ifirst, double *W);
/* Arguments:
 * ----------
 *
 * INPUTS: 
 * -------
 * comm              MPI communicator; commonly: MPI_COMM_WORLD.
 * nz                Number of eigenvalues local in W as returned 
 *                   from 'pmrrr'.
 * offset            Index, relative to the computed eigenvalues, of 
 *                   the smallest eigenvalue computed locally
 *                   (0-based indexing).
 *
 * INPUT + OUTPUT: 
 * ---------------
 * W (double[n])     Eigenvalues. 
 *                   On input the first nz elements of W contain 
 *                   the eigenvalues computed locally.
 *                   On output the first 'iu-il+1' of W contain 
 *                   all computed eigenvalues.
 *
 * The Fortran interface takes an additinal integer argument INFO 
 * to retrieve the return value. 
 * An example call in Fortran looks therefore like
 *
 * CALL PMR_COMM_EIGVALS(MPI_COMM_WORLD, NZ, OFFSET, W, INFO)
 *
 */

/* BLAS function prototypes
 * 
 * These macros were added by Jack Poulson for easing integration into 
 * Elemental */
#if defined(BLAS_POST)
# if defined(PREPEND_X)
#  define BLAS(name) x ## name ## _
# else
#  define BLAS(name) name ## _
# endif
#else
# if defined(PREPEND_X)
#  define BLAS(name) x ## name
# else
#  define BLAS(name) name
# endif
#endif
void pmrrr_dscal(int*, double*, double*, int*);
void BLAS(dscal)(int*, double*, double*, int*);

/* LAPACK function prototypes
 * Note: type specifier 'extern' does not matter in declaration
 * so here used to mark routines from LAPACK and BLAS libraries 
 * 
 * The macros for determining underscore convention were added by 
 * Jack Poulson to ease integration into Elemental */
#if defined(LAPACK_POST)
# if defined(PREPEND_X)
#  define LAPACK(name) x ## name ## _
# else
#  define LAPACK(name) name ## _
# endif
#else
# if defined(PREPEND_X)
#  define LAPACK(name) x ## name
# else
#  define LAPACK(name) name
# endif
#endif
double LAPACK(dlamch)(char*);

double LAPACK(dlanst)(char*, int*, double*, double*);

void LAPACK(dlarrr)(int*, double*, double*, int*);

void LAPACK(dlarra)
(int*, double*, double*, double*, double*,double*, int*, int*, int*);

void LAPACK(dlarrc)
(char*, int*, double*, double*, double*, double*, double*, int*, int*, int*,
 int*);

void LAPACK(dlarrd)
(char*, char*, int*, double*, double*, int*, int*, double*, double*, double*,
 double*, double*, double*, int*, int*, int*, double*, double*, double*, 
 double*, int*, int*, double*, int*, int*);

void LAPACK(dlarrb)
(int*, double*, double*, int*, int*, double*, double*, int*, double*, double*,
 double*, double*, int*, double*, double*, int*, int*);

void LAPACK(dlarrk)
(int*, int*, double*, double*, double*, double*, double*, double*, double*, 
 double*, int*);

void LAPACK(dlaebz)
(int*, int*, int*, int*, int*, int*, double*, double*, double*, double*, 
 double*, double*, int*, double*, double*, int*, int*, double*, int*, int*);

void LAPACK(dlarnv)(int*, int*, int*, double*);

void LAPACK(dlarrf)
(int*, double*, double*, double*, int*, int*, double*, double*, double*, 
 double*, double*, double*, double*, double*, double*, double*, double*, 
 int*);

void LAPACK(dlar1v)
(int*, int*, int*, double*, double*, double*, double*, double*, double*, 
 double*, double*, bool*, int*, double*, double*, int*, int*, double*, 
 double*, double*, double*);

void LAPACK(dlarrj)
(int*, double*, double*, int*, int*, double*, int*, double*, double*, double*,
 int*, double*, double*, int*);

void LAPACK(dstemr)
(char*, char*, int*, double*, double*, double*, double*, int*, int*, int*, 
 double*, double*, int*, int*, int*, int*, double*, int*, int*, int*, int*);

#endif /* End of header file */
