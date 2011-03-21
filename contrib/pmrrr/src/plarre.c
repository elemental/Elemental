/* Parallel computation of eigenvalues and symmetric tridiagonal 
 * matrix T, given by its diagonal elements D and its super-/sub-
 * diagonal elements E.
 *
 * Copyright (c) 2010, RWTH Aachen University
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include <pthread.h>
#include "mpi.h"
#include "pmrrr.h" 
#include "plarre.h"
#include "global.h"
#include "structs.h" 


#define ONE                1.0
#define HUNDRED          100.0
#define HALF               0.5
#define FOURTH             0.25


static void *eigval_subset_thread_a(void *argin);
static void *eigval_subset_thread_r(void *argin);
static void clean_up_plarre(double*, double*, int*, int*, int*);

static 
int eigval_subset_proc(proc_t*, char*, in_t*, double*, int, int,
		       tol_t*, val_t*, double*, int*);
static 
auxarg1_t *create_auxarg1(int, double*, double*, double*, int, int, 
			  int, int, int, int*, double,  double,
			  double*, double*, double*, int*, int*);
static 
void retrieve_auxarg1(auxarg1_t*, int*, double**, double**, double**,
		      int*, int*, int*, int*, int*, int**, double*, 
		      double*, double**, double**, double**, int**, 
		      int**);
static 
auxarg2_t *create_auxarg2(int, double*, double*, int, int, val_t*,
			  double, double, double, double);
static 
void retrieve_auxarg2(auxarg2_t*, int*, double**, double**, int*, 
		      int*, val_t**, double*, double*, double*, 
		      double*);





/* Routine to compute eigenvalues */
int plarre(proc_t *procinfo, char *jobz, char *range, in_t *Dstruct, 
	   val_t *Wstruct, tol_t *tolstruct, int *nzp, int *myfirstp)
{
  /* input variables */
  int              pid    = procinfo->pid;
  int              nproc  = procinfo->nproc;
  bool             wantZ  = (jobz[0]  == 'V' || jobz[0]  == 'v');
  bool             cntval = (jobz[0]  == 'C' || jobz[0]  == 'c');
  int              n      = Dstruct->n;
  double *restrict D      = Dstruct->D;
  double *restrict E      = Dstruct->E;
  int    *restrict isplit = Dstruct->isplit;
  double           *vl    = Wstruct->vl;
  double           *vu    = Wstruct->vu;
  int              *il    = Wstruct->il;
  int              *iu    = Wstruct->iu;
  double *restrict W      = Wstruct->W;
  double *restrict Werr   = Wstruct->Werr;
  double *restrict Wgap   = Wstruct->Wgap;
  int    *restrict Windex = Wstruct->Windex;
  int    *restrict iblock = Wstruct->iblock;
  double *restrict gersch = Wstruct->gersch;

  /* constants */
  int             IZERO = 0,   IONE = 1;
  double          DZERO = 0.0;

  /* work space */
  double          *E2, *work;
  int             *iwork;

  /* compute geschgorin disks and spectral diameter */
  double          gl, gu, bl_gu, eold, emax, eabs;

  /* compute splitting points */
  int             bl_begin, bl_end;

  /* distribute work among processes */
  int             ifirst, ilast, ifirst_tmp, ilast_tmp;
  int             chunk, isize, iil, iiu;

  /* gather results */
  int             *rcount, *rdispl;

  /* others */
  int             info, i, j, im, idummy, ind;
  double          tmp1, dummy;
  enum range_enum {allrng=1, valrng=2, indrng=3} irange;
  double          intervals[2];
  int             negcounts[2];
  double          sigma;

  if (range[0] == 'A' || range[0] == 'a') {
    irange = allrng;
  } else if (range[0] == 'V' || range[0] == 'v') {
    irange = valrng;
  } else if (range[0] == 'I' || range[0] == 'i') {
    irange = indrng;
  } else {
    return(1);
  }

  /* allocate work space */
  E2     = (double *) malloc(     n * sizeof(double) );
  assert(E2 != NULL);
  work   = (double *) malloc(   4*n * sizeof(double) );
  assert(work != NULL);
  iwork  = (int *)    malloc(   3*n * sizeof(int) );
  assert(iwork != NULL);
  rcount = (int *)    malloc( nproc * sizeof(int) );
  assert(rcount != NULL);
  rdispl = (int *)    malloc( nproc * sizeof(int) );
  assert(rdispl != NULL);

  /* Compute square of off-diagonal elements */
  for (i=0; i<n-1; i++) {
    E2[i] = E[i]*E[i];
  }

  /* compute geschgorin disks and spectral diameter */
  gl     = D[0];
  gu     = D[0];
  eold   =  0.0;
  emax   =  0.0;
  E[n-1] =  0.0;

  for (i=0; i<n; i++) {
    eabs = fabs(E[i]);
    if (eabs >= emax) emax = eabs;
    tmp1 = eabs + eold;
    gersch[2*i] = D[i] - tmp1;
    gl = fmin(gl, gersch[2*i]);
    gersch[2*i+1] = D[i] + tmp1;
    gu = fmax(gu, gersch[2*i+1]);
    eold = eabs;
  }
  /* min. pivot allowed in the Sturm sequence of T */
  tolstruct->pivmin = DBL_MIN * fmax(1.0, emax*emax);
  /* estimate of spectral diameter */
  Dstruct->spdiam = gu - gl;

  /* compute splitting points with threshold "split" */
  LAPACK(dlarra)
  (&n, D, E, E2, &tolstruct->split, &Dstruct->spdiam, &Dstruct->nsplit, 
   isplit, &info);
  assert(info == 0);

  if (irange == allrng || irange == indrng) {
    *vl = gl;
    *vu = gu;
  }

  /* set eigenvalue indices in case of all or subset by value has
   * to be computed; thereby convert all problem to subset by index
   * computation */
  if (irange == allrng) {
    *il = 1;
    *iu = n;
  } else if (irange == valrng) {
    intervals[0] = *vl; intervals[1] = *vu;
    
    /* find negcount at boundaries 'vl' and 'vu'; 
     * needs work of dim(n) and iwork of dim(n) */
    LAPACK(dlaebz)
    (&IONE, &IZERO, &n, &IONE, &IONE, &IZERO, &DZERO, &DZERO, 
     &tolstruct->pivmin, D, E, E2, &idummy, intervals, &dummy, &idummy, 
     negcounts, work, iwork, &info);
    assert(info == 0);
    
    /* update negcounts of whole matrix with negcounts found for block */
    *il = negcounts[0] + 1;
    *iu = negcounts[1];
  }

  if (cntval && irange == valrng) {
    /* clean up and return */
    *nzp = iceil(*iu-*il+1, nproc);
    clean_up_plarre(E2, work, iwork, rcount, rdispl);
    return(0);
  }

  /* in case only eigenvalues are desired compute eigenvalues 
   * "il" to "iu"; otherwise compute all */
  if (wantZ) {
    iil = 1;
    iiu = n;
  } else {
    iil = *il;
    iiu = *iu;
  }
  
  /* each process computes a subset of the eigenvalues */
  ifirst_tmp = iil;
  for (i=0; i<nproc; i++) {
    chunk  = (iiu-iil+1)/nproc + (i < (iiu-iil+1)%nproc);
    if (i == nproc-1) {
      ilast_tmp = iiu;
    } else {
      ilast_tmp = ifirst_tmp + chunk - 1;
      ilast_tmp = imin(ilast_tmp, iiu);
    }
    if (i == pid) {
      ifirst    = ifirst_tmp;
      ilast     = ilast_tmp;
      isize     = ilast - ifirst + 1;
      *myfirstp = ifirst - iil;;
      *nzp      = isize;
    }
    rcount[i]  = ilast_tmp - ifirst_tmp + 1;
    rdispl[i]  = ifirst_tmp - iil;
    ifirst_tmp = ilast_tmp + 1;
    ifirst_tmp = imin(ifirst_tmp, iiu + 1);
  }

  /* compute eigenvalues assigned to process */
  if (isize != 0) {
    info = eigval_subset_proc(procinfo, range, Dstruct, E2, ifirst, ilast,
			      tolstruct, Wstruct, work, iwork);
    assert(info == 0);
  }

  if (wantZ) {
    /* communicate results */
    memcpy(work, W, isize * sizeof(double) );
    MPI_Allgatherv(work, isize, MPI_DOUBLE, W, rcount, rdispl, 
		   MPI_DOUBLE, procinfo->comm);

    memcpy(work, Werr, isize * sizeof(double) );
    MPI_Allgatherv(work, isize, MPI_DOUBLE, Werr, rcount, rdispl, 
		   MPI_DOUBLE, procinfo->comm);
    
    memcpy(iwork, Windex, isize * sizeof(int) );
    MPI_Allgatherv(iwork, isize, MPI_INT, Windex, rcount, rdispl, 
		   MPI_INT, procinfo->comm);
    
    memcpy(iwork, iblock, isize * sizeof(int) );
    MPI_Allgatherv(iwork, isize, MPI_INT, iblock, rcount, rdispl, 
		   MPI_INT, procinfo->comm);

    /* sort by block */
    memcpy(&work[0],   W,      n*sizeof(double));
    memcpy(&work[n],   Werr,   n*sizeof(double));
    memcpy(&iwork[0],  Windex, n*sizeof(int));
    memcpy(&iwork[n],  iblock, n*sizeof(int));
    
    im = 0;
    for (i=1; i<=Dstruct->nsplit; i++) {
      for (j=0; j<n; j++) {
	if (iwork[j+n] == i) {    /* iblock == i */
	  W[im]      = work[j];
	  Werr[im]   = work[j+n];
	  Windex[im] = iwork[j];
	  iblock[im] = iwork[j+n];
	  im++;
	}
      }
    }
    
    /* recompute gap of blocks */
    bl_begin = 0;
    for (i=0; i < Dstruct->nsplit; i++) {
      bl_end  = isplit[i] - 1;
      sigma   = E[bl_end];
      
      /* find outer bounds GU for block used for last gap */
      bl_gu = D[bl_begin];
      for (j = bl_begin; j <= bl_end; j++) {
	bl_gu = fmax(bl_gu, gersch[2*j+1]);
      }
      
      /* recompute gaps within the blocks */
      for (j = bl_begin; j < bl_end; j++) {
	Wgap[j] = fmax(0.0, (W[j+1] - Werr[j+1]) - (W[j] + Werr[j]) );
      }
      Wgap[bl_end] = fmax(0.0, (bl_gu - sigma) - (W[bl_end] + Werr[bl_end]) );
      
      bl_begin = bl_end + 1;
    } /* end i */
 
  } else {

    /* compute UNSHIFTED eigenvalues */
    for (i=0; i < isize; i++) {
      ind     = iblock[i]   - 1;
      bl_end  = isplit[ind] - 1;
      sigma   = E[bl_end];
      W[i]   += sigma;
    }

  } /* if wantZ */

  /* free memory */
  clean_up_plarre(E2, work, iwork, rcount, rdispl);

  return(0);
}




/*
 * Free's on allocated memory of plarre routine
 */
static  
void clean_up_plarre(double *E2, double *work, int *iwork, 
		     int *rcount, int *rdispl)
{
  free(E2);
  free(work);
  free(iwork);
  free(rcount);
  free(rdispl);
}




static 
int eigval_subset_proc(proc_t *procinfo, char *range, in_t *Dstruct,
		       double *E2, int ifirst, int ilast, tol_t *tolstruct,
		       val_t *Wstruct, double *work,
		       int *iwork)
{
  /* Input parameter */
  int              max_nthreads = procinfo->nthreads;
  int              n            = Dstruct->n;
  double *restrict D            = Dstruct->D;
  double *restrict E            = Dstruct->E;
  int              nsplit       = Dstruct->nsplit;
  int    *restrict isplit       = Dstruct->isplit;
  int              isize        = ilast-ifirst+1;
  double *restrict W            = Wstruct->W;
  double *restrict Werr         = Wstruct->Werr;
  double *restrict Wgap         = Wstruct->Wgap;
  int    *restrict Windex       = Wstruct->Windex;
  int    *restrict iblock       = Wstruct->iblock;
  double *restrict gersch       = Wstruct->gersch;
  double           pivmin       = tolstruct->pivmin;

  double gl, gu, wl, wu;

  /* Tolerances */
  double bsrtol, rtl;

  /* Multithreading */
  int            nthreads;
  int            iifirst, iilast, chunk;
  pthread_t      *threads;
  pthread_attr_t attr;
  auxarg1_t      *auxarg1;
  auxarg2_t      *auxarg2;
  void           *status;

  /* Create random vector to perturb rrr, same seed */
  int    two_n = 2*n;
  int    iseed[4] = {1,1,1,1};
  double *randvec;

  /* loop over blocks */
  int    jbl, num_vals;
  int    bl_begin,  bl_end, bl_size;
  int    bl_Wbegin, bl_Wend;
  double isleft, isright, spdiam;
  int    i_low, i_upp, bl_m;
  double sigma, s1, s2;
  int    sgndef, cnt, negcnt_lft, negcnt_rgt;
  double tau;

  /* Compute RRR */
  int    jtry, off_L, off_invD;
  double Dpivot, Dmax;
  bool   noREP;

  /* Refine eigenvalues */
  int    off_DE2, offset;
  int    rf_begin, rf_end;

  /* Others */
  int    IONE = 1, ITWO = 2;
  int    info, m, i, j;
  double dummy, tmp, tmp1, tmp2;

  /* Allocate workspace */
  randvec = (double *) malloc( 2*n * sizeof(double) );
  assert(randvec != NULL);

  threads = (pthread_t *) malloc( max_nthreads * sizeof(pthread_t) );
  assert(threads != NULL);
  
  if (max_nthreads > 1) {
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM);
  }

  /* Set tolerance parameters */
  bsrtol = sqrt(DBL_EPSILON);
  rtl    = sqrt(DBL_EPSILON);

  /* create random vector to perturb rrr and broadcast it */
  LAPACK(dlarnv)(&ITWO, iseed, &two_n, randvec);

  /* compute approximations of the eigenvalues with muliple threads
   * equivalent to:
   * LAPACK(dlarrd)
   * ("I", "B", &n, &dummy, &dummy, &ifirst, &ilast, gersch,
   *  &bsrtol, D, E, E2, &pivmin, &nsplit, isplit, &m, W, Werr,
   *  &wl, &wu, iblock, Windex, work, iwork, &info);
   * assert(info == 0);
   * assert(m == ilast-ifirst+1); */

  nthreads = max_nthreads;
  while (nthreads > 1 && isize / nthreads < 2)
    nthreads--;

  if (nthreads > 1) {

    /* each threads computes W[iifirst:iilast] and places them in
     * work[0:n-1]; the corresponding errors in work[n:2*n-1];
     * the blocks they belong in iwork[0:n-1]; and their indices in
     * iwork[n:2*n-1]; */ 

    iifirst = ifirst;
    chunk = isize / nthreads;
    for (i=1; i<nthreads; i++) {

      iilast = iifirst + chunk - 1; 
      
      auxarg1 = create_auxarg1(n, D, E, E2, ifirst, ilast, iifirst, iilast,
			       nsplit, isplit, bsrtol, pivmin, gersch,
			       &work[0], &work[n], &iwork[n], &iwork[0]);

      info = pthread_create(&threads[i], &attr,
			    eigval_subset_thread_a,
			    (void *) auxarg1);
      assert(info == 0);

      iifirst = iilast + 1;
    }
    iilast = ilast;

    auxarg1 = create_auxarg1(n, D, E, E2, ifirst, ilast, iifirst, iilast,
			     nsplit, isplit, bsrtol, pivmin, gersch,
			     &work[0], &work[n], &iwork[n], &iwork[0]);

    status = eigval_subset_thread_a( (void *) auxarg1 );
    assert(status == NULL);

    /* join threads */
    for (i=1; i<nthreads; i++) {
      info = pthread_join(threads[i], &status);
      assert(info == 0 && status == NULL);
    }

    /* sort, m counts the numbers of eigenvalues computed by process */
    m = 0;
    for (i=1; i<=nsplit; i++) {
      for (j=0; j<isize; j++) {
	if (iwork[j] == i) {
	  W[m]      = work[j];
	  Werr[m]   = work[j+n];
	  iblock[m] = iwork[j];
	  Windex[m] = iwork[j+n];
	  m++;
	}
      }
    }

  } else {
    /* no multithreaded computation */
    LAPACK(dlarrd)
    ("I", "B", &n, &dummy, &dummy, &ifirst, &ilast, gersch, &bsrtol, D, E, E2,
     &pivmin, &nsplit, isplit, &m, W, Werr, &wl, &wu, iblock, Windex, work, 
     iwork, &info);
    assert(info == 0);
    assert(m == ilast-ifirst+1);
  }


  /* loop over unreduced blocks */  
  num_vals  = m;
  m         = 0; /* accumulates eigenvalues found or refined */
  bl_begin  = 0;
  bl_Wbegin = 0;
  
  for (jbl=0; jbl<nsplit; jbl++) {
    
    bl_end  = isplit[jbl] - 1;
    bl_size = bl_end - bl_begin + 1;

    /* deal with 1x1 block immediately */
    if (bl_size == 1) {
      E[bl_end] = 0.0;
      /* if eigenvalue part of process work */
      if (iblock[bl_Wbegin] == jbl+1) {
	W[m]      = D[bl_begin];
	Werr[m]   = 0.0;
	Werr[m]   = 0.0;
	iblock[m] = jbl + 1;
	Windex[m] = 1;
	m++;
	bl_Wbegin++;
      }
      bl_begin  = bl_end + 1;
      continue;
    }


    /* COMPUTE ROOT RRR */

    /* store shift of initial RRR, here set to zero */
    E[bl_end] = 0.0;

    /* find outer bounds GL, GU for block and spectral diameter */
    gl = D[bl_begin];
    gu = D[bl_begin];
    for (i = bl_begin; i <= bl_end; i++) {
      gl = fmin(gl, gersch[2*i]  );
      gu = fmax(gu, gersch[2*i+1]);
    }
    spdiam = gu - gl;

    /* find approximation of extremal eigenvalues of the block
     * dlarrk computes one eigenvalue of tridiagonal matrix T
     * tmp1 and tmp2 one hold the eigenvalue and error, respectively */
    LAPACK(dlarrk)
    (&bl_size, &IONE, &gl, &gu, &D[bl_begin], &E2[bl_begin], &pivmin, &rtl, 
     &tmp1, &tmp2, &info);
    assert(info == 0);  /* if info=-1 => eigenvalue did not converge */
    
    isleft = fmax(gl, tmp1-tmp2 - HUNDRED*DBL_EPSILON*fabs(tmp1-tmp2) );
    
    LAPACK(dlarrk)
    (&bl_size, &bl_size, &gl, &gu, &D[bl_begin], &E2[bl_begin], &pivmin, &rtl,     &tmp1, &tmp2, &info);
    assert(info == 0);  /* if info=-1 => eigenvalue did not converge */
    
    isright = fmin(gu, tmp1+tmp2 + HUNDRED*DBL_EPSILON*fabs(tmp1+tmp2) );
    
    spdiam = isright - isleft;

    /* compute negcount at points s1 and s2 */
    s1 = isleft  + HALF   * spdiam;
    s2 = isright - FOURTH * spdiam;  /* not needed currently */

    /* compute negcount at points s1 and s2 */
    /* cnt = number of eigenvalues in (s1,s2] = count_right - count_left
     * negcnt_lft = number of eigenvalues smaller equals than s1
     * negcnt_rgt = number of eigenvalues smaller equals than s2 */
    LAPACK(dlarrc)
    ("T", &bl_size, &s1, &s2, &D[bl_begin], &E[bl_begin], &pivmin, &cnt, 
     &negcnt_lft, &negcnt_rgt, &info);
    assert(info == 0);

    /* if more of the desired eigenvectors are in the left part shift left
     * and the other way around */
    if ( negcnt_lft >= bl_size - negcnt_lft ) {
      /* shift left */
      sigma = isleft;
      sgndef = ONE;
    } else {
      /* shift right */
      sigma = isright;
      sgndef = -ONE;
    }

    /* define increment to perturb initial shift to find RRR
     * with not too much element growth */
    tau = spdiam*DBL_EPSILON*n + 2.0*pivmin;


    /* try to find initial RRR of block:
     * need work space of 3*n here to store D, L, D^-1 of possible
     * representation:
     * D_try      = work[0  :  n-1] 
     * L_try      = work[n  :2*n-1]
     * inv(D_try) = work[2*n:3*n-1] */

    off_L    = n;
    off_invD = 2*n;
    
    for (jtry = 0; jtry < MAX_TRY_RRRR; jtry++) {

      Dpivot  = D[bl_begin] - sigma;
      work[0] = Dpivot;
      Dmax    = fabs( work[0] );
      j = bl_begin;

      for (i = 0; i < bl_size-1; i++) {
 	work[i+off_invD] = 1.0 / work[i];
	tmp = E[j] * work[i+off_invD];
	work[i+off_L] = tmp;
	Dpivot = (D[j+1] - sigma) - tmp*E[j];
	work[i+1] = Dpivot;
	Dmax = fmax(Dmax, fabs(Dpivot) );
	j++;
      }
      
      /* except representation only if not too much element growth */
      if (Dmax > MAX_GROWTH*spdiam) {
	noREP = true;
      } else {
	noREP = false;
      }
      
      if (noREP == true) {
	/* if all eigenvalues are desired shift is made definite to use DQDS
	 * so we should not end here */
	if (jtry == MAX_TRY_RRRR-2) {
	  if (sgndef == ONE) { /* floating point comparison okay here */
	    sigma = gl - FUDGE_FACTOR*spdiam*DBL_EPSILON*n 
	               - FUDGE_FACTOR*2.0*pivmin;
	  } else {
	    sigma = gu + FUDGE_FACTOR*spdiam*DBL_EPSILON*n 
                       + FUDGE_FACTOR*2.0*pivmin;
	  }
	} else if (jtry == MAX_TRY_RRRR-1) {
	  fprintf(stderr,"No initial representation could be found.\n");
	  exit(3);
	} else {
	  sigma -= sgndef*tau;
	  tau   *= 2.0;
	  continue;
	}
      } else {   /* found representation */
	break;
      }  
    }
    /* end trying to find initial RRR of block */


    /* save initial RRR and corresponding shift */
    E[bl_end] = sigma;
    memcpy(&D[bl_begin], &work[0],  bl_size    * sizeof(double) );
    memcpy(&E[bl_begin], &work[n], (bl_size-1) * sizeof(double) );
    /* work[0:4*n-1] can now be used again for anything */


    /* perturb root rrr by small relative amount, first make sure
     * that at least two values are actually disturbed enough,
     * which might not be necessary */
    while( fabs(randvec[bl_begin])*RAND_FACTOR < 1.0 )
      randvec[bl_begin] *= 2.0;
    while( fabs(randvec[bl_end])  *RAND_FACTOR < 1.0 )
      randvec[bl_end]   *= 2.0;

    for (i=bl_begin; i<bl_end; i++) {
      D[i] *= 1.0 + DBL_EPSILON*RAND_FACTOR*randvec[i];
      E[i] *= 1.0 + DBL_EPSILON*RAND_FACTOR*randvec[i+n];
    }
    D[bl_end] *= 1.0 + DBL_EPSILON*RAND_FACTOR*randvec[bl_end];


    /* REFINE EIGENVALUES WITH REPECT TO RRR */

    /* count number of eigenvalues in block and find smallest
     * and largest index of block */
    bl_m  = 0;
    i_low = n;
    i_upp = 1;
    for (i=bl_Wbegin; i<num_vals; i++) {
      if (iblock[i] == jbl+1)  {
	bl_m++;
	i_low = imin(i_low, Windex[i]);
	i_upp = imax(i_upp, Windex[i]);
      } else {
	break;
      }
    }

    if (bl_m == 0) {
      bl_begin  = bl_end + 1;
      continue; /* go to next block */
    }

    /* last index of W to store eigenvalues of block */
    bl_Wend = bl_Wbegin + bl_m - 1;

    /* calculate gaps */
    for (i=bl_Wbegin; i<bl_Wend; i++) {
      Wgap[i] = fmax(0.0, (W[i+1] - Werr[i+1]) - (W[i] + Werr[i]) );
    }	
    
    Wgap[bl_Wend] = fmax(0.0, gu - (W[bl_Wend] + Werr[bl_Wend]) );
    
    /* shift eigenvalues to be consistent with dqds 
     * and compute eigenvalues of SHIFTED matrix */
    for (i=bl_Wbegin; i<=bl_Wend; i++) {
      W[i]    -= sigma;
      Werr[i] += fabs(W[i])*DBL_EPSILON;
    }

    /* work  for sequential dlarrb = work[0:2*n-1]
     * iwork for sequential dlarrb = iwork[0:2*n-1]
     * DE2 = work[2*n:3*n-1] strting at bl_begin */
    off_DE2 = 2*n;
    
    /* compute DE2 at store it in work[bl_begin+2*n:bl_end-1+2*n] */
    for (i=bl_begin; i<bl_end; i++) {
      work[i+off_DE2] = D[i]*E[i]*E[i];
    }
    
    nthreads = max_nthreads;
    while (nthreads > 1 && bl_m/nthreads < 2) {
      nthreads--;
    }

    if (nthreads > 1) {

      rf_begin = bl_Wbegin;
      chunk    = bl_m / nthreads;
      for (i=1; i<nthreads; i++) {
	
	rf_end = rf_begin + chunk - 1; 

	auxarg2 = create_auxarg2(bl_size, &D[bl_begin],
				 &work[bl_begin+off_DE2],
				 rf_begin, rf_end, Wstruct,
				 tolstruct->rtol1, tolstruct->rtol2,
				 pivmin, spdiam);
	
	info = pthread_create(&threads[i], &attr,
			      eigval_subset_thread_r,
			      (void *) auxarg2);
	assert(info == 0);
	
	rf_begin = rf_end + 1;
      }
      rf_end = bl_Wend;

      auxarg2 = create_auxarg2(bl_size, &D[bl_begin],
			       &work[bl_begin+off_DE2],
			       rf_begin, rf_end, Wstruct,
			       tolstruct->rtol1, tolstruct->rtol2,
			       pivmin, spdiam);
      
      status = eigval_subset_thread_r( (void *) auxarg2 );
      assert(status == NULL);
    
      /* join threads */
      for (i=1; i<nthreads; i++) {
	info = pthread_join(threads[i], &status);
	assert(info == 0 && status == NULL);
      }
      /* should update gaps at splitting points here, but the gaps
       * will be recomputed anyway */
      
    } else {

      offset = i_low-1;
      
      /* refine eigenvalues found by dlarrd for i_low:i_upp */
      LAPACK(dlarrb)
      (&bl_size, &D[bl_begin], &work[bl_begin+off_DE2], &i_low, &i_upp, 
       &tolstruct->rtol1, &tolstruct->rtol2, &offset, &W[bl_Wbegin], 
       &Wgap[bl_Wbegin], &Werr[bl_Wbegin], work, iwork, &pivmin, &spdiam, 
       &bl_size, &info);
      assert(info == 0);
      /* needs work of dim(2*n) and iwork of dim(2*n) */
    }
    /* dlarrb computes gaps correctly, but not last one;
     * this is ignored since the gaps are recomputed anyway */
    
    /* this makes sure that the indices are in the right order */
    for (i=i_low; i<=i_upp; i++) {
      Windex[m] = i;
      m++;
    }

    /* proceed with next block */
    bl_begin  = bl_end  + 1;
    bl_Wbegin = bl_Wend + 1;
  }
  /* end of loop over unreduced blocks */  
  
  /* clean up */
  free(randvec);
  free(threads);

  if (max_nthreads > 1) {
    pthread_attr_destroy(&attr);
  }

  return(0);
}




static 
void *eigval_subset_thread_a(void *argin)
{
  /* from input argument */
  int    n, il, iu, my_il, my_iu;
  double *D, *E, *E2, *gersch;
  double bsrtol, pivmin;
  int    nsplit, *isplit;

  /* others */
  int    info;
  double dummy1, dummy2;
  int    num_vals;
  double *W_tmp, *Werr_tmp, *W, *Werr;
  int    *iblock_tmp, *Windex_tmp, *iblock, *Windex;
  double *work;
  int    *iwork;
  
  retrieve_auxarg1((auxarg1_t *) argin, &n, &D, &E, &E2,
		   &il, &iu, &my_il, &my_iu, &nsplit,
		   &isplit, &bsrtol, &pivmin, &gersch,
		   &W, &Werr, &Windex, &iblock);

  /* allocate memory needed for dlarrd */
  W_tmp = (double *) malloc( n * sizeof(double) );
  assert(W_tmp != NULL);
  
  Werr_tmp = (double *) malloc( n * sizeof(double) );
  assert(Werr_tmp != NULL);
  
  Windex_tmp = (int *) malloc( n * sizeof(int) );
  assert(Windex_tmp != NULL);

  iblock_tmp = (int *) malloc( n * sizeof(int) );
  assert(iblock_tmp != NULL);

  work  = (double *) malloc( 4*n * sizeof(double) );
  assert (work != NULL);

  iwork = (int *) malloc( 3*n * sizeof(int) );
  assert (iwork != NULL);

  /* compute eigenvalues 'my_il' to 'my_iu', put into temporary arrays */
  LAPACK(dlarrd)
  ("I", "B", &n, &dummy1, &dummy2, &my_il, &my_iu, gersch, &bsrtol, D, E, E2, 
   &pivmin, &nsplit, isplit, &num_vals, W_tmp, Werr_tmp, &dummy1, &dummy2, 
   iblock_tmp, Windex_tmp, work, iwork, &info);
  assert(info == 0);

  /* copy computed values in W, Werr, Windex, iblock (which are work space) */
  memcpy(&W[my_il-il],      W_tmp,      num_vals * sizeof(double) );
  memcpy(&Werr[my_il-il],   Werr_tmp,   num_vals * sizeof(double) );
  memcpy(&Windex[my_il-il], Windex_tmp, num_vals * sizeof(int)    );
  memcpy(&iblock[my_il-il], iblock_tmp, num_vals * sizeof(int)    );
  
  free(W_tmp);
  free(Werr_tmp);
  free(Windex_tmp);
  free(iblock_tmp);
  free(work);
  free(iwork);

  return(NULL);
}




static 
auxarg1_t *create_auxarg1(int n, double *D, double *E, double *E2,
			  int il, int iu, int my_il, int my_iu, 
			  int nsplit, int *isplit, double bsrtol, 
			  double pivmin, double *gersch, double *W, 
			  double *Werr, int *Windex, int *iblock)
{
  auxarg1_t *arg;

  arg = (auxarg1_t *) malloc( sizeof(auxarg1_t) );
  assert(arg != NULL);

  arg->n       = n;
  arg->D       = D;
  arg->E       = E;
  arg->E2      = E2;
  arg->il      = il;
  arg->iu      = iu;
  arg->my_il   = my_il;
  arg->my_iu   = my_iu;
  arg->nsplit  = nsplit;
  arg->isplit  = isplit;
  arg->bsrtol  = bsrtol;
  arg->pivmin  = pivmin;
  arg->gersch  = gersch;
  arg->W       = W;
  arg->Werr    = Werr;
  arg->Windex  = Windex;
  arg->iblock  = iblock;

  return(arg);
}




static 
void retrieve_auxarg1(auxarg1_t *arg, int *n, double **D, double **E,
		      double **E2, int *il, int *iu, int *my_il, 
		      int *my_iu, int *nsplit, int **isplit, 
		      double *bsrtol, double *pivmin, double **gersch, 
		      double **W, double **Werr, int **Windex, 
		      int **iblock)
{
  *n      = arg->n;
  *D      = arg->D;
  *E      = arg->E;
  *E2     = arg->E2;
  *il     = arg->il;
  *iu     = arg->iu;
  *my_il  = arg->my_il;
  *my_iu  = arg->my_iu;
  *nsplit = arg->nsplit;
  *isplit = arg->isplit;
  *bsrtol = arg->bsrtol;
  *pivmin = arg->pivmin;
  *gersch = arg->gersch;
  *W      = arg->W;
  *Werr   = arg->Werr;
  *Windex = arg->Windex;
  *iblock = arg->iblock;

  free(arg);
}




static 
void *eigval_subset_thread_r(void *argin)
{
  /* from input argument */
  int          bl_size, rf_begin, rf_end;
  double       *D, *DE2;
  double       rtol1, rtol2, pivmin;
  double       bl_spdiam;
  val_t        *Wstruct;

  /* others */
  int          info, offset;
  double       *W, *Werr, *Wgap;
  int          *Windex;
  double       *work;
  int          *iwork;

  retrieve_auxarg2((auxarg2_t *) argin, &bl_size, &D, &DE2,
		   &rf_begin, &rf_end, &Wstruct, &rtol1, &rtol2,
		   &pivmin, &bl_spdiam);

  /* malloc work space */
  work = (double *) malloc( 2*bl_size * sizeof(double) );
  assert(work != NULL);
  
  iwork = (int *)   malloc( 2*bl_size * sizeof(int) );
  assert(iwork != NULL);

  W      = Wstruct->W;
  Werr   = Wstruct->Werr;
  Wgap   = Wstruct->Wgap;
  Windex = Wstruct->Windex;

  /* special case of only one eigenvalue */
  if (rf_begin == rf_end)
    Wgap[rf_begin] = 0.0;
 
  offset = Windex[rf_begin] - 1;

  /* call bisection routine to refine the eigenvalues */
  LAPACK(dlarrb)
  (&bl_size, D, DE2, &Windex[rf_begin], &Windex[rf_end], &rtol1, &rtol2, 
   &offset, &W[rf_begin], &Wgap[rf_begin], &Werr[rf_begin], work, iwork, 
   &pivmin, &bl_spdiam, &bl_size, &info);
  assert(info == 0);

  /* clean up */
  free(work);
  free(iwork);

  return(NULL);
}




static 
auxarg2_t *create_auxarg2(int bl_size, double *D, double *DE2,
			  int rf_begin, int rf_end, val_t *Wstruct,
			  double rtol1, double rtol2, double pivmin, 
			  double bl_spdiam)
{
  auxarg2_t *arg;

  arg = (auxarg2_t *) malloc( sizeof(auxarg2_t) );
  assert(arg != NULL);

  arg->bl_size   = bl_size;
  arg->D         = D;
  arg->DE2       = DE2;
  arg->rf_begin  = rf_begin;
  arg->rf_end    = rf_end;
  arg->Wstruct   = Wstruct;
  arg->rtol1     = rtol1;
  arg->rtol2     = rtol2;
  arg->pivmin    = pivmin;
  arg->bl_spdiam = bl_spdiam;

  return(arg);
}




static 
void retrieve_auxarg2(auxarg2_t *arg, int *bl_size, double **D,
		      double **DE2, int *rf_begin, int *rf_end,
		      val_t **Wstruct, double *rtol1, double *rtol2,
		      double *pivmin, double *bl_spdiam)
{
  *bl_size   = arg->bl_size;
  *D         = arg->D;
  *DE2       = arg->DE2;
  *rf_begin  = arg->rf_begin;
  *rf_end    = arg->rf_end;
  *Wstruct   = arg->Wstruct;
  *rtol1     = arg->rtol1;
  *rtol2     = arg->rtol2;
  *pivmin    = arg->pivmin;
  *bl_spdiam = arg->bl_spdiam;

  free(arg);
}
