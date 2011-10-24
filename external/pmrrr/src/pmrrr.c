/* Computation of eigenvalues and eigenvectors of a symmetric
 * tridiagonal matrix T, given by its diagonal elements D
 * and its super-/subdiagonal elements E.
 *
 * See pmrrr.h for more information.
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
#include <math.h>
#include <string.h>
#include <assert.h>
#include <float.h>

#include "mpi.h"
#include "pmrrr.h"
#include "global.h"
#include "plarre.h"
#include "plarrv.h"
#include "structs.h"


static int handle_small_cases(char*, char*, int*, double*, double*,
			      double*, double*, int*, int*, int*,
			      MPI_Comm, int*, int*, double*, double*,
			      int*, int*);
static int cmp(const void*, const void*);
static int cmpb(const void*, const void*);
static double scale_matrix(in_t*, val_t*, bool);
static void invscale_eigenvalues(val_t*, double, int);
static int sort_eigenpairs(proc_t*, val_t*, vec_t*);
static void clean_up(MPI_Comm, double*, double*, double*, 
		     int*, int*, int*, int*, int*, proc_t*, 
		     in_t*, val_t*, vec_t*, tol_t*);
static int refine_to_highrac(proc_t*, char*, double*, double*,
			     in_t*, int*, val_t*, tol_t*);




/* 
 * Computation of eigenvalues and eigenvectors of a symmetric
 * tridiagonal matrix T, given by its diagonal elements D
 * and its super-/subdiagonal elements E.
 * See README or 'pmrrr.h' for details.
 */

int PMRRR(char *jobz, char *range, int *np, double  *D,
	  double *E, double *vl, double *vu, int *il,
	  int *iu, int *tryracp, MPI_Comm comm, int *nzp,
	  int *offsetp, double *W, double *Z, int *ldz,
	  int *Zsupp)
{
  /* Input parameter */
  int         n      = *np;
  bool        onlyW  = (jobz[0]  == 'N' || jobz[0]  == 'n');
  bool        wantZ  = (jobz[0]  == 'V' || jobz[0]  == 'v');
  bool        cntval = (jobz[0]  == 'C' || jobz[0]  == 'c');
  bool        alleig = (range[0] == 'A' || range[0] == 'a');
  bool        valeig = (range[0] == 'V' || range[0] == 'v');
  bool        indeig = (range[0] == 'I' || range[0] == 'i');

  /* Work space */
  double      *Werr,   *Wgap,  *gersch, *Dcopy, *E2copy;
  int         *iblock, *iproc, *Windex, *Zindex, *isplit;

  /* Structures to store variables */
  proc_t      *procinfo;
  in_t        *Dstruct;  
  val_t       *Wstruct;
  vec_t       *Zstruct;
  tol_t       *tolstruct;

  /* Multiprocessing and multithreading */
  int         nproc, pid, nthreads;             
  char        *ompvar;
  MPI_Comm    comm_dup;
  int         thread_support;

  /* Others */
  double      scale;              
  int         i, info;

  /* Check input parameters */
  if(!(onlyW  || wantZ  || cntval)) return(1);
  if(!(alleig || valeig || indeig)) return(1);
  if(n <= 0) return(1);
  if (valeig) {
    if(*vu<=*vl) return(1);
  } else if (indeig) {
    if (*il<1 || *il>n || *iu<*il || *iu>n) return(1);
  }
  
  /* MPI & multithreading info */
  MPI_Comm_dup(comm, &comm_dup);
  MPI_Comm_size(comm_dup, &nproc);
  MPI_Comm_rank(comm_dup, &pid);
  MPI_Query_thread(&thread_support);

  if ( !(thread_support == MPI_THREAD_MULTIPLE ||
         thread_support == MPI_THREAD_FUNNELED) ) {
    /* Disable multithreading; note: to support multithreading with 
     * MPI_THREAD_SERIALIZED the code must be changed slightly; this 
     * is not supported at the moment */
    nthreads = 1;
  } else {
    ompvar = getenv("PMR_NUM_THREADS");
    if (ompvar == NULL) {
      nthreads = DEFAULT_NUM_THREADS;
    } else {
      nthreads = atoi(ompvar);
    }
  }

  /* If only maximal number of local eigenvectors are queried
   * return if possible here */
  *nzp     = 0;
  *offsetp = 0;
  if (cntval) {
    if ( alleig || n < DSTEMR_IF_SMALLER ) {
      *nzp = iceil(n,nproc);
      MPI_Comm_free(&comm_dup);
      return(0);
    } else if (indeig) {
      *nzp = iceil(*iu-*il+1,nproc);
      MPI_Comm_free(&comm_dup);
      return(0);
    }
  }

  /* Check if computation should be done by multiple processes */
  if (n < DSTEMR_IF_SMALLER) {
    info = handle_small_cases(jobz, range, np, D, E, vl, vu, il,
			      iu, tryracp, comm, nzp, offsetp, W,
			      Z, ldz, Zsupp);
    MPI_Comm_free(&comm_dup);
    return(info);
  }

  /* Allocate memory */
  Werr    = (double *)   malloc( n * sizeof(double) );
  assert(Werr != NULL);
  Wgap      = (double *) malloc( n * sizeof(double) );
  assert(Wgap != NULL);
  gersch    = (double *) malloc( 2*n*sizeof(double) );
  assert(gersch != NULL);
  iblock    = (int *)    calloc( n , sizeof(int) );
  assert(iblock != NULL);
  iproc     = (int *)    malloc( n * sizeof(int) );
  assert(iproc != NULL);
  Windex    = (int *)    malloc( n * sizeof(int) );
  assert(Windex != NULL);
  isplit    = (int *)    malloc( n * sizeof(int) );
  assert(isplit != NULL);
  Zindex    = (int *)    malloc( n * sizeof(int) );
  assert(Zindex != NULL);
  procinfo  = (proc_t *) malloc( sizeof(proc_t) );
  assert(procinfo != NULL);
  Dstruct   = (in_t *)   malloc( sizeof(in_t) );
  assert(Dstruct != NULL);
  Wstruct   = (val_t *)  malloc( sizeof(val_t) );
  assert(Wstruct != NULL);
  Zstruct   = (vec_t *)  malloc( sizeof(vec_t) );
  assert(Zstruct != NULL);
  tolstruct = (tol_t *)  malloc( sizeof(tol_t) );
  assert(tolstruct != NULL);

  /* Bundle variables into a structures */
  procinfo->pid            = pid;
  procinfo->nproc          = nproc;
  procinfo->comm           = comm_dup;
  procinfo->nthreads       = nthreads;
  procinfo->thread_support = thread_support;

  Dstruct->n               = n;
  Dstruct->D               = D;
  Dstruct->E               = E;
  Dstruct->isplit          = isplit;

  Wstruct->n               = n;
  Wstruct->vl              = vl;
  Wstruct->vu              = vu;
  Wstruct->il              = il;
  Wstruct->iu              = iu;
  Wstruct->W               = W;
  Wstruct->Werr            = Werr;
  Wstruct->Wgap            = Wgap;
  Wstruct->Windex          = Windex;
  Wstruct->iblock          = iblock;
  Wstruct->iproc           = iproc;
  Wstruct->gersch          = gersch;

  Zstruct->ldz             = *ldz;
  Zstruct->nz              = 0;
  Zstruct->Z               = Z;
  Zstruct->Zsupp           = Zsupp;
  Zstruct->Zindex          = Zindex;

  /* Scale matrix to allowable range, returns 1.0 if not scaled */
  scale = scale_matrix(Dstruct, Wstruct, valeig);

  /*  Test if matrix warrants more expensive computations which
   *  guarantees high relative accuracy */
  if (*tryracp) {
    LAPACK(dlarrr)(&n, D, E, &info); /* 0 - rel acc */
  }
  else info = -1;

  if (info == 0) {
    /* This case is extremely rare in practice */ 
    tolstruct->split = DBL_EPSILON;
    /* Copy original data needed for refinement later */
    Dcopy  = (double *) malloc( n * sizeof(double) );
    assert(Dcopy != NULL);
    memcpy(Dcopy, D, n*sizeof(double));  
    E2copy = (double *) malloc( n * sizeof(double) );
    assert(E2copy != NULL);
    for (i=0; i<n-1; i++) E2copy[i] = E[i]*E[i];
  } else {
    /* Neg. threshold forces old splitting criterion */
    tolstruct->split = -DBL_EPSILON; 
    *tryracp = 0;
  }

  if (!wantZ) {
    /* Compute eigenvalues to full precision */
    tolstruct->rtol1 = 4.0 * DBL_EPSILON;
    tolstruct->rtol2 = 4.0 * DBL_EPSILON;
  } else {
    /* Do not compute to full accuracy first, but refine later */
    tolstruct->rtol1 = sqrt(DBL_EPSILON);
    tolstruct->rtol1 = fmin(1e-1*MIN_RELGAP, tolstruct->rtol1);
    tolstruct->rtol2 = sqrt(DBL_EPSILON)*5.0E-3;
    tolstruct->rtol2 = fmin(5e-6*MIN_RELGAP, tolstruct->rtol2);
    tolstruct->rtol2 = fmax(4.0 * DBL_EPSILON, tolstruct->rtol2);
  }

  /*  Compute eigenvalues (all in case vectors are desired too) */
  info = plarre(procinfo, jobz, range, Dstruct, Wstruct, tolstruct,
		nzp, offsetp);
  assert(info == 0);

  /* If just number of local eigenvectors are queried return, then 
   * value is set in 'plarre' and we can return */
  if (cntval & valeig) {
    clean_up(comm_dup, Werr, Wgap, gersch, iblock, iproc, Windex,
	     isplit, Zindex, procinfo, Dstruct, Wstruct, Zstruct,
	     tolstruct);
    return(0);
  }

  /* If only eigenvalues are to be computed, clean up and return */
  if (!wantZ) {

    /* Refine to high relative with respect to input T */
    if (*tryracp) {
      info = refine_to_highrac(procinfo, jobz, Dcopy, E2copy, 
			       Dstruct, nzp, Wstruct, tolstruct);
      assert(info == 0);
    }

    /* If matrix was scaled, rescale eigenvalues */
    invscale_eigenvalues(Wstruct, scale, *nzp);

    /* Sort eigenvalues */
    if (Dstruct->nsplit > 1) qsort(W, *nzp, sizeof(double), cmp);

    clean_up(comm_dup, Werr, Wgap, gersch, iblock, iproc, Windex,
	     isplit, Zindex, procinfo, Dstruct, Wstruct, Zstruct,
	     tolstruct);

    return(0);
  } /* end of only eigenvalues to compute */

  /* Compute eigenvectors */
  info = plarrv(procinfo, Dstruct, Wstruct, Zstruct, tolstruct, 
		nzp, offsetp);
  assert(info == 0);

  /* Refine to high relative with respect to input matrix */
  if (*tryracp) {
    info = refine_to_highrac(procinfo, jobz, Dcopy, E2copy, 
			     Dstruct, nzp, Wstruct, tolstruct);
    assert(info == 0);
  }

  /* If matrix was scaled, rescale eigenvalues */
  invscale_eigenvalues(Wstruct, scale, n);

  /* Sort eigenvalues and eigenvectors of process */
  sort_eigenpairs(procinfo, Wstruct, Zstruct);

  clean_up(comm_dup, Werr, Wgap, gersch, iblock, iproc, Windex,
	   isplit, Zindex, procinfo, Dstruct, Wstruct, Zstruct,
	   tolstruct);
  if (*tryracp) {
    free(Dcopy);
    free(E2copy);
  }

  return(0);
} /* end pmrrr */




/*
 * Free's on allocated memory of pmrrr routine
 */
static  
void clean_up(MPI_Comm comm, double *Werr, double *Wgap,
	      double *gersch, int *iblock, int *iproc,
	      int *Windex, int *isplit, int *Zindex,
	      proc_t *procinfo, in_t *Dstruct,
	      val_t *Wstruct, vec_t *Zstruct,
	      tol_t *tolstruct)
{
  MPI_Comm_free(&comm);
  free(Werr);
  free(Wgap);
  free(gersch);
  free(iblock);
  free(iproc);
  free(Windex);
  free(isplit);
  free(Zindex);
  free(procinfo);
  free(Dstruct);
  free(Wstruct);
  free(Zstruct);
  free(tolstruct);
}




/*
 * Wrapper to call LAPACKs DSTEMR for small matrices
 */
static
int handle_small_cases(char *jobz, char *range, int *np, double  *D,
		       double *E, double *vlp, double *vup, int *ilp,
		       int *iup, int *tryracp, MPI_Comm comm, int *nzp,
		       int *myfirstp, double *W, double *Z, int *ldzp,
		       int *Zsupp)
{
  bool   cntval  = (jobz[0]  == 'C' || jobz[0]  == 'c');
  bool   onlyW   = (jobz[0]  == 'N' || jobz[0]  == 'n');
  bool   wantZ   = (jobz[0]  == 'V' || jobz[0]  == 'v');
  bool   indeig  = (range[0] == 'I' || range[0] == 'i');
  int    n       = *np;
  int    ldz_tmp = *np;
  int    ldz     = *ldzp;

  int    nproc, pid;
  int    m, lwork, *iwork, liwork, info;
  double *Z_tmp, *work, cnt;
  int    i, itmp, MINUSONE=-1;
  int    chunk, myfirst, mylast, mysize;

  MPI_Comm_size(comm, &nproc);
  MPI_Comm_rank(comm, &pid);
  
  if (onlyW) {
    lwork  = 12*n;
    liwork =  8*n;
  } else if (cntval) {
    lwork  = 18*n;
    liwork = 10*n;
  } else if (wantZ) {
    lwork  = 18*n;
    liwork = 10*n;
    if (indeig) itmp = *iup-*ilp+1;
    else        itmp = n;
    Z_tmp = (double *) malloc(n*itmp * sizeof(double));
    assert(Z_tmp != NULL);
  } else {
    return(1);
  }

  work = (double *) malloc( lwork  * sizeof(double));
  assert(work != NULL);
  iwork = (int *)   malloc( liwork * sizeof(int));
  assert(iwork != NULL);

  if (cntval) {
    /* Note: at the moment, jobz="C" should never get here, since
     * it is blocked before. */
    LAPACK(dstemr)
    ("V", "V", np, D, E, vlp, vup, ilp, iup, &m, W, &cnt, &ldz_tmp, &MINUSONE, 
     Zsupp, tryracp, work, &lwork, iwork, &liwork, &info);
    assert(info == 0);
    
    *nzp = (int) ceil(cnt/nproc);
    free(work); free(iwork);
    return(0);
  }

  LAPACK(dstemr)
  (jobz, range, np, D, E, vlp, vup, ilp, iup, &m, W, Z_tmp, &ldz_tmp, np, Zsupp,
   tryracp, work, &lwork, iwork, &liwork, &info);
  assert(info == 0);

  chunk   = iceil(m,nproc);
  myfirst = imin(pid * chunk, m);
  mylast  = imin((pid+1)*chunk - 1, m - 1);
  mysize  = mylast - myfirst + 1;

  if (mysize > 0) {
    memmove(W, &W[myfirst], mysize*sizeof(double));
    if (wantZ) {
      if (ldz == ldz_tmp) {
	/* copy everything in one chunk */
	memcpy(Z, &Z_tmp[myfirst*ldz_tmp], n*mysize*sizeof(double));
      } else {
	/* copy each vector seperately */
	for (i=0; i<mysize; i++)
	  memcpy(&Z[i*ldz], &Z_tmp[(myfirst+i)*ldz_tmp], 
		 n*sizeof(double));
      } 
    } /* if (wantZ) */
  } 
  
  *myfirstp = myfirst;
  *nzp      = mysize;

  if (wantZ) free(Z_tmp);
  free(work);
  free(iwork);

  return(0);
}




/*
 * Scale matrix to allowable range, returns 1.0 if not scaled
 */
static 
double scale_matrix(in_t *Dstruct, val_t *Wstruct, bool valeig)
{
  int              n  = Dstruct->n;
  double *restrict D  = Dstruct->D;
  double *restrict E  = Dstruct->E;
  double          *vl = Wstruct->vl;
  double          *vu = Wstruct->vu;

  double           scale = 1.0;
  double           T_norm;              
  double           smlnum, bignum, rmin, rmax;
  int              IONE = 1, itmp;

  /* Set some machine dependent constants */
  smlnum = DBL_MIN / DBL_EPSILON;
  bignum = 1.0 / smlnum;
  rmin   = sqrt(smlnum);
  rmax   = fmin(sqrt(bignum), 1.0 / sqrt(sqrt(DBL_MIN)));

  /*  Scale matrix to allowable range */
  T_norm = LAPACK(dlanst)("M", &n, D, E);  /* returns max(|T(i,j)|) */
  if (T_norm > 0 && T_norm < rmin) {
    scale = rmin / T_norm;
  } else if (T_norm > rmax) {
    scale = rmax / T_norm;
  }

  if (scale != 1.0) {  /* FP cmp okay */
    /* Scale matrix and matrix norm */
    itmp = n-1;
    BLAS(dscal)(&n,    &scale, D, &IONE);
    BLAS(dscal)(&itmp, &scale, E, &IONE);
    if (valeig == true) {
      /* Scale eigenvalue bounds */
      *vl *= scale;
      *vu *= scale;
    }
  } /* end scaling */

  return(scale);
}




/*
 * If matrix scaled, rescale eigenvalues
 */
static 
void invscale_eigenvalues(val_t *Wstruct, double scale,
			  int size)
{
  double *vl = Wstruct->vl;
  double *vu = Wstruct->vu;
  double *W  = Wstruct->W;
  double invscale = 1.0 / scale;
  int    IONE = 1;

  if (scale != 1.0) {  /* FP cmp okay */
    *vl *= invscale;
    *vu *= invscale;
    BLAS(dscal)(&size, &invscale, W, &IONE);
  }
}




static 
int sort_eigenpairs_local(int m, val_t *Wstruct, vec_t *Zstruct)
{
  int              n        = Wstruct->n;
  double *restrict W        = Wstruct->W;
  double *restrict work     = Wstruct->gersch;
  int              ldz      = Zstruct->ldz;
  double *restrict Z        = Zstruct->Z;
  int    *restrict Zsupp    = Zstruct->Zsupp;
 
  bool             sorted;
  int              j;
  double           tmp;
  int              itmp1, itmp2;

  /* Make sure that sorted correctly; ineffective implementation,
   * but usually no or very little swapping should be done here */
  sorted = false;
  while (sorted == false) {
    sorted = true;
    for (j=0; j<m-1; j++) {
      if (W[j] > W[j+1]) {
	sorted = false;
	/* swap eigenvalue */
	tmp    = W[j];
	W[j]   = W[j+1];
	W[j+1] = tmp;
	/* swap eigenvalue support */
	itmp1 = Zsupp[2*j];
	Zsupp[2*j] = Zsupp[2*(j+1)];
	Zsupp[2*(j+1)] = itmp1;
	itmp2 = Zsupp[2*j + 1];
	Zsupp[2*j + 1] = Zsupp[2*(j+1) + 1];
	Zsupp[2*(j+1) +1 ] = itmp2;
	/* swap eigenvector */
	memcpy(work, &Z[j*ldz], n*sizeof(double));
	memcpy(&Z[j*ldz], &Z[(j+1)*ldz], n*sizeof(double));
	memcpy(&Z[(j+1)*ldz], work, n*sizeof(double));
      }
    }
  } /* end while */

  return(0);
}




static 
int sort_eigenpairs_global(int m, proc_t *procinfo, val_t *Wstruct, 
			   vec_t *Zstruct)
{
  int              pid   = procinfo->pid;
  int              nproc = procinfo->nproc;
  int              n     = Wstruct->n;
  double *restrict W     = Wstruct->W;
  double *restrict work  = Wstruct->gersch;
  int              ldz   = Zstruct->ldz;
  double *restrict Z     = Zstruct->Z;
  int    *restrict Zsupp = Zstruct->Zsupp;

  double           *minW, *maxW, *minmax; 
  int              i, p, lp, itmp[2];
  bool             sorted;
  MPI_Status       status;
  
  minW   = (double *) malloc(  nproc*sizeof(double));
  assert(minW != NULL);
  maxW   = (double *) malloc(  nproc*sizeof(double));
  assert(maxW != NULL);
  minmax = (double *) malloc(2*nproc*sizeof(double));
  assert(minmax != NULL);

  MPI_Allgather(&W[0], 1, MPI_DOUBLE, minW, 1, MPI_DOUBLE, 
		procinfo->comm); 
  MPI_Allgather(&W[m-1], 1, MPI_DOUBLE, maxW, 1, MPI_DOUBLE, 
		procinfo->comm); 
  
  for (i=0; i<nproc; i++) {
    minmax[2*i]   = minW[i];
    minmax[2*i+1] = maxW[i];
  }

  sorted = true;
  for (i=0; i<2*nproc-1; i++) {
    if (minmax[i] > minmax[i+1]) sorted = false;
  }

  /* Make sure that sorted correctly; ineffective implementation,
   * but usually no or very little swapping should be done here */
  while (sorted == false) {

    sorted = true;

    for (p=1; p<nproc; p++) {

      lp =  p - 1;

      /* swap one pair of eigenvalues and eigenvectors */
      if ((pid == lp || pid == p) && minW[p] < maxW[lp]) {
	if (pid == lp) {
	  W[m-1] = minW[p];
          MPI_Sendrecv(&Z[(m-1)*ldz], n, MPI_DOUBLE, p, lp, 
		       work, n, MPI_DOUBLE, p, p, 
		       procinfo->comm, &status);
	  memcpy(&Z[(m-1)*ldz], work, n*sizeof(double));
	}
	if (pid == p) {
	  W[0]   = maxW[p-1];
          MPI_Sendrecv(&Z[0], n, MPI_DOUBLE, lp, p, 
		       work,  n, MPI_DOUBLE, lp, lp, 
		       procinfo->comm, &status);
	  memcpy(&Z[0], work, n*sizeof(double));
	}
      }

      /* swap eigenvector support as well; 
       * (would better be recomputed here though) */
      if ((pid == lp || pid == p) && minW[p] < maxW[lp]) {
	if (pid == lp) {
          MPI_Sendrecv(&Zsupp[2*(m-1)], 2, MPI_INT, p, lp, 
		       itmp, 2, MPI_INT, p, p, 
		       procinfo->comm, &status);
	  Zsupp[2*(m-1)]     = itmp[0];
	  Zsupp[2*(m-1) + 1] = itmp[1];
	}
	if (pid == p) {
          MPI_Sendrecv(&Zsupp[0], 2, MPI_INT, lp, p, 
		       itmp,  2, MPI_INT, lp, lp, 
		       procinfo->comm, &status);
	  Zsupp[0] = itmp[0];
	  Zsupp[1] = itmp[1];
	}
      }
    }

    /* sort local again */
    sort_eigenpairs_local(m, Wstruct, Zstruct);
    
    /* check again if globally sorted */
    MPI_Allgather(&W[0], 1, MPI_DOUBLE, minW, 1, MPI_DOUBLE, 
		  procinfo->comm); 
    MPI_Allgather(&W[m-1], 1, MPI_DOUBLE, maxW, 1, MPI_DOUBLE, 
		  procinfo->comm); 
    
    for (i=0; i<nproc; i++) {
      minmax[2*i]   = minW[i];
      minmax[2*i+1] = maxW[i];
    }
    
    for (i=0; i<2*nproc-1; i++) {
      if (minmax[i] > minmax[i+1]) sorted = false;
    }
    
    if (pid == 0) printf("sorted = %d\n", sorted);

  } /* end while not sorted */

  free(minW);
  free(maxW);
  free(minmax);

  return(0);
}




/* Routine to sort the eigenpairs */
static 
int sort_eigenpairs(proc_t *procinfo, val_t *Wstruct, vec_t *Zstruct)
{
  /* From inputs */
  int              pid      = procinfo->pid;
  int              n        = Wstruct->n;
  double *restrict W        = Wstruct->W;
  int    *restrict Windex   = Wstruct->Windex;
  int    *restrict iproc    = Wstruct->iproc;
  int    *restrict Zindex   = Zstruct->Zindex;

  /* Others */
  int           im, j;
  sort_struct_t *sort_array;

  /* Make the first nz elements of W contains the eigenvalues
   * associated to the process */
  im = 0;
  for (j=0; j<n; j++) {
    if (iproc[j] == pid) {
      W[im]      = W[j];
      Windex[im] = Windex[j];
      Zindex[im] = Zindex[j];
      im++;
    }
  }

  sort_array = (sort_struct_t *) malloc(im*sizeof(sort_struct_t));

  for (j=0; j<im; j++) {
    sort_array[j].lambda    = W[j]; 
    sort_array[j].local_ind = Windex[j];
    sort_array[j].block_ind = 0;
    sort_array[j].ind       = Zindex[j];
  }

  /* Sort according to Zindex */
  qsort(sort_array, im, sizeof(sort_struct_t), cmpb);

  for (j=0; j<im; j++) {
    W[j]      = sort_array[j].lambda; 
    Windex[j] = sort_array[j].local_ind;
  }

  /* Make sure eigenpairs are sorted locally; this is a very 
   * inefficient way sorting, but in general no or very little 
   * swapping of eigenpairs is expected here */
  sort_eigenpairs_local(im, Wstruct, Zstruct);

  /* Make sure eigenpairs are sorted globally; this is a very 
   * inefficient way sorting, but in general no or very little 
   * swapping of eigenpairs is expected here */
  if (ASSERT_SORTED_EIGENPAIRS == true)
    sort_eigenpairs_global(im, procinfo, Wstruct, Zstruct);

  free(sort_array);

  return(0);
}





/*
 * Refines the eigenvalue to high relative accuracy with
 * respect to the input matrix;
 * Note: In principle this part could be multithreaded too,
 * but it will only rarely be called and not much work
 * is involved, if the eigenvalues are not small in magnitude
 * even no work at all is not uncommon. Therefore using
 * multiple threads seems not to be warranted.
 */
static 
int refine_to_highrac(proc_t *procinfo, char *jobz, double *D,
		      double *E2, in_t *Dstruct, int *nzp,
		      val_t *Wstruct, tol_t *tolstruct)
{
  int              pid    = procinfo->pid;
  bool             wantZ  = (jobz[0]  == 'V' || jobz[0]  == 'v');
  int              n      = Dstruct->n;
  int              nsplit = Dstruct->nsplit;
  int    *restrict isplit = Dstruct->isplit;
  double           spdiam = Dstruct->spdiam;
  int              nz     = *nzp;
  double *restrict W      = Wstruct->W;
  double *restrict Werr   = Wstruct->Werr;
  int    *restrict Windex = Wstruct->Windex;
  int    *restrict iblock = Wstruct->iblock;
  int    *restrict iproc  = Wstruct->iproc;
  double           pivmin = tolstruct->pivmin; 
  double           tol    = 4 * DBL_EPSILON; 
  
  double *work;
  int    *iwork;
  int    ifirst, ilast, offset, info;
  int    i, j, k;
  int    ibegin, iend, isize;
  int    iWbegin, iWend, nbl;

  work  = (double *) malloc( 2*n * sizeof(double) );
  assert (work != NULL);
  iwork = (int *)    malloc( 2*n * sizeof(int)    );
  assert (iwork != NULL);

  ibegin  = 0;
  for (j=0; j<nsplit; j++) {
    
    iend   = isplit[j] - 1;
    isize  = iend - ibegin + 1;
    nbl    = 0;
    
    /* Find eigenvalues in block */
    if (!wantZ) {
      /* Eigenvalues are in W[0...nbl] */
      iWbegin = 0;
      iWend   = 0;
      for (i=0; i<nz; i++) {
	if (nbl == 0 && iblock[i] == j+1) {
	  iWbegin = i;
	  iWend   = i;
	  nbl++;
	  k = i+1;
	  while (k <nz && iblock[k] == j+1) {
	    iWend++; nbl++; k++;
	  }
	}
      }
    } else {
      /* eigenvalues are in W[0...(n-1)] */
      iWbegin = iend   + 1;
      iWend   = ibegin - 1;
      for (i=ibegin; i<=iend; i++) {
	if (nbl == 0 && iproc[i] == pid) {
	  iWbegin = i;
	  iWend   = i;
	  nbl++;
	  k = i+1;
	  while (k <=iend && iproc[k] == pid) {
	    iWend++; nbl++; k++;
	  }
	}
      }
    }

    /* If no eigenvalues for process in block continue */
    if (nbl == 0) {
      ibegin = iend + 1;
      continue;
    }
    
    ifirst  = Windex[iWbegin];
    ilast   = Windex[iWend];
    offset  = Windex[iWbegin] - 1;

    LAPACK(dlarrj)
    (&isize, &D[ibegin], &E2[ibegin], &ifirst, &ilast, &tol, &offset, 
     &W[iWbegin], &Werr[iWbegin], work, iwork, &pivmin, &spdiam, &info);
    assert(info == 0);
    
    ibegin = iend + 1;
  } /* end j */
  
  free(work);
  free(iwork);
  return(0);
}




/* 
 * Compare function for using qsort() on an array of 
 * sort_structs
 */
static 
int cmpb(const void *a1, const void *a2)
{
  sort_struct_t *arg1, *arg2;

  arg1 = (sort_struct_t *) a1;
  arg2 = (sort_struct_t *) a2;

  /* Within block local index decides */
  if (arg1->ind < arg2->ind) 
    return(-1);
  else
    return(1);
}




/*
 * Compare function for using qsort() on an array
 * of doubles
 */
static 
int cmp(const void *a1, const void *a2)
{
  double arg1 = *(double *)a1;
  double arg2 = *(double *)a2;

  if (arg1 < arg2)
    return(-1);
  else
    return(1);
}




/*
 * Routine to communicate eigenvalues such that every process has
 * all computed eigenvalues (iu-il+1) in W; this routine is designed 
 * to be called right after 'PMRRR'.
 */
int PMR_comm_eigvals(MPI_Comm comm, int *nz, int *myfirstp, double *W)
{
  MPI_Comm comm_dup;
  int      nproc;
  double   *work;
  int      *rcount, *rdispl;

  MPI_Comm_dup(comm, &comm_dup);
  MPI_Comm_size(comm_dup, &nproc);

  rcount = (int *) malloc( nproc * sizeof(int) );
  assert(rcount != NULL);
  rdispl = (int *) malloc( nproc * sizeof(int) );
  assert(rdispl != NULL);
  work = (double *) malloc((*nz+1) * sizeof(double));
  assert(work != NULL);

  if (*nz > 0) {
    memcpy(work, W, (*nz) * sizeof(double) );
  }

  MPI_Allgather(nz, 1, MPI_INT, rcount, 1, MPI_INT, comm_dup);

  MPI_Allgather(myfirstp, 1, MPI_INT, rdispl, 1, MPI_INT, comm_dup);
  
  MPI_Allgatherv(work, *nz, MPI_DOUBLE, W, rcount, rdispl,
 		 MPI_DOUBLE, comm_dup);

  MPI_Comm_free(&comm_dup);
  free(rcount);
  free(rdispl);
  free(work);

  return(0);
}




/* Fortran function prototype */
void pmrrr_(char *jobz, char *range, int *n, double  *D,
	    double *E, double *vl, double *vu, int *il, int *iu,
	    int *tryracp, MPI_Fint *comm, int *nz, int *myfirst,
	    double *W, double *Z, int *ldz, int *Zsupp, int* info)
{
  MPI_Comm c_comm = MPI_Comm_f2c(*comm);

  *info = PMRRR(jobz, range, n, D, E, vl, vu, il, iu, tryracp, 
		c_comm, nz, myfirst, W, Z, ldz, Zsupp);
}

void pmr_comm_eigvals_(MPI_Fint *comm, int *nz, int *myfirstp, 
		       double *W, int *info)
{
  MPI_Comm c_comm = MPI_Comm_f2c(*comm);

  *info = PMR_comm_eigvals(c_comm, nz, myfirstp, W);
}
