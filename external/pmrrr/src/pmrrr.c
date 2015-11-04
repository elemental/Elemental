/* Computation of eigenvalues and eigenvectors of a symmetric
 * tridiagonal matrix T, given by its diagonal elements D
 * and its super-/subdiagonal elements E.
 *
 * See INCLUDE/pmrrr.h for more information.
 *
 * Copyright (c) 2010, RWTH Aachen University
 * All rights reserved.
 *
 * Copyright (c) 2015, Jack Poulson
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
#include "pmrrr.h"
#include "pmrrr/plarre.h"
#include "pmrrr/plarrv.h"
#include "pmrrr/structs.h"

static int handle_small_cases
(char*, char*, int*, double*, double*,
 double*, double*, int*, int*, int*,
 MPI_Comm, int*, int*, double*, double*,int*, int*);

static double scale_matrix(in_t*, val_t*, bool);
static void invscale_eigenvalues(val_t*, double, int);
static void clean_up
(MPI_Comm, double*, double*, double*, 
 int*, int*, int*, int*, int*, proc_t*, 
 in_t*, val_t*, vec_t*, tol_t*);
static int refine_to_highrac
(proc_t*, char*, double*, double*,in_t*, int*, val_t*, tol_t*);

static int cmp(const void*, const void*);
static int cmpb(const void*, const void*);

/* 
 * Computation of eigenvalues and eigenvectors of a symmetric
 * tridiagonal matrix T, given by its diagonal elements D
 * and its super-/subdiagonal elements E.
 * See README or 'pmrrr.h' for details.
 */

int pmrrr
(char *jobz, char *range, int *np, double  *D,
 double *E, double *vl, double *vu, int *il,
 int *iu, int *tryracp, MPI_Comm comm, int *nzp,
 int *offsetp, double *W, double *Z, int *ldz, int *Zsupp)
{
  /* Input parameter */
  int  n      = *np;
  bool onlyW = toupper(jobz[0]) == 'N';
  bool wantZ = toupper(jobz[0]) == 'V';
  bool cntval = toupper(jobz[0]) == 'C';
  bool alleig = toupper(range[0]) == 'A';
  bool valeig = toupper(range[0]) == 'V';
  bool indeig = toupper(range[0]) == 'I';

  /* Check input parameters */
  if(!(onlyW  || wantZ  || cntval)) return 1;
  if(!(alleig || valeig || indeig)) return 1;
  if(n <= 0) return 1;
  if (valeig) {
    if(*vu<=*vl) return 1;
  } else if (indeig) {
    if (*il<1 || *il>n || *iu<*il || *iu>n) return 1;
  }
  
  /* MPI & multithreading info */
  int is_init, is_final;
  MPI_Initialized(&is_init);
  MPI_Finalized(&is_final);
  if (is_init!=1 || is_final==1) {
    fprintf(stderr, "ERROR: MPI is not active! (init=%d, final=%d) \n", 
      is_init, is_final);
    return 1;
  }
  MPI_Comm comm_dup;
  MPI_Comm_dup(comm, &comm_dup);
  int nproc, pid, thread_support;
  MPI_Comm_size(comm_dup, &nproc);
  MPI_Comm_rank(comm_dup, &pid);
  MPI_Query_thread(&thread_support);

  int nthreads;
  if ( !(thread_support == MPI_THREAD_MULTIPLE ||
         thread_support == MPI_THREAD_FUNNELED) ) {
    /* Disable multithreading; note: to support multithreading with 
     * MPI_THREAD_SERIALIZED the code must be changed slightly; this 
     * is not supported at the moment */
    nthreads = 1;
  } else {
    char *ompvar = getenv("PMR_NUM_THREADS");
    if (ompvar == NULL) {
      nthreads = DEFAULT_NUM_THREADS;
    } else {
      nthreads = atoi(ompvar);
    }
  }

#if defined(MVAPICH2_VERSION)
  if (nthreads>1) {
    int mv2_affinity=1;
    char *mv2_string = getenv("MV2_ENABLE_AFFINITY");
    if (mv2_string != NULL) 
      mv2_affinity = atoi(mv2_string);
    if (mv2_affinity!=0) {
      nthreads = 1;
      if (pid==0) {
        fprintf(stderr, "WARNING: PMRRR incurs a significant performance penalty when multithreaded with MVAPICH2 with affinity enabled. The number of threads has been reduced to one; please rerun with MV2_ENABLE_AFFINITY=0 or PMR_NUM_THREADS=1 in the future.\n");
        fflush(stderr);
      }
    }
  }
#endif

  /* If only maximal number of local eigenvectors are queried
   * return if possible here */
  *nzp     = 0;
  *offsetp = 0;
  if (cntval) {
    if ( alleig || n < DSTEMR_IF_SMALLER ) {
      *nzp = iceil(n,nproc);
      MPI_Comm_free(&comm_dup);
      return 0;
    } else if (indeig) {
      *nzp = iceil(*iu-*il+1,nproc);
      MPI_Comm_free(&comm_dup);
      return 0;
    }
  }

  /* Check if computation should be done by multiple processes */
  int info;
  if (n < DSTEMR_IF_SMALLER) {
    info = handle_small_cases(jobz, range, np, D, E, vl, vu, il,
			      iu, tryracp, comm, nzp, offsetp, W,
			      Z, ldz, Zsupp);
    MPI_Comm_free(&comm_dup);
    return info;
  }

  /* Allocate memory */
  double *Werr = (double*)malloc(n*sizeof(double)); assert(Werr!=NULL);
  double *Wgap = (double*)malloc(n*sizeof(double)); assert(Wgap!=NULL);
  double *gersch = (double*)malloc(2*n*sizeof(double)); assert(gersch!=NULL);
  int *iblock = (int*)calloc(n,sizeof(int)); assert(iblock!=NULL);
  int *iproc  = (int*)malloc(n*sizeof(int)); assert(iproc!=NULL);
  int *Windex = (int*)malloc(n*sizeof(int)); assert(Windex!=NULL);
  int *isplit = (int*)malloc(n*sizeof(int)); assert(isplit!=NULL);
  int *Zindex = (int*)malloc(n*sizeof(int)); assert(Zindex!=NULL);
  proc_t *procinfo = (proc_t*)malloc(sizeof(proc_t)); assert(procinfo!=NULL);
  in_t *Dstruct = (in_t*)malloc(sizeof(in_t)); assert(Dstruct!=NULL);
  val_t *Wstruct = (val_t*)malloc(sizeof(val_t)); assert(Wstruct!=NULL);
  vec_t *Zstruct = (vec_t*)malloc(sizeof(vec_t)); assert(Zstruct!=NULL);
  tol_t *tolstruct = (tol_t*)malloc(sizeof(tol_t)); assert(tolstruct!=NULL);

  /* Bundle variables into a structures */
  procinfo->pid            = pid;
  procinfo->nproc          = nproc;
  procinfo->comm           = comm_dup;
  procinfo->nthreads       = nthreads;
  procinfo->thread_support = thread_support;

  Dstruct->n      = n;
  Dstruct->D      = D;
  Dstruct->E      = E;
  Dstruct->isplit = isplit;

  Wstruct->n      = n;
  Wstruct->vl     = vl;
  Wstruct->vu     = vu;
  Wstruct->il     = il;
  Wstruct->iu     = iu;
  Wstruct->W      = W;
  Wstruct->Werr   = Werr;
  Wstruct->Wgap   = Wgap;
  Wstruct->Windex = Windex;
  Wstruct->iblock = iblock;
  Wstruct->iproc  = iproc;
  Wstruct->gersch = gersch;

  Zstruct->ldz    = *ldz;
  Zstruct->nz     = 0;
  Zstruct->Z      = Z;
  Zstruct->Zsupp  = Zsupp;
  Zstruct->Zindex = Zindex;

  /* Scale matrix to allowable range, returns 1.0 if not scaled */
  double scale = scale_matrix(Dstruct, Wstruct, valeig);

  /*  Test if matrix warrants more expensive computations which
   *  guarantees high relative accuracy */
  if (*tryracp)
    odrrr(&n, D, E, &info); /* 0 - rel acc */
  else info = -1;

  int i;
  double *Dcopy, *E2copy;
  if (info == 0) {
    /* This case is extremely rare in practice */ 
    tolstruct->split = DBL_EPSILON;
    /* Copy original data needed for refinement later */
    Dcopy = (double*)malloc(n*sizeof(double)); assert(Dcopy!=NULL);
    memcpy(Dcopy, D, n*sizeof(double));  
    E2copy = (double*)malloc(n*sizeof(double)); assert(E2copy!=NULL);
    for (i=0; i<n-1; i++) 
      E2copy[i] = E[i]*E[i];
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
    tolstruct->rtol1 = fmin(1e-2*MIN_RELGAP, tolstruct->rtol1);
    tolstruct->rtol2 = sqrt(DBL_EPSILON)*5.0E-3;
    tolstruct->rtol2 = fmin(5e-6*MIN_RELGAP, tolstruct->rtol2);
    tolstruct->rtol2 = fmax(4.0 * DBL_EPSILON, tolstruct->rtol2);
  }

  /*  Compute all eigenvalues: sorted by block */
  info = plarre(procinfo,jobz,range,Dstruct,Wstruct,tolstruct,nzp,offsetp);
  assert(info == 0);

  /* If just number of local eigenvectors are queried */
  if (cntval & valeig) {    
    clean_up(comm_dup, Werr, Wgap, gersch, iblock, iproc, Windex,
	     isplit, Zindex, procinfo, Dstruct, Wstruct, Zstruct,
	     tolstruct);
    return 0;
  }

  /* If only eigenvalues are to be computed */
  if (!wantZ) {

    /* Refine to high relative with respect to input T */
    if (*tryracp) {
      info = 
        refine_to_highrac
        (procinfo, jobz, Dcopy, E2copy, Dstruct, nzp, Wstruct, tolstruct);
      assert(info == 0);
    }

    /* Sort eigenvalues */
    qsort(W, n, sizeof(double), cmp);

    /* Only keep subset ifirst:ilast */
    int ifirst, ilast, isize;
    int iil = *il;
    int iiu = *iu;
    int ifirst_tmp=iil;
    for (i=0; i<nproc; i++) {
      int chunk  = (iiu-iil+1)/nproc + (i < (iiu-iil+1)%nproc);
      int ilast_tmp;
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
	*offsetp = ifirst - iil;
	*nzp      = isize;
      }
      ifirst_tmp = ilast_tmp + 1;
      ifirst_tmp = imin(ifirst_tmp, iiu + 1);
    }
    if (isize > 0) {
      memmove(W, &W[ifirst-1], *nzp * sizeof(double));
    }

    /* If matrix was scaled, rescale eigenvalues */
    invscale_eigenvalues(Wstruct, scale, *nzp);

    clean_up
    (comm_dup, Werr, Wgap, gersch, iblock, iproc, Windex,
     isplit, Zindex, procinfo, Dstruct, Wstruct, Zstruct, tolstruct);

    return 0;
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

  /* Make the first nz elements of W contains the eigenvalues
   * associated to the process */
  int j, im=0;
  for (j=0; j<n; j++) {
    if (iproc[j] == pid) {
      W[im]      = W[j];
      Windex[im] = Windex[j];
      Zindex[im] = Zindex[j];
      im++;
    }
  }
  
  sort_struct_t *sort_array = (sort_struct_t*)malloc(im*sizeof(sort_struct_t));
  for (j=0; j<im; j++) {
    sort_array[j].lambda = W[j];
    sort_array[j].local_ind = Windex[j];
    sort_array[j].block_ind = 0;
    sort_array[j].ind = Zindex[j];
  }
  qsort(sort_array, im, sizeof(sort_struct_t), cmpb);
  for (j=0; j<im; j++) {
    W[j] = sort_array[j].lambda;
    Windex[j] = sort_array[j].local_ind;
  }
  free(sort_array);

  clean_up(comm_dup, Werr, Wgap, gersch, iblock, iproc, Windex,
	   isplit, Zindex, procinfo, Dstruct, Wstruct, Zstruct,
	   tolstruct);
  if (*tryracp) {
    free(Dcopy);
    free(E2copy);
  }

  return 0;
} /* end pmrrr */

static int cmpb(const void *a1, const void *a2)
{
  sort_struct_t *arg1, *arg2;
  
  arg1 = (sort_struct_t*) a1;
  arg2 = (sort_struct_t*) a2;
  
  if (arg1->ind < arg2->ind)
    return -1;
  else
    return 1;
}

/*
 * Free's on allocated memory of pmrrr routine
 */
static  
void clean_up
(MPI_Comm comm, double *Werr, double *Wgap,
 double *gersch, int *iblock, int *iproc,
 int *Windex, int *isplit, int *Zindex,
 proc_t *procinfo, in_t *Dstruct,
 val_t *Wstruct, vec_t *Zstruct, tol_t *tolstruct)
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
int handle_small_cases
(char *jobz, char *range, int *np, double  *D,
 double *E, double *vlp, double *vup, int *ilp,
 int *iup, int *tryracp, MPI_Comm comm, int *nzp,
 int *myfirstp, double *W, double *Z, int *ldzp, int *Zsupp)
{
  bool cntval = toupper(jobz[0]) == 'C';
  bool onlyW = toupper(jobz[0]) == 'N';
  bool wantZ = toupper(jobz[0]) == 'V';
  bool indeig = toupper(range[0]) == 'I';
  int n       = *np;
  int ldz_tmp = *np;
  int ldz     = *ldzp;

  int nproc, pid;
  MPI_Comm_size(comm, &nproc);
  MPI_Comm_rank(comm, &pid);
  
  int lwork, liwork;
  double *Z_tmp;
  if (onlyW) {
    lwork  = 12*n;
    liwork =  8*n;
  } else if (cntval) {
    lwork  = 18*n;
    liwork = 10*n;
  } else if (wantZ) {
    lwork  = 18*n;
    liwork = 10*n;
    int itmp;
    if (indeig) itmp = *iup-*ilp+1;
    else        itmp = n;
    Z_tmp = (double*)malloc(n*itmp*sizeof(double)); assert(Z_tmp!=NULL);
  } else {
    return 1;
  }

  double *work = (double*)malloc(lwork*sizeof(double)); assert(work != NULL);
  int *iwork = (int*)malloc(liwork*sizeof(int)); assert(iwork!=NULL);

  if (cntval) {
    /* Note: at the moment, jobz="C" should never get here, since
     * it is blocked before. */
    int m, info, MINUSONE=-1;
    double cnt;
    odstmr("V", "V", np, D, E, vlp, vup, ilp, iup, &m, W, &cnt,
	   &ldz_tmp, &MINUSONE, Zsupp, tryracp, work, &lwork, iwork,
	   &liwork, &info);
    assert(info == 0);
    
    *nzp = (int) ceil(cnt/nproc);
    free(work); free(iwork);
    return 0;
  }

  int m, info;
  odstmr(jobz, range, np, D, E, vlp, vup, ilp, iup, &m, W, Z_tmp,
	 &ldz_tmp, np, Zsupp, tryracp, work, &lwork, iwork,
	 &liwork, &info);
  assert(info == 0);

  int chunk   = iceil(m,nproc);
  int myfirst = imin(pid * chunk, m);
  int mylast  = imin((pid+1)*chunk - 1, m - 1);
  int mysize  = mylast - myfirst + 1;

  if (mysize > 0) {
    memmove(W, &W[myfirst], mysize*sizeof(double));
    if (wantZ) {
      if (ldz == ldz_tmp) {
	/* copy everything in one chunk */
	memcpy(Z, &Z_tmp[myfirst*ldz_tmp], n*mysize*sizeof(double));
      } else {
	/* copy each vector seperately */
        int i;
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

  return 0;
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

  /* Set some machine dependent constants */
  double smlnum = DBL_MIN / DBL_EPSILON;
  double bignum = 1.0 / smlnum;
  double rmin   = sqrt(smlnum);
  double rmax   = fmin(sqrt(bignum), 1.0 / sqrt(sqrt(DBL_MIN)));

  /*  Scale matrix to allowable range */
  double scale = 1.0;
  double T_norm = odnst("M", &n, D, E);  /* returns max(|T(i,j)|) */
  if (T_norm > 0 && T_norm < rmin) {
    scale = rmin / T_norm;
  } else if (T_norm > rmax) {
    scale = rmax / T_norm;
  }

  if (scale != 1.0) {  /* FP cmp okay */
    /* Scale matrix and matrix norm */
    int itmp = n-1;
    int IONE = 1;
    odscal(&n,    &scale, D, &IONE);
    odscal(&itmp, &scale, E, &IONE);
    if (valeig == true) {
      /* Scale eigenvalue bounds */
      *vl *= scale;
      *vu *= scale;
    }
  } /* end scaling */

  return scale;
}

/*
 * If matrix scaled, rescale eigenvalues
 */
static 
void invscale_eigenvalues(val_t *Wstruct, double scale, int size)
{
  if (scale != 1.0) {  /* FP cmp okay */
    double *vl = Wstruct->vl;
    double *vu = Wstruct->vu;
    double *W  = Wstruct->W;

    double invscale = 1.0 / scale;
    *vl *= invscale;
    *vu *= invscale;
    int IONE = 1;
    odscal(&size, &invscale, W, &IONE);
  }
}

/*
 * Refines the eigenvalue to high relative accuracy with
 * respect to the input matrix;
 * Note: In principle this part could be fully parallel too,
 * but it will only rarely be called and not much work
 * is involved, if the eigenvalues are not small in magnitude
 * even no work at all is not uncommon. 
 */
static 
int refine_to_highrac
(proc_t *procinfo, char *jobz, double *D, double *E2, in_t *Dstruct, int *nzp,
 val_t *Wstruct, tol_t *tolstruct)
{
  int              n      = Dstruct->n;
  int              nsplit = Dstruct->nsplit;
  int    *restrict isplit = Dstruct->isplit;
  double           spdiam = Dstruct->spdiam;
  double *restrict W      = Wstruct->W;
  double *restrict Werr   = Wstruct->Werr;
  
  double *work = (double*)malloc(2*n*sizeof(double)); assert(work!=NULL);
  int *iwork = (int*)malloc(2*n*sizeof(int)); assert(iwork!=NULL);

  int j, ibegin=0;
  for (j=0; j<nsplit; j++) {
    int iend = isplit[j] - 1;
    int isize = iend - ibegin + 1;
    int nbl = isize;
    if (nbl == 1) {
      ibegin = iend + 1;
      continue;
    }
    
    int ifirst=1, ilast=nbl, offset=0, info;
    double tol = 4*DBL_EPSILON;
    double pivmin = tolstruct->pivmin; 
    odrrj(&isize, &D[ibegin], &E2[ibegin], &ifirst, &ilast, &tol,
	  &offset, &W[ibegin], &Werr[ibegin], work, iwork, &pivmin,
	  &spdiam, &info);
    assert(info == 0);
    
    ibegin = iend + 1;
  } /* end j */
  
  free(work);
  free(iwork);
  return 0;
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
    return -1;
  else
    return 1;
}

/*
 * Routine to communicate eigenvalues such that every process has
 * all computed eigenvalues (iu-il+1) in W; this routine is designed 
 * to be called right after 'pmrrr'.
 */
int PMR_comm_eigvals(MPI_Comm comm, int *nz, int *myfirstp, double *W)
{
  MPI_Comm comm_dup;
  MPI_Comm_dup(comm, &comm_dup);
  int nproc;
  MPI_Comm_size(comm_dup, &nproc);

  int *rcount = (int*)malloc(nproc*sizeof(int)); assert(rcount!=NULL);
  int *rdispl = (int*)malloc(nproc*sizeof(int)); assert(rdispl!=NULL);
  double *work = (double*)malloc((*nz+1)*sizeof(double)); assert(work!=NULL);

  if (*nz > 0)
    memcpy(work, W, (*nz)*sizeof(double) );

  MPI_Allgather(nz, 1, MPI_INT, rcount, 1, MPI_INT, comm_dup);

  MPI_Allgather(myfirstp, 1, MPI_INT, rdispl, 1, MPI_INT, comm_dup);
  
  MPI_Allgatherv
  (work, *nz, MPI_DOUBLE, W, rcount, rdispl, MPI_DOUBLE, comm_dup);

  MPI_Comm_free(&comm_dup);
  free(rcount);
  free(rdispl);
  free(work);

  return 0;
}

void pmr_comm_eigvals_
(MPI_Fint *comm, int *nz, int *myfirstp, double *W, int *info)
{
  MPI_Comm c_comm = MPI_Comm_f2c(*comm);
  *info = PMR_comm_eigvals(c_comm, nz, myfirstp, W);
}
