/* Parallel computation of eigenvectors and symmetric tridiagonal 
 * matrix T, which is preprocessed by the routine 'plarre'.
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
#include <assert.h>
#include <pthread.h>
#include "mpi.h"
#include "pmrrr.h"
#include "plarrv.h"
#include "process_task.h"
#include "global.h"
#include "rrr.h"
#include "queue.h"
#include "structs.h"
#include "counter.h"


static int assign_to_proc(proc_t *procinfo, in_t *Dstruct,
			  val_t *Wstruct, vec_t *Zstruct, int *nzp,
			  int *myfirstp);
static int cmpa(const void*, const void*);
static int init_workQ(proc_t *procinfo, in_t *Dstruct,
			   val_t *Wstruct, int *nzp,
			   workQ_t *workQ);
static void *empty_workQ(void*);
static workQ_t *create_workQ();
static void destroy_workQ(workQ_t*);
static auxarg3_t *create_auxarg3(int, proc_t*, val_t*, vec_t*,
				 tol_t*, workQ_t*, counter_t*);
static void retrieve_auxarg3(auxarg3_t*, int*, proc_t**, val_t**,
			     vec_t**, tol_t**, workQ_t**, 
			     counter_t**);



/*
 * Computation of eigenvectors of a symmetric tridiagonal
 */
int plarrv(proc_t *procinfo, in_t *Dstruct, val_t *Wstruct,
	   vec_t *Zstruct, tol_t *tolstruct, int *nzp,
	   int *myfirstp)
{
  /* Input variables */
  int            nthreads = procinfo->nthreads;
  int            n        = Dstruct->n;
  double         *W       = Wstruct->W;

  /* Work space */
  double         *Wshifted;
 
  /* Multi-threading */
  pthread_t      *threads;   
  pthread_attr_t attr;
  void           *status;
  auxarg3_t      *auxarg;
  counter_t      *num_left;
  
  /* Others */
  int            info, i;
  workQ_t        *workQ;

  /* Allocate work space and copy eigenvalues */
  Wshifted = (double *) malloc( n * sizeof(double) );
  assert(Wshifted != NULL);

  memcpy(Wshifted, W, n*sizeof(double));
  Wstruct->Wshifted = Wshifted;

  threads = (pthread_t *) malloc(nthreads * sizeof(pthread_t));
  assert(threads != NULL);

  /* Assign eigenvectors to processes */
  assign_to_proc(procinfo, Dstruct, Wstruct, Zstruct, nzp,
		 myfirstp);

  /* Create work queue Q, counter, threads to empty Q */
  workQ    = create_workQ( );
  num_left = PMR_create_counter(*nzp);

  threads[0] = pthread_self();
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
  pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM); 

  for (i=1; i<nthreads; i++) {
    auxarg = create_auxarg3(i, procinfo, Wstruct, Zstruct, tolstruct,
			    workQ, num_left);
    info = pthread_create(&threads[i], &attr, empty_workQ, 
			  (void *) auxarg);
    assert(info == 0);
  }

  /* Initialize work queue of process */
  info = init_workQ(procinfo, Dstruct, Wstruct, nzp, workQ);
  assert(info == 0);

  /* Empty the work queue */
  auxarg = create_auxarg3(0, procinfo, Wstruct, Zstruct, tolstruct,
			  workQ, num_left);
  status = empty_workQ((void *) auxarg);
  assert(status == NULL);

  /* Join all the worker thread */
  for (i=1; i<nthreads; i++) {
    info = pthread_join(threads[i], &status);
    assert(info == 0 && status == NULL);
  }
  
  /* Clean up */
  free(Wshifted);
  free(threads);
  pthread_attr_destroy(&attr);
  destroy_workQ(workQ);
  PMR_destroy_counter(num_left);

  return(0);
}




/* 
 * Assign the computation of eigenvectors to the processes
 */
static  
int assign_to_proc(proc_t *procinfo, in_t *Dstruct, val_t *Wstruct,
		   vec_t *Zstruct, int *nzp, int *myfirstp)
{
  /* From inputs */
  int              pid     = procinfo->pid;
  int              nproc   = procinfo->nproc;
  double *restrict L       = Dstruct->E;
  int    *restrict isplit  = Dstruct->isplit;
  int              n       = Wstruct->n;
  int              il      = *(Wstruct->il);
  int              iu      = *(Wstruct->iu);
  double *restrict W       = Wstruct->W;
  int    *restrict Windex  = Wstruct->Windex;
  int    *restrict iblock  = Wstruct->iblock;
  int    *restrict iproc   = Wstruct->iproc;
  int    *restrict Zindex  = Zstruct->Zindex;
  
  /* Others */
  int               i, id, j, k, isize, iblk, ishift,
                    ibegin, iend, chunk, ind;
  double            sigma;
  sort_struct_t     *array;

  array = (sort_struct_t *) malloc(n*sizeof(sort_struct_t));

  for (i=0; i<n; i++) {
    /* Find shift of block */
    iblk                    = iblock[i];
    ishift                  = isplit[iblk-1] - 1;
    sigma                   = L[ishift];
    /* Apply shift so that unshifted eigenvalues can be sorted */
    array[i].lambda    = W[i] + sigma; 
    array[i].local_ind = Windex[i];
    array[i].block_ind = iblk;
    array[i].ind       = i;
  }
  /* Alternatively could sort list of pointers, which
     requires less data copying */

  qsort(array, n, sizeof(sort_struct_t), cmpa);

  for (j = 0; j < il-1; j++ ) {
    ind = array[j].ind;
    iproc[ind]  = -1;
    Zindex[ind] = -1;
  }

  isize = iu - il + 1;

  ibegin = il - 1;
  for (id=0; id<nproc; id++) {

    chunk = imax(1, isize/nproc + (id < isize%nproc));
    
    if (id==nproc-1) {
      iend = iu - 1;
    } else {
      iend = ibegin + chunk - 1;
      iend = imin(iend, iu -1);
    }

    k = 0;
    for (j=ibegin; j<=iend; j++) {
      ind = array[j].ind;
      iproc[ind]  = id;
      Zindex[ind] = k;
      k++;
    }

    if (id == pid) {
      *myfirstp   = ibegin - il + 1;
      *nzp        = iend - ibegin + 1;
      Zstruct->nz = *nzp; 
    }

    ibegin = iend + 1;
    ibegin = imin(ibegin, iu);
  } /* end id */

  for (j = iend+1; j < n; j++ ) {
    ind = array[j].ind;
    iproc[ind]  = -1;
    Zindex[ind] = -1;
  }

  free(array);
  return(0);
}




/* 
 * Compare function for using qsort() on an array of 
 * sort_structs
 */
static 
int cmpa(const void *a1, const void *a2)
{
  sort_struct_t *arg1, *arg2;

  arg1 = (sort_struct_t *) a1;
  arg2 = (sort_struct_t *) a2;

  /* Within block local index decides */
  if (arg1->block_ind == arg2->block_ind) {
    return (arg1->local_ind - arg2->local_ind);
  } else {
    if (arg1->lambda < arg2->lambda) {
      return(-1);
    } else if (arg1->lambda > arg2->lambda) {
      return(1);
    } else {
      if (arg1->local_ind < arg2->local_ind)
	return(-1);
      else
	return(1);
    }
  }
}




/*
 * Initialize work queue by putting all tasks for the process
 * into the work queue.
 */
static 
int init_workQ(proc_t *procinfo, in_t *Dstruct, val_t *Wstruct,
	       int *nzp, workQ_t *workQ)
{
  /* Input arguments */
  int              pid      = procinfo->pid;
  int              nproc    = procinfo->nproc;
  int              nthreads = procinfo->nthreads;
  double *restrict D        = Dstruct->D;
  double *restrict L        = Dstruct->E;
  int              nsplit   = Dstruct->nsplit;
  int    *restrict isplit   = Dstruct->isplit;
  double *restrict W        = Wstruct->W;
  double *restrict Werr     = Wstruct->Werr;
  double *restrict Wgap     = Wstruct->Wgap;
  int    *restrict iproc    = Wstruct->iproc;
  double *restrict Wshifted = Wstruct->Wshifted;
  double *restrict gersch   = Wstruct->gersch;
  int              nz       = *nzp;

  /* Loop over blocks */
  int              ibegin, iend, isize, iWbegin, iWend, nbl;
  double           sigma, gl, gu, avggap, spdiam;
  double *restrict DL;
  double *restrict DLL;
  rrr_t            *RRR, *RRR_parent;

  /* Splitting into singletons and cluster */
  int              new_first, new_last, new_size;
  int              sn_first,  sn_last,  sn_size;
  int              cl_first,  cl_last,  cl_size;
  bool             task_inserted;
  int              max_size, left_pid, right_pid;
  double           lgap;
 
  /* Others */
  int              i, j, k, l;
  double           tmp;
  task_t           *task;

  /* Loop over blocks */
  ibegin  = 0;
  for ( j=0; j<nsplit; j++ ) {

    iend   = isplit[j] - 1;
    isize  = iend - ibegin + 1;
    sigma  = L[iend];

    /* Use Gerschgorin disks to find spectral diameter */
    gl = gersch[2*ibegin    ];
    gu = gersch[2*ibegin + 1];
    for (i=ibegin+1; i<iend; i++) {
      gl = fmin(gl, gersch[2*i    ]);
      gu = fmax(gu, gersch[2*i + 1]);
    }
    spdiam = gu - gl;
    avggap = spdiam / (isize-1);

    /* Find eigenvalues in block */
    nbl = 0;
    iWbegin = iend   + 1;
    iWend   = ibegin - 1;
    for (i=ibegin; i<=iend; i++) {
      if (nbl == 0 && iproc[i] == pid) {
	iWbegin = i;
	iWend   = i;
	nbl++;
	k = i+1;
	while (k <=iend && iproc[k] == pid) {
	  iWend++;
	  nbl++;
	  k++;
	}
	/* iWend = iWbegin + nbl - 1; instead of incrementing in loop */
      }
    }

    /* If no eigenvalues for process in block continue */
    if (nbl == 0) {
      ibegin = iend + 1;
      continue;
    }

    /* Compute DL and DLL for later computation of singletons
     * (freed when all singletons of root are computed) */
    DL  = (double *) malloc(isize * sizeof(double));
    assert(DL != NULL);
    
    DLL = (double *) malloc(isize * sizeof(double));
    assert(DLL != NULL);

    for (i=0; i<isize-1; i++) {
      tmp    = D[i+ibegin]*L[i+ibegin];
      DL[i]  = tmp;
      DLL[i] = tmp*L[i+ibegin];
    }

    RRR = PMR_create_rrr(&D[ibegin], &L[ibegin], DL, DLL, isize, 0);
    PMR_increment_rrr_dependencies(RRR);
    
    /* In W apply shift of current block to eigenvalues
     * to get unshifted values w.r.t. T */
    for (i=ibegin; i<=iend; i++) {
      W[i] += sigma;
    }

    /* Split eigenvalues of block into singletons and clusters
     * and add them to process work queue */
    max_size = imax(1, nz/nthreads);
    task_inserted = false;
    new_first = ibegin;
    for (i=ibegin; i<=iend; i++) {    

      if (i == iend)
	new_last = i;
      else if (Wgap[i] >= MIN_RELGAP*fabs(Wshifted[i]))
	new_last = i;
      else
	continue;

      /* Skip rest if no eigenvalues of process */
      if (new_first > iWend || new_last < iWbegin) {
	new_first = i + 1;
	continue;
      }

      new_size = new_last - new_first + 1;
      
      if (new_size == 1) {
	/* Singleton was found */
	
	if (new_first < iWbegin || new_first > iWend) {
	  new_first = i + 1;
	  continue;
	} else {
	  if (new_first==iWbegin || task_inserted==true) {
	    /* Initialize new singleton task */
	    sn_first = new_first;
	    sn_last  = new_first;
	    sn_size  = 1;
	  } else {
	    /* Extend singleton task by one */
	    sn_last++;
	    sn_size++;
	  }
	}

	/* Insert task if ... */
	if (i==iWend || sn_size>=max_size ||
	    Wgap[i+1] < MIN_RELGAP*fabs(Wshifted[i+1])) {

	  if (sn_first == ibegin) {
	    lgap = fmax(0.0, W[ibegin] - Werr[ibegin] - gl );
	  } else {
	    lgap = Wgap[sn_first-1];
	  }

	  PMR_increment_rrr_dependencies(RRR);

	  task = PMR_create_s_task(sn_first, sn_last, 1, ibegin, 
				   iend, spdiam, lgap, RRR);
	  
 	  PMR_insert_task_at_back(workQ->s_queue, task);
	  
	  task_inserted = true;
	} else {
	  task_inserted = false;
	}

      } else {
	/* Cluster was found */

	cl_first = new_first;
	cl_last  = new_last;
	cl_size  = new_size;

	/* Split cluster into clusters by absolut criterion */
	if (cl_size > 3) {

	  /* Split cluster to smaller clusters [cl_first:cl_last] */
	  for (k=new_first+1; k<new_last; k++) {

	    if (k == new_last-1)
	      cl_last = new_last;
             else if (k != cl_first && Wgap[k] > 0.8*avggap) 
	      cl_last = k;
	    else
	      continue;

	    /* Skip cluster if no eigenvalues of process in it */
	    if (cl_last < iWbegin || cl_first > iWend) {
	      cl_first = k + 1;
	      continue;
	    }

	    /* Record left gap of cluster */
	    if (cl_first == ibegin) {
	      lgap = fmax(0.0, W[ibegin] - Werr[ibegin] - gl);
	    } else {
	      lgap = Wgap[cl_first-1];
	    }

	    /* Determine processes involved in processing the cluster */
	    left_pid  = nproc-1;
	    right_pid = 0;
	    for (l=cl_first; l<=cl_last; l++) {
	      if (iproc[l] != -1) {
		left_pid  = imin(left_pid,  iproc[l]);
		right_pid = imax(right_pid, iproc[l]);
	      }
	    }

	    RRR_parent = PMR_create_rrr(&D[ibegin], &L[ibegin], 
					NULL, NULL, isize, 0);

	    task = PMR_create_c_task(cl_first, cl_last, 1, ibegin, 
				     iend, spdiam, lgap, iWbegin, 
				     iWend, left_pid, right_pid, 
				     RRR_parent);

	    /* Insert task into queue, depending if cluster need
	     * communication with other processes */
	    if (left_pid != right_pid)
	      PMR_insert_task_at_back(workQ->r_queue, task);
	    else
	      PMR_insert_task_at_back(workQ->c_queue, task);
	    
	    cl_first = k + 1;
	  } /* end k */

	} else {
	  /* Cluster is too small to split, so insert it to queue */

	  /* Record left gap of cluster */
	  if (cl_first == ibegin) {
	    lgap = fmax(0.0, W[ibegin] - Werr[ibegin] - gl );
	  } else {
	    lgap = Wgap[cl_first-1];
	  }

	  /* Determine processes involved */
	  left_pid  = nproc-1;
	  right_pid = 0;
	  for (l=cl_first; l<=cl_last; l++) {
	    if (iproc[l] != -1) {
	      left_pid  = imin(left_pid,  iproc[l]);
	      right_pid = imax(right_pid, iproc[l]);
	    }
	  }

	  RRR_parent = PMR_create_rrr(&D[ibegin], &L[ibegin], 
				      NULL, NULL, isize, 0);

	  task = PMR_create_c_task(cl_first, cl_last, 1, ibegin, 
				   iend, spdiam, lgap, iWbegin, iWend,
				   left_pid, right_pid, RRR_parent);

	  /* Insert task into queue, depending if cluster need
	   * communication with other processes */
	  if (left_pid != right_pid)
	    PMR_insert_task_at_back(workQ->r_queue, task);
	  else
	    PMR_insert_task_at_back(workQ->c_queue, task);
	  
	}
	task_inserted = true;

      } /* end new_size */

      new_first = i + 1;
    } /* end of splitting eigenvalues into tasks */

    /* Set flag in RRR that last singleton is created */
    PMR_set_parent_processed_flag(RRR);
    PMR_try_destroy_rrr(RRR);

    ibegin = iend + 1;
  } /* end loop over blocks */ 

  return(0);
}




/*
 * Processes all the tasks put in the work queue.
 */
static 
void *empty_workQ(void *argin)
{
  /* input arguments */
  int          tid;
  proc_t       *procinfo;
  val_t        *Wstruct;
  vec_t        *Zstruct;
  tol_t        *tolstruct;
  workQ_t *workQ;
  counter_t    *num_left;
  int          n;

  /* others */
  task_t       *task;
  double       *work;
  int          *iwork;

  /* retrieve necessary arguments from structures */
  retrieve_auxarg3((auxarg3_t *) argin, &tid, &procinfo, &Wstruct,
		   &Zstruct, &tolstruct, &workQ, &num_left);

  n        = Wstruct->n;

  /* max. needed double precision work space: dlar1v */
  work      = (double *) malloc(4*n * sizeof(double));
  assert(work != NULL);

  /* max. needed double precision work space: dlarrb */
  iwork     = (int *)    malloc(2*n * sizeof(int)   );
  assert(iwork != NULL);


  /* while loop to empty the work queue */
  while (PMR_get_counter_value(num_left) > 0) {

    /* empty r-queue before processing other tasks */
    PMR_process_r_queue(tid, procinfo, Wstruct, Zstruct, tolstruct,
			workQ, num_left, work, iwork);

    task = PMR_remove_task_at_front(workQ->s_queue);
    if ( task != NULL ) {
      assert(task->flag == SINGLETON_TASK_FLAG);

      PMR_process_s_task((singleton_t *) task->data, tid, procinfo,
			 Wstruct, Zstruct, tolstruct, num_left, 
			 work, iwork);
      free(task);
      continue;
    }
    
    task = PMR_remove_task_at_front(workQ->c_queue);
    if ( task != NULL ) {
      assert(task->flag == CLUSTER_TASK_FLAG);

      PMR_process_c_task((cluster_t *) task->data, tid, procinfo,
			 Wstruct, Zstruct, tolstruct, workQ,
			 num_left, work, iwork);
      free(task);
      continue;
    }
    
  } /* end while */

  free(work);
  free(iwork);

  return(NULL);
}




static workQ_t *create_workQ()
{
  workQ_t *wq;

  wq = (workQ_t *) malloc(sizeof(workQ_t));

  wq->r_queue = PMR_create_empty_queue();
  wq->s_queue = PMR_create_empty_queue();
  wq->c_queue = PMR_create_empty_queue();

  return(wq);
}




static void destroy_workQ(workQ_t *wq)
{
  PMR_destroy_queue(wq->r_queue);
  PMR_destroy_queue(wq->s_queue);
  PMR_destroy_queue(wq->c_queue);
  free(wq);
}




static auxarg3_t*
create_auxarg3(int tid, proc_t *procinfo, val_t *Wstruct,
	       vec_t *Zstruct, tol_t *tolstruct,
	       workQ_t *workQ, counter_t *num_left)
{
  auxarg3_t *arg;

  arg = (auxarg3_t *) malloc( sizeof(auxarg3_t) );
  assert(arg != NULL);

  arg->tid         = tid;
  arg->procinfo    = procinfo;
  arg->Wstruct     = Wstruct;
  arg->Zstruct     = Zstruct;
  arg->tolstruct   = tolstruct; 
  arg->workQ  = workQ;
  arg->num_left    = num_left;

  return(arg);
}




static void 
retrieve_auxarg3(auxarg3_t *arg, int *tid, proc_t **procinfo,
		 val_t **Wstruct, vec_t **Zstruct,
		 tol_t **tolstruct, workQ_t **workQ,
		 counter_t **num_left)
{
  *tid         = arg->tid;
  *procinfo    = arg->procinfo;
  *Wstruct     = arg->Wstruct;
  *Zstruct     = arg->Zstruct;
  *tolstruct   = arg->tolstruct;
  *workQ  = arg->workQ;
  *num_left    = arg->num_left;

  free(arg);
}
