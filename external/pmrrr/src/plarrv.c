/* Parallel computation of eigenvectors and symmetric tridiagonal 
 * matrix T, which is preprocessed by the routine 'plarre'.
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
#include "pmrrr/plarrv.h"
#include "pmrrr/process_task.h"
#include "pmrrr/rrr.h"
#include "pmrrr/queue.h"
#include "pmrrr/structs.h"
#include "pmrrr/counter.h"

static int assign_to_proc
(proc_t *procinfo, in_t *Dstruct, val_t *Wstruct, vec_t *Zstruct, 
 int *nzp, int *myfirstp);
static int cmpa(const void*, const void*);
static int init_workQ
(proc_t *procinfo, in_t *Dstruct, val_t *Wstruct, int *nzp, workQ_t *workQ);
static void *empty_workQ(void*);
static workQ_t *create_workQ();
static void destroy_workQ(workQ_t*);
static auxarg3_t *create_auxarg3
(int, proc_t*, val_t*, vec_t*, tol_t*, workQ_t*, counter_t*);
static void retrieve_auxarg3
(auxarg3_t*, int*, proc_t**, val_t**, vec_t**, tol_t**, workQ_t**, counter_t**);

/*
 * Computation of eigenvectors of a symmetric tridiagonal
 */
#ifndef DISABLE_PTHREADS
int plarrv
(proc_t *procinfo, in_t *Dstruct, val_t *Wstruct,
 vec_t *Zstruct, tol_t *tolstruct, int *nzp, int *myfirstp)
{
  int     nthreads = procinfo->nthreads;
  int     n        = Dstruct->n;
  double  *W       = Wstruct->W;

  /* Allocate work space and copy eigenvalues */
  double *Wshifted = (double*)malloc(n*sizeof(double));
  assert(Wshifted != NULL);

  memcpy(Wshifted, W, n*sizeof(double));
  Wstruct->Wshifted = Wshifted;

  pthread_t *threads = (pthread_t*)malloc(nthreads*sizeof(pthread_t));
  assert(threads != NULL);

  /* Assign eigenvectors to processes */
  assign_to_proc(procinfo, Dstruct, Wstruct, Zstruct, nzp, myfirstp);

  /* Create work queue Q, counter, threads to empty Q */
  workQ_t *workQ = create_workQ();
  counter_t *num_left = PMR_create_counter(*nzp);

  threads[0] = pthread_self();
  pthread_attr_t attr;
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
  pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM); 

  int i;
  for (i=1; i<nthreads; i++) {
    auxarg3_t *auxarg = 
      create_auxarg3(i, procinfo, Wstruct, Zstruct, tolstruct, workQ, num_left);
    int info = pthread_create(&threads[i], &attr, empty_workQ, (void*)auxarg);
    assert(info == 0);
  }

  /* Initialize work queue of process */
  int info = init_workQ(procinfo, Dstruct, Wstruct, nzp, workQ);
  assert(info == 0);

  /* Empty the work queue */
  auxarg3_t *auxarg = 
    create_auxarg3(0, procinfo, Wstruct, Zstruct, tolstruct, workQ, num_left);
  void *status = empty_workQ((void*)auxarg);
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

  return 0;
}
#else
int plarrv
(proc_t *procinfo, in_t *Dstruct, val_t *Wstruct,
 vec_t *Zstruct, tol_t *tolstruct, int *nzp, int *myfirstp)
{
  int     n  = Dstruct->n;
  double  *W = Wstruct->W;

  /* Allocate work space and copy eigenvalues */
  double *Wshifted = (double*)malloc(n*sizeof(double));
  assert(Wshifted != NULL);

  memcpy(Wshifted, W, n*sizeof(double));
  Wstruct->Wshifted = Wshifted;

  /* Assign eigenvectors to processes */
  assign_to_proc(procinfo, Dstruct, Wstruct, Zstruct, nzp, myfirstp);

  /* Create work queue Q, counter, threads to empty Q */
  workQ_t *workQ = create_workQ();
  counter_t *num_left = PMR_create_counter(*nzp);

  /* Initialize work queue of process */
  int info = init_workQ(procinfo, Dstruct, Wstruct, nzp, workQ);
  assert(info == 0);

  /* Empty the work queue */
  auxarg3_t *auxarg = 
    create_auxarg3(0, procinfo, Wstruct, Zstruct, tolstruct, workQ, num_left);
  void *status = empty_workQ((void*)auxarg);
  assert(status == NULL);

  /* Clean up */
  free(Wshifted);
  destroy_workQ(workQ);
  PMR_destroy_counter(num_left);

  return 0;
}
#endif

/* 
 * Assign the computation of eigenvectors to the processes
 */
static  
int assign_to_proc
(proc_t *procinfo, in_t *Dstruct, val_t *Wstruct,
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
  
  sort_struct_t *array = (sort_struct_t *) malloc(n*sizeof(sort_struct_t));
  
  int i;
  for (i=0; i<n; i++) {
    /* Find shift of block */
    int iblk = iblock[i];
    int ishift  = isplit[iblk-1] - 1;
    /* Apply shift so that unshifted eigenvalues can be sorted */
    array[i].lambda    = W[i] + L[ishift]; 
    array[i].local_ind = Windex[i];
    array[i].block_ind = iblk;
    array[i].ind       = i;
  }
  /* Alternatively could sort list of pointers, which
     requires less data copying */

  qsort(array, n, sizeof(sort_struct_t), cmpa);

  /* Mark eigenvectors that do not need to be computed */
  int j;
  for (j = 0; j < il-1; j++ ) {
    iproc[array[j].ind]  = -1;
    Zindex[array[j].ind] = -1;
  }

  int isize = iu - il + 1;

  int ibegin=il-1, iend;
  int id;
  for (id=0; id<nproc; id++) {

    int chunk = imax(1, isize/nproc + (id < isize%nproc));
    
    if (id==nproc-1) {
      iend = iu - 1;
    } else {
      iend = ibegin + chunk - 1;
      iend = imin(iend, iu -1);
    }

    int k = 0;
    for (j=ibegin; j<=iend; j++) {
      iproc[array[j].ind]  = id;
      Zindex[array[j].ind] = k;
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
    iproc[array[j].ind]  = -1;
    Zindex[array[j].ind] = -1;
  }
  
  free(array);
  return 0;
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
      return -1;
    } else if (arg1->lambda > arg2->lambda) {
      return 1;
    } else {
      if (arg1->local_ind < arg2->local_ind)
        return -1;
      else
        return 1;
    }
  }
}

/*
 * Initialize work queue by putting all tasks for the process
 * into the work queue.
 */
static 
int init_workQ
(proc_t *procinfo, in_t *Dstruct, val_t *Wstruct, int *nzp, workQ_t *workQ)
{
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
  int i, j, k, l;
  int ibegin  = 0;
  for ( j=0; j<nsplit; j++ ) {
    int iend = isplit[j] - 1;
    int isize = iend - ibegin + 1;
    double sigma = L[iend];

    /* Use Gerschgorin disks to find spectral diameter */
    double gl = gersch[2*ibegin    ];
    double gu = gersch[2*ibegin + 1];
    for (i=ibegin+1; i<iend; i++) {
      gl = fmin(gl, gersch[2*i    ]);
      gu = fmax(gu, gersch[2*i + 1]);
    }
    double spdiam = gu - gl;
    double avggap = spdiam / (isize-1);

    /* Find eigenvalues in block */
    int nbl = 0;
    int iWbegin = iend   + 1;
    int iWend   = ibegin - 1;
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
    double *DL = (double*)malloc(isize*sizeof(double));
    assert(DL != NULL);
    
    double *DLL = (double*)malloc(isize*sizeof(double));
    assert(DLL != NULL);

    for (i=0; i<isize-1; i++) {
      double tmp = D[i+ibegin]*L[i+ibegin];
      DL[i]      = tmp;
      DLL[i]     = tmp*L[i+ibegin];
    }

    rrr_t *RRR = PMR_create_rrr(&D[ibegin], &L[ibegin], DL, DLL, isize, 0);
    PMR_increment_rrr_dependencies(RRR);
    
    /* In W apply shift of current block to eigenvalues
     * to get unshifted values w.r.t. T */
    for (i=ibegin; i<=iend; i++)
      W[i] += sigma;

    /* Split eigenvalues of block into singletons and clusters
     * and add them to process work queue */
    int max_size = imax(1, nz/nthreads);
    bool task_inserted = false;
    int new_first=ibegin, new_last;
    int sn_first, sn_last, sn_size;
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

      int new_size = new_last - new_first + 1;
      
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

          double lgap;
          if (sn_first == ibegin) {
            lgap = fmax(0.0, W[ibegin] - Werr[ibegin] - gl );
          } else {
            lgap = Wgap[sn_first-1];
          }

          PMR_increment_rrr_dependencies(RRR);

          task_t *task = 
            PMR_create_s_task
            (sn_first, sn_last, 1, ibegin, iend, spdiam, lgap, RRR);
          
          PMR_insert_task_at_back(workQ->s_queue, task);
          
          task_inserted = true;
        } else {
          task_inserted = false;
        }
      } else {
        /* Cluster was found */
        int cl_first = new_first;
        int cl_last  = new_last;
        int cl_size  = new_size;

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
            double lgap;
            if (cl_first == ibegin) {
              lgap = fmax(0.0, W[ibegin] - Werr[ibegin] - gl);
            } else {
              lgap = Wgap[cl_first-1];
            }

            /* Determine processes involved in processing the cluster */
            int left_pid  = nproc-1;
            int right_pid = 0;
            for (l=cl_first; l<=cl_last; l++) {
              if (iproc[l] != -1) {
                left_pid  = imin(left_pid,  iproc[l]);
                right_pid = imax(right_pid, iproc[l]);
              }
            }

            rrr_t *RRR_parent = 
              PMR_create_rrr(&D[ibegin], &L[ibegin], NULL, NULL, isize, 0);

            task_t *task = 
              PMR_create_c_task
              (cl_first, cl_last, 1, ibegin, iend, spdiam, lgap, iWbegin, 
               iWend, left_pid, right_pid, RRR_parent);

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
          double lgap;
          if (cl_first == ibegin) {
            lgap = fmax(0.0, W[ibegin] - Werr[ibegin] - gl );
          } else {
            lgap = Wgap[cl_first-1];
          }

          /* Determine processes involved */
          int left_pid  = nproc-1;
          int right_pid = 0;
          for (l=cl_first; l<=cl_last; l++) {
            if (iproc[l] != -1) {
              left_pid  = imin(left_pid,  iproc[l]);
              right_pid = imax(right_pid, iproc[l]);
            }
          }

          rrr_t *RRR_parent = 
            PMR_create_rrr
            (&D[ibegin], &L[ibegin], NULL, NULL, isize, 0);

          task_t *task = 
            PMR_create_c_task
            (cl_first, cl_last, 1, ibegin, 
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

  return 0;
}

/*
 * Processes all the tasks put in the work queue.
 */
static 
void *empty_workQ(void *argin)
{
  int        tid;
  proc_t    *procinfo;
  val_t     *Wstruct;
  vec_t     *Zstruct;
  tol_t     *tolstruct;
  workQ_t   *workQ;
  counter_t *num_left;
  retrieve_auxarg3
  ((auxarg3_t*)argin, &tid, &procinfo, &Wstruct,
   &Zstruct, &tolstruct, &workQ, &num_left);

  int n = Wstruct->n;

  /* max. needed double precision work space: odr1v */
  double *work = (double*)malloc(4*n*sizeof(double));
  assert(work != NULL);

  /* max. needed double precision work space: odrrb */
  int *iwork = (int*)malloc(2*n*sizeof(int));
  assert(iwork != NULL);

  /* while loop to empty the work queue */
  while (PMR_get_counter_value(num_left) > 0) {
    /* empty r-queue before processing other tasks */
    PMR_process_r_queue
    (tid, procinfo, Wstruct, Zstruct, tolstruct, workQ, num_left, work, iwork);

    task_t *task = PMR_remove_task_at_front(workQ->s_queue);
    if ( task != NULL ) {
      assert(task->flag == SINGLETON_TASK_FLAG);

      PMR_process_s_task
      ((singleton_t*)task->data, tid, procinfo,
       Wstruct, Zstruct, tolstruct, num_left, work, iwork);
      free(task);
      continue;
    }
    
    task = PMR_remove_task_at_front(workQ->c_queue);
    if ( task != NULL ) {
      assert(task->flag == CLUSTER_TASK_FLAG);

      PMR_process_c_task
      ((cluster_t*)task->data, tid, procinfo,
       Wstruct, Zstruct, tolstruct, workQ, num_left, work, iwork);
      free(task);
      continue;
    }
  } /* end while */

  free(work);
  free(iwork);

  return NULL;
}

static workQ_t *create_workQ()
{
  workQ_t *wq = (workQ_t*)malloc(sizeof(workQ_t));

  wq->r_queue = PMR_create_empty_queue();
  wq->s_queue = PMR_create_empty_queue();
  wq->c_queue = PMR_create_empty_queue();

  return wq;
}

static void destroy_workQ(workQ_t *wq)
{
  PMR_destroy_queue(wq->r_queue);
  PMR_destroy_queue(wq->s_queue);
  PMR_destroy_queue(wq->c_queue);
  free(wq);
}

static auxarg3_t*
create_auxarg3
(int tid, proc_t *procinfo, val_t *Wstruct, vec_t *Zstruct, 
 tol_t *tolstruct, workQ_t *workQ, counter_t *num_left)
{
  auxarg3_t *arg = (auxarg3_t*)malloc(sizeof(auxarg3_t));
  assert(arg != NULL);

  arg->tid       = tid;
  arg->procinfo  = procinfo;
  arg->Wstruct   = Wstruct;
  arg->Zstruct   = Zstruct;
  arg->tolstruct = tolstruct; 
  arg->workQ     = workQ;
  arg->num_left  = num_left;

  return arg;
}

static void 
retrieve_auxarg3
(auxarg3_t *arg, int *tid, proc_t **procinfo, 
 val_t **Wstruct, vec_t **Zstruct, tol_t **tolstruct, 
 workQ_t **workQ, counter_t **num_left)
{
  *tid       = arg->tid;
  *procinfo  = arg->procinfo;
  *Wstruct   = arg->Wstruct;
  *Zstruct   = arg->Zstruct;
  *tolstruct = arg->tolstruct;
  *workQ     = arg->workQ;
  *num_left  = arg->num_left;

  free(arg);
}
