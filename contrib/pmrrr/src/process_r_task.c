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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <semaphore.h>
#include <assert.h>
#include "mpi.h"
#include "pmrrr.h"
#include "global.h"
#include "plarrv.h"
#include "queue.h"
#include "counter.h"
#include "structs.h"
#include "tasks.h"
#include "process_task.h"


int PMR_process_r_task(refine_t *rf, proc_t *procinfo, val_t *Wstruct,
		       tol_t *tolstruct, double *work, int *iwork);


/*
 * Executes all tasks which are in the r-queue at the moment of the 
 * call. This routine is called to make sure that all tasks in the 
 * queue are dequeued before continueing with other tasks.
 */
void PMR_process_r_queue(int tid, proc_t *procinfo, val_t *Wstruct, 
			 vec_t *Zstruct, tol_t *tolstruct, 
			 workQ_t *workQ, counter_t *num_left, 
			 double *work, int *iwork)
{
  int        thread_support = procinfo->thread_support;
  int        t, num_tasks;
  int        status;
  task_t     *task;

  num_tasks = PMR_get_num_tasks(workQ->r_queue);

  for (t=0; t<num_tasks; t++) {
    
    task = PMR_remove_task_at_front(workQ->r_queue);

    if ( task != NULL ) {
    
      if (task->flag == CLUSTER_TASK_FLAG) {

	if (thread_support != MPI_THREAD_FUNNELED || tid == 0) {
	  /* if MPI_THREAD_FUNNELED only tid==0 should process 
           * these tasks, otherwise any thread can do it */
	  status = PMR_process_c_task((cluster_t *) task->data,
				      tid, procinfo, Wstruct,
				      Zstruct, tolstruct, workQ,
				      num_left, work, iwork);
	  
	  if (status == C_TASK_PROCESSED) {
	    free(task);
	  } else {
	    PMR_insert_task_at_back(workQ->r_queue, task);
	  }
	} else {
	    PMR_insert_task_at_back(workQ->r_queue, task);
	}

      } /* end if cluster task */

      if (task->flag == REFINE_TASK_FLAG) {
	PMR_process_r_task((refine_t *) task->data, procinfo,
			   Wstruct, tolstruct, work, iwork);
	free(task);
      }
 
    } /* end if task removed */
  } /* end for t */
} /* end process_entire_r_queue */
  



/*
 * Process the task of refining a subset of eigenvalues.
 */
int PMR_process_r_task(refine_t *rf, proc_t *procinfo, 
		       val_t *Wstruct, tol_t *tolstruct, 
		       double *work, int *iwork)
{
  /* From inputs */
  int              ts_begin  = rf->begin;
  double *restrict D         = rf->D;
  double *restrict DLL       = rf->DLL;
  int              p         = rf->p;
  int              q         = rf->q;
  int              bl_size   = rf->bl_size;
  double           bl_spdiam = rf->bl_spdiam;
  sem_t            *sem      = rf->sem;

  double *restrict Werr      = Wstruct->Werr;
  double *restrict Wgap      = Wstruct->Wgap;
  int    *restrict Windex    = Wstruct->Windex;
  double *restrict Wshifted  = Wstruct->Wshifted;
  
  double           rtol1     = tolstruct->rtol1;
  double           rtol2     = tolstruct->rtol2;
  double           pivmin    = tolstruct->pivmin;

  /* Others */
  int    info, offset;
  double savegap;

  offset = Windex[ts_begin] - 1;

  if (p == q) {
    savegap = Wgap[ts_begin];
    Wgap[ts_begin] = 0.0;
  }  

  LAPACK(dlarrb)
  (&bl_size, D, DLL, &p, &q, &rtol1, &rtol2, &offset, &Wshifted[ts_begin], 
   &Wgap[ts_begin], &Werr[ts_begin], work, iwork, &pivmin, &bl_spdiam, 
   &bl_size, &info);
  assert(info == 0);

  if (p == q) {
    Wgap[ts_begin] = savegap;
  }  

  sem_post(sem);
  free(rf);

  return(0);
}
