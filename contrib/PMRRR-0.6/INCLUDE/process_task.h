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

#ifndef PPROCESS_TASK_H
#define PPROCESS_TASK_H

#include "global.h"
#include "structs.h"
#include "tasks.h"
#include "counter.h"
#include "queue.h"

/* Function prototypes */
int PMR_process_c_task(cluster_t *cl, int tid, proc_t *procinfo,
		       val_t *Wstruct, vec_t *Zstruct, 
		       tol_t *tolstruct, workQ_t *workQ, 
		       counter_t *num_left, double *work, int *iwork);

int PMR_process_s_task(singleton_t *sng, int tid, proc_t *procinfo,
		       val_t *Wstruct, vec_t *Zstruct, 
		       tol_t *tolstruct, counter_t *num_left, 
		       double *work, int *iwork);

int PMR_process_r_task(refine_t *rf, proc_t *procinfo, val_t *Wstruct,
		       tol_t *tolstruct, double *work, int *iwork);

void PMR_process_r_queue(int tid, proc_t *procinfo, val_t *Wstruct,
			 vec_t *Zstruct, tol_t *tolstruct,
			 workQ_t *workQ, 
			 counter_t *num_left, double *work, 
			 int *iwork);

#endif /* End of header file */
