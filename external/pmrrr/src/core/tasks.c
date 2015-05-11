/* Copyright (c) 2010, RWTH Aachen University
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
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "pmrrr/global.h"
#include "pmrrr/tasks.h"
#include "pmrrr/rrr.h"

task_t *PMR_create_s_task
(int first, int last, int depth, int bl_begin, int bl_end, 
 double spdiam, double lgap, rrr_t *RRR)
{
  task_t *t = (task_t*)malloc(sizeof(task_t)); assert(t!=NULL);
  singleton_t *s = (singleton_t*)malloc(sizeof(singleton_t)); assert(s!=NULL);

  s->begin     = first;
  s->end       = last;
  s->depth     = depth;  
  s->bl_begin  = bl_begin;
  s->bl_end    = bl_end;
  s->bl_spdiam = spdiam;
  s->lgap      = lgap;
  s->RRR       = RRR;

  t->data = (void*)s;
  t->flag = SINGLETON_TASK_FLAG;
  t->next = NULL;
  t->prev = NULL;

  return t;
}

task_t *PMR_create_c_task
(int first, int last, int depth, int bl_begin, int bl_end, double spdiam,
 double lgap, int proc_W_begin, int proc_W_end, int lpid, int rpid, rrr_t *RRR)
{
  task_t *t = (task_t*)malloc(sizeof(task_t)); assert(t!=NULL);
  cluster_t *c = (cluster_t*)malloc(sizeof(cluster_t)); assert(c!=NULL);

  c->begin              = first;
  c->end                = last;
  c->depth              = depth;  
  c->bl_begin           = bl_begin;
  c->bl_end             = bl_end;
  c->bl_spdiam          = spdiam;
  c->lgap               = lgap;
  c->proc_W_begin       = proc_W_begin;
  c->proc_W_end         = proc_W_end;
  c->left_pid           = lpid;
  c->right_pid          = rpid;
  c->RRR                = RRR;
  c->wait_until_refined = false;

  t->data = (void*)c;
  t->flag = CLUSTER_TASK_FLAG;
  t->next = NULL;
  t->prev = NULL;

  return t;
}

int PMR_refine_sem_init(refine_t *refine)
{
#ifndef DISABLE_PTHREADS
  int info = sem_init(refine->sem, 0, 0);
  assert(info == 0);
  return info;
#else
  return 0;
#endif
}

int PMR_refine_sem_destroy(refine_t *refine)
{
#ifndef DISABLE_PTHREADS
  int info = sem_destroy(refine->sem);
  assert(info == 0);
  return info;
#else
  return 0;
#endif
}

int PMR_refine_sem_wait(refine_t *refine)
{
#ifndef DISABLE_PTHREADS
  int info = sem_wait(refine->sem);
  assert(info == 0);
  return info;
#else
  return 0;
#endif
}

int PMR_refine_sem_post(refine_t *refine)
{
#ifndef DISABLE_PTHREADS
  int info = sem_post(refine->sem);
  assert(info == 0);
  return info;
#else
  return 0;
#endif
}

task_t *PMR_create_r_task
(int begin, int end, double *D, double *DLL, 
 int p, int q, int bl_size, double bl_spdiam, int tid)
{
  task_t *t = (task_t*)malloc(sizeof(task_t)); assert(t!=NULL);
  refine_t *r = (refine_t*)malloc(sizeof(refine_t)); assert(r!=NULL); 

  r->begin        = begin;
  r->end          = end;
  r->D            = D;
  r->DLL          = DLL;
  r->p            = p;
  r->q            = q;
  r->bl_size      = bl_size;
  r->bl_spdiam    = bl_spdiam;
  r->producer_tid = tid;

  PMR_refine_sem_init(r);

  t->data = (void*)r;
  t->flag = REFINE_TASK_FLAG;
  t->next = NULL;
  t->prev = NULL;

  return t;
}
