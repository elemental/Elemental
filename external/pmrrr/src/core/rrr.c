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
#include "pmrrr/rrr.h"
#include "pmrrr/global.h"

#ifndef DISABLE_PTHREADS
# include <errno.h>
#endif

int PMR_rrr_init_lock(rrr_t *RRR)
{
#ifndef DISABLE_PTHREADS
  int info = pthread_mutex_init(&RRR->mutex, NULL);
  assert(info == 0);
  return info;
#else
  return 0;
#endif
}

void PMR_rrr_destroy_lock(rrr_t *RRR)
{
#ifndef DISABLE_PTHREADS
  pthread_mutex_destroy(&RRR->mutex);
#endif
}

int PMR_rrr_lock(rrr_t *RRR)
{
#ifndef DISABLE_PTHREADS
  int info = pthread_mutex_lock(&RRR->mutex);
  if( info == EINVAL )
    fprintf(stderr,"pthread_mutex_lock returned EINVAL\n");
  else if( info == EAGAIN )
    fprintf(stderr,"pthread_mutex_lock returned EAGAIN\n");
  else if( info == EDEADLK )
    fprintf(stderr,"pthread_mutex_lock returned EDEADLK\n");
  else if( info == EPERM )
    fprintf(stderr,"pthread_mutex_lock returned EPERM\n");
  else
    fprintf(stderr,"pthread_mutex_lock returned %d\n",info);
  assert(info == 0);
  return info;
#else
  return 0;
#endif
}

int PMR_rrr_unlock(rrr_t *RRR)
{
#ifndef DISABLE_PTHREADS
  int info = pthread_mutex_unlock(&RRR->mutex);
  if( info == EINVAL )
    fprintf(stderr,"pthread_mutex_unlock returned EINVAL\n");
  else if( info == EAGAIN )
    fprintf(stderr,"pthread_mutex_unlock returned EAGAIN\n");
  else if( info == EDEADLK )
    fprintf(stderr,"pthread_mutex_unlock returned EDEADLK\n");
  else if( info == EPERM )
    fprintf(stderr,"pthread_mutex_unlock returned EPERM\n");
  else
    fprintf(stderr,"pthread_mutex_unlock returned %d\n",info);
  assert(info == 0);
  return info;
#else
  return 0;
#endif
}

rrr_t *PMR_create_rrr
(double *restrict D, double *restrict L,
 double *restrict DL, double *restrict DLL, int size, int depth)
{
  rrr_t* RRR = (rrr_t*)malloc(sizeof(rrr_t)); assert(RRR!=NULL);

  RRR->D                 = D;
  RRR->L                 = L;
  RRR->DL                = DL;
  RRR->DLL               = DLL;
  RRR->size              = size;
  RRR->depth             = depth;
  RRR->parent_processed  = false;
  RRR->copied_parent_rrr = false;
  RRR->ndepend           = 0;

  int info = PMR_rrr_init_lock(RRR);
  return RRR;
}
 
rrr_t *PMR_reset_rrr
(rrr_t *RRR, double *restrict D, 
 double *restrict L, double *restrict DL, 
 double *restrict DLL, int size, int depth)
{
  RRR->D                = D;
  RRR->L                = L;
  RRR->DL               = DL;
  RRR->DLL              = DLL;
  RRR->size             = size;
  RRR->depth            = depth;
  RRR->parent_processed = false;

  return RRR;
}

int PMR_increment_rrr_dependencies(rrr_t *RRR)
{
  int info = PMR_rrr_lock(RRR);
  RRR->ndepend++;
  int i = RRR->ndepend;
  info |= PMR_rrr_unlock(RRR);
  return i;
}

int PMR_set_parent_processed_flag(rrr_t *RRR)
{
  int info = PMR_rrr_lock(RRR);
  RRR->parent_processed = true;
  info |= PMR_rrr_unlock(RRR);
  return info;
}

int PMR_set_copied_parent_rrr_flag(rrr_t *RRR, bool val)
{
  int info = PMR_rrr_lock(RRR);
  RRR->copied_parent_rrr = val;
  info |= PMR_rrr_unlock(RRR);
  return info;
}

/* return 0 on success, otherwise 1 */
int PMR_try_destroy_rrr(rrr_t *RRR)
{
  int info = PMR_rrr_lock(RRR);

  RRR->ndepend--;
  int tmp=0;
  if (RRR->ndepend == 0 && RRR->parent_processed == true) {
    if (RRR->depth >0) {
      free(RRR->D);
      free(RRR->L);
    }
    if (RRR->depth >=0) {
      free(RRR->DL);
      free(RRR->DLL);
    }
    tmp = 1;
  }
  
  info |= PMR_rrr_unlock(RRR);

  if (tmp == 1) {
    PMR_rrr_destroy_lock(RRR);
    free(RRR);
    return 0;
  } else {
    return 1;
  }
}
