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
#include "pmrrr/counter.h"

#ifndef DISABLE_PTHREADS
# include <errno.h>
#endif

int PMR_counter_init_lock(counter_t *counter)
{
#ifndef DISABLE_PTHREADS
 #ifdef NOSPINLOCKS
  int info = pthread_mutex_init(&counter->lock, NULL);
 #else
  int info = pthread_spin_init(&counter->lock, PTHREAD_PROCESS_PRIVATE);
 #endif
  assert(info == 0);
  return info;
#else
  return 0;
#endif
}

void PMR_counter_destroy_lock(counter_t *counter)
{
#ifndef DISABLE_PTHREADS
 #ifdef NOSPINLOCKS
  pthread_mutex_destroy(&counter->lock);
 #else
  pthread_spin_destroy(&counter->lock);
 #endif
#endif
}

int PMR_counter_lock(counter_t *counter)
{
#ifndef DISABLE_PTHREADS
 #ifdef NOSPINLOCKS
  int info = pthread_mutex_lock(&counter->lock);
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
 #else
  int info = pthread_spin_lock(&counter->lock);
 #endif
  assert(info == 0);
 return info;
#else
  return 0;
#endif
}

int PMR_counter_unlock(counter_t *counter)
{
#ifndef DISABLE_PTHREADS
 #ifdef NOSPINLOCKS
  int info = pthread_mutex_unlock(&counter->lock);
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
 #else
  int info = pthread_spin_unlock(&counter->lock);
 #endif
  assert(info == 0);
  return info;
#else
  return 0;
#endif
}

counter_t *PMR_create_counter(int init_value)
{
  counter_t *counter = (counter_t *) malloc( sizeof(counter_t) );
  counter->value = init_value;
  int info = PMR_counter_init_lock(counter);
  return counter;
}

void PMR_destroy_counter(counter_t *counter)
{
  PMR_counter_destroy_lock(counter);
  free(counter);
}

int PMR_get_counter_value(counter_t *counter)
{
  int info = PMR_counter_lock(counter);
  int value = counter->value;
  info |= PMR_counter_unlock(counter);
  return value;
}

inline int PMR_set_counter_value(counter_t *counter, int value)
{
  int info = PMR_counter_lock(counter);
  counter->value = value;
  info |= PMR_counter_unlock(counter);
  return value;
}

int PMR_decrement_counter(counter_t *counter, int amount)
{
  int info = PMR_counter_lock(counter);
  counter->value -= amount;
  int value = counter->value;
  info |= PMR_counter_unlock(counter);
  return value;
}

int PMR_increment_counter(counter_t *counter, int amount)
{
  int info = PMR_counter_lock(counter);
  counter->value += amount;
  int value = counter->value;
  info |= PMR_counter_unlock(counter);
  return value;
}
