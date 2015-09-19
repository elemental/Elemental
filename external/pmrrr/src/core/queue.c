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
#include "pmrrr/queue.h"

#ifndef DISABLE_PTHREADS
# include <errno.h>
#endif

int PMR_queue_init_lock(queue_t *queue)
{
#ifndef DISABLE_PTHREADS
 #ifdef NOSPINLOCKS
  int info = pthread_mutex_init(&queue->lock, NULL);
 #else
  int info = pthread_spin_init(&queue->lock, PTHREAD_PROCESS_PRIVATE);
 #endif
  assert(info == 0);
  return info;
#else
  return 0;
#endif
}

void PMR_queue_destroy_lock(queue_t *queue)
{
#ifndef DISABLE_PTHREADS
 #ifdef NOSPINLOCKS
  pthread_mutex_destroy(&queue->lock);
 #else
  pthread_spin_destroy(&queue->lock);
 #endif
#endif
}

int PMR_queue_lock(queue_t *queue)
{
#ifndef DISABLE_PTHREADS
 #ifdef NOSPINLOCKS
  int info = pthread_mutex_lock(&queue->lock);
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
  int info = pthread_spin_lock(&queue->lock);
 #endif
  assert(info == 0);
  return info;
#else
  return 0;
#endif
}

int PMR_queue_unlock(queue_t *queue)
{
#ifndef DISABLE_PTHREADS
 #ifdef NOSPINLOCKS
  int info = pthread_mutex_unlock(&queue->lock);
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
  int info = pthread_spin_unlock(&queue->lock);
 #endif
  assert(info == 0);
  return info;
#else
  return 0;
#endif
}

queue_t *PMR_create_empty_queue(void)
{
  queue_t* queue = (queue_t*)malloc(sizeof(queue_t)); assert(queue!=NULL);
  
  queue->num_tasks = 0;
  queue->head      = NULL;
  queue->back      = NULL;
  
  PMR_queue_init_lock(queue);
  return queue;
}

void PMR_destroy_queue(queue_t *queue)
{
  PMR_queue_destroy_lock(queue);
  free(queue);
}

int PMR_insert_task_at_front(queue_t *queue, task_t *task)
{
  int info = PMR_queue_lock(queue);

  queue->num_tasks++;
  task->next = queue->head;
  if (queue->head == NULL)
    queue->back = task;
  else
    queue->head->prev = task;
  queue->head = task;

  info |= PMR_queue_unlock(queue);
  return info;
}

int PMR_insert_task_at_back(queue_t *queue, task_t *task)
{
  int info = PMR_queue_lock(queue);

  queue->num_tasks++;
  task->prev = queue->back;
  task->next = NULL;
  if (queue->head == NULL)
    queue->head = task;
  else
    queue->back->next = task;
  queue->back = task;

  info |= PMR_queue_unlock(queue);
  return info;
}

/* returns NULL when empty */
task_t *PMR_remove_task_at_front(queue_t *queue)
{
  int info = PMR_queue_lock(queue);
  
  task_t *task = queue->head;
  if (queue->head != NULL) {
    /* at least one element */
    queue->num_tasks--;
    if (queue->head->next == NULL) {
      /* last task removed */
      queue->head = NULL;
      queue->back = NULL;
    } else {
      /* at least two tasks */
      queue->head->next->prev = NULL;
      queue->head = queue->head->next;
    }
  }
  
  info |= PMR_queue_unlock(queue);
  return task;
}

/* returns NULL when empty */
task_t *PMR_remove_task_at_back (queue_t *queue)
{
  int info = PMR_queue_lock(queue);
  
  task_t *task = queue->back;
  if (queue->back != NULL) {
    /* at least one element */
    queue->num_tasks--;
    if (queue->back->prev == NULL) {
      /* last task removed */
      queue->head = NULL;
      queue->back = NULL;
    } else {
      /* at least two tasks */
      queue->back->prev->next = NULL;
      queue->back = queue->back->prev;
    }
  }
  
  info |= PMR_queue_unlock(queue);
  return task;
}

int PMR_get_num_tasks(queue_t *queue)
{
  int info = PMR_queue_lock(queue);
  int num_tasks = queue->num_tasks;
  info |= PMR_queue_unlock(queue);
  return num_tasks;
}
