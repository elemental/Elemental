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

#ifndef QQUEUE_H
#define QQUEUE_H

#include <pthread.h>
#include "global.h"

typedef struct task_aux task_t;
struct task_aux {
  void       *data;       /* ptr to data, has to be casted */
  int         flag;       /* flag specifying the task */
  task_t     *next;       /* ptr to next  task; NULL if non-existing; */
  task_t     *prev;       /* ptr to prev. task; NULL if non-existing; */
};

typedef struct {
  int                num_tasks;
  task_t            *head;
  task_t            *back;
#ifdef NOSPINLOCKS
  pthread_mutex_t    lock;
#else
  pthread_spinlock_t lock;
#endif
} queue_t;


/* functionality of the queue */
queue_t *PMR_create_empty_queue  (void);
int     PMR_insert_task_at_front (queue_t *queue, task_t *task);
int     PMR_insert_task_at_back  (queue_t *queue, task_t *task);
task_t  *PMR_remove_task_at_front(queue_t *queue);
task_t  *PMR_remove_task_at_back (queue_t *queue);
int     PMR_get_num_tasks(queue_t *queue);
void    PMR_destroy_queue(queue_t *queue);

#endif
