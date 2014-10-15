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

#ifndef RRR_H
#define RRR_H

#include <pthread.h>
#include "global.h"

typedef struct {
  double          *restrict D;
  double          *restrict L;
  double          *restrict DL;
  double          *restrict DLL;
  int             size;
  int             depth;
  bool            parent_processed;
  bool            copied_parent_rrr;
  int             ndepend;
  pthread_mutex_t mutex;
} rrr_t;


rrr_t *PMR_create_rrr(double *restrict D, double *restrict L,
		      double *restrict DL, double *restrict DLL,
		      int size, int depth);

rrr_t *PMR_reset_rrr (rrr_t *restrict RRR, double *restrict D,
		      double *restrict L, double *restrict DL,
		      double *restrict DLL, int size, int depth);

int  PMR_increment_rrr_dependencies(rrr_t *RRR);
int  PMR_set_parent_processed_flag (rrr_t *RRR);
int  PMR_set_copied_parent_rrr_flag(rrr_t *RRR, bool val);
int  PMR_try_destroy_rrr(rrr_t *RRR);

#endif
