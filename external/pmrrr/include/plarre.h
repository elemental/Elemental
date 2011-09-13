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

#ifndef PPLARRE_H
#define PPLARRE_H

#include "global.h"
#include "structs.h"

/*
 * Computing the eigenvalues and of a symmetric tridiagonal matrix 
 * by bisection implemented in LAPACK routines 'dlarrd' and 'dlarrb'.
 *
 * Note: this implementation is here is not really optimized in
 * terms of performance and memory usage and could definitely be
 * optimized.
 */
int plarre(proc_t *procinfo, char *jobz, char *range, in_t *Dstruct, 
	   val_t *Wstruct, tol_t *tolstruct, int *nzp, int *myfirstp);

/* Perturb the initial root representation by "1 + eps*RAND_FACTOR*rand"; 
 * default: 8.0 */ 
#define RAND_FACTOR         8.0 

/* Fudge of shift, to get root rrr; default: 2.0 */
#define FUDGE_FACTOR        2.0 

#endif

