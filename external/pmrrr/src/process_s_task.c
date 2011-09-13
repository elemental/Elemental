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
#include <math.h>
#include <float.h>
#include <assert.h>
#include "pmrrr.h"
#include "global.h"
#include "rrr.h"
#include "counter.h"
#include "structs.h"
#include "tasks.h"
#include "process_task.h"


int PMR_process_s_task(singleton_t *sng, int tid, proc_t *procinfo,
		       val_t *Wstruct, vec_t *Zstruct, 
		       tol_t *tolstruct, counter_t *num_left, 
		       double *work, int *iwork)
{
  /* Inputs */
  int    begin         = sng->begin; 
  int    end           = sng->end;
  int    bl_begin      = sng->bl_begin;
  int    bl_end        = sng->bl_end;
  int    bl_size       = bl_end - bl_begin + 1;
  double bl_spdiam     = sng->bl_spdiam; 
  rrr_t  *RRR          = sng->RRR;
  double *restrict D   = RRR->D; 
  double *restrict L   = RRR->L; 
  double *restrict DL  = RRR->DL;
  double *restrict DLL = RRR->DLL;

  int              pid      = procinfo->pid;  
  int              n        = Wstruct->n;
  double *restrict W        = Wstruct->W;
  double *restrict Werr     = Wstruct->Werr;
  double *restrict Wgap     = Wstruct->Wgap;
  int    *restrict Windex   = Wstruct->Windex;  
  int    *restrict iproc    = Wstruct->iproc;  
  double *restrict Wshifted = Wstruct->Wshifted;
  int              ldz      = Zstruct->ldz;
  double *restrict Z        = Zstruct->Z;
  int    *restrict isuppZ   = Zstruct->Zsupp;;
  int    *restrict Zindex   = Zstruct->Zindex;
  double           pivmin   = tolstruct->pivmin;

  /* others */
  int              info, i, k, itmp, num_decrement=0;
  int              IONE = 1;
  double           DZERO = 0.0;
  double           tol, lambda, left, right;
  int              i_local, zind;
  double           gap, lgap, rgap, gaptol, savedgap, tmp;
  bool             usedBS, usedRQ, needBS, wantNC, step2II;
  int              r, offset;
  double           twoeps = 2*DBL_EPSILON, RQtol = 2*DBL_EPSILON;
  double           residual, bstres, bstw; 
  int              i_supmn, i_supmx;
  double           RQcorr;
  int              negcount;
  int              sgndef, suppsize;
  double           sigma;
  int              i_Zfrom, i_Zto;
  double           ztz, norminv, mingma;


  /* set tolerance parameter */
  tol  = 4.0 * log( (double) bl_size ) * DBL_EPSILON;

  /* loop over all singletons in the task */
  for (i=begin; i<=end; i++) {

    /* check if eigenvector is supposed to be computed by
     * the process */
    if (iproc[i] != pid)
      continue;
    num_decrement++;

    if (bl_size == 1) {
      /* set eigenvector to column of identity matrix */
      zind = Zindex[i];
      memset(&Z[zind*ldz], 0.0, n*sizeof(double) );
      Z[zind*ldz + bl_begin] = 1.0;
      isuppZ[2*zind    ]     = bl_begin;
      isuppZ[2*zind + 1]     = bl_begin;
      continue;
    }

    lambda  = Wshifted[i];  
    left    = Wshifted[i] - Werr[i];
    right   = Wshifted[i] + Werr[i];
    i_local = Windex[i];
    r       = 0;
    
    /* compute left and right gap */
    if (i == bl_begin)
      lgap = DBL_EPSILON * fmax( fabs(left), fabs(right) );
    else if (i == begin)
      lgap = sng->lgap;
    else
      lgap = Wgap[i-1];

    if (i == bl_end) {
      rgap = DBL_EPSILON * fmax( fabs(left), fabs(right) );
    } else {
      rgap = Wgap[i];
    }

    gap = fmin(lgap, rgap);

    if ( i == bl_begin || i == bl_end ) {
      gaptol = 0.0;
    } else {
      gaptol = gap * DBL_EPSILON;
    }

    /* initialize lower and upper value of support */
    i_supmn = bl_size;
    i_supmx = 1;

    /* update Wgap so that it holds minimum gap and save the 
     * old value */
    savedgap  = Wgap[i];
    Wgap[i]   = gap;
    
    /* initialize flags indicating if bisection or Rayleigh-Quotient
     * correction was used */
    usedBS = false;
    usedRQ = false;
  
    /* the need for bisection is initially turned off */
    needBS = !TRY_RQC;

    /* IEEE floating point is assumed, so that all 0 bits are 0.0 */
    zind = Zindex[i];
    memset(&Z[zind*ldz], 0.0, n*sizeof(double));

    /* inverse iteration with twisted factorization */
    for (k=1; k<=MAXITER; k++) {

      if (needBS == true) {
	usedBS = true;
	itmp   = r;
	
	offset  = Windex[i] - 1;
	tmp     = Wgap[i]; 
	Wgap[i] = 0.0;
	
	LAPACK(dlarrb)
        (&bl_size, D, DLL, &i_local, &i_local, &DZERO, &twoeps, &offset, 
         &Wshifted[i], &Wgap[i], &Werr[i], work, iwork, &pivmin, &bl_spdiam,
         &itmp, &info);
	assert(info == 0);
	
	Wgap[i] = tmp;
	lambda = Wshifted[i];
	r = 0;
      }
      wantNC = (usedBS == true) ? false : true;

      /* compute the eigenvector corresponding to lambda */
      LAPACK(dlar1v)
      (&bl_size, &IONE, &bl_size, &lambda, D, L, DL, DLL, &pivmin, &gaptol, 
       &Z[zind*ldz+bl_begin], &wantNC, &negcount, &ztz, &mingma, &r, 
       &isuppZ[2*zind], &norminv, &residual, &RQcorr, work);

      if (k == 1) {
	bstres = residual;
	bstw   = lambda;
      } else if (residual < bstres) {
	bstres = residual;
	bstw   = lambda;
      }
      
      /* update support held */
      i_supmn = imin(i_supmn, isuppZ[2*zind    ]);
      i_supmx = imax(i_supmx, isuppZ[2*zind + 1]);

      /* Convergence test for Rayleigh Quotient Iteration
       * not done if bisection was used */
      if ( !usedBS && residual > tol*gap 
	   && fabs(RQcorr) > RQtol*fabs(lambda) ) {
      
	if (i_local <= negcount) {
	  sgndef = -1;    /* wanted eigenvalue lies to the left  */
	} else {
	  sgndef =  1;    /* wanted eigenvalue lies to the right */
	}
	
	if ( RQcorr*sgndef >= 0.0
	     && lambda+RQcorr <= right 
	     && lambda+RQcorr >= left ) {
	  usedRQ = true;
	  if ( sgndef == 1 )
	    left  = lambda;
	  else
	    right = lambda;
	  Wshifted[i] = 0.5*(left + right);
	  lambda     += RQcorr;
	} else { /* bisection is needed */
	  needBS = true;
	}
	
	if ( right-left < RQtol*fabs(lambda) ) {
	  /* eigenvalue computed to bisection accuracy
	   * => compute eigenvector */
	  usedBS = true;
	} else if ( k == MAXITER-1 ) {
	  /* for last iteration use bisection */
	  needBS = true;
	}
      } else {
	/* go to next iteration */
	break;
      }

    } /* end k */

    /* if necessary call dlar1v to improve error angle by 2nd step */
    step2II = false;
    if ( usedRQ && usedBS && (bstres <= residual) ) {
      lambda = bstw;
      step2II = true;
    }
    if ( step2II == true ) {
      LAPACK(dlar1v)
      (&bl_size, &IONE, &bl_size, &lambda, D, L, DL, DLL, &pivmin, &gaptol, 
       &Z[zind*ldz+bl_begin], &wantNC, &negcount, &ztz, &mingma, &r, 
       &isuppZ[2*zind], &norminv, &residual, &RQcorr, work);
    }
    Wshifted[i] = lambda;

    /* compute support w.r.t. whole matrix
     * block beginning is offset for each support */
    isuppZ[2*zind    ] += bl_begin;
    isuppZ[2*zind + 1] += bl_begin;
  
    /* ensure vector is okay if support changed in RQI 
     * minus ones because of indices starting from zero */
    i_Zfrom    = isuppZ[2*zind    ] - 1;
    i_Zto      = isuppZ[2*zind + 1] - 1;
    i_supmn   += bl_begin - 1;
    i_supmx   += bl_begin - 1;
    if ( i_supmn < i_Zfrom ) {
      for ( k=i_supmn; k < i_Zfrom; k++ ) {
	Z[k + zind*ldz] = 0.0;
      }
    }
    if ( i_supmx > i_Zto ) {
      for ( k=i_Zto+1; k <= i_supmx; k++ ) {
	Z[k + zind*ldz] = 0.0;
      }
    }
    
    /* normalize eigenvector */
    suppsize = i_Zto - i_Zfrom + 1;
    pmrrr_dscal(&suppsize, &norminv, &Z[i_Zfrom + zind*ldz], &IONE);

    sigma = L[bl_size-1];
    W[i]  = lambda + sigma;
    
    if (i < end)
      Wgap[i] = fmax(savedgap, W[i+1]-Werr[i+1] - W[i]-Werr[i]);

  } /* end i */

  /* decrement counter */
  PMR_decrement_counter(num_left, num_decrement);

  /* clean up */
  free(sng);
  PMR_try_destroy_rrr(RRR);

  return(0);
}
