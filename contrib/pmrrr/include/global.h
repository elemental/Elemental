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

#ifndef GGLOBAL_H
#define GGLOBAL_H

#ifdef __STDC__
 /* Some version of Standard C */
#if defined (__STDC_VERSION__) && __STDC_VERSION__>=199901L
 /* C99 */

#include <stdbool.h>
/* #pragma STDC FP_CONTRACT ON|OFF|DEFAULT */

#else
 /* C89 with or without Amendment 1, see H&S p.53 */

#define restrict /*nothing*/
#define inline   /*nothing*/

#define fmax(a,b) ( (a) > (b) ? (a) : (b) )
#define fmin(a,b) ( (a) < (b) ? (a) : (b) )

typedef int    bool;
#define false  0
#define true   1

#endif
#else  /* __STDC__ not defined */
 /* Not Standard C */

#define inline /*nothing*/

 /* in case compiler does not support type qualifiers
  * see Harbison and Steele p. 89*/
#define const /*nothing*/
#define volatile /*nothing*/
#define restrict /*nothing*/

#define fmax(a,b) ( (a) > (b) ? (a) : (b) )
#define fmin(a,b) ( (a) < (b) ? (a) : (b) )

typedef int    bool;
#define false  0
#define true   1
 
 /* need also remove all // style comments */
 /* inline ... */

#endif
/* see H&S p.70 to include C++ in the list above */

/* max/min/ceil for integers */
#define imax(a,b) ( (a) > (b) ? (a) : (b) )
#define imin(a,b) ( (a) < (b) ? (a) : (b) )
#define iceil(a,b) ( (((a)%(b))!=0) ? (((a)/(b))+1) : ((a)/(b)) )

/* End of header file */
#endif
