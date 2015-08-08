/* ========================================================================== */
/* === ldl.h:  include file for the LDL package ============================= */
/* ========================================================================== */

/* 
 * Copyright (c) Timothy A Davis, http://www.suitesparse.com.
 * All Rights Reserved. 
 *
 * Your use or distribution of LDL or any modified version of
 * LDL implies that you agree to this License.

 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301
 * USA
 *
 * Permission is hereby granted to use or copy this program under the
 * terms of the GNU LGPL, provided that the Copyright, this License,
 * and the Availability of the original version is retained on all copies.
 * User documentation of any code that uses this code or any modified
 * version of this code must cite the Copyright, this License, the
 * Availability note, and "Used by permission." Permission to modify
 * the code and to distribute modified code is granted, provided the
 * Copyright, this License, and the Availability note are retained,
 * and a notice that the code was modified is included.
 */
#include "SuiteSparse_config.h"

#ifdef LDL_LONG
#define LDL_int SuiteSparse_long
#define LDL_ID SuiteSparse_long_id

#define LDL_symbolic ldl_l_symbolic
#define LDL_numeric ldl_l_numeric
#define LDL_lsolve ldl_l_lsolve
#define LDL_dsolve ldl_l_dsolve
#define LDL_ltsolve ldl_l_ltsolve
#define LDL_perm ldl_l_perm
#define LDL_permt ldl_l_permt
#define LDL_valid_perm ldl_l_valid_perm
#define LDL_valid_matrix ldl_l_valid_matrix

#else
#define LDL_int int
#define LDL_ID "%d"

#define LDL_symbolic ldl_symbolic
#define LDL_numeric ldl_numeric
#define LDL_lsolve ldl_lsolve
#define LDL_dsolve ldl_dsolve
#define LDL_ltsolve ldl_ltsolve
#define LDL_perm ldl_perm
#define LDL_permt ldl_permt
#define LDL_valid_perm ldl_valid_perm
#define LDL_valid_matrix ldl_valid_matrix

#endif

/* ========================================================================== */
/* === int version ========================================================== */
/* ========================================================================== */

void ldl_symbolic (int n, int Ap [ ], int Ai [ ], int Lp [ ],
    int Parent [ ], int Lnz [ ], int Flag [ ], int P [ ], int Pinv [ ]) ;

int ldl_numeric (int n, int Ap [ ], int Ai [ ], double Ax [ ],
    int Lp [ ], int Parent [ ], int Lnz [ ], int Li [ ], double Lx [ ],
    double D [ ], double Y [ ], int Pattern [ ], int Flag [ ],
    int P [ ], int Pinv [ ]) ;

void ldl_lsolve (int n, double X [ ], int Lp [ ], int Li [ ],
    double Lx [ ]) ;

void ldl_dsolve (int n, double X [ ], double D [ ]) ;

void ldl_ltsolve (int n, double X [ ], int Lp [ ], int Li [ ],
    double Lx [ ]) ;

void ldl_perm  (int n, double X [ ], double B [ ], int P [ ]) ;
void ldl_permt (int n, double X [ ], double B [ ], int P [ ]) ;

int ldl_valid_perm (int n, int P [ ], int Flag [ ]) ;
int ldl_valid_matrix ( int n, int Ap [ ], int Ai [ ]) ;

/* ========================================================================== */
/* === long version ========================================================= */
/* ========================================================================== */

void ldl_l_symbolic (SuiteSparse_long n, SuiteSparse_long Ap [ ],
    SuiteSparse_long Ai [ ], SuiteSparse_long Lp [ ],
    SuiteSparse_long Parent [ ], SuiteSparse_long Lnz [ ],
    SuiteSparse_long Flag [ ], SuiteSparse_long P [ ],
    SuiteSparse_long Pinv [ ]) ;

SuiteSparse_long ldl_l_numeric (SuiteSparse_long n, SuiteSparse_long Ap [ ],
    SuiteSparse_long Ai [ ], double Ax [ ], SuiteSparse_long Lp [ ],
    SuiteSparse_long Parent [ ], SuiteSparse_long Lnz [ ],
    SuiteSparse_long Li [ ], double Lx [ ], double D [ ], double Y [ ],
    SuiteSparse_long Pattern [ ], SuiteSparse_long Flag [ ],
    SuiteSparse_long P [ ], SuiteSparse_long Pinv [ ]) ;

void ldl_l_lsolve (SuiteSparse_long n, double X [ ], SuiteSparse_long Lp [ ],
    SuiteSparse_long Li [ ], double Lx [ ]) ;

void ldl_l_dsolve (SuiteSparse_long n, double X [ ], double D [ ]) ;

void ldl_l_ltsolve (SuiteSparse_long n, double X [ ], SuiteSparse_long Lp [ ],
    SuiteSparse_long Li [ ], double Lx [ ]) ;

void ldl_l_perm  (SuiteSparse_long n, double X [ ], double B [ ],
    SuiteSparse_long P [ ]) ;
void ldl_l_permt (SuiteSparse_long n, double X [ ], double B [ ],
    SuiteSparse_long P [ ]) ;

SuiteSparse_long ldl_l_valid_perm (SuiteSparse_long n, SuiteSparse_long P [ ],
    SuiteSparse_long Flag [ ]) ;
SuiteSparse_long ldl_l_valid_matrix ( SuiteSparse_long n,
    SuiteSparse_long Ap [ ], SuiteSparse_long Ai [ ]) ;

/* ========================================================================== */
/* === LDL version ========================================================== */
/* ========================================================================== */

#define LDL_DATE "Oct 10, 2014"
#define LDL_VERSION_CODE(main,sub) ((main) * 1000 + (sub))
#define LDL_MAIN_VERSION 2
#define LDL_SUB_VERSION 2
#define LDL_SUBSUB_VERSION 1
#define LDL_VERSION LDL_VERSION_CODE(LDL_MAIN_VERSION,LDL_SUB_VERSION)

