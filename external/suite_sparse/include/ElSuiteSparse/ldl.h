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
#include "ElSuiteSparse/config.h"

#ifdef EL_LDL_LONG
#define ElLDL_int ElSuiteSparse_long
#define ElLDL_ID ElSuiteSparse_long_id

#define ElLDL_symbolic El_ldl_l_symbolic
#define ElLDL_numeric El_ldl_l_numeric
#define ElLDL_lsolve El_ldl_l_lsolve
#define ElLDL_dsolve El_ldl_l_dsolve
#define ElLDL_ltsolve El_ldl_l_ltsolve
#define ElLDL_perm El_ldl_l_perm
#define ElLDL_permt El_ldl_l_permt
#define ElLDL_valid_perm El_ldl_l_valid_perm
#define ElLDL_valid_matrix El_ldl_l_valid_matrix

#else
#define ElLDL_int int
#define ElLDL_ID "%d"

#define ElLDL_symbolic El_ldl_symbolic
#define ElLDL_numeric El_ldl_numeric
#define ElLDL_lsolve El_ldl_lsolve
#define ElLDL_dsolve El_ldl_dsolve
#define ElLDL_ltsolve El_ldl_ltsolve
#define ElLDL_perm El_ldl_perm
#define ElLDL_permt El_ldl_permt
#define ElLDL_valid_perm El_ldl_valid_perm
#define ElLDL_valid_matrix El_ldl_valid_matrix

#endif

/* ========================================================================== */
/* === int version ========================================================== */
/* ========================================================================== */

void El_ldl_symbolic (int n, int Ap [ ], int Ai [ ], int Lp [ ],
    int Parent [ ], int Lnz [ ], int Flag [ ], int P [ ], int Pinv [ ]) ;

int El_ldl_numeric (int n, int Ap [ ], int Ai [ ], double Ax [ ],
    int Lp [ ], int Parent [ ], int Lnz [ ], int Li [ ], double Lx [ ],
    double D [ ], double Y [ ], int Pattern [ ], int Flag [ ],
    int P [ ], int Pinv [ ]) ;

void El_ldl_lsolve (int n, double X [ ], int Lp [ ], int Li [ ],
    double Lx [ ]) ;

void El_ldl_dsolve (int n, double X [ ], double D [ ]) ;

void El_ldl_ltsolve (int n, double X [ ], int Lp [ ], int Li [ ],
    double Lx [ ]) ;

void El_ldl_perm  (int n, double X [ ], double B [ ], int P [ ]) ;
void El_ldl_permt (int n, double X [ ], double B [ ], int P [ ]) ;

int El_ldl_valid_perm (int n, int P [ ], int Flag [ ]) ;
int El_ldl_valid_matrix ( int n, int Ap [ ], int Ai [ ]) ;

/* ========================================================================== */
/* === long version ========================================================= */
/* ========================================================================== */

void El_ldl_l_symbolic (ElSuiteSparse_long n, ElSuiteSparse_long Ap [ ],
    ElSuiteSparse_long Ai [ ], ElSuiteSparse_long Lp [ ],
    ElSuiteSparse_long Parent [ ], ElSuiteSparse_long Lnz [ ],
    ElSuiteSparse_long Flag [ ], ElSuiteSparse_long P [ ],
    ElSuiteSparse_long Pinv [ ]) ;

ElSuiteSparse_long ldl_l_numeric (
    ElSuiteSparse_long n, ElSuiteSparse_long Ap [ ],
    ElSuiteSparse_long Ai [ ], double Ax [ ], ElSuiteSparse_long Lp [ ],
    ElSuiteSparse_long Parent [ ], ElSuiteSparse_long Lnz [ ],
    ElSuiteSparse_long Li [ ], double Lx [ ], double D [ ], double Y [ ],
    ElSuiteSparse_long Pattern [ ], ElSuiteSparse_long Flag [ ],
    ElSuiteSparse_long P [ ], ElSuiteSparse_long Pinv [ ]) ;

void El_ldl_l_lsolve (
    ElSuiteSparse_long n, double X [ ], ElSuiteSparse_long Lp [ ],
    ElSuiteSparse_long Li [ ], double Lx [ ]) ;

void El_ldl_l_dsolve (ElSuiteSparse_long n, double X [ ], double D [ ]) ;

void El_ldl_l_ltsolve (
    ElSuiteSparse_long n, double X [ ], ElSuiteSparse_long Lp [ ],
    ElSuiteSparse_long Li [ ], double Lx [ ]) ;

void El_ldl_l_perm  (ElSuiteSparse_long n, double X [ ], double B [ ],
    ElSuiteSparse_long P [ ]) ;
void El_ldl_l_permt (ElSuiteSparse_long n, double X [ ], double B [ ],
    ElSuiteSparse_long P [ ]) ;

ElSuiteSparse_long ldl_l_valid_perm (
    ElSuiteSparse_long n, ElSuiteSparse_long P [ ],
    ElSuiteSparse_long Flag [ ]) ;
ElSuiteSparse_long ldl_l_valid_matrix (
    ElSuiteSparse_long n,
    ElSuiteSparse_long Ap [ ], ElSuiteSparse_long Ai [ ]) ;

/* ========================================================================== */
/* === LDL version ========================================================== */
/* ========================================================================== */

#define EL_LDL_DATE "Oct 10, 2014"
#define EL_LDL_VERSION_CODE(main,sub) ((main) * 1000 + (sub))
#define EL_LDL_MAIN_VERSION 2
#define EL_LDL_SUB_VERSION 2
#define EL_LDL_SUBSUB_VERSION 1
#define EL_LDL_VERSION EL_LDL_VERSION_CODE( \
 EL_LDL_MAIN_VERSION,EL_LDL_SUB_VERSION)

