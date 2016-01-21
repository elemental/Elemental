/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_FLAMEPART_SLIDEPARTITION_C_H
#define EL_FLAMEPART_SLIDEPARTITION_C_H

#ifdef __cplusplus
extern "C" {
#endif

/* Downward
   ======== */
EL_EXPORT ElError ElSlidePartitionDown_i
( ElMatrix_i AT, ElMatrix_i A0,
                 ElMatrix_i A1,
  ElMatrix_i AB, ElMatrix_i A2 );
EL_EXPORT ElError ElSlidePartitionDown_s
( ElMatrix_s AT, ElMatrix_s A0,
                 ElMatrix_s A1,
  ElMatrix_s AB, ElMatrix_s A2 );
EL_EXPORT ElError ElSlidePartitionDown_d
( ElMatrix_d AT, ElMatrix_d A0,
                 ElMatrix_d A1,
  ElMatrix_d AB, ElMatrix_d A2 );
EL_EXPORT ElError ElSlidePartitionDown_c
( ElMatrix_c AT, ElMatrix_c A0,
                 ElMatrix_c A1,
  ElMatrix_c AB, ElMatrix_c A2 );
EL_EXPORT ElError ElSlidePartitionDown_z
( ElMatrix_z AT, ElMatrix_z A0,
                 ElMatrix_z A1,
  ElMatrix_z AB, ElMatrix_z A2 );
EL_EXPORT ElError ElSlidePartitionDownDist_i
( ElDistMatrix_i AT, ElDistMatrix_i A0,
                     ElDistMatrix_i A1,
  ElDistMatrix_i AB, ElDistMatrix_i A2 );
EL_EXPORT ElError ElSlidePartitionDownDist_s
( ElDistMatrix_s AT, ElDistMatrix_s A0,
                     ElDistMatrix_s A1,
  ElDistMatrix_s AB, ElDistMatrix_s A2 );
EL_EXPORT ElError ElSlidePartitionDownDist_d
( ElDistMatrix_d AT, ElDistMatrix_d A0,
                     ElDistMatrix_d A1,
  ElDistMatrix_d AB, ElDistMatrix_d A2 );
EL_EXPORT ElError ElSlidePartitionDownDist_c
( ElDistMatrix_c AT, ElDistMatrix_c A0,
                     ElDistMatrix_c A1,
  ElDistMatrix_c AB, ElDistMatrix_c A2 );
EL_EXPORT ElError ElSlidePartitionDownDist_z
( ElDistMatrix_z AT, ElDistMatrix_z A0,
                     ElDistMatrix_z A1,
  ElDistMatrix_z AB, ElDistMatrix_z A2 );

EL_EXPORT ElError ElSlideLockedPartitionDown_i
( ElMatrix_i AT, ElConstMatrix_i A0,
                 ElConstMatrix_i A1,
  ElMatrix_i AB, ElConstMatrix_i A2 );
EL_EXPORT ElError ElSlideLockedPartitionDown_s
( ElMatrix_s AT, ElConstMatrix_s A0,
                 ElConstMatrix_s A1,
  ElMatrix_s AB, ElConstMatrix_s A2 );
EL_EXPORT ElError ElSlideLockedPartitionDown_d
( ElMatrix_d AT, ElConstMatrix_d A0,
                 ElConstMatrix_d A1,
  ElMatrix_d AB, ElConstMatrix_d A2 );
EL_EXPORT ElError ElSlideLockedPartitionDown_c
( ElMatrix_c AT, ElConstMatrix_c A0,
                 ElConstMatrix_c A1,
  ElMatrix_c AB, ElConstMatrix_c A2 );
EL_EXPORT ElError ElSlideLockedPartitionDown_z
( ElMatrix_z AT, ElConstMatrix_z A0,
                 ElConstMatrix_z A1,
  ElMatrix_z AB, ElConstMatrix_z A2 );
EL_EXPORT ElError ElSlideLockedPartitionDownDist_i
( ElDistMatrix_i AT, ElConstDistMatrix_i A0,
                     ElConstDistMatrix_i A1,
  ElDistMatrix_i AB, ElConstDistMatrix_i A2 );
EL_EXPORT ElError ElSlideLockedPartitionDownDist_s
( ElDistMatrix_s AT, ElConstDistMatrix_s A0,
                     ElConstDistMatrix_s A1,
  ElDistMatrix_s AB, ElConstDistMatrix_s A2 );
EL_EXPORT ElError ElSlideLockedPartitionDownDist_d
( ElDistMatrix_d AT, ElConstDistMatrix_d A0,
                     ElConstDistMatrix_d A1,
  ElDistMatrix_d AB, ElConstDistMatrix_d A2 );
EL_EXPORT ElError ElSlideLockedPartitionDownDist_c
( ElDistMatrix_c AT, ElConstDistMatrix_c A0,
                     ElConstDistMatrix_c A1,
  ElDistMatrix_c AB, ElConstDistMatrix_c A2 );
EL_EXPORT ElError ElSlideLockedPartitionDownDist_z
( ElDistMatrix_z AT, ElConstDistMatrix_z A0,
                     ElConstDistMatrix_z A1,
  ElDistMatrix_z AB, ElConstDistMatrix_z A2 );

/* Upward
   ====== */
EL_EXPORT ElError ElSlidePartitionUp_i
( ElMatrix_i AT, ElMatrix_i A0,
                 ElMatrix_i A1,
  ElMatrix_i AB, ElMatrix_i A2 );
EL_EXPORT ElError ElSlidePartitionUp_s
( ElMatrix_s AT, ElMatrix_s A0,
                 ElMatrix_s A1,
  ElMatrix_s AB, ElMatrix_s A2 );
EL_EXPORT ElError ElSlidePartitionUp_d
( ElMatrix_d AT, ElMatrix_d A0,
                 ElMatrix_d A1,
  ElMatrix_d AB, ElMatrix_d A2 );
EL_EXPORT ElError ElSlidePartitionUp_c
( ElMatrix_c AT, ElMatrix_c A0,
                 ElMatrix_c A1,
  ElMatrix_c AB, ElMatrix_c A2 );
EL_EXPORT ElError ElSlidePartitionUp_z
( ElMatrix_z AT, ElMatrix_z A0,
                 ElMatrix_z A1,
  ElMatrix_z AB, ElMatrix_z A2 );
EL_EXPORT ElError ElSlidePartitionUpDist_i
( ElDistMatrix_i AT, ElDistMatrix_i A0,
                     ElDistMatrix_i A1,
  ElDistMatrix_i AB, ElDistMatrix_i A2 );
EL_EXPORT ElError ElSlidePartitionUpDist_s
( ElDistMatrix_s AT, ElDistMatrix_s A0,
                     ElDistMatrix_s A1,
  ElDistMatrix_s AB, ElDistMatrix_s A2 );
EL_EXPORT ElError ElSlidePartitionUpDist_d
( ElDistMatrix_d AT, ElDistMatrix_d A0,
                     ElDistMatrix_d A1,
  ElDistMatrix_d AB, ElDistMatrix_d A2 );
EL_EXPORT ElError ElSlidePartitionUpDist_c
( ElDistMatrix_c AT, ElDistMatrix_c A0,
                     ElDistMatrix_c A1,
  ElDistMatrix_c AB, ElDistMatrix_c A2 );
EL_EXPORT ElError ElSlidePartitionUpDist_z
( ElDistMatrix_z AT, ElDistMatrix_z A0,
                     ElDistMatrix_z A1,
  ElDistMatrix_z AB, ElDistMatrix_z A2 );

EL_EXPORT ElError ElSlideLockedPartitionUp_i
( ElMatrix_i AT, ElConstMatrix_i A0,
                 ElConstMatrix_i A1,
  ElMatrix_i AB, ElConstMatrix_i A2 );
EL_EXPORT ElError ElSlideLockedPartitionUp_s
( ElMatrix_s AT, ElConstMatrix_s A0,
                 ElConstMatrix_s A1,
  ElMatrix_s AB, ElConstMatrix_s A2 );
EL_EXPORT ElError ElSlideLockedPartitionUp_d
( ElMatrix_d AT, ElConstMatrix_d A0,
                 ElConstMatrix_d A1,
  ElMatrix_d AB, ElConstMatrix_d A2 );
EL_EXPORT ElError ElSlideLockedPartitionUp_c
( ElMatrix_c AT, ElConstMatrix_c A0,
                 ElConstMatrix_c A1,
  ElMatrix_c AB, ElConstMatrix_c A2 );
EL_EXPORT ElError ElSlideLockedPartitionUp_z
( ElMatrix_z AT, ElConstMatrix_z A0,
                 ElConstMatrix_z A1,
  ElMatrix_z AB, ElConstMatrix_z A2 );
EL_EXPORT ElError ElSlideLockedPartitionUpDist_i
( ElDistMatrix_i AT, ElConstDistMatrix_i A0,
                     ElConstDistMatrix_i A1,
  ElDistMatrix_i AB, ElConstDistMatrix_i A2 );
EL_EXPORT ElError ElSlideLockedPartitionUpDist_s
( ElDistMatrix_s AT, ElConstDistMatrix_s A0,
                     ElConstDistMatrix_s A1,
  ElDistMatrix_s AB, ElConstDistMatrix_s A2 );
EL_EXPORT ElError ElSlideLockedPartitionUpDist_d
( ElDistMatrix_d AT, ElConstDistMatrix_d A0,
                     ElConstDistMatrix_d A1,
  ElDistMatrix_d AB, ElConstDistMatrix_d A2 );
EL_EXPORT ElError ElSlideLockedPartitionUpDist_c
( ElDistMatrix_c AT, ElConstDistMatrix_c A0,
                     ElConstDistMatrix_c A1,
  ElDistMatrix_c AB, ElConstDistMatrix_c A2 );
EL_EXPORT ElError ElSlideLockedPartitionUpDist_z
( ElDistMatrix_z AT, ElConstDistMatrix_z A0,
                     ElConstDistMatrix_z A1,
  ElDistMatrix_z AB, ElConstDistMatrix_z A2 );

/* Rightward
   ========= */
EL_EXPORT ElError ElSlidePartitionRight_i
( ElMatrix_i AL, ElMatrix_i AR, 
  ElMatrix_i A0, ElMatrix_i A1, ElMatrix_i A2 );
EL_EXPORT ElError ElSlidePartitionRight_s
( ElMatrix_s AL, ElMatrix_s AR, 
  ElMatrix_s A0, ElMatrix_s A1, ElMatrix_s A2 );
EL_EXPORT ElError ElSlidePartitionRight_d
( ElMatrix_d AL, ElMatrix_d AR, 
  ElMatrix_d A0, ElMatrix_d A1, ElMatrix_d A2 );
EL_EXPORT ElError ElSlidePartitionRight_c
( ElMatrix_c AL, ElMatrix_c AR, 
  ElMatrix_c A0, ElMatrix_c A1, ElMatrix_c A2 );
EL_EXPORT ElError ElSlidePartitionRight_z
( ElMatrix_z AL, ElMatrix_z AR, 
  ElMatrix_z A0, ElMatrix_z A1, ElMatrix_z A2 );
EL_EXPORT ElError ElSlidePartitionRightDist_i
( ElDistMatrix_i AL, ElDistMatrix_i AR, 
  ElDistMatrix_i A0, ElDistMatrix_i A1, ElDistMatrix_i A2 );
EL_EXPORT ElError ElSlidePartitionRightDist_s
( ElDistMatrix_s AL, ElDistMatrix_s AR, 
  ElDistMatrix_s A0, ElDistMatrix_s A1, ElDistMatrix_s A2 );
EL_EXPORT ElError ElSlidePartitionRightDist_d
( ElDistMatrix_d AL, ElDistMatrix_d AR, 
  ElDistMatrix_d A0, ElDistMatrix_d A1, ElDistMatrix_d A2 );
EL_EXPORT ElError ElSlidePartitionRightDist_c
( ElDistMatrix_c AL, ElDistMatrix_c AR, 
  ElDistMatrix_c A0, ElDistMatrix_c A1, ElDistMatrix_c A2 );
EL_EXPORT ElError ElSlidePartitionRightDist_z
( ElDistMatrix_z AL, ElDistMatrix_z AR, 
  ElDistMatrix_z A0, ElDistMatrix_z A1, ElDistMatrix_z A2 );

EL_EXPORT ElError ElSlideLockedPartitionRight_i
( ElMatrix_i AL, ElMatrix_i AR, 
  ElConstMatrix_i A0, ElConstMatrix_i A1, ElConstMatrix_i A2 );
EL_EXPORT ElError ElSlideLockedPartitionRight_s
( ElMatrix_s AL, ElMatrix_s AR, 
  ElConstMatrix_s A0, ElConstMatrix_s A1, ElConstMatrix_s A2 );
EL_EXPORT ElError ElSlideLockedPartitionRight_d
( ElMatrix_d AL, ElMatrix_d AR, 
  ElConstMatrix_d A0, ElConstMatrix_d A1, ElConstMatrix_d A2 );
EL_EXPORT ElError ElSlideLockedPartitionRight_c
( ElMatrix_c AL, ElMatrix_c AR, 
  ElConstMatrix_c A0, ElConstMatrix_c A1, ElConstMatrix_c A2 );
EL_EXPORT ElError ElSlideLockedPartitionRight_z
( ElMatrix_z AL, ElMatrix_z AR, 
  ElConstMatrix_z A0, ElConstMatrix_z A1, ElConstMatrix_z A2 );
EL_EXPORT ElError ElSlideLockedPartitionRightDist_i
( ElDistMatrix_i AL, ElDistMatrix_i AR, 
  ElConstDistMatrix_i A0, ElConstDistMatrix_i A1, ElConstDistMatrix_i A2 );
EL_EXPORT ElError ElSlideLockedPartitionRightDist_s
( ElDistMatrix_s AL, ElDistMatrix_s AR, 
  ElConstDistMatrix_s A0, ElConstDistMatrix_s A1, ElConstDistMatrix_s A2 );
EL_EXPORT ElError ElSlideLockedPartitionRightDist_d
( ElDistMatrix_d AL, ElDistMatrix_d AR, 
  ElConstDistMatrix_d A0, ElConstDistMatrix_d A1, ElConstDistMatrix_d A2 );
EL_EXPORT ElError ElSlideLockedPartitionRightDist_c
( ElDistMatrix_c AL, ElDistMatrix_c AR, 
  ElConstDistMatrix_c A0, ElConstDistMatrix_c A1, ElConstDistMatrix_c A2 );
EL_EXPORT ElError ElSlideLockedPartitionRightDist_z
( ElDistMatrix_z AL, ElDistMatrix_z AR, 
  ElConstDistMatrix_z A0, ElConstDistMatrix_z A1, ElConstDistMatrix_z A2 );

/* Leftward
   ======== */
EL_EXPORT ElError ElSlidePartitionLeft_i
( ElMatrix_i AL, ElMatrix_i AR, 
  ElMatrix_i A0, ElMatrix_i A1, ElMatrix_i A2 );
EL_EXPORT ElError ElSlidePartitionLeft_s
( ElMatrix_s AL, ElMatrix_s AR, 
  ElMatrix_s A0, ElMatrix_s A1, ElMatrix_s A2 );
EL_EXPORT ElError ElSlidePartitionLeft_d
( ElMatrix_d AL, ElMatrix_d AR, 
  ElMatrix_d A0, ElMatrix_d A1, ElMatrix_d A2 );
EL_EXPORT ElError ElSlidePartitionLeft_c
( ElMatrix_c AL, ElMatrix_c AR, 
  ElMatrix_c A0, ElMatrix_c A1, ElMatrix_c A2 );
EL_EXPORT ElError ElSlidePartitionLeft_z
( ElMatrix_z AL, ElMatrix_z AR, 
  ElMatrix_z A0, ElMatrix_z A1, ElMatrix_z A2 );
EL_EXPORT ElError ElSlidePartitionLeftDist_i
( ElDistMatrix_i AL, ElDistMatrix_i AR, 
  ElDistMatrix_i A0, ElDistMatrix_i A1, ElDistMatrix_i A2 );
EL_EXPORT ElError ElSlidePartitionLeftDist_s
( ElDistMatrix_s AL, ElDistMatrix_s AR, 
  ElDistMatrix_s A0, ElDistMatrix_s A1, ElDistMatrix_s A2 );
EL_EXPORT ElError ElSlidePartitionLeftDist_d
( ElDistMatrix_d AL, ElDistMatrix_d AR, 
  ElDistMatrix_d A0, ElDistMatrix_d A1, ElDistMatrix_d A2 );
EL_EXPORT ElError ElSlidePartitionLeftDist_c
( ElDistMatrix_c AL, ElDistMatrix_c AR, 
  ElDistMatrix_c A0, ElDistMatrix_c A1, ElDistMatrix_c A2 );
EL_EXPORT ElError ElSlidePartitionLeftDist_z
( ElDistMatrix_z AL, ElDistMatrix_z AR, 
  ElDistMatrix_z A0, ElDistMatrix_z A1, ElDistMatrix_z A2 );

EL_EXPORT ElError ElSlideLockedPartitionLeft_i
( ElMatrix_i AL, ElMatrix_i AR, 
  ElConstMatrix_i A0, ElConstMatrix_i A1, ElConstMatrix_i A2 );
EL_EXPORT ElError ElSlideLockedPartitionLeft_s
( ElMatrix_s AL, ElMatrix_s AR, 
  ElConstMatrix_s A0, ElConstMatrix_s A1, ElConstMatrix_s A2 );
EL_EXPORT ElError ElSlideLockedPartitionLeft_d
( ElMatrix_d AL, ElMatrix_d AR, 
  ElConstMatrix_d A0, ElConstMatrix_d A1, ElConstMatrix_d A2 );
EL_EXPORT ElError ElSlideLockedPartitionLeft_c
( ElMatrix_c AL, ElMatrix_c AR, 
  ElConstMatrix_c A0, ElConstMatrix_c A1, ElConstMatrix_c A2 );
EL_EXPORT ElError ElSlideLockedPartitionLeft_z
( ElMatrix_z AL, ElMatrix_z AR, 
  ElConstMatrix_z A0, ElConstMatrix_z A1, ElConstMatrix_z A2 );
EL_EXPORT ElError ElSlideLockedPartitionLeftDist_i
( ElDistMatrix_i AL, ElDistMatrix_i AR, 
  ElConstDistMatrix_i A0, ElConstDistMatrix_i A1, ElConstDistMatrix_i A2 );
EL_EXPORT ElError ElSlideLockedPartitionLeftDist_s
( ElDistMatrix_s AL, ElDistMatrix_s AR, 
  ElConstDistMatrix_s A0, ElConstDistMatrix_s A1, ElConstDistMatrix_s A2 );
EL_EXPORT ElError ElSlideLockedPartitionLeftDist_d
( ElDistMatrix_d AL, ElDistMatrix_d AR, 
  ElConstDistMatrix_d A0, ElConstDistMatrix_d A1, ElConstDistMatrix_d A2 );
EL_EXPORT ElError ElSlideLockedPartitionLeftDist_c
( ElDistMatrix_c AL, ElDistMatrix_c AR, 
  ElConstDistMatrix_c A0, ElConstDistMatrix_c A1, ElConstDistMatrix_c A2 );
EL_EXPORT ElError ElSlideLockedPartitionLeftDist_z
( ElDistMatrix_z AL, ElDistMatrix_z AR, 
  ElConstDistMatrix_z A0, ElConstDistMatrix_z A1, ElConstDistMatrix_z A2 );

/* Down a diagonal
   =============== */
EL_EXPORT ElError ElSlidePartitionDownDiagonal_i
( ElMatrix_i ATL, ElMatrix_i ATR, ElMatrix_i A00, ElMatrix_i A01, ElMatrix_i A02,
                                  ElMatrix_i A10, ElMatrix_i A11, ElMatrix_i A12,
  ElMatrix_i ABL, ElMatrix_i ABR, ElMatrix_i A20, ElMatrix_i A21, ElMatrix_i A22 );
EL_EXPORT ElError ElSlidePartitionDownDiagonal_s
( ElMatrix_s ATL, ElMatrix_s ATR, ElMatrix_s A00, ElMatrix_s A01, ElMatrix_s A02,
                                  ElMatrix_s A10, ElMatrix_s A11, ElMatrix_s A12,
  ElMatrix_s ABL, ElMatrix_s ABR, ElMatrix_s A20, ElMatrix_s A21, ElMatrix_s A22 );
EL_EXPORT ElError ElSlidePartitionDownDiagonal_d
( ElMatrix_d ATL, ElMatrix_d ATR, ElMatrix_d A00, ElMatrix_d A01, ElMatrix_d A02,
                                  ElMatrix_d A10, ElMatrix_d A11, ElMatrix_d A12,
  ElMatrix_d ABL, ElMatrix_d ABR, ElMatrix_d A20, ElMatrix_d A21, ElMatrix_d A22 );
EL_EXPORT ElError ElSlidePartitionDownDiagonal_c
( ElMatrix_c ATL, ElMatrix_c ATR, ElMatrix_c A00, ElMatrix_c A01, ElMatrix_c A02,
                                  ElMatrix_c A10, ElMatrix_c A11, ElMatrix_c A12,
  ElMatrix_c ABL, ElMatrix_c ABR, ElMatrix_c A20, ElMatrix_c A21, ElMatrix_c A22 );
EL_EXPORT ElError ElSlidePartitionDownDiagonal_z
( ElMatrix_z ATL, ElMatrix_z ATR, ElMatrix_z A00, ElMatrix_z A01, ElMatrix_z A02,
                                  ElMatrix_z A10, ElMatrix_z A11, ElMatrix_z A12,
  ElMatrix_z ABL, ElMatrix_z ABR, ElMatrix_z A20, ElMatrix_z A21, ElMatrix_z A22 );
EL_EXPORT ElError ElSlidePartitionDownDiagonalDist_i
( ElDistMatrix_i ATL, ElDistMatrix_i ATR, ElDistMatrix_i A00, ElDistMatrix_i A01, ElDistMatrix_i A02,
                                          ElDistMatrix_i A10, ElDistMatrix_i A11, ElDistMatrix_i A12,
  ElDistMatrix_i ABL, ElDistMatrix_i ABR, ElDistMatrix_i A20, ElDistMatrix_i A21, ElDistMatrix_i A22 );
EL_EXPORT ElError ElSlidePartitionDownDiagonalDist_s
( ElDistMatrix_s ATL, ElDistMatrix_s ATR, ElDistMatrix_s A00, ElDistMatrix_s A01, ElDistMatrix_s A02,
                                          ElDistMatrix_s A10, ElDistMatrix_s A11, ElDistMatrix_s A12,
  ElDistMatrix_s ABL, ElDistMatrix_s ABR, ElDistMatrix_s A20, ElDistMatrix_s A21, ElDistMatrix_s A22 );
EL_EXPORT ElError ElSlidePartitionDownDiagonalDist_d
( ElDistMatrix_d ATL, ElDistMatrix_d ATR, ElDistMatrix_d A00, ElDistMatrix_d A01, ElDistMatrix_d A02,
                                          ElDistMatrix_d A10, ElDistMatrix_d A11, ElDistMatrix_d A12,
  ElDistMatrix_d ABL, ElDistMatrix_d ABR, ElDistMatrix_d A20, ElDistMatrix_d A21, ElDistMatrix_d A22 );
EL_EXPORT ElError ElSlidePartitionDownDiagonalDist_c
( ElDistMatrix_c ATL, ElDistMatrix_c ATR, ElDistMatrix_c A00, ElDistMatrix_c A01, ElDistMatrix_c A02,
                                          ElDistMatrix_c A10, ElDistMatrix_c A11, ElDistMatrix_c A12,
  ElDistMatrix_c ABL, ElDistMatrix_c ABR, ElDistMatrix_c A20, ElDistMatrix_c A21, ElDistMatrix_c A22 );
EL_EXPORT ElError ElSlidePartitionDownDiagonalDist_z
( ElDistMatrix_z ATL, ElDistMatrix_z ATR, ElDistMatrix_z A00, ElDistMatrix_z A01, ElDistMatrix_z A02,
                                          ElDistMatrix_z A10, ElDistMatrix_z A11, ElDistMatrix_z A12,
  ElDistMatrix_z ABL, ElDistMatrix_z ABR, ElDistMatrix_z A20, ElDistMatrix_z A21, ElDistMatrix_z A22 );

EL_EXPORT ElError ElSlideLockedPartitionDownDiagonal_i
( ElMatrix_i ATL, ElMatrix_i ATR, ElConstMatrix_i A00, ElConstMatrix_i A01, ElConstMatrix_i A02,
                                  ElConstMatrix_i A10, ElConstMatrix_i A11, ElConstMatrix_i A12,
  ElMatrix_i ABL, ElMatrix_i ABR, ElConstMatrix_i A20, ElConstMatrix_i A21, ElConstMatrix_i A22 );
EL_EXPORT ElError ElSlideLockedPartitionDownDiagonal_s
( ElMatrix_s ATL, ElMatrix_s ATR, ElConstMatrix_s A00, ElConstMatrix_s A01, ElConstMatrix_s A02,
                                  ElConstMatrix_s A10, ElConstMatrix_s A11, ElConstMatrix_s A12,
  ElMatrix_s ABL, ElMatrix_s ABR, ElConstMatrix_s A20, ElConstMatrix_s A21, ElConstMatrix_s A22 );
EL_EXPORT ElError ElSlideLockedPartitionDownDiagonal_d
( ElMatrix_d ATL, ElMatrix_d ATR, ElConstMatrix_d A00, ElConstMatrix_d A01, ElConstMatrix_d A02,
                                  ElConstMatrix_d A10, ElConstMatrix_d A11, ElConstMatrix_d A12,
  ElMatrix_d ABL, ElMatrix_d ABR, ElConstMatrix_d A20, ElConstMatrix_d A21, ElConstMatrix_d A22 );
EL_EXPORT ElError ElSlideLockedPartitionDownDiagonal_c
( ElMatrix_c ATL, ElMatrix_c ATR, ElConstMatrix_c A00, ElConstMatrix_c A01, ElConstMatrix_c A02,
                                  ElConstMatrix_c A10, ElConstMatrix_c A11, ElConstMatrix_c A12,
  ElMatrix_c ABL, ElMatrix_c ABR, ElConstMatrix_c A20, ElConstMatrix_c A21, ElConstMatrix_c A22 );
EL_EXPORT ElError ElSlideLockedPartitionDownDiagonal_z
( ElMatrix_z ATL, ElMatrix_z ATR, ElConstMatrix_z A00, ElConstMatrix_z A01, ElConstMatrix_z A02,
                                  ElConstMatrix_z A10, ElConstMatrix_z A11, ElConstMatrix_z A12,
  ElMatrix_z ABL, ElMatrix_z ABR, ElConstMatrix_z A20, ElConstMatrix_z A21, ElConstMatrix_z A22 );
EL_EXPORT ElError ElSlideLockedPartitionDownDiagonalDist_i
( ElDistMatrix_i ATL, ElDistMatrix_i ATR, ElConstDistMatrix_i A00, ElConstDistMatrix_i A01, ElConstDistMatrix_i A02,
                                          ElConstDistMatrix_i A10, ElConstDistMatrix_i A11, ElConstDistMatrix_i A12,
  ElDistMatrix_i ABL, ElDistMatrix_i ABR, ElConstDistMatrix_i A20, ElConstDistMatrix_i A21, ElConstDistMatrix_i A22 );
EL_EXPORT ElError ElSlideLockedPartitionDownDiagonalDist_s
( ElDistMatrix_s ATL, ElDistMatrix_s ATR, ElConstDistMatrix_s A00, ElConstDistMatrix_s A01, ElConstDistMatrix_s A02,
                                          ElConstDistMatrix_s A10, ElConstDistMatrix_s A11, ElConstDistMatrix_s A12,
  ElDistMatrix_s ABL, ElDistMatrix_s ABR, ElConstDistMatrix_s A20, ElConstDistMatrix_s A21, ElConstDistMatrix_s A22 );
EL_EXPORT ElError ElSlideLockedPartitionDownDiagonalDist_d
( ElDistMatrix_d ATL, ElDistMatrix_d ATR, ElConstDistMatrix_d A00, ElConstDistMatrix_d A01, ElConstDistMatrix_d A02,
                                          ElConstDistMatrix_d A10, ElConstDistMatrix_d A11, ElConstDistMatrix_d A12,
  ElDistMatrix_d ABL, ElDistMatrix_d ABR, ElConstDistMatrix_d A20, ElConstDistMatrix_d A21, ElConstDistMatrix_d A22 );
EL_EXPORT ElError ElSlideLockedPartitionDownDiagonalDist_c
( ElDistMatrix_c ATL, ElDistMatrix_c ATR, ElConstDistMatrix_c A00, ElConstDistMatrix_c A01, ElConstDistMatrix_c A02,
                                          ElConstDistMatrix_c A10, ElConstDistMatrix_c A11, ElConstDistMatrix_c A12,
  ElDistMatrix_c ABL, ElDistMatrix_c ABR, ElConstDistMatrix_c A20, ElConstDistMatrix_c A21, ElConstDistMatrix_c A22 );
EL_EXPORT ElError ElSlideLockedPartitionDownDiagonalDist_z
( ElDistMatrix_z ATL, ElDistMatrix_z ATR, ElConstDistMatrix_z A00, ElConstDistMatrix_z A01, ElConstDistMatrix_z A02,
                                          ElConstDistMatrix_z A10, ElConstDistMatrix_z A11, ElConstDistMatrix_z A12,
  ElDistMatrix_z ABL, ElDistMatrix_z ABR, ElConstDistMatrix_z A20, ElConstDistMatrix_z A21, ElConstDistMatrix_z A22 );

/* Up a diagonal
   ============= */
EL_EXPORT ElError ElSlidePartitionUpDiagonal_i
( ElMatrix_i ATL, ElMatrix_i ATR, ElMatrix_i A00, ElMatrix_i A01, ElMatrix_i A02,
                                  ElMatrix_i A10, ElMatrix_i A11, ElMatrix_i A12,
  ElMatrix_i ABL, ElMatrix_i ABR, ElMatrix_i A20, ElMatrix_i A21, ElMatrix_i A22 );
EL_EXPORT ElError ElSlidePartitionUpDiagonal_s
( ElMatrix_s ATL, ElMatrix_s ATR, ElMatrix_s A00, ElMatrix_s A01, ElMatrix_s A02,
                                  ElMatrix_s A10, ElMatrix_s A11, ElMatrix_s A12,
  ElMatrix_s ABL, ElMatrix_s ABR, ElMatrix_s A20, ElMatrix_s A21, ElMatrix_s A22 );
EL_EXPORT ElError ElSlidePartitionUpDiagonal_d
( ElMatrix_d ATL, ElMatrix_d ATR, ElMatrix_d A00, ElMatrix_d A01, ElMatrix_d A02,
                                  ElMatrix_d A10, ElMatrix_d A11, ElMatrix_d A12,
  ElMatrix_d ABL, ElMatrix_d ABR, ElMatrix_d A20, ElMatrix_d A21, ElMatrix_d A22 );
EL_EXPORT ElError ElSlidePartitionUpDiagonal_c
( ElMatrix_c ATL, ElMatrix_c ATR, ElMatrix_c A00, ElMatrix_c A01, ElMatrix_c A02,
                                  ElMatrix_c A10, ElMatrix_c A11, ElMatrix_c A12,
  ElMatrix_c ABL, ElMatrix_c ABR, ElMatrix_c A20, ElMatrix_c A21, ElMatrix_c A22 );
EL_EXPORT ElError ElSlidePartitionUpDiagonal_z
( ElMatrix_z ATL, ElMatrix_z ATR, ElMatrix_z A00, ElMatrix_z A01, ElMatrix_z A02,
                                  ElMatrix_z A10, ElMatrix_z A11, ElMatrix_z A12,
  ElMatrix_z ABL, ElMatrix_z ABR, ElMatrix_z A20, ElMatrix_z A21, ElMatrix_z A22 );
EL_EXPORT ElError ElSlidePartitionUpDiagonalDist_i
( ElDistMatrix_i ATL, ElDistMatrix_i ATR, ElDistMatrix_i A00, ElDistMatrix_i A01, ElDistMatrix_i A02,
                                          ElDistMatrix_i A10, ElDistMatrix_i A11, ElDistMatrix_i A12,
  ElDistMatrix_i ABL, ElDistMatrix_i ABR, ElDistMatrix_i A20, ElDistMatrix_i A21, ElDistMatrix_i A22 );
EL_EXPORT ElError ElSlidePartitionUpDiagonalDist_s
( ElDistMatrix_s ATL, ElDistMatrix_s ATR, ElDistMatrix_s A00, ElDistMatrix_s A01, ElDistMatrix_s A02,
                                          ElDistMatrix_s A10, ElDistMatrix_s A11, ElDistMatrix_s A12,
  ElDistMatrix_s ABL, ElDistMatrix_s ABR, ElDistMatrix_s A20, ElDistMatrix_s A21, ElDistMatrix_s A22 );
EL_EXPORT ElError ElSlidePartitionUpDiagonalDist_d
( ElDistMatrix_d ATL, ElDistMatrix_d ATR, ElDistMatrix_d A00, ElDistMatrix_d A01, ElDistMatrix_d A02,
                                          ElDistMatrix_d A10, ElDistMatrix_d A11, ElDistMatrix_d A12,
  ElDistMatrix_d ABL, ElDistMatrix_d ABR, ElDistMatrix_d A20, ElDistMatrix_d A21, ElDistMatrix_d A22 );
EL_EXPORT ElError ElSlidePartitionUpDiagonalDist_c
( ElDistMatrix_c ATL, ElDistMatrix_c ATR, ElDistMatrix_c A00, ElDistMatrix_c A01, ElDistMatrix_c A02,
                                          ElDistMatrix_c A10, ElDistMatrix_c A11, ElDistMatrix_c A12,
  ElDistMatrix_c ABL, ElDistMatrix_c ABR, ElDistMatrix_c A20, ElDistMatrix_c A21, ElDistMatrix_c A22 );
EL_EXPORT ElError ElSlidePartitionUpDiagonalDist_z
( ElDistMatrix_z ATL, ElDistMatrix_z ATR, ElDistMatrix_z A00, ElDistMatrix_z A01, ElDistMatrix_z A02,
                                          ElDistMatrix_z A10, ElDistMatrix_z A11, ElDistMatrix_z A12,
  ElDistMatrix_z ABL, ElDistMatrix_z ABR, ElDistMatrix_z A20, ElDistMatrix_z A21, ElDistMatrix_z A22 );

EL_EXPORT ElError ElSlideLockedPartitionUpDiagonal_i
( ElMatrix_i ATL, ElMatrix_i ATR, ElConstMatrix_i A00, ElConstMatrix_i A01, ElConstMatrix_i A02,
                                  ElConstMatrix_i A10, ElConstMatrix_i A11, ElConstMatrix_i A12,
  ElMatrix_i ABL, ElMatrix_i ABR, ElConstMatrix_i A20, ElConstMatrix_i A21, ElConstMatrix_i A22 );
EL_EXPORT ElError ElSlideLockedPartitionUpDiagonal_s
( ElMatrix_s ATL, ElMatrix_s ATR, ElConstMatrix_s A00, ElConstMatrix_s A01, ElConstMatrix_s A02,
                                  ElConstMatrix_s A10, ElConstMatrix_s A11, ElConstMatrix_s A12,
  ElMatrix_s ABL, ElMatrix_s ABR, ElConstMatrix_s A20, ElConstMatrix_s A21, ElConstMatrix_s A22 );
EL_EXPORT ElError ElSlideLockedPartitionUpDiagonal_d
( ElMatrix_d ATL, ElMatrix_d ATR, ElConstMatrix_d A00, ElConstMatrix_d A01, ElConstMatrix_d A02,
                                  ElConstMatrix_d A10, ElConstMatrix_d A11, ElConstMatrix_d A12,
  ElMatrix_d ABL, ElMatrix_d ABR, ElConstMatrix_d A20, ElConstMatrix_d A21, ElConstMatrix_d A22 );
EL_EXPORT ElError ElSlideLockedPartitionUpDiagonal_c
( ElMatrix_c ATL, ElMatrix_c ATR, ElConstMatrix_c A00, ElConstMatrix_c A01, ElConstMatrix_c A02,
                                  ElConstMatrix_c A10, ElConstMatrix_c A11, ElConstMatrix_c A12,
  ElMatrix_c ABL, ElMatrix_c ABR, ElConstMatrix_c A20, ElConstMatrix_c A21, ElConstMatrix_c A22 );
EL_EXPORT ElError ElSlideLockedPartitionUpDiagonal_z
( ElMatrix_z ATL, ElMatrix_z ATR, ElConstMatrix_z A00, ElConstMatrix_z A01, ElConstMatrix_z A02,
                                  ElConstMatrix_z A10, ElConstMatrix_z A11, ElConstMatrix_z A12,
  ElMatrix_z ABL, ElMatrix_z ABR, ElConstMatrix_z A20, ElConstMatrix_z A21, ElConstMatrix_z A22 );
EL_EXPORT ElError ElSlideLockedPartitionUpDiagonalDist_i
( ElDistMatrix_i ATL, ElDistMatrix_i ATR, ElConstDistMatrix_i A00, ElConstDistMatrix_i A01, ElConstDistMatrix_i A02,
                                          ElConstDistMatrix_i A10, ElConstDistMatrix_i A11, ElConstDistMatrix_i A12,
  ElDistMatrix_i ABL, ElDistMatrix_i ABR, ElConstDistMatrix_i A20, ElConstDistMatrix_i A21, ElConstDistMatrix_i A22 );
EL_EXPORT ElError ElSlideLockedPartitionUpDiagonalDist_s
( ElDistMatrix_s ATL, ElDistMatrix_s ATR, ElConstDistMatrix_s A00, ElConstDistMatrix_s A01, ElConstDistMatrix_s A02,
                                          ElConstDistMatrix_s A10, ElConstDistMatrix_s A11, ElConstDistMatrix_s A12,
  ElDistMatrix_s ABL, ElDistMatrix_s ABR, ElConstDistMatrix_s A20, ElConstDistMatrix_s A21, ElConstDistMatrix_s A22 );
EL_EXPORT ElError ElSlideLockedPartitionUpDiagonalDist_d
( ElDistMatrix_d ATL, ElDistMatrix_d ATR, ElConstDistMatrix_d A00, ElConstDistMatrix_d A01, ElConstDistMatrix_d A02,
                                          ElConstDistMatrix_d A10, ElConstDistMatrix_d A11, ElConstDistMatrix_d A12,
  ElDistMatrix_d ABL, ElDistMatrix_d ABR, ElConstDistMatrix_d A20, ElConstDistMatrix_d A21, ElConstDistMatrix_d A22 );
EL_EXPORT ElError ElSlideLockedPartitionUpDiagonalDist_c
( ElDistMatrix_c ATL, ElDistMatrix_c ATR, ElConstDistMatrix_c A00, ElConstDistMatrix_c A01, ElConstDistMatrix_c A02,
                                          ElConstDistMatrix_c A10, ElConstDistMatrix_c A11, ElConstDistMatrix_c A12,
  ElDistMatrix_c ABL, ElDistMatrix_c ABR, ElConstDistMatrix_c A20, ElConstDistMatrix_c A21, ElConstDistMatrix_c A22 );
EL_EXPORT ElError ElSlideLockedPartitionUpDiagonalDist_z
( ElDistMatrix_z ATL, ElDistMatrix_z ATR, ElConstDistMatrix_z A00, ElConstDistMatrix_z A01, ElConstDistMatrix_z A02,
                                          ElConstDistMatrix_z A10, ElConstDistMatrix_z A11, ElConstDistMatrix_z A12,
  ElDistMatrix_z ABL, ElDistMatrix_z ABR, ElConstDistMatrix_z A20, ElConstDistMatrix_z A21, ElConstDistMatrix_z A22 );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_FLAMEPART_SLIDEPARTITION_C_H */
