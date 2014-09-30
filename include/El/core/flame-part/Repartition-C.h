/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_FLAMEPART_REPARTITION_C_H
#define EL_FLAMEPART_REPARTITION_C_H

#ifdef __cplusplus
extern "C" {
#endif

/* Downward
   ======== */
ElError ElRepartitionDown_i
( ElMatrix_i AT, ElMatrix_i A0,
                 ElMatrix_i A1,
  ElMatrix_i AB, ElMatrix_i A2, ElInt bsize );
ElError ElRepartitionDown_s
( ElMatrix_s AT, ElMatrix_s A0,
                 ElMatrix_s A1,
  ElMatrix_s AB, ElMatrix_s A2, ElInt bsize );
ElError ElRepartitionDown_d
( ElMatrix_d AT, ElMatrix_d A0,
                 ElMatrix_d A1,
  ElMatrix_d AB, ElMatrix_d A2, ElInt bsize );
ElError ElRepartitionDown_c
( ElMatrix_c AT, ElMatrix_c A0,
                 ElMatrix_c A1,
  ElMatrix_c AB, ElMatrix_c A2, ElInt bsize );
ElError ElRepartitionDown_z
( ElMatrix_z AT, ElMatrix_z A0,
                 ElMatrix_z A1,
  ElMatrix_z AB, ElMatrix_z A2, ElInt bsize );
ElError ElRepartitionDownDist_i
( ElDistMatrix_i AT, ElDistMatrix_i A0,
                     ElDistMatrix_i A1,
  ElDistMatrix_i AB, ElDistMatrix_i A2, ElInt bsize );
ElError ElRepartitionDownDist_s
( ElDistMatrix_s AT, ElDistMatrix_s A0,
                     ElDistMatrix_s A1,
  ElDistMatrix_s AB, ElDistMatrix_s A2, ElInt bsize );
ElError ElRepartitionDownDist_d
( ElDistMatrix_d AT, ElDistMatrix_d A0,
                     ElDistMatrix_d A1,
  ElDistMatrix_d AB, ElDistMatrix_d A2, ElInt bsize );
ElError ElRepartitionDownDist_c
( ElDistMatrix_c AT, ElDistMatrix_c A0,
                     ElDistMatrix_c A1,
  ElDistMatrix_c AB, ElDistMatrix_c A2, ElInt bsize );
ElError ElRepartitionDownDist_z
( ElDistMatrix_z AT, ElDistMatrix_z A0,
                     ElDistMatrix_z A1,
  ElDistMatrix_z AB, ElDistMatrix_z A2, ElInt bsize );

ElError ElLockedRepartitionDown_i
( ElConstMatrix_i AT, ElMatrix_i A0,
                      ElMatrix_i A1,
  ElConstMatrix_i AB, ElMatrix_i A2, ElInt bsize );
ElError ElLockedRepartitionDown_s
( ElConstMatrix_s AT, ElMatrix_s A0,
                      ElMatrix_s A1,
  ElConstMatrix_s AB, ElMatrix_s A2, ElInt bsize );
ElError ElLockedRepartitionDown_d
( ElConstMatrix_d AT, ElMatrix_d A0,
                      ElMatrix_d A1,
  ElConstMatrix_d AB, ElMatrix_d A2, ElInt bsize );
ElError ElLockedRepartitionDown_c
( ElConstMatrix_c AT, ElMatrix_c A0,
                      ElMatrix_c A1,
  ElConstMatrix_c AB, ElMatrix_c A2, ElInt bsize );
ElError ElLockedRepartitionDown_z
( ElConstMatrix_z AT, ElMatrix_z A0,
                      ElMatrix_z A1,
  ElConstMatrix_z AB, ElMatrix_z A2, ElInt bsize );
ElError ElLockedRepartitionDownDist_i
( ElConstDistMatrix_i AT, ElDistMatrix_i A0,
                          ElDistMatrix_i A1,
  ElConstDistMatrix_i AB, ElDistMatrix_i A2, ElInt bsize );
ElError ElLockedRepartitionDownDist_s
( ElConstDistMatrix_s AT, ElDistMatrix_s A0,
                          ElDistMatrix_s A1,
  ElConstDistMatrix_s AB, ElDistMatrix_s A2, ElInt bsize );
ElError ElLockedRepartitionDownDist_d
( ElConstDistMatrix_d AT, ElDistMatrix_d A0,
                          ElDistMatrix_d A1,
  ElConstDistMatrix_d AB, ElDistMatrix_d A2, ElInt bsize );
ElError ElLockedRepartitionDownDist_c
( ElConstDistMatrix_c AT, ElDistMatrix_c A0,
                          ElDistMatrix_c A1,
  ElConstDistMatrix_c AB, ElDistMatrix_c A2, ElInt bsize );
ElError ElLockedRepartitionDownDist_z
( ElConstDistMatrix_z AT, ElDistMatrix_z A0,
                          ElDistMatrix_z A1,
  ElConstDistMatrix_z AB, ElDistMatrix_z A2, ElInt bsize );

/* Upward
   ====== */
ElError ElRepartitionUp_i
( ElMatrix_i AT, ElMatrix_i A0,
                 ElMatrix_i A1,
  ElMatrix_i AB, ElMatrix_i A2, ElInt bsize );
ElError ElRepartitionUp_s
( ElMatrix_s AT, ElMatrix_s A0,
                 ElMatrix_s A1,
  ElMatrix_s AB, ElMatrix_s A2, ElInt bsize );
ElError ElRepartitionUp_d
( ElMatrix_d AT, ElMatrix_d A0,
                 ElMatrix_d A1,
  ElMatrix_d AB, ElMatrix_d A2, ElInt bsize );
ElError ElRepartitionUp_c
( ElMatrix_c AT, ElMatrix_c A0,
                 ElMatrix_c A1,
  ElMatrix_c AB, ElMatrix_c A2, ElInt bsize );
ElError ElRepartitionUp_z
( ElMatrix_z AT, ElMatrix_z A0,
                 ElMatrix_z A1,
  ElMatrix_z AB, ElMatrix_z A2, ElInt bsize );
ElError ElRepartitionUpDist_i
( ElDistMatrix_i AT, ElDistMatrix_i A0,
                     ElDistMatrix_i A1,
  ElDistMatrix_i AB, ElDistMatrix_i A2, ElInt bsize );
ElError ElRepartitionUpDist_s
( ElDistMatrix_s AT, ElDistMatrix_s A0,
                     ElDistMatrix_s A1,
  ElDistMatrix_s AB, ElDistMatrix_s A2, ElInt bsize );
ElError ElRepartitionUpDist_d
( ElDistMatrix_d AT, ElDistMatrix_d A0,
                     ElDistMatrix_d A1,
  ElDistMatrix_d AB, ElDistMatrix_d A2, ElInt bsize );
ElError ElRepartitionUpDist_c
( ElDistMatrix_c AT, ElDistMatrix_c A0,
                     ElDistMatrix_c A1,
  ElDistMatrix_c AB, ElDistMatrix_c A2, ElInt bsize );
ElError ElRepartitionUpDist_z
( ElDistMatrix_z AT, ElDistMatrix_z A0,
                     ElDistMatrix_z A1,
  ElDistMatrix_z AB, ElDistMatrix_z A2, ElInt bsize );

ElError ElLockedRepartitionUp_i
( ElConstMatrix_i AT, ElMatrix_i A0,
                      ElMatrix_i A1,
  ElConstMatrix_i AB, ElMatrix_i A2, ElInt bsize );
ElError ElLockedRepartitionUp_s
( ElConstMatrix_s AT, ElMatrix_s A0,
                      ElMatrix_s A1,
  ElConstMatrix_s AB, ElMatrix_s A2, ElInt bsize );
ElError ElLockedRepartitionUp_d
( ElConstMatrix_d AT, ElMatrix_d A0,
                      ElMatrix_d A1,
  ElConstMatrix_d AB, ElMatrix_d A2, ElInt bsize );
ElError ElLockedRepartitionUp_c
( ElConstMatrix_c AT, ElMatrix_c A0,
                      ElMatrix_c A1,
  ElConstMatrix_c AB, ElMatrix_c A2, ElInt bsize );
ElError ElLockedRepartitionUp_z
( ElConstMatrix_z AT, ElMatrix_z A0,
                      ElMatrix_z A1,
  ElConstMatrix_z AB, ElMatrix_z A2, ElInt bsize );
ElError ElLockedRepartitionUpDist_i
( ElConstDistMatrix_i AT, ElDistMatrix_i A0,
                          ElDistMatrix_i A1,
  ElConstDistMatrix_i AB, ElDistMatrix_i A2, ElInt bsize );
ElError ElLockedRepartitionUpDist_s
( ElConstDistMatrix_s AT, ElDistMatrix_s A0,
                          ElDistMatrix_s A1,
  ElConstDistMatrix_s AB, ElDistMatrix_s A2, ElInt bsize );
ElError ElLockedRepartitionUpDist_d
( ElConstDistMatrix_d AT, ElDistMatrix_d A0,
                          ElDistMatrix_d A1,
  ElConstDistMatrix_d AB, ElDistMatrix_d A2, ElInt bsize );
ElError ElLockedRepartitionUpDist_c
( ElConstDistMatrix_c AT, ElDistMatrix_c A0,
                          ElDistMatrix_c A1,
  ElConstDistMatrix_c AB, ElDistMatrix_c A2, ElInt bsize );
ElError ElLockedRepartitionUpDist_z
( ElConstDistMatrix_z AT, ElDistMatrix_z A0,
                          ElDistMatrix_z A1,
  ElConstDistMatrix_z AB, ElDistMatrix_z A2, ElInt bsize );

/* Rightward
   ========= */
ElError ElRepartitionRight_i
( ElMatrix_i AL, ElMatrix_i AR, 
  ElMatrix_i A0, ElMatrix_i A1, ElMatrix_i A2, ElInt bsize );
ElError ElRepartitionRight_s
( ElMatrix_s AL, ElMatrix_s AR, 
  ElMatrix_s A0, ElMatrix_s A1, ElMatrix_s A2, ElInt bsize );
ElError ElRepartitionRight_d
( ElMatrix_d AL, ElMatrix_d AR, 
  ElMatrix_d A0, ElMatrix_d A1, ElMatrix_d A2, ElInt bsize );
ElError ElRepartitionRight_c
( ElMatrix_c AL, ElMatrix_c AR, 
  ElMatrix_c A0, ElMatrix_c A1, ElMatrix_c A2, ElInt bsize );
ElError ElRepartitionRight_z
( ElMatrix_z AL, ElMatrix_z AR, 
  ElMatrix_z A0, ElMatrix_z A1, ElMatrix_z A2, ElInt bsize );
ElError ElRepartitionRightDist_i
( ElDistMatrix_i AL, ElDistMatrix_i AR, 
  ElDistMatrix_i A0, ElDistMatrix_i A1, ElDistMatrix_i A2, ElInt bsize );
ElError ElRepartitionRightDist_s
( ElDistMatrix_s AL, ElDistMatrix_s AR, 
  ElDistMatrix_s A0, ElDistMatrix_s A1, ElDistMatrix_s A2, ElInt bsize );
ElError ElRepartitionRightDist_d
( ElDistMatrix_d AL, ElDistMatrix_d AR, 
  ElDistMatrix_d A0, ElDistMatrix_d A1, ElDistMatrix_d A2, ElInt bsize );
ElError ElRepartitionRightDist_c
( ElDistMatrix_c AL, ElDistMatrix_c AR, 
  ElDistMatrix_c A0, ElDistMatrix_c A1, ElDistMatrix_c A2, ElInt bsize );
ElError ElRepartitionRightDist_z
( ElDistMatrix_z AL, ElDistMatrix_z AR, 
  ElDistMatrix_z A0, ElDistMatrix_z A1, ElDistMatrix_z A2, ElInt bsize );

ElError ElLockedRepartitionRight_i
( ElConstMatrix_i AL, ElConstMatrix_i AR, 
  ElMatrix_i A0, ElMatrix_i A1, ElMatrix_i A2, ElInt bsize );
ElError ElLockedRepartitionRight_s
( ElConstMatrix_s AL, ElConstMatrix_s AR, 
  ElMatrix_s A0, ElMatrix_s A1, ElMatrix_s A2, ElInt bsize );
ElError ElLockedRepartitionRight_d
( ElConstMatrix_d AL, ElConstMatrix_d AR, 
  ElMatrix_d A0, ElMatrix_d A1, ElMatrix_d A2, ElInt bsize );
ElError ElLockedRepartitionRight_c
( ElConstMatrix_c AL, ElConstMatrix_c AR, 
  ElMatrix_c A0, ElMatrix_c A1, ElMatrix_c A2, ElInt bsize );
ElError ElLockedRepartitionRight_z
( ElConstMatrix_z AL, ElConstMatrix_z AR, 
  ElMatrix_z A0, ElMatrix_z A1, ElMatrix_z A2, ElInt bsize );
ElError ElLockedRepartitionRightDist_i
( ElConstDistMatrix_i AL, ElConstDistMatrix_i AR, 
  ElDistMatrix_i A0, ElDistMatrix_i A1, ElDistMatrix_i A2, ElInt bsize );
ElError ElLockedRepartitionRightDist_s
( ElConstDistMatrix_s AL, ElConstDistMatrix_s AR, 
  ElDistMatrix_s A0, ElDistMatrix_s A1, ElDistMatrix_s A2, ElInt bsize );
ElError ElLockedRepartitionRightDist_d
( ElConstDistMatrix_d AL, ElConstDistMatrix_d AR, 
  ElDistMatrix_d A0, ElDistMatrix_d A1, ElDistMatrix_d A2, ElInt bsize );
ElError ElLockedRepartitionRightDist_c
( ElConstDistMatrix_c AL, ElConstDistMatrix_c AR, 
  ElDistMatrix_c A0, ElDistMatrix_c A1, ElDistMatrix_c A2, ElInt bsize );
ElError ElLockedRepartitionRightDist_z
( ElConstDistMatrix_z AL, ElConstDistMatrix_z AR, 
  ElDistMatrix_z A0, ElDistMatrix_z A1, ElDistMatrix_z A2, ElInt bsize );

/* Leftward
   ======== */
ElError ElRepartitionLeft_i
( ElMatrix_i AL, ElMatrix_i AR, 
  ElMatrix_i A0, ElMatrix_i A1, ElMatrix_i A2, ElInt bsize );
ElError ElRepartitionLeft_s
( ElMatrix_s AL, ElMatrix_s AR, 
  ElMatrix_s A0, ElMatrix_s A1, ElMatrix_s A2, ElInt bsize );
ElError ElRepartitionLeft_d
( ElMatrix_d AL, ElMatrix_d AR, 
  ElMatrix_d A0, ElMatrix_d A1, ElMatrix_d A2, ElInt bsize );
ElError ElRepartitionLeft_c
( ElMatrix_c AL, ElMatrix_c AR, 
  ElMatrix_c A0, ElMatrix_c A1, ElMatrix_c A2, ElInt bsize );
ElError ElRepartitionLeft_z
( ElMatrix_z AL, ElMatrix_z AR, 
  ElMatrix_z A0, ElMatrix_z A1, ElMatrix_z A2, ElInt bsize );
ElError ElRepartitionLeftDist_i
( ElDistMatrix_i AL, ElDistMatrix_i AR, 
  ElDistMatrix_i A0, ElDistMatrix_i A1, ElDistMatrix_i A2, ElInt bsize );
ElError ElRepartitionLeftDist_s
( ElDistMatrix_s AL, ElDistMatrix_s AR, 
  ElDistMatrix_s A0, ElDistMatrix_s A1, ElDistMatrix_s A2, ElInt bsize );
ElError ElRepartitionLeftDist_d
( ElDistMatrix_d AL, ElDistMatrix_d AR, 
  ElDistMatrix_d A0, ElDistMatrix_d A1, ElDistMatrix_d A2, ElInt bsize );
ElError ElRepartitionLeftDist_c
( ElDistMatrix_c AL, ElDistMatrix_c AR, 
  ElDistMatrix_c A0, ElDistMatrix_c A1, ElDistMatrix_c A2, ElInt bsize );
ElError ElRepartitionLeftDist_z
( ElDistMatrix_z AL, ElDistMatrix_z AR, 
  ElDistMatrix_z A0, ElDistMatrix_z A1, ElDistMatrix_z A2, ElInt bsize );

ElError ElLockedRepartitionLeft_i
( ElConstMatrix_i AL, ElConstMatrix_i AR, 
  ElMatrix_i A0, ElMatrix_i A1, ElMatrix_i A2, ElInt bsize );
ElError ElLockedRepartitionLeft_s
( ElConstMatrix_s AL, ElConstMatrix_s AR, 
  ElMatrix_s A0, ElMatrix_s A1, ElMatrix_s A2, ElInt bsize );
ElError ElLockedRepartitionLeft_d
( ElConstMatrix_d AL, ElConstMatrix_d AR, 
  ElMatrix_d A0, ElMatrix_d A1, ElMatrix_d A2, ElInt bsize );
ElError ElLockedRepartitionLeft_c
( ElConstMatrix_c AL, ElConstMatrix_c AR, 
  ElMatrix_c A0, ElMatrix_c A1, ElMatrix_c A2, ElInt bsize );
ElError ElLockedRepartitionLeft_z
( ElConstMatrix_z AL, ElConstMatrix_z AR, 
  ElMatrix_z A0, ElMatrix_z A1, ElMatrix_z A2, ElInt bsize );
ElError ElLockedRepartitionLeftDist_i
( ElConstDistMatrix_i AL, ElConstDistMatrix_i AR, 
  ElDistMatrix_i A0, ElDistMatrix_i A1, ElDistMatrix_i A2, ElInt bsize );
ElError ElLockedRepartitionLeftDist_s
( ElConstDistMatrix_s AL, ElConstDistMatrix_s AR, 
  ElDistMatrix_s A0, ElDistMatrix_s A1, ElDistMatrix_s A2, ElInt bsize );
ElError ElLockedRepartitionLeftDist_d
( ElConstDistMatrix_d AL, ElConstDistMatrix_d AR, 
  ElDistMatrix_d A0, ElDistMatrix_d A1, ElDistMatrix_d A2, ElInt bsize );
ElError ElLockedRepartitionLeftDist_c
( ElConstDistMatrix_c AL, ElConstDistMatrix_c AR, 
  ElDistMatrix_c A0, ElDistMatrix_c A1, ElDistMatrix_c A2, ElInt bsize );
ElError ElLockedRepartitionLeftDist_z
( ElConstDistMatrix_z AL, ElConstDistMatrix_z AR, 
  ElDistMatrix_z A0, ElDistMatrix_z A1, ElDistMatrix_z A2, ElInt bsize );

/* Down a diagonal
   =============== */
ElError ElRepartitionDownDiagonal_i
( ElMatrix_i ATL, ElMatrix_i ATR, ElMatrix_i A00, ElMatrix_i A01, ElMatrix_i A02,
                                  ElMatrix_i A10, ElMatrix_i A11, ElMatrix_i A12,
  ElMatrix_i ABL, ElMatrix_i ABR, ElMatrix_i A20, ElMatrix_i A21, ElMatrix_i A22, 
  ElInt bsize );
ElError ElRepartitionDownDiagonal_s
( ElMatrix_s ATL, ElMatrix_s ATR, ElMatrix_s A00, ElMatrix_s A01, ElMatrix_s A02,
                                  ElMatrix_s A10, ElMatrix_s A11, ElMatrix_s A12,
  ElMatrix_s ABL, ElMatrix_s ABR, ElMatrix_s A20, ElMatrix_s A21, ElMatrix_s A22, 
  ElInt bsize );
ElError ElRepartitionDownDiagonal_d
( ElMatrix_d ATL, ElMatrix_d ATR, ElMatrix_d A00, ElMatrix_d A01, ElMatrix_d A02,
                                  ElMatrix_d A10, ElMatrix_d A11, ElMatrix_d A12,
  ElMatrix_d ABL, ElMatrix_d ABR, ElMatrix_d A20, ElMatrix_d A21, ElMatrix_d A22, 
  ElInt bsize );
ElError ElRepartitionDownDiagonal_c
( ElMatrix_c ATL, ElMatrix_c ATR, ElMatrix_c A00, ElMatrix_c A01, ElMatrix_c A02,
                                  ElMatrix_c A10, ElMatrix_c A11, ElMatrix_c A12,
  ElMatrix_c ABL, ElMatrix_c ABR, ElMatrix_c A20, ElMatrix_c A21, ElMatrix_c A22, 
  ElInt bsize );
ElError ElRepartitionDownDiagonal_z
( ElMatrix_z ATL, ElMatrix_z ATR, ElMatrix_z A00, ElMatrix_z A01, ElMatrix_z A02,
                                  ElMatrix_z A10, ElMatrix_z A11, ElMatrix_z A12,
  ElMatrix_z ABL, ElMatrix_z ABR, ElMatrix_z A20, ElMatrix_z A21, ElMatrix_z A22, 
  ElInt bsize );
ElError ElRepartitionDownDiagonalDist_i
( ElDistMatrix_i ATL, ElDistMatrix_i ATR, ElDistMatrix_i A00, ElDistMatrix_i A01, ElDistMatrix_i A02,
                                          ElDistMatrix_i A10, ElDistMatrix_i A11, ElDistMatrix_i A12,
  ElDistMatrix_i ABL, ElDistMatrix_i ABR, ElDistMatrix_i A20, ElDistMatrix_i A21, ElDistMatrix_i A22, 
  ElInt bsize );
ElError ElRepartitionDownDiagonalDist_s
( ElDistMatrix_s ATL, ElDistMatrix_s ATR, ElDistMatrix_s A00, ElDistMatrix_s A01, ElDistMatrix_s A02,
                                          ElDistMatrix_s A10, ElDistMatrix_s A11, ElDistMatrix_s A12,
  ElDistMatrix_s ABL, ElDistMatrix_s ABR, ElDistMatrix_s A20, ElDistMatrix_s A21, ElDistMatrix_s A22, 
  ElInt bsize );
ElError ElRepartitionDownDiagonalDist_d
( ElDistMatrix_d ATL, ElDistMatrix_d ATR, ElDistMatrix_d A00, ElDistMatrix_d A01, ElDistMatrix_d A02,
                                          ElDistMatrix_d A10, ElDistMatrix_d A11, ElDistMatrix_d A12,
  ElDistMatrix_d ABL, ElDistMatrix_d ABR, ElDistMatrix_d A20, ElDistMatrix_d A21, ElDistMatrix_d A22, 
  ElInt bsize );
ElError ElRepartitionDownDiagonalDist_c
( ElDistMatrix_c ATL, ElDistMatrix_c ATR, ElDistMatrix_c A00, ElDistMatrix_c A01, ElDistMatrix_c A02,
                                          ElDistMatrix_c A10, ElDistMatrix_c A11, ElDistMatrix_c A12,
  ElDistMatrix_c ABL, ElDistMatrix_c ABR, ElDistMatrix_c A20, ElDistMatrix_c A21, ElDistMatrix_c A22, 
  ElInt bsize );
ElError ElRepartitionDownDiagonalDist_z
( ElDistMatrix_z ATL, ElDistMatrix_z ATR, ElDistMatrix_z A00, ElDistMatrix_z A01, ElDistMatrix_z A02,
                                          ElDistMatrix_z A10, ElDistMatrix_z A11, ElDistMatrix_z A12,
  ElDistMatrix_z ABL, ElDistMatrix_z ABR, ElDistMatrix_z A20, ElDistMatrix_z A21, ElDistMatrix_z A22, 
  ElInt bsize );

ElError ElLockedRepartitionDownDiagonal_i
( ElConstMatrix_i ATL, ElConstMatrix_i ATR, ElMatrix_i A00, ElMatrix_i A01, ElMatrix_i A02,
                                            ElMatrix_i A10, ElMatrix_i A11, ElMatrix_i A12,
  ElConstMatrix_i ABL, ElConstMatrix_i ABR, ElMatrix_i A20, ElMatrix_i A21, ElMatrix_i A22, 
  ElInt bsize );
ElError ElLockedRepartitionDownDiagonal_s
( ElConstMatrix_s ATL, ElConstMatrix_s ATR, ElMatrix_s A00, ElMatrix_s A01, ElMatrix_s A02,
                                            ElMatrix_s A10, ElMatrix_s A11, ElMatrix_s A12,
  ElConstMatrix_s ABL, ElConstMatrix_s ABR, ElMatrix_s A20, ElMatrix_s A21, ElMatrix_s A22, 
  ElInt bsize );
ElError ElLockedRepartitionDownDiagonal_d
( ElConstMatrix_d ATL, ElConstMatrix_d ATR, ElMatrix_d A00, ElMatrix_d A01, ElMatrix_d A02,
                                            ElMatrix_d A10, ElMatrix_d A11, ElMatrix_d A12,
  ElConstMatrix_d ABL, ElConstMatrix_d ABR, ElMatrix_d A20, ElMatrix_d A21, ElMatrix_d A22, 
  ElInt bsize );
ElError ElLockedRepartitionDownDiagonal_c
( ElConstMatrix_c ATL, ElConstMatrix_c ATR, ElMatrix_c A00, ElMatrix_c A01, ElMatrix_c A02,
                                            ElMatrix_c A10, ElMatrix_c A11, ElMatrix_c A12,
  ElConstMatrix_c ABL, ElConstMatrix_c ABR, ElMatrix_c A20, ElMatrix_c A21, ElMatrix_c A22, 
  ElInt bsize );
ElError ElLockedRepartitionDownDiagonal_z
( ElConstMatrix_z ATL, ElConstMatrix_z ATR, ElMatrix_z A00, ElMatrix_z A01, ElMatrix_z A02,
                                            ElMatrix_z A10, ElMatrix_z A11, ElMatrix_z A12,
  ElConstMatrix_z ABL, ElConstMatrix_z ABR, ElMatrix_z A20, ElMatrix_z A21, ElMatrix_z A22, 
  ElInt bsize );
ElError ElLockedRepartitionDownDiagonalDist_i
( ElConstDistMatrix_i ATL, ElConstDistMatrix_i ATR, ElDistMatrix_i A00, ElDistMatrix_i A01, ElDistMatrix_i A02,
                                                    ElDistMatrix_i A10, ElDistMatrix_i A11, ElDistMatrix_i A12,
  ElConstDistMatrix_i ABL, ElConstDistMatrix_i ABR, ElDistMatrix_i A20, ElDistMatrix_i A21, ElDistMatrix_i A22, 
  ElInt bsize );
ElError ElLockedRepartitionDownDiagonalDist_s
( ElConstDistMatrix_s ATL, ElConstDistMatrix_s ATR, ElDistMatrix_s A00, ElDistMatrix_s A01, ElDistMatrix_s A02,
                                                    ElDistMatrix_s A10, ElDistMatrix_s A11, ElDistMatrix_s A12,
  ElConstDistMatrix_s ABL, ElConstDistMatrix_s ABR, ElDistMatrix_s A20, ElDistMatrix_s A21, ElDistMatrix_s A22, 
  ElInt bsize );
ElError ElLockedRepartitionDownDiagonalDist_d
( ElConstDistMatrix_d ATL, ElConstDistMatrix_d ATR, ElDistMatrix_d A00, ElDistMatrix_d A01, ElDistMatrix_d A02,
                                                    ElDistMatrix_d A10, ElDistMatrix_d A11, ElDistMatrix_d A12,
  ElConstDistMatrix_d ABL, ElConstDistMatrix_d ABR, ElDistMatrix_d A20, ElDistMatrix_d A21, ElDistMatrix_d A22, 
  ElInt bsize );
ElError ElLockedRepartitionDownDiagonalDist_c
( ElConstDistMatrix_c ATL, ElConstDistMatrix_c ATR, ElDistMatrix_c A00, ElDistMatrix_c A01, ElDistMatrix_c A02,
                                                    ElDistMatrix_c A10, ElDistMatrix_c A11, ElDistMatrix_c A12,
  ElConstDistMatrix_c ABL, ElConstDistMatrix_c ABR, ElDistMatrix_c A20, ElDistMatrix_c A21, ElDistMatrix_c A22, 
  ElInt bsize );
ElError ElLockedRepartitionDownDiagonalDist_z
( ElConstDistMatrix_z ATL, ElConstDistMatrix_z ATR, ElDistMatrix_z A00, ElDistMatrix_z A01, ElDistMatrix_z A02,
                                                    ElDistMatrix_z A10, ElDistMatrix_z A11, ElDistMatrix_z A12,
  ElConstDistMatrix_z ABL, ElConstDistMatrix_z ABR, ElDistMatrix_z A20, ElDistMatrix_z A21, ElDistMatrix_z A22, 
  ElInt bsize );

/* Up a diagonal
   ============= */
ElError ElRepartitionUpDiagonal_i
( ElMatrix_i ATL, ElMatrix_i ATR, ElMatrix_i A00, ElMatrix_i A01, ElMatrix_i A02,
                                  ElMatrix_i A10, ElMatrix_i A11, ElMatrix_i A12,
  ElMatrix_i ABL, ElMatrix_i ABR, ElMatrix_i A20, ElMatrix_i A21, ElMatrix_i A22, 
  ElInt bsize );
ElError ElRepartitionUpDiagonal_s
( ElMatrix_s ATL, ElMatrix_s ATR, ElMatrix_s A00, ElMatrix_s A01, ElMatrix_s A02,
                                  ElMatrix_s A10, ElMatrix_s A11, ElMatrix_s A12,
  ElMatrix_s ABL, ElMatrix_s ABR, ElMatrix_s A20, ElMatrix_s A21, ElMatrix_s A22, 
  ElInt bsize );
ElError ElRepartitionUpDiagonal_d
( ElMatrix_d ATL, ElMatrix_d ATR, ElMatrix_d A00, ElMatrix_d A01, ElMatrix_d A02,
                                  ElMatrix_d A10, ElMatrix_d A11, ElMatrix_d A12,
  ElMatrix_d ABL, ElMatrix_d ABR, ElMatrix_d A20, ElMatrix_d A21, ElMatrix_d A22, 
  ElInt bsize );
ElError ElRepartitionUpDiagonal_c
( ElMatrix_c ATL, ElMatrix_c ATR, ElMatrix_c A00, ElMatrix_c A01, ElMatrix_c A02,
                                  ElMatrix_c A10, ElMatrix_c A11, ElMatrix_c A12,
  ElMatrix_c ABL, ElMatrix_c ABR, ElMatrix_c A20, ElMatrix_c A21, ElMatrix_c A22, 
  ElInt bsize );
ElError ElRepartitionUpDiagonal_z
( ElMatrix_z ATL, ElMatrix_z ATR, ElMatrix_z A00, ElMatrix_z A01, ElMatrix_z A02,
                                  ElMatrix_z A10, ElMatrix_z A11, ElMatrix_z A12,
  ElMatrix_z ABL, ElMatrix_z ABR, ElMatrix_z A20, ElMatrix_z A21, ElMatrix_z A22, 
  ElInt bsize );
ElError ElRepartitionUpDiagonalDist_i
( ElDistMatrix_i ATL, ElDistMatrix_i ATR, ElDistMatrix_i A00, ElDistMatrix_i A01, ElDistMatrix_i A02,
                                          ElDistMatrix_i A10, ElDistMatrix_i A11, ElDistMatrix_i A12,
  ElDistMatrix_i ABL, ElDistMatrix_i ABR, ElDistMatrix_i A20, ElDistMatrix_i A21, ElDistMatrix_i A22, 
  ElInt bsize );
ElError ElRepartitionUpDiagonalDist_s
( ElDistMatrix_s ATL, ElDistMatrix_s ATR, ElDistMatrix_s A00, ElDistMatrix_s A01, ElDistMatrix_s A02,
                                          ElDistMatrix_s A10, ElDistMatrix_s A11, ElDistMatrix_s A12,
  ElDistMatrix_s ABL, ElDistMatrix_s ABR, ElDistMatrix_s A20, ElDistMatrix_s A21, ElDistMatrix_s A22, 
  ElInt bsize );
ElError ElRepartitionUpDiagonalDist_d
( ElDistMatrix_d ATL, ElDistMatrix_d ATR, ElDistMatrix_d A00, ElDistMatrix_d A01, ElDistMatrix_d A02,
                                          ElDistMatrix_d A10, ElDistMatrix_d A11, ElDistMatrix_d A12,
  ElDistMatrix_d ABL, ElDistMatrix_d ABR, ElDistMatrix_d A20, ElDistMatrix_d A21, ElDistMatrix_d A22, 
  ElInt bsize );
ElError ElRepartitionUpDiagonalDist_c
( ElDistMatrix_c ATL, ElDistMatrix_c ATR, ElDistMatrix_c A00, ElDistMatrix_c A01, ElDistMatrix_c A02,
                                          ElDistMatrix_c A10, ElDistMatrix_c A11, ElDistMatrix_c A12,
  ElDistMatrix_c ABL, ElDistMatrix_c ABR, ElDistMatrix_c A20, ElDistMatrix_c A21, ElDistMatrix_c A22, 
  ElInt bsize );
ElError ElRepartitionUpDiagonalDist_z
( ElDistMatrix_z ATL, ElDistMatrix_z ATR, ElDistMatrix_z A00, ElDistMatrix_z A01, ElDistMatrix_z A02,
                                          ElDistMatrix_z A10, ElDistMatrix_z A11, ElDistMatrix_z A12,
  ElDistMatrix_z ABL, ElDistMatrix_z ABR, ElDistMatrix_z A20, ElDistMatrix_z A21, ElDistMatrix_z A22, 
  ElInt bsize );

ElError ElLockedRepartitionUpDiagonal_i
( ElConstMatrix_i ATL, ElConstMatrix_i ATR, ElMatrix_i A00, ElMatrix_i A01, ElMatrix_i A02,
                                            ElMatrix_i A10, ElMatrix_i A11, ElMatrix_i A12,
  ElConstMatrix_i ABL, ElConstMatrix_i ABR, ElMatrix_i A20, ElMatrix_i A21, ElMatrix_i A22, 
  ElInt bsize );
ElError ElLockedRepartitionUpDiagonal_s
( ElConstMatrix_s ATL, ElConstMatrix_s ATR, ElMatrix_s A00, ElMatrix_s A01, ElMatrix_s A02,
                                            ElMatrix_s A10, ElMatrix_s A11, ElMatrix_s A12,
  ElConstMatrix_s ABL, ElConstMatrix_s ABR, ElMatrix_s A20, ElMatrix_s A21, ElMatrix_s A22, 
  ElInt bsize );
ElError ElLockedRepartitionUpDiagonal_d
( ElConstMatrix_d ATL, ElConstMatrix_d ATR, ElMatrix_d A00, ElMatrix_d A01, ElMatrix_d A02,
                                            ElMatrix_d A10, ElMatrix_d A11, ElMatrix_d A12,
  ElConstMatrix_d ABL, ElConstMatrix_d ABR, ElMatrix_d A20, ElMatrix_d A21, ElMatrix_d A22, 
  ElInt bsize );
ElError ElLockedRepartitionUpDiagonal_c
( ElConstMatrix_c ATL, ElConstMatrix_c ATR, ElMatrix_c A00, ElMatrix_c A01, ElMatrix_c A02,
                                            ElMatrix_c A10, ElMatrix_c A11, ElMatrix_c A12,
  ElConstMatrix_c ABL, ElConstMatrix_c ABR, ElMatrix_c A20, ElMatrix_c A21, ElMatrix_c A22, 
  ElInt bsize );
ElError ElLockedRepartitionUpDiagonal_z
( ElConstMatrix_z ATL, ElConstMatrix_z ATR, ElMatrix_z A00, ElMatrix_z A01, ElMatrix_z A02,
                                            ElMatrix_z A10, ElMatrix_z A11, ElMatrix_z A12,
  ElConstMatrix_z ABL, ElConstMatrix_z ABR, ElMatrix_z A20, ElMatrix_z A21, ElMatrix_z A22, 
  ElInt bsize );
ElError ElLockedRepartitionUpDiagonalDist_i
( ElConstDistMatrix_i ATL, ElConstDistMatrix_i ATR, ElDistMatrix_i A00, ElDistMatrix_i A01, ElDistMatrix_i A02,
                                                    ElDistMatrix_i A10, ElDistMatrix_i A11, ElDistMatrix_i A12,
  ElConstDistMatrix_i ABL, ElConstDistMatrix_i ABR, ElDistMatrix_i A20, ElDistMatrix_i A21, ElDistMatrix_i A22, 
  ElInt bsize );
ElError ElLockedRepartitionUpDiagonalDist_s
( ElConstDistMatrix_s ATL, ElConstDistMatrix_s ATR, ElDistMatrix_s A00, ElDistMatrix_s A01, ElDistMatrix_s A02,
                                                    ElDistMatrix_s A10, ElDistMatrix_s A11, ElDistMatrix_s A12,
  ElConstDistMatrix_s ABL, ElConstDistMatrix_s ABR, ElDistMatrix_s A20, ElDistMatrix_s A21, ElDistMatrix_s A22, 
  ElInt bsize );
ElError ElLockedRepartitionUpDiagonalDist_d
( ElConstDistMatrix_d ATL, ElConstDistMatrix_d ATR, ElDistMatrix_d A00, ElDistMatrix_d A01, ElDistMatrix_d A02,
                                                    ElDistMatrix_d A10, ElDistMatrix_d A11, ElDistMatrix_d A12,
  ElConstDistMatrix_d ABL, ElConstDistMatrix_d ABR, ElDistMatrix_d A20, ElDistMatrix_d A21, ElDistMatrix_d A22, 
  ElInt bsize );
ElError ElLockedRepartitionUpDiagonalDist_c
( ElConstDistMatrix_c ATL, ElConstDistMatrix_c ATR, ElDistMatrix_c A00, ElDistMatrix_c A01, ElDistMatrix_c A02,
                                                    ElDistMatrix_c A10, ElDistMatrix_c A11, ElDistMatrix_c A12,
  ElConstDistMatrix_c ABL, ElConstDistMatrix_c ABR, ElDistMatrix_c A20, ElDistMatrix_c A21, ElDistMatrix_c A22, 
  ElInt bsize );
ElError ElLockedRepartitionUpDiagonalDist_z
( ElConstDistMatrix_z ATL, ElConstDistMatrix_z ATR, ElDistMatrix_z A00, ElDistMatrix_z A01, ElDistMatrix_z A02,
                                                    ElDistMatrix_z A10, ElDistMatrix_z A11, ElDistMatrix_z A12,
  ElConstDistMatrix_z ABL, ElConstDistMatrix_z ABR, ElDistMatrix_z A20, ElDistMatrix_z A21, ElDistMatrix_z A22, 
  ElInt bsize );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_FLAMEPART_REPARTITION_C_H */
