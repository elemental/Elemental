/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_FLAMEPART_PARTITION_C_H
#define EL_FLAMEPART_PARTITION_C_H

#ifdef __cplusplus
extern "C" {
#endif

/* Partition downwards from the top
   ================================ */
EL_EXPORT ElError ElPartitionDown_i
( ElMatrix_i A, ElMatrix_i AT, ElMatrix_i AB, ElInt heightAT );
EL_EXPORT ElError ElPartitionDown_s
( ElMatrix_s A, ElMatrix_s AT, ElMatrix_s AB, ElInt heightAT );
EL_EXPORT ElError ElPartitionDown_d
( ElMatrix_d A, ElMatrix_d AT, ElMatrix_d AB, ElInt heightAT );
EL_EXPORT ElError ElPartitionDown_c
( ElMatrix_c A, ElMatrix_c AT, ElMatrix_c AB, ElInt heightAT );
EL_EXPORT ElError ElPartitionDown_z
( ElMatrix_z A, ElMatrix_z AT, ElMatrix_z AB, ElInt heightAT );
EL_EXPORT ElError ElPartitionDownDist_i
( ElDistMatrix_i A, ElDistMatrix_i AT, ElDistMatrix_i AB, ElInt heightAT );
EL_EXPORT ElError ElPartitionDownDist_s
( ElDistMatrix_s A, ElDistMatrix_s AT, ElDistMatrix_s AB, ElInt heightAT );
EL_EXPORT ElError ElPartitionDownDist_d
( ElDistMatrix_d A, ElDistMatrix_d AT, ElDistMatrix_d AB, ElInt heightAT );
EL_EXPORT ElError ElPartitionDownDist_c
( ElDistMatrix_c A, ElDistMatrix_c AT, ElDistMatrix_c AB, ElInt heightAT );
EL_EXPORT ElError ElPartitionDownDist_z
( ElDistMatrix_z A, ElDistMatrix_z AT, ElDistMatrix_z AB, ElInt heightAT );

EL_EXPORT ElError ElLockedPartitionDown_i
( ElConstMatrix_i A, ElMatrix_i AT, ElMatrix_i AB, ElInt heightAT );
EL_EXPORT ElError ElLockedPartitionDown_s
( ElConstMatrix_s A, ElMatrix_s AT, ElMatrix_s AB, ElInt heightAT );
EL_EXPORT ElError ElLockedPartitionDown_d
( ElConstMatrix_d A, ElMatrix_d AT, ElMatrix_d AB, ElInt heightAT );
EL_EXPORT ElError ElLockedPartitionDown_c
( ElConstMatrix_c A, ElMatrix_c AT, ElMatrix_c AB, ElInt heightAT );
EL_EXPORT ElError ElLockedPartitionDown_z
( ElConstMatrix_z A, ElMatrix_z AT, ElMatrix_z AB, ElInt heightAT );
EL_EXPORT ElError ElLockedPartitionDownDist_i
( ElConstDistMatrix_i A, ElDistMatrix_i AT, ElDistMatrix_i AB, ElInt heightAT );
EL_EXPORT ElError ElLockedPartitionDownDist_s
( ElConstDistMatrix_s A, ElDistMatrix_s AT, ElDistMatrix_s AB, ElInt heightAT );
EL_EXPORT ElError ElLockedPartitionDownDist_d
( ElConstDistMatrix_d A, ElDistMatrix_d AT, ElDistMatrix_d AB, ElInt heightAT );
EL_EXPORT ElError ElLockedPartitionDownDist_c
( ElConstDistMatrix_c A, ElDistMatrix_c AT, ElDistMatrix_c AB, ElInt heightAT );
EL_EXPORT ElError ElLockedPartitionDownDist_z
( ElConstDistMatrix_z A, ElDistMatrix_z AT, ElDistMatrix_z AB, ElInt heightAT );

/* Partition upwards from the bottom
   ================================= */
EL_EXPORT ElError ElPartitionUp_i
( ElMatrix_i A, ElMatrix_i AT, ElMatrix_i AB, ElInt heightAB );
EL_EXPORT ElError ElPartitionUp_s
( ElMatrix_s A, ElMatrix_s AT, ElMatrix_s AB, ElInt heightAB );
EL_EXPORT ElError ElPartitionUp_d
( ElMatrix_d A, ElMatrix_d AT, ElMatrix_d AB, ElInt heightAB );
EL_EXPORT ElError ElPartitionUp_c
( ElMatrix_c A, ElMatrix_c AT, ElMatrix_c AB, ElInt heightAB );
EL_EXPORT ElError ElPartitionUp_z
( ElMatrix_z A, ElMatrix_z AT, ElMatrix_z AB, ElInt heightAB );
EL_EXPORT ElError ElPartitionUpDist_i
( ElDistMatrix_i A, ElDistMatrix_i AT, ElDistMatrix_i AB, ElInt heightAB );
EL_EXPORT ElError ElPartitionUpDist_s
( ElDistMatrix_s A, ElDistMatrix_s AT, ElDistMatrix_s AB, ElInt heightAB );
EL_EXPORT ElError ElPartitionUpDist_d
( ElDistMatrix_d A, ElDistMatrix_d AT, ElDistMatrix_d AB, ElInt heightAB );
EL_EXPORT ElError ElPartitionUpDist_c
( ElDistMatrix_c A, ElDistMatrix_c AT, ElDistMatrix_c AB, ElInt heightAB );
EL_EXPORT ElError ElPartitionUpDist_z
( ElDistMatrix_z A, ElDistMatrix_z AT, ElDistMatrix_z AB, ElInt heightAB );

EL_EXPORT ElError ElLockedPartitionUp_i
( ElConstMatrix_i A, ElMatrix_i AT, ElMatrix_i AB, ElInt heightAB );
EL_EXPORT ElError ElLockedPartitionUp_s
( ElConstMatrix_s A, ElMatrix_s AT, ElMatrix_s AB, ElInt heightAB );
EL_EXPORT ElError ElLockedPartitionUp_d
( ElConstMatrix_d A, ElMatrix_d AT, ElMatrix_d AB, ElInt heightAB );
EL_EXPORT ElError ElLockedPartitionUp_c
( ElConstMatrix_c A, ElMatrix_c AT, ElMatrix_c AB, ElInt heightAB );
EL_EXPORT ElError ElLockedPartitionUp_z
( ElConstMatrix_z A, ElMatrix_z AT, ElMatrix_z AB, ElInt heightAB );
EL_EXPORT ElError ElLockedPartitionUpDist_i
( ElConstDistMatrix_i A, ElDistMatrix_i AT, ElDistMatrix_i AB, ElInt heightAB );
EL_EXPORT ElError ElLockedPartitionUpDist_s
( ElConstDistMatrix_s A, ElDistMatrix_s AT, ElDistMatrix_s AB, ElInt heightAB );
EL_EXPORT ElError ElLockedPartitionUpDist_d
( ElConstDistMatrix_d A, ElDistMatrix_d AT, ElDistMatrix_d AB, ElInt heightAB );
EL_EXPORT ElError ElLockedPartitionUpDist_c
( ElConstDistMatrix_c A, ElDistMatrix_c AT, ElDistMatrix_c AB, ElInt heightAB );
EL_EXPORT ElError ElLockedPartitionUpDist_z
( ElConstDistMatrix_z A, ElDistMatrix_z AT, ElDistMatrix_z AB, ElInt heightAB );

/* Partition leftward from the right
   ================================= */
EL_EXPORT ElError ElPartitionLeft_i
( ElMatrix_i A, ElMatrix_i AL, ElMatrix_i AR, ElInt widthAR );
EL_EXPORT ElError ElPartitionLeft_s
( ElMatrix_s A, ElMatrix_s AL, ElMatrix_s AR, ElInt widthAR );
EL_EXPORT ElError ElPartitionLeft_d
( ElMatrix_d A, ElMatrix_d AL, ElMatrix_d AR, ElInt widthAR );
EL_EXPORT ElError ElPartitionLeft_c
( ElMatrix_c A, ElMatrix_c AL, ElMatrix_c AR, ElInt widthAR );
EL_EXPORT ElError ElPartitionLeft_z
( ElMatrix_z A, ElMatrix_z AL, ElMatrix_z AR, ElInt widthAR );
EL_EXPORT ElError ElPartitionLeftDist_i
( ElDistMatrix_i A, ElDistMatrix_i AL, ElDistMatrix_i AR, ElInt widthAR );
EL_EXPORT ElError ElPartitionLeftDist_s
( ElDistMatrix_s A, ElDistMatrix_s AL, ElDistMatrix_s AR, ElInt widthAR );
EL_EXPORT ElError ElPartitionLeftDist_d
( ElDistMatrix_d A, ElDistMatrix_d AL, ElDistMatrix_d AR, ElInt widthAR );
EL_EXPORT ElError ElPartitionLeftDist_c
( ElDistMatrix_c A, ElDistMatrix_c AL, ElDistMatrix_c AR, ElInt widthAR );
EL_EXPORT ElError ElPartitionLeftDist_z
( ElDistMatrix_z A, ElDistMatrix_z AL, ElDistMatrix_z AR, ElInt widthAR );

EL_EXPORT ElError ElLockedPartitionLeft_i
( ElConstMatrix_i A, ElMatrix_i AL, ElMatrix_i AR, ElInt widthAR );
EL_EXPORT ElError ElLockedPartitionLeft_s
( ElConstMatrix_s A, ElMatrix_s AL, ElMatrix_s AR, ElInt widthAR );
EL_EXPORT ElError ElLockedPartitionLeft_d
( ElConstMatrix_d A, ElMatrix_d AL, ElMatrix_d AR, ElInt widthAR );
EL_EXPORT ElError ElLockedPartitionLeft_c
( ElConstMatrix_c A, ElMatrix_c AL, ElMatrix_c AR, ElInt widthAR );
EL_EXPORT ElError ElLockedPartitionLeft_z
( ElConstMatrix_z A, ElMatrix_z AL, ElMatrix_z AR, ElInt widthAR );
EL_EXPORT ElError ElLockedPartitionLeftDist_i
( ElConstDistMatrix_i A, ElDistMatrix_i AL, ElDistMatrix_i AR, ElInt widthAR );
EL_EXPORT ElError ElLockedPartitionLeftDist_s
( ElConstDistMatrix_s A, ElDistMatrix_s AL, ElDistMatrix_s AR, ElInt widthAR );
EL_EXPORT ElError ElLockedPartitionLeftDist_d
( ElConstDistMatrix_d A, ElDistMatrix_d AL, ElDistMatrix_d AR, ElInt widthAR );
EL_EXPORT ElError ElLockedPartitionLeftDist_c
( ElConstDistMatrix_c A, ElDistMatrix_c AL, ElDistMatrix_c AR, ElInt widthAR );
EL_EXPORT ElError ElLockedPartitionLeftDist_z
( ElConstDistMatrix_z A, ElDistMatrix_z AL, ElDistMatrix_z AR, ElInt widthAR );

/* Partition rightward from the left
   ================================= */
EL_EXPORT ElError ElPartitionRight_i
( ElMatrix_i A, ElMatrix_i AL, ElMatrix_i AR, ElInt widthAL );
EL_EXPORT ElError ElPartitionRight_s
( ElMatrix_s A, ElMatrix_s AL, ElMatrix_s AR, ElInt widthAL );
EL_EXPORT ElError ElPartitionRight_d
( ElMatrix_d A, ElMatrix_d AL, ElMatrix_d AR, ElInt widthAL );
EL_EXPORT ElError ElPartitionRight_c
( ElMatrix_c A, ElMatrix_c AL, ElMatrix_c AR, ElInt widthAL );
EL_EXPORT ElError ElPartitionRight_z
( ElMatrix_z A, ElMatrix_z AL, ElMatrix_z AR, ElInt widthAL );
EL_EXPORT ElError ElPartitionRightDist_i
( ElDistMatrix_i A, ElDistMatrix_i AL, ElDistMatrix_i AR, ElInt widthAL );
EL_EXPORT ElError ElPartitionRightDist_s
( ElDistMatrix_s A, ElDistMatrix_s AL, ElDistMatrix_s AR, ElInt widthAL );
EL_EXPORT ElError ElPartitionRightDist_d
( ElDistMatrix_d A, ElDistMatrix_d AL, ElDistMatrix_d AR, ElInt widthAL );
EL_EXPORT ElError ElPartitionRightDist_c
( ElDistMatrix_c A, ElDistMatrix_c AL, ElDistMatrix_c AR, ElInt widthAL );
EL_EXPORT ElError ElPartitionRightDist_z
( ElDistMatrix_z A, ElDistMatrix_z AL, ElDistMatrix_z AR, ElInt widthAL );

EL_EXPORT ElError ElLockedPartitionRight_i
( ElConstMatrix_i A, ElMatrix_i AL, ElMatrix_i AR, ElInt widthAL );
EL_EXPORT ElError ElLockedPartitionRight_s
( ElConstMatrix_s A, ElMatrix_s AL, ElMatrix_s AR, ElInt widthAL );
EL_EXPORT ElError ElLockedPartitionRight_d
( ElConstMatrix_d A, ElMatrix_d AL, ElMatrix_d AR, ElInt widthAL );
EL_EXPORT ElError ElLockedPartitionRight_c
( ElConstMatrix_c A, ElMatrix_c AL, ElMatrix_c AR, ElInt widthAL );
EL_EXPORT ElError ElLockedPartitionRight_z
( ElConstMatrix_z A, ElMatrix_z AL, ElMatrix_z AR, ElInt widthAL );
EL_EXPORT ElError ElLockedPartitionRightDist_i
( ElConstDistMatrix_i A, ElDistMatrix_i AL, ElDistMatrix_i AR, ElInt widthAL );
EL_EXPORT ElError ElLockedPartitionRightDist_s
( ElConstDistMatrix_s A, ElDistMatrix_s AL, ElDistMatrix_s AR, ElInt widthAL );
EL_EXPORT ElError ElLockedPartitionRightDist_d
( ElConstDistMatrix_d A, ElDistMatrix_d AL, ElDistMatrix_d AR, ElInt widthAL );
EL_EXPORT ElError ElLockedPartitionRightDist_c
( ElConstDistMatrix_c A, ElDistMatrix_c AL, ElDistMatrix_c AR, ElInt widthAL );
EL_EXPORT ElError ElLockedPartitionRightDist_z
( ElConstDistMatrix_z A, ElDistMatrix_z AL, ElDistMatrix_z AR, ElInt widthAL );

/* Down the main diagonal
   ====================== */
EL_EXPORT ElError ElPartitionDownDiagonal_i
( ElMatrix_i A, 
  ElMatrix_i ATL, ElMatrix_i ATR, 
  ElMatrix_i ABL, ElMatrix_i ABR, ElInt diagDist );
EL_EXPORT ElError ElPartitionDownDiagonal_s
( ElMatrix_s A, 
  ElMatrix_s ATL, ElMatrix_s ATR, 
  ElMatrix_s ABL, ElMatrix_s ABR, ElInt diagDist );
EL_EXPORT ElError ElPartitionDownDiagonal_d
( ElMatrix_d A, 
  ElMatrix_d ATL, ElMatrix_d ATR, 
  ElMatrix_d ABL, ElMatrix_d ABR, ElInt diagDist );
EL_EXPORT ElError ElPartitionDownDiagonal_c
( ElMatrix_c A, 
  ElMatrix_c ATL, ElMatrix_c ATR, 
  ElMatrix_c ABL, ElMatrix_c ABR, ElInt diagDist );
EL_EXPORT ElError ElPartitionDownDiagonal_z
( ElMatrix_z A, 
  ElMatrix_z ATL, ElMatrix_z ATR, 
  ElMatrix_z ABL, ElMatrix_z ABR, ElInt diagDist );
EL_EXPORT ElError ElPartitionDownDiagonalDist_i
( ElDistMatrix_i A, 
  ElDistMatrix_i ATL, ElDistMatrix_i ATR, 
  ElDistMatrix_i ABL, ElDistMatrix_i ABR, ElInt diagDist );
EL_EXPORT ElError ElPartitionDownDiagonalDist_s
( ElDistMatrix_s A, 
  ElDistMatrix_s ATL, ElDistMatrix_s ATR, 
  ElDistMatrix_s ABL, ElDistMatrix_s ABR, ElInt diagDist );
EL_EXPORT ElError ElPartitionDownDiagonalDist_d
( ElDistMatrix_d A, 
  ElDistMatrix_d ATL, ElDistMatrix_d ATR, 
  ElDistMatrix_d ABL, ElDistMatrix_d ABR, ElInt diagDist );
EL_EXPORT ElError ElPartitionDownDiagonalDist_c
( ElDistMatrix_c A, 
  ElDistMatrix_c ATL, ElDistMatrix_c ATR, 
  ElDistMatrix_c ABL, ElDistMatrix_c ABR, ElInt diagDist );
EL_EXPORT ElError ElPartitionDownDiagonalDist_z
( ElDistMatrix_z A, 
  ElDistMatrix_z ATL, ElDistMatrix_z ATR, 
  ElDistMatrix_z ABL, ElDistMatrix_z ABR, ElInt diagDist );

EL_EXPORT ElError ElLockedPartitionDownDiagonal_i
( ElConstMatrix_i A, 
  ElMatrix_i ATL, ElMatrix_i ATR, 
  ElMatrix_i ABL, ElMatrix_i ABR, ElInt diagDist );
EL_EXPORT ElError ElLockedPartitionDownDiagonal_s
( ElConstMatrix_s A, 
  ElMatrix_s ATL, ElMatrix_s ATR, 
  ElMatrix_s ABL, ElMatrix_s ABR, ElInt diagDist );
EL_EXPORT ElError ElLockedPartitionDownDiagonal_d
( ElConstMatrix_d A, 
  ElMatrix_d ATL, ElMatrix_d ATR, 
  ElMatrix_d ABL, ElMatrix_d ABR, ElInt diagDist );
EL_EXPORT ElError ElLockedPartitionDownDiagonal_c
( ElConstMatrix_c A, 
  ElMatrix_c ATL, ElMatrix_c ATR, 
  ElMatrix_c ABL, ElMatrix_c ABR, ElInt diagDist );
EL_EXPORT ElError ElLockedPartitionDownDiagonal_z
( ElConstMatrix_z A, 
  ElMatrix_z ATL, ElMatrix_z ATR, 
  ElMatrix_z ABL, ElMatrix_z ABR, ElInt diagDist );
EL_EXPORT ElError ElLockedPartitionDownDiagonalDist_i
( ElConstDistMatrix_i A, 
  ElDistMatrix_i ATL, ElDistMatrix_i ATR, 
  ElDistMatrix_i ABL, ElDistMatrix_i ABR, ElInt diagDist );
EL_EXPORT ElError ElLockedPartitionDownDiagonalDist_s
( ElConstDistMatrix_s A, 
  ElDistMatrix_s ATL, ElDistMatrix_s ATR, 
  ElDistMatrix_s ABL, ElDistMatrix_s ABR, ElInt diagDist );
EL_EXPORT ElError ElLockedPartitionDownDiagonalDist_d
( ElConstDistMatrix_d A, 
  ElDistMatrix_d ATL, ElDistMatrix_d ATR, 
  ElDistMatrix_d ABL, ElDistMatrix_d ABR, ElInt diagDist );
EL_EXPORT ElError ElLockedPartitionDownDiagonalDist_c
( ElConstDistMatrix_c A, 
  ElDistMatrix_c ATL, ElDistMatrix_c ATR, 
  ElDistMatrix_c ABL, ElDistMatrix_c ABR, ElInt diagDist );
EL_EXPORT ElError ElLockedPartitionDownDiagonalDist_z
( ElConstDistMatrix_z A, 
  ElDistMatrix_z ATL, ElDistMatrix_z ATR, 
  ElDistMatrix_z ABL, ElDistMatrix_z ABR, ElInt diagDist );

/* Down an offset diagonal
   ======================= */
EL_EXPORT ElError ElPartitionDownOffsetDiagonal_i
( ElInt offset, ElMatrix_i A, 
  ElMatrix_i ATL, ElMatrix_i ATR, 
  ElMatrix_i ABL, ElMatrix_i ABR, ElInt diagDist );
EL_EXPORT ElError ElPartitionDownOffsetDiagonal_s
( ElInt offset, ElMatrix_s A, 
  ElMatrix_s ATL, ElMatrix_s ATR, 
  ElMatrix_s ABL, ElMatrix_s ABR, ElInt diagDist );
EL_EXPORT ElError ElPartitionDownOffsetDiagonal_d
( ElInt offset, ElMatrix_d A, 
  ElMatrix_d ATL, ElMatrix_d ATR, 
  ElMatrix_d ABL, ElMatrix_d ABR, ElInt diagDist );
EL_EXPORT ElError ElPartitionDownOffsetDiagonal_c
( ElInt offset, ElMatrix_c A, 
  ElMatrix_c ATL, ElMatrix_c ATR, 
  ElMatrix_c ABL, ElMatrix_c ABR, ElInt diagDist );
EL_EXPORT ElError ElPartitionDownOffsetDiagonal_z
( ElInt offset, ElMatrix_z A, 
  ElMatrix_z ATL, ElMatrix_z ATR, 
  ElMatrix_z ABL, ElMatrix_z ABR, ElInt diagDist );
EL_EXPORT ElError ElPartitionDownOffsetDiagonalDist_i
( ElInt offset, ElDistMatrix_i A, 
  ElDistMatrix_i ATL, ElDistMatrix_i ATR, 
  ElDistMatrix_i ABL, ElDistMatrix_i ABR, ElInt diagDist );
EL_EXPORT ElError ElPartitionDownOffsetDiagonalDist_s
( ElInt offset, ElDistMatrix_s A, 
  ElDistMatrix_s ATL, ElDistMatrix_s ATR, 
  ElDistMatrix_s ABL, ElDistMatrix_s ABR, ElInt diagDist );
EL_EXPORT ElError ElPartitionDownOffsetDiagonalDist_d
( ElInt offset, ElDistMatrix_d A, 
  ElDistMatrix_d ATL, ElDistMatrix_d ATR, 
  ElDistMatrix_d ABL, ElDistMatrix_d ABR, ElInt diagDist );
EL_EXPORT ElError ElPartitionDownOffsetDiagonalDist_c
( ElInt offset, ElDistMatrix_c A, 
  ElDistMatrix_c ATL, ElDistMatrix_c ATR, 
  ElDistMatrix_c ABL, ElDistMatrix_c ABR, ElInt diagDist );
EL_EXPORT ElError ElPartitionDownOffsetDiagonalDist_z
( ElInt offset, ElDistMatrix_z A, 
  ElDistMatrix_z ATL, ElDistMatrix_z ATR, 
  ElDistMatrix_z ABL, ElDistMatrix_z ABR, ElInt diagDist );

EL_EXPORT ElError ElLockedPartitionDownOffsetDiagonal_i
( ElInt offset, ElConstMatrix_i A, 
  ElMatrix_i ATL, ElMatrix_i ATR, 
  ElMatrix_i ABL, ElMatrix_i ABR, ElInt diagDist );
EL_EXPORT ElError ElLockedPartitionDownOffsetDiagonal_s
( ElInt offset, ElConstMatrix_s A, 
  ElMatrix_s ATL, ElMatrix_s ATR, 
  ElMatrix_s ABL, ElMatrix_s ABR, ElInt diagDist );
EL_EXPORT ElError ElLockedPartitionDownOffsetDiagonal_d
( ElInt offset, ElConstMatrix_d A, 
  ElMatrix_d ATL, ElMatrix_d ATR, 
  ElMatrix_d ABL, ElMatrix_d ABR, ElInt diagDist );
EL_EXPORT ElError ElLockedPartitionDownOffsetDiagonal_c
( ElInt offset, ElConstMatrix_c A, 
  ElMatrix_c ATL, ElMatrix_c ATR, 
  ElMatrix_c ABL, ElMatrix_c ABR, ElInt diagDist );
EL_EXPORT ElError ElLockedPartitionDownOffsetDiagonal_z
( ElInt offset, ElConstMatrix_z A, 
  ElMatrix_z ATL, ElMatrix_z ATR, 
  ElMatrix_z ABL, ElMatrix_z ABR, ElInt diagDist );
EL_EXPORT ElError ElLockedPartitionDownOffsetDiagonalDist_i
( ElInt offset, ElConstDistMatrix_i A, 
  ElDistMatrix_i ATL, ElDistMatrix_i ATR, 
  ElDistMatrix_i ABL, ElDistMatrix_i ABR, ElInt diagDist );
EL_EXPORT ElError ElLockedPartitionDownOffsetDiagonalDist_s
( ElInt offset, ElConstDistMatrix_s A, 
  ElDistMatrix_s ATL, ElDistMatrix_s ATR, 
  ElDistMatrix_s ABL, ElDistMatrix_s ABR, ElInt diagDist );
EL_EXPORT ElError ElLockedPartitionDownOffsetDiagonalDist_d
( ElInt offset, ElConstDistMatrix_d A, 
  ElDistMatrix_d ATL, ElDistMatrix_d ATR, 
  ElDistMatrix_d ABL, ElDistMatrix_d ABR, ElInt diagDist );
EL_EXPORT ElError ElLockedPartitionDownOffsetDiagonalDist_c
( ElInt offset, ElConstDistMatrix_c A, 
  ElDistMatrix_c ATL, ElDistMatrix_c ATR, 
  ElDistMatrix_c ABL, ElDistMatrix_c ABR, ElInt diagDist );
EL_EXPORT ElError ElLockedPartitionDownOffsetDiagonalDist_z
( ElInt offset, ElConstDistMatrix_z A, 
  ElDistMatrix_z ATL, ElDistMatrix_z ATR, 
  ElDistMatrix_z ABL, ElDistMatrix_z ABR, ElInt diagDist );

/* Up the main diagonal
   ==================== */
EL_EXPORT ElError ElPartitionUpDiagonal_i
( ElMatrix_i A, 
  ElMatrix_i ATL, ElMatrix_i ATR, 
  ElMatrix_i ABL, ElMatrix_i ABR, ElInt diagDist );
EL_EXPORT ElError ElPartitionUpDiagonal_s
( ElMatrix_s A, 
  ElMatrix_s ATL, ElMatrix_s ATR, 
  ElMatrix_s ABL, ElMatrix_s ABR, ElInt diagDist );
EL_EXPORT ElError ElPartitionUpDiagonal_d
( ElMatrix_d A, 
  ElMatrix_d ATL, ElMatrix_d ATR, 
  ElMatrix_d ABL, ElMatrix_d ABR, ElInt diagDist );
EL_EXPORT ElError ElPartitionUpDiagonal_c
( ElMatrix_c A, 
  ElMatrix_c ATL, ElMatrix_c ATR, 
  ElMatrix_c ABL, ElMatrix_c ABR, ElInt diagDist );
EL_EXPORT ElError ElPartitionUpDiagonal_z
( ElMatrix_z A, 
  ElMatrix_z ATL, ElMatrix_z ATR, 
  ElMatrix_z ABL, ElMatrix_z ABR, ElInt diagDist );
EL_EXPORT ElError ElPartitionUpDiagonalDist_i
( ElDistMatrix_i A, 
  ElDistMatrix_i ATL, ElDistMatrix_i ATR, 
  ElDistMatrix_i ABL, ElDistMatrix_i ABR, ElInt diagDist );
EL_EXPORT ElError ElPartitionUpDiagonalDist_s
( ElDistMatrix_s A, 
  ElDistMatrix_s ATL, ElDistMatrix_s ATR, 
  ElDistMatrix_s ABL, ElDistMatrix_s ABR, ElInt diagDist );
EL_EXPORT ElError ElPartitionUpDiagonalDist_d
( ElDistMatrix_d A, 
  ElDistMatrix_d ATL, ElDistMatrix_d ATR, 
  ElDistMatrix_d ABL, ElDistMatrix_d ABR, ElInt diagDist );
EL_EXPORT ElError ElPartitionUpDiagonalDist_c
( ElDistMatrix_c A, 
  ElDistMatrix_c ATL, ElDistMatrix_c ATR, 
  ElDistMatrix_c ABL, ElDistMatrix_c ABR, ElInt diagDist );
EL_EXPORT ElError ElPartitionUpDiagonalDist_z
( ElDistMatrix_z A, 
  ElDistMatrix_z ATL, ElDistMatrix_z ATR, 
  ElDistMatrix_z ABL, ElDistMatrix_z ABR, ElInt diagDist );

EL_EXPORT ElError ElLockedPartitionUpDiagonal_i
( ElConstMatrix_i A, 
  ElMatrix_i ATL, ElMatrix_i ATR, 
  ElMatrix_i ABL, ElMatrix_i ABR, ElInt diagDist );
EL_EXPORT ElError ElLockedPartitionUpDiagonal_s
( ElConstMatrix_s A, 
  ElMatrix_s ATL, ElMatrix_s ATR, 
  ElMatrix_s ABL, ElMatrix_s ABR, ElInt diagDist );
EL_EXPORT ElError ElLockedPartitionUpDiagonal_d
( ElConstMatrix_d A, 
  ElMatrix_d ATL, ElMatrix_d ATR, 
  ElMatrix_d ABL, ElMatrix_d ABR, ElInt diagDist );
EL_EXPORT ElError ElLockedPartitionUpDiagonal_c
( ElConstMatrix_c A, 
  ElMatrix_c ATL, ElMatrix_c ATR, 
  ElMatrix_c ABL, ElMatrix_c ABR, ElInt diagDist );
EL_EXPORT ElError ElLockedPartitionUpDiagonal_z
( ElConstMatrix_z A, 
  ElMatrix_z ATL, ElMatrix_z ATR, 
  ElMatrix_z ABL, ElMatrix_z ABR, ElInt diagDist );
EL_EXPORT ElError ElLockedPartitionUpDiagonalDist_i
( ElConstDistMatrix_i A, 
  ElDistMatrix_i ATL, ElDistMatrix_i ATR, 
  ElDistMatrix_i ABL, ElDistMatrix_i ABR, ElInt diagDist );
EL_EXPORT ElError ElLockedPartitionUpDiagonalDist_s
( ElConstDistMatrix_s A, 
  ElDistMatrix_s ATL, ElDistMatrix_s ATR, 
  ElDistMatrix_s ABL, ElDistMatrix_s ABR, ElInt diagDist );
EL_EXPORT ElError ElLockedPartitionUpDiagonalDist_d
( ElConstDistMatrix_d A, 
  ElDistMatrix_d ATL, ElDistMatrix_d ATR, 
  ElDistMatrix_d ABL, ElDistMatrix_d ABR, ElInt diagDist );
EL_EXPORT ElError ElLockedPartitionUpDiagonalDist_c
( ElConstDistMatrix_c A, 
  ElDistMatrix_c ATL, ElDistMatrix_c ATR, 
  ElDistMatrix_c ABL, ElDistMatrix_c ABR, ElInt diagDist );
EL_EXPORT ElError ElLockedPartitionUpDiagonalDist_z
( ElConstDistMatrix_z A, 
  ElDistMatrix_z ATL, ElDistMatrix_z ATR, 
  ElDistMatrix_z ABL, ElDistMatrix_z ABR, ElInt diagDist );

/* Up an offset diagonal
   ===================== */
EL_EXPORT ElError ElPartitionUpOffsetDiagonal_i
( ElInt offset, ElMatrix_i A, 
  ElMatrix_i ATL, ElMatrix_i ATR, 
  ElMatrix_i ABL, ElMatrix_i ABR, ElInt diagDist );
EL_EXPORT ElError ElPartitionUpOffsetDiagonal_s
( ElInt offset, ElMatrix_s A, 
  ElMatrix_s ATL, ElMatrix_s ATR, 
  ElMatrix_s ABL, ElMatrix_s ABR, ElInt diagDist );
EL_EXPORT ElError ElPartitionUpOffsetDiagonal_d
( ElInt offset, ElMatrix_d A, 
  ElMatrix_d ATL, ElMatrix_d ATR, 
  ElMatrix_d ABL, ElMatrix_d ABR, ElInt diagDist );
EL_EXPORT ElError ElPartitionUpOffsetDiagonal_c
( ElInt offset, ElMatrix_c A, 
  ElMatrix_c ATL, ElMatrix_c ATR, 
  ElMatrix_c ABL, ElMatrix_c ABR, ElInt diagDist );
EL_EXPORT ElError ElPartitionUpOffsetDiagonal_z
( ElInt offset, ElMatrix_z A, 
  ElMatrix_z ATL, ElMatrix_z ATR, 
  ElMatrix_z ABL, ElMatrix_z ABR, ElInt diagDist );
EL_EXPORT ElError ElPartitionUpOffsetDiagonalDist_i
( ElInt offset, ElDistMatrix_i A, 
  ElDistMatrix_i ATL, ElDistMatrix_i ATR, 
  ElDistMatrix_i ABL, ElDistMatrix_i ABR, ElInt diagDist );
EL_EXPORT ElError ElPartitionUpOffsetDiagonalDist_s
( ElInt offset, ElDistMatrix_s A, 
  ElDistMatrix_s ATL, ElDistMatrix_s ATR, 
  ElDistMatrix_s ABL, ElDistMatrix_s ABR, ElInt diagDist );
EL_EXPORT ElError ElPartitionUpOffsetDiagonalDist_d
( ElInt offset, ElDistMatrix_d A, 
  ElDistMatrix_d ATL, ElDistMatrix_d ATR, 
  ElDistMatrix_d ABL, ElDistMatrix_d ABR, ElInt diagDist );
EL_EXPORT ElError ElPartitionUpOffsetDiagonalDist_c
( ElInt offset, ElDistMatrix_c A, 
  ElDistMatrix_c ATL, ElDistMatrix_c ATR, 
  ElDistMatrix_c ABL, ElDistMatrix_c ABR, ElInt diagDist );
EL_EXPORT ElError ElPartitionUpOffsetDiagonalDist_z
( ElInt offset, ElDistMatrix_z A, 
  ElDistMatrix_z ATL, ElDistMatrix_z ATR, 
  ElDistMatrix_z ABL, ElDistMatrix_z ABR, ElInt diagDist );

EL_EXPORT ElError ElLockedPartitionUpOffsetDiagonal_i
( ElInt offset, ElConstMatrix_i A, 
  ElMatrix_i ATL, ElMatrix_i ATR, 
  ElMatrix_i ABL, ElMatrix_i ABR, ElInt diagDist );
EL_EXPORT ElError ElLockedPartitionUpOffsetDiagonal_s
( ElInt offset, ElConstMatrix_s A, 
  ElMatrix_s ATL, ElMatrix_s ATR, 
  ElMatrix_s ABL, ElMatrix_s ABR, ElInt diagDist );
EL_EXPORT ElError ElLockedPartitionUpOffsetDiagonal_d
( ElInt offset, ElConstMatrix_d A, 
  ElMatrix_d ATL, ElMatrix_d ATR, 
  ElMatrix_d ABL, ElMatrix_d ABR, ElInt diagDist );
EL_EXPORT ElError ElLockedPartitionUpOffsetDiagonal_c
( ElInt offset, ElConstMatrix_c A, 
  ElMatrix_c ATL, ElMatrix_c ATR, 
  ElMatrix_c ABL, ElMatrix_c ABR, ElInt diagDist );
EL_EXPORT ElError ElLockedPartitionUpOffsetDiagonal_z
( ElInt offset, ElConstMatrix_z A, 
  ElMatrix_z ATL, ElMatrix_z ATR, 
  ElMatrix_z ABL, ElMatrix_z ABR, ElInt diagDist );
EL_EXPORT ElError ElLockedPartitionUpOffsetDiagonalDist_i
( ElInt offset, ElConstDistMatrix_i A, 
  ElDistMatrix_i ATL, ElDistMatrix_i ATR, 
  ElDistMatrix_i ABL, ElDistMatrix_i ABR, ElInt diagDist );
EL_EXPORT ElError ElLockedPartitionUpOffsetDiagonalDist_s
( ElInt offset, ElConstDistMatrix_s A, 
  ElDistMatrix_s ATL, ElDistMatrix_s ATR, 
  ElDistMatrix_s ABL, ElDistMatrix_s ABR, ElInt diagDist );
EL_EXPORT ElError ElLockedPartitionUpOffsetDiagonalDist_d
( ElInt offset, ElConstDistMatrix_d A, 
  ElDistMatrix_d ATL, ElDistMatrix_d ATR, 
  ElDistMatrix_d ABL, ElDistMatrix_d ABR, ElInt diagDist );
EL_EXPORT ElError ElLockedPartitionUpOffsetDiagonalDist_c
( ElInt offset, ElConstDistMatrix_c A, 
  ElDistMatrix_c ATL, ElDistMatrix_c ATR, 
  ElDistMatrix_c ABL, ElDistMatrix_c ABR, ElInt diagDist );
EL_EXPORT ElError ElLockedPartitionUpOffsetDiagonalDist_z
( ElInt offset, ElConstDistMatrix_z A, 
  ElDistMatrix_z ATL, ElDistMatrix_z ATR, 
  ElDistMatrix_z ABL, ElDistMatrix_z ABR, ElInt diagDist );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_FLAMEPART_PARTITION_C_H */
