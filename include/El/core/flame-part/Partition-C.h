/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_FLAMEPART_PARTITION_C_H
#define EL_FLAMEPART_PARTITION_C_H

#ifdef __cplusplus
extern "C" {
#endif

/* Partition downwards from the top
   ================================ */
ElError ElPartitionDown_i
( ElMatrix_i A, ElMatrix_i AT, ElMatrix_i AB, ElInt heightAT );
ElError ElPartitionDown_s
( ElMatrix_s A, ElMatrix_s AT, ElMatrix_s AB, ElInt heightAT );
ElError ElPartitionDown_d
( ElMatrix_d A, ElMatrix_d AT, ElMatrix_d AB, ElInt heightAT );
ElError ElPartitionDown_c
( ElMatrix_c A, ElMatrix_c AT, ElMatrix_c AB, ElInt heightAT );
ElError ElPartitionDown_z
( ElMatrix_z A, ElMatrix_z AT, ElMatrix_z AB, ElInt heightAT );
ElError ElPartitionDownDist_i
( ElDistMatrix_i A, ElDistMatrix_i AT, ElDistMatrix_i AB, ElInt heightAT );
ElError ElPartitionDownDist_s
( ElDistMatrix_s A, ElDistMatrix_s AT, ElDistMatrix_s AB, ElInt heightAT );
ElError ElPartitionDownDist_d
( ElDistMatrix_d A, ElDistMatrix_d AT, ElDistMatrix_d AB, ElInt heightAT );
ElError ElPartitionDownDist_c
( ElDistMatrix_c A, ElDistMatrix_c AT, ElDistMatrix_c AB, ElInt heightAT );
ElError ElPartitionDownDist_z
( ElDistMatrix_z A, ElDistMatrix_z AT, ElDistMatrix_z AB, ElInt heightAT );

ElError ElLockedPartitionDown_i
( ElConstMatrix_i A, ElMatrix_i AT, ElMatrix_i AB, ElInt heightAT );
ElError ElLockedPartitionDown_s
( ElConstMatrix_s A, ElMatrix_s AT, ElMatrix_s AB, ElInt heightAT );
ElError ElLockedPartitionDown_d
( ElConstMatrix_d A, ElMatrix_d AT, ElMatrix_d AB, ElInt heightAT );
ElError ElLockedPartitionDown_c
( ElConstMatrix_c A, ElMatrix_c AT, ElMatrix_c AB, ElInt heightAT );
ElError ElLockedPartitionDown_z
( ElConstMatrix_z A, ElMatrix_z AT, ElMatrix_z AB, ElInt heightAT );
ElError ElLockedPartitionDownDist_i
( ElConstDistMatrix_i A, ElDistMatrix_i AT, ElDistMatrix_i AB, ElInt heightAT );
ElError ElLockedPartitionDownDist_s
( ElConstDistMatrix_s A, ElDistMatrix_s AT, ElDistMatrix_s AB, ElInt heightAT );
ElError ElLockedPartitionDownDist_d
( ElConstDistMatrix_d A, ElDistMatrix_d AT, ElDistMatrix_d AB, ElInt heightAT );
ElError ElLockedPartitionDownDist_c
( ElConstDistMatrix_c A, ElDistMatrix_c AT, ElDistMatrix_c AB, ElInt heightAT );
ElError ElLockedPartitionDownDist_z
( ElConstDistMatrix_z A, ElDistMatrix_z AT, ElDistMatrix_z AB, ElInt heightAT );

/* Partition upwards from the bottom
   ================================= */
ElError ElPartitionUp_i
( ElMatrix_i A, ElMatrix_i AT, ElMatrix_i AB, ElInt heightAB );
ElError ElPartitionUp_s
( ElMatrix_s A, ElMatrix_s AT, ElMatrix_s AB, ElInt heightAB );
ElError ElPartitionUp_d
( ElMatrix_d A, ElMatrix_d AT, ElMatrix_d AB, ElInt heightAB );
ElError ElPartitionUp_c
( ElMatrix_c A, ElMatrix_c AT, ElMatrix_c AB, ElInt heightAB );
ElError ElPartitionUp_z
( ElMatrix_z A, ElMatrix_z AT, ElMatrix_z AB, ElInt heightAB );
ElError ElPartitionUpDist_i
( ElDistMatrix_i A, ElDistMatrix_i AT, ElDistMatrix_i AB, ElInt heightAB );
ElError ElPartitionUpDist_s
( ElDistMatrix_s A, ElDistMatrix_s AT, ElDistMatrix_s AB, ElInt heightAB );
ElError ElPartitionUpDist_d
( ElDistMatrix_d A, ElDistMatrix_d AT, ElDistMatrix_d AB, ElInt heightAB );
ElError ElPartitionUpDist_c
( ElDistMatrix_c A, ElDistMatrix_c AT, ElDistMatrix_c AB, ElInt heightAB );
ElError ElPartitionUpDist_z
( ElDistMatrix_z A, ElDistMatrix_z AT, ElDistMatrix_z AB, ElInt heightAB );

ElError ElLockedPartitionUp_i
( ElConstMatrix_i A, ElMatrix_i AT, ElMatrix_i AB, ElInt heightAB );
ElError ElLockedPartitionUp_s
( ElConstMatrix_s A, ElMatrix_s AT, ElMatrix_s AB, ElInt heightAB );
ElError ElLockedPartitionUp_d
( ElConstMatrix_d A, ElMatrix_d AT, ElMatrix_d AB, ElInt heightAB );
ElError ElLockedPartitionUp_c
( ElConstMatrix_c A, ElMatrix_c AT, ElMatrix_c AB, ElInt heightAB );
ElError ElLockedPartitionUp_z
( ElConstMatrix_z A, ElMatrix_z AT, ElMatrix_z AB, ElInt heightAB );
ElError ElLockedPartitionUpDist_i
( ElConstDistMatrix_i A, ElDistMatrix_i AT, ElDistMatrix_i AB, ElInt heightAB );
ElError ElLockedPartitionUpDist_s
( ElConstDistMatrix_s A, ElDistMatrix_s AT, ElDistMatrix_s AB, ElInt heightAB );
ElError ElLockedPartitionUpDist_d
( ElConstDistMatrix_d A, ElDistMatrix_d AT, ElDistMatrix_d AB, ElInt heightAB );
ElError ElLockedPartitionUpDist_c
( ElConstDistMatrix_c A, ElDistMatrix_c AT, ElDistMatrix_c AB, ElInt heightAB );
ElError ElLockedPartitionUpDist_z
( ElConstDistMatrix_z A, ElDistMatrix_z AT, ElDistMatrix_z AB, ElInt heightAB );

/* Partition leftward from the right
   ================================= */
ElError ElPartitionLeft_i
( ElMatrix_i A, ElMatrix_i AL, ElMatrix_i AR, ElInt widthAR );
ElError ElPartitionLeft_s
( ElMatrix_s A, ElMatrix_s AL, ElMatrix_s AR, ElInt widthAR );
ElError ElPartitionLeft_d
( ElMatrix_d A, ElMatrix_d AL, ElMatrix_d AR, ElInt widthAR );
ElError ElPartitionLeft_c
( ElMatrix_c A, ElMatrix_c AL, ElMatrix_c AR, ElInt widthAR );
ElError ElPartitionLeft_z
( ElMatrix_z A, ElMatrix_z AL, ElMatrix_z AR, ElInt widthAR );
ElError ElPartitionLeftDist_i
( ElDistMatrix_i A, ElDistMatrix_i AL, ElDistMatrix_i AR, ElInt widthAR );
ElError ElPartitionLeftDist_s
( ElDistMatrix_s A, ElDistMatrix_s AL, ElDistMatrix_s AR, ElInt widthAR );
ElError ElPartitionLeftDist_d
( ElDistMatrix_d A, ElDistMatrix_d AL, ElDistMatrix_d AR, ElInt widthAR );
ElError ElPartitionLeftDist_c
( ElDistMatrix_c A, ElDistMatrix_c AL, ElDistMatrix_c AR, ElInt widthAR );
ElError ElPartitionLeftDist_z
( ElDistMatrix_z A, ElDistMatrix_z AL, ElDistMatrix_z AR, ElInt widthAR );

ElError ElLockedPartitionLeft_i
( ElConstMatrix_i A, ElMatrix_i AL, ElMatrix_i AR, ElInt widthAR );
ElError ElLockedPartitionLeft_s
( ElConstMatrix_s A, ElMatrix_s AL, ElMatrix_s AR, ElInt widthAR );
ElError ElLockedPartitionLeft_d
( ElConstMatrix_d A, ElMatrix_d AL, ElMatrix_d AR, ElInt widthAR );
ElError ElLockedPartitionLeft_c
( ElConstMatrix_c A, ElMatrix_c AL, ElMatrix_c AR, ElInt widthAR );
ElError ElLockedPartitionLeft_z
( ElConstMatrix_z A, ElMatrix_z AL, ElMatrix_z AR, ElInt widthAR );
ElError ElLockedPartitionLeftDist_i
( ElConstDistMatrix_i A, ElDistMatrix_i AL, ElDistMatrix_i AR, ElInt widthAR );
ElError ElLockedPartitionLeftDist_s
( ElConstDistMatrix_s A, ElDistMatrix_s AL, ElDistMatrix_s AR, ElInt widthAR );
ElError ElLockedPartitionLeftDist_d
( ElConstDistMatrix_d A, ElDistMatrix_d AL, ElDistMatrix_d AR, ElInt widthAR );
ElError ElLockedPartitionLeftDist_c
( ElConstDistMatrix_c A, ElDistMatrix_c AL, ElDistMatrix_c AR, ElInt widthAR );
ElError ElLockedPartitionLeftDist_z
( ElConstDistMatrix_z A, ElDistMatrix_z AL, ElDistMatrix_z AR, ElInt widthAR );

/* Partition rightward from the left
   ================================= */
ElError ElPartitionRight_i
( ElMatrix_i A, ElMatrix_i AL, ElMatrix_i AR, ElInt widthAL );
ElError ElPartitionRight_s
( ElMatrix_s A, ElMatrix_s AL, ElMatrix_s AR, ElInt widthAL );
ElError ElPartitionRight_d
( ElMatrix_d A, ElMatrix_d AL, ElMatrix_d AR, ElInt widthAL );
ElError ElPartitionRight_c
( ElMatrix_c A, ElMatrix_c AL, ElMatrix_c AR, ElInt widthAL );
ElError ElPartitionRight_z
( ElMatrix_z A, ElMatrix_z AL, ElMatrix_z AR, ElInt widthAL );
ElError ElPartitionRightDist_i
( ElDistMatrix_i A, ElDistMatrix_i AL, ElDistMatrix_i AR, ElInt widthAL );
ElError ElPartitionRightDist_s
( ElDistMatrix_s A, ElDistMatrix_s AL, ElDistMatrix_s AR, ElInt widthAL );
ElError ElPartitionRightDist_d
( ElDistMatrix_d A, ElDistMatrix_d AL, ElDistMatrix_d AR, ElInt widthAL );
ElError ElPartitionRightDist_c
( ElDistMatrix_c A, ElDistMatrix_c AL, ElDistMatrix_c AR, ElInt widthAL );
ElError ElPartitionRightDist_z
( ElDistMatrix_z A, ElDistMatrix_z AL, ElDistMatrix_z AR, ElInt widthAL );

ElError ElLockedPartitionRight_i
( ElConstMatrix_i A, ElMatrix_i AL, ElMatrix_i AR, ElInt widthAL );
ElError ElLockedPartitionRight_s
( ElConstMatrix_s A, ElMatrix_s AL, ElMatrix_s AR, ElInt widthAL );
ElError ElLockedPartitionRight_d
( ElConstMatrix_d A, ElMatrix_d AL, ElMatrix_d AR, ElInt widthAL );
ElError ElLockedPartitionRight_c
( ElConstMatrix_c A, ElMatrix_c AL, ElMatrix_c AR, ElInt widthAL );
ElError ElLockedPartitionRight_z
( ElConstMatrix_z A, ElMatrix_z AL, ElMatrix_z AR, ElInt widthAL );
ElError ElLockedPartitionRightDist_i
( ElConstDistMatrix_i A, ElDistMatrix_i AL, ElDistMatrix_i AR, ElInt widthAL );
ElError ElLockedPartitionRightDist_s
( ElConstDistMatrix_s A, ElDistMatrix_s AL, ElDistMatrix_s AR, ElInt widthAL );
ElError ElLockedPartitionRightDist_d
( ElConstDistMatrix_d A, ElDistMatrix_d AL, ElDistMatrix_d AR, ElInt widthAL );
ElError ElLockedPartitionRightDist_c
( ElConstDistMatrix_c A, ElDistMatrix_c AL, ElDistMatrix_c AR, ElInt widthAL );
ElError ElLockedPartitionRightDist_z
( ElConstDistMatrix_z A, ElDistMatrix_z AL, ElDistMatrix_z AR, ElInt widthAL );

/* Down the main diagonal
   ====================== */
ElError ElPartitionDownDiagonal_i
( ElMatrix_i A, 
  ElMatrix_i ATL, ElMatrix_i ATR, 
  ElMatrix_i ABL, ElMatrix_i ABR, ElInt diagDist );
ElError ElPartitionDownDiagonal_s
( ElMatrix_s A, 
  ElMatrix_s ATL, ElMatrix_s ATR, 
  ElMatrix_s ABL, ElMatrix_s ABR, ElInt diagDist );
ElError ElPartitionDownDiagonal_d
( ElMatrix_d A, 
  ElMatrix_d ATL, ElMatrix_d ATR, 
  ElMatrix_d ABL, ElMatrix_d ABR, ElInt diagDist );
ElError ElPartitionDownDiagonal_c
( ElMatrix_c A, 
  ElMatrix_c ATL, ElMatrix_c ATR, 
  ElMatrix_c ABL, ElMatrix_c ABR, ElInt diagDist );
ElError ElPartitionDownDiagonal_z
( ElMatrix_z A, 
  ElMatrix_z ATL, ElMatrix_z ATR, 
  ElMatrix_z ABL, ElMatrix_z ABR, ElInt diagDist );
ElError ElPartitionDownDiagonalDist_i
( ElDistMatrix_i A, 
  ElDistMatrix_i ATL, ElDistMatrix_i ATR, 
  ElDistMatrix_i ABL, ElDistMatrix_i ABR, ElInt diagDist );
ElError ElPartitionDownDiagonalDist_s
( ElDistMatrix_s A, 
  ElDistMatrix_s ATL, ElDistMatrix_s ATR, 
  ElDistMatrix_s ABL, ElDistMatrix_s ABR, ElInt diagDist );
ElError ElPartitionDownDiagonalDist_d
( ElDistMatrix_d A, 
  ElDistMatrix_d ATL, ElDistMatrix_d ATR, 
  ElDistMatrix_d ABL, ElDistMatrix_d ABR, ElInt diagDist );
ElError ElPartitionDownDiagonalDist_c
( ElDistMatrix_c A, 
  ElDistMatrix_c ATL, ElDistMatrix_c ATR, 
  ElDistMatrix_c ABL, ElDistMatrix_c ABR, ElInt diagDist );
ElError ElPartitionDownDiagonalDist_z
( ElDistMatrix_z A, 
  ElDistMatrix_z ATL, ElDistMatrix_z ATR, 
  ElDistMatrix_z ABL, ElDistMatrix_z ABR, ElInt diagDist );

ElError ElLockedPartitionDownDiagonal_i
( ElConstMatrix_i A, 
  ElMatrix_i ATL, ElMatrix_i ATR, 
  ElMatrix_i ABL, ElMatrix_i ABR, ElInt diagDist );
ElError ElLockedPartitionDownDiagonal_s
( ElConstMatrix_s A, 
  ElMatrix_s ATL, ElMatrix_s ATR, 
  ElMatrix_s ABL, ElMatrix_s ABR, ElInt diagDist );
ElError ElLockedPartitionDownDiagonal_d
( ElConstMatrix_d A, 
  ElMatrix_d ATL, ElMatrix_d ATR, 
  ElMatrix_d ABL, ElMatrix_d ABR, ElInt diagDist );
ElError ElLockedPartitionDownDiagonal_c
( ElConstMatrix_c A, 
  ElMatrix_c ATL, ElMatrix_c ATR, 
  ElMatrix_c ABL, ElMatrix_c ABR, ElInt diagDist );
ElError ElLockedPartitionDownDiagonal_z
( ElConstMatrix_z A, 
  ElMatrix_z ATL, ElMatrix_z ATR, 
  ElMatrix_z ABL, ElMatrix_z ABR, ElInt diagDist );
ElError ElLockedPartitionDownDiagonalDist_i
( ElConstDistMatrix_i A, 
  ElDistMatrix_i ATL, ElDistMatrix_i ATR, 
  ElDistMatrix_i ABL, ElDistMatrix_i ABR, ElInt diagDist );
ElError ElLockedPartitionDownDiagonalDist_s
( ElConstDistMatrix_s A, 
  ElDistMatrix_s ATL, ElDistMatrix_s ATR, 
  ElDistMatrix_s ABL, ElDistMatrix_s ABR, ElInt diagDist );
ElError ElLockedPartitionDownDiagonalDist_d
( ElConstDistMatrix_d A, 
  ElDistMatrix_d ATL, ElDistMatrix_d ATR, 
  ElDistMatrix_d ABL, ElDistMatrix_d ABR, ElInt diagDist );
ElError ElLockedPartitionDownDiagonalDist_c
( ElConstDistMatrix_c A, 
  ElDistMatrix_c ATL, ElDistMatrix_c ATR, 
  ElDistMatrix_c ABL, ElDistMatrix_c ABR, ElInt diagDist );
ElError ElLockedPartitionDownDiagonalDist_z
( ElConstDistMatrix_z A, 
  ElDistMatrix_z ATL, ElDistMatrix_z ATR, 
  ElDistMatrix_z ABL, ElDistMatrix_z ABR, ElInt diagDist );

/* Down an offset diagonal
   ======================= */
ElError ElPartitionDownOffsetDiagonal_i
( ElInt offset, ElMatrix_i A, 
  ElMatrix_i ATL, ElMatrix_i ATR, 
  ElMatrix_i ABL, ElMatrix_i ABR, ElInt diagDist );
ElError ElPartitionDownOffsetDiagonal_s
( ElInt offset, ElMatrix_s A, 
  ElMatrix_s ATL, ElMatrix_s ATR, 
  ElMatrix_s ABL, ElMatrix_s ABR, ElInt diagDist );
ElError ElPartitionDownOffsetDiagonal_d
( ElInt offset, ElMatrix_d A, 
  ElMatrix_d ATL, ElMatrix_d ATR, 
  ElMatrix_d ABL, ElMatrix_d ABR, ElInt diagDist );
ElError ElPartitionDownOffsetDiagonal_c
( ElInt offset, ElMatrix_c A, 
  ElMatrix_c ATL, ElMatrix_c ATR, 
  ElMatrix_c ABL, ElMatrix_c ABR, ElInt diagDist );
ElError ElPartitionDownOffsetDiagonal_z
( ElInt offset, ElMatrix_z A, 
  ElMatrix_z ATL, ElMatrix_z ATR, 
  ElMatrix_z ABL, ElMatrix_z ABR, ElInt diagDist );
ElError ElPartitionDownOffsetDiagonalDist_i
( ElInt offset, ElDistMatrix_i A, 
  ElDistMatrix_i ATL, ElDistMatrix_i ATR, 
  ElDistMatrix_i ABL, ElDistMatrix_i ABR, ElInt diagDist );
ElError ElPartitionDownOffsetDiagonalDist_s
( ElInt offset, ElDistMatrix_s A, 
  ElDistMatrix_s ATL, ElDistMatrix_s ATR, 
  ElDistMatrix_s ABL, ElDistMatrix_s ABR, ElInt diagDist );
ElError ElPartitionDownOffsetDiagonalDist_d
( ElInt offset, ElDistMatrix_d A, 
  ElDistMatrix_d ATL, ElDistMatrix_d ATR, 
  ElDistMatrix_d ABL, ElDistMatrix_d ABR, ElInt diagDist );
ElError ElPartitionDownOffsetDiagonalDist_c
( ElInt offset, ElDistMatrix_c A, 
  ElDistMatrix_c ATL, ElDistMatrix_c ATR, 
  ElDistMatrix_c ABL, ElDistMatrix_c ABR, ElInt diagDist );
ElError ElPartitionDownOffsetDiagonalDist_z
( ElInt offset, ElDistMatrix_z A, 
  ElDistMatrix_z ATL, ElDistMatrix_z ATR, 
  ElDistMatrix_z ABL, ElDistMatrix_z ABR, ElInt diagDist );

ElError ElLockedPartitionDownOffsetDiagonal_i
( ElInt offset, ElConstMatrix_i A, 
  ElMatrix_i ATL, ElMatrix_i ATR, 
  ElMatrix_i ABL, ElMatrix_i ABR, ElInt diagDist );
ElError ElLockedPartitionDownOffsetDiagonal_s
( ElInt offset, ElConstMatrix_s A, 
  ElMatrix_s ATL, ElMatrix_s ATR, 
  ElMatrix_s ABL, ElMatrix_s ABR, ElInt diagDist );
ElError ElLockedPartitionDownOffsetDiagonal_d
( ElInt offset, ElConstMatrix_d A, 
  ElMatrix_d ATL, ElMatrix_d ATR, 
  ElMatrix_d ABL, ElMatrix_d ABR, ElInt diagDist );
ElError ElLockedPartitionDownOffsetDiagonal_c
( ElInt offset, ElConstMatrix_c A, 
  ElMatrix_c ATL, ElMatrix_c ATR, 
  ElMatrix_c ABL, ElMatrix_c ABR, ElInt diagDist );
ElError ElLockedPartitionDownOffsetDiagonal_z
( ElInt offset, ElConstMatrix_z A, 
  ElMatrix_z ATL, ElMatrix_z ATR, 
  ElMatrix_z ABL, ElMatrix_z ABR, ElInt diagDist );
ElError ElLockedPartitionDownOffsetDiagonalDist_i
( ElInt offset, ElConstDistMatrix_i A, 
  ElDistMatrix_i ATL, ElDistMatrix_i ATR, 
  ElDistMatrix_i ABL, ElDistMatrix_i ABR, ElInt diagDist );
ElError ElLockedPartitionDownOffsetDiagonalDist_s
( ElInt offset, ElConstDistMatrix_s A, 
  ElDistMatrix_s ATL, ElDistMatrix_s ATR, 
  ElDistMatrix_s ABL, ElDistMatrix_s ABR, ElInt diagDist );
ElError ElLockedPartitionDownOffsetDiagonalDist_d
( ElInt offset, ElConstDistMatrix_d A, 
  ElDistMatrix_d ATL, ElDistMatrix_d ATR, 
  ElDistMatrix_d ABL, ElDistMatrix_d ABR, ElInt diagDist );
ElError ElLockedPartitionDownOffsetDiagonalDist_c
( ElInt offset, ElConstDistMatrix_c A, 
  ElDistMatrix_c ATL, ElDistMatrix_c ATR, 
  ElDistMatrix_c ABL, ElDistMatrix_c ABR, ElInt diagDist );
ElError ElLockedPartitionDownOffsetDiagonalDist_z
( ElInt offset, ElConstDistMatrix_z A, 
  ElDistMatrix_z ATL, ElDistMatrix_z ATR, 
  ElDistMatrix_z ABL, ElDistMatrix_z ABR, ElInt diagDist );

/* Up the main diagonal
   ==================== */
ElError ElPartitionUpDiagonal_i
( ElMatrix_i A, 
  ElMatrix_i ATL, ElMatrix_i ATR, 
  ElMatrix_i ABL, ElMatrix_i ABR, ElInt diagDist );
ElError ElPartitionUpDiagonal_s
( ElMatrix_s A, 
  ElMatrix_s ATL, ElMatrix_s ATR, 
  ElMatrix_s ABL, ElMatrix_s ABR, ElInt diagDist );
ElError ElPartitionUpDiagonal_d
( ElMatrix_d A, 
  ElMatrix_d ATL, ElMatrix_d ATR, 
  ElMatrix_d ABL, ElMatrix_d ABR, ElInt diagDist );
ElError ElPartitionUpDiagonal_c
( ElMatrix_c A, 
  ElMatrix_c ATL, ElMatrix_c ATR, 
  ElMatrix_c ABL, ElMatrix_c ABR, ElInt diagDist );
ElError ElPartitionUpDiagonal_z
( ElMatrix_z A, 
  ElMatrix_z ATL, ElMatrix_z ATR, 
  ElMatrix_z ABL, ElMatrix_z ABR, ElInt diagDist );
ElError ElPartitionUpDiagonalDist_i
( ElDistMatrix_i A, 
  ElDistMatrix_i ATL, ElDistMatrix_i ATR, 
  ElDistMatrix_i ABL, ElDistMatrix_i ABR, ElInt diagDist );
ElError ElPartitionUpDiagonalDist_s
( ElDistMatrix_s A, 
  ElDistMatrix_s ATL, ElDistMatrix_s ATR, 
  ElDistMatrix_s ABL, ElDistMatrix_s ABR, ElInt diagDist );
ElError ElPartitionUpDiagonalDist_d
( ElDistMatrix_d A, 
  ElDistMatrix_d ATL, ElDistMatrix_d ATR, 
  ElDistMatrix_d ABL, ElDistMatrix_d ABR, ElInt diagDist );
ElError ElPartitionUpDiagonalDist_c
( ElDistMatrix_c A, 
  ElDistMatrix_c ATL, ElDistMatrix_c ATR, 
  ElDistMatrix_c ABL, ElDistMatrix_c ABR, ElInt diagDist );
ElError ElPartitionUpDiagonalDist_z
( ElDistMatrix_z A, 
  ElDistMatrix_z ATL, ElDistMatrix_z ATR, 
  ElDistMatrix_z ABL, ElDistMatrix_z ABR, ElInt diagDist );

ElError ElLockedPartitionUpDiagonal_i
( ElConstMatrix_i A, 
  ElMatrix_i ATL, ElMatrix_i ATR, 
  ElMatrix_i ABL, ElMatrix_i ABR, ElInt diagDist );
ElError ElLockedPartitionUpDiagonal_s
( ElConstMatrix_s A, 
  ElMatrix_s ATL, ElMatrix_s ATR, 
  ElMatrix_s ABL, ElMatrix_s ABR, ElInt diagDist );
ElError ElLockedPartitionUpDiagonal_d
( ElConstMatrix_d A, 
  ElMatrix_d ATL, ElMatrix_d ATR, 
  ElMatrix_d ABL, ElMatrix_d ABR, ElInt diagDist );
ElError ElLockedPartitionUpDiagonal_c
( ElConstMatrix_c A, 
  ElMatrix_c ATL, ElMatrix_c ATR, 
  ElMatrix_c ABL, ElMatrix_c ABR, ElInt diagDist );
ElError ElLockedPartitionUpDiagonal_z
( ElConstMatrix_z A, 
  ElMatrix_z ATL, ElMatrix_z ATR, 
  ElMatrix_z ABL, ElMatrix_z ABR, ElInt diagDist );
ElError ElLockedPartitionUpDiagonalDist_i
( ElConstDistMatrix_i A, 
  ElDistMatrix_i ATL, ElDistMatrix_i ATR, 
  ElDistMatrix_i ABL, ElDistMatrix_i ABR, ElInt diagDist );
ElError ElLockedPartitionUpDiagonalDist_s
( ElConstDistMatrix_s A, 
  ElDistMatrix_s ATL, ElDistMatrix_s ATR, 
  ElDistMatrix_s ABL, ElDistMatrix_s ABR, ElInt diagDist );
ElError ElLockedPartitionUpDiagonalDist_d
( ElConstDistMatrix_d A, 
  ElDistMatrix_d ATL, ElDistMatrix_d ATR, 
  ElDistMatrix_d ABL, ElDistMatrix_d ABR, ElInt diagDist );
ElError ElLockedPartitionUpDiagonalDist_c
( ElConstDistMatrix_c A, 
  ElDistMatrix_c ATL, ElDistMatrix_c ATR, 
  ElDistMatrix_c ABL, ElDistMatrix_c ABR, ElInt diagDist );
ElError ElLockedPartitionUpDiagonalDist_z
( ElConstDistMatrix_z A, 
  ElDistMatrix_z ATL, ElDistMatrix_z ATR, 
  ElDistMatrix_z ABL, ElDistMatrix_z ABR, ElInt diagDist );

/* Up an offset diagonal
   ===================== */
ElError ElPartitionUpOffsetDiagonal_i
( ElInt offset, ElMatrix_i A, 
  ElMatrix_i ATL, ElMatrix_i ATR, 
  ElMatrix_i ABL, ElMatrix_i ABR, ElInt diagDist );
ElError ElPartitionUpOffsetDiagonal_s
( ElInt offset, ElMatrix_s A, 
  ElMatrix_s ATL, ElMatrix_s ATR, 
  ElMatrix_s ABL, ElMatrix_s ABR, ElInt diagDist );
ElError ElPartitionUpOffsetDiagonal_d
( ElInt offset, ElMatrix_d A, 
  ElMatrix_d ATL, ElMatrix_d ATR, 
  ElMatrix_d ABL, ElMatrix_d ABR, ElInt diagDist );
ElError ElPartitionUpOffsetDiagonal_c
( ElInt offset, ElMatrix_c A, 
  ElMatrix_c ATL, ElMatrix_c ATR, 
  ElMatrix_c ABL, ElMatrix_c ABR, ElInt diagDist );
ElError ElPartitionUpOffsetDiagonal_z
( ElInt offset, ElMatrix_z A, 
  ElMatrix_z ATL, ElMatrix_z ATR, 
  ElMatrix_z ABL, ElMatrix_z ABR, ElInt diagDist );
ElError ElPartitionUpOffsetDiagonalDist_i
( ElInt offset, ElDistMatrix_i A, 
  ElDistMatrix_i ATL, ElDistMatrix_i ATR, 
  ElDistMatrix_i ABL, ElDistMatrix_i ABR, ElInt diagDist );
ElError ElPartitionUpOffsetDiagonalDist_s
( ElInt offset, ElDistMatrix_s A, 
  ElDistMatrix_s ATL, ElDistMatrix_s ATR, 
  ElDistMatrix_s ABL, ElDistMatrix_s ABR, ElInt diagDist );
ElError ElPartitionUpOffsetDiagonalDist_d
( ElInt offset, ElDistMatrix_d A, 
  ElDistMatrix_d ATL, ElDistMatrix_d ATR, 
  ElDistMatrix_d ABL, ElDistMatrix_d ABR, ElInt diagDist );
ElError ElPartitionUpOffsetDiagonalDist_c
( ElInt offset, ElDistMatrix_c A, 
  ElDistMatrix_c ATL, ElDistMatrix_c ATR, 
  ElDistMatrix_c ABL, ElDistMatrix_c ABR, ElInt diagDist );
ElError ElPartitionUpOffsetDiagonalDist_z
( ElInt offset, ElDistMatrix_z A, 
  ElDistMatrix_z ATL, ElDistMatrix_z ATR, 
  ElDistMatrix_z ABL, ElDistMatrix_z ABR, ElInt diagDist );

ElError ElLockedPartitionUpOffsetDiagonal_i
( ElInt offset, ElConstMatrix_i A, 
  ElMatrix_i ATL, ElMatrix_i ATR, 
  ElMatrix_i ABL, ElMatrix_i ABR, ElInt diagDist );
ElError ElLockedPartitionUpOffsetDiagonal_s
( ElInt offset, ElConstMatrix_s A, 
  ElMatrix_s ATL, ElMatrix_s ATR, 
  ElMatrix_s ABL, ElMatrix_s ABR, ElInt diagDist );
ElError ElLockedPartitionUpOffsetDiagonal_d
( ElInt offset, ElConstMatrix_d A, 
  ElMatrix_d ATL, ElMatrix_d ATR, 
  ElMatrix_d ABL, ElMatrix_d ABR, ElInt diagDist );
ElError ElLockedPartitionUpOffsetDiagonal_c
( ElInt offset, ElConstMatrix_c A, 
  ElMatrix_c ATL, ElMatrix_c ATR, 
  ElMatrix_c ABL, ElMatrix_c ABR, ElInt diagDist );
ElError ElLockedPartitionUpOffsetDiagonal_z
( ElInt offset, ElConstMatrix_z A, 
  ElMatrix_z ATL, ElMatrix_z ATR, 
  ElMatrix_z ABL, ElMatrix_z ABR, ElInt diagDist );
ElError ElLockedPartitionUpOffsetDiagonalDist_i
( ElInt offset, ElConstDistMatrix_i A, 
  ElDistMatrix_i ATL, ElDistMatrix_i ATR, 
  ElDistMatrix_i ABL, ElDistMatrix_i ABR, ElInt diagDist );
ElError ElLockedPartitionUpOffsetDiagonalDist_s
( ElInt offset, ElConstDistMatrix_s A, 
  ElDistMatrix_s ATL, ElDistMatrix_s ATR, 
  ElDistMatrix_s ABL, ElDistMatrix_s ABR, ElInt diagDist );
ElError ElLockedPartitionUpOffsetDiagonalDist_d
( ElInt offset, ElConstDistMatrix_d A, 
  ElDistMatrix_d ATL, ElDistMatrix_d ATR, 
  ElDistMatrix_d ABL, ElDistMatrix_d ABR, ElInt diagDist );
ElError ElLockedPartitionUpOffsetDiagonalDist_c
( ElInt offset, ElConstDistMatrix_c A, 
  ElDistMatrix_c ATL, ElDistMatrix_c ATR, 
  ElDistMatrix_c ABL, ElDistMatrix_c ABR, ElInt diagDist );
ElError ElLockedPartitionUpOffsetDiagonalDist_z
( ElInt offset, ElConstDistMatrix_z A, 
  ElDistMatrix_z ATL, ElDistMatrix_z ATR, 
  ElDistMatrix_z ABL, ElDistMatrix_z ABR, ElInt diagDist );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_FLAMEPART_PARTITION_C_H */
