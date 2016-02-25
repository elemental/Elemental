/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_VIEW_C_H
#define EL_VIEW_C_H

#ifdef __cplusplus
extern "C" {
#endif

/* UNEXPLAINED COMPILER BUG: 
  'cc' 4.8.2 in Ubuntu fails if "ElRange_i iRange" is replaced with 
  "ElRange_i I" */

EL_EXPORT ElError ElView_i
( ElMatrix_i A, ElMatrix_i B, ElRange_i iRange, ElRange_i jRange );
EL_EXPORT ElError ElView_s
( ElMatrix_s A, ElMatrix_s B, ElRange_i iRange, ElRange_i jRange );
EL_EXPORT ElError ElView_d
( ElMatrix_d A, ElMatrix_d B, ElRange_i iRange, ElRange_i jRange );
EL_EXPORT ElError ElView_c
( ElMatrix_c A, ElMatrix_c B, ElRange_i iRange, ElRange_i jRange );
EL_EXPORT ElError ElView_z
( ElMatrix_z A, ElMatrix_z B, ElRange_i iRange, ElRange_i jRange );
EL_EXPORT ElError ElViewDist_i
( ElDistMatrix_i A, ElDistMatrix_i B, ElRange_i iRange, ElRange_i jRange );
EL_EXPORT ElError ElViewDist_s
( ElDistMatrix_s A, ElDistMatrix_s B, ElRange_i iRange, ElRange_i jRange );
EL_EXPORT ElError ElViewDist_d
( ElDistMatrix_d A, ElDistMatrix_d B, ElRange_i iRange, ElRange_i jRange );
EL_EXPORT ElError ElViewDist_c
( ElDistMatrix_c A, ElDistMatrix_c B, ElRange_i iRange, ElRange_i jRange );
EL_EXPORT ElError ElViewDist_z
( ElDistMatrix_z A, ElDistMatrix_z B, ElRange_i iRange, ElRange_i jRange );

EL_EXPORT ElError ElLockedView_i
( ElMatrix_i A, ElConstMatrix_i B, ElRange_i iRange, ElRange_i jRange );
EL_EXPORT ElError ElLockedView_s
( ElMatrix_s A, ElConstMatrix_s B, ElRange_i iRange, ElRange_i jRange );
EL_EXPORT ElError ElLockedView_d
( ElMatrix_d A, ElConstMatrix_d B, ElRange_i iRange, ElRange_i jRange );
EL_EXPORT ElError ElLockedView_c
( ElMatrix_c A, ElConstMatrix_c B, ElRange_i iRange, ElRange_i jRange );
EL_EXPORT ElError ElLockedView_z
( ElMatrix_z A, ElConstMatrix_z B, ElRange_i iRange, ElRange_i jRange );
EL_EXPORT ElError ElLockedViewDist_i
( ElDistMatrix_i A, ElConstDistMatrix_i B, ElRange_i iRange, ElRange_i jRange );
EL_EXPORT ElError ElLockedViewDist_s
( ElDistMatrix_s A, ElConstDistMatrix_s B, ElRange_i iRange, ElRange_i jRange );
EL_EXPORT ElError ElLockedViewDist_d
( ElDistMatrix_d A, ElConstDistMatrix_d B, ElRange_i iRange, ElRange_i jRange );
EL_EXPORT ElError ElLockedViewDist_c
( ElDistMatrix_c A, ElConstDistMatrix_c B, ElRange_i iRange, ElRange_i jRange );
EL_EXPORT ElError ElLockedViewDist_z
( ElDistMatrix_z A, ElConstDistMatrix_z B, ElRange_i iRange, ElRange_i jRange );

EL_EXPORT ElError ElViewOffset_i
( ElMatrix_i A, ElMatrix_i B, 
  ElInt i, ElInt j, ElInt height, ElInt width );
EL_EXPORT ElError ElViewOffset_s
( ElMatrix_s A, ElMatrix_s B, 
  ElInt i, ElInt j, ElInt height, ElInt width );
EL_EXPORT ElError ElViewOffset_d
( ElMatrix_d A, ElMatrix_d B, 
  ElInt i, ElInt j, ElInt height, ElInt width );
EL_EXPORT ElError ElViewOffset_c
( ElMatrix_c A, ElMatrix_c B, 
  ElInt i, ElInt j, ElInt height, ElInt width );
EL_EXPORT ElError ElViewOffset_z
( ElMatrix_z A, ElMatrix_z B, 
  ElInt i, ElInt j, ElInt height, ElInt width );
EL_EXPORT ElError ElViewOffsetDist_i
( ElDistMatrix_i A, ElDistMatrix_i B, 
  ElInt i, ElInt j, ElInt height, ElInt width );
EL_EXPORT ElError ElViewOffsetDist_s
( ElDistMatrix_s A, ElDistMatrix_s B, 
  ElInt i, ElInt j, ElInt height, ElInt width );
EL_EXPORT ElError ElViewOffsetDist_d
( ElDistMatrix_d A, ElDistMatrix_d B, 
  ElInt i, ElInt j, ElInt height, ElInt width );
EL_EXPORT ElError ElViewOffsetDist_c
( ElDistMatrix_c A, ElDistMatrix_c B, 
  ElInt i, ElInt j, ElInt height, ElInt width );
EL_EXPORT ElError ElViewOffsetDist_z
( ElDistMatrix_z A, ElDistMatrix_z B, 
  ElInt i, ElInt j, ElInt height, ElInt width );

EL_EXPORT ElError ElLockedViewOffset_i
( ElMatrix_i A, ElConstMatrix_i B, 
  ElInt i, ElInt j, ElInt height, ElInt width );
EL_EXPORT ElError ElLockedViewOffset_s
( ElMatrix_s A, ElConstMatrix_s B, 
  ElInt i, ElInt j, ElInt height, ElInt width );
EL_EXPORT ElError ElLockedViewOffset_d
( ElMatrix_d A, ElConstMatrix_d B, 
  ElInt i, ElInt j, ElInt height, ElInt width );
EL_EXPORT ElError ElLockedViewOffset_c
( ElMatrix_c A, ElConstMatrix_c B, 
  ElInt i, ElInt j, ElInt height, ElInt width );
EL_EXPORT ElError ElLockedViewOffset_z
( ElMatrix_z A, ElConstMatrix_z B, 
  ElInt i, ElInt j, ElInt height, ElInt width );
EL_EXPORT ElError ElLockedViewOffsetDist_i
( ElDistMatrix_i A, ElConstDistMatrix_i B, 
  ElInt i, ElInt j, ElInt height, ElInt width );
EL_EXPORT ElError ElLockedViewOffsetDist_s
( ElDistMatrix_s A, ElConstDistMatrix_s B, 
  ElInt i, ElInt j, ElInt height, ElInt width );
EL_EXPORT ElError ElLockedViewOffsetDist_d
( ElDistMatrix_d A, ElConstDistMatrix_d B, 
  ElInt i, ElInt j, ElInt height, ElInt width );
EL_EXPORT ElError ElLockedViewOffsetDist_c
( ElDistMatrix_c A, ElConstDistMatrix_c B, 
  ElInt i, ElInt j, ElInt height, ElInt width );
EL_EXPORT ElError ElLockedViewOffsetDist_z
( ElDistMatrix_z A, ElConstDistMatrix_z B, 
  ElInt i, ElInt j, ElInt height, ElInt width );

EL_EXPORT ElError ElViewFull_i( ElMatrix_i A, ElMatrix_i B );
EL_EXPORT ElError ElViewFull_s( ElMatrix_s A, ElMatrix_s B );
EL_EXPORT ElError ElViewFull_d( ElMatrix_d A, ElMatrix_d B );
EL_EXPORT ElError ElViewFull_c( ElMatrix_c A, ElMatrix_c B );
EL_EXPORT ElError ElViewFull_z( ElMatrix_z A, ElMatrix_z B );
EL_EXPORT ElError ElViewFullDist_i( ElDistMatrix_i A, ElDistMatrix_i B );
EL_EXPORT ElError ElViewFullDist_s( ElDistMatrix_s A, ElDistMatrix_s B );
EL_EXPORT ElError ElViewFullDist_d( ElDistMatrix_d A, ElDistMatrix_d B );
EL_EXPORT ElError ElViewFullDist_c( ElDistMatrix_c A, ElDistMatrix_c B );
EL_EXPORT ElError ElViewFullDist_z( ElDistMatrix_z A, ElDistMatrix_z B );

EL_EXPORT ElError ElLockedViewFull_i( ElMatrix_i A, ElConstMatrix_i B );
EL_EXPORT ElError ElLockedViewFull_s( ElMatrix_s A, ElConstMatrix_s B );
EL_EXPORT ElError ElLockedViewFull_d( ElMatrix_d A, ElConstMatrix_d B );
EL_EXPORT ElError ElLockedViewFull_c( ElMatrix_c A, ElConstMatrix_c B );
EL_EXPORT ElError ElLockedViewFull_z( ElMatrix_z A, ElConstMatrix_z B );
EL_EXPORT ElError ElLockedViewFullDist_i
( ElDistMatrix_i A, ElConstDistMatrix_i B );
EL_EXPORT ElError ElLockedViewFullDist_s
( ElDistMatrix_s A, ElConstDistMatrix_s B );
EL_EXPORT ElError ElLockedViewFullDist_d
( ElDistMatrix_d A, ElConstDistMatrix_d B );
EL_EXPORT ElError ElLockedViewFullDist_c
( ElDistMatrix_c A, ElConstDistMatrix_c B );
EL_EXPORT ElError ElLockedViewFullDist_z
( ElDistMatrix_z A, ElConstDistMatrix_z B );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_VIEW_C_H */
