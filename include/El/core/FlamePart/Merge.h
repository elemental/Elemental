/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_FLAMEPART_MERGE_C_H
#define EL_FLAMEPART_MERGE_C_H

#ifdef __cplusplus
extern "C" {
#endif

/* Horizontally merge two contiguous matrices
   ========================================== */

EL_EXPORT ElError ElMerge1x2_i( ElMatrix_i A, ElMatrix_i BL, ElMatrix_i BR );
EL_EXPORT ElError ElMerge1x2_s( ElMatrix_s A, ElMatrix_s BL, ElMatrix_s BR );
EL_EXPORT ElError ElMerge1x2_d( ElMatrix_d A, ElMatrix_d BL, ElMatrix_d BR );
EL_EXPORT ElError ElMerge1x2_c( ElMatrix_c A, ElMatrix_c BL, ElMatrix_c BR );
EL_EXPORT ElError ElMerge1x2_z( ElMatrix_z A, ElMatrix_z BL, ElMatrix_z BR );
EL_EXPORT ElError ElMerge1x2Dist_i
( ElDistMatrix_i A, ElDistMatrix_i BL, ElDistMatrix_i BR );
EL_EXPORT ElError ElMerge1x2Dist_s
( ElDistMatrix_s A, ElDistMatrix_s BL, ElDistMatrix_s BR );
EL_EXPORT ElError ElMerge1x2Dist_d
( ElDistMatrix_d A, ElDistMatrix_d BL, ElDistMatrix_d BR );
EL_EXPORT ElError ElMerge1x2Dist_c
( ElDistMatrix_c A, ElDistMatrix_c BL, ElDistMatrix_c BR );
EL_EXPORT ElError ElMerge1x2Dist_z
( ElDistMatrix_z A, ElDistMatrix_z BL, ElDistMatrix_z BR );

EL_EXPORT ElError ElLockedMerge1x2_i
( ElMatrix_i A, ElConstMatrix_i BL, ElConstMatrix_i BR );
EL_EXPORT ElError ElLockedMerge1x2_s
( ElMatrix_s A, ElConstMatrix_s BL, ElConstMatrix_s BR );
EL_EXPORT ElError ElLockedMerge1x2_d
( ElMatrix_d A, ElConstMatrix_d BL, ElConstMatrix_d BR );
EL_EXPORT ElError ElLockedMerge1x2_c
( ElMatrix_c A, ElConstMatrix_c BL, ElConstMatrix_c BR );
EL_EXPORT ElError ElLockedMerge1x2_z
( ElMatrix_z A, ElConstMatrix_z BL, ElConstMatrix_z BR );
EL_EXPORT ElError ElLockedMerge1x2Dist_i
( ElDistMatrix_i A, ElConstDistMatrix_i BL, ElConstDistMatrix_i BR );
EL_EXPORT ElError ElLockedMerge1x2Dist_s
( ElDistMatrix_s A, ElConstDistMatrix_s BL, ElConstDistMatrix_s BR );
EL_EXPORT ElError ElLockedMerge1x2Dist_d
( ElDistMatrix_d A, ElConstDistMatrix_d BL, ElConstDistMatrix_d BR );
EL_EXPORT ElError ElLockedMerge1x2Dist_c
( ElDistMatrix_c A, ElConstDistMatrix_c BL, ElConstDistMatrix_c BR );
EL_EXPORT ElError ElLockedMerge1x2Dist_z
( ElDistMatrix_z A, ElConstDistMatrix_z BL, ElConstDistMatrix_z BR );

/* Vertically merge two contiguous matrices
   ======================================== */

EL_EXPORT ElError ElMerge2x1_i( ElMatrix_i A, ElMatrix_i BT, ElMatrix_i BB );
EL_EXPORT ElError ElMerge2x1_s( ElMatrix_s A, ElMatrix_s BT, ElMatrix_s BB );
EL_EXPORT ElError ElMerge2x1_d( ElMatrix_d A, ElMatrix_d BT, ElMatrix_d BB );
EL_EXPORT ElError ElMerge2x1_c( ElMatrix_c A, ElMatrix_c BT, ElMatrix_c BB );
EL_EXPORT ElError ElMerge2x1_z( ElMatrix_z A, ElMatrix_z BT, ElMatrix_z BB );
EL_EXPORT ElError ElMerge2x1Dist_i
( ElDistMatrix_i A, ElDistMatrix_i BT, ElDistMatrix_i BB );
EL_EXPORT ElError ElMerge2x1Dist_s
( ElDistMatrix_s A, ElDistMatrix_s BT, ElDistMatrix_s BB );
EL_EXPORT ElError ElMerge2x1Dist_d
( ElDistMatrix_d A, ElDistMatrix_d BT, ElDistMatrix_d BB );
EL_EXPORT ElError ElMerge2x1Dist_c
( ElDistMatrix_c A, ElDistMatrix_c BT, ElDistMatrix_c BB );
EL_EXPORT ElError ElMerge2x1Dist_z
( ElDistMatrix_z A, ElDistMatrix_z BT, ElDistMatrix_z BB );

EL_EXPORT ElError ElLockedMerge2x1_i
( ElMatrix_i A, ElConstMatrix_i BT, ElConstMatrix_i BB );
EL_EXPORT ElError ElLockedMerge2x1_s
( ElMatrix_s A, ElConstMatrix_s BT, ElConstMatrix_s BB );
EL_EXPORT ElError ElLockedMerge2x1_d
( ElMatrix_d A, ElConstMatrix_d BT, ElConstMatrix_d BB );
EL_EXPORT ElError ElLockedMerge2x1_c
( ElMatrix_c A, ElConstMatrix_c BT, ElConstMatrix_c BB );
EL_EXPORT ElError ElLockedMerge2x1_z
( ElMatrix_z A, ElConstMatrix_z BT, ElConstMatrix_z BB );
EL_EXPORT ElError ElLockedMerge2x1Dist_i
( ElDistMatrix_i A, ElConstDistMatrix_i BT, ElConstDistMatrix_i BB );
EL_EXPORT ElError ElLockedMerge2x1Dist_s
( ElDistMatrix_s A, ElConstDistMatrix_s BT, ElConstDistMatrix_s BB );
EL_EXPORT ElError ElLockedMerge2x1Dist_d
( ElDistMatrix_d A, ElConstDistMatrix_d BT, ElConstDistMatrix_d BB );
EL_EXPORT ElError ElLockedMerge2x1Dist_c
( ElDistMatrix_c A, ElConstDistMatrix_c BT, ElConstDistMatrix_c BB );
EL_EXPORT ElError ElLockedMerge2x1Dist_z
( ElDistMatrix_z A, ElConstDistMatrix_z BT, ElConstDistMatrix_z BB );

/* Merge a contiguous 2x2 block of matrices
   ======================================== */

EL_EXPORT ElError ElMerge2x2_i
( ElMatrix_i A, 
  ElMatrix_i BTL, ElMatrix_i BTR, 
  ElMatrix_i BBL, ElMatrix_i BBR );
EL_EXPORT ElError ElMerge2x2_s
( ElMatrix_s A, 
  ElMatrix_s BTL, ElMatrix_s BTR, 
  ElMatrix_s BBL, ElMatrix_s BBR );
EL_EXPORT ElError ElMerge2x2_d
( ElMatrix_d A, 
  ElMatrix_d BTL, ElMatrix_d BTR, 
  ElMatrix_d BBL, ElMatrix_d BBR );
EL_EXPORT ElError ElMerge2x2_c
( ElMatrix_c A, 
  ElMatrix_c BTL, ElMatrix_c BTR, 
  ElMatrix_c BBL, ElMatrix_c BBR );
EL_EXPORT ElError ElMerge2x2_z
( ElMatrix_z A, 
  ElMatrix_z BTL, ElMatrix_z BTR, 
  ElMatrix_z BBL, ElMatrix_z BBR );
EL_EXPORT ElError ElMerge2x2Dist_i
( ElDistMatrix_i A, 
  ElDistMatrix_i BTL, ElDistMatrix_i BTR, 
  ElDistMatrix_i BBL, ElDistMatrix_i BBR );
EL_EXPORT ElError ElMerge2x2Dist_s
( ElDistMatrix_s A, 
  ElDistMatrix_s BTL, ElDistMatrix_s BTR, 
  ElDistMatrix_s BBL, ElDistMatrix_s BBR );
EL_EXPORT ElError ElMerge2x2Dist_d
( ElDistMatrix_d A, 
  ElDistMatrix_d BTL, ElDistMatrix_d BTR, 
  ElDistMatrix_d BBL, ElDistMatrix_d BBR );
EL_EXPORT ElError ElMerge2x2Dist_c
( ElDistMatrix_c A, 
  ElDistMatrix_c BTL, ElDistMatrix_c BTR, 
  ElDistMatrix_c BBL, ElDistMatrix_c BBR );
EL_EXPORT ElError ElMerge2x2Dist_z
( ElDistMatrix_z A, 
  ElDistMatrix_z BTL, ElDistMatrix_z BTR, 
  ElDistMatrix_z BBL, ElDistMatrix_z BBR );

EL_EXPORT ElError ElLockedMerge2x2_i
( ElMatrix_i A, 
  ElConstMatrix_i BTL, ElConstMatrix_i BTR, 
  ElConstMatrix_i BBL, ElConstMatrix_i BBR );
EL_EXPORT ElError ElLockedMerge2x2_s
( ElMatrix_s A, 
  ElConstMatrix_s BTL, ElConstMatrix_s BTR, 
  ElConstMatrix_s BBL, ElConstMatrix_s BBR );
EL_EXPORT ElError ElLockedMerge2x2_d
( ElMatrix_d A, 
  ElConstMatrix_d BTL, ElConstMatrix_d BTR, 
  ElConstMatrix_d BBL, ElConstMatrix_d BBR );
EL_EXPORT ElError ElLockedMerge2x2_c
( ElMatrix_c A, 
  ElConstMatrix_c BTL, ElConstMatrix_c BTR, 
  ElConstMatrix_c BBL, ElConstMatrix_c BBR );
EL_EXPORT ElError ElLockedMerge2x2_z
( ElMatrix_z A, 
  ElConstMatrix_z BTL, ElConstMatrix_z BTR, 
  ElConstMatrix_z BBL, ElConstMatrix_z BBR );
EL_EXPORT ElError ElLockedMerge2x2Dist_i
( ElDistMatrix_i A, 
  ElConstDistMatrix_i BTL, ElConstDistMatrix_i BTR, 
  ElConstDistMatrix_i BBL, ElConstDistMatrix_i BBR );
EL_EXPORT ElError ElLockedMerge2x2Dist_s
( ElDistMatrix_s A, 
  ElConstDistMatrix_s BTL, ElConstDistMatrix_s BTR, 
  ElConstDistMatrix_s BBL, ElConstDistMatrix_s BBR );
EL_EXPORT ElError ElLockedMerge2x2Dist_d
( ElDistMatrix_d A, 
  ElConstDistMatrix_d BTL, ElConstDistMatrix_d BTR, 
  ElConstDistMatrix_d BBL, ElConstDistMatrix_d BBR );
EL_EXPORT ElError ElLockedMerge2x2Dist_c
( ElDistMatrix_c A, 
  ElConstDistMatrix_c BTL, ElConstDistMatrix_c BTR, 
  ElConstDistMatrix_c BBL, ElConstDistMatrix_c BBR );
EL_EXPORT ElError ElLockedMerge2x2Dist_z
( ElDistMatrix_z A, 
  ElConstDistMatrix_z BTL, ElConstDistMatrix_z BTR, 
  ElConstDistMatrix_z BBL, ElConstDistMatrix_z BBR );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_FLAMEPART_MERGE_C_H */
