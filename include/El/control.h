/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_CONTROL_C_H
#define EL_CONTROL_C_H

#ifdef __cplusplus
extern "C" {
#endif

/* Lyapunov
   ======== */
EL_EXPORT ElError ElLyapunov_s
( ElConstMatrix_s A, ElConstMatrix_s C, ElMatrix_s X );
EL_EXPORT ElError ElLyapunov_d
( ElConstMatrix_d A, ElConstMatrix_d C, ElMatrix_d X );
EL_EXPORT ElError ElLyapunov_c
( ElConstMatrix_c A, ElConstMatrix_c C, ElMatrix_c X );
EL_EXPORT ElError ElLyapunov_z
( ElConstMatrix_z A, ElConstMatrix_z C, ElMatrix_z X );

EL_EXPORT ElError ElLyapunovDist_s
( ElConstDistMatrix_s A, ElConstDistMatrix_s C, ElDistMatrix_s X );
EL_EXPORT ElError ElLyapunovDist_d
( ElConstDistMatrix_d A, ElConstDistMatrix_d C, ElDistMatrix_d X );
EL_EXPORT ElError ElLyapunovDist_c
( ElConstDistMatrix_c A, ElConstDistMatrix_c C, ElDistMatrix_c X );
EL_EXPORT ElError ElLyapunovDist_z
( ElConstDistMatrix_z A, ElConstDistMatrix_z C, ElDistMatrix_z X );

/* TODO: Expert versions */

/* Ricatti
   ======= */
EL_EXPORT ElError ElRicatti_s
( ElUpperOrLower uplo, 
  ElConstMatrix_s A, ElConstMatrix_s K, ElConstMatrix_s L, 
  ElMatrix_s X );
EL_EXPORT ElError ElRicatti_d
( ElUpperOrLower uplo, 
  ElConstMatrix_d A, ElConstMatrix_d K, ElConstMatrix_d L, 
  ElMatrix_d X );
EL_EXPORT ElError ElRicatti_c
( ElUpperOrLower uplo, 
  ElConstMatrix_c A, ElConstMatrix_c K, ElConstMatrix_c L, 
  ElMatrix_c X );
EL_EXPORT ElError ElRicatti_z
( ElUpperOrLower uplo, 
  ElConstMatrix_z A, ElConstMatrix_z K, ElConstMatrix_z L, 
  ElMatrix_z X );

EL_EXPORT ElError ElRicattiDist_s
( ElUpperOrLower uplo, 
  ElConstDistMatrix_s A, ElConstDistMatrix_s K, ElConstDistMatrix_s L, 
  ElDistMatrix_s X );
EL_EXPORT ElError ElRicattiDist_d
( ElUpperOrLower uplo, 
  ElConstDistMatrix_d A, ElConstDistMatrix_d K, ElConstDistMatrix_d L, 
  ElDistMatrix_d X );
EL_EXPORT ElError ElRicattiDist_c
( ElUpperOrLower uplo, 
  ElConstDistMatrix_c A, ElConstDistMatrix_c K, ElConstDistMatrix_c L, 
  ElDistMatrix_c X );
EL_EXPORT ElError ElRicattiDist_z
( ElUpperOrLower uplo, 
  ElConstDistMatrix_z A, ElConstDistMatrix_z K, ElConstDistMatrix_z L, 
  ElDistMatrix_z X );

EL_EXPORT ElError ElRicattiPreformed_s( ElMatrix_s W, ElMatrix_s X );
EL_EXPORT ElError ElRicattiPreformed_d( ElMatrix_d W, ElMatrix_d X );
EL_EXPORT ElError ElRicattiPreformed_c( ElMatrix_c W, ElMatrix_c X );
EL_EXPORT ElError ElRicattiPreformed_z( ElMatrix_z W, ElMatrix_z X );

EL_EXPORT ElError ElRicattiPreformedDist_s
( ElDistMatrix_s W, ElDistMatrix_s X );
EL_EXPORT ElError ElRicattiPreformedDist_d
( ElDistMatrix_d W, ElDistMatrix_d X );
EL_EXPORT ElError ElRicattiPreformedDist_c
( ElDistMatrix_c W, ElDistMatrix_c X );
EL_EXPORT ElError ElRicattiPreformedDist_z
( ElDistMatrix_z W, ElDistMatrix_z X );

/* TODO: Expert versions */

/* Sylvester
   ========= */
EL_EXPORT ElError ElSylvester_s
( ElConstMatrix_s A, ElConstMatrix_s B, ElConstMatrix_s C,
  ElMatrix_s X );
EL_EXPORT ElError ElSylvester_d
( ElConstMatrix_d A, ElConstMatrix_d B, ElConstMatrix_d C,
  ElMatrix_d X );
EL_EXPORT ElError ElSylvester_c
( ElConstMatrix_c A, ElConstMatrix_c B, ElConstMatrix_c C,
  ElMatrix_c X );
EL_EXPORT ElError ElSylvester_z
( ElConstMatrix_z A, ElConstMatrix_z B, ElConstMatrix_z C,
  ElMatrix_z X );

EL_EXPORT ElError ElSylvesterDist_s
( ElConstDistMatrix_s A, ElConstDistMatrix_s B, ElConstDistMatrix_s C,
  ElDistMatrix_s X );
EL_EXPORT ElError ElSylvesterDist_d
( ElConstDistMatrix_d A, ElConstDistMatrix_d B, ElConstDistMatrix_d C,
  ElDistMatrix_d X );
EL_EXPORT ElError ElSylvesterDist_c
( ElConstDistMatrix_c A, ElConstDistMatrix_c B, ElConstDistMatrix_c C,
  ElDistMatrix_c X );
EL_EXPORT ElError ElSylvesterDist_z
( ElConstDistMatrix_z A, ElConstDistMatrix_z B, ElConstDistMatrix_z C,
  ElDistMatrix_z X );

EL_EXPORT ElError ElSylvesterPreformed_s( ElInt m, ElMatrix_s W, ElMatrix_s X );
EL_EXPORT ElError ElSylvesterPreformed_d( ElInt m, ElMatrix_d W, ElMatrix_d X );
EL_EXPORT ElError ElSylvesterPreformed_c( ElInt m, ElMatrix_c W, ElMatrix_c X );
EL_EXPORT ElError ElSylvesterPreformed_z( ElInt m, ElMatrix_z W, ElMatrix_z X );

EL_EXPORT ElError ElSylvesterPreformedDist_s
( ElInt m, ElDistMatrix_s W, ElDistMatrix_s X );
EL_EXPORT ElError ElSylvesterPreformedDist_d
( ElInt m, ElDistMatrix_d W, ElDistMatrix_d X );
EL_EXPORT ElError ElSylvesterPreformedDist_c
( ElInt m, ElDistMatrix_c W, ElDistMatrix_c X );
EL_EXPORT ElError ElSylvesterPreformedDist_z
( ElInt m, ElDistMatrix_z W, ElDistMatrix_z X );

/* TODO: Expert versions */

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_CONTROL_C_H */
