/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_VIEWS_C_H
#define EL_VIEWS_C_H

#ifdef __cplusplus
extern "C" {
#endif

/* NOTE: 'cc' 4.8.2 in Ubuntu fails if "ElRange_i iRange" is replaced with
         "ElRange_i iRange" */

ElError ElView_i
( ElMatrix_i A, ElMatrix_i B, ElRange_i iRange, ElRange_i jRange );
ElError ElView_s
( ElMatrix_s A, ElMatrix_s B, ElRange_i iRange, ElRange_i jRange );
ElError ElView_d
( ElMatrix_d A, ElMatrix_d B, ElRange_i iRange, ElRange_i jRange );
ElError ElView_c
( ElMatrix_c A, ElMatrix_c B, ElRange_i iRange, ElRange_i jRange );
ElError ElView_z
( ElMatrix_z A, ElMatrix_z B, ElRange_i iRange, ElRange_i jRange );
ElError ElViewDist_i
( ElDistMatrix_i A, ElDistMatrix_i B, ElRange_i iRange, ElRange_i jRange );
ElError ElViewDist_s
( ElDistMatrix_s A, ElDistMatrix_s B, ElRange_i iRange, ElRange_i jRange );
ElError ElViewDist_d
( ElDistMatrix_d A, ElDistMatrix_d B, ElRange_i iRange, ElRange_i jRange );
ElError ElViewDist_c
( ElDistMatrix_c A, ElDistMatrix_c B, ElRange_i iRange, ElRange_i jRange );
ElError ElViewDist_z
( ElDistMatrix_z A, ElDistMatrix_z B, ElRange_i iRange, ElRange_i jRange );

ElError ElLockedView_i
( ElMatrix_i A, ElConstMatrix_i B, ElRange_i iRange, ElRange_i jRange );
ElError ElLockedView_s
( ElMatrix_s A, ElConstMatrix_s B, ElRange_i iRange, ElRange_i jRange );
ElError ElLockedView_d
( ElMatrix_d A, ElConstMatrix_d B, ElRange_i iRange, ElRange_i jRange );
ElError ElLockedView_c
( ElMatrix_c A, ElConstMatrix_c B, ElRange_i iRange, ElRange_i jRange );
ElError ElLockedView_z
( ElMatrix_z A, ElConstMatrix_z B, ElRange_i iRange, ElRange_i jRange );
ElError ElLockedViewDist_i
( ElDistMatrix_i A, ElConstDistMatrix_i B, ElRange_i iRange, ElRange_i jRange );
ElError ElLockedViewDist_s
( ElDistMatrix_s A, ElConstDistMatrix_s B, ElRange_i iRange, ElRange_i jRange );
ElError ElLockedViewDist_d
( ElDistMatrix_d A, ElConstDistMatrix_d B, ElRange_i iRange, ElRange_i jRange );
ElError ElLockedViewDist_c
( ElDistMatrix_c A, ElConstDistMatrix_c B, ElRange_i iRange, ElRange_i jRange );
ElError ElLockedViewDist_z
( ElDistMatrix_z A, ElConstDistMatrix_z B, ElRange_i iRange, ElRange_i jRange );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_VIEWS_C_H */
