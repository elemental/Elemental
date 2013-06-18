/*
   Copyright (c) 2009-2013, Jack Poulson
                      2013, Michael C. Grant
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

%module elem_matrices

%include "common.swg"
%import "elem.i"

/*
 * SPECIAL MATRICES
 */

// The equivalent of elemental/matrices.hpp
%include "elemental/matrices/Cauchy.hpp"
%include "elemental/matrices/CauchyLike.hpp"
%include "elemental/matrices/Circulant.hpp"
%include "elemental/matrices/Diagonal.hpp"
%include "elemental/matrices/Egorov.hpp"
%include "elemental/matrices/ExtendedKahan.hpp"
%include "elemental/matrices/Fiedler.hpp"
%include "elemental/matrices/Forsythe.hpp"
%include "elemental/matrices/Fourier.hpp"
%include "elemental/matrices/GCDMatrix.hpp"
%include "elemental/matrices/Gear.hpp"
%include "elemental/matrices/GKS.hpp"
%include "elemental/matrices/Grcar.hpp"
%include "elemental/matrices/Hankel.hpp"
%include "elemental/matrices/Hanowa.hpp"
%include "elemental/matrices/Helmholtz.hpp"
%include "elemental/matrices/HermitianUniformSpectrum.hpp"
%include "elemental/matrices/Hilbert.hpp"
%include "elemental/matrices/Identity.hpp"
%include "elemental/matrices/Jordan.hpp"
%include "elemental/matrices/Kahan.hpp"
%include "elemental/matrices/KMS.hpp"
%include "elemental/matrices/Laplacian.hpp"
%include "elemental/matrices/Lauchli.hpp"
%include "elemental/matrices/Legendre.hpp"
%include "elemental/matrices/Lehmer.hpp"
%include "elemental/matrices/Lotkin.hpp"
%include "elemental/matrices/MinIJ.hpp"
%include "elemental/matrices/NormalUniformSpectrum.hpp"
%include "elemental/matrices/Ones.hpp"
%include "elemental/matrices/OneTwoOne.hpp"
%include "elemental/matrices/Parter.hpp"
%include "elemental/matrices/Pei.hpp"
%include "elemental/matrices/Redheffer.hpp"
%include "elemental/matrices/Riemann.hpp"
%include "elemental/matrices/Ris.hpp"
%include "elemental/matrices/Toeplitz.hpp"
%include "elemental/matrices/TriW.hpp"
%include "elemental/matrices/Uniform.hpp"
%include "elemental/matrices/Walsh.hpp"
%include "elemental/matrices/Wilkinson.hpp"
%include "elemental/matrices/Zeros.hpp"

namespace elem {
OVERLOAD01(Cauchy)
OVERLOAD01(CauchyLike)
OVERLOAD01_int(Circulant)
OVERLOAD01_int(Diagonal)
// Not sure how to handle Egorov yet...
OVERLOAD01(ExtendedKahan)
OVERLOAD01_int(Fiedler)
OVERLOAD01_int(Forsythe)
OVERLOAD01_cpx(Fourier)
OVERLOAD01_int(GCDMatrix)
OVERLOAD01_int(Gear)
OVERLOAD01(GKS)
OVERLOAD01_int(Grcar)
OVERLOAD01_int(Hankel)
OVERLOAD01_int(Hanowa)
OVERLOAD01(Helmholtz)
OVERLOAD01(HermitianUniformSpectrum)
OVERLOAD01(Hilbert)
OVERLOAD01_int(Identity)
OVERLOAD01_int(Jordan)
OVERLOAD01(Kahan)
OVERLOAD01_int(KMS)
OVERLOAD01(Laplacian)
OVERLOAD01_int(Lauchli)
OVERLOAD01(Legendre)
OVERLOAD01(Lehmer)
OVERLOAD01(Lotkin)
OVERLOAD01_int(MinIJ)
OVERLOAD01_cpx(NormalUniformSpectrum)
OVERLOAD01_int(Ones)
OVERLOAD01_int(OneTwoOne)
OVERLOAD01(Parter)
OVERLOAD01_int(Pei)
OVERLOAD01_int(Redheffer)
OVERLOAD01_int(Riemann)
OVERLOAD01(Ris)
OVERLOAD01_int(Toeplitz)
OVERLOAD01_int(TriW)
OVERLOAD01_int(Uniform)
OVERLOAD1_int_UV(Uniform,CIRC,CIRC)
OVERLOAD1_int_UV(Uniform,MC,STAR)
OVERLOAD1_int_UV(Uniform,MR,MC)
OVERLOAD1_int_UV(Uniform,MR,STAR)
OVERLOAD1_int_UV(Uniform,STAR,MC)
OVERLOAD1_int_UV(Uniform,STAR,MR)
OVERLOAD1_int_UV(Uniform,STAR,STAR)
OVERLOAD1_int_UV(Uniform,STAR,VC)
OVERLOAD1_int_UV(Uniform,STAR,VR)
OVERLOAD1_int_UV(Uniform,VC,STAR)
OVERLOAD1_int_UV(Uniform,VR,STAR)
OVERLOAD01_int(Walsh)
OVERLOAD01_int(Wilkinson)
OVERLOAD01_int(Zeros)
OVERLOAD1_int_UV(Zeros,CIRC,CIRC)
OVERLOAD1_int_UV(Zeros,MC,STAR)
OVERLOAD1_int_UV(Zeros,MR,MC)
OVERLOAD1_int_UV(Zeros,MR,STAR)
OVERLOAD1_int_UV(Zeros,STAR,MC)
OVERLOAD1_int_UV(Zeros,STAR,MR)
OVERLOAD1_int_UV(Zeros,STAR,STAR)
OVERLOAD1_int_UV(Zeros,STAR,VC)
OVERLOAD1_int_UV(Zeros,STAR,VR)
OVERLOAD1_int_UV(Zeros,VC,STAR)
OVERLOAD1_int_UV(Zeros,VR,STAR)

};
