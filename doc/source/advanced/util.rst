Utilities
=========

Householder reflectors
----------------------
**TODO:** Describe major difference from LAPACK's conventions (i.e., we do not 
treat the identity matrix as a Householder transform since it requires the 
:math:`u` in :math:`H=I-2uu'` to have norm zero rather than one). 

Applying packed Householder transforms
--------------------------------------
**TODO:** Describe ``advanced::ApplyPackedReflectors`` here.

Applying pivots
---------------
**TODO**

Reduction of Hermitian generalized-definite EVPs
------------------------------------------------
**TODO:** Describe the reduction steps of :math:`ABx=\lambda x`, 
:math:`BAx=\lambda x`, and :math:`Ax=\lambda Bx` using the operations 
:math:`A := L^H A L` and :math:`A := L^{-1} A L^{-H}`.

