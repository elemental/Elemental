Special matrices
****************

It is frequently useful to test algorithms on well-known, trivial, and random 
matrices, such as::

1. matrices with entries sampled from a uniform distribution,
2. matrices with spectrum sampled from a uniform distribution,
3. Wilkinson matrices,
4. identity matrices,
5. matrices of all ones, and
6. matrices of all zeros.

Elemental therefore provides utilities for generating many such matrices.

.. toctree::
   :maxdepth: 2

   special_matrices/deterministic
   special_matrices/random
