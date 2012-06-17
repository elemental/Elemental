Tuning parameters
=================

Hermitian to tridiagonal
------------------------
Two different basic strategies are available for the reduction to tridiagonal
form:

1. Run a pipelined algorithm designed for general (rectangular) process grids.
2. Redistribute the matrix so that it is owned by a perfect square number of
   processes, perform a fast reduction to tridiaogal form, and redistribute
   the data back to the original process grid. This algorithm is essentially
   an evolution of the HJS tridiagonalization approach
   (see "*Towards an efficient parallel eigensolver for dense symmetric 
   matrices*" by Bruce Hendrickson, Elizabeth Jessup, and Christopher Smith)
   which is described in detail in Ken Stanley's dissertation, "*Execution time 
   of symmetric eigensolvers*". 

There is clearly a small penalty associated with the extra redistributions
necessary for the second approach, but the benefit from using a square process
grid is usually quite signficant. By default, ``HermitianTridiag`` will run the
standard algorithm (approach 1) unless the matrix is already distributed over a
square process grid. The reasoning is that good performance depends upon a
"good" ordering of the square (say, :math:`\hat p \times \hat p`) subgrid,
though usually either a row-major or column-major ordering of the first
:math:`\hat p^2` processes suffices.


.. cpp:type:: HermitianTridiagApproach

   * ``HERMITIAN_TRIDIAG_NORMAL``: Run the pipelined rectangular algorithm.
   * ``HERMITIAN_TRIDIAG_SQUARE``: Run the square grid algorithm on the largest
     possible square process grid.
   * ``HERMITIAN_TRIDIAG_DEFAULT``: If the given process grid is already square,
     run the square grid algorithm, otherwise use the pipelined non-square
     approach.

   .. note::

      A properly tuned ``HERMITIAN_TRIDIAG_SQUARE`` approach is almost always 
      fastest, so it is worthwhile to test it with both the ``COLUMN_MAJOR`` and 
      ``ROW_MAJOR`` subgrid orderings, as described below.


   .. note::
   
      The first algorithm heavily depends upon the performance of distributed 
      ``Symv`` and ``Hemv`` (for real and complex data, 
      respectively), so users interested in maximizing the performance of the 
      first algorithm will likely want to investigate different values for the 
      local blocksizes through the routines 
      ``SetLocalSymvBlocksize<T>( int blocksize )`` and 
      ``SetLocalHemvBlocksize<T>( int blocksize )``; the default values 
      are both 64.

.. cpp:function:: void SetHermitianTridiagApproach( HermitianTridiagApproach approach )

   Sets the algorithm used by subsequent calls to
   ``HermitianTridiag``.

.. cpp:function:: HermitianTridiagApproach GetHermitianTridiagApproach()

   Queries the currently set approach for the reduction of a Hermitian matrix
   to tridiagonal form.

.. cpp:function:: void SetHermitianTridiagGridOrder( GridOrder order )

   Sets the ordering to use for the first :math:`\hat p^2` processes in the
   construction of the :math:`\hat p \times \hat p` subgrid. This is only
   relevant to the ``HERMITIAN_TRIDIAG_SQUARE`` approach.

.. cpp:function:: GridOrder GetHermitianTridiagGridOrder()

   Queries the currently set approach for the ordering of the square subgrid
   needed by the ``HERMITIAN_TRIDIAG_SQUARE`` approach to the
   tridiagonalization of a Hermitian matrix.

