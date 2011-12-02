PMRRR
=====
Rather than directly using Petschow and Bientinesi's parallel implementation of 
the Multiple Relatively Robust Representations (MRRR) algorithm, several 
simplified interfaces have been exposed.

The prototypes can be found in
`include/elemental/imports/pmrrr.hpp <../../../../include/elemental/imports/pmrrr.hpp>`_,
while the implementations are in the folder
`external/pmrrr/ <../../../../external/pmrrr/>`_.

Data structures
---------------

.. cpp:type:: struct pmrrr::Estimate

   For returning upper bounds on the number of local and global eigenvalues
   with eigenvalues lying in the specified interval, :math:`(a,b]`.

   .. cpp:member:: int numLocalEigenvalues

      The upper bound on the number of eigenvalues in the specified interval 
      that our process stores locally.

   .. cpp:member:: int numGlobalEigenvalues

      The upper bound on the number of eigenvalues in the specified interval.

.. cpp:type:: struct pmrrr::Info

   For returning information about the computed eigenvalues.

   .. cpp:member:: int numLocalEigenvalues

      The number of computed eigenvalues that our process locally stores.

   .. cpp:member:: int numGlobalEigenvalues

      The number of computed eigenvalues.

   .. cpp:member:: int firstLocalEigenvalue

      The index of the first eigenvalue stored locally on our process.

Compute all eigenvalues
-----------------------

.. cpp:function:: pmrrr::Info pmrrr::Eig( int n, const double* d, const double* e, double* w, mpi::Comm comm )

   Compute all of the eigenvalues of the real symmetric tridiagonal matrix with 
   diagonal `d` and subdiagonal `e`: the eigenvalues will be stored in 
   `w` and the work will be divided among the processors in `comm`.

.. cpp:function:: pmrrr::Info pmrrr::Eig( int n, const double* d, const double* e, double* w, double* Z, int ldz, mpi::Comm comm )

   Same as above, but also compute the corresponding eigenvectors.

Compute eigenvalues within interval
-----------------------------------

.. cpp:function:: pmrrr::Info pmrrr::Eig( int n, const double* d, const double* e, double* w, mpi::Comm comm, double a, double b )

   Only compute the eigenvalues lying within the interval :math:`(a,b]`.

.. cpp:function:: pmrrr::Info pmrrr::Eig( int n, const double* d, const double* e, double* w, double* Z, int ldz, mpi::Comm comm, double a, double b )

   Same as above, but also compute the corresponding eigenvectors.

.. cpp:function:: pmrrr::Estimate pmrrr::EigEstimate( int n, const double* d, const double* w, mpi::Comm comm, double a, double b )

   Return upper bounds on the number of local and global eigenvalues lying 
   within the specified interval.

Compute eigenvalues in index range
----------------------------------

.. cpp:function:: pmrrr::Info pmrrr::Eig( int n, const double* d, const double* e, double* w, mpi::Comm comm, int a, int b )

   Only compute the eigenvalues with indices ranging from ``a`` to ``b``, where 
   :math:`0 \le a \le b < n`.

.. cpp:function:: pmrrr::Info pmrrr::Eig( int n, const double* d, const double* e, double* w, double* Z, int ldz, mpi::Comm comm, int a, int b )

   Same as above, but also compute the corresponding eigenvectors.
