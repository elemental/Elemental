libFLAME
--------

.. cpp:function:: int FLA_Bsvd_v_opd_var1( int k, int mU, int mV, int nGH, int nIterMax, double* d, int dInc, double* e, int eInc, Complex<double>* G, int rsG, int csG, Complex<double>* H, int rsH, int csH, double* U, int rsU, int csU, double* V, int rsV, int csV, int nb )
.. cpp:function:: int FLA_Bsvd_v_opd_var1( int k, int mU, int mV, int nGH, int nIterMax, double* d, int dInc, double* e, int eInc, Complex<double>* G, int rsG, int csG, Complex<double>* H, int rsH, int csH, Complex<double>* U, int rsU, int csU, Complex<double>* V, int rsV, int csV, int nb )

Optional high-performance implementations of the bidiagonal QR algorithm.
This can lead to substantial improvements in Elemental's distributed-memory SVD
on supported architectures (as of now, modern Intel architectures).
