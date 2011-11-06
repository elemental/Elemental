Parallel LCG
============
Since it is often necessary to generate a large matrix with pseudo-random 
entries in parallel, a method for ensuring that a large set of processes can 
each generate independent uniformly random samples is required. The purpose of
Parallel LCG (PLCG) is to provide a provably independent generalization of a
simple (but well-studied) Linear Congruential Generator. Knuth's constants from
The Art of Computer Programming Vol. 2 are used.

Datatypes
---------

.. cpp:type:: plcg::UInt32

   Since the vast majority of modern systems make use of ``unsigned`` for
   storing 32-bit unsigned integers, we simply hardcode the type. If your 
   system does not follow this convention, then this typedef will need to be
   changed!

.. cpp:type:: struct plcg::UInt64

   A custom 64-bit unsigned integer which is simply the concatenation of two 
   32-bit unsigned integers (``UInt32``).

.. cpp:type:: struct plcg::ExpandedUInt64

   A custom 64-bit unsigned integer which is stores each of the four 16-bit
   pieces within the first 16 bits of a 32-bit unsigned integer. This is done
   so that two such expanded 16-bit numbers can be multiplied without any 
   chance of overflow.

LCG primitives
--------------

.. cpp:function:: plcg::UInt32 plcg::Lower16Bits( plcg::UInt32 a )

   Return the lower 16 bits of ``a`` in the lower 16 bits of the returned 
   32-bit unsigned integer.

.. cpp:function:: plcg::UInt32 plcg::Upper16Bits( plcg::UInt32 a )

   Return the upper 16 bits of ``a`` in the lower 16 bits of the returned
   32-bit unsigned integer.

.. cpp:function:: plcg::ExpandedUInt64 plcg::Expand( plcg::UInt32 a )

   Expand a 32-bit unsigned integer into a 64-bit expanded representation.

.. cpp:function:: plcg::ExpandedUInt64 plcg::Expand( plcg::UInt64 a )

   Expand a 64-bit unsigned integer into a 64-bit expanded representation.

.. cpp:function:: plcg::UInt64 plcg::Deflate( plcg::ExpandedUInt64 a )

   Deflate an expanded 64-bit unsigned integer into the standard 64-bit form.

.. cpp:function:: void plcg::CarryUpper16Bits( plcg::ExpandedUInt64& a )

   Carry the results stored in the upper 16-bits of each of the four pieces 
   into the next lower 16 bits.

.. cpp:function:: plcg::ExpandedUInt64 plcg::AddWith64BitMod( plcg::ExpandedUInt64 a, plcg::ExpandedUInt64 b )

   Return :math:`a+b \mod 2^{64}`.

.. cpp:function:: plcg::ExpandedUInt64 plcg::MultiplyWith64BitMod( plcg::ExpandedUInt64 a, plcg::ExpandedUInt64 b )

   Return :math:`ab \mod 2^{64}`.

.. cpp:function:: plcg::ExpandedUInt64 plcg::IntegerPowerWith64BitMod( plcg::ExpandedUInt64 x, plcg::ExpandedUInt64 n )

   Return :math:`x^n \mod 2^{64}`.

.. cpp:function:: void plcg::Halve( plcg::ExpandedUInt64& a )

   :math:`a := a/2`.

.. cpp:function:: void plcg::SeedSerialLcg( plcg::UInt64 globalSeed )

   Set the initial state of the serial Linear Congruential Generator.

.. cpp:function:: void plcg::SeedParallelLcg( plcg::UInt32 rank, plcg::UInt32 commSize, plcg::UInt64 globalSeed )

   Have our process seed a separate LCG meant for parallel computation, where 
   the calling process has the given rank within a communicator of the 
   specified size.

.. cpp:function:: plcg::UInt64 plcg::SerialLcg()

   Return the current state of the serial LCG, and then advance to the next one.

.. cpp:function:: plcg::UInt64 plcg::ParallelLcg()

   Return the current state of our process's portion of the parallel LCG, 
   and then advance to our next local state.

.. cpp:function:: void plcg::ManualLcg( plcg::ExpandedUInt64 a, plcg::ExpandedUInt64 c, plcg::ExpandedUInt64& X )

   :math:`X := a X + c \mod 2^{64}`.

Sampling
--------

.. cpp:function:: R plcg::SerialUniform()

   Return a uniform sample from :math:`(0,1]` using the serial LCG.

.. cpp:function:: R plcg::ParallelUniform()

   Return a uniform sample from :math:`(0,1]` using the parallel LCG.

.. cpp:function:: void plcg::SerialBoxMuller( R& X, R& Y )

   Return two samples from a normal distribution with mean 0 and standard 
   deviation of 1 using the serial LCG.

.. cpp:function:: void plcg::ParallelBoxMuller( R& X, R& Y )

   Return two samples from a normal distribution with mean 0 and standard
   deviation 1, but using the parallel LCG.

.. cpp:function:: void plcg::SerialGaussianRandomVariable( R& X )

   Return a single sample from a normal distribution with mean 0 and 
   standard deviation 1 using the serial LCG.

.. cpp:function:: void plcg::ParallelGaussianRandomVariable( R& X )

   Return a single sample from a normal distribution with mean 0 and 
   standard deviation 1, but using the parallel LCG.
   
.. cpp:function:: void plcg::SerialGaussianRandomVariable( std::complex<R>& X )

   Return a single complex sample from a normal distribution with mean 0 and 
   standard deviation 1 using the serial LCG.

.. cpp:function:: void plcg::ParallelGaussianRandomVariable( std::complex<R>& X )

   Return a single complex sample from a normal distribution with mean 0 and 
   standard deviation 1, but using the parallel LCG.
