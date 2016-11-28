#
#  Copyright (c) 2009-2016, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
import math, time, El
if El.havePyPlot:
  El.plt.set_cmap('bone')  
  # pyplot.set_cmap seems to open an empty figure (which can be detected by
  # running pyplot.get_fignums()), and so we manually close it. Unfortunately,
  # some versions of pyplot generate a spurious warning of:
  #
  #   can't invoke "event" command: application has been destroyed
  #
  # when calling pyplot.close().
  El.plt.close(1)

realRes = imagRes = 100 # grid resolution

# Display an instance of the Fox-Li/Landau matrix
A = El.DistMatrix(El.zTag)

nList = (50,100,300)
for n in nList:
  El.FoxLi(A,n)

  # Show the Real part of the matrix
  AReal = El.DistMatrix(El.dTag)
  El.RealPart(A,AReal)
  El.Display(AReal,'Real part of Fox-Li matrix (n={})'.format(n))

  # Compute the Schur decomposition (overwriting A with the Schur factor)
  schurStart = time.time()
  w = El.Schur(A)
  schurStop = time.time()
  if A.Grid().Rank() is 0:
    print('Schur decomp for n={}: {} [sec]'.format(n,schurStop-schurStart,))
  
  # Compute the spectral portrait of the Schur factor
  portraitStart = time.time()
  portrait, box = El.TriangularSpectralPortrait(A,realRes,imagRes)
  portraitStop = time.time()
  if A.Grid().Rank() is 0:
    print('Portrait for n={}: {} [sec]'.format(n,portraitStop-portraitStart,))

  # Display the eigenvalues on top of log10 plot of portrait
  El.DisplayPortrait(portrait,box,title='Fox-Li portrait (n={})'.format(n),
                     eigvals=w)

El.Finalize()
