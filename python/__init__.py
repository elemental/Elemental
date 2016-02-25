#
#  Copyright (c) 2009-2016, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#

# Attempt to import matplotlib.pyplot and save whether or not this succeeded
try:
  import numpy as np
  import matplotlib as mpl
  import matplotlib.pyplot as plt
  havePyPlot=True
except:
  havePyPlot=False 
  print 'Could not import matplotlib.pyplot'

if havePyPlot:
  try:
    import networkx as nx
    haveNetworkX = True
  except:
    haveNetworkX = False 
    print 'Could not import networkx'
else:
  haveNetworkX = False

from core          import *
from blas_like     import *
from io            import *
from lapack_like   import *
from matrices      import *
from optimization  import *
from control       import *
from lattice       import *
