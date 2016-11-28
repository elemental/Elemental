#
#  Copyright (c) 2009-2016, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
from core          import *
from blas_like     import *
from io            import *
from lapack_like   import *
from matrices      import *
from optimization  import *
from control       import *
from lattice       import *

lib.ElFinalize.argtypes = []
def Finalize():
  if havePyPlot:
    if len(plt.get_fignums()) > 0 and not isInlinePyPlot:
      plt.show()
  lib.ElFinalize()
