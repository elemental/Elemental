#
#  Copyright (c) 2009-2015, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
from environment import *
from imports     import mpi

import DistGraph as DG

# DistSparseMatrix
# ================

lib.ElDistSparseMatrixCreate_i.argtypes = [POINTER(c_void_p),mpi.Comm]
lib.ElDistSparseMatrixCreate_i.restype = c_uint
lib.ElDistSparseMatrixCreate_s.argtypes = [POINTER(c_void_p),mpi.Comm]
lib.ElDistSparseMatrixCreate_s.restype = c_uint
lib.ElDistSparseMatrixCreate_d.argtypes = [POINTER(c_void_p),mpi.Comm]
lib.ElDistSparseMatrixCreate_d.restype = c_uint
lib.ElDistSparseMatrixCreate_c.argtypes = [POINTER(c_void_p),mpi.Comm]
lib.ElDistSparseMatrixCreate_c.restype = c_uint
lib.ElDistSparseMatrixCreate_z.argtypes = [POINTER(c_void_p),mpi.Comm]
lib.ElDistSparseMatrixCreate_z.restype = c_uint

lib.ElDistSparseMatrixDestroy_i.argtypes = [c_void_p]
lib.ElDistSparseMatrixDestroy_i.restype = c_uint
lib.ElDistSparseMatrixDestroy_s.argtypes = [c_void_p]
lib.ElDistSparseMatrixDestroy_s.restype = c_uint
lib.ElDistSparseMatrixDestroy_d.argtypes = [c_void_p]
lib.ElDistSparseMatrixDestroy_d.restype = c_uint
lib.ElDistSparseMatrixDestroy_c.argtypes = [c_void_p]
lib.ElDistSparseMatrixDestroy_c.restype = c_uint
lib.ElDistSparseMatrixDestroy_z.argtypes = [c_void_p]
lib.ElDistSparseMatrixDestroy_z.restype = c_uint

lib.ElDistSparseMatrixEmpty_i.argtypes = [c_void_p]
lib.ElDistSparseMatrixEmpty_i.restype = c_uint
lib.ElDistSparseMatrixEmpty_s.argtypes = [c_void_p]
lib.ElDistSparseMatrixEmpty_s.restype = c_uint
lib.ElDistSparseMatrixEmpty_d.argtypes = [c_void_p]
lib.ElDistSparseMatrixEmpty_d.restype = c_uint
lib.ElDistSparseMatrixEmpty_c.argtypes = [c_void_p]
lib.ElDistSparseMatrixEmpty_c.restype = c_uint
lib.ElDistSparseMatrixEmpty_z.argtypes = [c_void_p]
lib.ElDistSparseMatrixEmpty_z.restype = c_uint

lib.ElDistSparseMatrixResize_i.argtypes = [c_void_p,iType,iType]
lib.ElDistSparseMatrixResize_i.restype = c_uint
lib.ElDistSparseMatrixResize_s.argtypes = [c_void_p,iType,iType]
lib.ElDistSparseMatrixResize_s.restype = c_uint
lib.ElDistSparseMatrixResize_d.argtypes = [c_void_p,iType,iType]
lib.ElDistSparseMatrixResize_d.restype = c_uint
lib.ElDistSparseMatrixResize_c.argtypes = [c_void_p,iType,iType]
lib.ElDistSparseMatrixResize_c.restype = c_uint
lib.ElDistSparseMatrixResize_z.argtypes = [c_void_p,iType,iType]
lib.ElDistSparseMatrixResize_z.restype = c_uint

lib.ElDistSparseMatrixSetComm_i.argtypes = [c_void_p,mpi.Comm]
lib.ElDistSparseMatrixSetComm_i.restype = c_uint
lib.ElDistSparseMatrixSetComm_s.argtypes = [c_void_p,mpi.Comm]
lib.ElDistSparseMatrixSetComm_s.restype = c_uint
lib.ElDistSparseMatrixSetComm_d.argtypes = [c_void_p,mpi.Comm]
lib.ElDistSparseMatrixSetComm_d.restype = c_uint
lib.ElDistSparseMatrixSetComm_c.argtypes = [c_void_p,mpi.Comm]
lib.ElDistSparseMatrixSetComm_c.restype = c_uint
lib.ElDistSparseMatrixSetComm_z.argtypes = [c_void_p,mpi.Comm]
lib.ElDistSparseMatrixSetComm_z.restype = c_uint

lib.ElDistSparseMatrixReserve_i.argtypes = [c_void_p,iType]
lib.ElDistSparseMatrixReserve_i.restype = c_uint
lib.ElDistSparseMatrixReserve_s.argtypes = [c_void_p,iType]
lib.ElDistSparseMatrixReserve_s.restype = c_uint
lib.ElDistSparseMatrixReserve_d.argtypes = [c_void_p,iType]
lib.ElDistSparseMatrixReserve_d.restype = c_uint
lib.ElDistSparseMatrixReserve_c.argtypes = [c_void_p,iType]
lib.ElDistSparseMatrixReserve_c.restype = c_uint
lib.ElDistSparseMatrixReserve_z.argtypes = [c_void_p,iType]
lib.ElDistSparseMatrixReserve_z.restype = c_uint

lib.ElDistSparseMatrixUpdate_i.argtypes = [c_void_p,iType,iType,iType]
lib.ElDistSparseMatrixUpdate_i.restype = c_uint
lib.ElDistSparseMatrixUpdate_s.argtypes = [c_void_p,iType,iType,sType]
lib.ElDistSparseMatrixUpdate_s.restype = c_uint
lib.ElDistSparseMatrixUpdate_d.argtypes = [c_void_p,iType,iType,dType]
lib.ElDistSparseMatrixUpdate_d.restype = c_uint
lib.ElDistSparseMatrixUpdate_c.argtypes = [c_void_p,iType,iType,cType]
lib.ElDistSparseMatrixUpdate_c.restype = c_uint
lib.ElDistSparseMatrixUpdate_z.argtypes = [c_void_p,iType,iType,zType]
lib.ElDistSparseMatrixUpdate_z.restype = c_uint

lib.ElDistSparseMatrixUpdateLocal_i.argtypes = [c_void_p,iType,iType,iType]
lib.ElDistSparseMatrixUpdateLocal_i.restype = c_uint
lib.ElDistSparseMatrixUpdateLocal_s.argtypes = [c_void_p,iType,iType,sType]
lib.ElDistSparseMatrixUpdateLocal_s.restype = c_uint
lib.ElDistSparseMatrixUpdateLocal_d.argtypes = [c_void_p,iType,iType,dType]
lib.ElDistSparseMatrixUpdateLocal_d.restype = c_uint
lib.ElDistSparseMatrixUpdateLocal_c.argtypes = [c_void_p,iType,iType,cType]
lib.ElDistSparseMatrixUpdateLocal_c.restype = c_uint
lib.ElDistSparseMatrixUpdateLocal_z.argtypes = [c_void_p,iType,iType,zType]
lib.ElDistSparseMatrixUpdateLocal_z.restype = c_uint

lib.ElDistSparseMatrixZero_i.argtypes = [c_void_p,iType,iType]
lib.ElDistSparseMatrixZero_i.restype = c_uint
lib.ElDistSparseMatrixZero_s.argtypes = [c_void_p,iType,iType]
lib.ElDistSparseMatrixZero_s.restype = c_uint
lib.ElDistSparseMatrixZero_d.argtypes = [c_void_p,iType,iType]
lib.ElDistSparseMatrixZero_d.restype = c_uint
lib.ElDistSparseMatrixZero_c.argtypes = [c_void_p,iType,iType]
lib.ElDistSparseMatrixZero_c.restype = c_uint
lib.ElDistSparseMatrixZero_z.argtypes = [c_void_p,iType,iType]
lib.ElDistSparseMatrixZero_z.restype = c_uint

lib.ElDistSparseMatrixZeroLocal_i.argtypes = [c_void_p,iType,iType]
lib.ElDistSparseMatrixZeroLocal_i.restype = c_uint
lib.ElDistSparseMatrixZeroLocal_s.argtypes = [c_void_p,iType,iType]
lib.ElDistSparseMatrixZeroLocal_s.restype = c_uint
lib.ElDistSparseMatrixZeroLocal_d.argtypes = [c_void_p,iType,iType]
lib.ElDistSparseMatrixZeroLocal_d.restype = c_uint
lib.ElDistSparseMatrixZeroLocal_c.argtypes = [c_void_p,iType,iType]
lib.ElDistSparseMatrixZeroLocal_c.restype = c_uint
lib.ElDistSparseMatrixZeroLocal_z.argtypes = [c_void_p,iType,iType]
lib.ElDistSparseMatrixZeroLocal_z.restype = c_uint

lib.ElDistSparseMatrixQueueUpdate_i.argtypes = [c_void_p,iType,iType,iType]
lib.ElDistSparseMatrixQueueUpdate_i.restype = c_uint
lib.ElDistSparseMatrixQueueUpdate_s.argtypes = [c_void_p,iType,iType,sType]
lib.ElDistSparseMatrixQueueUpdate_s.restype = c_uint
lib.ElDistSparseMatrixQueueUpdate_d.argtypes = [c_void_p,iType,iType,dType]
lib.ElDistSparseMatrixQueueUpdate_d.restype = c_uint
lib.ElDistSparseMatrixQueueUpdate_c.argtypes = [c_void_p,iType,iType,cType]
lib.ElDistSparseMatrixQueueUpdate_c.restype = c_uint
lib.ElDistSparseMatrixQueueUpdate_z.argtypes = [c_void_p,iType,iType,zType]
lib.ElDistSparseMatrixQueueUpdate_z.restype = c_uint

lib.ElDistSparseMatrixQueueLocalUpdate_i.argtypes = [c_void_p,iType,iType,iType]
lib.ElDistSparseMatrixQueueLocalUpdate_i.restype = c_uint
lib.ElDistSparseMatrixQueueLocalUpdate_s.argtypes = [c_void_p,iType,iType,sType]
lib.ElDistSparseMatrixQueueLocalUpdate_s.restype = c_uint
lib.ElDistSparseMatrixQueueLocalUpdate_d.argtypes = [c_void_p,iType,iType,dType]
lib.ElDistSparseMatrixQueueLocalUpdate_d.restype = c_uint
lib.ElDistSparseMatrixQueueLocalUpdate_c.argtypes = [c_void_p,iType,iType,cType]
lib.ElDistSparseMatrixQueueLocalUpdate_c.restype = c_uint
lib.ElDistSparseMatrixQueueLocalUpdate_z.argtypes = [c_void_p,iType,iType,zType]
lib.ElDistSparseMatrixQueueLocalUpdate_z.restype = c_uint

lib.ElDistSparseMatrixQueueZero_i.argtypes = [c_void_p,iType,iType]
lib.ElDistSparseMatrixQueueZero_i.restype = c_uint
lib.ElDistSparseMatrixQueueZero_s.argtypes = [c_void_p,iType,iType]
lib.ElDistSparseMatrixQueueZero_s.restype = c_uint
lib.ElDistSparseMatrixQueueZero_d.argtypes = [c_void_p,iType,iType]
lib.ElDistSparseMatrixQueueZero_d.restype = c_uint
lib.ElDistSparseMatrixQueueZero_c.argtypes = [c_void_p,iType,iType]
lib.ElDistSparseMatrixQueueZero_c.restype = c_uint
lib.ElDistSparseMatrixQueueZero_z.argtypes = [c_void_p,iType,iType]
lib.ElDistSparseMatrixQueueZero_z.restype = c_uint

lib.ElDistSparseMatrixQueueLocalZero_i.argtypes = [c_void_p,iType,iType]
lib.ElDistSparseMatrixQueueLocalZero_i.restype = c_uint
lib.ElDistSparseMatrixQueueLocalZero_s.argtypes = [c_void_p,iType,iType]
lib.ElDistSparseMatrixQueueLocalZero_s.restype = c_uint
lib.ElDistSparseMatrixQueueLocalZero_d.argtypes = [c_void_p,iType,iType]
lib.ElDistSparseMatrixQueueLocalZero_d.restype = c_uint
lib.ElDistSparseMatrixQueueLocalZero_c.argtypes = [c_void_p,iType,iType]
lib.ElDistSparseMatrixQueueLocalZero_c.restype = c_uint
lib.ElDistSparseMatrixQueueLocalZero_z.argtypes = [c_void_p,iType,iType]
lib.ElDistSparseMatrixQueueLocalZero_z.restype = c_uint

lib.ElDistSparseMatrixMakeConsistent_i.argtypes = [c_void_p]
lib.ElDistSparseMatrixMakeConsistent_i.restype = c_uint
lib.ElDistSparseMatrixMakeConsistent_s.argtypes = [c_void_p]
lib.ElDistSparseMatrixMakeConsistent_s.restype = c_uint
lib.ElDistSparseMatrixMakeConsistent_d.argtypes = [c_void_p]
lib.ElDistSparseMatrixMakeConsistent_d.restype = c_uint
lib.ElDistSparseMatrixMakeConsistent_c.argtypes = [c_void_p]
lib.ElDistSparseMatrixMakeConsistent_c.restype = c_uint
lib.ElDistSparseMatrixMakeConsistent_z.argtypes = [c_void_p]
lib.ElDistSparseMatrixMakeConsistent_z.restype = c_uint

lib.ElDistSparseMatrixHeight_i.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistSparseMatrixHeight_i.restype = c_uint
lib.ElDistSparseMatrixHeight_s.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistSparseMatrixHeight_s.restype = c_uint
lib.ElDistSparseMatrixHeight_d.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistSparseMatrixHeight_d.restype = c_uint
lib.ElDistSparseMatrixHeight_c.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistSparseMatrixHeight_c.restype = c_uint
lib.ElDistSparseMatrixHeight_z.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistSparseMatrixHeight_z.restype = c_uint

lib.ElDistSparseMatrixWidth_i.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistSparseMatrixWidth_i.restype = c_uint
lib.ElDistSparseMatrixWidth_s.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistSparseMatrixWidth_s.restype = c_uint
lib.ElDistSparseMatrixWidth_d.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistSparseMatrixWidth_d.restype = c_uint
lib.ElDistSparseMatrixWidth_c.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistSparseMatrixWidth_c.restype = c_uint
lib.ElDistSparseMatrixWidth_z.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistSparseMatrixWidth_z.restype = c_uint

lib.ElDistSparseMatrixDistGraph_i.argtypes = [c_void_p,POINTER(c_void_p)]
lib.ElDistSparseMatrixDistGraph_i.restype = c_uint
lib.ElDistSparseMatrixDistGraph_s.argtypes = [c_void_p,POINTER(c_void_p)]
lib.ElDistSparseMatrixDistGraph_s.restype = c_uint
lib.ElDistSparseMatrixDistGraph_d.argtypes = [c_void_p,POINTER(c_void_p)]
lib.ElDistSparseMatrixDistGraph_d.restype = c_uint
lib.ElDistSparseMatrixDistGraph_c.argtypes = [c_void_p,POINTER(c_void_p)]
lib.ElDistSparseMatrixDistGraph_c.restype = c_uint
lib.ElDistSparseMatrixDistGraph_z.argtypes = [c_void_p,POINTER(c_void_p)]
lib.ElDistSparseMatrixDistGraph_z.restype = c_uint

lib.ElDistSparseMatrixLockedDistGraph_i.argtypes = [c_void_p,POINTER(c_void_p)]
lib.ElDistSparseMatrixLockedDistGraph_i.restype = c_uint
lib.ElDistSparseMatrixLockedDistGraph_s.argtypes = [c_void_p,POINTER(c_void_p)]
lib.ElDistSparseMatrixLockedDistGraph_s.restype = c_uint
lib.ElDistSparseMatrixLockedDistGraph_d.argtypes = [c_void_p,POINTER(c_void_p)]
lib.ElDistSparseMatrixLockedDistGraph_d.restype = c_uint
lib.ElDistSparseMatrixLockedDistGraph_c.argtypes = [c_void_p,POINTER(c_void_p)]
lib.ElDistSparseMatrixLockedDistGraph_c.restype = c_uint
lib.ElDistSparseMatrixLockedDistGraph_z.argtypes = [c_void_p,POINTER(c_void_p)]
lib.ElDistSparseMatrixLockedDistGraph_z.restype = c_uint

lib.ElDistSparseMatrixFirstLocalRow_i.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistSparseMatrixFirstLocalRow_i.restype = c_uint
lib.ElDistSparseMatrixFirstLocalRow_s.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistSparseMatrixFirstLocalRow_s.restype = c_uint
lib.ElDistSparseMatrixFirstLocalRow_d.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistSparseMatrixFirstLocalRow_d.restype = c_uint
lib.ElDistSparseMatrixFirstLocalRow_c.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistSparseMatrixFirstLocalRow_c.restype = c_uint
lib.ElDistSparseMatrixFirstLocalRow_z.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistSparseMatrixFirstLocalRow_z.restype = c_uint

lib.ElDistSparseMatrixLocalHeight_i.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistSparseMatrixLocalHeight_i.restype = c_uint
lib.ElDistSparseMatrixLocalHeight_s.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistSparseMatrixLocalHeight_s.restype = c_uint
lib.ElDistSparseMatrixLocalHeight_d.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistSparseMatrixLocalHeight_d.restype = c_uint
lib.ElDistSparseMatrixLocalHeight_c.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistSparseMatrixLocalHeight_c.restype = c_uint
lib.ElDistSparseMatrixLocalHeight_z.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistSparseMatrixLocalHeight_z.restype = c_uint

lib.ElDistSparseMatrixNumLocalEntries_i.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistSparseMatrixNumLocalEntries_i.restype = c_uint
lib.ElDistSparseMatrixNumLocalEntries_s.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistSparseMatrixNumLocalEntries_s.restype = c_uint
lib.ElDistSparseMatrixNumLocalEntries_d.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistSparseMatrixNumLocalEntries_d.restype = c_uint
lib.ElDistSparseMatrixNumLocalEntries_c.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistSparseMatrixNumLocalEntries_c.restype = c_uint
lib.ElDistSparseMatrixNumLocalEntries_z.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistSparseMatrixNumLocalEntries_z.restype = c_uint

lib.ElDistSparseMatrixCapacity_i.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistSparseMatrixCapacity_i.restype = c_uint
lib.ElDistSparseMatrixCapacity_s.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistSparseMatrixCapacity_s.restype = c_uint
lib.ElDistSparseMatrixCapacity_d.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistSparseMatrixCapacity_d.restype = c_uint
lib.ElDistSparseMatrixCapacity_c.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistSparseMatrixCapacity_c.restype = c_uint
lib.ElDistSparseMatrixCapacity_z.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistSparseMatrixCapacity_z.restype = c_uint

lib.ElDistSparseMatrixConsistent_i.argtypes = [c_void_p,POINTER(bType)]
lib.ElDistSparseMatrixConsistent_i.restype = c_uint
lib.ElDistSparseMatrixConsistent_s.argtypes = [c_void_p,POINTER(bType)]
lib.ElDistSparseMatrixConsistent_s.restype = c_uint
lib.ElDistSparseMatrixConsistent_d.argtypes = [c_void_p,POINTER(bType)]
lib.ElDistSparseMatrixConsistent_d.restype = c_uint
lib.ElDistSparseMatrixConsistent_c.argtypes = [c_void_p,POINTER(bType)]
lib.ElDistSparseMatrixConsistent_c.restype = c_uint
lib.ElDistSparseMatrixConsistent_z.argtypes = [c_void_p,POINTER(bType)]
lib.ElDistSparseMatrixConsistent_z.restype = c_uint

lib.ElDistSparseMatrixComm_i.argtypes = [c_void_p,POINTER(mpi.Comm)]
lib.ElDistSparseMatrixComm_i.restype = c_uint
lib.ElDistSparseMatrixComm_s.argtypes = [c_void_p,POINTER(mpi.Comm)]
lib.ElDistSparseMatrixComm_s.restype = c_uint
lib.ElDistSparseMatrixComm_d.argtypes = [c_void_p,POINTER(mpi.Comm)]
lib.ElDistSparseMatrixComm_d.restype = c_uint
lib.ElDistSparseMatrixComm_c.argtypes = [c_void_p,POINTER(mpi.Comm)]
lib.ElDistSparseMatrixComm_c.restype = c_uint
lib.ElDistSparseMatrixComm_z.argtypes = [c_void_p,POINTER(mpi.Comm)]
lib.ElDistSparseMatrixComm_z.restype = c_uint

lib.ElDistSparseMatrixBlocksize_i.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistSparseMatrixBlocksize_i.restype = c_uint
lib.ElDistSparseMatrixBlocksize_s.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistSparseMatrixBlocksize_s.restype = c_uint
lib.ElDistSparseMatrixBlocksize_d.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistSparseMatrixBlocksize_d.restype = c_uint
lib.ElDistSparseMatrixBlocksize_c.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistSparseMatrixBlocksize_c.restype = c_uint
lib.ElDistSparseMatrixBlocksize_z.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistSparseMatrixBlocksize_z.restype = c_uint

lib.ElDistSparseMatrixRowOwner_i.argtypes = [c_void_p,iType,POINTER(c_int)]
lib.ElDistSparseMatrixRowOwner_i.restype = c_uint
lib.ElDistSparseMatrixRowOwner_s.argtypes = [c_void_p,iType,POINTER(c_int)]
lib.ElDistSparseMatrixRowOwner_s.restype = c_uint
lib.ElDistSparseMatrixRowOwner_d.argtypes = [c_void_p,iType,POINTER(c_int)]
lib.ElDistSparseMatrixRowOwner_d.restype = c_uint
lib.ElDistSparseMatrixRowOwner_c.argtypes = [c_void_p,iType,POINTER(c_int)]
lib.ElDistSparseMatrixRowOwner_c.restype = c_uint
lib.ElDistSparseMatrixRowOwner_z.argtypes = [c_void_p,iType,POINTER(c_int)]
lib.ElDistSparseMatrixRowOwner_z.restype = c_uint

lib.ElDistSparseMatrixRow_i.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDistSparseMatrixRow_i.restype = c_uint
lib.ElDistSparseMatrixRow_s.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDistSparseMatrixRow_s.restype = c_uint
lib.ElDistSparseMatrixRow_d.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDistSparseMatrixRow_d.restype = c_uint
lib.ElDistSparseMatrixRow_c.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDistSparseMatrixRow_c.restype = c_uint
lib.ElDistSparseMatrixRow_z.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDistSparseMatrixRow_z.restype = c_uint

lib.ElDistSparseMatrixCol_i.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDistSparseMatrixCol_i.restype = c_uint
lib.ElDistSparseMatrixCol_s.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDistSparseMatrixCol_s.restype = c_uint
lib.ElDistSparseMatrixCol_d.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDistSparseMatrixCol_d.restype = c_uint
lib.ElDistSparseMatrixCol_c.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDistSparseMatrixCol_c.restype = c_uint
lib.ElDistSparseMatrixCol_z.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDistSparseMatrixCol_z.restype = c_uint

lib.ElDistSparseMatrixValue_i.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDistSparseMatrixValue_i.restype = c_uint
lib.ElDistSparseMatrixValue_s.argtypes = [c_void_p,iType,POINTER(sType)]
lib.ElDistSparseMatrixValue_s.restype = c_uint
lib.ElDistSparseMatrixValue_d.argtypes = [c_void_p,iType,POINTER(dType)]
lib.ElDistSparseMatrixValue_d.restype = c_uint
lib.ElDistSparseMatrixValue_c.argtypes = [c_void_p,iType,POINTER(cType)]
lib.ElDistSparseMatrixValue_c.restype = c_uint
lib.ElDistSparseMatrixValue_z.argtypes = [c_void_p,iType,POINTER(zType)]
lib.ElDistSparseMatrixValue_z.restype = c_uint

lib.ElDistSparseMatrixEntryOffset_i.argtypes = \
  [c_void_p,iType,POINTER(iType)]
lib.ElDistSparseMatrixEntryOffset_i.restype = c_uint
lib.ElDistSparseMatrixEntryOffset_s.argtypes = \
  [c_void_p,iType,POINTER(iType)]
lib.ElDistSparseMatrixEntryOffset_s.restype = c_uint
lib.ElDistSparseMatrixEntryOffset_d.argtypes = \
  [c_void_p,iType,POINTER(iType)]
lib.ElDistSparseMatrixEntryOffset_d.restype = c_uint
lib.ElDistSparseMatrixEntryOffset_c.argtypes = \
  [c_void_p,iType,POINTER(iType)]
lib.ElDistSparseMatrixEntryOffset_c.restype = c_uint
lib.ElDistSparseMatrixEntryOffset_z.argtypes = \
  [c_void_p,iType,POINTER(iType)]
lib.ElDistSparseMatrixEntryOffset_z.restype = c_uint

lib.ElDistSparseMatrixNumConnections_i.argtypes = \
  [c_void_p,iType,POINTER(iType)]
lib.ElDistSparseMatrixNumConnections_i.restype = c_uint
lib.ElDistSparseMatrixNumConnections_s.argtypes = \
  [c_void_p,iType,POINTER(iType)]
lib.ElDistSparseMatrixNumConnections_s.restype = c_uint
lib.ElDistSparseMatrixNumConnections_d.argtypes = \
  [c_void_p,iType,POINTER(iType)]
lib.ElDistSparseMatrixNumConnections_d.restype = c_uint
lib.ElDistSparseMatrixNumConnections_c.argtypes = \
  [c_void_p,iType,POINTER(iType)]
lib.ElDistSparseMatrixNumConnections_c.restype = c_uint
lib.ElDistSparseMatrixNumConnections_z.argtypes = \
  [c_void_p,iType,POINTER(iType)]
lib.ElDistSparseMatrixNumConnections_z.restype = c_uint

lib.ElDistSparseMatrixSourceBuffer_i.argtypes = \
  [c_void_p,POINTER(POINTER(iType))]
lib.ElDistSparseMatrixSourceBuffer_i.restype = c_uint
lib.ElDistSparseMatrixSourceBuffer_s.argtypes = \
  [c_void_p,POINTER(POINTER(iType))]
lib.ElDistSparseMatrixSourceBuffer_s.restype = c_uint
lib.ElDistSparseMatrixSourceBuffer_d.argtypes = \
  [c_void_p,POINTER(POINTER(iType))]
lib.ElDistSparseMatrixSourceBuffer_d.restype = c_uint
lib.ElDistSparseMatrixSourceBuffer_c.argtypes = \
  [c_void_p,POINTER(POINTER(iType))]
lib.ElDistSparseMatrixSourceBuffer_c.restype = c_uint
lib.ElDistSparseMatrixSourceBuffer_z.argtypes = \
  [c_void_p,POINTER(POINTER(iType))]
lib.ElDistSparseMatrixSourceBuffer_z.restype = c_uint

lib.ElDistSparseMatrixLockedSourceBuffer_i.argtypes = \
  [c_void_p,POINTER(POINTER(iType))]
lib.ElDistSparseMatrixLockedSourceBuffer_i.restype = c_uint
lib.ElDistSparseMatrixLockedSourceBuffer_s.argtypes = \
  [c_void_p,POINTER(POINTER(iType))]
lib.ElDistSparseMatrixLockedSourceBuffer_s.restype = c_uint
lib.ElDistSparseMatrixLockedSourceBuffer_d.argtypes = \
  [c_void_p,POINTER(POINTER(iType))]
lib.ElDistSparseMatrixLockedSourceBuffer_d.restype = c_uint
lib.ElDistSparseMatrixLockedSourceBuffer_c.argtypes = \
  [c_void_p,POINTER(POINTER(iType))]
lib.ElDistSparseMatrixLockedSourceBuffer_c.restype = c_uint
lib.ElDistSparseMatrixLockedSourceBuffer_z.argtypes = \
  [c_void_p,POINTER(POINTER(iType))]
lib.ElDistSparseMatrixLockedSourceBuffer_z.restype = c_uint

lib.ElDistSparseMatrixTargetBuffer_i.argtypes = \
  [c_void_p,POINTER(POINTER(iType))]
lib.ElDistSparseMatrixTargetBuffer_i.restype = c_uint
lib.ElDistSparseMatrixTargetBuffer_s.argtypes = \
  [c_void_p,POINTER(POINTER(iType))]
lib.ElDistSparseMatrixTargetBuffer_s.restype = c_uint
lib.ElDistSparseMatrixTargetBuffer_d.argtypes = \
  [c_void_p,POINTER(POINTER(iType))]
lib.ElDistSparseMatrixTargetBuffer_d.restype = c_uint
lib.ElDistSparseMatrixTargetBuffer_c.argtypes = \
  [c_void_p,POINTER(POINTER(iType))]
lib.ElDistSparseMatrixTargetBuffer_c.restype = c_uint
lib.ElDistSparseMatrixTargetBuffer_z.argtypes = \
  [c_void_p,POINTER(POINTER(iType))]
lib.ElDistSparseMatrixTargetBuffer_z.restype = c_uint

lib.ElDistSparseMatrixLockedTargetBuffer_i.argtypes = \
  [c_void_p,POINTER(POINTER(iType))]
lib.ElDistSparseMatrixLockedTargetBuffer_i.restype = c_uint
lib.ElDistSparseMatrixLockedTargetBuffer_s.argtypes = \
  [c_void_p,POINTER(POINTER(iType))]
lib.ElDistSparseMatrixLockedTargetBuffer_s.restype = c_uint
lib.ElDistSparseMatrixLockedTargetBuffer_d.argtypes = \
  [c_void_p,POINTER(POINTER(iType))]
lib.ElDistSparseMatrixLockedTargetBuffer_d.restype = c_uint
lib.ElDistSparseMatrixLockedTargetBuffer_c.argtypes = \
  [c_void_p,POINTER(POINTER(iType))]
lib.ElDistSparseMatrixLockedTargetBuffer_c.restype = c_uint
lib.ElDistSparseMatrixLockedTargetBuffer_z.argtypes = \
  [c_void_p,POINTER(POINTER(iType))]
lib.ElDistSparseMatrixLockedTargetBuffer_z.restype = c_uint

lib.ElDistSparseMatrixValueBuffer_i.argtypes = \
  [c_void_p,POINTER(POINTER(iType))]
lib.ElDistSparseMatrixValueBuffer_i.restype = c_uint
lib.ElDistSparseMatrixValueBuffer_s.argtypes = \
  [c_void_p,POINTER(POINTER(sType))]
lib.ElDistSparseMatrixValueBuffer_s.restype = c_uint
lib.ElDistSparseMatrixValueBuffer_d.argtypes = \
  [c_void_p,POINTER(POINTER(dType))]
lib.ElDistSparseMatrixValueBuffer_d.restype = c_uint
lib.ElDistSparseMatrixValueBuffer_c.argtypes = \
  [c_void_p,POINTER(POINTER(cType))]
lib.ElDistSparseMatrixValueBuffer_c.restype = c_uint
lib.ElDistSparseMatrixValueBuffer_z.argtypes = \
  [c_void_p,POINTER(POINTER(zType))]
lib.ElDistSparseMatrixValueBuffer_z.restype = c_uint

lib.ElDistSparseMatrixLockedValueBuffer_i.argtypes = \
  [c_void_p,POINTER(POINTER(iType))]
lib.ElDistSparseMatrixLockedValueBuffer_i.restype = c_uint
lib.ElDistSparseMatrixLockedValueBuffer_s.argtypes = \
  [c_void_p,POINTER(POINTER(sType))]
lib.ElDistSparseMatrixLockedValueBuffer_s.restype = c_uint
lib.ElDistSparseMatrixLockedValueBuffer_d.argtypes = \
  [c_void_p,POINTER(POINTER(dType))]
lib.ElDistSparseMatrixLockedValueBuffer_d.restype = c_uint
lib.ElDistSparseMatrixLockedValueBuffer_c.argtypes = \
  [c_void_p,POINTER(POINTER(cType))]
lib.ElDistSparseMatrixLockedValueBuffer_c.restype = c_uint
lib.ElDistSparseMatrixLockedValueBuffer_z.argtypes = \
  [c_void_p,POINTER(POINTER(zType))]
lib.ElDistSparseMatrixLockedValueBuffer_z.restype = c_uint

class DistSparseMatrix(object):
  # Constructors and destructors
  # ============================
  def __init__(self,tag=dTag,comm=mpi.COMM_WORLD(),create=True):
    self.obj = c_void_p()
    self.tag = tag
    CheckTag(tag)
    if create:
      args = [pointer(self.obj),comm]
      if   tag == iTag: lib.ElDistSparseMatrixCreate_i(*args)
      elif tag == sTag: lib.ElDistSparseMatrixCreate_s(*args)
      elif tag == dTag: lib.ElDistSparseMatrixCreate_d(*args)
      elif tag == cTag: lib.ElDistSparseMatrixCreate_c(*args)
      elif tag == zTag: lib.ElDistSparseMatrixCreate_z(*args)
      else: DataExcept()
  def Destroy(self):
    args = [self.obj]
    if   self.tag == iTag: lib.ElDistSparseMatrixDestroy_i(*args)
    elif self.tag == sTag: lib.ElDistSparseMatrixDestroy_s(*args)
    elif self.tag == dTag: lib.ElDistSparseMatrixDestroy_d(*args)
    elif self.tag == cTag: lib.ElDistSparseMatrixDestroy_c(*args)
    elif self.tag == zTag: lib.ElDistSparseMatrixDestroy_z(*args)
    else: DataExcept()
  # Assignment and reconfiguration
  # ==============================
  def Empty(self):
    args = [self.obj]
    if   self.tag == iTag: lib.ElDistSparseMatrixEmpty_i(*args)
    elif self.tag == sTag: lib.ElDistSparseMatrixEmpty_s(*args)
    elif self.tag == dTag: lib.ElDistSparseMatrixEmpty_d(*args)
    elif self.tag == cTag: lib.ElDistSparseMatrixEmpty_c(*args)
    elif self.tag == zTag: lib.ElDistSparseMatrixEmpty_z(*args)
    else: DataExcept()
  def Resize(self,height,width):
    args = [self.obj,height,width]
    if   self.tag == iTag: lib.ElDistSparseMatrixResize_i(*args)
    elif self.tag == sTag: lib.ElDistSparseMatrixResize_s(*args)
    elif self.tag == dTag: lib.ElDistSparseMatrixResize_d(*args)
    elif self.tag == cTag: lib.ElDistSparseMatrixResize_c(*args)
    elif self.tag == zTag: lib.ElDistSparseMatrixResize_z(*args)
    else: DataExcept()
  def SetComm(self,comm):
    args = [self.obj,comm]
    if   self.tag == iTag: lib.ElDistSparseMatrixSetComm_i(*args)
    elif self.tag == sTag: lib.ElDistSparseMatrixSetComm_s(*args)
    elif self.tag == dTag: lib.ElDistSparseMatrixSetComm_d(*args)
    elif self.tag == cTag: lib.ElDistSparseMatrixSetComm_c(*args)
    elif self.tag == zTag: lib.ElDistSparseMatrixSetComm_z(*args)
    else: DataExcept()
  def Reserve(self,numEntries):
    args = [self.obj,numEntries]
    if   self.tag == iTag: lib.ElDistSparseMatrixReserve_i(*args)
    elif self.tag == sTag: lib.ElDistSparseMatrixReserve_s(*args)
    elif self.tag == dTag: lib.ElDistSparseMatrixReserve_d(*args)
    elif self.tag == cTag: lib.ElDistSparseMatrixReserve_c(*args)
    elif self.tag == zTag: lib.ElDistSparseMatrixReserve_z(*args)
    else: DataExcept()
  def Update(self,row,col,value):
    args = [self.obj,row,col,value]
    if   self.tag == iTag: lib.ElDistSparseMatrixUpdate_i(*args)
    elif self.tag == sTag: lib.ElDistSparseMatrixUpdate_s(*args)
    elif self.tag == dTag: lib.ElDistSparseMatrixUpdate_d(*args)
    elif self.tag == cTag: lib.ElDistSparseMatrixUpdate_c(*args)
    elif self.tag == zTag: lib.ElDistSparseMatrixUpdate_z(*args)
    else: DataExcept()
  def UpdateLocal(self,localRow,col,value):
    args = [self.obj,localRow,col,value]
    if   self.tag == iTag: lib.ElDistSparseMatrixUpdateLocal_i(*args)
    elif self.tag == sTag: lib.ElDistSparseMatrixUpdateLocal_s(*args)
    elif self.tag == dTag: lib.ElDistSparseMatrixUpdateLocal_d(*args)
    elif self.tag == cTag: lib.ElDistSparseMatrixUpdateLocal_c(*args)
    elif self.tag == zTag: lib.ElDistSparseMatrixUpdateLocal_z(*args)
    else: DataExcept()
  def Zero(self,row,col):
    args = [self.obj,row,col]
    if   self.tag == iTag: lib.ElDistSparseMatrixZero_i(*args)
    elif self.tag == sTag: lib.ElDistSparseMatrixZero_s(*args)
    elif self.tag == dTag: lib.ElDistSparseMatrixZero_d(*args)
    elif self.tag == cTag: lib.ElDistSparseMatrixZero_c(*args)
    elif self.tag == zTag: lib.ElDistSparseMatrixZero_z(*args)
    else: DataExcept()
  def ZeroLocal(self,localRow,col):
    args = [self.obj,localRow,col]
    if   self.tag == iTag: lib.ElDistSparseMatrixZeroLocal_i(*args)
    elif self.tag == sTag: lib.ElDistSparseMatrixZeroLocal_s(*args)
    elif self.tag == dTag: lib.ElDistSparseMatrixZeroLocal_d(*args)
    elif self.tag == cTag: lib.ElDistSparseMatrixZeroLocal_c(*args)
    elif self.tag == zTag: lib.ElDistSparseMatrixZeroLocal_z(*args)
    else: DataExcept()
  def QueueUpdate(self,row,col,value):
    args = [self.obj,row,col,value]
    if   self.tag == iTag: lib.ElDistSparseMatrixQueueUpdate_i(*args)
    elif self.tag == sTag: lib.ElDistSparseMatrixQueueUpdate_s(*args)
    elif self.tag == dTag: lib.ElDistSparseMatrixQueueUpdate_d(*args)
    elif self.tag == cTag: lib.ElDistSparseMatrixQueueUpdate_c(*args)
    elif self.tag == zTag: lib.ElDistSparseMatrixQueueUpdate_z(*args)
    else: DataExcept()
  def QueueLocalUpdate(self,localRow,col,value):
    args = [self.obj,localRow,col,value]
    if   self.tag == iTag: lib.ElDistSparseMatrixQueueLocalUpdate_i(*args)
    elif self.tag == sTag: lib.ElDistSparseMatrixQueueLocalUpdate_s(*args)
    elif self.tag == dTag: lib.ElDistSparseMatrixQueueLocalUpdate_d(*args)
    elif self.tag == cTag: lib.ElDistSparseMatrixQueueLocalUpdate_c(*args)
    elif self.tag == zTag: lib.ElDistSparseMatrixQueueLocalUpdate_z(*args)
    else: DataExcept()
  def QueueZero(self,row,col):
    args = [self.obj,row,col]
    if   self.tag == iTag: lib.ElDistSparseMatrixQueueZero_i(*args)
    elif self.tag == sTag: lib.ElDistSparseMatrixQueueZero_s(*args)
    elif self.tag == dTag: lib.ElDistSparseMatrixQueueZero_d(*args)
    elif self.tag == cTag: lib.ElDistSparseMatrixQueueZero_c(*args)
    elif self.tag == zTag: lib.ElDistSparseMatrixQueueZero_z(*args)
    else: DataExcept()
  def QueueLocalZero(self,localRow,col):
    args = [self.obj,localRow,col]
    if   self.tag == iTag: lib.ElDistSparseMatrixQueueLocalZero_i(*args)
    elif self.tag == sTag: lib.ElDistSparseMatrixQueueLocalZero_s(*args)
    elif self.tag == dTag: lib.ElDistSparseMatrixQueueLocalZero_d(*args)
    elif self.tag == cTag: lib.ElDistSparseMatrixQueueLocalZero_c(*args)
    elif self.tag == zTag: lib.ElDistSparseMatrixQueueLocalZero_z(*args)
    else: DataExcept()
  def MakeConsistent(self):
    args = [self.obj]
    if   self.tag == iTag: lib.ElDistSparseMatrixMakeConsistent_i(*args)
    elif self.tag == sTag: lib.ElDistSparseMatrixMakeConsistent_s(*args)
    elif self.tag == dTag: lib.ElDistSparseMatrixMakeConsistent_d(*args)
    elif self.tag == cTag: lib.ElDistSparseMatrixMakeConsistent_c(*args)
    elif self.tag == zTag: lib.ElDistSparseMatrixMakeConsistent_z(*args)
    else: DataExcept()
  # Queries
  # =======
  def Height(self):
    height = iType()
    args = [self.obj,pointer(height)]
    if   self.tag == iTag: lib.ElDistSparseMatrixHeight_i(*args)
    elif self.tag == sTag: lib.ElDistSparseMatrixHeight_s(*args)
    elif self.tag == dTag: lib.ElDistSparseMatrixHeight_d(*args)
    elif self.tag == cTag: lib.ElDistSparseMatrixHeight_c(*args)
    elif self.tag == zTag: lib.ElDistSparseMatrixHeight_z(*args)
    else: DataExcept()
    return height.value
  def Width(self):
    width = iType()
    args = [self.obj,pointer(width)]
    if   self.tag == iTag: lib.ElDistSparseMatrixWidth_i(*args)
    elif self.tag == sTag: lib.ElDistSparseMatrixWidth_s(*args)
    elif self.tag == dTag: lib.ElDistSparseMatrixWidth_d(*args)
    elif self.tag == cTag: lib.ElDistSparseMatrixWidth_c(*args)
    elif self.tag == zTag: lib.ElDistSparseMatrixWidth_z(*args)
    else: DataExcept()
    return width.value
  def DistGraph(self):
    graph = DG.DistGraph(mpi.COMM_WORLD(),False)
    args = [self.obj,pointer(graph.obj)]
    if   self.tag == iTag: lib.ElDistSparseMatrixDistGraph_i(*args)  
    elif self.tag == sTag: lib.ElDistSparseMatrixDistGraph_s(*args)
    elif self.tag == dTag: lib.ElDistSparseMatrixDistGraph_d(*args)
    elif self.tag == cTag: lib.ElDistSparseMatrixDistGraph_c(*args)
    elif self.tag == zTag: lib.ElDistSparseMatrixDistGraph_z(*args)
    else: DataExcept()
    return graph
  def LockedDistGraph(self):
    graph = DG.DistGraph(mpi.COMM_WORLD(),False)
    args = [self.obj,pointer(graph.obj)]
    if   self.tag == iTag: lib.ElDistSparseMatrixLockedDistGraph_i(*args)  
    elif self.tag == sTag: lib.ElDistSparseMatrixLockedDistGraph_s(*args)
    elif self.tag == dTag: lib.ElDistSparseMatrixLockedDistGraph_d(*args)
    elif self.tag == cTag: lib.ElDistSparseMatrixLockedDistGraph_c(*args)
    elif self.tag == zTag: lib.ElDistSparseMatrixLockedDistGraph_z(*args)
    else: DataExcept()
    return graph
  def FirstLocalRow(self):
    firstLocalRow = iType()
    args = [self.obj,pointer(firstLocalRow)]
    if   self.tag == iTag: lib.ElDistSparseMatrixFirstLocalRow_i(*args)
    elif self.tag == sTag: lib.ElDistSparseMatrixFirstLocalRow_s(*args)
    elif self.tag == dTag: lib.ElDistSparseMatrixFirstLocalRow_d(*args)
    elif self.tag == cTag: lib.ElDistSparseMatrixFirstLocalRow_c(*args)
    elif self.tag == zTag: lib.ElDistSparseMatrixFirstLocalRow_z(*args)
    else: DataExcept()
    return firstLocalRow.value
  def LocalHeight(self):
    localHeight = iType()
    args = [self.obj,pointer(localHeight)]
    if   self.tag == iTag: lib.ElDistSparseMatrixLocalHeight_i(*args)
    elif self.tag == sTag: lib.ElDistSparseMatrixLocalHeight_s(*args)
    elif self.tag == dTag: lib.ElDistSparseMatrixLocalHeight_d(*args)
    elif self.tag == cTag: lib.ElDistSparseMatrixLocalHeight_c(*args)
    elif self.tag == zTag: lib.ElDistSparseMatrixLocalHeight_z(*args)
    else: DataExcept()
    return localHeight.value
  def NumLocalEntries(self):
    numLocalEntries = iType()
    args = [self.obj,pointer(numLocalEntries)]
    if   self.tag == iTag: lib.ElDistSparseMatrixNumLocalEntries_i(*args)
    elif self.tag == sTag: lib.ElDistSparseMatrixNumLocalEntries_s(*args)
    elif self.tag == dTag: lib.ElDistSparseMatrixNumLocalEntries_d(*args)
    elif self.tag == cTag: lib.ElDistSparseMatrixNumLocalEntries_c(*args)
    elif self.tag == zTag: lib.ElDistSparseMatrixNumLocalEntries_z(*args)
    else: DataExcept()
    return numLocalEntries.value
  def Capacity(self):
    capacity = iType()
    args = [self.obj,pointer(capacity)]
    if   self.tag == iTag: lib.ElDistSparseMatrixCapacity_i(*args)
    elif self.tag == sTag: lib.ElDistSparseMatrixCapacity_s(*args)
    elif self.tag == dTag: lib.ElDistSparseMatrixCapacity_d(*args)
    elif self.tag == cTag: lib.ElDistSparseMatrixCapacity_c(*args)
    elif self.tag == zTag: lib.ElDistSparseMatrixCapacity_z(*args)
    else: DataExcept()
    return capacity.value
  def Consistent(self):
    consistent = bType()
    args = [self.obj,pointer(consistent)]
    if   self.tag == iTag: lib.ElDistSparseMatrixConsistent_i(*args)
    elif self.tag == sTag: lib.ElDistSparseMatrixConsistent_s(*args)
    elif self.tag == dTag: lib.ElDistSparseMatrixConsistent_d(*args)
    elif self.tag == cTag: lib.ElDistSparseMatrixConsistent_c(*args)
    elif self.tag == zTag: lib.ElDistSparseMatrixConsistent_z(*args)
    else: DataExcept()
    return consistent.value
  def Comm(self):
    comm = mpi.Comm()
    args = [self.obj,pointer(comm)]
    if   self.tag == iTag: lib.ElDistSparseMatrixComm_i(*args)
    elif self.tag == sTag: lib.ElDistSparseMatrixComm_s(*args)
    elif self.tag == dTag: lib.ElDistSparseMatrixComm_d(*args)
    elif self.tag == cTag: lib.ElDistSparseMatrixComm_c(*args)
    elif self.tag == zTag: lib.ElDistSparseMatrixComm_z(*args)
    else: DataExcept()
    return comm
  def Blocksize(self):
    blocksize = iType()
    args = [self.obj,pointer(blocksize)]
    if   self.tag == iTag: lib.ElDistSparseMatrixBlocksize_i(*args)
    elif self.tag == sTag: lib.ElDistSparseMatrixBlocksize_s(*args)
    elif self.tag == dTag: lib.ElDistSparseMatrixBlocksize_d(*args)
    elif self.tag == cTag: lib.ElDistSparseMatrixBlocksize_c(*args)
    elif self.tag == zTag: lib.ElDistSparseMatrixBlocksize_z(*args)
    else: DataExcept()
    return blocksize.value
  def RowOwner(self,i):
    owner = c_int()
    args = [self.obj,i,pointer(owner)]
    if   self.tag == iTag: lib.ElDistSparseMatrixRowOwner_i(*args)
    elif self.tag == sTag: lib.ElDistSparseMatrixRowOwner_s(*args)
    elif self.tag == dTag: lib.ElDistSparseMatrixRowOwner_d(*args)
    elif self.tag == cTag: lib.ElDistSparseMatrixRowOwner_c(*args)
    elif self.tag == zTag: lib.ElDistSparseMatrixRowOwner_z(*args)
    else: DataExcept()
    return owner.value
  def Row(self,localInd):
    row = iType() 
    args = [self.obj,localInd,pointer(row)]
    if   self.tag == iTag: lib.ElDistSparseMatrixRow_i(*args)
    elif self.tag == sTag: lib.ElDistSparseMatrixRow_s(*args)
    elif self.tag == dTag: lib.ElDistSparseMatrixRow_d(*args)
    elif self.tag == cTag: lib.ElDistSparseMatrixRow_c(*args)
    elif self.tag == zTag: lib.ElDistSparseMatrixRow_z(*args)
    else: DataExcept()
    return row.value
  def Col(self,localInd):
    col = iType() 
    args = [self.obj,localInd,pointer(col)]
    if   self.tag == iTag: lib.ElDistSparseMatrixCol_i(*args)
    elif self.tag == sTag: lib.ElDistSparseMatrixCol_s(*args)
    elif self.tag == dTag: lib.ElDistSparseMatrixCol_d(*args)
    elif self.tag == cTag: lib.ElDistSparseMatrixCol_c(*args)
    elif self.tag == zTag: lib.ElDistSparseMatrixCol_z(*args)
    else: DataExcept()
    return col.value
  def Value(self,localInd):
    value =  TagToType(self.tag)()
    args = [self.obj,localInd,pointer(value)]
    if   self.tag == iTag: lib.ElDistSparseMatrixValue_i(*args)
    elif self.tag == sTag: lib.ElDistSparseMatrixValue_s(*args)
    elif self.tag == dTag: lib.ElDistSparseMatrixValue_d(*args)
    elif self.tag == cTag: lib.ElDistSparseMatrixValue_c(*args)
    elif self.tag == zTag: lib.ElDistSparseMatrixValue_z(*args)
    else: DataExcept()
    return value.value
  def EntryOffset(self,localRow):
    localOffset = iType()
    args = [self.obj,localRow,pointer(localOffset)]
    if   self.tag == iTag: lib.ElDistSparseMatrixEntryOffset_i(*args)
    elif self.tag == sTag: lib.ElDistSparseMatrixEntryOffset_s(*args)
    elif self.tag == dTag: lib.ElDistSparseMatrixEntryOffset_d(*args)
    elif self.tag == cTag: lib.ElDistSparseMatrixEntryOffset_c(*args)
    elif self.tag == zTag: lib.ElDistSparseMatrixEntryOffset_z(*args)
    else: DataExcept()
    return localOffset.value
  def NumConnections(self,localRow):
    numConnections = iType()
    args = [self.obj,localRow,pointer(numConnections)]
    if   self.tag == iTag: lib.ElDistSparseMatrixNumConnections_i(*args)
    elif self.tag == sTag: lib.ElDistSparseMatrixNumConnections_s(*args)
    elif self.tag == dTag: lib.ElDistSparseMatrixNumConnections_d(*args)
    elif self.tag == cTag: lib.ElDistSparseMatrixNumConnections_c(*args)
    elif self.tag == zTag: lib.ElDistSparseMatrixNumConnections_z(*args)
    else: DataExcept()
    return numConnections.value
  def SourceBuffer(self):
    sourceBuf = POINTER(iType)()
    args = [self.obj,pointer(sourceBuf)]    
    if   self.tag == iTag: lib.ElDistSparseMatrixSourceBuffer_i(*args)
    elif self.tag == sTag: lib.ElDistSparseMatrixSourceBuffer_s(*args)
    elif self.tag == dTag: lib.ElDistSparseMatrixSourceBuffer_d(*args)
    elif self.tag == cTag: lib.ElDistSparseMatrixSourceBuffer_c(*args)
    elif self.tag == zTag: lib.ElDistSparseMatrixSourceBuffer_z(*args)
    else: DataExcept()
    return sourceBuf
  def LockedSourceBuffer(self):
    sourceBuf = POINTER(iType)()
    args = [self.obj,pointer(sourceBuf)]    
    if   self.tag == iTag: lib.ElDistSparseMatrixLockedSourceBuffer_i(*args)
    elif self.tag == sTag: lib.ElDistSparseMatrixLockedSourceBuffer_s(*args)
    elif self.tag == dTag: lib.ElDistSparseMatrixLockedSourceBuffer_d(*args)
    elif self.tag == cTag: lib.ElDistSparseMatrixLockedSourceBuffer_c(*args)
    elif self.tag == zTag: lib.ElDistSparseMatrixLockedSourceBuffer_z(*args)
    else: DataExcept()
    return sourceBuf
  def TargetBuffer(self):
    targetBuf = POINTER(iType)()
    args = [self.obj,pointer(targetBuf)]    
    if   self.tag == iTag: lib.ElDistSparseMatrixTargetBuffer_i(*args)
    elif self.tag == sTag: lib.ElDistSparseMatrixTargetBuffer_s(*args)
    elif self.tag == dTag: lib.ElDistSparseMatrixTargetBuffer_d(*args)
    elif self.tag == cTag: lib.ElDistSparseMatrixTargetBuffer_c(*args)
    elif self.tag == zTag: lib.ElDistSparseMatrixTargetBuffer_z(*args)
    else: DataExcept()
    return targetBuf
  def LockedTargetBuffer(self):
    targetBuf = POINTER(iType)()
    args = [self.obj,pointer(targetBuf)]    
    if   self.tag == iTag: lib.ElDistSparseMatrixLockedTargetBuffer_i(*args)
    elif self.tag == sTag: lib.ElDistSparseMatrixLockedTargetBuffer_s(*args)
    elif self.tag == dTag: lib.ElDistSparseMatrixLockedTargetBuffer_d(*args)
    elif self.tag == cTag: lib.ElDistSparseMatrixLockedTargetBuffer_c(*args)
    elif self.tag == zTag: lib.ElDistSparseMatrixLockedTargetBuffer_z(*args)
    else: DataExcept()
    return targetBuf
  def ValueBuffer(self):
    valueBuf = POINTER(TagToType(self.tag))()
    args = [self.obj,pointer(valueBuf)]
    if   self.tag == iTag: lib.ElDistSparseMatrixValueBuffer_i(*args)
    elif self.tag == sTag: lib.ElDistSparseMatrixValueBuffer_s(*args)
    elif self.tag == dTag: lib.ElDistSparseMatrixValueBuffer_d(*args)
    elif self.tag == cTag: lib.ElDistSparseMatrixValueBuffer_c(*args)
    elif self.tag == zTag: lib.ElDistSparseMatrixValueBuffer_z(*args)
    else: DataExcept()
    return valueBuf
  def LockedValueBuffer(self):
    valueBuf = POINTER(TagToType(self.tag))()
    args = [self.obj,pointer(valueBuf)]
    if   self.tag == iTag: lib.ElDistSparseMatrixLockedValueBuffer_i(*args)
    elif self.tag == sTag: lib.ElDistSparseMatrixLockedValueBuffer_s(*args)
    elif self.tag == dTag: lib.ElDistSparseMatrixLockedValueBuffer_d(*args)
    elif self.tag == cTag: lib.ElDistSparseMatrixLockedValueBuffer_c(*args)
    elif self.tag == zTag: lib.ElDistSparseMatrixLockedValueBuffer_z(*args)
    else: DataExcept()
    return valueBuf
