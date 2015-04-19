#include "El/config.h"

#ifdef EL_HAVE_PARMETIS
/*
 * Copyright 1997, Regents of the University of Minnesota
 * Modified by Jack Poulson, 2012
 */
#include <parmetislib.h>

void ElParallelLabelVertices
( ctrl_t *ctrl, graph_t *graph, idx_t *order, idx_t *sizes )
{ 
  idx_t i, j, nvtxs, id; 
  idx_t *where, *lpwgts, *gpwgts;
  idx_t sizescan[3];

  nvtxs  = graph->nvtxs;
  where  = graph->where;
  lpwgts = graph->lpwgts;
  gpwgts = graph->gpwgts;

  /* Compute the local sizes of the left side, right side, and separator */
  iset(3, 0, lpwgts);
  for (i=0; i<nvtxs; i++) 
      lpwgts[where[i]]++;

  /* Perform a Prefix scan of the separator size to determine the boundaries */
  gkMPI_Scan((void *)lpwgts, (void *)sizescan, 3, IDX_T, MPI_SUM, ctrl->comm);
  gkMPI_Allreduce
  ((void *)lpwgts, (void *)gpwgts, 3, IDX_T, MPI_SUM, ctrl->comm);

  /* Fill in the size of the partition */
  sizes[0] = gpwgts[0];
  sizes[1] = gpwgts[1];
  sizes[2] = gpwgts[2];

  for( i=2; i>=0; --i )
      for( j=i+1; j<3; ++j )
          sizescan[i] += gpwgts[j];
  for( i=0; i<3; i++ )
      sizescan[i] -= lpwgts[i];

  for( i=0; i<nvtxs; i++ ) 
  {
      id = where[i];
      PASSERT(ctrl, id <= 2);
      sizescan[id]++;
      PASSERT(ctrl, order[i] == -1);
      order[i] = graph->gnvtxs - sizescan[id];
  }
}

void ElParallelOrder
( ctrl_t *ctrl, graph_t *graph, idx_t *order, idx_t *sizes )
{
  idx_t i, nvtxs;

  nvtxs = graph->nvtxs;
  iset(nvtxs, -1, order);

  /* graph->where = ismalloc(nvtxs, 0, "ElOrder: graph->where"); */
  /* If we computed an initial partition with Global_Partition, then we 
     should run the following instead of the above ismalloc of graph->where*/
  iset(nvtxs, 0, graph->where); 
  gk_free((void **)&graph->match, 
          (void **)&graph->cmap, 
          (void **)&graph->rlens, 
          (void **)&graph->slens, 
          (void **)&graph->rcand, LTERM);

  Order_Partition_Multiple(ctrl, graph);

  ElParallelLabelVertices(ctrl, graph, order, sizes);
}

void ElParallelBisect
( idx_t *vtxdist, idx_t *xadj, idx_t *adjncy, 
  idx_t *p_nseps, idx_t *s_nseps, 
  real_t *ubfrac, idx_t *idbglvl, idx_t *order, idx_t *sizes, 
  MPI_Comm *comm )
{
  idx_t i, j, npes, npesNonzero, mype, mypeNonzero, dbglvl, status, haveData;
  ctrl_t *ctrl;
  graph_t *graph;
  MPI_Comm nonzeroComm, nullComm;
  size_t curmem;

  gkMPI_Comm_size(*comm, &npes);
  gkMPI_Comm_rank(*comm, &mype);

  if( vtxdist[npes] == 0 )
  {
      sizes[0] = 0;
      sizes[1] = 0;
      sizes[2] = 0;
      return;
  }

  haveData = ( vtxdist[mype+1]-vtxdist[mype] != 0 );
  if( haveData )
      gkMPI_Comm_split(*comm, 1, mype, &nonzeroComm);
  else
      gkMPI_Comm_split(*comm, MPI_UNDEFINED, 0, &nullComm);

  if( !haveData )
  {
      sizes[0] = sizes[1] = sizes[2] = 0;
      gkMPI_Allreduce(MPI_IN_PLACE, (void *)sizes, 3, IDX_T, MPI_SUM, *comm);
      return;
  }

  gkMPI_Comm_size(nonzeroComm, &npesNonzero);
  gkMPI_Comm_rank(nonzeroComm, &mypeNonzero);

  /* Compress the vtxdist data to make it match the new communicator */
  j=0;
  for( i=1; i<npes+1; ++i )
      if( vtxdist[i] != vtxdist[j] )
          vtxdist[++j] = vtxdist[i];

  status = METIS_OK;
  gk_malloc_init();
  curmem = gk_GetCurMemoryUsed();

  ctrl = SetupCtrl(PARMETIS_OP_KMETIS, NULL, 1, 2, NULL, NULL, nonzeroComm);

  dbglvl = (idbglvl == NULL ? 0 : *idbglvl);
  ctrl->dbglvl = dbglvl;

  graph = SetupGraph(ctrl, 1, vtxdist, xadj, NULL, NULL, adjncy, NULL, 0);
  AllocateWSpace(ctrl, 10*graph->nvtxs);

  /* Compute an initial partition: for some reason this improves the quality */
  ctrl->CoarsenTo = gk_min(vtxdist[npesNonzero]+1, 
                           200*gk_max(npesNonzero,ctrl->nparts));
  Global_Partition(ctrl, graph); 

  /* Compute an ordering */
  ctrl->optype    = PARMETIS_OP_OMETIS;
  ctrl->partType  = ORDER_PARTITION;
  ctrl->mtype     = PARMETIS_MTYPE_GLOBAL;
  ctrl->rtype     = PARMETIS_SRTYPE_2PHASE;
  ctrl->p_nseps   = (p_nseps  == NULL ? 1 : *p_nseps);
  ctrl->s_nseps   = (s_nseps  == NULL ? 1 : *s_nseps);
  ctrl->ubfrac    = (ubfrac == NULL ? ORDER_UNBALANCE_FRACTION : *ubfrac);
  ctrl->dbglvl    = dbglvl;
  ctrl->ipart     = ISEP_NODE;
  ctrl->CoarsenTo = gk_min(graph->gnvtxs-1,1500*npesNonzero); 
  ElParallelOrder(ctrl, graph, order, sizes);

  FreeInitialGraphAndRemap(graph);

  /* Pass the data to the early-exiting processes with an allreduce */
  if( mypeNonzero != 0 )
      sizes[0] = sizes[1] = sizes[2] = 0;
  gkMPI_Allreduce(MPI_IN_PLACE, (void*)sizes, 3, IDX_T, MPI_SUM, *comm);

  MPI_Comm_free( &nonzeroComm );

  goto DONE;

DONE:
  FreeCtrl(&ctrl);
  if (gk_GetCurMemoryUsed() - curmem > 0) {
    printf("ParMETIS appears to have a memory leak of %zdbytes. Report this.\n",
        (ssize_t)(gk_GetCurMemoryUsed() - curmem));
  }
  gk_malloc_cleanup(0);
}

#else

/* This is only here to avoid warnings */
void ElFoo( int* bar ) { *bar = 42; }

#endif /* ifdef EL_HAVE_PARMETIS */
