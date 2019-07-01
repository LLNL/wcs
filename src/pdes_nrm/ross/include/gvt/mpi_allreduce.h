#ifndef INC_gvt_mpi_allreduce_h
#define INC_gvt_mpi_allreduce_h

extern double gvt_print_interval;
extern double percent_complete;

static inline int 
tw_gvt_inprogress(tw_pe * pe)
{
	return pe->gvt_status;
}

static inline void 
gvt_print(tw_stime gvt)
{
  double ttreal;
  {
    static int first = 1;
    static struct timeval t0;
    struct timeval tt;
    
    gettimeofday(&tt,NULL);
    if(first) {
      t0 = tt;
      first = 0;
    }
    tt.tv_sec -= t0.tv_sec;
    tt.tv_usec -= t0.tv_usec;
    ttreal = tt.tv_sec + 1e-6*tt.tv_usec;
  }

	if(gvt_print_interval == 1.0)
		return;

	if(percent_complete == 0.0)
	{
		percent_complete = gvt_print_interval;
		return;
	}

	printf("GVT #%d: simulation %d%% complete, max event queue size %u (",
		g_tw_gvt_done,
	       (int) ROSS_MIN(100, floor(100 * (gvt.real()/g_tw_ts_end))), 
           tw_pq_max_size(g_tw_pe[0]->pq));

	if (gvt == DBL_MAX)
		printf("GVT = %s", "MAX");
	else
	  printf("GVT = %.4f", gvt.real());

	//printf(").\n");
	printf(").  Wall time = %.3fs\n",ttreal);
    
#ifdef AVL_TREE
    printf("AVL tree size: %d\n", g_tw_pe[0]->avl_tree_size);
#endif
    
	percent_complete += gvt_print_interval;
}

#endif
