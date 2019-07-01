#include <ross.h>

#include "simtime.hh"
#include "matlpcancel.hh"

matlp::global_stat_t matlp::global_stat;

/* Global variables */
static unsigned int offset_lpid = 0;
static unsigned int ttl_lps = 0;
static unsigned int nlp_per_pe = 8;

// rate for timestamp exponential distribution
static double mean = 1.0;

//static char run_id[1024] = "undefined";

tw_peid ctr_map(tw_lpid gid) {
  return (tw_peid) gid / g_tw_nlp;
}

void do_nothing(void *a,void *b,void *c) {}

tw_lptype mylps[] = {
  { (init_f) matlp::init,
    (pre_run_f) do_nothing,
    (event_f) matlp::forward_event,
    (revent_f) matlp::backward_event,
    (commit_f) matlp::commit_event,
    (final_f) matlp::finalize,
    (map_f) ctr_map,
    //NULL,
    sizeof(matlp)
  },
  {0}
};

const tw_optdef app_opt[] =
  {
    TWOPT_GROUP("Particle spammer"),
    TWOPT_UINT("nlp", nlp_per_pe, "number of LPs per processor"),
    TWOPT_DOUBLE("mean", mean, "exponential distribution mean for timestamps"),
    TWOPT_END()
  };

int
main(int argc, char **argv, char **env)
{
  int i;

  //printf("Process starting main...\n");

  // get rid of error if compiled w/ MEMORY queues
  g_tw_memory_nqueues=1;

  //printf("Calling tw_opt_add\n");
  tw_opt_add(app_opt);

  //printf("Calling tw_init\n");
  tw_init(&argc, &argv);
  //printf("tw_init returned.\n");

  g_tw_memory_nqueues = 16; // give at least 16 memory queue event

  offset_lpid = g_tw_mynode * nlp_per_pe;
  ttl_lps = tw_nnodes() * g_tw_npe * nlp_per_pe;
  g_tw_events_per_pe = nlp_per_pe * 10;

  /* Initialize statistics data... */ {
    matlp::global_stat.nforward = 0;
    matlp::global_stat.nbackward = 0;
    matlp::global_stat.ncommit = 0;

    matlp::global_stat.ncancelforward = 0;
    matlp::global_stat.ncancelbackward = 0;
    matlp::global_stat.ncancelcommit = 0;

    matlp::global_stat.nretract = 0;

    for(int i = 0; i<nmat; i++)
      for(int j = 0; j<nmat; j++)
	matlp::global_stat.A[i*nmat+j] = (i == j);
  }

  tw_define_lps(nlp_per_pe, sizeof(Event), 0);

  if(0) {
    for(int i = 0; i<(int) tw_nnodes(); i++) {
      MPI_Barrier(MPI_COMM_WORLD);
      if(g_tw_mynode == i) {
	printf("I am proc. %d, g_tw_lpoffset = %d\n",
	       i,(int) g_tw_lp_offset);
	fflush(stdout);
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  for(i = 0; i < (int) g_tw_nlp; i++)
    tw_lp_settype(i, &mylps[0]);

  if(0) {
    for(int i = 0; i<(int) tw_nnodes(); i++) {
      MPI_Barrier(MPI_COMM_WORLD);
      if(g_tw_mynode == i) {
	printf("I am proc. %d, g_tw_lpoffset = %d\n",
	       i,(int) g_tw_lp_offset);
	fflush(stdout);
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  if(g_tw_mynode == 0) {
      printf("========================================\n");
      printf("Particle Spammer Configuration\n");
      printf("   Mean...................%lf\n", mean);
      printf("   nlp_per_pe = %lld\n", (long long int) nlp_per_pe);
      printf("   g_tw_nlp = %lld\n", (long long int) g_tw_nlp);
      printf("========================================\n\n");
    }


  MPI_Barrier(MPI_COMM_WORLD);
  if(tw_ismaster()) printf("@@@ Calling run...\n");
  MPI_Barrier(MPI_COMM_WORLD);

  tw_run();

  MPI_Barrier(MPI_COMM_WORLD);
  if(tw_ismaster()) printf("@@@ Run returned...\n");
  MPI_Barrier(MPI_COMM_WORLD);


  if(tw_ismaster()) {
    
  }


#if 1
  {
    int pid,np;
    MPI_Comm_rank(MPI_COMM_WORLD,&pid);
    MPI_Comm_size(MPI_COMM_WORLD,&np);

    MPI_Barrier(MPI_COMM_WORLD);

    for(int bit = 1; bit<np; bit<<=1) {
      if((pid & (bit-1)) == 0 && (pid | bit) < np) {
	if((pid & bit)) {
	  //printf("   %d sending to %d\n",pid,pid-bit);
	  MPI_Send(&matlp::global_stat,sizeof(matlp::global_stat),MPI_BYTE,pid-bit,142,MPI_COMM_WORLD);
	} else {
	  matlp::global_stat_t tmp,tmp2;
	  //printf("   %d receiving from %d\n",pid,pid+bit);
	  MPI_Recv(&tmp,sizeof(tmp),MPI_BYTE,pid+bit,142,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	  matlp::global_stat.nforward += tmp.nforward;
	  matlp::global_stat.nbackward += tmp.nbackward;
	  matlp::global_stat.ncommit += tmp.ncommit;

	  matlp::global_stat.ncancelforward += tmp.ncancelforward;
	  matlp::global_stat.ncancelbackward += tmp.ncancelbackward;
	  matlp::global_stat.ncancelcommit += tmp.ncancelcommit;

	  matlp::global_stat.nretract += tmp.nretract;

	  matmul(nmat,matlp::global_stat.A,tmp.A,tmp2.A);
	  memcpy(matlp::global_stat.A,tmp2.A,nmat*nmat*sizeof(*(tmp2.A)));

#if 0
	  printf("  A =\n");
	  for(int i = 0; i<nmat; i++) {
	    printf("    ");
	    for(int j = 0; j<nmat; j++) {
	      printf("  %12llu",(unsigned long long int) matlp::global_stat.A[i*nmat+j]);
	    }
	    printf("\n");
	  }
	  printf("\n");
#endif
	}
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if(pid == 0) {
      printf("FINAL STATE:\n");
      printf("  nforward  = %12lld  ncancelforward         = %12lld  nretract = %12lld\n",
	     matlp::global_stat.nforward,matlp::global_stat.ncancelforward,
	     matlp::global_stat.nretract);
      printf("  nbackward = %12lld  ncancelbackward        = %12lld\n",
	     matlp::global_stat.nbackward,matlp::global_stat.ncancelbackward);
      printf("  ncommit   = %12lld  ncancelcommit+nretract = %12lld\n",
	     matlp::global_stat.ncommit,
	     matlp::global_stat.ncancelcommit + matlp::global_stat.nretract);
      printf("  nf-nb-nc  = %12lld  ncf-ncb-ncc          = %12lld\n",
	     matlp::global_stat.nforward -
	     matlp::global_stat.nbackward -
	     matlp::global_stat.ncommit,
	     matlp::global_stat.ncancelforward -
	     matlp::global_stat.ncancelbackward -
	     matlp::global_stat.ncancelcommit
	     );
      printf("  A =\n");
      for(int i = 0; i<nmat; i++) {
	printf("    ");
	for(int j = 0; j<nmat; j++) {
	  printf("  %12llu",(unsigned long long int) matlp::global_stat.A[i*nmat+j]);
	}
	printf("\n");
      }
    }
  }
#endif

  tw_end();

  return 0;
}
