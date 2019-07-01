#ifndef MATLP__
#define MATLP__

#include <cstring>
#include <list>

#include "zinv.hh"

static const int nmat = 5;
static const double mean_delay = 1.0;

template<typename uint>
void matprint(int n,uint A[],const char *label) {
  if(label)
    printf("  %s =\n",label);
  for(int i = 0; i<nmat; i++) {
    printf("    ");
    for(int j = 0; j<nmat; j++)
      printf("  %12llu",(unsigned long long int) A[i*nmat+j]);
    printf("\n");
  }
  printf("\n");
}


DECLARE_LP(Matrix) {

};

struct AllData {
   uint A[nmat*nmat];
   int cancel_count;
   int nretract,nforward,nbackward;

   void print(const char *label = 0) {
      matprint(nmat,A,label);
   }
};


DECLARE_EVENT(Ping, Matrix, AllData) {
   virtual void forward();
   virtual void backward();
   virtual void commit();
};

DECLARE_EVENT(Self, Matrix, AllData) {
   virtual void forward();
   virtual void backward();
   virtual void commit();
};



class LP {
 public:
   virtual Event* create(EventData* evt) = 0;
   virtual ~LP() {}
   static void forward_event(void *state, tw_bf *bf, void *evt, tw_lp *thislp) {
      Time tnow = tw_now(thislp);
      EventData *e = reinterpret_cast<EventData *>(evt);
      LP *lpptr = reinterpret_cast<LP *>(state);
      e->obj = lpptr->createEvent(e);
      e->obj->forward(tnow,bf,thislp);
   }

   static void backward_event(void *state, tw_bf *bf, void *evt, tw_lp *thislp) {
      Time tnow = tw_now(thislp);
      EventData *e = reinterpret_cast<EventData *>(evt);
      e->obj->backward(tnow,bf,thislp);
      delete e->obj;
   }
   static void commit_event(void *state, tw_bf *bf, void *evt, tw_lp *thislp) {
      Time tnow = tw_now(thislp);
      EventData *e = reinterpret_cast<EventData *>(evt);
      e->obj->commit(tnow,bf,thislp);
      delete e->obj;
   }
};




template<typename myuint>
bool invcheck(int n,myuint A[]) {
  myuint Atmp[n*n],Ai[n*n];
  memcpy(Atmp,A,n*n*sizeof(*A));
  return matinv(n,Atmp,Ai) == 1;
}

template<typename myuint>
void rand_matrix(tw_lp *lpptr,int n,myuint A[],int imax) {
  do {
    for(int i = 0; i<n; i++)
      for(int j = 0; j<n; j++) {
	myuint x = tw_rand_integer(lpptr->rng,0,imax);
	A[i*n+j] = x;
      }
  } while(!invcheck(n,A));
}


struct matlp : public LP {
  typedef unsigned int uint;
  uint *M,*T; 
  const int n;

  Event *next_event;

  struct global_stat_t {
    long long int nforward,nbackward,ncommit;
    long long int ncancelforward,ncancelbackward,ncancelcommit;
    long long int nretract;
    uint A[nmat*nmat];
  };
  static global_stat_t global_stat;

  long long int nforward,nbackward,ncommit;
  long long int ncancelforward,ncancelbackward,ncancelcommit;
  long long int nretract;

  void apply(uint X[]) {
    uint tmp[n*n];
    assert(invcheck(n,X));
    assert(invcheck(n,M));

    matmul(n,M,X,tmp);
    memcpy(M,tmp,n*n*sizeof(*M));
  }
  void unapply(uint X[]) {
    uint Xtmp[n*n],Xi[n*n];

    assert(invcheck(n,M));

    memcpy(Xtmp,X,n*n*sizeof(*Xtmp));
    assert(matinv(n,Xtmp,Xi) == 1);
    matmul(n,M,Xi,Xtmp);
    memcpy(M,Xtmp,n*n*sizeof(*M));

    assert(invcheck(n,M));
  }


  matlp(tw_lp *lpptr,int nin = 2) :
    n(nin),nforward(0),nbackward(0),ncommit(0),
    ncancelforward(0),ncancelbackward(0),ncancelcommit(0),
    nretract(0) {

    M = new uint[n*n];
    T = new uint[n*n];
    rand_matrix(lpptr,n,M,1000);
    rand_matrix(lpptr,n,T,1000);
    next_event = 0;

    if(0) {
      printf("Creating lp %d, matrix size %d x %d:\n",(int) lpptr->gid,n,n);
      matprint(n,M,"M");
      matprint(n,T,"T");
    }
    
    {
      tw_stime now = tw_now(lpptr);
      tw_stime trecv;

      double dt = tw_rand_exponential(lpptr->rng,mean_delay);
      trecv.t = now.t + dt;
      trecv.bits[0] = lpptr->gid;
      assert(trecv > now);
      if(trecv.t < g_tw_ts_end) {
	tw_event *evtptr = tw_event_new(lpptr->gid,trecv,lpptr);
	Event *evt = static_cast<Event *>(tw_event_data(evtptr));
	memcpy(evt->A,T,n*n*sizeof(*T));
	evt->type = Event::self;
	evt->cancel_count = 0;
	evt->canceller = 0;
	evt->cancelled_by_me = 0;
	evt->nretract = 0;
	evt->nforward = 0;
	evt->nbackward = 0;

	next_event = evt;
	assert(invcheck(n,evt->A));

	if(0) {
	  printf("Making event at time %.5f, to arrive at time %.5f, "
		 "from %d to self; evtptr = %llu\n",
		 now.t,trecv.t,(int) lpptr->gid,(unsigned long long int) evtptr);
	  matprint(n,evt->A,"A");
	}
	
	tw_event_send(evtptr);
      }
    }
  }
  
  ~matlp() {
    // Aggregate local statistics to global counters
    global_stat.nforward += nforward;
    global_stat.nbackward += nbackward;
    global_stat.ncommit += ncommit;

    global_stat.ncancelforward += ncancelforward;
    global_stat.ncancelbackward += ncancelbackward;
    global_stat.ncancelcommit += ncancelcommit;

    global_stat.nretract += nretract;

    {
      uint tmp[n*n];
      matmul(n,global_stat.A,M,tmp);

      if(0) {
	matprint(n,global_stat.A,"A");
	matprint(n,M,"M");
      }

      memcpy(global_stat.A,tmp,n*n*sizeof(*tmp));
    }

    delete[] T;
    delete[] M;
  }

  static void init(void *state,tw_lp *lpptr) {
    (void) new(state) matlp(lpptr,nmat);

    if(0)
      printf("Initializing state %llu for lp %llu with id %lld\n",
	     (unsigned long long int) state,
	     (unsigned long long int) lpptr,
	     (long long int) (lpptr->id));
  }
  
  static void finalize(void *state,tw_lp *lpptr) {
    static_cast<lp *>(state)->~lp();
  }

   
};


struct matlp : public lp {


  void forward(tw_lp *twlp,Event *evt,const Time& now) {
    evt->nforward++;
    if(evt->cancel_count != 0) { ncancelforward++; return; }
    apply(evt->A);
    nforward++;

    evt->cancelled_by_me = new std::list<Event *>();

    /* Make new event to send... */ {
      long long int nlp_total =
	((long long int) tw_nnodes()) *
	((long long int) g_tw_npe) *
	((long long int) g_tw_nlp);
      tw_stime now = tw_now(twlp);
      tw_stime trecv;
      int nsend;

      if(evt->type == Event::ping) {
	if(next_event != 0) {
	  next_event->cancel_count++;
	  if(next_event->cancel_count == 1)
	    next_event->canceller = evt;
	  evt->cancelled_by_me->push_back(next_event);
	  next_event = 0;
	}
	nsend = 1;
      } else if(evt->type == Event::self) {
	assert(evt == next_event);
	next_event = 0;
	nsend = 2;
      } else {
	assert(0);
      }

      for(int isend = 0; isend<nsend; isend++) {
	double dt = tw_rand_exponential(twlp->rng,mean_delay);
	tw_lpid dest = tw_rand_integer(twlp->rng,0,nlp_total-1);
	trecv.t = now.t + dt;
	trecv.bits[0] = twlp->gid;
	assert(trecv > now);

	if(isend == 0) dest = twlp->gid;

	if(trecv.t < g_tw_ts_end) {
	  tw_event *evtptr = tw_event_new(dest,trecv,twlp);
	  Event *evt = static_cast<Event *>(tw_event_data(evtptr));
	  memcpy(evt->A,T,n*n*sizeof(*T));
	  evt->cancel_count = 0;
	  evt->canceller = 0;
	  evt->cancelled_by_me = 0;

	  evt->nretract = 0;
	  evt->nforward = 0;
	  evt->nbackward = 0;
	  
	  if(isend == 0) {
	    next_event = evt;
	    evt->type = Event::self;
	  } else {
	    evt->type = Event::ping;
	  }
	  assert(invcheck(n,evt->A));
	  tw_event_send(evtptr);
	}
      }
    }

  }

  void backward(tw_lp *twlp,Event *evt,const Time& now) {
    evt->nbackward++;
    if(evt->cancel_count != 0) { ncancelbackward++; return; }

    int nsend;
    if(evt->type == Event::ping) {
      if(!evt->cancelled_by_me->empty()) {
	next_event = evt->cancelled_by_me->back();
	evt->cancelled_by_me->pop_back();
	if(next_event->cancel_count == 1)
	  next_event->canceller = 0;
	next_event->cancel_count--;
      } else
	next_event = 0;

      nsend = 1;
    } else if(evt->type == Event::self) {
      next_event = evt;
      nsend = 2;
    } else {
      assert(0);
    }
    
    for(int isend = 0; isend<nsend; isend++) {
      tw_rand_reverse_unif(twlp->rng);
      tw_rand_reverse_unif(twlp->rng);
    }

    assert(evt->cancelled_by_me->empty());
    delete evt->cancelled_by_me;
    unapply(evt->A);
    nbackward++;
  }

  tw_event *get_raw_evt_ptr(Event *evt) {
    size_t offset = (size_t) tw_event_data((tw_event *) (size_t) 0);
    return (tw_event *) (((size_t) evt) - offset);
  }
  void commit(tw_lp *twlp,Event *evt,const Time& now) {
    assert(evt->nforward - evt->nbackward == 1);
    assert(evt->cancel_count == 0 || evt->cancel_count == 1);
    if(evt->cancel_count > 0) { ncancelcommit++; return; }

    for(std::list<Event *>::iterator ii = evt->cancelled_by_me->begin();
	ii != evt->cancelled_by_me->end(); ++ii) {

      assert((*ii)->cancel_count == 1);
      assert((*ii)->canceller == evt);

      {
	tw_event *raw = get_raw_evt_ptr(*ii);
	tw_stime cmptime = tw_now(twlp);
	if(cmptime < twlp->pe->GVT) cmptime = twlp->pe->GVT;
	if(raw->recv_ts > cmptime) {
	  assert((*ii)->nforward == (*ii)->nbackward);
	  assert((*ii)->nretract == 0);
	  (*ii)->nretract++;

	  assert(raw->state.owner != TW_pe_free_q);
	  nretract++;
	  tw_event_cancel(get_raw_evt_ptr(*ii));
	}
      }
    }

    delete evt->cancelled_by_me;
    ncommit++;
  }
  
};


#endif
