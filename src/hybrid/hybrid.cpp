/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#if defined(WCS_HAS_CONFIG)
#include "wcs_config.hpp"
#else
#error "no config"
#endif

#if defined(WCS_HAS_ROSS)

#include <string>
#include <cstring> // memset
#include <iostream>
#include <vector>
#include <getopt.h>
#include "utils/write_graphviz.hpp"
#include "utils/timer.hpp"
#include "utils/to_string.hpp"
#include "reaction_network/network.hpp"
#include "hybrid.hpp"
#include "hybrid/simtime.hh"
#include "hybrid/Funnel.hh"

#define OPTIONS "hi:s:t:"
static const struct option longopts[] = {
    {"help",    no_argument,  0, 'h'},
    {"interval", required_argument, 0, 'i'},
    {"seedbase", required_argument,  0, 's'},
    {"trust_region", required_argument,  0, 't'},
    { 0, 0, 0, 0 },
};

void do_nothing(void *,void *,void *) {}

tw_peid ctr_map(tw_lpid gid) {
  return (tw_peid) gid / g_tw_nlp;
}

template <typename TTT>
class RegisterBufferAndOffset {
 public:
   static std::vector<TTT*> buffer;
   static int offset;

   static void init_f(void *state, tw_lp* me)
   {
      LP<TTT>::init_lp(state, me);
      buffer[me->gid-offset] = reinterpret_cast<TTT*>(state);
   }

   static void init(const int bufferSize, const int _offset)
   {
      buffer.resize(bufferSize);
      offset = _offset;
   }
   static TTT** getBuffer() {
      return &buffer[0];
   }
};

tw_lptype Funnellps[] = {
   { (init_f) RegisterBufferAndOffset<Funnel>::init_f,
     (pre_run_f) do_nothing,
     (event_f) Funnel::forward_event,
     (revent_f) Funnel::backward_event,
     (commit_f) Funnel::commit_event,
     (final_f) Funnel::final_lp,
     (map_f) ctr_map,
    //NULL,
    sizeof(Funnel)
  },
  {0}
};
template<>
std::vector<Funnel*> RegisterBufferAndOffset<Funnel>::buffer = std::vector<Funnel*>();
template<>
int RegisterBufferAndOffset<Funnel>::offset = 0;

tw_lptype NRMlps[] = {
   { (init_f) RegisterBufferAndOffset<NextReactionMethodCrucible>::init_f,
     (pre_run_f) do_nothing,
     (event_f) NextReactionMethodCrucible::forward_event,
     (revent_f) NextReactionMethodCrucible::backward_event,
     (commit_f) NextReactionMethodCrucible::commit_event,
     (final_f) NextReactionMethodCrucible::final_lp,
     (map_f) ctr_map,
    //NULL,
    sizeof(NextReactionMethodCrucible)
  },
  {0}
};

template<>
std::vector<NextReactionMethodCrucible*> RegisterBufferAndOffset<NextReactionMethodCrucible>::buffer = std::vector<NextReactionMethodCrucible*>();
template<>
int RegisterBufferAndOffset<NextReactionMethodCrucible>::offset = 0;


tw_lptype Heartbeatlps[] = {
   { (init_f) RegisterBufferAndOffset<Heartbeat>::init_f,
     (pre_run_f) do_nothing,
     (event_f) Heartbeat::forward_event,
     (revent_f) Heartbeat::backward_event,
     (commit_f) Heartbeat::commit_event,
     (final_f) Heartbeat::final_lp,
     (map_f) ctr_map,
    //NULL,
    sizeof(Heartbeat)
  },
  {0}
};


template<>
std::vector<Heartbeat*> RegisterBufferAndOffset<Heartbeat>::buffer = std::vector<Heartbeat*>();
template<>
int RegisterBufferAndOffset<Heartbeat>::offset = 0;

// void post_lpinit_setup(tw_pe* pe);

tw_petype MyPEType = {
   NULL,
   post_lpinit_setup,
   NULL,
   NULL
};



//FIXME, use graph ids
enum ReactionEnum {
   A_TO_X,
   Y_2X_TO_3X,
   B_X_TO_Y_D,
   X_TO_E,
   numReactions
};

enum SpeciesEnum {
   A_idx,
   B_idx,
   D_idx,
   E_idx,
   X_idx,
   Y_idx,
   numSpecies
};



// FIXME, comment out, every stoic reference has to check the graph
std::vector<std::unordered_map<SpeciesTag, int>> stoic(numReactions);


//FIXME, all of these functions are linked in with dlsym trick

int factor = 200;
Real k_A_TO_X = 1;
Real k_Y_2X_TO_3X = 2.7574E-06/factor/factor;
Real k_B_X_TO_Y_D = 0.001660538/factor;
Real k_X_TO_E = 1;

Real function_A_TO_X(const std::vector<SpeciesValue>& species) {
   // for (auto specie : species){
   //    std::cout << specie << std::endl;
   // }
   return k_A_TO_X*species[0];
    
}

Real function_Y_2X_TO_3X(const std::vector<SpeciesValue>& species) {
   return k_Y_2X_TO_3X*species[0]*species[0]*species[1];
}

Real function_B_X_TO_Y_D(const std::vector<SpeciesValue>& species) {
   return k_B_X_TO_Y_D*species[0]*species[1];
}

Real function_X_TO_E(const std::vector<SpeciesValue>& species) {
   return k_X_TO_E*species[0];
}


void print_usage(const std::string exec, int code)
{
  std::cerr <<
    "Usage: " << exec << " <filename>.xml\n"
    "    Run the stochastic simulation algorithm (SSA) on a reaction\n"
    "    network graph (<filename>.xml).\n"
    "\n"
    "    OPTIONS:\n"
    "    -h, --help\n"
    "            Display this usage information\n"
    "\n"
    "    -i, --interval\n"
    "            Specify the interval\n"
    "\n"
    "    -s, --seedbase\n"
    "            Specify the seed for random number generator\n"
    "\n"
    "    -t, --trust_region\n"
    "            Specify the trust region\n"
    "\n";
  exit(code);
}


int main(int argc, char **argv)
{
   // get rid of error if compiled w/ MEMORY queues
   //g_tw_memory_nqueues=1;
   //g_tw_memory_nqueues = 16; // give at least 16 memory queue event
   
   //printf("Calling tw_opt_add\n");
   //tw_opt_add(app_opt);


  int c;
  int rc = EXIT_SUCCESS;
  std::string outfile;
  double interval = 0.2;
  int seedBase = 137;
  int trust_region = 20000;
  while ((c = getopt_long(argc, argv, OPTIONS, longopts, NULL)) != -1) {
    switch (c) {
      case 'h': /* --help */
        print_usage(argv[0], 0);
        break;
      case 'i': /* --interval */
        interval = atof(optarg);
        break;
      case 's': /* --seedbase */
        seedBase = static_cast<int>(atoi(optarg));
        break;
      case 't': /* --trust_region */
        trust_region = static_cast<int>(atoi(optarg));
        break;
      default:
        print_usage(argv[0], 1);
        break;
    }
  }
  if (optind != (argc - 1)) {
    print_usage (argv[0], 1);
  }
  std::string fn(argv[optind]);
  std::shared_ptr<wcs::Network> rnet_ptr = std::make_shared<wcs::Network>();
  wcs::Network& rnet = *rnet_ptr;
  rnet.load(fn);
  rnet.init();
  const wcs::Network::graph_t& g = rnet.graph();


  //FIXME, needs to come from graph => DONE
  int numReactions = rnet.reaction_list().size();
   
  //printf("Calling tw_init\n");
  tw_init(&argc, &argv);
  //printf("tw_init returned.\n");
  // g_tw_nlp = numReactions+numReactions+numReactions;
  g_tw_nlp = numReactions+numReactions+1; 
  g_tw_events_per_pe = g_tw_nlp * 100;

  //FIXME, yeom2 .  I need to fill out these arrays with pointers to the actual objects.
  RegisterBufferAndOffset<Funnel>::init(numReactions,0);
  RegisterBufferAndOffset<NextReactionMethodCrucible>::init(numReactions,numReactions);
  RegisterBufferAndOffset<Heartbeat>::init(1,2*numReactions);

  tw_define_lps(g_tw_nlp, messageSize());

  for(int ii = 0; ii < numReactions ; ii++) {
    tw_lp_settype(ii, &Funnellps[0]);
  }
  for(int ii = numReactions; ii < numReactions+numReactions; ii++) {
    tw_lp_settype(ii, &NRMlps[0]);
  }
  tw_lp_settype(2*numReactions, &Heartbeatlps[0]);

  tw_pe_settype(&MyPEType);
  // std::cout << "sampling end: " << g_st_sampling_end <<std::endl; 
  if (g_st_sampling_end > 0) {
    g_tw_ts_end = g_st_sampling_end; // default 100000
  } else {
    g_tw_ts_end = 1; // default 100000
  }
   
  // g_tw_events_per_pe =10; //default 800
  if( g_tw_mynode == 0 )
  {
    std::cout << "=========================================" << std::endl;
    std::cout << "WCS ROSS Configuration.............." << std::endl;
    // std::cout << "   run_id:\t" + std::string(run_id) << std::endl;
    std::cout << "   nlp_per_pe:\t" << g_tw_nlp << std::endl;
    std::cout << "   g_tw_ts_end:\t" << g_tw_ts_end << std::endl;
    std::cout << "   g_tw_events_per_pe:\t" << g_tw_events_per_pe << std::endl;
    
    std::cout << "   gvt-interval:\t" << g_tw_gvt_interval << std::endl;;
    std::cout << "   extramem:\t" << g_tw_events_per_pe_extra << std::endl;
    std::cout << "   ......................................" << std::endl;
    std::cout << "   Num nodes:\t" << tw_nnodes() << std::endl;
    // std::cout << "   Message size:\t" << sizeof(WCS_Message) << std::endl;
    std::cout << "========================================="
              << std::endl << std::endl;
  }

  
  tw_run();

  tw_end();
  
//  return 0;
return EXIT_SUCCESS;
}

void post_lpinit_setup(tw_pe* pe) {
  
  Funnel** funnels = RegisterBufferAndOffset<Funnel>::getBuffer();
  NextReactionMethodCrucible** crucibles = RegisterBufferAndOffset<NextReactionMethodCrucible>::getBuffer();
  Heartbeat** heartbeats = RegisterBufferAndOffset<Heartbeat>::getBuffer();


  //FIXME, output interval comes from options?
  double interval = 0.2;
  heartbeats[0]->setInterval(interval);

  //FIXME, seed comes from options
  int seedBase = 137;
  for (int ireaction=0; ireaction<numReactions; ireaction++) {
    crucibles[ireaction]->setTag(ireaction);
    crucibles[ireaction]->setSeed(seedBase+ireaction);

    funnels[ireaction]->setTag(ireaction);
    funnels[ireaction]->setSeed(seedBase+numReactions+ireaction);
  }


  //FIXME, inititial values are looked up in graph
  //setup initial conditions
  SpeciesValue initial_value[numSpecies];
  initial_value[A_idx] = 301*factor;
  initial_value[B_idx] = 1806*factor;
  initial_value[D_idx] = 0*factor;
  initial_value[E_idx] = 0*factor;
  initial_value[X_idx] = 1806*factor;
  initial_value[Y_idx] = 1806*factor;

  //FIXME, numerical parameter, comes from options
  int trust_region = 20000;


  //FIXME, run addSpecies for each dep species, add the rate function.
  
  //stoic[A_TO_X][A_idx] = -1;
  stoic[A_TO_X][X_idx] = 1;
  funnels[A_TO_X]->_rateFunction = &function_A_TO_X;
  funnels[A_TO_X]->addSpecies(A_idx, initial_value[A_idx], trust_region);
  //funnels[A_TO_X]->addSpecies(X_idx, initial_value[X_idx]);
  
  stoic[Y_2X_TO_3X][X_idx] = -2+3;
  stoic[Y_2X_TO_3X][Y_idx] = -1;
  funnels[Y_2X_TO_3X]->_rateFunction = function_Y_2X_TO_3X;
  funnels[Y_2X_TO_3X]->addSpecies(X_idx, initial_value[X_idx], trust_region);
  funnels[Y_2X_TO_3X]->addSpecies(Y_idx, initial_value[Y_idx], trust_region);

  //stoic[B_X_TO_Y_D][B_idx] = -1;
  stoic[B_X_TO_Y_D][X_idx] = -1;
  stoic[B_X_TO_Y_D][Y_idx] = 1;
  //stoic[B_X_TO_Y_D][D_idx] = 1;
  funnels[B_X_TO_Y_D]->_rateFunction = function_B_X_TO_Y_D;
  funnels[B_X_TO_Y_D]->addSpecies(B_idx, initial_value[B_idx], trust_region);
  //funnels[B_X_TO_Y_D]->addSpecies(D_idx, initial_value[D_idx], trust_region);
  funnels[B_X_TO_Y_D]->addSpecies(X_idx, initial_value[X_idx], trust_region);
  //funnels[B_X_TO_Y_D]->addSpecies(Y_idx, initial_value[Y_idx], trust_region);

  stoic[X_TO_E][X_idx] = -1;
  //stoic[X_TO_E][E_idx] = 1;
  funnels[X_TO_E]->_rateFunction = function_X_TO_E;
  //funnels[X_TO_E]->addSpecies(E_idx, initial_value[E_idx], trust_region);
  funnels[X_TO_E]->addSpecies(X_idx, initial_value[X_idx], trust_region);

  //Connect up the funnels to the crucibles they own.
  for( int ii=0; ii<numReactions; ii++) {
    CONNECT(*(funnels[ii]),Funnel::updateRateNRM,
            *(crucibles[ii]),NextReactionMethodCrucible::updatedRate);
  }

  //FIXME, use graph instead of bruteforcing all possible connections
  //connect the reactions to each other.  Do this by scanning the reaction products against
  //the dependent reactions.
  //using O(n^2) loop here, fix with proper bipartite graph later.
  for (int ireaction=0; ireaction<numReactions; ireaction++)
  {
    std::vector<int> speciesForThisRxn;
    //replace stoic instance, possibly do better than a n^2 loop because we have a graph
    for (auto iter : stoic[ireaction]) {
        speciesForThisRxn.push_back(iter.first);
    }
    funnels[ireaction]->setupRecvFromReaction(ireaction);
    for (int jreaction=0; jreaction<numReactions; jreaction++)
    {
        for (auto speciesTag : speciesForThisRxn) {
          if (funnels[jreaction]->dependsOnSpecies(speciesTag)) {
              CONNECT(*(crucibles[ireaction]),NextReactionMethodCrucible::fire,
                      *(funnels[jreaction]),Funnel::firedReaction);
              if (ireaction != jreaction) {
                funnels[ireaction]->setupSendToReaction(jreaction, funnels[jreaction]->getID());
                funnels[jreaction]->setupRecvFromReaction(ireaction);
              }
              break;
          }
        }
    }
  }

  //connect the printing to the heartbeat for output
  CONNECT((*heartbeats[0]),Heartbeat::activate,(*heartbeats[0]),Heartbeat::activated);
  for (int ireaction=0; ireaction<numReactions; ireaction++) {
    CONNECT(*(heartbeats[0]),Heartbeat::activate,
            *(funnels[ireaction]),Funnel::printData);
    
  }
  
  //setup the initial messages
  Time tnow = 0;
  for (int ireaction=0; ireaction<numReactions; ireaction++) {
    funnels[ireaction]->beginReactions(tnow);

  }
  NoopMsg msg;
  heartbeats[0]->activate(tnow, msg);
}




#endif // defined(WCS_HAS_ROSS)
