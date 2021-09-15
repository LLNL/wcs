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
#include "reaction_network/species.hpp"
#include "reaction_network/reaction.hpp"
#include "reaction_network/vertex.hpp"
#include "hybrid.hpp"
#include "hybrid/simtime.hh"
#include "hybrid/Funnel.hh"

#define OPTIONS "hi:s:r:t:"
static const struct option longopts[] = {
    {"help",    no_argument,  0, 'h'},
    {"interval", required_argument, 0, 'i'},
    {"seedbase", required_argument,  0, 's'},
    {"trust_region", required_argument,  0, 'r'},
    {"simulation_time", required_argument,  0, 't'},
    { 0, 0, 0, 0},
};

using r_prop_t = wcs::Network::r_prop_t;
// /// The type of BGL vertex descriptor for graph_t
using v_desc_t = boost::graph_traits<wcs::Network::graph_t>::vertex_descriptor;
double interval = 0.2;
int seedBase = 137;
int trust_region = 20000;
int simulation_time = 100;
int numReactions = 0;
int numSpecies = 0;
const wcs::Network::reaction_list_t* reaction_list;
const wcs::Network::species_list_t* species_list;
const LIBSBML_CPP_NAMESPACE::Model* model;
std::string filename = "";
wcs::Network rnet;

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
     (event_f) Funnel::FORWARD_event,
     (revent_f) Funnel::BACKWARD_event,
     (commit_f) Funnel::COMMIT_event,
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
     (event_f) NextReactionMethodCrucible::FORWARD_event,
     (revent_f) NextReactionMethodCrucible::BACKWARD_event,
     (commit_f) NextReactionMethodCrucible::COMMIT_event,
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
     (event_f) Heartbeat::FORWARD_event,
     (revent_f) Heartbeat::BACKWARD_event,
     (commit_f) Heartbeat::COMMIT_event,
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


std::vector<std::unordered_map<SpeciesTag, int>> stoic;


void print_usage(const std::string exec, int code)
{
  std::cerr <<
    "Usage: " << exec << " <filename>.xml\n"
    "    Run the simulation using the hybrid simulation mode which combines \n"
    "    Stochastic differential equations (SDEs) with stochastic simulation \n" 
    "    algorithm (SSA) on a reaction network graph (<filename>.xml).\n"
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
      case 'r': /* --trust_region */
        trust_region = static_cast<int>(atoi(optarg));
        break;
      case 't': /* --simulation_time */
        simulation_time = static_cast<int>(atoi(optarg));
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
  filename = fn;
  rnet.load(fn);
  rnet.init();
  const wcs::Network::graph_t& g = rnet.graph();

  LIBSBML_CPP_NAMESPACE::SBMLReader reader;
  const LIBSBML_CPP_NAMESPACE::SBMLDocument* document
    = reader.readSBML(fn);

  model = document->getModel();


  //FIXME, needs to come from graph => DONE

  reaction_list = &rnet.reaction_list();
  species_list = &rnet.species_list();  
  numReactions = reaction_list->size();
  numSpecies = species_list->size();
   
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
  g_st_sampling_end = simulation_time; 
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


  heartbeats[0]->setInterval(interval);

  for (int ireaction=0; ireaction<numReactions; ireaction++) {
    crucibles[ireaction]->setTag(ireaction);
    crucibles[ireaction]->setSeed(seedBase+ireaction);

    funnels[ireaction]->setTag(ireaction);
    funnels[ireaction]->setSeed(seedBase+numReactions+ireaction);
  }


  const wcs::Network::graph_t& g = rnet.graph();

  const wcs::Network::species_list_t& species_list = rnet.species_list(); 
  const wcs::Network::reaction_list_t& reaction_list = rnet.reaction_list(); 
  const LIBSBML_CPP_NAMESPACE::ListOfReactions* model_reactions
    = model->getListOfReactions();
  const LIBSBML_CPP_NAMESPACE::ListOfSpecies* model_species
    = model->getListOfSpecies();

  
  for (size_t i = 0u; i < reaction_list.size(); ++i) {
    // std::cout << "reaction" << i << std::endl; 
    const auto& vd = reaction_list[i];

    const auto &reaction = *(model_reactions->get(g[vd].get_label()));
    const LIBSBML_CPP_NAMESPACE::ListOfSpeciesReferences* reaction_reactants
    = reaction.getListOfReactants(); 
    const LIBSBML_CPP_NAMESPACE::ListOfSpeciesReferences* reaction_modifiers
    = reaction.getListOfModifiers();

    // std::cout << "Reaction label " << g[vd].get_label() << std::endl;
    const auto& rp = g[vd].property<r_prop_t>();
    funnels[i]->_rateFunction = rp.get_calc_rate_fn();
    // reactant species
    // std::cout << "reactants  " << std::endl;
    for (const auto ei_in : boost::make_iterator_range(boost::in_edges(vd, g))) {
      const auto vd_reactant = boost::source(ei_in, g);
      const auto& sv_reactant = g[vd_reactant];

      const auto& sp_reactant = sv_reactant.property<wcs::Species>();
      const auto stoichio = g[ei_in].get_stoichiometry_ratio();
      // Include only reactants without the constant reactants/modifiers which have been converted into modifiers
      if (reaction_reactants->get(sv_reactant.get_label()) != NULL || reaction_modifiers->get(sv_reactant.get_label()) != NULL) {
        funnels[i]->addSpecies(vd_reactant, sp_reactant.get_count(), trust_region);
      }
    }
  }
   
  //Connect up the funnels to the crucibles they own.
  for( int ii=0; ii<numReactions; ii++) {
    CONNECT(*(funnels[ii]),Funnel::updateRateNRM,
            *(crucibles[ii]),NextReactionMethodCrucible::updatedRate);
  }

  for (size_t ireaction = 0u; ireaction < reaction_list.size(); ++ireaction) {
    funnels[ireaction]->setupRecvFromReaction(ireaction);
    std::unordered_map<SpeciesTag, wcs::stoic_t> dependentSpecies;
    std::unordered_map<SpeciesTag, wcs::stoic_t>::iterator speciesit;

    const auto& vd = reaction_list[ireaction];

    // reactant species
    for (const auto ei_in : boost::make_iterator_range(boost::in_edges(vd, g))) {
      const auto vd_reactant = boost::source(ei_in, g);
      const auto& sv_reactant = g[vd_reactant];

      const auto stoichio = g[ei_in].get_stoichiometry_ratio();
      // Include only reactants without modifiers
      if (stoichio > 0) {
        dependentSpecies.insert(std::make_pair(vd_reactant, stoichio * -1));
      }
    }
    // product species
    for (const auto ei_out : boost::make_iterator_range(boost::out_edges(vd, g))) {
      const auto vd_product = boost::target(ei_out, g);
      const auto& sv_product = g[vd_product];

      const auto stoichio = g[ei_out].get_stoichiometry_ratio();
      speciesit = dependentSpecies.find(vd_product);
      if (speciesit == dependentSpecies.cend()) {
        dependentSpecies.insert(std::make_pair(vd_product, stoichio)); 
      } else {
        speciesit->second += stoichio;
      }
    }
    stoic.push_back(dependentSpecies);
    for (size_t jreaction = 0u; jreaction < reaction_list.size(); ++jreaction) {
      // const auto& vdj = reaction_list[jreaction];
      for (auto speciesTag : dependentSpecies) {
        if (funnels[jreaction]->dependsOnSpecies(speciesTag.first)) {
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
