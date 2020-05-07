/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#include <ross.h>
#include <vector>

#include "simtime.hh"
#include "reactions.hh"

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
      LP::init_lp<TTT>(state, me);
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

tw_lptype Specieslps[] = {
   { (init_f) RegisterBufferAndOffset<Species>::init_f,
     (pre_run_f) do_nothing,
     (event_f) LP::forward_event<Species>,
     (revent_f) LP::backward_event,
     (commit_f) LP::commit_event,
     (final_f) LP::final_lp<Species>,
     (map_f) ctr_map,
    //NULL,
    sizeof(Species)
  },
  {0}
};
template<>
std::vector<Species*> RegisterBufferAndOffset<Species>::buffer = std::vector<Species*>();
template<>
int RegisterBufferAndOffset<Species>::offset = 0;

tw_lptype Reactionlps[] = {
   { (init_f) RegisterBufferAndOffset<Reaction>::init_f,
     (pre_run_f) do_nothing,
     (event_f) LP::forward_event<Reaction>,
     (revent_f) LP::backward_event,
     (commit_f) LP::commit_event,
     (final_f) LP::final_lp<Reaction>,
     (map_f) ctr_map,
    //NULL,
    sizeof(Reaction)
  },
  {0}
};
template<>
std::vector<Reaction*> RegisterBufferAndOffset<Reaction>::buffer = std::vector<Reaction*>();
template<>
int RegisterBufferAndOffset<Reaction>::offset = 0;

enum SpeciesNames {
   L,
   R,
   LR,
   pLR,
   RR,
   NUM_SPECIES
};

enum ReactionNames {
   L_R__LR,
   LR__L_R,
   LR__pLR,
   pLR__LR,
   R_R__RR,
   RR__R_R,
   NUM_REACTIONS
};

double massAction(const int numSpecies, const SpeciesValue species[], const int numExpr, const ExpressionValue expr[])
{
   assert(numExpr == 2);
   double prod = expr[0]/expr[1];
   for (int ii=0; ii<numSpecies; ii++)
   {
      prod *= species[ii];
   }
   return prod;
}


int main(int argc, char **argv)
{
   // get rid of error if compiled w/ MEMORY queues
   g_tw_memory_nqueues=1;
   
   //printf("Calling tw_opt_add\n");
   //tw_opt_add(app_opt);
   
   //printf("Calling tw_init\n");
   tw_init(&argc, &argv);
   //printf("tw_init returned.\n");
   
   g_tw_memory_nqueues = 16; // give at least 16 memory queue event

   
   g_tw_events_per_pe = g_tw_nlp * 10;

   //FIXME, yeom2 .  I need to fill out these arrays with pointers to the actual objects.
   RegisterBufferAndOffset<Species>::init(NUM_SPECIES,0);
   RegisterBufferAndOffset<Reaction>::init(NUM_REACTIONS,NUM_SPECIES);
   g_tw_nlp = NUM_SPECIES+NUM_REACTIONS; 

   Species** species = RegisterBufferAndOffset<Species>::getBuffer();
   Reaction** rxn = RegisterBufferAndOffset<Reaction>::getBuffer();
   
   tw_define_lps(g_tw_nlp, sizeof(EventData), 0);
   
   for(int ii = 0; ii < NUM_SPECIES ; ii++) {
      tw_lp_settype(ii, &Specieslps[0]);
   }
   for(int ii = NUM_SPECIES; ii < NUM_SPECIES+NUM_REACTIONS; ii++) {
      tw_lp_settype(ii, &Reactionlps[0]);
   }
   
   if(g_tw_mynode == 0) {
      printf("========================================\n");
      printf("Particle Spammer Configuration\n");
      printf("   g_tw_nlp = %lld\n", (long long int) g_tw_nlp);
      printf("========================================\n\n");
   }

   //Set up tags
   for (int ii=0; ii<NUM_SPECIES; ii++) {
      species[ii]->tag = species[ii]->twlp->gid; 
   }
   for (int ii=0; ii<NUM_REACTIONS; ii++) {
      rxn[ii]->tag = rxn[ii]->twlp->gid; 
   }
   
   
   //Set up valences
   rxn[L_R__LR]->rxnValences[species[L]->tag] = -1;
   rxn[L_R__LR]->rxnValences[species[R]->tag] = -1;
   rxn[L_R__LR]->rxnValences[species[LR]->tag] = 1;

   rxn[LR__L_R]->rxnValences[species[LR]->tag] = -1;
   rxn[LR__L_R]->rxnValences[species[L]->tag] = 1;
   rxn[LR__L_R]->rxnValences[species[R]->tag] = 1;

   rxn[LR__pLR]->rxnValences[species[LR]->tag] = -1;
   rxn[LR__pLR]->rxnValences[species[pLR]->tag] = 1;

   rxn[pLR__LR]->rxnValences[species[pLR]->tag] = -1;
   rxn[pLR__LR]->rxnValences[species[LR]->tag] = 1;

   rxn[R_R__RR]->rxnValences[species[R]->tag] = -2;
   rxn[R_R__RR]->rxnValences[species[RR]->tag] = 1;

   rxn[RR__R_R]->rxnValences[species[RR]->tag] = -1;
   rxn[RR__R_R]->rxnValences[species[R]->tag] = 2;

   //connect the models to the species
   for (int ireaction=0; ireaction<NUM_REACTIONS; ireaction++)
   {
      //counter + unordered map here could cause chaos if position
      //dependence matters for the evaluation function.  For mass
      //action though it's not a problem.  When we have arbitrary
      //functions we'll need to fix this.
      int counter=0;
      for (auto iter : rxn[ireaction]->rxnValences) 
      {
         if (iter.second > 0)
         {
            rxn[ireaction]->dependentSpecies[iter.first] = counter++;
            CONNECT(*species[iter.first],valueChanged,rxn[ireaction]->tag,Reaction::UpdateSpeciesSlot);
         }
      }
      rxn[ireaction]->speciesValues.resize(counter);
   }

   //connect the reactions to each other.  Do this by scanning the reaction products against
   //the dependent reactions.
   //using O(n^2) loop here, fix with proper bipartite graph later.
   for (int ireaction=0; ireaction<NUM_REACTIONS; ireaction++)
   {
      std::vector<int> speciesForThisRxn;
      for (auto iter : rxn[ireaction]->rxnValences) {
         speciesForThisRxn.push_back(iter.first);
      }
      for (int jreaction=0; jreaction<NUM_REACTIONS; jreaction++)
      {
         if (ireaction == jreaction)
         {
            CONNECT(*rxn[ireaction],updateDependentRxns,rxn[ireaction]->tag,Reaction::UpdateReactionRateSlot);
         }
         else
         {
            for (auto speciesTag : speciesForThisRxn)
            {
               auto jter=rxn[jreaction]->dependentSpecies.find(speciesTag);
               if (jter != rxn[jreaction]->dependentSpecies.end()) {
                  CONNECT(*rxn[ireaction],updateDependentRxns,rxn[jreaction]->tag,Reaction::UpdateReactionRateSlot);
                  break;
               }
            }

         }
      }
   }
   
   //set up volumes
   double Volume=1;
   for (int ireaction=0; ireaction<NUM_REACTIONS; ireaction++)
   {
      rxn[ireaction]->exprValues.resize(2);
      rxn[ireaction]->exprValues[1] = Volume;
   }

   //set up rates
   rxn[L_R__LR]->exprValues[0] = 1;
   rxn[LR__L_R]->exprValues[0] = 1;
   rxn[LR__pLR]->exprValues[0] = 1;
   rxn[pLR__LR]->exprValues[0] = 1;
   rxn[R_R__RR]->exprValues[0] = 1;
   rxn[RR__R_R]->exprValues[0] = 1;
   
   //set up the functions
   for (int ireaction=0; ireaction<NUM_REACTIONS; ireaction++)
   {
      rxn[ireaction]->rateFunction = massAction;
   }

   //declare the initial conditions
   int initialConditions[NUM_SPECIES] = {0};
   initialConditions[L] = 1000;
   initialConditions[R] = 1000;

   //fire off the initial events
   for (int ispecies=0; ispecies<NUM_SPECIES; ispecies++)
   {
      species[ispecies]->value = initialConditions[ispecies];
   }
   for (int ireaction=0; ireaction<NUM_REACTIONS; ireaction++)
   {
      Reaction& thisReaction=*rxn[ireaction];
      for (auto iter : thisReaction.dependentSpecies)
      {
         //assumes iter.first == tag for the species.
         thisReaction.speciesValues[iter.second] = initialConditions[iter.first];
      }
      thisReaction.nextReaction = NULL;
      thisReaction.rate=1; //just need something other than 0 here to avoid a NaN.
      Time tBegin;
      tBegin.t=0;
      tBegin.bits[0]=tBegin.bits[1]=0;
      thisReaction.updateCountdown(tBegin);
      thisReaction.updateReactionRate(tBegin);
   }
   
   tw_run();

   tw_end();
   
   return 0;
}
