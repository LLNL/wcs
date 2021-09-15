#pragma once

#include "actor.hh"
#include "reaction_network/network.hpp"
// #include "utils/graph_factory.hpp"
// #include <boost/graph/graphml.hpp>
#include <iostream>
#include <random>
#include <functional>
#include "wcs_types.hpp"

#define MAX_TIME 100 //FIXME!

typedef wcs::Network::v_desc_t SpeciesTag;
typedef std::size_t ReactionTag;
typedef wcs::reaction_rate_t  Real;
typedef std::size_t SpeciesValue;


typedef std::function<wcs::reaction_rate_t (const std::vector<wcs::reaction_rate_t>&)> RateFunction;


struct NoopMsg {
};
struct RateMsg {
   Real rate;
};
struct ReactionMsg {
   ReactionTag reaction;
};
struct FunnelRateExchangeMsg {
   ReactionTag reaction;
   bool nrmMode;
   Real rate;
   Real normalContrib;
};
struct SwitchModeMsg {
   ReactionTag reaction;
   bool nrmMode;
};

class NextReactionMethodCrucible : public LP<NextReactionMethodCrucible> {
 public:

   ACTOR_DEF();
      
   SLOT(firedMyself, NoopMsg);
   void FORWARD(firedMyself)(Time tnow, NoopMsg& msg, tw_bf*) {
      ReactionMsg r;
      r.reaction = _rxntag;
      fire(0, r);
      _countdown = newCountdown();
      fireNextReaction(tnow);
   }
   void COMMIT(firedMyself)(Time tnow, NoopMsg& msg, tw_bf*) {
      //std::cout << tnow << ": Fired Reaction " << _rxntag << "!\n";
   }
   template<typename Checkpointer>
   void CHECKPOINT(firedMyself)(Checkpointer& ar, Time tnow, NoopMsg& msg, tw_bf*) {
      ar & _countdown & _messageOutstanding & _lastUpdatedTime & _simpleRand;
   }

   SLOT(updatedRate, RateMsg);
   void FORWARD(updatedRate)(Time tnow, RateMsg& msg, tw_bf*) {
      if (_rate != 0) {
         EVENT_CANCEL(_messageOutstanding);
         _countdown -= countdownChange(tnow);
      }
      _rate = msg.rate;
      if (_rate != 0) {
         fireNextReaction(tnow);
      }
   }
   template<typename Checkpointer>
   void CHECKPOINT(updatedRate)(Checkpointer& ar, Time tnow, RateMsg& msg, tw_bf*) {
      ar & _messageOutstanding & _countdown & _rate & _lastUpdatedTime & _simpleRand;
   }

   NextReactionMethodCrucible() {
      _rate = 0;
      _countdown = newCountdown();
      _messageOutstanding = PDES_NULL_EVENT;
   }
   
   SIGNAL(fire,ReactionMsg);

   void setTag(const ReactionTag mytag) {
      _rxntag = mytag;
   }

   void setSeed(const int seed) {
      _simpleRand.seed(seed);
   }
   
 private:
   ReactionTag _rxntag;
   Real _countdown;
   Real _rate;
   std:: mt19937 _simpleRand;
   PdesEvent _messageOutstanding;
   Real _lastUpdatedTime;
   // PdesLpid _self_lpid;
   
   void fireNextReaction(Time tnow) {
      _lastUpdatedTime = tnow;
      Time timeToFire = firingTime(_countdown,_rate);
      NoopMsg msg;
      _messageOutstanding = EVENT_SEND(twlp->gid, firedMyself, timeToFire, msg);
   }

   Real newCountdown() {
      //From Anderson, Algo 3
      //choose a uniform random number 0-1
      // std::random_device rd;  //Will be used to obtain a seed for the random number engine
      // std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
      std::uniform_real_distribution<> dis(0.0, 1.0);
      Real uniform_random_number;
      uniform_random_number = dis(_simpleRand);

      // std::cout << "Uniform: " << uniform_random_number << '\n';
      return log(1/uniform_random_number);
   }

   Real countdownChange(Time tnow) {
      return _rate*(tnow-_lastUpdatedTime);
   }

   Real firingTime(Real countdown, Real rate) {
      return countdown/rate;
   }
};

extern
std::vector<std::unordered_map<SpeciesTag, int>> stoic;

class Funnel : public LP<Funnel> {
 public:

   ACTOR_DEF();

   Funnel() {
      //Things needed to initiate the model
      _rateNRM = 0;
      _messageOutstanding = PDES_NULL_EVENT;
      _numReactionsInBox = 0;
         

      //tolerances the users can adjust
      _tolerance = 0.003; 
      _sigma_max = 6;
      _numReactionCutoff = 100;
      _maxTimestep = MAX_TIME;
   }


   SLOT(firedReaction, ReactionMsg);
   void FORWARD(firedReaction)(Time tnow, ReactionMsg& msg, tw_bf* twbf) {
      //update species
      for (const auto& kv : stoic[msg.reaction]) {
         auto iter = _indexFromSpecTag.find(kv.first);
         if (iter != _indexFromSpecTag.end()) {
            _species[iter->second] += kv.second;
         }
      }
      updateSpeciesBasedOnRates(tnow);
      bool allInBounds = checkRateBounds(tnow);
      if (!allInBounds) {
         interruptNextWakeup(tnow);
         setupFutureEvents(tnow, true);
      }
   }
   void BACKWARD(firedReaction)(Time tnow, ReactionMsg& msg, tw_bf*) {
   }
   void COMMIT(firedReaction)(Time tnow, ReactionMsg& msg, tw_bf*) {
      //std::cout << tnow << ": NRM " << msg.reaction << "->" << _rxntag << std::endl;
      //printSpecies(tnow);
   }
   template<typename Checkpointer>
   void CHECKPOINT(firedReaction)(Checkpointer& ar, Time tnow, ReactionMsg& msg, tw_bf*) {
   }

   
   SLOT(changedRate, FunnelRateExchangeMsg);
   void FORWARD(changedRate)(Time tnow, FunnelRateExchangeMsg& msg, tw_bf* twbf) {
      updateSpeciesBasedOnRates(tnow);
      updateSavedRates(msg.reaction, msg.nrmMode, msg.rate, msg.normalContrib);
      interruptNextWakeup(tnow);
      //std::cout << tnow << ": %%% " << ((tnow-_lastRecomputeTime)>=_maxTimestep) << " " << checkRateBounds(tnow) << std::endl;
      setupFutureEvents(tnow, (tnow-_lastRecomputeTime)>=_maxTimestep || !checkRateBounds(tnow));
   }
   void BACKWARD(changedRate)(Time tnow, FunnelRateExchangeMsg& msg, tw_bf*) {
   }
   void COMMIT(changedRate)(Time tnow, FunnelRateExchangeMsg& msg, tw_bf*) {
      // std::cout << tnow << ": Changed rate for " << msg.reaction << "->"
      //          << _rxntag << ": "
      //          << (msg.nrmMode ? "NRM " : "SDE ") 
      //          << msg.rate << " "
      //          << msg.normalContrib << std::endl;
      printSpecies(tnow);
   }
   template<typename Checkpointer>
   void CHECKPOINT(changedRate)(Checkpointer& ar, Time tnow, FunnelRateExchangeMsg& msg, tw_bf*) {
   }

   SLOT(wakeup, NoopMsg);
   void FORWARD(wakeup)(Time tnow, NoopMsg& msg, tw_bf* twbf) {
      // This gets called when we're sure that we're going to go
      // out-of-bounds due to our tau condition, based on what we
      // currently know about the SDE rates.
      updateSpeciesBasedOnRates(tnow);
      setupFutureEvents(tnow, (tnow-_lastRecomputeTime)>=_maxTimestep || !checkRateBounds(tnow));
   }
   void COMMIT(wakeup)(Time tnow, NoopMsg& msg, tw_bf*) {
      //std::cout << tnow << ": Woke Reaction " << _rxntag << "!\n";
      //printSpecies(tnow);
   }
   template <typename Checkpointer>
   void CHECKPOINT(wakeup)(Checkpointer& ar, Time tnow, NoopMsg& msg, tw_bf*) {
   }

   SLOT(changedMode, SwitchModeMsg);
   void FORWARD(changedMode)(Time tnow, SwitchModeMsg& msg, tw_bf*) {
      updateSpeciesBasedOnRates(tnow);
      int irxn = _indexFromRxnTag[msg.reaction];
      updateSavedRates(msg.reaction, msg.nrmMode, _rxnRate[irxn], _normalContrib[irxn]);
      bool allInBounds = checkRateBounds(tnow);
      if (!allInBounds) {
         interruptNextWakeup(tnow);
         setupFutureEvents(tnow, true);
      }
   }
   void COMMIT(changedMode)(Time tnow, SwitchModeMsg& msg, tw_bf*) {
      //std::cout << tnow << ": Changed mode " << msg.reaction << "->" << _rxntag << ": " << msg.nrmMode << std::endl;
   }
   template <typename Checkpointer>
   void CHECKPOINT(changedMode)(Checkpointer& ar, Time tnow, SwitchModeMsg& msg, tw_bf*) {
   }

   SLOT(printData, NoopMsg);
   void FORWARD(printData)(Time tnow, NoopMsg& msg, tw_bf*) {
   }
   void COMMIT(printData)(Time tnow, NoopMsg& msg, tw_bf*) {
      printSpecies(tnow);
   }
   
   SIGNAL(updateRateNRM, RateMsg);
   SIGNAL(updateRateFunnel, FunnelRateExchangeMsg);
   SIGNAL(updateModeFunnel, SwitchModeMsg);

   void addSpecies(SpeciesTag tag, SpeciesValue initialValue, int trustRegion) {
      int oldSize = _indexFromSpecTag.size();
      _indexFromSpecTag[tag] = oldSize;
      _tagFromSpecIndex.push_back(tag);
      _species.push_back(initialValue);
      _lowerBound.push_back(0);
      _upperBound.push_back(0);
      _speciesRate.push_back(0);
      _trustRegion.push_back(trustRegion);
   }

   void setTag(ReactionTag tag) {
      _rxntag = tag;
   }

   void setSeed(const int seed) {
      _simpleRand.seed(seed);
   }

   void setupSendToReaction(ReactionTag otherRxnTag, tw_lpid destID) {
      if (otherRxnTag != _rxntag) {
         CONNECT(*this,Funnel::updateRateFunnel,destID,Funnel::changedRate);
         CONNECT(*this,Funnel::updateModeFunnel,destID,Funnel::changedMode);
      }
   }

   void setupRecvFromReaction(ReactionTag otherRxnTag) {
      int index = _tagFromRxnIndex.size();
      _tagFromRxnIndex.push_back(otherRxnTag);
      _indexFromRxnTag[otherRxnTag] = index;
      _sqrtTimestep.push_back(0);
      _normalContrib.push_back(0);
      _rxnRate.push_back(0);
      _rxnOverflow.push_back(0);
      _rxnNrmMode.push_back(0);
      if (otherRxnTag == _rxntag) {
         _thisrxnindex = index;
      }
   }
      
   bool dependsOnSpecies(int specTag) {
      return _indexFromSpecTag.find(specTag) != _indexFromSpecTag.end();
   }
   
   void beginReactions(Time tnow) {
      _lastRecomputeTime = tnow;
      _lastUpdatedTime = tnow;

      //recalculate rate
      setupFutureEvents(tnow, true);
      printSpecies(tnow);
   }

   void printSpecies(Time tnow) {
      // std::cout << tnow << ": "
      //           << "Reaction " << _rxntag << ": (wakeup " << _wakeupTime << ") (nrm=" << _rxnNrmMode[_thisrxnindex] << " rate=" << _rxnRate[_thisrxnindex] << ")" ;
      // for (std::size_t ii=0; ii<_tagFromSpecIndex.size(); ii++) {
      //    std::cout << " {" << _tagFromSpecIndex[ii] << ": "
      //              << "(" << _lowerBound[ii] << "<"
      //              << _species[ii]
      //              << "<" << _upperBound[ii] << ") "
      //              << _speciesRate[ii]
      //              << "}";
      // }
      // std::cout << std::endl;
      int both_species = 0;
      for (int ii=0; ii<_tagFromSpecIndex.size(); ii++) {
         if (_tagFromSpecIndex[ii] == 2 || _tagFromSpecIndex[ii] == 4){
            both_species = both_species +1;
         }
      }
      if (both_species == 2 ) {
         std::cout << tnow << " ";
         for (int ii=0; ii<_tagFromSpecIndex.size(); ii++) {
            if (_tagFromSpecIndex[ii] == 2 || _tagFromSpecIndex[ii] == 4){
               std::cout << _tagFromSpecIndex[ii] << ": " << _species[ii] << " ";	
            }
         }
         std::cout << std::endl;
      }

   }
   
 public:
   
   RateFunction _rateFunction; //function to use for rate calculations
   Real _tolerance; // tolerance for the tau independence condition
   Real _sigma_max; // tolerance for switching between SDE and NRM
   int _numReactionCutoff; // min number of reactions before we consider SDE
   Real _maxTimestep; // maximum timestep allowed in SDE mode.


 private:
   std::unordered_map<SpeciesTag, int> _indexFromSpecTag; // internal index from tag
   std::vector<SpeciesTag> _tagFromSpecIndex; // tag from internal index
   std::vector<Real> _species; //species counts
   std::vector<Real> _speciesRate; //first derivative of the species
   std::vector<Real> _lowerBound; // lb of the tau independence region
   std::vector<Real> _upperBound; // ub of the tau independence region
   std::vector<Real> _trustRegion; //how far should we trust this species for a trust region?  Think of this as a max for the ub/lb, independent of the derivative

   std::unordered_map<ReactionTag, int> _indexFromRxnTag; //internal index from tag
   std::vector<ReactionTag> _tagFromRxnIndex; // tag from internal index
   std::vector<Real> _sqrtTimestep; // sqrt of the timestep, used for SDE
   std::vector<Real> _normalContrib; // sqrt(rate)*normal, used for SDE
   std::vector<Real> _rxnOverflow; // Accumulator for SDE updates
   std::vector<Real> _rxnRate; // rate for each reaction
   std::vector<int> _rxnNrmMode; // boolean, is the rxn in nrmMode
   
   Real _lastUpdatedTime;
   Real _lastRecomputeTime;
   Real _wakeupTime;

   ReactionTag _rxntag;
   int _thisrxnindex;

   Real _rateNRM;
   PdesEvent _messageOutstanding;
   Real _numReactionsInBox;

   std:: mt19937 _simpleRand; 
   
   Real recalculateRate(Time tnow) {
      //recalculate rate
      Real newRate = _rateFunction(_species);
      _numReactionsInBox = _numReactionCutoff+1;
      //compute gradient
      std::vector<Real> gradient(_species.size());
      for (unsigned int ii=0; ii<_species.size(); ii++) {
         std::vector<Real> modSpecies(_species);
         modSpecies[ii]++;
         Real gradRate = _rateFunction(modSpecies);
         gradient[ii] = gradRate-newRate;
      }

      int nonZeroGradCount = 0;
      for (auto grad : gradient) {
         if (grad!=0) {
            nonZeroGradCount++;
         }
      }
      //compute bounding box
      for (unsigned int ii=0; ii<_species.size(); ii++) {
         Real maxDelta = _trustRegion[ii];
         Real delta;
         if (gradient[ii] == 0) {
            delta = maxDelta;
         } else {
            delta = _tolerance*newRate/fabs(gradient[ii])/sqrt(nonZeroGradCount)/2;  
            if (delta > maxDelta) {
               delta = maxDelta;
            }
         }
         _lowerBound[ii] = _species[ii]-delta;
         _upperBound[ii] = _species[ii]+delta;

         SpeciesTag specTag = _tagFromSpecIndex[ii];
         auto iter = stoic[_rxntag].find(specTag);
         if (iter != stoic[_rxntag].end()) {
            int sss = ((iter->second >= 0) ? iter->second : -iter->second);
            if (sss) {
               _numReactionsInBox = std::min(_numReactionsInBox, delta/sss);
            }
         }
      }
      _lastRecomputeTime = tnow;
      return newRate;
   }

   bool checkRateBounds(Time tnow) const {
      //are we within our bounds
      bool allWithinBounds=true;
      for (unsigned int ii=0; ii<_species.size(); ii++) {
         if (!(_lowerBound[ii] <= _species[ii] && _species[ii] <= _upperBound[ii])) {
            allWithinBounds = false;
            break;
         }
      }
      return allWithinBounds; 
   }

   void updateSpeciesBasedOnRates(Time tnow) {

      // Use SDE euler method to update the rates for things in SDE
      // mode.
      if (tnow <= _lastUpdatedTime) {
         return;
      }
      
      for (unsigned int irxn=0; irxn<_tagFromRxnIndex.size(); irxn++) {
         if (_rxnNrmMode[irxn]==false && _rxnRate[irxn] != 0) {
            _rxnOverflow[irxn] += _rxnRate[irxn]*(tnow-_lastUpdatedTime);
            if (_normalContrib[irxn] != 0) {
               Real c_old = _sqrtTimestep[irxn];
               Real c_now = std::sqrt(tnow - _lastUpdatedTime + c_old * c_old);
               _rxnOverflow[irxn] += _normalContrib[irxn]*(c_now-c_old);
               _sqrtTimestep[irxn] = c_now;
            }
            if (_rxnOverflow[irxn] >= 1) {
               int rxnCount = int(_rxnOverflow[irxn]);
               for (const auto& kv : stoic[_tagFromRxnIndex[irxn]]) {
                  auto iter = _indexFromSpecTag.find(kv.first);
                  if (iter != _indexFromSpecTag.end()) {
                     _species[iter->second] += rxnCount*kv.second;
                  }
               }
               _rxnOverflow[irxn] -= rxnCount;
            }
         }
      }
      _lastUpdatedTime = tnow;
   }

   void updateSavedRates(ReactionTag otherRxnTag, bool nrmMode, Real rate, Real normalContrib) {
      int irxn = _indexFromRxnTag[otherRxnTag];
      _normalContrib[irxn] = normalContrib; 
      _sqrtTimestep[irxn] = 0;
      Real delta_rate = rate-_rxnRate[irxn];
      _rxnRate[irxn] = rate;
      _rxnNrmMode[irxn] = nrmMode;
      if (delta_rate!=0) {
         for (const auto& kv : stoic[otherRxnTag]) {
            auto iter = _indexFromSpecTag.find(kv.first);
            if (iter != _indexFromSpecTag.end()) {
               _speciesRate[iter->second] += kv.second*delta_rate;
            }
         }
      }
   }
   
   void interruptNextWakeup(Time tnow) {
      EVENT_CANCEL(_messageOutstanding);
   }      
   
   void scheduleNextWakeup(Time tnow, Real delta_time) {
      NoopMsg msg;
      _messageOutstanding = EVENT_SEND(twlp->gid, wakeup, delta_time, msg);
      _wakeupTime = tnow+delta_time;
   }

   void setupFutureEvents(Time tnow, bool recompute) {

      bool nrmMode;
      Real delta_time;
      Real nextRate;
      if (recompute) {
         nextRate = recalculateRate(tnow);
      } else {
         nextRate = _rxnRate[_thisrxnindex];
      }         
      getNextMode(nextRate, delta_time, nrmMode);
      if (recompute) {
         FunnelRateExchangeMsg newRateMsg;
         newRateMsg.reaction = _rxntag;
         newRateMsg.rate = nextRate;
         newRateMsg.nrmMode = nrmMode;

         std::normal_distribution<> dis(0,1);
         Real unit_normal_random = dis(_simpleRand);
         newRateMsg.normalContrib = unit_normal_random * std::sqrt(nextRate);
         
         updateSavedRates(newRateMsg.reaction, newRateMsg.nrmMode, newRateMsg.rate, newRateMsg.normalContrib);
         updateRateFunnel(0, newRateMsg);
      }
      Real newRateNRM = (nrmMode==true ? nextRate : 0);
      if (newRateNRM != _rateNRM) {
         RateMsg newRateMsg;
         newRateMsg.rate = newRateNRM;
         updateRateNRM(0, newRateMsg);
         _rateNRM = newRateNRM;

         if (!recompute) {
            SwitchModeMsg newSwitchMsg;
            newSwitchMsg.reaction = _rxntag;
            newSwitchMsg.nrmMode = nrmMode;
            updateSavedRates(_rxntag, nrmMode, nextRate, _normalContrib[_thisrxnindex]);
            updateModeFunnel(0, newSwitchMsg);
         }
      }
      scheduleNextWakeup(tnow, delta_time);
   }
   
   void getNextMode(Real nextRate, Real& delta_time, bool& nextNrmMode) const {

      delta_time = findWakeupInterval(nextRate);
      if (_numReactionsInBox > _numReactionCutoff && nextRate*(delta_time+_sqrtTimestep[_thisrxnindex]*_sqrtTimestep[_thisrxnindex]) > _sigma_max*_sigma_max) {
         nextNrmMode = false;
      } else {
         nextNrmMode = true;
      }      
   }

   Real findWakeupInterval(Real nextRate) const {
      // find when species will cross boundary
      Real minTime = _maxTimestep;
      Real delta_rate = nextRate-(_rxnRate[_thisrxnindex]);
      
      for (unsigned int ispec=0; ispec<_species.size(); ispec++) {
         Real thisRate = _speciesRate[ispec];
         if (delta_rate != 0) {
            auto iter=stoic[_rxntag].find(_tagFromSpecIndex[ispec]);
            if (iter != stoic[_rxntag].end()) {
               thisRate += delta_rate*iter->second;
            }
         }
         if (thisRate > 0) {
            minTime = std::min(minTime,(_upperBound[ispec]+1 - _species[ispec]) / thisRate);
         } else if (thisRate < 0) {
            minTime = std::min(minTime,(_lowerBound[ispec]-1 - _species[ispec]) / thisRate);
         } else {
            //do nothing
         }
         //std::cout << "%%% " << thisRate << " " << minTime << std::endl;
      }
      return minTime;
   }
};


class Heartbeat : public LP<Heartbeat> {
 public:

   ACTOR_DEF();
   
   void setInterval(double interval) { _interval = interval; }
   
   SLOT(activated, NoopMsg);
   void FORWARD(activated)(Time tnow, NoopMsg& msg, tw_bf*) {
      activate(_interval, msg);
   }

   SIGNAL(activate, NoopMsg);
   
 private:
   double _interval;
};
