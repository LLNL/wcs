#pragma once

#include <vector>
#include <unordered_map>
#include <list>
#include <ross.h>
// #include "wcs-ross-bf.hpp"
#include "simtime.hh"
#include "Checkpoint.hh"

//// Classes that make the whole API work

struct EventData; //will be defined by a perl script. 
typedef std::size_t Tag;
typedef tw_lpid PdesLpid;
typedef tw_event* PdesEvent;
#define PDES_NULL_EVENT NULL

class BaseLP {
 public:
   virtual ~BaseLP() {}

   struct ObserverData {
      tw_lpid dest;
      Tag slot;
   };
   tw_lp* twlp; //this has to be initialized properly! TODO

   tw_lpid getID() const;
};

extern CheckpointQueue checkpointQueue;

template <typename Dependent>
class LP : public BaseLP {
 public:
   typedef void (Dependent::*MethodPtr)(Time tnow, void*, tw_bf*);
   static void init_lp(void *state, tw_lp *me);
   static void forward_event(void *state, tw_bf *bf, void *evt, tw_lp *thislp);
   static void backward_event(void *state, tw_bf *bf, void *evt, tw_lp *thislp);
   static void commit_event(void *state, tw_bf *bf, void *evt, tw_lp *thislp);
   static void final_lp(void *state, tw_lp *me);
   static MethodPtr forwardDispatch[];
   static MethodPtr backwardDispatch[];
   static MethodPtr commitDispatch[];

   template <typename SlotType>
   void forwardWithCheckpoint(Time tnow, typename SlotType::MsgType& msg, tw_bf* twbf) {
      InputCheckpoint ic;
      static_cast<Dependent*>(this)->Dependent::template checkpointFull<SlotType>(ic, tnow, msg, twbf);
      static_cast<Dependent*>(this)->Dependent::template forwardFull<SlotType>(tnow, msg, twbf);
      checkpointQueue.push_back(ic.str());
   }
   template <typename SlotType>
   void backwardWithCheckpoint(Time tnow, typename SlotType::MsgType& msg, tw_bf* twbf) {
      static_cast<Dependent*>(this)->Dependent::template backwardFull<SlotType>(tnow, msg, twbf);
      OutputCheckpoint oc(checkpointQueue.back());
      static_cast<Dependent*>(this)->Dependent::template checkpointFull<SlotType>(oc, tnow, msg, twbf);
      checkpointQueue.pop_back();
   }
   template <typename SlotType>
   void commitWithCheckpoint(Time tnow, typename SlotType::MsgType& msg, tw_bf* twbf) {
      static_cast<Dependent*>(this)->Dependent::template commitFull<SlotType>(tnow, msg, twbf);
      checkpointQueue.pop_front();
   }
};

#define ACTOR_DEF()                                                     \
   template<typename SlotType>                                          \
   void forwardFull(Time tnow, typename SlotType::MsgType& msg, tw_bf*);    \
   template<typename SlotType>                                          \
   void backwardFull(Time tnow, typename SlotType::MsgType& msg, tw_bf*) {} \
   template<typename SlotType>                                          \
   void commitFull(Time tnow, typename SlotType::MsgType& msg, tw_bf*) {} \
   template <typename SlotType, typename Archiver>                                        \
   inline void checkpointFull(Archiver& ar, Time tnow, typename SlotType::MsgType& msg, tw_bf*);


#define TAGIFY(SlotType) SlotType::tag

#define forward(SlotType) forward_##SlotType
#define backward(SlotType) backward_##SlotType
#define commit(SlotType) commit_##SlotType
#define checkpoint(SlotType) checkpoint_##SlotType

#define SLOT(name, msgType)                     \
class name {                                    \
 public:                                        \
   typedef msgType MsgType;                     \
   static Tag tag;                              \
}

#define SIGNAL(name, type)                                              \
   class name##_signal_class { public: typedef type MsgType; };         \
   void name(Time trecv, type& msg) {                                   \
      /*for all LPs attached to this signal*/                           \
      for (auto&& observerData : name##_listeners) {                    \
         my_event_send(observerData.dest, trecv, twlp, observerData.slot, msg);      \
      }                                                                 \
   }                                                                    \
   std::list<BaseLP::ObserverData> name##_listeners

#define CONNECT(sendingLP, signal, dest, slot) signal_add_listener<signal##_signal_class,slot>((sendingLP).signal##_listeners,dest)

template <typename SignalType, typename SlotType>
inline
typename std::enable_if<std::is_same<typename SignalType::MsgType,typename SlotType::MsgType>::value>::type
signal_add_listener(std::list<BaseLP::ObserverData>& listeners, const BaseLP& dest) {
   tw_lpid destID = dest.getID();
   signal_add_listener<SignalType,SlotType>(listeners, destID);
}

template <typename SignalType, typename SlotType>
inline
typename std::enable_if<std::is_same<typename SignalType::MsgType,typename SlotType::MsgType>::value>::type
signal_add_listener(std::list<BaseLP::ObserverData>& listeners, const tw_lpid destID) {
   listeners.push_back({destID,TAGIFY(SlotType)});
}

#define EVENT_SEND(lpid, slotClass, msgTime, msg) my_event_send<slotClass>(lpid, msgTime, twlp, msg)

#define EVENT_CANCEL(event) my_event_cancel(event)

template <typename SlotType>
PdesEvent my_event_send(tw_lpid dest, Time trecv, tw_lp* twlp, typename SlotType::MsgType& msg)
{
   return my_event_send(dest, trecv, twlp, TAGIFY(SlotType), msg);
}

template <typename TTT>
PdesEvent my_event_send(tw_lpid dest, Time trecv, tw_lp* twlp, Tag tag, TTT& msg);

void my_event_cancel(PdesEvent evtptr);

std::size_t messageSize();
