#!/usr/bin/env python

import re
from collections import OrderedDict

class Slot:
    def __init__(self, name, msg):
        self.name = name
        self.msg = msg
        self.modes = set()
    def hasMode(self, mode):
        return mode in self.modes
    def addMode(self, mode):
        self.modes.add(mode)
        
class Actor:
    def __init__(self, name):
        self.name = name
        self.slots = {}
    def addSlot(self,slot):
        self.slots[slot.name] = slot
    def getSlot(self,slotName):
        return self.slots[slotName]
    def getSlots(self):
        return self.slots.values()

def main():
    import sys

    filenames = sys.argv[1:]

    actors = OrderedDict()
    messages = set()
    for filename in filenames:
        for line in open(filename, "r"):
            mmm = re.search(r'public LP<([A-Za-z0-9_]+)>',line)
            if mmm:
                actorName = mmm.group(1)
                actors[actorName] = Actor(actorName)
            mmm = re.search(r'SLOT\(\s*([A-Za-z0-9_]+)\s*,\s*([A-Za-z0-9_]+)\s*\)',line)
            if mmm:
                slotName = mmm.group(1)
                msgName = mmm.group(2)
                actors[actorName].addSlot(Slot(slotName,msgName))
            for eventType in ['FORWARD','BACKWARD','COMMIT','CHECKPOINT']:
                mmm = re.search("void "+eventType+r'\(([A-Za-z0-9_]+)\)',line)
                if mmm:
                    slotName = mmm.group(1)
                    actors[actorName].getSlot(slotName).addMode(eventType)

    print('#include "actor.hh"')
    print('#include "Checkpoint.hh"')
    
    for filename in filenames:
        print('#include "%s"' % filename)

    print('''
struct EventData {
   Tag tag;
   union {
      char msg;''')
    allMessages = set()
    for actor in actors.values():
        for slot in actor.getSlots():
            allMessages.add(slot.msg)
    ii=0
    for msg in allMessages:
        print("      %s msg%d;" % (msg, ii))
        ii += 1
    print('''
   };
};

std::size_t messageSize() {
   return sizeof(EventData);
}

CHECKPOINTQueue checkpointQueue;

template <typename TTT>
PdesEvent my_event_send(tw_lpid dest, Time trecv, tw_lp* twlp, Tag tag, TTT& msg)
{
   tw_stime recv_ts = TW_STIME_ADD(tw_now(twlp), trecv);
   if (recv_ts >= g_tw_ts_end) {
      return NULL;
   }
   tw_event* evtptr = tw_event_new(dest, trecv, twlp);
   EventData* data = static_cast<EventData *>(tw_event_data(evtptr));
   data->tag = tag;
   memcpy(&data->msg, &msg, sizeof(msg));
   tw_event_send(evtptr);
   return evtptr;
}''')
    for msg in allMessages:
        print("template PdesEvent my_event_send(tw_lpid dest, Time trecv, tw_lp* twlp, Tag tag, %s& msg);" % msg)
    for actor in actors.values():
        for slot in actor.getSlots():
            ttt = {
                "actor" : actor.name,
                "slot" : slot.name,
                "msg" : slot.msg,
            }
            for eventType in "FORWARD","BACKWARD","COMMIT":
                ttt['event'] = eventType
                if slot.hasMode(eventType):
                    print("""template <>
inline void %(actor)s::%(event)sFull<%(actor)s::%(slot)s>(Time tnow, %(msg)s& msg, tw_bf* twbf) { %(event)s(%(slot)s)(tnow,msg,twbf); }""" % ttt)
            if slot.hasMode('CHECKPOINT'):
                print("""template <>
inline void %(actor)s::CHECKPOINTFull<%(actor)s::%(slot)s>(InputCHECKPOINT& ar, Time tnow, %(msg)s& msg, tw_bf* twbf) { CHECKPOINT(%(slot)s)(ar,tnow,msg,twbf); }
template <>
inline void %(actor)s::CHECKPOINTFull<%(actor)s::%(slot)s>(OutputCHECKPOINT& ar, Time tnow, %(msg)s& msg, tw_bf* twbf) { CHECKPOINT(%(slot)s)(ar,tnow,msg,twbf); }
""" % ttt)
    for actor in actors.values():
        for (ii,slot) in enumerate(actor.getSlots()):
            print("Tag %s::%s::tag = %d;" % (actor.name,slot.name,ii))
            
    for actor in actors.values():
        for eventType in ["FORWARD","BACKWARD","COMMIT"]:
            print("template<> %s::MethodPtr LP<%s>::%sDispatch[] = {" % (actor.name, actor.name,eventType))
            for slot in actor.getSlots():
                ttt = {
                    "actor" : actor.name,
                    "slot" : slot.name,
                    "event" : eventType,
                }
                if slot.hasMode('CHECKPOINT'):
                    method = '%(actor)s::%(event)sWithCHECKPOINT<%(actor)s::%(slot)s>' % ttt
                else:
                    method = '%(actor)s::%(event)sFull<%(actor)s::%(slot)s>' % ttt
                ttt['method'] = method
                print("  reinterpret_cast<%(actor)s::MethodPtr>(&%(method)s)," % ttt)
            print("0};")

    print("""
///This stuff needs EventData to be defined, but doesn't need specific generation

template <typename Dependent>
void LP<Dependent>::init_lp(void *state, tw_lp *thislp) {
   Dependent* derivedObj = new (reinterpret_cast<Dependent *>(state)) Dependent();
   BaseLP *lpptr = static_cast<BaseLP*>(derivedObj);
   lpptr->twlp = thislp;
}

template <typename Dependent>
void LP<Dependent>::FORWARD_event(void *state, tw_bf *bf, void *evt, tw_lp *thislp) {
   Time tnow = tw_now(thislp);
   EventData *e = reinterpret_cast<EventData *>(evt);
   Dependent* derivedObj = reinterpret_cast<Dependent*>(state);
   (derivedObj->*FORWARDDispatch[e->tag])(tnow, &e->msg, bf);
}

template <typename Dependent>
void LP<Dependent>::BACKWARD_event(void *state, tw_bf *bf, void *evt, tw_lp *thislp) {
   Time tnow = tw_now(thislp);
   EventData *e = reinterpret_cast<EventData *>(evt);
   Dependent* derivedObj = reinterpret_cast<Dependent*>(state);
   (derivedObj->*BACKWARDDispatch[e->tag])(tnow, &e->msg, bf);
}

template <typename Dependent>
void LP<Dependent>::COMMIT_event(void *state, tw_bf *bf, void *evt, tw_lp *thislp) {
   Time tnow = tw_now(thislp);
   EventData *e = reinterpret_cast<EventData *>(evt);
   Dependent* derivedObj = reinterpret_cast<Dependent*>(state);
   (derivedObj->*COMMITDispatch[e->tag])(tnow, &e->msg, bf);
}  

template <typename Dependent>
void LP<Dependent>::final_lp(void *state, tw_lp *) {
   Dependent* derivedObj = reinterpret_cast<Dependent *>(state);
   derivedObj->~Dependent();
}

""")
    for actor in actors.values():
        print("""
template void LP<%(actor)s>::FORWARD_event(void *state, tw_bf *bf, void *evt, tw_lp *thislp);
template void LP<%(actor)s>::BACKWARD_event(void *state, tw_bf *bf, void *evt, tw_lp *thislp);
template void LP<%(actor)s>::COMMIT_event(void *state, tw_bf *bf, void *evt, tw_lp *thislp);
template void LP<%(actor)s>::final_lp(void *state, tw_lp *me);
template void LP<%(actor)s>::init_lp(void *state, tw_lp *me);

""" % {"actor" : actor.name})
                

if __name__=="__main__":
    main()
