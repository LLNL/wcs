
#include "actor.hh"
#include <iostream>

void my_event_cancel(PdesEvent evtptr) {
   if (evtptr) {
      tw_event_rescind(evtptr);
   }
}

tw_lpid BaseLP::getID() const {
   return twlp->gid;
}
