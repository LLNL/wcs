#ifndef SIMTIME__
#define SIMTIME__

#ifndef KMCSIM_TYPES_ONLY
#include <cassert>
#include <ross.h>
#endif

struct Time {
  double t0;
  long long int bits[2];
  
  Time(double x = 0.0) : t0(x) {
    for(int i = 0; i<(int) (sizeof(bits)/sizeof(*bits)); i++)
      bits[i] = 0;
  }
  Time(const Time &t) : t0(t.t0) {
    for(int i = 0; i<(int) (sizeof(bits)/sizeof(*bits)); i++)
      bits[i] = t.bits[i];
  }
#ifndef KMCSIM_TYPES_ONLY
  Time(const tw_stime &t) : t0(t.t) {
    assert(sizeof(*this) == sizeof(t));
    assert(sizeof(bits)/sizeof(*bits) == sizeof(t.bits)/sizeof(*(t.bits)));
    for(int i = 0; i<(int) (sizeof(bits)/sizeof(*bits)); i++)
      bits[i] = 0;
  }
#endif

  int cmp(const Time &b) const {
    if(t0 < b.t0) return -1;
    if(t0 > b.t0) return  1;
    else {
      for(int i = 0; i<(int) (sizeof(bits)/sizeof(*bits)); i++)
	if(bits[i] < b.bits[i]) return -1;
	else if(bits[i] > b.bits[i]) return 1;
    }
    return 0;
  }
  bool operator<  (const Time &b) const { return (cmp(b) <  0); }
  bool operator<= (const Time &b) const { return (cmp(b) <= 0); }
  bool operator== (const Time &b) const { return (cmp(b) == 0); }
  bool operator>= (const Time &b) const { return (cmp(b) >= 0); }
  bool operator>  (const Time &b) const { return (cmp(b) >  0); }
  bool operator!= (const Time &b) const { return (cmp(b) != 0); }
};

#endif
