#include "utils/rngen.hpp"
#include <iostream>
#include <limits>

int main(int argc, char** argv)
{
  using namespace std;
  if (argc > 2) {
    cout << "usage: > " << argv[0] << " [seed]" << endl;
    return 0;
  }

  using rng1_t = wcs::RNGen<std::uniform_real_distribution>;
  rng1_t r1;
  
  if (argc == 1)
    r1.set_seed();
  else
    r1.set_seed(static_cast<unsigned>(atoi(argv[1])));

  r1.param(typename rng1_t::param_type(0.0, 1.0));
  for(int i = 0; i < 10 ; ++i) {
    cout << ' ' << r1();
  }
  cout << endl;

  using rng2_t = wcs::RNGen<std::uniform_int_distribution, unsigned>;
  rng2_t r2;
  
  if (argc == 1)
    r2.set_seed();
  else
    r2.set_seed(static_cast<unsigned>(atoi(argv[1])));

  constexpr unsigned unsigned_max = std::numeric_limits<unsigned>::max();
  r2.param(typename rng2_t::param_type(1, unsigned_max-1));
  for(int i = 0; i < 10 ; ++i) {
    cout << ' ' << static_cast<double>(r2())/unsigned_max;
  }
  cout << endl;

  if (static_cast<double>(r2.distribution().min())/unsigned_max > 0.0)
    cout << "min of the distribution is greater than 0.0" << endl;
  if (static_cast<double>(r2.distribution().max())/unsigned_max < 1.0)
    cout << "max of the distribution is less than 1.0" << endl;

  return 0;
}
