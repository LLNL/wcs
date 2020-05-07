/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#include <vector>
#include <random>
#include <cstring>
#include <sstream>

#include <chrono>
inline double get_time() {
  using namespace std::chrono;
  return duration_cast<duration<double>>(
           steady_clock::now().time_since_epoch()).count();
}

#include "utils/state_io.hpp"

using namespace wcs;

struct params_t {
  params_t()
  : num_gen(1u), idx_sr(0u), num_sr(0u),
    seedval(std::chrono::system_clock::now().time_since_epoch().count()) {}

  unsigned num_gen;
  unsigned idx_sr;
  unsigned num_sr;
  int seedval;
};


params_t get_params(int &argc, char** &argv)
{
  params_t p;

  if (argc > 1) {
    p.num_gen = static_cast<unsigned>(atoi(argv[1]));
  }
  if (argc > 2) {
    p.seedval = atoi(argv[2]);
  }
  if (argc > 3) {
    p.idx_sr = static_cast<unsigned>(atoi(argv[3]));
    p.num_sr = 1u;
  }
  if (argc > 4) {
    p.num_sr = static_cast<unsigned>(atoi(argv[4]));
  }
  std::cout << argv[0] << ' ' << p.num_gen << ' '
            << p.seedval << ' ' << p.idx_sr << ' '
            << p.num_sr << std::endl;
  return p;
}


void test_rng_state(const params_t& p)
{
  using generator_t = std::mt19937;
  //using generator_t = std::minstd_rand;
  //using generator_t = std::ranlux48_base;
  using char_t = char;

  std::seed_seq sseq{p.seedval};
  generator_t gen;
  gen.seed(sseq);

  std::vector<char_t> buffer(32, 'c');

  // generate RNs
  for (unsigned i = 0u; (i < p.idx_sr) && (i < p.num_gen); i++) {
    std::cout << '\t' << i << '\t' << gen() << std::endl;
  }

  const auto t1 = get_time();
  if ((p.idx_sr < p.num_gen) && ( 0u < p.num_sr)) {
    std::cout << "\trepeat save/load " << p.num_sr << " times" << std::endl;
    for (unsigned int j = 0u; j < p.num_sr; j++) {
      save_state(gen, buffer); // save RNG state
      load_state(gen, buffer); // load RNG state
    }
  }
  const auto t2 = get_time();

  // generate RNs
  for (unsigned i = p.idx_sr; i < p.num_gen; i++) {
    std::cout << '\t' << i << '\t' << gen() << std::endl;
  }

  if (p.num_sr > 0u) {
    std::cout << "Time for state save and load: " << t2 - t1 << std::endl;
    std::cout << "state size: " << buffer.size() * sizeof(char_t) << std::endl;
  }
}


void test_streamvec()
{
  using char_t = char;
  std::vector<char_t> buf;
#if 1
  uint32_t c = 0x01020304;
  uint32_t d = 0;
#else
  uint64_t c = 0x0102030405060708;
  uint64_t d = 0;
#endif

  std::stringstream ss;
  ss << "Saving state: 0x" << std::hex << c << std::endl;

  bool ok1 = save_state(c, buf); // save state
  bool ok2 = ok1 && load_state(d, buf); // load state

  if (!ok1) {
    std::cout << "Fail to save state!" << std::endl;
  } else if (!ok2) {
    std::cout << "Fail to load state!" << std::endl;
  }

  ss << "Restored state: 0x" << std::hex << d << std::endl;
  std::cout << ss.str() << std::endl
            << "State size: " << buf.size() * sizeof(char_t) << std::endl;
}


template<typename ObjT,
         typename CharT = char,
         typename Traits = std::char_traits<CharT> >
void test_iostremavec(ObjT& obj)
{
  const int cnt = 10;
  std::vector<CharT> buffer;
  wcs::streamvec<CharT, Traits> strmbuf(buffer);
  std::basic_iostream<CharT, Traits> ss(&strmbuf);

  // state size grows as data accumulate
  for (int i = 1; i <= cnt; i++) {
    ObjT obj2;
    ss << bits(obj);
    ss >> bits(obj2);
    strmbuf.print(std::cout, true);
  }
}


void test_streambuff()
{
  using char_t = char;
  std::vector<char_t> buf;
#if 1
  uint32_t c = 0x01020304;
  uint32_t d = 0;
#else
  uint64_t c = 0x0102030405060708;
  uint64_t d = 0;
#endif
  buf.resize(sizeof(c));

  std::stringstream ss;
  ss << "Saving state: 0x" << std::hex << c << std::endl;

  bool ok1 = save_state(c, buf.data()); // save state
  bool ok2 = ok1 && load_state(d, buf.data()); // load state

  if (!ok1) {
    std::cout << "Fail to save state!" << std::endl;
  } else if (!ok2) {
    std::cout << "Fail to load state!" << std::endl;
  }

  ss << "Restored state: 0x" << std::hex << d << std::endl;
  std::cout << ss.str() << std::endl
            << "State size: " << buf.size() * sizeof(char_t) << std::endl;
}


template<typename ObjT,
         typename CharT = char,
         typename Traits = std::char_traits<CharT> >
void test_iostremabuff(ObjT& obj)
{
  const int cnt = 10;
  std::vector<CharT> buffer;
  buffer.resize(cnt*sizeof(obj));
  wcs::streambuff<CharT, Traits> strmbuf(buffer.data(), buffer.size());
  std::basic_iostream<CharT, Traits> ss(&strmbuf);


  // state size grows as data accumulate
  for (int i = 1; i <= cnt; i++) {
    ObjT obj2;
    ss << bits(obj);
    ss >> bits(obj2);
    strmbuf.print(std::cout, true);
  }
}


int main(int argc, char** argv)
{
  if (argc == 1) {
    std::cout << "Usage: > " << argv[0] << " num_gen [seed [iter_ckpt [num_save/load]]]" << std::endl;
    return 0;
  }

  { // compilation test
    int a = 5;
    std::vector<int> b(3);
    bits_t<int&> ba(bits(a)); // works
    (void) ba;
    //bits_t bb(bits(b)); // does not work because bits() is not defined for std::vector
  }

  const params_t p = get_params(argc, argv);
  std::cout << "-------- test binary I/O for rng state --------" << std::endl;
  test_rng_state(p);

  std::cout << "-------- test streamvec (separate for input and output) --------" << std::endl;
  test_streamvec();

  std::cout << "-------- test streamvec (shared by input and output) --------" << std::endl;
  int c = 0x04030201;
  test_iostremavec(c);

  std::cout << "-------- test streambuff (separate for input and output) --------" << std::endl;
  test_streambuff();

  std::cout << "-------- test streambuff (shared by input and output) --------" << std::endl;
  //int c = 0x04030201;
  test_iostremabuff(c);

  return 0;
}
