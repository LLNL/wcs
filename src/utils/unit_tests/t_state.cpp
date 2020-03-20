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

namespace wcs {

template<typename ObjT,
         typename CharT = char,
         typename Traits = std::char_traits<CharT> >
void test_iostremabuf(ObjT& obj)
{
  std::cout << "TEST iostream" << std::endl;
  std::vector<CharT> buffer;
  streamvec<CharT, Traits> strmbuf(buffer);
  std::basic_iostream<CharT, Traits> ss(&strmbuf);
  for (int i = 1; i < 11; i++) {
    ObjT obj2;
    ss << bits(obj);
    ss >> bits(obj2);
    strmbuf.print(std::cout, true);
  }
}

} // end of namespace wcs

using namespace wcs;


int main(int argc, char** argv)
{
  using generator_t = std::mt19937;
  //using generator_t = std::minstd_rand;
  //using generator_t = std::ranlux48_base;
  using char_t = char;

  unsigned num_gen = 1u;
  unsigned idx_sr = 0u;
  unsigned num_sr = 0u;
  auto seedval = std::chrono::system_clock::now().time_since_epoch().count();

  if (argc == 1) {
    std::cout << "Usage: > " << argv[0] << " num_gen [seed [iter_ckpt [num_save/load]]]" << std::endl;
    return 0;
  }
  if (argc > 1) {
    num_gen = static_cast<unsigned>(atoi(argv[1]));
  }
  if (argc > 2) {
    seedval = atoi(argv[2]);
  }
  if (argc > 3) {
    idx_sr = static_cast<unsigned>(atoi(argv[3]));
    num_sr = 1u;
  }
  if (argc > 4) {
    num_sr = static_cast<unsigned>(atoi(argv[4]));
  }
  std::cout << argv[0] << ' ' << num_gen << ' '  << seedval << ' ' << idx_sr << ' ' << num_sr << std::endl;

  std::seed_seq sseq{seedval};
  generator_t gen;
  gen.seed(sseq);

  std::vector<char_t> buffer(32, 'c');

  for (unsigned i = 0u; (i < idx_sr) && (i < num_gen); i++) {
    std::cout << '\t' << i << '\t' << gen() << std::endl;
  }

  const auto t1 = get_time();
  if ((idx_sr < num_gen) && ( 0u < num_sr)) {
    std::cout << "\trepeat save/load " << num_sr << " times" << std::endl;
    for (unsigned int j = 0u; j < num_sr; j++) {
      save_state(gen, buffer);
      load_state(gen, buffer);
    }
  }
  const auto t2 = get_time();

  for (unsigned i = idx_sr; i < num_gen; i++) {
    std::cout << '\t' << i << '\t' << gen() << std::endl;
  }

  if (num_sr > 0u) {
    std::cout << "Time for state save and load: " << t2 - t1 << std::endl;
    std::cout << "state size: " << buffer.size() * sizeof(char_t) << std::endl;
  }

  {
    int a = 5;
    std::vector<int> b(3);
    bits_t<int&> ba(bits(a)); // works
    (void) ba;
    //bits_t bb(bits(b)); // does not work because bits() is not defined for std::vector
  }


  std::cout << "----------------" << std::endl;
  {
    std::vector<char_t> buf;
  #if 1
    uint32_t c = 0x01020304;
    uint32_t d = 0;
  #else
    uint64_t c = 0x0102030405060708;
    uint64_t d = 0;
  #endif

    std::stringstream ss;
    ss << "saving state: 0x" << std::hex << c << std::endl;
    bool ok1 = save_state(c, buf);
    bool ok2 = ok1 && load_state(d, buf);
    if (!ok1) {
      std::cout << "Fail to save state!" << std::endl;
    } else if (!ok2) {
      std::cout << "Fail to load state!" << std::endl;
    }
    ss << "restored state: 0x" << std::hex << d << std::endl;
    std::cout << ss.str() << std::endl;
    std::cout << "state size: " << buf.size() * sizeof(char_t) << std::endl;
  }

  std::cout << "----------------" << std::endl;
  {
    uint32_t c = 5;
    test_iostremabuf(c);
  }
  return 0;
}
