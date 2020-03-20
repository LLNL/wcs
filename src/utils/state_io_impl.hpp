#ifndef __STATE_IO_IMPL_HPP__
#define __STATE_IO_IMPL_HPP__

namespace wcs {

template<typename ObjT,
         typename CharT,
         typename Traits>
bool export_state(const ObjT& obj, std::vector<CharT>& buffer) {
/*
  if (buffer.size()*sizeof(CharT) < sizeof(obj)) {
    buffer.resize((sizeof(obj) + sizeof(CharT) - 1u)/sizeof(CharT));
  }
*/

  ostreamvec<CharT, Traits> ostrmbuf(buffer);
  //streamvec<CharT, Traits> ostrmbuf(buffer);
  std::basic_ostream<CharT, Traits> oss(&ostrmbuf);

  oss << bits(obj);

  return oss.good();
}

template<typename ObjT,
         typename CharT,
         typename Traits>
bool import_state(ObjT& obj, std::vector<CharT>& buffer) {
  istreamvec<CharT, Traits> istrmbuf(buffer);
  //streamvec<CharT, Traits> istrmbuf(buffer, true);
  std::basic_istream<CharT, Traits> iss(&istrmbuf);

  if (istrmbuf.size()*sizeof(CharT) < sizeof(ObjT)) {
    return false;
  }

  iss >> bits(obj);

  return iss.good();
}

} // end of namespace wcs

#endif // __STATE_IO_IMPL_HPP__
