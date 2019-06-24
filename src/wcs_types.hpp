#ifndef __WCS_TYPES_HPP__
#define __WCS_TYPES_HPP__
#include <iostream>
#include <utility> // std::move
#include <memory>  // std::unique_ptr

namespace wcs {
using species_cnt_t = unsigned int;
using reaction_rate_t = double;

template <typename G>
std::ostream& write_graphviz_of_any_vertex_list(std::ostream& os, const G& g);

class GraphFactory;

constexpr size_t num_in_edges_to_reserve = 8ul;

} // end of namespace wcs
#endif // __WCS_TYPES_HPP__
