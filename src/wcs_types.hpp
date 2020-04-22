#ifndef __WCS_TYPES_HPP__
#define __WCS_TYPES_HPP__
#include <iostream>
#include <utility> // std::move
#include <memory>  // std::unique_ptr

namespace wcs {

using species_cnt_t = unsigned int;
using species_cnt_diff_t = int;
using reaction_rate_t = double;
using sim_time_t = double;
using sim_iter_t = unsigned;
using stoic_t = int;

template <typename G>
std::ostream& write_graphviz_of_any_vertex_list(std::ostream& os, const G& g);

class GraphFactory;

constexpr size_t num_in_edges_to_reserve = 8ul;
/// ceil(log10(2^sizeof(species_cnt_t))) + sizeof('\t')
constexpr size_t cnt_digits = 21ul;

} // end of namespace wcs
#endif // __WCS_TYPES_HPP__
