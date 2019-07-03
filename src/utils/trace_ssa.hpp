#ifndef	 __WCS_UTILS_TRACE_SSA_HPP__
#define	 __WCS_UTILS_TRACE_SSA_HPP__
#include "utils/trace.hpp"

namespace wcs {
/** \addtogroup wcs_utils
 *  *  @{ */

class TraceSSA : public Trace {

protected:
  void build_index_maps() override;
  std::ostream& write_header(std::ostream& os) const override;
  void count_reaction(r_desc_t r) override;
  size_t estimate_tmpstr_size() const override;
  std::ostream& print_stats(const sim_time_t sim_time,
                            const std::vector<species_cnt_t>& species,
                            const std::string rlabel,
                            std::string& tmpstr, std::ostream& os) const override;

protected:
  /// Show how many times each reaction fires
  std::vector<size_t> m_reaction_counts;
  /// Map a BGL vertex descriptor to the reaction index
  std::unordered_map<r_desc_t, size_t> m_r_id_map;
};

/**@}*/
} // end of namespace wcs
#endif // __WCS_UTILS_TRACE_SSA_HPP__
