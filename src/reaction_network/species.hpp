#ifndef __WCS_SPECIES_HPP__
#define __WCS_SPECIES_HPP__

#include "reaction_network/vertex_property_base.hpp"

namespace wcs {
/** \addtogroup wcs_reaction_network
 *  *  @{ */

class Species : public VertexPropertyBase {
 public:
  Species();
  Species(const Species& rhs);
  Species(Species&& rhs) noexcept;
  Species& operator=(const Species& rhs);
  Species& operator=(Species&& rhs) noexcept;
  ~Species() override;
  std::unique_ptr<Species> clone() const;

  /// Increase the count by one. Return false under over/under-flow condition.
  bool inc_count();
  /// Decrease the count by one. Return false under over/under-flow condition.
  bool dec_count();
  /// Increase the count by the given amount c.
  bool inc_count(const species_cnt_t c);
  /// Decrease the count by the given amount c.
  bool dec_count(const species_cnt_t c);
  /// Set the count to the given amount.
  bool set_count(const species_cnt_t c);
  /// Return the current count.
  species_cnt_t get_count() const;

 protected:
  void reset(Species& obj);

 private:
  Species* clone_impl() const override;

 protected:
  species_cnt_t m_count; ///< copy number of the species
  /// The maximum count for a species allowed
  static species_cnt_t m_max_count;
};

/**@}*/
} // end of namespace wcs
#endif // __WCS_SPECIES_HPP__
