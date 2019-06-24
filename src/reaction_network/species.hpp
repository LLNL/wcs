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

  bool inc_count();
  bool dec_count();
  bool set_count(species_cnt_t c);
  species_cnt_t get_count() const;

 protected:
  void reset(Species& obj);

 private:
  Species* clone_impl() const override;

 protected:
  unsigned m_count; ///< copy number of the species
};

/**@}*/
} // end of namespace wcs
#endif // __WCS_SPECIES_HPP__
