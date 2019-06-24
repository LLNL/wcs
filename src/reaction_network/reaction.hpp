#ifndef __WCS_REACTION_HPP__
#define __WCS_REACTION_HPP__
#include "reaction_network/reaction_base.hpp"

namespace wcs {
/** \addtogroup wcs_reaction_network
 *  *  @{ */

template <typename VD>
class Reaction : public ReactionBase {
 public:
  using rate_input_t = std::vector<VD>;

  Reaction();
  Reaction(const Reaction& rhs);
  Reaction(Reaction&& rhs) noexcept;
  Reaction& operator=(const Reaction& rhs);
  Reaction& operator=(Reaction&& rhs) noexcept;
  ~Reaction() override;
  std::unique_ptr<Reaction> clone() const;

  void set_rate_inputs(const std::map<std::string, VD>& species_linked);
  const rate_input_t& get_rate_inputs() const;

 protected:
  void reset(Reaction& obj);

 private:
  Reaction* clone_impl() const override;

 protected:
  /// The BGL descriptors of input vertices to reaction rate formula
  rate_input_t m_rate_inputs;
};

/**@}*/
} // end of namespace wcs

#include "reaction_network/reaction_impl.hpp"
#endif // __WCS_REACTION_HPP__
