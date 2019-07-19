#ifndef __WCS_REACTION_NETWORK_REACTION_HPP__
#define __WCS_REACTION_NETWORK_REACTION_HPP__
#include "reaction_network/reaction_base.hpp"
#include "include/exprtk.hpp"

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

  void set_rate_inputs(const std::map<std::string, VD>& reactants);
  void set_outputs(const std::map<std::string, VD>& products);
  const rate_input_t& get_rate_inputs() const;
  reaction_rate_t calc_rate(std::vector<reaction_rate_t> params) override;

 protected:
  void reset(Reaction& obj);

 private:
  Reaction* clone_impl() const override;

  std::vector<reaction_rate_t> m_params;
  exprtk::symbol_table<reaction_rate_t> m_sym_table;
  exprtk::parser<reaction_rate_t> m_parser;
  exprtk::expression<reaction_rate_t> m_expr;

 protected:
  /// The BGL descriptors of input vertices to reaction rate formula
  rate_input_t m_rate_inputs;
  std::vector<VD> m_outputs;
};

/**@}*/
} // end of namespace wcs

#include "reaction_network/reaction_impl.hpp"
#endif // __WCS_REACTION_NETWORK_REACTION_HPP__
