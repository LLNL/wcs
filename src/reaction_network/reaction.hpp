/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#ifndef __WCS_REACTION_NETWORK_REACTION_HPP__
#define __WCS_REACTION_NETWORK_REACTION_HPP__
#if defined(WCS_HAS_CONFIG)
#include "wcs_config.hpp"
#endif

#include "reaction_network/reaction_base.hpp"

#if defined(WCS_HAS_EXPRTK)
#include "exprtk.hpp"
#endif

namespace wcs {
/** \addtogroup wcs_reaction_network
 *  @{ */

template <typename VD>
class Reaction : public ReactionBase {
 public:
  using rdriver_t = std::template pair<VD, stoic_t>;
  using involved_species_t = std::template vector<rdriver_t>;

  Reaction();
  Reaction(const Reaction& rhs);
  Reaction(Reaction&& rhs) noexcept;
  Reaction& operator=(const Reaction& rhs);
  Reaction& operator=(Reaction&& rhs) noexcept;
  ~Reaction() override;
  std::unique_ptr<Reaction> clone() const;

  const involved_species_t& get_rate_inputs() const;
  void set_products(const std::map<std::string, rdriver_t>& products);
  reaction_rate_t calc_rate(std::vector<reaction_rate_t>&& params) override;

#if defined(WCS_HAS_EXPRTK)
  void set_rate_inputs(const std::map<std::string, rdriver_t>& species_involved);
  void show_compile_error() const;
  bool detect_composite() const;
#else
  void set_rate_inputs(const std::map<std::string, rdriver_t>& species_involved,
    std::vector<std::string>& dep_params_f,
    std::vector<std::string>& dep_params_nf);
#endif // defined(WCS_HAS_EXPRTK)

 protected:
  void reset(Reaction& obj);

 private:
  Reaction* clone_impl() const override;

  std::vector<reaction_rate_t> m_params;
  bool m_is_composite;
#if defined(WCS_HAS_EXPRTK)
  exprtk::symbol_table<reaction_rate_t> m_sym_table;
  exprtk::parser<reaction_rate_t> m_parser;
  exprtk::expression<reaction_rate_t> m_expr;
#endif // defined(WCS_HAS_EXPRTK)

 protected:
  /// The BGL descriptors of input vertices to reaction rate formula
  involved_species_t m_rate_inputs;
  involved_species_t m_products;
};

/**@}*/
} // end of namespace wcs

#include "reaction_network/reaction_impl.hpp"
#endif // __WCS_REACTION_NETWORK_REACTION_HPP__
