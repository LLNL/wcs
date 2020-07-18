/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#ifndef __WCS_REACTION_NETWORK_REACTION_BASE_HPP__
#define __WCS_REACTION_NETWORK_REACTION_BASE_HPP__
#include "reaction_network/vertex_property_base.hpp"
#include "utils/exception.hpp"
#include <vector>
#include <string>
#include <functional>
#include <regex>
#include <algorithm>
#include <cctype>
#include <map>

namespace wcs {
/** \addtogroup wcs_reaction_network
 *  @{ */

class ReactionBase : public VertexPropertyBase {
 public:
  ReactionBase();
  ReactionBase(const ReactionBase& rhs);
  ReactionBase(ReactionBase&& rhs) noexcept;
  ReactionBase& operator=(const ReactionBase& rhs);
  ReactionBase& operator=(ReactionBase&& rhs) noexcept;
  ~ReactionBase() override = 0;
  std::unique_ptr<ReactionBase> clone() const;

  void set_rate_constant(reaction_rate_t k);
  reaction_rate_t get_rate_constant() const;
  /**
   * Set a default method to calculate the rate of reaction, which is the
   * product of all the reactant populations and the reaction constant.
   */
  void set_calc_rate_fn();
  /**
   * Allow setting rate calculation method specific to each reaction.
   * This is a fallback mechanism where no parsing method is available.
   * Optinoal parsing methods include ExprTk and SBML.
   */
  void set_calc_rate_fn(const std::function<
                          reaction_rate_t (const std::vector<reaction_rate_t>&)
                        >& calc_rate);
  reaction_rate_t get_rate() const;
  void set_rate_formula(const std::string& f);
  const std::string& get_rate_formula() const;

  virtual reaction_rate_t calc_rate(std::vector<reaction_rate_t> params);

 protected:
  void reset(ReactionBase& obj);
  /**
    * Parse out the species names from the given reaction formula string,
    * and translate those into the information that can identify the
    * corresponding vertices of the graph data structure using the map given.
    * The template parameter RD represents the reaction driving species.
    * i.e., std::pair<v_desc_t, stoic_t>.
    */
  template <typename RD>
  static std::vector<RD> interpret_species_name(
                           const std::string& formula,
                           const std::map<std::string, RD>& species_linked);

 private:
  ReactionBase* clone_impl() const override = 0;

 protected:
  reaction_rate_t m_rate; ///< reaction rate
  reaction_rate_t m_rate_const; ///< rate constant
  std::string m_rate_formula; ///< reaction rate formula
  std::function<reaction_rate_t (const std::vector<reaction_rate_t>&)> m_calc_rate;
};


template <typename RD>
std::vector<RD> ReactionBase::interpret_species_name(
  const std::string& formula,
  const std::map<std::string, RD>& species_linked)
{
  // concentration pattern, e.g., [A] and [B] and [C] in [A] + [B] -> [C]
  std::regex conc_pattern("\\[([\\s\\n\\r\\t]*\\w+[\\s\\n\\r\\t]*)]");
  std::regex space_filter("\\S");

  // find species names in the formula by matching against the concentration pattern
  std::sregex_iterator i
   = std::sregex_iterator(formula.begin(), formula.end(), conc_pattern);

  typename std::vector<RD> vertices;

  for( ; i != std::sregex_iterator(); ++i) {
    std::smatch m = *i;
    std::string species = m[1].str();
    // remove whitespace in each match such that, for example,
    // it gets A even with [A ] or [ A ].
    species.erase(std::remove_if(species.begin(), species.end(), ::isspace),
                                 species.end());
    typename std::map<std::string, RD>::const_iterator it
      = species_linked.find(species);
    if (it == species_linked.cend()) {
      WCS_THROW("Cannot interpret " + species + " in " + formula);
    }
    vertices.push_back(it->second);
  }
  return vertices;
}

/**@}*/
} // end of namespace wcs
#endif // __WCS_REACTION_NETWORK_REACTION_BASE_HPP__
