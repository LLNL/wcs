/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#include <iostream>
#include <string>

namespace wcs {
/** \addtogroup wcs_reaction_network
 *  @{ */

template <typename VD>
inline Reaction<VD>::Reaction()
: ReactionBase()
{}

template <typename VD>
inline Reaction<VD>::Reaction(const Reaction<VD>& rhs)
: ReactionBase(rhs),
  m_rate_inputs(rhs.m_rate_inputs)
{}

template <typename VD>
inline Reaction<VD>::Reaction(Reaction<VD>&& rhs) noexcept
: ReactionBase(std::move(rhs))
{
  if (this != &rhs) {
    m_rate_inputs = std::move(rhs.m_rate_inputs);
    reset(rhs);
  }
}

template <typename VD>
inline Reaction<VD>& Reaction<VD>::operator=(const Reaction<VD>& rhs)
{
  if (this != &rhs) {
    ReactionBase::operator=(rhs);
    m_rate_inputs = rhs.m_rate_inputs;
  }
  return *this;
}

template <typename VD>
inline Reaction<VD>& Reaction<VD>::operator=(Reaction<VD>&& rhs) noexcept
{
  if (this != &rhs) {
    ReactionBase::operator=(std::move(rhs));
    m_rate_inputs = std::move(rhs.m_rate_inputs);
    reset(rhs);
  }
  return *this;
}

template <typename VD>
inline Reaction<VD>::~Reaction() {}

template <typename VD>
inline std::unique_ptr< Reaction<VD> > Reaction<VD>::clone() const
{
  return std::unique_ptr< Reaction<VD> >(this->clone_impl());
}

template <typename VD>
inline Reaction<VD>* Reaction<VD>::clone_impl() const
{
  return (new Reaction<VD>(*this));
}

template <typename VD>
inline void Reaction<VD>::reset(Reaction& obj)
{
  obj.m_rate_inputs.clear();
  ReactionBase::reset(obj);
}

#if defined(WCS_HAS_EXPRTK)
template <typename VD>
inline void Reaction<VD>::set_rate_inputs(const std::map<std::string, rdriver_t>& species_involved)
{
  const auto num_inputs = species_involved.size();
  m_rate_inputs = std::vector<rdriver_t>( num_inputs );
  m_params.resize(num_inputs);
  size_t i = 0ul;

  std::set<std::string> var_names;
  using new_species_involvedt = std::map<std::string, rdriver_t> ;
  typename new_species_involvedt::const_iterator it_species;

  std::string formula = this->get_rate_formula();
  std::string str2 ("m_rate := ");
  std::size_t found = formula.find(str2);
  formula.replace(formula.begin(), formula.begin() + found + 10,""); //remove vars before formula
  formula.replace(formula.end()-1, formula.end(),""); //remove ; from the end of the formula
  //using std::operator<<;
  for (const auto& x: species_involved) {
    var_names.insert(x.first);
    //std::cout << " " << x.first << " ,";
  }
  //std::cout << '\n';

  // put the function parameters with the order met in the formula
  static const std::regex symbol_name("[\\w_]+"); // sequence of alnum or '_'
  std::unordered_map<std::string, size_t> input_map;
  std::sregex_token_iterator p(formula.cbegin(), formula.cend(), symbol_name);
  std::sregex_token_iterator e;
  for ( ; p!=e ; ++p) {
    for (const auto& vn : var_names) {
      if (vn == *p) {
        if (input_map.count(vn) == 0u) {
          input_map.insert(std::make_pair(vn, i++));
        }
        if (input_map.size() == var_names.size()) {
          break;
        }
      }
    }
  }

  std::vector<std::string> var_names_ord (input_map.size());
  for (const auto& x: input_map) {
    var_names_ord[x.second] = x.first;
  }

  size_t j = 0ul;
  //put only species_involved are used in formula
  std::vector<std::string>::const_iterator it;
  for (it = var_names_ord.cbegin(); it < var_names_ord.cend(); it++) {
    it_species = species_involved.find(*it);
    if (it_species != species_involved.cend()) {
      std::string var_str = it_species->first;
      m_rate_inputs[j] = it_species->second;
      // Need to remap if the m_params reallocates
      m_sym_table.add_variable( var_str, m_params[j] );
      j++;
    }
  }

  //resize both
  m_rate_inputs.resize(j);
  m_params.resize(j);

  m_sym_table.add_variable("m_rate", m_rate);
  m_sym_table.add_constant("r_const", m_rate_const);
  m_sym_table.add_constants();

  m_is_composite = detect_composite();

  m_expr.register_symbol_table(m_sym_table);
  if (!m_parser.compile(m_rate_formula, m_expr)) {
    show_compile_error();
    return;
  }
}

template <typename VD>
reaction_rate_t Reaction<VD>::calc_rate(std::vector<reaction_rate_t>&& params)
{
  // GG: This copy could be avoided by directly linking species count to sym_table
  if (m_params.size() != params.size()) {
    using std::operator<<;
    WCS_THROW("The number of involved species differs from what is expected");
    // If params is larger than m_params, m_params will reallocate its data.
    // Then, the symbol table needs to reset, and re-registered.
    // If it is smaller, some values are missing, and the symbol table may
    // not make sense at all.
  }
  // The order of parameters is the same as the one in the return of
  // get_rate_inputs()
  m_params.assign(params.begin(), params.end());

  if (!m_is_composite) {
    m_rate = m_expr.value();
  } else {
    m_rate = static_cast<reaction_rate_t>(0.0);
    m_expr.value();
  }
  // Depending on the species population, reaction rate formula can evaluate
  // to a negative value. Rather than having a formula include a conditional
  // logic, we deal with it here.
  if (m_rate < static_cast<reaction_rate_t>(0.0)) {
    m_rate = static_cast<reaction_rate_t>(0.0);
  }

  return m_rate;
}

template <typename VD>
inline void Reaction<VD>::set_products(const std::map<std::string, rdriver_t>& products)
{
  for( auto &e: products )
      m_products.push_back( e.second );
}

template <typename VD>
inline void Reaction<VD>::show_compile_error () const
{
  using std::operator<<;

  std::string err =
    "Error: " + m_parser.error() + "\tExpression: " + this->m_rate_formula;
  std::cerr << err << std::endl;

  for (size_t i = 0u; i < m_parser.error_count(); ++i)
  {
     exprtk::parser_error::type error = m_parser.get_error(i);
     std::string errmsg
       = "Error: " + std::to_string(i)
       + "\tPosition: " + std::to_string(error.token.position)
       + "\tType: [" + exprtk::parser_error::to_str(error.mode) + ']'
       + "\tMsg: " + error.diagnostic
       + "\tExpression: " + m_rate_formula;
     std::cerr << errmsg << std::endl;
  }
}

template <typename VD>
inline bool Reaction<VD>::detect_composite() const
{
  static const std::string whitespaces(" \t\f\v\n\r");
  const auto pos1 = m_rate_formula.find_last_of (";");
  const auto head = m_rate_formula.substr(0, pos1);
  const auto tail = m_rate_formula.substr(pos1+1);
  const auto pos2 = tail.find_last_not_of(whitespaces);
  const auto pos3 = head.find_first_of(";");
  const auto pos4 = head.find_last_not_of(whitespaces);
  return ( pos1 != std::string::npos ) &&
         ( ( pos2 != std::string::npos ) ||
           ( ( pos3 != std::string::npos ) &&
             ( pos4 != std::string::npos ) &&
             ( pos3 < pos4)
           )
         );
}
#elif defined(WCS_HAS_SBML)
template <typename VD>
inline void Reaction<VD>::set_rate_inputs(const std::map<std::string, rdriver_t>& species_involved,
std::vector<std::string>& dep_params_rf,
std::vector<std::string>& dep_params_nrf)
{
  const auto num_inputs = species_involved.size();
  m_rate_inputs = std::vector<rdriver_t>( num_inputs );
  m_params.resize(num_inputs);
  size_t i = 0ul;
  /*for(auto &e : species_involved ) {
      std::string var_str = e.first;
      m_rate_inputs[i] = e.second;
      i++;
  } */

  std::set<std::string> var_names;
  using new_species_involvedt = std::map<std::string, rdriver_t> ;
  typename new_species_involvedt::const_iterator it_species;
  new_species_involvedt new_species_involved;

  std::string formula = this->get_rate_formula();
  std::string str2 ("m_rate := ");
  std::size_t found = formula.find(str2);
  formula.replace(formula.begin(), formula.begin() + found + 10,""); //remove vars before formula
  formula.replace(formula.end()-1, formula.end(),""); //remove ; from the end of the formula

  using std::operator<<;
  for (const auto& x: species_involved) {
    //std::cout << " " << x.first << " ,";
    var_names.insert(x.first);
  }
  //std::cout << '\n';

  // put the function parameters with the order met in the formula
  static const std::regex symbol_name("[\\w_]+"); // sequence of alnum or '_'
  std::unordered_map<std::string, size_t> input_map;
  std::sregex_token_iterator p(formula.cbegin(), formula.cend(), symbol_name);
  std::sregex_token_iterator e;
  i = 0ul;
  for ( ; p!=e ; ++p) {
    for (const auto& vn : var_names) {
      if (vn == *p) {
        if (input_map.count(vn) == 0u) {
          input_map.insert(std::make_pair(vn, i++));
        }
        if (input_map.size() == var_names.size()) {
          break;
        }
      }
    }
  }

  std::vector<std::string> var_names_ord (input_map.size());
  for (const auto& x: input_map) {
    var_names_ord[x.second] = x.first;
  }

  //std::cout << std::endl << "\nReceived:";
  //put first parameters in formula taking input
  std::vector<std::string>::const_iterator it;
  size_t j = 0ul;
  for (it = var_names_ord.cbegin(); it < var_names_ord.cend(); it++) {
    it_species = species_involved.find(*it);
    if (it_species != species_involved.end()){
      // include only species from species_involved expected as input
      if ( std::find(dep_params_rf.cbegin(), dep_params_rf.cend(), *it)
      != dep_params_rf.cend()) {
        std::string var_str = it_species->first;
        m_rate_inputs[j] = it_species->second;
        j++;
        var_names.erase(*it);
      }
    }
  }

  //put rest parameters (not in formula) taking input
  for (const auto& x: var_names) {
    it_species = species_involved.find(x);
    if (it_species != species_involved.cend()) {
      // include only species from species_involved expected as input
      if ( std::find(dep_params_nrf.cbegin(), dep_params_nrf.cend(), x)
      != dep_params_nrf.cend()) {
        std::string var_str = it_species->first;
        m_rate_inputs[j] = it_species->second;
        j++;
      }
    }
  }
  //resize both
  m_rate_inputs.resize(j);
  m_params.resize(j);

}

template <typename VD>
reaction_rate_t Reaction<VD>::calc_rate(std::vector<reaction_rate_t>&& params)
{
  // GG: This copy could be avoided by directly linking species count to sym_table
  if (m_params.size() != params.size()) {
    using std::operator<<;
    WCS_THROW("The number of involved species differs from what is expected");
    // If params is larger than m_params, m_params will reallocate its data.
    // Then, the symbol table needs to reset, and re-registered.
    // If it is smaller, some values are missing, and the symbol table may
    // not make sense at all.
  }
  // The order of parameters is the same as the one in the return of
  // get_rate_inputs()
  m_params.assign(params.begin(), params.end());

  m_rate = m_calc_rate(m_params);

  // Depending on the species population, reaction rate formula can evaluate
  // to a negative value. Rather than having a formula include a conditional
  // logic, we deal with it here.
  if (m_rate < static_cast<reaction_rate_t>(0.0)) {
    m_rate = static_cast<reaction_rate_t>(0.0);
  }

  return m_rate;
}

template <typename VD>
inline void Reaction<VD>::set_products(const std::map<std::string, rdriver_t>& products)
{
  for( auto &e: products )
      m_products.push_back( e.second );
}
#else
#error "Must enable either ExprTk or SBML"
#endif

template <typename VD>
inline const typename Reaction<VD>::involved_species_t& Reaction<VD>::get_rate_inputs() const
{
  return m_rate_inputs;
}

#ifdef WCS_CACHE_DEPENDENT
template <typename VD>
void Reaction<VD>::set_dependent_reactions(const std::set<VD>& dependent)
{
  m_dependent_reactions.assign(dependent.begin(), dependent.end());
}

template <typename VD>
const std::vector<VD>& Reaction<VD>::get_dependent_reactions() const
{
  return m_dependent_reactions;
}

template <typename VD>
void Reaction<VD>::clear_dependent_reactions()
{
  m_dependent_reactions.clear();
}
#endif // WCS_CACHE_DEPENDENT

/**@}*/
} // end of namespace wcs
