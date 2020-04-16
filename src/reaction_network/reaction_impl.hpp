#include <iostream>

namespace wcs {
/** \addtogroup wcs_reaction_network
 *  *  @{ */

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

#if defined(WCS_HAS_SBML)
#pragma message ("libSBML is not supported yet.")
template <typename VD>
inline void Reaction<VD>::set_rate_inputs(const std::map<std::string, rdriver_t>& species_involved)
{
  m_rate_inputs = ReactionBase::interpret_species_name(this->get_rate_formula(), species_involved);
}
#elif defined(WCS_HAS_EXPRTK)
template <typename VD>
inline void Reaction<VD>::set_rate_inputs(const std::map<std::string, rdriver_t>& species_involved)
{
  m_rate_inputs = std::vector<rdriver_t>( species_involved.size() );
  m_params = std::vector<reaction_rate_t>( species_involved.size() );
  size_t i = 0;
  for(auto &e : species_involved ) {
      std::string var_str = e.first;
      m_rate_inputs[i] = e.second;
      // Need to remap if the m_params reallocates
      m_sym_table.add_variable( var_str, m_params[i] );
      i++;
  }
  m_sym_table.add_constant("r_const", m_rate_const);
  m_sym_table.add_constants();

  m_expr.register_symbol_table(m_sym_table);
  m_parser.compile(m_rate_formula, m_expr);
}

template <typename VD>
reaction_rate_t Reaction<VD>::calc_rate(std::vector<reaction_rate_t> params)
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
  m_rate = m_expr.value();
  return m_rate;
}

template <typename VD>
inline void Reaction<VD>::set_products(const std::map<std::string, rdriver_t>& products)
{
  for( auto &e: products )
      m_products.push_back( e.second );
}
#else
template <typename VD>
inline void Reaction<VD>::set_rate_inputs(const std::map<std::string, rdriver_t>& species_involved)
{
  m_rate_inputs = ReactionBase::interpret_species_name(this->get_rate_formula(), species_involved);
}
#endif

template <typename VD>
inline const typename Reaction<VD>::involved_species_t& Reaction<VD>::get_rate_inputs() const
{
  return m_rate_inputs;
}

/**@}*/
} // end of namespace wcs
