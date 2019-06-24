//#include "reaction_network/reaction.hpp"

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

template <typename VD>
inline void Reaction<VD>::set_rate_inputs(const std::map<std::string, VD>& species_linked)
{
  m_rate_inputs = ReactionBase::interpret_species_name(this->get_rate_formula(), species_linked);
}

template <typename VD>
inline const typename Reaction<VD>::rate_input_t& Reaction<VD>::get_rate_inputs() const
{
  return m_rate_inputs;
}

/**@}*/
} // end of namespace wcs
